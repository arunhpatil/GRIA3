import os
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch


# ================================
# SETUP ENVIRONMENT & GPU
# ================================
scvi.settings.device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {scvi.settings.device}")

# Set paths

REF_PATH = "/dcs04/lieber/hwanglab/Arun/bmaher/Annotations/HNOCA_annotation/archmap/treutlein_neural_organoid_atlas/data.h5ad"
REF_PATH_dataOnly = "/dcs04/lieber/hwanglab/Arun/bmaher/Annotations/HNOCA_annotation/archmap/treutlein_neural_organoid_atlas/data_only_count.h5ad"
QUERY_PATH = "/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/3_analysis/02_cell_annotations/cellAtlasBased/1_0_markerBased/hnoca_adata_with_markerGenesAnnotation_090826.h5ad"
REF_MODEL_DIR = "/dcs04/lieber/hwanglab/Arun/bmaher/Annotations/HNOCA_annotation/archmap/saved_scanvi_reference_model"
# MODEL_PATH = "/dcs04/lieber/hwanglab/Arun/bmaher/Annotations/HNOCA_annotation/archmap/treutlein_neural_organoid_atlas/"
# import shutil; shutil.rmtree(REF_MODEL_DIR, ignore_errors=True) # Only to remove the wrong model. Dont execute again!

OUTPUT_DIR = "/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/3_analysis/02_cell_annotations/cellAtlasBased/1_1_refBased_neuralOrganoid/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

labels_key_ref = "annot_level_1"
# labels_key_query = "batch" # Change to your actual query label column

# ================================
# 1. LOAD DATA
# ================================
print("Loading reference and query data...")
ref_adata = sc.read_h5ad(REF_PATH)
count_ref_adata = sc.read_h5ad(REF_PATH_dataOnly)
query_adata = sc.read_h5ad(QUERY_PATH)

# Check if reference .X is empty and restore it from your counts layer
if hasattr(ref_adata.X, "nnz") and ref_adata.X.nnz == 0:
    print("Reference .X is empty. Restoring from layers['counts']...")
    # Note: Change "counts" if your raw counts are stored under a different layer name (e.g., "raw_counts")
    # ref_adata.X = ref_adata.layers["counts"].copy()

#assert ref_adata.X.nnz > 0, "Reference counts matrix (.X) is still empty! Check your layers."

ref_adata.X = count_ref_adata.X.tocsr().copy()  #The ref_adata doesn't come with counts, they should be loaded separately (also provided by the authors).

# ================================
# 2. PRESERVE FULL-GENE BACKUP & NORM COUNTS
# ================================
# Save current normalized data to a layer
query_adata.layers["norm_counts"] = query_adata.X.copy()

# Keep a full-gene copy (all 34,584 genes) for downstream biology
query_full = query_adata.copy()

# Ensure the working query uses raw counts in .X for scANVI
if "counts" in query_adata.layers:
  query_adata.X = query_adata.layers["counts"].copy()


# # 2. Strict filtering to remove empty droplets, dying cells, and doublets
# sc.pp.filter_cells(query_adata, min_genes=500)
# sc.pp.filter_cells(query_adata, min_counts=1000)
# Optional: filter out cells with exceptionally high counts (likely doublets)
# max_counts = np.percentile(query_adata.X.sum(axis=1), 99)
# query_adata = query_adata[query_adata.X.sum(axis=1) < max_counts].copy()
# 
# # Convert the 2D sum matrix into a flat 1D NumPy array
# library_sizes = query_adata.X.sum(axis=1)
# library_sizes = library_sizes.A1 if hasattr(library_sizes, "A1") else library_sizes.ravel()
# 
# # Calculate max_counts threshold
# max_counts = np.percentile(library_sizes, 99)
# 
# # Now boolean indexing will work correctly!
# query_adata = query_adata[library_sizes < max_counts].copy()

query_adata.obs["batch"] = query_adata.obs["sample"].astype(str)

# 1. Inspect your label column for hidden NaNs or weird types
# labels_key = "batch" # Replace with your actual label column name
# print("Unique labels before cleanup:", query_adata.obs[labels_key].unique())
# 
# # 2. Check for any rogue NaNs, Infs, or negative numbers inside the sparse matrix data
# if np.isnan(query_adata.X.data).any() or np.isinf(query_adata.X.data).any():
#     print("Found NaNs or Infs in sparse matrix! Cleaning...")
#     query_adata.X.data[np.isnan(query_adata.X.data)] = 0.0
#     query_adata.X.data[np.isinf(query_adata.X.data)] = 0.0
# 
# if (query_adata.X.data < 0).any():
#     print("Warning: Negative values found in count matrix. Clipping to 0...")
#     query_adata.X.data = np.clip(query_adata.X.data, 0, None)

# ================================
# 3. REFERENCE MODEL SETUP & TRAINING
# ================================
# Setup reference AnnData for scANVI
scvi.model.SCANVI.setup_anndata(
    ref_adata,
    layer=None,
    # layer="counts",  # Adjust if counts are stored in a different layer or .X
    batch_key="batch",  # Column name for reference batches
    labels_key="annot_level_1",  # Column name for reference cell type annotations
    unlabeled_category="Unknown",
)

# Create or load reference model matching current scvi-tools version
if not os.path.exists(REF_MODEL_DIR):
  print("Training reference scANVI model from scratch...")
  ref_scanvi_model = scvi.model.SCANVI(ref_adata)
  ref_scanvi_model.train(max_epochs=25)  # Adjust epochs depending on dataset size
  ref_scanvi_model.save(REF_MODEL_DIR, overwrite=True)
else:
  print("Loading existing reference scANVI model...")
  # ref_scanvi_model = scvi.model.SCANVI.load(REF_MODEL_DIR)

# ================================
# 4. PREPARE QUERY DATA (AUTO ZERO-PADDING)
# ================================
print("Preparing query data for scArches alignment...")
# Label query cells as 'Unknown' for semi-supervised transfer
query_adata.obs["annot_level_1"] = "Unknown"

# scvi.model.SCVI.setup_anndata(query_adata, layer="counts", batch_key="batch")

# Automatically aligns genes and zero-pads missing reference features
scvi.model.SCANVI.prepare_query_anndata(query_adata, REF_MODEL_DIR)


# ================================
# 5. SURGICAL QUERY TRAINING & PREDICTION
# ================================
print("Initializing query model and starting surgical training...")
query_scanvi_model = scvi.model.SCANVI.load_query_data(query_adata, REF_MODEL_DIR)

# weight_decay=0.0 prevents reference space distortion
query_scanvi_model.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0}, early_stopping=True)
# query_scanvi_model.train(
#     max_epochs=100, 
#     plan_kwargs={"weight_decay": 0.0, "lr": 1e-4}, 
#     early_stopping=True
# )
# ref = scvi.model.SCANVI.load(REF_MODEL_DIR, adata=ref_adata)
print(any(torch.isnan(p).any().item() for p in query_scanvi_model.module.parameters()))

##########
# query_adata.obs["annot_level_1"] = "Unknown"
# 
# # Align to reference genes FIRST (in place)
# scvi.model.SCANVI.prepare_query_anndata(query_adata, REF_MODEL_DIR)
# 
# # Now check/filter on the genes the model actually sees
# lib = np.asarray(query_adata.X.sum(axis=1)).ravel()
# print("zero-lib cells:", (lib == 0).sum(), "| <100 counts:", (lib < 100).sum())
# query_adata = query_adata[lib >= 200].copy()   # pick a sensible threshold
# 
# # Sanity checks
# X = query_adata.X
# vals = X.data if hasattr(X, "data") else np.asarray(X).ravel()
# assert np.isfinite(vals).all() and (vals >= 0).all()
# assert np.allclose(vals, np.round(vals)), "X is not integer counts"
# 
# query_scanvi_model = scvi.model.SCANVI.load_query_data(query_adata, REF_MODEL_DIR)
# query_scanvi_model.train(
#     max_epochs=100,
#     plan_kwargs={"weight_decay": 0.0},
#     early_stopping=True,
# )

##########

# print("X min/max:", query_adata.X.min(), query_adata.X.max())
# print("X mean/std:", query_adata.X.mean(), query_adata.X.std())
# print("Counts min/max:", query_adata.layers["counts"].min(), query_adata.layers["counts"].max())
# print("Counts mean/std:", query_adata.layers["counts"].mean(), query_adata.layers["counts"].std())
# Train with stabilized parameters
# query_scanvi_model.train(
#     max_epochs=100, 
#     plan_kwargs={
#         "weight_decay": 0.0, 
#         "lr": 1e-5
#     }, 
#     early_stopping=True,
#     precision="32",              # Passed directly, not inside a dict
#     gradient_clip_val=0.5        # Passed directly
# )

print("Extracting latent space and predicting labels...")
query_adata.obsm["X_scANVI"] = query_scanvi_model.get_latent_representation(
    query_adata
)
query_adata.obs["predicted_label"] = query_scanvi_model.predict(query_adata)

# Get soft prediction probabilities for confidence scoring
soft_preds = query_scanvi_model.predict(query_adata, soft=True)
query_adata.obs["prediction_confidence"] = soft_preds.max(axis=1)

# ================================
# 6. TRANSFER RESULTS TO FULL-GENE BACKUP
# ================================
print("Transferring outputs back to the full-gene backup object...")
query_full.obsm["X_scANVI"] = query_adata.obsm["X_scANVI"]
query_full.obs["predicted_label"] = query_adata.obs["predicted_label"]
query_full.obs["prediction_confidence"] = query_adata.obs[
    "prediction_confidence"
]

# Save final annotated full-transcriptome object
output_file = os.path.join(
    OUTPUT_DIR, "harmonized_full_genes_scanvi_annotated.h5ad"
)
query_full.write(output_file)
print(f"Successfully saved annotated object to: {output_file}")

# ================================
# END Here!!!
# ================================

# 
# 
# 
# from scvi_colab import install
# 
# install()
# 
# import os
# import tempfile
# 
# import anndata
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# import scanpy as sc
# import scvi
# import seaborn as sns
# import torch
# 
# # import anndata
# # import hnoca.map as hmap
# import hnoca.mapping.mapper as hmap
# import scipy.sparse as sp
# 
# scvi.settings.seed = 0
# print("Last run with scvi-tools version:", scvi.__version__)
# 
# # Ensure GPU is available
# print("Using GPU:", torch.cuda.is_available())
# scvi.settings.device = "cuda" if torch.cuda.is_available() else "cpu"
# 
# try:
#   # Programmatically apply IPython magic configurations if available
#   ipython = get_ipython()
#   if ipython is not None:
#     ipython.run_line_magic(
#         "config", 'InlineBackend.print_figure_kwargs={"facecolor": "w"}'
#     )
#     ipython.run_line_magic("config", 'InlineBackend.figure_format="retina"')
# except NameError:
#   # Fallback for standard Python scripts using native matplotlib
#   import matplotlib.pyplot as plt
#   plt.rcParams["figure.facecolor"] = "w"
# 
# 
# sc.set_figure_params(figsize=(6, 6), frameon=False)
# sns.set_theme()
# torch.set_float32_matmul_precision("high")
# # save_dir = tempfile.TemporaryDirectory()
# # 
# # REF_PATH = "/dcs04/lieber/hwanglab/Arun/bmaher/Annotations/HNOCA_annotation/archmap/treutlein_neural_organoid_atlas/data.h5ad"
# # QUERY_PATH = "/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/3_analysis/02_cell_annotations/cellAtlasBased/1_0_markerBased/hnoca_adata_with_markerGenesAnnotation_090826.h5ad"
# # MODEL_PATH = "/dcs04/lieber/hwanglab/Arun/bmaher/Annotations/HNOCA_annotation/archmap/treutlein_neural_organoid_atlas/"
# # 
# # save_dir = "/dcs04/lieber/hwanglab/Arun/bmaher/GRIA3_fai/3_analysis/02_cell_annotations/cellAtlasBased/1_1_refBased_neuralOrganoid/"
# # os.makedirs(save_dir, exist_ok=True)
# # 
# # MODEL_PATH = "/dcs04/lieber/hwanglab/Arun/bmaher/Annotations/HNOCA_annotation/archmap/treutlein_neural_organoid_atlas"
# # CONVERTED_MODEL_PATH = MODEL_PATH + "_converted"
# 
# # Convert the legacy model if it hasn't been converted yet
# if not os.path.exists(CONVERTED_MODEL_PATH):
#   print("Converting legacy scvi-tools model to updated format...")
#   scvi.model.SCANVI.convert_legacy_save(
#       dir_path=MODEL_PATH,
#       output_dir_path=CONVERTED_MODEL_PATH,
#       overwrite=True,
#   )
#   print("Conversion complete!")
# else:
#   print("Converted model directory already exists.")
# 
# print("Loading data...")
# ref_adata = sc.read(REF_PATH)
# query_adata = sc.read(QUERY_PATH)
# 
# print(ref_adata.obs.columns)
# print(query_adata.obs.columns)
# 
# #####
# 
# # Save current normalized/log1p data into a layer
# query_adata.layers["norm_counts"] = query_adata.X.copy()
# 
# # CRITICAL: scvi-tools requires raw counts. Ensure query uses the counts layer.
# if "counts" in query_adata.layers:
#   query_adata.X = query_adata.layers["counts"].copy()
# 
# # Keep this intact with all 34,584 genes and your "norm_counts" layer
# query_full = query_adata.copy()
# 
# # Subset query data strictly to reference genes
# # query_adata = query_adata[:, ref_adata.var_names].copy()
# 
# # ================================
# # 2. PREPARE QUERY FOR scANVI MAPPING
# # ================================
# # Assign an 'Unknown' (unlabeled) category to query cells for semi-supervised transfer
# query_adata.obs["scanvi_label"] = "Unknown"
# 
# 
# import anndata
# anndata.read = anndata.read_h5ad
# import scarches as sca
# 
# print("Loading reference scPoli model...")
# hnoca_model = sca.models.scPoli.load(
#     MODEL_PATH,
#     adata=ref_adata,
# )
# 
# 
# # ========================================================
# # ALIGN QUERY GENES TO REFERENCE (WITH ZERO-PADDING)
# # ========================================================
# print("Aligning query genes to reference (3000 genes)...")
# present_genes = np.intersect1d(query_adata.var_names, ref_adata.var_names)
# missing_genes = np.setdiff1d(ref_adata.var_names, query_adata.var_names)
# 
# # Subset query to present reference genes
# q_sub = query_adata[:, present_genes].copy()
# 
# if len(missing_genes) > 0:
#   print(
#       f"Padding {len(missing_genes)} missing reference genes with zeros..."
#   )
#   if sp.issparse(q_sub.X):
#     zero_mat = sp.csr_matrix((q_sub.n_obs, len(missing_genes)), dtype=q_sub.X.dtype)
#     new_X = sp.hstack([q_sub.X, zero_mat], format="csr")
#   else:
#     zero_mat = np.zeros((q_sub.n_obs, len(missing_genes)), dtype=q_sub.X.dtype)
#     new_X = np.hstack([q_sub.X, zero_mat])
#     
#   # Match reference var dataframe order exactly
#   temp_var = pd.concat([
#       q_sub.var.loc[present_genes],
#       pd.DataFrame(index=missing_genes, columns=q_sub.var.columns),
#   ])
#   current_order = np.concatenate([present_genes, missing_genes])
#   reorder_idx = [np.where(current_order == g)[0][0] for g in ref_adata.var_names]
#   if sp.issparse(new_X):
#     new_X = new_X.tocsr()[:, reorder_idx]
#   else:
#     new_X = new_X[:, reorder_idx]
#   query_aligned = anndata.AnnData(
#       X=new_X, obs=query_adata.obs.copy(), var=ref_adata.var.copy()
#   )
# else:
#   query_aligned = q_sub[:, ref_adata.var_names].copy()
# 
# # ================================
# # scPoli ATLAS MAPPING
# # ================================
# print("Mapping query using scPoli AtlasMapper...")
# mapper = hmap.AtlasMapper(hnoca_model)
# mapper.map_query(
#     query_aligned,
#     retrain="partial",
#     n_epochs=30,
#     pretraining_epochs=20,
#     eta=10,
#     batch_size=1024,
# )
# 
# mapper_dir = os.path.join(OUTPUT_DIR, "HNOCA_FXS_mapper/")
# mapper.save(mapper_dir)
# mapper = hmap.AtlasMapper.load(mapper_dir)
# 
# mapper.compute_wknn(k=100)
# ct_transfer = mapper.transfer_labels(label_key="annot_level_1")
# 
# # Extract latent representation from the aligned mapping
# X_latent = mapper.get_latent_representation(query_aligned)
# 
# # ========================================================
# # TRANSFER RESULTS BACK TO FULL-GENE OBJECT (query_full)
# # ========================================================
# query_full.obsm["X_latent"] = X_latent
# query_full.obs["best_label"] = ct_transfer["best_label"].values
# query_full.obs["best_score"] = ct_transfer["best_score"].values
# 
# # Save the final harmonized, full-transcriptome object
# output_h5ad = os.path.join(
#     OUTPUT_DIR, "harmonized_combined_HNOAC_annotated_090726.h5ad"
# )
# query_full.write(output_h5ad)
# print(f"Successfully saved annotated full-gene object to: {output_h5ad}")
# 
# # Load the pretrained reference scANVI model
# print("Loading reference scANVI model...")
# ref_scanvi_model = scvi.model.SCANVI.load(CONVERTED_MODEL_PATH, adata=ref_adata)
# 
# # Prepare query AnnData structures using scArches alignment
# scvi.model.SCANVI.prepare_query_anndata(query_adata, ref_scanvi_model)
# 
# # Initialize query model from reference model
# query_scanvi_model = scvi.model.SCANVI.load_query_data(
#     query_adata, ref_scanvi_model
# )
# 
# # ================================
# # 3. SURGICAL TRAINING & PREDICTION
# # ================================
# print("Training query model (surgical mapping)...")
# # weight_decay=0.0 prevents shifting/drifting the pre-trained reference space
# query_scanvi_model.train(
#     max_epochs=100,
#     plan_kwargs={"weight_decay": 0.0},
#     early_stopping=True,
# )
# # 
# # # Extract latent representation and predicted labels
# # query_adata.obsm["X_scANVI"] = query_scanvi_model.get_latent_representation(
# #     query_adata
# # )
# # query_adata.obs["predicted_label"] = query_scanvi_model.predict(query_adata)
# # 
# # # Optional: Get soft prediction probabilities for confidence scoring
# # soft_preds = query_scanvi_model.predict(query_adata, soft=True)
# # query_adata.obs["prediction_confidence"] = soft_preds.max(axis=1)
# # 
# # print("Label transfer complete!")
# 
# # Extract mapping outputs
# X_scanvi = query_scanvi_model.get_latent_representation(query_subset)
# predicted_labels = query_scanvi_model.predict(query_subset)
# soft_preds = query_scanvi_model.predict(query_subset, soft=True)
# confidence = soft_preds.max(axis=1)
# 
# # ==========================================
# # STEP C: TRANSFER RESULTS BACK TO QUERY_FULL
# # ==========================================
# # Now inject the labels, confidence scores, and latent coordinates
# # back into the full-gene AnnData object using cell alignment
# query_full.obsm["X_scANVI"] = X_scanvi
# query_full.obs["predicted_label"] = predicted_labels
# query_full.obs["prediction_confidence"] = confidence
# 
# # Save the harmonized full-gene object for your downstream analysis
# query_full.write(
#     os.path.join(OUTPUT_DIR, "harmonized_full_genes_annotated.h5ad")
# )


