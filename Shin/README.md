# GRIA3
GRIA3 knockout organoids - single cell sequencing analysis pipeline 

```
adata.obs.head()
                                    sample cell_line     genome genotype         pair   batch sample_batch source_batch
Batch1_Koala_2_AAACCCAAGAGAGAAC-1  Koala_2   Koala_2  old_Koala       WT  stem_1017BM  Batch1      Koala_2       Batch1
Batch1_Koala_2_AAACCCAAGAGCCCAA-1  Koala_2   Koala_2  old_Koala       WT  stem_1017BM  Batch1      Koala_2       Batch1
Batch1_Koala_2_AAACCCAAGATGCTGG-1  Koala_2   Koala_2  old_Koala       WT  stem_1017BM  Batch1      Koala_2       Batch1
Batch1_Koala_2_AAACCCAAGCATCAAA-1  Koala_2   Koala_2  old_Koala       WT  stem_1017BM  Batch1      Koala_2       Batch1
Batch1_Koala_2_AAACCCAAGGAGAATG-1  Koala_2   Koala_2  old_Koala       WT  stem_1017BM  Batch1      Koala_2       Batch1
```

```
adata.obs['cell_line'].value_counts()
cell_line
1017_G3-8    20544
4001_NC-1    19103
4002_NC-1    15285
4002_G3-2    13967
1017_NC-1    13741
4001_G3-1    12822
Lemur        12448
Sloth        12412
Whale        10397
Koala_2      10186
Otter         9997
Zebra         8952
Name: count, dtype: int64
```

```
adata.obs[['batch']].value_counts()
batch 
Batch2    95462
Batch1    64392
Name: count, dtype: int64
```

```
print(adata.obs[['sample', 'genome', 'batch','sample_batch', 'source_batch', 'pair', 'cell_line', 'genotype']].drop_duplicates().sort_values(['p\
air', 'genotype']).to_markdown(index=False,tablefmt="github"))


| sample    | genome    | batch   | sample_batch   | source_batch   | pair        | cell_line   | genotype   |
|-----------|-----------|---------|----------------|----------------|-------------|-------------|------------|
| Otter     | old_otter | Batch1  | Otter          | Batch1         | stem_1017BM | Otter       | KO         |
| Koala_2   | old_Koala | Batch1  | Koala_2        | Batch1         | stem_1017BM | Koala_2     | WT         |
| 1017_G3-8 | otter     | Batch2  | 1017_G3-8      | Batch2         | stem_1017SM | 1017_G3-8   | KO         |
| 1017_NC-1 | Koala_2   | Batch2  | 1017_NC-1      | Batch2         | stem_1017SM | 1017_NC-1   | WT         |
| Lemur     | old_lemur | Batch1  | Lemur          | Batch1         | stem_4001BM | Lemur       | KO         |
| Whale     | old_whale | Batch1  | Whale          | Batch1         | stem_4001BM | Whale       | WT         |
| 4001_G3-1 | lemur     | Batch2  | 4001_G3-1      | Batch2         | stem_4001SM | 4001_G3-1   | KO         |
| 4001_NC-1 | whale     | Batch2  | 4001_NC-1      | Batch2         | stem_4001SM | 4001_NC-1   | WT         |
| Zebra     | old_zebra | Batch1  | Zebra          | Batch1         | stem_4002BM | Zebra       | KO         |
| Sloth     | old_sloth | Batch1  | Sloth          | Batch1         | stem_4002BM | Sloth       | WT         |
| 4002_G3-2 | zebra     | Batch2  | 4002_G3-2      | Batch2         | stem_4002SM | 4002_G3-2   | KO         |
| 4002_NC-1 | sloth     | Batch2  | 4002_NC-1      | Batch2         | stem_4002SM | 4002_NC-1   | WT         |
```
