#  Polygenic Risk Score (PRS) Calculation



## Prerequisite

* PLINK

* R

* R packages:

  ```R
  # Install packages in R
  
  install.packages('BEDMatrix')
  install.packages('data.table')
  install.packages('parallel')
  ```

  


## Test

* Clone the repository

```bash
git clone https://github.com/now2014/B2PRS.git
```

* Test in R

Here, the genotype data in the `B2PRS/test-bed/` directory, which comes from the 1000 Genomes Project.

```R
bed.file.pattern <- 'B2PRS/test-bed/1kg-eur-chr@.bed' # @ is replaced by chromosome number
source('B2PRS/B2PRS.r')

# max number of SNPs to be read from bed file in each chunk
# for bed file with large number of samples, reduce chunk.size to avoid memory overflow, e.g., 1e3
chunk.size <- 1e5 # test for small bed file with 5 samples


PRS.results <- PRS.wg(
  EAF.b.rds='B2PRS/PRS-b/EAF.b.SumHEM.FEMALE_n16600_xgb_pred_proba.regenie.rds',
  snp.info.rds='B2PRS/PRS-b/snp.info.SumHEM.LDpred.rds',
  bed.file.pattern=bed.file.pattern, fill.missing=c('mean'), chunk.size=chunk.size,
  mc.cores=1, PRS.name='PRS.SumHEM.FEMALE_n16600'
)

print(PRS.results)
###### Output for the test #######
#         PRS.SumHEM.FEMALE_n16600
# HG00096               0.02322992
# HG00097               0.01947996
# HG00099               0.02164128
# HG00101               0.02291205
# HG00102               0.02969029
```



## Parameters of the `PRS.wg` function

| Param            | Default                | Meaning                                                      |
| ---------------- | ---------------------- | ------------------------------------------------------------ |
| EAF.b.rds        | './PRS-b/EAF.b.rds'    | path to the RDS file with `'EAF'` and `'b'` for PRS calculation |
| snp.info.rds     | './PRS-b/snp.info.rds' | path to the RDS file with SNP information                    |
| bed.file.pattern | './chr@.bed'           | a pattern of bed file names, where '@' is replaced by chromosome number |
| fill.missing     | 'mean'                 | a string of `'EAF'` or `'mean'` to fill missing genotype with `'EAF'` column or `'mean'` frequency derived from the input genotype data |
| chunk.size       | 1000                   | max number of SNPs to be read from bed file in each chunk    |
| mc.cores         | 1                      | number of cores for parallel computing                       |
| PRS.name         | 'PRS'                  | name of the PRS column in the output                         |



## Input RDS files

* Please **make sure** that each **row** of the data.frame in the `snp.info.rds` **matches** the corresponding **row** of the data.frame in the `EAF.b.rds`.

* Columns in `snp.info.rds`

  ```R
  head(readRDS('B2PRS/PRS-b/PRS-b/snp.info.SumHEM.LDpred.rds'))
  #   CHROM    POS       rsid EA OA
  # 1     1 752721  rs3131972  G  A
  # 2     1 754182  rs3131969  G  A
  # 3     1 760912  rs1048488  T  C
  # 4     1 768448 rs12562034  A  G
  # 5     1 779322  rs4040617  G  A
  # 6     1 838555  rs4970383  A  C
  ```

  **CHROM: Chromosome number (1, 2, 3, ..., 22, 23)**

  **POS: GRCh37 position**

  **EA: Effect Allele**

  **OA: Other allele**

* Columns in the `EAF.b.rds`

  ```R
  head(readRDS('B2PRS/PRS-b/PRS-b/EAF.b.SumHEM.FEMALE_n16600_xgb_pred_proba.regenie.rds'))
  #        EAF             b
  # 1 0.838825  2.734815e-10
  # 2 0.868619  1.987269e-06
  # 3 0.838086 -4.706651e-10
  # 4 0.106446  1.674726e-06
  # 5 0.129488 -2.342625e-06
  # 6 0.244398 -7.139709e-06
  ```

  **EAF: Effect Allele Frequency**

  **b: Independent genetic effect of EA**

