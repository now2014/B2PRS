#  Polygenic Risk Score (PRS) Calculation



## Prerequisite

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
chunk.size <- 1e3 # test for small bed file with 5 samples


PRS.results <- PRS.wg(
  EAF.b.rds='B2PRS/PRS-b/MaleScore/EAF.b.SumHEM.FEMALE_n16600_xgb_pred_proba.regenie.rds',
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
  head(readRDS('B2PRS/PRS-b/snp.info.SumHEM.LDpred.rds'))
  #   CHROM    POS       rsid EA OA
  # 1     1 752721  rs3131972  G  A
  # 2     1 754182  rs3131969  G  A
  # 3     1 760912  rs1048488  T  C
  # 4     1 768448 rs12562034  A  G
  # 5     1 779322  rs4040617  G  A
  # 6     1 838555  rs4970383  A  C
  ```

  CHROM (**required**): Chromosome number (1, 2, 3, ..., 22, 23)

  POS (optional): GRCh37 position

  rsid  (**required**): RSID

  EA (**required**): Effect Allele

  OA (**required**): Other allele

* Columns in the `EAF.b.rds`

  ```R
  head(readRDS('B2PRS/PRS-b/MaleScore/EAF.b.SumHEM.FEMALE_n16600_xgb_pred_proba.regenie.rds'))
  #        EAF             b
  # 1 0.838825  2.734815e-10
  # 2 0.868619  1.987269e-06
  # 3 0.838086 -4.706651e-10
  # 4 0.106446  1.674726e-06
  # 5 0.129488 -2.342625e-06
  # 6 0.244398 -7.139709e-06
  ```

  EAF (**required**): Effect Allele Frequency

  b (**required**): Independent genetic effect of EA



## Pre-computed RDS files

This repository contains several pre-computed RDS files stored in the `./PRS-b` directory:

1. **TopSNPs** method (LD pruned ($r^2$ < 0.1) variants with GWAS p-value < 5e-8)

* `snp.info.rds`

| GWAS sumstats            | snp.info.rds                                                 |
| ------------------------ | ------------------------------------------------------------ |
| MaleScore/ALL_n30677     | MaleScore/snp.info.TopSNPs.ALL_n30677_xgb_pred_proba.regenie.rds |
| MaleScore/FEMALE_n16600  | MaleScore/snp.info.TopSNPs.FEMALE_n16600_xgb_pred_proba.regenie.rds |
| MaleScore/MALE_n14077    | MaleScore/snp.info.TopSNPs.MALE_n14077_xgb_pred_proba.regenie.rds |
| ProtSexGap/FEMALE_n16600 | ProtSexGap/snp.info.TopSNPs.FEMALE_n16600_ProtSexGap.rds     |
| ProtSexGap/MALE_n14077   | ProtSexGap/snp.info.TopSNPs.MALE_n14077_ProtSexGap.rds       |

* `EAF.b.rds`

| GWAS sumstats            | EAF.b.rds                                                    |
| ------------------------ | ------------------------------------------------------------ |
| MaleScore/ALL_n30677     | MaleScore/EAF.b.TopSNPs.ALL_n30677_xgb_pred_proba.regenie.rds |
| MaleScore/FEMALE_n16600  | MaleScore/EAF.b.TopSNPs.FEMALE_n16600_xgb_pred_proba.regenie.rds |
| MaleScore/MALE_n14077    | MaleScore/EAF.b.TopSNPs.MALE_n14077_xgb_pred_proba.regenie.rds |
| ProtSexGap/FEMALE_n16600 | ProtSexGap/EAF.b.TopSNPs.FEMALE_n16600_ProtSexGap.rds        |
| ProtSexGap/MALE_n14077   | ProtSexGap/EAF.b.TopSNPs.MALE_n14077_ProtSexGap.rds          |



2. **SumHEM** method

* `snp.info.rds='PRS-b/snp.info.SumHEM.LDpred.rds' `
* `EAF.b.rds`

| GWAS sumstats            | EAF.b.rds                                                    |
| ------------------------ | ------------------------------------------------------------ |
| MaleScore/ALL_n30677     | MaleScore/EAF.b.SumHEM.ALL_n30677_xgb_pred_proba.regenie.rds |
| MaleScore/FEMALE_n16600  | MaleScore/EAF.b.SumHEM.FEMALE_n16600_xgb_pred_proba.regenie.rds |
| MaleScore/MALE_n14077    | MaleScore/EAF.b.SumHEM.MALE_n14077_xgb_pred_proba.regenie.rds |
| ProtSexGap/META_n30677   | ProtSexGap/EAF.b.SumHEM.META_n30677_ProtSexGap.rds           |
| ProtSexGap/FEMALE_n16600 | ProtSexGap/EAF.b.SumHEM.FEMALE_n16600_ProtSexGap.rds         |
| ProtSexGap/MALE_n14077   | ProtSexGap/EAF.b.SumHEM.MALE_n14077_ProtSexGap.rds           |



3. **LDpred** method

* `snp.info.rds='PRS-b/snp.info.SumHEM.LDpred.rds' `
* `EAF.b.rds`

| Model       | GWAS sumstats            | EAF.b.rds                                                    |
| ----------- | ------------------------ | ------------------------------------------------------------ |
| LDpred-auto | MaleScore/ALL_n30677     | MaleScore/EAF.b.LDpred-Auto.ALL_n30677_xgb_pred_proba.regenie.rds |
| LDpred-auto | MaleScore/FEMALE_n16600  | MaleScore/EAF.b.LDpred-Auto.FEMALE_n16600_xgb_pred_proba.regenie.rds |
| LDpred-auto | MaleScore/MALE_n14077    | MaleScore/EAF.b.LDpred-Auto.MALE_n14077_xgb_pred_proba.regenie.rds |
| LDpred-inf  | MaleScore/ALL_n30677     | MaleScore/EAF.b.LDpred-Inf.ALL_n30677_xgb_pred_proba.regenie.rds |
| LDpred-inf  | MaleScore/FEMALE_n16600  | MaleScore/EAF.b.LDpred-Inf.FEMALE_n16600_xgb_pred_proba.regenie.rds |
| LDpred-inf  | MaleScore/MALE_n14077    | MaleScore/EAF.b.LDpred-Inf.MALE_n14077_xgb_pred_proba.regenie.rds |
| LDpred-auto | ProtSexGap/META_n30677   | ProtSexGap/EAF.b.LDpred-Auto.META_n30677_ProtSexGap.rds      |
| LDpred-auto | ProtSexGap/FEMALE_n16600 | ProtSexGap/EAF.b.LDpred-Auto.FEMALE_n16600_ProtSexGap.rds    |
| LDpred-auto | ProtSexGap/MALE_n14077   | ProtSexGap/EAF.b.LDpred-Auto.MALE_n14077_ProtSexGap.rds      |
| LDpred-inf  | ProtSexGap/META_n30677   | ProtSexGap/EAF.b.LDpred-Inf.META_n30677_ProtSexGap.rds       |
| LDpred-inf  | ProtSexGap/FEMALE_n16600 | ProtSexGap/EAF.b.LDpred-Inf.FEMALE_n16600_ProtSexGap.rds     |
| LDpred-inf  | ProtSexGap/MALE_n14077   | ProtSexGap/EAF.b.LDpred-Inf.MALE_n14077_ProtSexGap.rds       |

