# Parkinson's disease age-stratified polygenic risck score



### ==================== 1. PRS calculations ====================

#### PRSice-2 package

```
Rscript PRSice.R

--prsice PRSice_linux

--target /home/charlise/runs/charlise/extract/MCGILL/MCGILL_imp_soft_1800

--base PD_GWAS_2019.PRS.1800.txt

--beta --stat b

--snp SNP

--A1 A1

--A2 A2

--score avg

--binary-target T

--prevalence 0.005

--pvalue p

--out MCGILL_imp_soft_1800_PRS

--print-snp

--fastscore

--no-regress
```

### ==================== 2. Logistic regression in R ====================
```R

library("data.table")

library("dplyr")

#strat by age

study <- fread('MCGILL_imp_soft_1800_Z.all_score')

#convert age and sex to 0/1

study$sex <- study$sex - 1

study$phenotype <- study$phenotype - 1

#stratify by age

age_50_below <- subset(study, age < 50)

age_70_above <- subset(study, age > 70)

age50.60 <- subset(study, study$age <= 60 & study$age >= 50)

age60.70 <- subset(study, study$age <= 70 & study$age > 60)

#function for sum stat

fit_sum <- function(k)

{

fit <- glm(phenotype ~ Z + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = k, family = binomial())

summary(fit)

}

fit_sum(age_50_below)

fit_sum(age50.60)

fit_sum(age60.70)

fit_sum(age_70_above)
```

### ==== 3. Meta-analysis =====

#### METAL

```
MARKER age_grp

PVALUE p

EFFECT beta

WEIGHT ntotal

SCHEME STDERR

STDERR error

PROCESS  IPDGC.txt

PROCESS  APDGC.txt

PROCESS  NIND.txt

PROCESS  NGRC.txt

PROCESS  IPDGC.txt

PROCESS  mcgill.txt

PROCESS  Oslo.txt

OUTFILE  METAL-stratified_PRS .tbl

ANALYZE HETEROGENEITY

QUIT
```
