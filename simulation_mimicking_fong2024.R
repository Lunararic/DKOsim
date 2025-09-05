# Title: Double-CRISPR Knockout Simulation (DKOsim) - Mimicking Fong-2024 A549 Data
# Author: Yue Gu, Luis Leon Novelo, John Paul Shen, Traver Hart
# Fong-2024 Data Reference: Fong, S.H., Kuenzi, B.M., Mattson, N.M. et al. A multilineage screen identifies actionable synthetic lethal interactions in human cancers. Nat Genet 57, 154–164 (2025). https://doi.org/10.1038/s41588-024-01971-9

# Load necessary libraries
if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("xtable", quietly = TRUE)) install.packages("xtable")
if (!requireNamespace("MCMCpack", quietly = TRUE)) install.packages("MCMCpack")
if (!requireNamespace("entropy", quietly = TRUE)) install.packages("entropy")
if (!requireNamespace("truncnorm", quietly = TRUE)) install.packages("truncnorm")
if (!requireNamespace("ggtern", quietly = TRUE)) install.packages("ggtern")
if (!requireNamespace("gtools", quietly = TRUE)) install.packages("gtools")
if (!requireNamespace("plyr", quietly = TRUE)) install.packages("plyr")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("dplyr")

library(xtable)
library(MCMCpack)
library(entropy)
library(truncnorm)
library(ggtern)
library(gtools)
library(plyr)
library(doParallel)
library(foreach)
library(extraDistr)
library(tidyverse)
library(data.table)
source(dkosim.R)
### Major Modifications Before Each Run
### 1.sample_name
### 2.parameter specifications: simulation settings
### 3.codes: (1)proportion of each gene class; (2)gene class distribution w/ means and std. dev

## initialize parameters
sample_name = "DKO_simulation_hpc_v8_guide_250x3x3x2_1000x_run139"
set.seed(888) # tune the random seeds to ensure the initialized guide-level counts are exceeding the same initialized bottleneck size if needed

## initialized library parameters
coverage = 1000 # coverage: cell representation per guide
n = 246 # number of unique single gene
n_guide_g = 3 # number of guide per gene
n_gene_pairs = n * (n-1) / 2 + n  # number of unique gene pairs (both SKO and DKO)
n_construct = (n*n_guide_g) * ((n-1)*n_guide_g) / 2 + n*n_guide_g  # total number of constructs
library_size = n_construct * coverage # number of total cells in the initialized gene-level library
sd_freq0 = 1/(2*qnorm(0.90)) # std dev of initialized counts distribution

## Run Simulations
# Setup parallel computing background
if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
if (!requireNamespace("gtools", quietly = TRUE)) install.packages("gtools")  # For rdirichlet()
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("dplyr")
## Load libraries
library(foreach)
library(doParallel)
library(gtools)  # For rdirichlet()
library(dplyr)
library(data.table)

## Use all available cores
num_cores <- parallel::detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
start_time <- proc.time()

# Gene-level Cell Library
##initialize shared gene-level cell library, same across replicates
cell_lib0 <- initialize_gene_cell_lib0()

# Guide-level Cell Library
## initialize independent guide-level cell library by replicates
cell_lib_guide00 <- foreach(replicate = c("repA", "repB"), .packages = c("dplyr", "data.table", "truncnorm"), .combine = "list") %dopar% {
  initialize_guide_cell_lib0(replicate, cell_lib0)
}
## extract from the list into separated initialized guide-level cell library
cell_lib_guide0_A <- cell_lib_guide00[[1]]
cell_lib_guide0_B <- cell_lib_guide00[[2]]

# Define phenotypes and GI
cell_lib_guide0 <- foreach(replicate = list(cell_lib_guide0_A, cell_lib_guide0_B), .packages = c("gtools", "dplyr", "foreach"), .combine = "list") %dopar% {
  define_phenotype_gi(replicate)
}
## extract initialized guide-level cell library with all elements
cell_lib_guide0_A <- cell_lib_guide0[[1]]
cell_lib_guide0_B <- cell_lib_guide0[[2]]
## check the initialized guide-level data
print(cell_lib_guide0_A)
print(cell_lib_guide0_B)

# Run Simulations
## Optimization: Run guide-level simulation for both replicates in parallel
results <- foreach(replicate = c("repA", "repB"), .packages = c("dplyr", "foreach", "doParallel", "tidyr", "data.table", "extraDistr"), .combine = "c") %dopar% {
  if (replicate == "repA") {
    run_replicate("repA", cell_lib_guide0_A)
  } else {
    run_replicate("repB", cell_lib_guide0_B)
  }
  return(paste(replicate, "completed"))
}
## Print the collected results
print(results)

# Stop the cluster and collect running time in seconds
stopCluster(cl)
end_time <- proc.time()
elapsed_time <- end_time - start_time; elapsed_time
# check used cores and running time
print(paste("number of cores", num_cores))
print(paste("Run Time (hrs): ", elapsed_time["elapsed"]/3600))

moi = 0.3 # moi
moi_pois = dpois(1, moi) # get the number of viral particles delivered per cell during transfection from Poisson(moi) to calculate resampling size
p_gi = 0.03 # % of genetic interaction presence
sd_gi = 1.5 # std dev of re-sampled phenoytpes w/ gi presence
## guide parameters
p_high = 0.75 # % of high-efficacy guides
mode = "CRISPRn" # CRISPR mode: use CRISPRn-100%Eff for need 100% efficient guides; use CRISPRn for high-efficient guides draw from distribution

## gene class parameters: 
### % of theoretical phenotype to each gene class - add up to 1
pt_neg = 64/246 # % negative phenotype (essential)
pt_pos = 0 # % positive phenotype (non-essential)
pt_wt = 178/246 # % wt
pt_ctrl = 4/246 # % non-targeting control
### mean and std dev of theoretical phenotypes
mu_neg = -0.03 # mean: negative phenotype
sd_neg = 0.25 # std dev: negative phenotype

mu_pos = 0.75 # mean: positive phenotype
sd_pos = 0.1 # std dev: positive phenotype

sd_wt = 0.2 # std dev: wildtype phenotype

## bottleneck parameters
bottleneck = 2 * library_size # bottleneck size
n.bottlenecks = 1 # how many times do we encountering bottlenecks?
n.iterations = 30 # assuming a maximum of 30 doubling cycles if we didnt encounter bottleneck
resampling = round(moi_pois * bottleneck)# determine resampling size based on moi and bottleneck size


