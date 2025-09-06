# wrapped up all inputting paramters into one function
dkosim <- function(sample_name, 
                   coverage, n, n_guide_g, n_gene_pairs, n_construct, library_size, sd_freq0,
                   moi, moi_pois, p_gi, sd_gi, p_high, mode,
                   pt_neg, pt_pos, pt_wt, pt_ctrl,
                   mu_neg, sd_neg, mu_pos, sd_pos, sd_wt,
                   bottleneck, n.bottlenecks, n.iterations, resampling){

# print out initialized parameters for this run
cat("
# ------------------------------------------------------------
# Simulation Settings Summary:
# ------------------------------------------------------------
# Sample Name:", sample_name, "
# Number of Genes:", n, "
# Cell Library Size (Initial):", library_size, "
# Coverage:", coverage,"x
# Number of Single Knockout(SKO):", n, "
# Number of Double Knockout(DKO):", n_gene_pairs - n, "
# Number of Guides per Gene:", n_guide_g, "
# Number of Constructs:", n_construct, "
# Variance of Initialized Counts:", round(sd_freq0^2,2), "

# Genetic Interactions (GI):
## Proportion of GI(%):", p_gi * 100, "
## Number of Interacting Gene Pairs:", round(p_gi * (n_gene_pairs-n)), " 
## Variance of re-sampled phenotypes w/ GI:", sd_gi^2, "

# Proportion of Each Initialized Gene Class (by theoretical phenotypes):
## Negative(%):", pt_neg * 100, "~ TN(", mu_neg, "," , sd_neg^2, ",-1,-0.025)
## Positive(%):", pt_pos * 100, "~ TN(", mu_pos, "," ,sd_pos^2, ",0.025, 1)
## Wild-Type(%):", pt_wt * 100, "~ TN(0,", sd_wt^2, ",-0.025, 0.025)
## Non-Targeting Control(%):", pt_ctrl * 100, "~ Delta(0)               

# Proportion of Guides (by efficacy):
## High-efficacy(%):", p_high * 100, "~ 1
## Low-efficacy(%):", (1-p_high) * 100, "~ TN(0.05, 0.0049, 0, 1)

# Multiplicity of Infection (MOI):", moi , "
# Percentage of viral particles delivered in cells during transfection(%):", round(dpois(1, lambda = moi) * 100, 2) , "~ Poisson(",moi,")
# Resampling Size based on MOI (Passage Size):", resampling, "
# Bottleneck Size (", bottleneck / library_size, "x Initial Guide-Level Library Size):", bottleneck, "
# Number of Bottleneck Encounters (Number of Passages):", n.bottlenecks, "
# Total Available Doublings:", n.iterations, "
# Number of Replicates:", 2, "
# Pseudo-count:", 5 * 10^(-floor(log10(bottleneck))-1), "

# ------------------------------------------------------------
\n")


  # %% [markdown]
  # ## PART 1: Parameters Specifications
  
  # %%
  ################### STEP 1: define functions ######################################
  # FUNCTION 1: define p_y calculation function
  p_y<-function(p1,p2){
    fp<-function(x){max(c(0,x))}
    prob2=c(fp(-p2),1-abs(p2),fp(p2))
    prob1=c(fp(-p1),1-abs(p1),fp(p1))  
    a=outer(prob1,prob2)
    # define the p^y vectors
    return(c(a[1,1]+a[1,2]+a[2,1]+a[1,3]+a[3,1],
             a[2,2],
             a[3,2]+a[2,3],
             a[3,3])
    )
  }
  
  # FUNCTION 2: define E(Cy|p_y) - update from E(y|.) to E(C_y|.)
  EC_y<-function(p0_0,p1_0,p2_0,p3_0, p0,p1,p2,p3){
    return(0*(p0 - p0_0) + 2*(p1 - p1_0) + 4*(p2 - p2_0) + 8*(p3 - p3_0))
  }
  
  # FUNCTION 3: simulate p1 by pre-specify gene class
  p1_sample <- function(n, pt_neg, pt_pos, pt_wt, pt_ctrl) {
    # sample counts per class
    n_neg <- round(pt_neg * n)
    n_pos <- round(pt_pos * n)
    n_wt  <- round(pt_wt * n)
    n_ctrl <- round(pt_ctrl * n)
    
    # sample values
    p_neg <- rtruncnorm(n_neg, a = -1, b = -0.025, mean = mu_neg, sd = sd_neg)
    p_pos <- rtruncnorm(n_pos, a = 0.025, b = 1, mean = mu_pos, sd = sd_pos)
    p_wt  <- rtruncnorm(n_wt, a = -0.025, b = 0.025, mean = 0, sd = sd_wt)
    p_ctrl <- rep(0, n_ctrl)
    
    # combine into one vector
    p1 <- c(p_neg, p_pos, p_wt, p_ctrl)
    
    # create class labels
    class <- c(rep("negative", n_neg),
               rep("positive", n_pos),
               rep("WT", n_wt),
               rep("Non-Targeting Control", n_ctrl))
    
    # return data frame
    return(data.frame(p1 = p1, class = class))
  }
  
  # FUNCTION 4: grow the cell from sampling by relative frequency
  ## need to decide when the cell does not divide, its' yielding 0 or 1 cells? ko essential gene lead to cell death - use 0 instead
  cell_grow <- function(n0, a){
    n1 = 0
    if (n0 == 0) {
      return(0)
    }
    else{
      multi_dist = rmultinom(n0, size = 1, prob = a)
      for (i in 1:n0){
        n1[i] = multi_dist[,i][1]*0 + multi_dist[,i][2]*2 + multi_dist[,i][3]*4 + multi_dist[,i][4]*8
      }
      return(sum(n1))
    } 
  }
  
  # FUNCTION 5: guides per gene generation
  guide_sample <- function(n, mode, high_prop) {
    # print error messages if input is not correct
    if (!(mode %in% c("CRISPRn", "CRISPRn-100%Eff"))) stop("mode must be 'CRISPRn' or 'CRISPRn-100%Eff'")
    # determine the gudie type
    guide_type <- sample(c(TRUE, FALSE), size = n, replace = TRUE, prob = c(high_prop, 1 - high_prop))
    # specify different distribution based on CRISPRi or CRISPRn
    if (mode == "CRISPRn") {
      guide_eff <- ifelse(
        guide_type,
        rtruncnorm(n, a = 0.6, b = 1, mean = 0.9, sd = 0.1),
        rtruncnorm(n, a = 0, b = 0.6, mean = 0.05, sd = 0.07)
      )
    } else if (mode == "CRISPRn-100%Eff") {
      guide_eff <- ifelse(
        guide_type,
        1,
        rtruncnorm(n, a = 0, b = 1, mean = 0.05, sd = 0.07)
      )
    }
    
    return(cbind(guide_type, guide_eff))
  }
  
  # FUNCTION 6: genetic interaction values generation - increase sd to enlarge the range of GI
  gi_sample <- function(p){
    p_i = rtruncnorm(1, a = -1, b = 1, mean = p, sd = sd_gi)
    return(p_i)
  }
  
  # %% [markdown]
  # ## PART 2: Gene-Level Cell Library Initialization
  library(dplyr)
  # Optimization: Wrapped up into function
  initialize_gene_cell_lib0 <- function(){
    ################### STEP 2: initialize p1 & p2 ######################################
    # gene1: simulate first n unique gene theoretical phenotype using defined function
    gene1_p1 = cbind(gene1 = 1:n, p1_sample(n, pt_neg, pt_pos, pt_wt, pt_ctrl)) %>% 
      as.data.frame() %>% 
      dplyr::mutate(
        # override or keep the original class from p1_sample
        gene1_behavior = class,  # use this to keep true class labels
        freq_t0_g1 = 10^(rnorm(n, mean = 0, sd = sd_freq0)) # using 80% for the expected z-score to freq_t0 - borrow from CRISPulator, z = 2 * qnorm(0.975)
      ) %>%
      dplyr::select(-class)
    # p2: simulate same phenotype as p1 to make sure the total number of unique genes is n
    gene2_p2 = gene1_p1 %>% 
      dplyr::rename(gene2 = gene1,
                    p2 = p1,
                    gene2_behavior = gene1_behavior,
                    freq_t0_g2 = freq_t0_g1)
    
    ############## STEP 3: initialize the interaction index between pre-specified genes##############
    ############## Generate all unique combinations of two elements from the sequence########
    combinations <- t(combn(1:(n+1), 2))
    ## Convert to a data frame and add index
    combinations_df <- as.data.frame(combinations)
    combinations_df$gene_pair_id = 1:nrow(combinations_df)
    combinations_df = combinations_df %>% 
      mutate(V2 = ifelse(V2-1 == V1, NA, V2-1)) %>% 
      dplyr::rename(gene1_id = V1,
                    gene2_id = V2)
    
    ########### STEP 4: initialize gene-level cell library ###############################
    # initialize cell library columns
    cell_lib0<-bind_cols(gene1=combinations_df$gene1_id,
                         gene2=combinations_df$gene2_id,
                         gene_pair_id=combinations_df$gene_pair_id,
                         p1i = NA,
                         p2i = NA,
                         p0y_0 = NA,
                         p1y_0 = NA,
                         p2y_0 = NA,
                         p3y_0 = NA,
                         p0y = NA,
                         p1y = NA,
                         p2y = NA,
                         p3y = NA,
                         gene1_behavior = NA,
                         gene2_behavior = NA,
                         #KL_divergence = NA,
                         interaction_gene = NA,
                         interaction_gene_type = NA) %>% 
      inner_join(gene1_p1, by = "gene1") %>% 
      left_join(gene2_p2, by = "gene2") %>% 
      # counts are the counts at later timept, set it equal to initial timept at the first stage
      mutate(gene2 = ifelse(is.na(gene2), NA, gene2),
             p2 = ifelse(is.na(p2), 0, p2),
             gene1_behavior = gene1_behavior.y,
             gene2_behavior = gene2_behavior.y,
             KO_type = ifelse(is.na(gene2), "SKO", "DKO"), 
             freq_t0 = ifelse(is.na(gene2), freq_t0_g1,
                              (freq_t0_g1 + freq_t0_g2)/2),
             rel_freq_t0 = freq_t0/sum(freq_t0),
             counts_t0 = round(rel_freq_t0 * library_size),
             counts_t1 = counts_t0
      )
    # sampling for genetic interaction index from pre-specified GI proportion w/o replacement
    dko_indices = which(!is.na(cell_lib0$gene2))
    gi_sampled_indices = sample(which(!is.na(cell_lib0$gene2) &
                                        cell_lib0$gene1_behavior != "Non-Targeting Control" &
                                        cell_lib0$gene2_behavior != "Non-Targeting Control"), 
                                size = round(length(dko_indices) * p_gi), replace = FALSE)
    # update initialized cell library with sampled genetic interaction index
    cell_lib0 <- cell_lib0 %>%
      mutate(i_interaction = as.integer(seq_len(nrow(cell_lib0)) %in% gi_sampled_indices),
             p1i = if_else(i_interaction == 0, p1, gi_sample(p1)),
             p2i = if_else(i_interaction == 0, p2, gi_sample(p2))) %>%
      dplyr::select(
        gene_pair_id, gene1, gene2, p1, p2, p1i, p2i, KO_type,
        gene1_behavior, gene2_behavior, i_interaction,
        interaction_gene, interaction_gene_type,
        p0y_0, p1y_0, p2y_0, p3y_0, p0y, p1y, p2y, p3y,
        freq_t0, rel_freq_t0, counts_t0, counts_t1)
    return(cell_lib0)
  }
  
  # %% [markdown]
  # ## PART 3: Guide-Level Cell Library Initialization
  # Optimization: Wrap up into function
  initialize_guide_cell_lib0 <- function(rep_name, cell_lib0){
    library(dplyr)
    library(purrr)
    library(gtools)
    set.seed(ifelse(rep_name == "repA", 666, 777))  # Different seed for Rep A & Rep B
    ########### STEP 5: initialize guide1 & guide2 w/ type and efficiency ################################################################
    # initialize guides
    guide1 = cbind(1:(n*n_guide_g), guide_sample(n*n_guide_g, mode, p_high)) %>% 
      as.data.frame() %>% 
      dplyr::rename(guide1 = V1,
                    guide1_type = guide_type,
                    guide1_eff = guide_eff) %>% 
      mutate(gene1=rep(c(1:n), n_guide_g)) %>% 
      arrange(gene1) %>% 
      mutate(guide1_id = paste0(gene1, "_", rep(c(1:n_guide_g), n)))
    
    guide2 = guide1 %>% 
      dplyr::rename(guide2 = guide1,
                    guide2_type = guide1_type,
                    guide2_eff = guide1_eff,
                    gene2 = gene1,
                    guide2_id = guide1_id)
    
    ########### STEP 6: initialize guide-level cell library based on initialized gene-level cell library ###############################
    # initialize guide-level data by expanding the gene-level data
    cell_lib_guide0 <- cell_lib0 %>%
      dplyr::mutate(
        guides_number = if_else(KO_type == "SKO", n_guide_g, n_guide_g^2),
        clone_idx = row_number()
      ) %>%
      slice(rep(1:n(), guides_number)) %>%
      group_by(clone_idx) %>%
      dplyr::mutate(guide = row_number()) %>%
      ungroup() %>%
      group_split(clone_idx) %>%
      map_dfr(function(df) { # assign the initialized guide-level freq by random draws from Dirichlet with large variance
        g_n <- nrow(df)
        # use large variance
        alpha <- rep(100, g_n)
        # Sample Dirichlet weights
        dirich_weights <- as.vector(rdirichlet(1, alpha))
        freq_raw <- df$freq_t0[1] * dirich_weights
        df$freq_dirich <- dirich_weights
        df$freq_guide_t0 <- freq_raw
        return(df)
      }) %>%
      dplyr::mutate(
        freq_guide_t0 = ifelse(KO_type == "DKO", freq_guide_t0 * n_guide_g, freq_guide_t0), # multiply initialized guide-level freq from Dirichlet by number of guides if it's a DKO
        rel_freq_guide_t0 = freq_guide_t0 / sum(freq_guide_t0),
        counts_guide_t0 = round(rel_freq_guide_t0 * library_size),
        counts_guide_t1 = counts_guide_t0,
        guide1_id = paste0(gene1, "_", rep(1:n_guide_g, times = nrow(.) / n_guide_g)),
        guide2_id = ifelse(is.na(gene2), NA,
                           paste0(gene2, "_", rep(1:n_guide_g, each = n_guide_g, times = nrow(.) / (n_guide_g^2)))),
        construct_id = paste0(guide1_id, "__", ifelse(is.na(gene2), "0", guide2_id))
      ) %>%
      inner_join(guide1, by = "guide1_id") %>%
      left_join(guide2, by = "guide2_id") %>%
      dplyr::mutate(
        gene1 = gene1.x,
        gene2 = gene2.x,
        p1ip = guide1_eff * p1i,
        p2ip = ifelse(is.na(gene2), 0, guide2_eff * p2i),
        p1p = guide1_eff * p1,
        p2p = ifelse(is.na(gene2), 0, guide2_eff * p2),
        gene1_gene2_id = paste(gene1, ifelse(is.na(gene2), "0", gene2), sep = "_"),
        p0y_p0 = NA, p1y_p0 = NA, p2y_p0 = NA, p3y_p0 = NA,
        p0y_p = NA, p1y_p = NA, p2y_p = NA, p3y_p = NA,
        interaction_guide = NA,
        interaction_guide_type = NA,
        guide_id = row_number()
      ) %>%
      dplyr::select(
        guide_id, gene_pair_id, gene1_gene2_id,
        guide1_id, guide2_id, construct_id,
        gene1, gene2, guide1_eff, guide2_eff,
        p1, p2, p1i, p2i, p1p, p2p, p1ip, p2ip,
        KO_type, gene1_behavior, gene2_behavior,
        i_interaction, interaction_gene, interaction_guide,
        interaction_gene_type, interaction_guide_type,
        p0y_0, p1y_0, p2y_0, p3y_0,
        p0y, p1y, p2y, p3y,
        p0y_p0, p1y_p0, p2y_p0, p3y_p0,
        p0y_p, p1y_p, p2y_p, p3y_p,
        freq_t0, rel_freq_t0,
        guide1_type, guide2_type, freq_guide_t0,
        counts_guide_t0, counts_guide_t1, rel_freq_guide_t0
      )
    
    return(cell_lib_guide0)
  }
  
  # %%
  ########### STEP 7: calculate the theoretical phenotype with/without genetic interaction##############################
  # Optimization: Wrap up into function
  define_phenotype_gi <- function(cell_lib_guide0){
    # Parallel guide-level data generation: replace gene-level initialization cell_lib0 to guide-level initialization cell_lib_guide0
    n <- nrow(cell_lib_guide0)  # Number of rows
    results0 <- foreach(i = 1:n, .combine = rbind, .packages = c("gtools", "dplyr")) %dopar% {
      # Initialize the output row for this iteration
      row_result <- cell_lib_guide0[i, ]
      
      # Calculate gene-level theoretical phenotype w/o genetic interaction
      p_vals_0 <- p_y(cell_lib_guide0$p1[i], cell_lib_guide0$p2[i])
      row_result$p0y_0 <- p_vals_0[1]
      row_result$p1y_0 <- p_vals_0[2]
      row_result$p2y_0 <- p_vals_0[3]
      row_result$p3y_0 <- p_vals_0[4]
      # Calculate gene-level theoretical phenotype with genetic interaction
      p_vals_1 <- p_y(cell_lib_guide0$p1i[i], cell_lib_guide0$p2i[i])
      row_result$p0y <- p_vals_1[1]
      row_result$p1y <- p_vals_1[2]
      row_result$p2y <- p_vals_1[3]
      row_result$p3y <- p_vals_1[4]
      # Calculate guide-level theoretical phenotype with guide efficacy but w/o genetic interaction
      p_vals_2 <- p_y(cell_lib_guide0$p1p[i], cell_lib_guide0$p2p[i])
      row_result$p0y_p0 <- p_vals_2[1]
      row_result$p1y_p0 <- p_vals_2[2]
      row_result$p2y_p0 <- p_vals_2[3]
      row_result$p3y_p0 <- p_vals_2[4]  
      # Calculate guide-level theoretical phenotype with both guide efficacy and genetic interaction
      p_vals_3 <- p_y(cell_lib_guide0$p1ip[i], cell_lib_guide0$p2ip[i])
      row_result$p0y_p <- p_vals_3[1]
      row_result$p1y_p <- p_vals_3[2]
      row_result$p2y_p <- p_vals_3[3]
      row_result$p3y_p <- p_vals_3[4]
      
      # Calculate genetic interaction term
      ## gene-level
      interaction_value1 <- EC_y(row_result$p0y_0, row_result$p1y_0, row_result$p2y_0, row_result$p3y_0, 
                                 row_result$p0y, row_result$p1y, row_result$p2y, row_result$p3y)
      ## guide-level
      interaction_value2 <- EC_y(row_result$p0y_p0, row_result$p1y_p0, row_result$p2y_p0, row_result$p3y_p0, 
                                 row_result$p0y_p, row_result$p1y_p, row_result$p2y_p, row_result$p3y_p)
      
      # update both gene-level interaction and guide-level interaction values
      row_result$interaction_gene <- interaction_value1
      row_result$interaction_guide <- interaction_value2
      # add categories based on interaction values
      row_result$interaction_gene_type <- ifelse(interaction_value1 < 0, "negative",
                                                 ifelse(interaction_value1 > 0, "positive", "none"))
      row_result$interaction_guide_type <- ifelse(interaction_value2 < 0, "negative",
                                                  ifelse(interaction_value2 > 0, "positive", "none"))
      
      # Return the modified row as the result of this iteration
      row_result
    }
    return(results0)
  }
  
  
  # %% [markdown]
  # ## PART 4: Build Guide-level Simulation Schemes
  # %%
  ########### STEP 8: Sampling by relative frequency of each construct and update the storing output ################
  # Define a function to run Replicates in parallel
  run_replicate <- function(replicate_name, cell_lib_guide0) {
    # create and write on a log file to track the iterations and bottleneck encountering
    log_file <- paste0("./logs/", sample_name, "_", replicate_name, "_log.txt")
    if (!dir.exists("logs")) {
      dir.create("logs")
    }
    write(paste0(Sys.time(), " - ", replicate_name, " started execution\n"), 
          file = log_file, append = TRUE)
    
    i.iterations <- 0
    i.bottleneck <- 0
    cell_lib_guide1 <- copy(cell_lib_guide0)
    
    while (i.bottleneck < n.bottlenecks & i.iterations < n.iterations) {
      param_matrix <- as.matrix(cell_lib_guide1[, c("p0y_p", "p1y_p", "p2y_p", "p3y_p")])     # define para_matrix for vectorized cell growth
      # Check for bottleneck condition
      if (sum(cell_lib_guide1$counts_guide_t1) > bottleneck) {
        # update bottleneck counters
        i.bottleneck <- i.bottleneck + 1
        write(paste0(Sys.time(), " - Bottleneck encountered: ", i.bottleneck, "\n"), file = log_file, append = TRUE)
        # use hypergeometric distribution to draw w/o replacement
        cell_lib_guide1$counts_guide_t1 = as.numeric(rmvhyper(nn = 1, k = resampling, n = cell_lib_guide1$counts_guide_t1))
      }
      
      #### cell population growth
      cell_lib_guide1$counts_guide_t1 <- vapply(seq_len(nrow(cell_lib_guide1)), function(i) {
        cell_grow(cell_lib_guide1$counts_guide_t1[i], param_matrix[i, ])
      }, numeric(1))
      
      # Update iteration counter to log file
      i.iterations <- i.iterations + 1
      write(paste0(Sys.time(), "-", replicate_name, "- Iteration:", i.iterations),
            file = log_file, append = TRUE)
      
    }
    
    # update the log file for completing execution
    write(paste0(Sys.time(), "-", replicate_name, "- completed execution"), 
          file = log_file, append = TRUE)
    
    # Normalization, Log Fold Change(LFC) Calculation and Save outputs
    ## adjust pseudo_counts according to the bottleneck number
    pseudo_counts = 5 * 10^(-floor(log10(bottleneck))-1)
    if (!dir.exists("data")) {
      dir.create("data")
    }
    ## stored updated cell library and calculate relative frequency
    cell_lib_guide2 = cell_lib_guide1 %>%
      mutate(counts_guide_t2 = counts_guide_t1,
             rel_freq_guide_t2 = counts_guide_t2 / sum(counts_guide_t2), 
             LFC = log2(((rel_freq_guide_t2 + pseudo_counts) / (rel_freq_guide_t0 + pseudo_counts)))) # add pseudocounts to calculate log fold change to avoid infinity
    write.csv(cell_lib_guide2, paste0("./data/", sample_name, "_", replicate_name, ".csv"))
  }
  
  # ## PART 5: Run Simulations
  ########### STEP 9: Utilize all written functions to actually run simulations ##############################
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
  
}
