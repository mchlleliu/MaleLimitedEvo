###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                        Permutation & Bootstrap functions
# 
# 
###################################

# packages required
#########
library(tidyr)
library(plyr)
library(dplyr)
library(broom)
########

# \\||// #
# Permutation tests
##########

# OnePerm
# One sample permutation test
# args: perm_dat = the dataset containing a numeric column to permute
#       x_col = string denoting the column containing the data to permute
#       myfun = test statistic function
#       n_perm = number of permutations. default = 10,000
#       alternative = alt. for hypothesis testing
OnePerm <- function(perm_dat, x_col, 
                      myfun = mean, mu = 0, n_perm = 10000,
                      alternative = c("two.sided","less","greater")){
  set.seed(1000)
  x = perm_dat[[x_col]] # get data to permute
  x = x - mu
  n = length(x) # number of data points
  obs.stat = myfun(x) # get observed test statistic
  
  # make 10,000 matrices of g vectors denoting observations that are above (g=1) or below (g=-1) the null
  gmat = replicate(n_perm, sample(x = c(1,-1),
                                 size = n,
                                 replace = TRUE))
  
  # randomly flip/keep signs for each observation using gmat, and
  # find permuted test stat for each permutation
  perm.stat = apply(gmat*abs(x), 2, myfun) 
  
  # assign p-value
  if(alternative[1]=="less"){
    pval = (sum(perm.stat <= obs.stat) + 1) / (n_perm + 1)
  } else if(alternative[1]=="greater"){
    pval = (sum(perm.stat >= obs.stat) + 1) / (n_perm + 1)
  } else{
    pval = (sum(abs(perm.stat) > abs(obs.stat)) + 1) /(n_perm + 1)
  }
  
  # store results in data frame object and return
  summary.test <- c(obs.stat, n, pval, pval < 0.05)
  return(summary.test)
}


# Function to permute based on SBGE category
# args: perm_dat = the dataset containing a numeric column to permute
#       x_col =string denoting the column containing the data to permute
#       SBGE_cat = string denoting the column which specifies SBGE bins
#       n_perm = optional argument for number of permutations. default is 10,000
OnePerm_SBGE <- function(perm_dat, x_col, SBGE_cat, n_perm = 10000){
  dat <- data.frame() # initialize data.frame object to store results
  # make sure the column of categories is set to factor
  perm_dat[[SBGE_cat]] <- factor(perm_dat[[SBGE_cat]])
  
  # go through each level (SBGE categories)
  for(i in levels(perm_dat[[SBGE_cat]])){
    # for each SBGE category, call the permutation function
    perm_cat <- c(i, OnePerm(perm_dat[perm_dat[[SBGE_cat]] == i,], 
                               x_col = x_col, alternative = "two.sided"))
    
    # concatenate results in the data.frame to return
    dat <- rbind(dat, perm_cat)
  }
  
  colnames(dat) <- c(SBGE_cat, "obs.stat", "n", "pval", "Sig")
  dat$obs.stat <- as.numeric(dat$obs.stat)
  dat$n <- as.numeric(dat$n)
  return(dat)
}



# TwoPerm
# Two sample permutation test for difference in means
# args: perm_dat = the dataset containing a numeric column to permute
#       x_col = string denoting the column containing the data to permute
#       groupBy = logical column in the dataset which separates group 1 and group 2
#       n_perm = number of permutations. default = 10,000
#       alternative = alt. for hypothesis testing
TwoPerm <- function(perm_dat, x_col, 
                      groupBy, n_perm = 10000,
                      alternative = c("two.sided","less","greater")){
  set.seed(1000)
  # separate data for groups 1 and 2
  trt1 <- perm_dat[perm_dat[[groupBy]],][[x_col]]
  trt2 <- perm_dat[!perm_dat[[groupBy]],][[x_col]]
  # number of observations
  n1 <- length(trt1)
  n2 <- length(trt2)
  
  # observed difference in mean
  obs.diff <- mean(trt1) - mean(trt2)  
  # shuffle data from trt1 and trt2
  rand.samples <- sample(c(trt1, trt2), n1 + n2, replace = F) 
  
  # initialize dataframe for permuted groups
  permuted.trt1 <- c()
  permuted.trt2 <- c()
  diff.permuted <- c()
  for (i in 1:n_perm){
    # shuffle and build permuted test stat for each group
    permuted.trt1 <- sample(rand.samples, length(n1), replace = F)
    permuted.trt2 <- sample(rand.samples, length(n2), replace = F)
    
    # calculated difference in means for each permutation
    diff.permuted <- append(diff.permuted, 
                            mean(permuted.trt1) - mean(permuted.trt2)) # permutation test stat
  }
  
  # assign p-value
  if(alternative[1]=="less"){
    pval = (sum(diff.permuted <= obs.diff) + 1) / (n_perm + 1)
  } else if(alternative[1]=="greater"){
    pval = (sum(diff.permuted >= obs.diff) + 1) / (n_perm + 1)
  } else{
    pval = (sum(abs(diff.permuted) >= abs(obs.diff)) + 1) / (n_perm + 1)
  }
  
  
  summary.test <- c(obs.diff, n1, n2, pval, pval < 0.05)
  return(summary.test)
}



# Function to permute based on SBGE category
# args: perm_dat = the dataset containing a numeric column to permute
#       x_col = string denoting denoting the column containing the data to permute
#       groupBy = logical column which separates group 1 and group 2
#       SBGE_cat = column specifying SBGE bins
TwoPerm_SBGE <- function(perm_dat, x_col, groupBy, SBGE_cat){
  dat <- data.frame() # initialize data.frame object to store results
  # make sure the column of categories is set to factor
  perm_dat[[SBGE_cat]] <- factor(perm_dat[[SBGE_cat]])
  
  # loop through each sex-bias category
  for(i in levels(perm_dat[[SBGE_cat]])){
    # for each SBGE category, call the permutation function
    comp <- c(i, TwoPerm(perm_dat[perm_dat[[SBGE_cat]] == i,], 
                           x_col = x_col, 
                           groupBy = groupBy, 
                           alternative = "two.sided"))
    
    # concatenate results in the data.frame to return
    dat <- rbind(dat, comp)
  }
  colnames(dat) <- c(SBGE_cat, "obs.diff", "n_TRUE", "n_FALSE", "pval", "Sig")
  dat$obs.diff <- as.numeric(dat$obs.diff)
  dat$n_TRUE <- as.numeric(dat$n_TRUE)
  dat$n_FALSE <- as.numeric(dat$n_FALSE)
  return(dat)
}



##########


# Bootstrapping functions
##########
# Bootstrapped 95% confidence intervals for two groups
# args: boot_dat = dataset to bootstrap
#       x_col = string denoting the column containing the data to bootstrap
#       groupBy = logical column which separates group 1 and group 2
#       boot_n = number of bootstrap replicates. Default is 1,000
#       myfun = test stat function. Default is mean
TwoBoot <- function(boot_dat, x_col, groupBy, 
                    boot_n = 1000, 
                    myfun = mean){
  boot_tabs <- as_tibble(data_frame(bs = 1:boot_n) %>% # make boot_n bootstrap replicates
                           group_by(bs) %>% # for each bootstrap,
                           # sample randomly x boot_n times from each group
                           mutate(data = list(boot_dat %>% 
                                                group_by(.[[groupBy]]) %>% # separate group 1 and group 2
                                                sample_frac(size = 1, replace = T)))) %>% 
    unnest(c(bs, data)) %>% # create separate rows for each bootstrap replicate in the list
    group_by(bs, .[[groupBy]]) %>% # group the data by bootstrap replicate and sig/non-sig
    # for each bootstrap replicate, calculate the test stat
    do(tidy(myfun(.[[x_col]]))) %>%
    rename({{groupBy}} := 2)
  
  # summarise bootstrap replicates
  boot_SE <- boot_tabs %>% 
    ungroup() %>%
    group_by(.[[groupBy]]) %>% 
    summarise(q05 = quantile(x, 0.025),
              q50 = quantile(x, 0.5),
              q95 = quantile(x, 0.975)) %>%
    rename({{groupBy}} := 1)
  return(boot_SE)
}

TwoBoot_SBGE <- function(boot_dat, x_col, groupBy, SBGE_cat){
  dat <- data.frame() # initialize data.frame object to store results
  # make sure the column of categories is set to factor
  boot_dat[[SBGE_cat]] <- factor(boot_dat[[SBGE_cat]])
  comp <- c()
  # loop through each sex-bias category
  for(i in levels(boot_dat[[SBGE_cat]])){
    # for each SBGE category, call the permutation function
    comp <- TwoBoot(boot_dat[boot_dat[[SBGE_cat]] == i,], 
                         x_col = x_col, 
                         groupBy = groupBy)
    comp[[SBGE_cat]] <- i
    # concatenate results in the data.frame to return
    dat <- rbind(dat, comp)
  }
  colnames(dat) <- c("Sig", "q05", "q50", "q95", SBGE_cat)
  return(dat)
}

##########


# testing
test_dat <- data.frame(dat = c(1,2,3,4,5,6,7,8,9, 10), group = rep(c(TRUE, FALSE)))

