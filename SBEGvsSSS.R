# Looking at SSS vs SBGE

# More SSS exons are female-biased?
# seems to not be true here...
MF.expressed <- jseq.ASE[!is.na(jseq.ASE$expr_M) & !is.na(jseq.ASE$expr_F) &
                           jseq.ASE$FlyBaseID %in% filter.low.exp.genes.q25,]
dim(MF.expressed)
hist(MF.expressed[MF.expressed$SSS,]$log2FCvst.M.F., breaks = 100)
length(MF.expressed[MF.expressed$SSS,]$FlyBaseID[MF.expressed$log2FCvst.M.F. < -0.5])
length(MF.expressed[MF.expressed$SSS,]$FlyBaseID[MF.expressed$log2FCvst.M.F. > 0.5])


test <- merge(male.exp_ASE[,c("FlyBaseID", "countBin", "glob.exp", "frac.exp.per.gene")], 
              fem.exp_ASE[,c("FlyBaseID", "countBin", "glob.exp", "frac.exp.per.gene")], 
              by = c("FlyBaseID", "countBin"))
colnames(test)[3:6] <- c("glob.exp.M", "frac.exp.M", "glob.exp.F", "frac.exp.F")
t.test(test$glob.exp.M, test$glob.exp.F, paired = T)
sum(test$frac.exp.M < test$frac.exp.F, na.rm = T)
sum(test$frac.exp.M > test$frac.exp.F, na.rm = T)

test <- merge(male.exp_SinghAgrw[,c("FlyBaseID", "countBin", "glob.exp", "frac.exp.per.gene")], 
              fem.exp_SinghAgrw[,c("FlyBaseID", "countBin", "glob.exp", "frac.exp.per.gene")], 
              by = c("FlyBaseID", "countBin"))
colnames(test)[3:6] <- c("glob.exp.M", "frac.exp.M", "glob.exp.F", "frac.exp.F")
t.test(test[test$FlyBaseID %in% subset.sss,]$frac.exp.M, test[test$FlyBaseID %in% subset.sss,]$frac.exp.F, paired = T)
sum(test[test$FlyBaseID %in% subset.sss,]$frac.exp.M < test[test$FlyBaseID %in% subset.sss,]$frac.exp.F, na.rm = T)
sum(test[test$FlyBaseID %in% subset.sss,]$frac.exp.M > test[test$FlyBaseID %in% subset.sss,]$frac.exp.F , na.rm = T)
hist(test$frac.exp.M - test$frac.exp.F, breaks = 100)

ReadCountFilesCalcNumExons <- function(data.files, decoder.file){
  count.file <- NULL
  for(i in data.files){
    tmp <- read.delim(i, header=F, sep = "\t")
    count.file <- cbind(count.file, tmp[,2])
  }
  
  countBins <- tmp %>% separate(V1, into = c("FlyBaseID", "countBin"), sep = ":")
  count.file <- cbind(countBins[,1:2], count.file)
  colnames(count.file) <- c("FlyBaseID", "countBin", decoder.file$unique.ID)
  print(dim(count.file)) # check number of splice sites
  
  # exclude novel splice sites
  count.file <- count.file %>% filter(!str_detect(countBin, "N") & !str_detect(countBin, "A"))
  print(dim(count.file)) # check number of remaining splice sites
  
  count.file <- count.file %>% group_by(FlyBaseID) %>%
    mutate(numExons = n())
  
  return(count.file)
}

numExons <- jseq.A.f.geno[,c(2,3)] %>% group_by(FlyBaseID) %>%
  summarise(numExon = n())

test <- merge(jseq.A.m.geno, numExons, by = c("FlyBaseID"))
summary(glm(sig.hit ~ numExon, data= test, family = "binomial"))
ggplot(test %>% mutate(prob=ifelse(sig.hit, 1, 0)), aes(numExon, prob)) + 
  stat_smooth(formula = y ~ x, method = "glm", 
              method.args = list(family = "binomial"))

test <- merge(A.m.NR.compare_ase, numExons, by = c("FlyBaseID"))
summary(glm(M.dis.F ~ numExon, data= test))
ggplot(test, aes(numExon, M.sub.F)) + 
  stat_smooth(formula = y ~ x, method = "glm", 
              method.args = list(family = "gaussian"))


# genes that are significant in A females tend to be also consistently dimorphic
# genes that are significant in A males not associated with the subset of consistently dimorphic genes
test <- numExons %>% 
  dplyr::mutate(subset = ifelse(FlyBaseID %in% subset.sss, TRUE, FALSE),
                sig.AS= ifelse(FlyBaseID %in% unique(c(A.m.sig.AS, A.f.sig.AS)), TRUE, FALSE))
fisher.test(test$sig.AS, test$subset)
# genes that are longer are more likely to be called significant, but they are also more likely to be
# consistently dimorphic between different lab populations.
summary(glm(subset.sss ~ numExon, data = test, family = "binomial"))


