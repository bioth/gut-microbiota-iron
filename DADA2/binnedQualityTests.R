library("dada2"); packageVersion("dada2")
library(ggplot2)
library(dplyr)

#weights and span and enforce monotonicity
loessErrfun_mod1 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

#weights and span and degree, also enforce monotonicity
loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

#testing my own parameters
#weights and span and degree, also enforce monotonicity
loessErrfun_modP <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

#microbiota 17 - testing stuff regarding binned quality scores
#loading filtered files of interest
setwd("~/Documents/CHUM_git/Microbiota_17/dada2_filtered_and_trimmed/")

#Listing files in the current working directory
list.files()

#Forward and reverse filtered fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
filtR1 <- sort(list.files(pattern="_R1_filt.fastq.gz", full.names = TRUE))
filtR2 <- sort(list.files(pattern="_R2_filt.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub("^(.*)_R1_filt\\.fastq\\.gz$", "\\1", basename(filtR1))

#dereplication of the fastq files
derepR1 <- derepFastq(filtR1)
derepR2 <- derepFastq(filtR2)
names(derepR1) <- sample_names
names(derepR2) <- sample_names

set.seed(100) # set seed to ensure that randomized steps are replicatable


#error rates estimation
errR1 <- learnErrors(derepR1, multithread=TRUE, randomize = TRUE, errorEstimationFunction = loessErrfun_modP)
errR2 <- learnErrors(derepR2, multithread=TRUE, randomize = TRUE, errorEstimationFunction = loessErrfun_mod1)
plotErrors(errR1, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/error_plotR1_m1.png", height = 10, width = 10)
plotErrors(errR2, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/error_plotR2_m1.png", height = 10, width = 10)


#running the inference algorithm => based on trimmed/filtered fastq and error rates
dadaR1 <- dada(derepR1, err=errR1, multithread=TRUE, pool = "pseudo", errorEstimationFunction = loessErrfun_mod1)
dadaR2 <- dada(derepR2, err=errR2, multithread=TRUE, pool = "pseudo", errorEstimationFunction = loessErrfun_mod1)
plotErrors(dadaR1, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/dada_plotR1_m1.png", height = 10, width = 10)
plotErrors(dadaR2, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/dada_plotR2_m1.png", height = 10, width = 10)

#merge paired reads 
mergers <- mergePairs(dadaR1, derepR1, dadaR2, derepR2, verbose=TRUE)

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#what percentage of chimears over the total dataset
sum(seqtab.nochim)/sum(seqtab)


#trying to export the data
setwd("../asv_table/")
write.table(seqtab.nochim, sep = ";", file = "seqtab.nochim_run_m1.csv", col.names = TRUE)



#tracking what reads made it through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaR1, getN), sapply(dadaR2, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
head(track)
track

mean(track[,4]/track[,1])*100
max(track[,4]/track[,1])*100
min(track[,4]/track[,1])*100



#extracting taxonomical information
#Loading ASV table previously generated
setwd("~/Documents/CHUM_git/Microbiota_17/asv_table/")
asv_table2 <- read.csv("seqtab.nochim_run_m1.csv", sep = ";")

#transforming asv_table into matrix so that it can be used by dada2 taxonomic assignment algorithm
asv_table <- as.matrix(asv_table2)
taxa <- assignTaxonomy(asv_table2, "../../training_set/silva_nr99_v138.1_wSpecies_train_set.fa.gz")
taxa.print <- taxa #removing rownames for display
rownames(taxa.print) <- NULL
head(taxa.print)


#save taxa matrix so that we can use it later
write.table(taxa, sep = ";", "../taxonomy/taxa_annotation_m1.csv", col.names = TRUE)


sink()


