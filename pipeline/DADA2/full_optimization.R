library("dada2"); packageVersion("dada2")
library(ggplot2)
library(dplyr)
library(ShortRead)
library(data.table)

# Function checking if a dir exists and creating it otherwise
existingDirCheck <- function(path){
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Directory created: ", path)
  } else {
    message("Directory already exists: ", path)
  }
  
}

# Custom functions for error modeling 
{
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
}

# Microbiota 18
# Setting working directory and selecting loading files of interest
setwd("~/Documents/CHUM_git/Microbiota_18/trimmed_fastq/")
setwd("~/Documents/CHUM_git/Microbiota_18/2106590411/")

# Loading metdata
meta <- read.csv("../../gut-microbiota-iron/experiments/finished exp/young-DSS-exp3/young48_dss_dissection.csv", header = TRUE, sep = ";")
meta$id <- substring(meta$id, 1, 5)
meta$gg_group <- paste(meta$diet,meta$treatment,sep=":")
samples <- c()
samples <- append(samples, meta[meta$gg_group == "50:water","id"][1:2])
samples <- append(samples, meta[meta$gg_group == "500:water","id"][1:2])
samples <- append(samples, meta[meta$gg_group == "50:dss","id"][1:2])
samples <- append(samples, meta[meta$gg_group == "500:dss","id"][1:2])
timepoints <- c("T0","T35","T49","T54","Tfinal")
samples <- as.vector(outer(samples, timepoints, paste, sep = "_"))

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
r1_fastq <- sort(list.files(pattern="R1_trimmed.fastq", full.names = TRUE)) 
r2_fastq <- sort(list.files(pattern="R2_trimmed.fastq", full.names = TRUE))

r1_fastq <- sort(list.files(pattern="R1.fastq", full.names = TRUE)) 
r2_fastq <- sort(list.files(pattern="R2.fastq", full.names = TRUE))

r1_fastq <- r1_fastq[grepl(paste0(paste0(samples, "s*"), collapse = "|"), r1_fastq)]
r2_fastq <- r2_fastq[grepl(paste0(paste0(samples, "s*"), collapse = "|"), r2_fastq)]

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub(".*\\.(.*?)_R1_trimmed\\.fastq$", "\\1", basename(r1_fastq))

# Creating a dir for r console output and pictures
existingDirCheck("~/Documents/CHUM_git/full_optimization/r_console_output/")

# Printing aggregate quality profile plot for all of the R1 files
plotQualityProfile(r1_fastq, aggregate = TRUE)
ggsave(dpi = 300, filename = "~/Documents/CHUM_git/full_optimization/r_console_output/quality_profileR1.png", height = 6, width = 10)

# Printing aggregate quality profile plot for all of the R2 files
plotQualityProfile(r2_fastq, aggregate = TRUE)
ggsave(dpi = 300, filename = "~/Documents/CHUM_git/full_optimization/r_console_output/quality_profileR2.png", height = 6, width = 10)

# Quality filtering of the reads
# Place filtered files in filtered/ subdirectory
existingDirCheck("~/Documents/CHUM_git/full_optimization/dada2_filtered_and_trimmed/")
path = "~/Documents/CHUM_git/full_optimization/dada2_filtered_and_trimmed"

filtR1 <- file.path(path, paste0(sample_names, "_R1_filt.fastq.gz"))
filtR2 <- file.path(path, paste0(sample_names, "_R2_filt.fastq.gz"))
names(filtR1) <- sample_names
names(filtR2) <- sample_names

# Each read is 300 bp. R1 - 18 bases primers, 282
# R2 - 15 bases primer, 285. truncLen=c(270,250)
# Start redirecting output to a file
run <- stringr::str_replace_all(as.character(Sys.time()), "[ :]", "-")
dataset <- "optimization"
run <- paste0("~/Documents/CHUM_git/full_optimization/r_console_output/", run, dataset,".txt")
sink(run)

out <- filterAndTrim(r1_fastq, filtR1, r2_fastq, filtR2, truncLen=c(270,250),trimLeft = c(18,15),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE) # On Windows set multithread=FALSE

#Check how many reads made it out
print(out)
print("Mean number of reads that passed filter:")
print(mean(out[,2]/out[,1])*100)
print("Min:")
print(min(out[,2]/out[,1])*100)
print("Max:")
print(max(out[,2]/out[,1])*100)

# Dereplication of the fastq files
derepR1 <- derepFastq(filtR1)
derepR2 <- derepFastq(filtR2)
names(derepR1) <- sample_names
names(derepR2) <- sample_names

# Error rates estimation
errR1 <- learnErrors(derepR1, multithread=TRUE, randomize = TRUE, errorEstimationFunction = loessErrfun_mod1) # , errorEstimationFunction = loessErrfun_mod4
errR2 <- learnErrors(derepR2, multithread=TRUE, randomize = TRUE, errorEstimationFunction = loessErrfun_mod1) # , errorEstimationFunction = loessErrfun_mod4
plotErrors(errR1, nominalQ=TRUE)
ggsave(dpi = 300, filename = "~/Documents/CHUM_git/full_optimization/r_console_output/error_plotR1_m1.png", height = 10, width = 10)
plotErrors(errR2, nominalQ=TRUE)
ggsave(dpi = 300, filename = "~/Documents/CHUM_git/full_optimization/r_console_output/error_plotR2_m1.png", height = 10, width = 10)

# Running the inference algorithm => based on trimmed/filtered fastq and error rates
dadaR1 <- dada(derepR1, err=errR1, multithread=TRUE, pool = TRUE, errorEstimationFunction = loessErrfun_mod1) # , errorEstimationFunction = loessErrfun_mod4
dadaR2 <- dada(derepR2, err=errR2, multithread=TRUE, pool = TRUE, errorEstimationFunction = loessErrfun_mod1) # , errorEstimationFunction = loessErrfun_mod4
plotErrors(dadaR1, nominalQ=TRUE)
ggsave(dpi = 300, filename = "~/Documents/CHUM_git/full_optimization/r_console_output/dada_plotR1_m1.png", height = 10, width = 10)
plotErrors(dadaR2, nominalQ=TRUE)
ggsave(dpi = 300, filename = "~/Documents/CHUM_git/full_optimization/r_console_output/dada_plotR2_m1.png", height = 10, width = 10)

# Merge paired reads 
mergers <- mergePairs(dadaR1, derepR1, dadaR2, derepR2, verbose=TRUE, maxMismatch = 0, justConcatenate = FALSE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
print("Dim of the sequence table:")
dim(seqtab)

# Inspect distribution of sequence lengths
print("Distribution of sequence lengths:")
table(nchar(getSequences(seqtab)))

# Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE,
                                    allowOneOff=FALSE, minFoldParentOverAbundance=8)

# What percentage of chimeras over the total dataset
print("Percentage of chimeras over the total dataset:")
1-sum(seqtab.nochim)/sum(seqtab)

table(nchar(getSequences(seqtab.nochim)))

# Tracking what reads made it through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(r1_fastq, getN), sapply(r2_fastq, getN), sapply(filtR1, getN), sapply(filtR2, getN), sapply(dadaR1, getN), sapply(dadaR2, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("Input", "InputR", "FilteredF", "FilteredR", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
print("Tracking of reads that made it through the pipeline:")
track

print("Mean %:")
mean(track[,4]/track[,3])*100
print("Max %:")
max(track[,4]/track[,3])*100
print("Min %:")
min(track[,4]/track[,3])*100

# Stop writing things in the output file
sink()