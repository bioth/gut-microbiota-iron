require("dada2"); packageVersion("dada2")
require(ggplot2)
require(dplyr)

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

set.seed(100) # Set seed to ensure that randomized steps are repeatable

# Loading filtered files of interest
setwd("path/to/trimmed_fastq/files")

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
r1_fastq <- sort(list.files(pattern="R1_trimmed.fastq", full.names = TRUE)) 
r2_fastq <- sort(list.files(pattern="R2_trimmed.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub(".*\\.(.*?)_R1_trimmed\\.fastq$", "\\1", basename(r1_fastq))

# Creating a dir for r console output and pictures
existingDirCheck("../r_console_output/")

# Printing aggregate quality profile plot for all of the R1 files
plotQualityProfile(r1_fastq, aggregate = TRUE)
ggsave(dpi = 300, filename = "../r_console_output/quality_profileR1.png", height = 6, width = 10)

# Printing aggregate quality profile plot for all of the R2 files
plotQualityProfile(r2_fastq, aggregate = TRUE)
ggsave(dpi = 300, filename = "../r_console_output/quality_profileR2.png", height = 6, width = 10)

# Quality filtering of the reads
# Place filtered files in filtered/ subdirectory
existingDirCheck("../dada2_filtered_and_trimmed/")
path = "../dada2_filtered_and_trimmed"

filtR1 <- file.path(path, paste0(sample_names, "_R1_filt.fastq.gz"))
filtR2 <- file.path(path, paste0(sample_names, "_R2_filt.fastq.gz"))
names(filtR1) <- sample_names
names(filtR2) <- sample_names

# Start redirecting output to a file
run <- stringr::str_replace_all(as.character(Sys.time()), "[ :]", "-")
dataset <- "optimization"
run <- paste0("../r_console_output/", run, dataset,".txt")
sink(run)

out <- filterAndTrim(r1_fastq, filtR1, r2_fastq, filtR2, truncLen=c(270,250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE)

#Check how many reads made it out
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
errR1 <- learnErrors(derepR1, multithread=TRUE, randomize = TRUE, errorEstimationFunction = loessErrfun_mod1) 
errR2 <- learnErrors(derepR2, multithread=TRUE, randomize = TRUE, errorEstimationFunction = loessErrfun_mod1)
plotErrors(errR1, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/error_plotR1.png", height = 10, width = 10)
plotErrors(errR2, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/error_plotR2.png", height = 10, width = 10)

# Running the inference algorithm => based on trimmed/filtered fastq and error rates
dadaR1 <- dada(derepR1, err=errR1, multithread=TRUE, pool = TRUE, errorEstimationFunction = loessErrfun_mod1) 
dadaR2 <- dada(derepR2, err=errR2, multithread=TRUE, pool = TRUE, errorEstimationFunction = loessErrfun_mod1)
plotErrors(dadaR1, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/dada_plotR1.png", height = 10, width = 10)
plotErrors(dadaR2, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/dada_plotR2.png", height = 10, width = 10)

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
colnames(track) <- c("InputF", "InputR", "FilteredF", "FilteredR", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
print("Tracking of reads that made it through the pipeline:")
print(track)

# Stop writing content in the output file
sink()