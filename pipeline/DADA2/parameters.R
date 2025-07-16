mergers <- mergePairs(dadaR1, derepR1, dadaR2, derepR2, verbose=TRUE, maxMismatch = 0, justConcatenate = FALSE)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, allowOneOff = FALSE) # Much more conservative, not so good when having batch effects and timepoints data
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, allowOneOff = FALSE)
# M18
260,240

# M19
250,230
                                                
                                                
                                              