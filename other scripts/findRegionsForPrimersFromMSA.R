library(Biostrings)

# Convert an aligned coordinate (with gaps) to the true ungapped sequence position
aligned_to_ungapped_pos <- function(aligned_seq, aln_pos) {
  aligned_seq <- paste0(aligned_seq, collapse = "")
  ungapped_pos <- 0L
  for (i in seq_len(aln_pos)) {
    if (substr(aligned_seq, i, i) != "-") ungapped_pos <- ungapped_pos + 1L
  }
  return(ungapped_pos)
}

# Identify windows where the target is specific (as before)
find_specific_windows <- function(fasta_file, window_size = 20, step = 1, max_mismatch = 3) {
  aln <- readDNAStringSet(fasta_file)
  aln_mat <- as.matrix(aln)
  target <- aln_mat[1, ]
  n_seqs <- nrow(aln_mat)
  aln_length <- ncol(aln_mat)
  starts <- integer()
  
  for (s in seq(1, aln_length - window_size + 1, by = step)) {
    w <- s:(s + window_size - 1)
    tw <- target[w]
    if ("-" %in% tw) next
    # check mismatches
    for (i in 2:n_seqs) {
      cw <- aln_mat[i, w]
      if ("-" %in% cw) next
      if (sum(tw != cw) <= max_mismatch) {
        tw <- NULL
        break
      }
    }
    if (!is.null(tw)) starts <- c(starts, s)
  }
  aligned_seq <- paste0(as.character(readDNAStringSet(fasta_file)[1]), collapse = "")
  list(starts = starts,
       window_size = window_size,
       aligned_seq = aligned_seq)
}

# Now build excluded regions from the complement of those starts
exclude_non_specific_regions <- function(fasta_file,
                                         window_size = 20,
                                         step = 1,
                                         max_mismatch = 3,
                                         primer_len = 20) {
  spec <- find_specific_windows(fasta_file, window_size, step, max_mismatch)
  aln_seq <- spec$aligned_seq
  aln_len <- nchar(aln_seq)
  # all possible window starts
  all_starts <- seq(1, aln_len - window_size + 1, by = step)
  forbidden_starts <- setdiff(all_starts, spec$starts)
  if (length(forbidden_starts) == 0) {
    message("No forbidden regions: everything is specific.")
    return(invisible(NULL))
  }
  
  # Group into contiguous runs
  groups <- split(forbidden_starts, cumsum(c(1, diff(forbidden_starts) != step)))
  
  excluded_blocks <- lapply(groups, function(g) {
    aln_start <- min(g)                             # first aligned pos
    aln_end   <- max(g) + window_size - 1L          # last aligned pos of that window
    # expand on both ends so no primer of length primer_len can overlap
    aln_start <- max(1, aln_start - (primer_len - 1L))
    aln_end   <- min(aln_len, aln_end   + (primer_len - 1L))
    # convert to ungapped coordinates
    start_true <- aligned_to_ungapped_pos(aln_seq, aln_start)
    end_true   <- aligned_to_ungapped_pos(aln_seq, aln_end)
    length_true <- end_true - start_true + 1L
    c(start_true, length_true)
  })
  
  # collapse for Primer3
  reg_str <- paste(vapply(excluded_blocks,
                          function(x) paste0(x[1], ",", x[2]),
                          character(1)),
                   collapse = " ")
  
  cat("SEQUENCE_PRIMER_EXCLUDED_REGION=", reg_str, "\n\n", sep = "")
  
  # annotate with <...>
  ungapped <- strsplit(gsub("-", "", aln_seq), "")[[1]]
  mask <- rep(FALSE, length(ungapped))
  for (blk in excluded_blocks) {
    idx <- blk[1]:(blk[1] + blk[2] - 1L)
    mask[idx] <- TRUE
  }
  
  annotated <- character(0)
  i <- 1L
  while (i <= length(ungapped)) {
    if (mask[i]) {
      annotated <- c(annotated, "<")
      while (i <= length(ungapped) && mask[i]) {
        annotated <- c(annotated, ungapped[i])
        i <- i + 1L
      }
      annotated <- c(annotated, ">")
    } else {
      annotated <- c(annotated, ungapped[i])
      i <- i + 1L
    }
  }
  
  cat("Annotated sequence:\n")
  cat(paste0(annotated, collapse = ""), "\n")
  
  invisible(list(excluded = excluded_blocks,
                 regions_str = reg_str,
                 annotated = paste0(annotated, collapse = "")))
}

find_specific_windows(fasta_file =  "~/Downloads/clustalo-I20250513-192439-0258-24961506-p1m - Copie.fasta", window_size = 20)

# — Example usage —
exclude_non_specific_regions(
  fasta_file  = "~/Downloads/clustalo-I20250513-192439-0258-24961506-p1m - Copie.fasta",
  window_size = 20,
  step        = 1,
  max_mismatch= 3,
  primer_len  = 20
)







sequence = "AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCATGCCTAATACATGCAAGTCGAACGAGAGACCTT 
CGGGTCTCTAGTGGCGAACGGGTGAGTAACACGTAGGGAACCTGCCCGCGCACCGGGAATACGCTCTGGA 
AACGGAGAACAAATCCCGATGTACAGGAAGGAGGCATCTTCTTTCTGTGAAACATCCTTTAGGGGATGGG 
GCGCGGATGGACCTGCGGTGCATTAGTTGGTTGGCGGGGTAAAGGCCCACCAAGACGATGATGCATAGCC 
GGCCTGAGAGGGCGGACGGCCACATTGGGACTGAGACACGGCCCAGACTCCTGCGGGAGGCAGCAGTAGG 
GAATTTTCGTCAATGGGCGCAAGCCTGAACGAGCGATGCCGCGTGAGTGAAGAAGGCCTTCGGGTCGTAA 
AGCTCTGTTGCGGGGGAAAAAAGGCAGCATCAGGAAATGGGTGCTGACTGATGGTGCCCCGCCAGAAAGT 
CACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTATCCGGAATGATTGGGCGT 
AAAGGGTGCGCAGGCGGTCCTGCAAGTCTGGAGTGAAACGCATGAGCTCAACTCATGCATGGCTTTGGAA 
ACTGGAGGACTGGAGAGCAGGAGAGGGCGGTGGAACTCCATGTGTAGCGGTAAAATGCGTAGATATATGG 
AAGAACACCAGTGGCGAAGGCGGCCGCCTGGCCTGTTGCTGACGCTGAGGCACGAAAGCGTGGGGAGCAA 
ATAGGATTAGATACCCTAGTAGTCCACGCCGTAAACGATGAGGACCAAGTGTTGGGGGTGAAACCTCAGT 
GCTGAAGTTAACGCAGTGAGTCCTCCGCCTGGGGAGTATGCACGCAAGTGTGAAACTCAAAGGAATTGAC 
GGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGCCTTGA 
CATGGGATGCGAAGATGCAGAGATGCATCGGAGGTCAACATCCACACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTCAAGTCCCGCAACGAGCGCAACCCTTGTGGCATGTTGCTAACAGGAAA 
AGCTGAGGACTCATGCCAGACTGCCGGTGACAAACCGGAGGAAGGCGGGGATGACGTCAAATCATCATGC 
CCCTTATGGCCTGGGCTACACACGTACTACAATGGCGGCTACAAAGAGCAGCGAGACAGGGATGTCGAGC 
GAATCTCATAAAAGCCGTCCCAGTTCGGATTGGAGGCTGCAACCCGCCTCCATGAAGTTGGAATCGCTAG 
TAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAACCATGGG 
AGTCGGTAATGCCCGAAGCCGGTGGCATGACCTCATAAGAGGAGTGAGCCGTCGAAGGCAGGATCGATGA 
CTGGGGTTAAGTCGTAA"
sequence = gsub(" \n", "", sequence)
substring(sequence, first = 65, last = 84)
