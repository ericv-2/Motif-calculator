##Motifs exercise by Èric Vargas Riera

#Select your motifs
motifs <- c("TTGCAAAT",
            "ATGCATAT",
            "ATTCAAAT",
            "ATGCAAAA",
            "ATGCAAAT",
            "ATGCAATT",
            "ATGCAATT",
            "AAGCAAAG")

NUC <- c("A", "C", "G", "T")
motif_length <- length(strsplit(motifs[1], split = "")[[1]])

#Build the PPSM matrix
PPSM <- matrix(nrow = length(NUC), ncol = motif_length)
rownames(PPSM) <- NUC
colnames(PPSM) <- c({1:motif_length})
PPSM[,] <- 0

for (i in motifs){
  change = 0
  j <- strsplit(i, split = "")
  for (nucleotide in j[[1]]){
    change = change + 1
    if (nucleotide == "A")
      PPSM[1, change] = PPSM[1, change] + 1
    else if (nucleotide == "C")
      PPSM[2, change] = PPSM[2, change] + 1
    else if (nucleotide == "G")
      PPSM[3, change] = PPSM[3, change] + 1
    else if (nucleotide == "T")
      PPSM[4, change] = PPSM[4, change] + 1
    else
      break
  }
}

#Build PWM matrix
PWM <- matrix(nrow = length(NUC), ncol = motif_length)
rownames(PWM) <- c(NUC)
colnames(PWM) <- c({1:8})
PWM[,] <- 0

for (k in 1:ncol(PPSM)){
  total <- sum(PPSM[1:4, k])
  PWM[1:4, k] <- PPSM[1:4, k]/total
}

#Pseudocount matrix
PSE <- matrix(nrow = length(NUC), ncol = motif_length)
rownames(PSE) <- c(NUC)
colnames(PSE) <- c({1:8})
PSE[,] <- 0 
pseudocounts = 1

PSE <- PPSM + pseudocounts/length(NUC)
for (l in 1:ncol(PPSM)){
  total <- sum(PPSM[1:4, l])
  PSE[1:4, l] <- PSE[1:4, l]/(total + pseudocounts)
}

#log-likelihood matrix
LOG <- matrix(nrow = length(NUC), ncol = motif_length)
rownames(LOG) <- c(NUC)
colnames(LOG) <- c({1:8})
LOG[,] <- 0 
val <- 1/length(NUC)

for (m in 1:ncol(PSE)){
  for (n in 1:nrow(PSE)){
    LOG[n,m] <- log(PSE[n,m]/val)
  }
}

#Information content matrix
IC_total = log2(length(NUC))
IC <- matrix(nrow = length(NUC), ncol = motif_length)
rownames(IC) <- c(NUC)
colnames(IC) <- c({1:8})
IC[,] <- 0 

for (o in 1:ncol(IC)){
  U = 0
  for (p in 1:nrow(IC)){
    log2_ <- log2(PWM[p,o])
    if (log2_ == -Inf){
      log2_ = 0
      U = ((PWM[p,o] * log2_) + U)
    }
    else {
      U = ((PWM[p,o] * log2_) + U)
    }
  }
  U = - U
  IC_final = IC_total - U
  for (q in 1:nrow(IC)){
    IC[q, o] <- PWM[q, o] * IC_final
  }
}

cat("The PPSM matrix is:\n") 
print(PPSM)

cat("The PWM matrix is:\n") 
print(round(PWM, digits = 2))

cat("The pseudocount matrix is:\n") 
print(round(PSE, digits = 2))

cat("The log-likelihood matrix is:\n") 
print(round(LOG, digits = 2))

cat("The Information content matrix is:\n") 
print(round(IC, digits = 2))


