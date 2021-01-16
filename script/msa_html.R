library(Biostrings)
library(DECIPHER)

argv <- commandArgs(T)
protein_seq_msa_file <- argv[1]
html_file <- gsub('.msa.fa','.html',protein_seq_msa_file)
protein_sequences <- readAAStringSet(protein_seq_msa_file, format="fasta")

patterns = c("-", alphabet(protein_sequences))
BrowseSeqs(protein_sequences, colorPatterns=T,colors = rainbow(length(patterns)),htmlFile = html_file,openURL = F)
