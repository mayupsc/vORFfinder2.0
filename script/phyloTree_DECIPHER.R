suppressPackageStartupMessages(library(DECIPHER))

argv <- commandArgs(T)
nucleotide_seq_msa_file <- argv[1] 
constree_file <- argv[2]

nucleotide_sequences <- readDNAStringSet(nucleotide_seq_msa_file, format="fasta")
d <- DistanceMatrix(nucleotide_sequences,correction = "Jukes-Cantor",verbose=FALSE)
dend <- IdClusters(d,method="ML",type = "dendrogram",myXStringSet=nucleotide_sequences)
WriteDendrogram(dend,file=constree_file)

