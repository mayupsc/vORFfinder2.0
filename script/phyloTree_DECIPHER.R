library(DECIPHER)
nucleotide_seq_msa_file <- 'database/virus_msa.fasta'

nucleotide_sequences <- readDNAStringSet(nucleotide_seq_msa_file, format="fasta")
d <- DistanceMatrix(nucleotide_sequences,correction = "Jukes-Cantor",verbose=FALSE)
dend <- IdClusters(d,method="ML",type = "dendrogram",myXStringSet=nucleotide_sequences)
WriteDendrogram(dend,file="database/constree")

