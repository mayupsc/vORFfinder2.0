library(phylotools)

argv <- commandArgs(T)

virus_fasta_file <- argv[1]
virus_length_file <- argv[2] 

virus_fasta <- read.fasta(virus_fasta_file)

virus_fasta$len <- apply(virus_fasta,1,function(x){nchar(x[2])})
colnames(virus_fasta) <- c('virus','seq','length')

write.table(virus_fasta,virus_length_file,sep = "\t",quote = F,row.names = F,col.names = T)
