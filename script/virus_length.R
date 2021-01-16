library(phylotools)

virus_fasta_file <- 'database/virus_db.fasta'
virus_fasta <- read.fasta(virus_fasta_file)

virus_fasta$len <- apply(virus_fasta,1,function(x){nchar(x[2])})
colnames(virus_fasta) <- c('virus','seq','length')

write.table(virus_fasta,'database/virus_length.txt',sep = "\t",quote = F,row.names = F,col.names = T)
