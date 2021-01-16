library(phylotools)

aligned_fa_file <- 'database/virus_msa.fasta'
out_phy_file <- 'database/virus.phy'
phy_fa <- read.fasta(aligned_fa_file)
ID_convert <- data.frame('seq.name'=phy_fa$seq.name,'ID'=paste0('spe_',1:nrow(phy_fa)))
write.table(ID_convert,'database/ID_convert.txt',sep = "\t",quote = F,row.names = F,col.names = T)

phy_fa$seq.name <- ID_convert$ID[match(phy_fa$seq.name,ID_convert$seq.name)]
#phy_fa$seq.name <- paste0(phy_fa$seq.name,rep(" ",10))
dat2phylip(phy_fa,outfile = out_phy_file)
system(paste0("awk '{printf(\"\\%-10s\\%s\\n\",$1,$2)}' ", aligned_fa_file," > database/tmp"))
system(paste0("mv ",aligned_fa_file, " database/tmp"))