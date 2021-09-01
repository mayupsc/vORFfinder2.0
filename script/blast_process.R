suppressPackageStartupMessages(library(phylotools))

argv <- commandArgs(T)
usr_interest_blast <- argv[1]
taskID <- argv[2]

blast_result <- read.table(usr_interest_blast,sep = ",",stringsAsFactors = F)
colnames(blast_result) <- c('queryID','subjectID','identity','Alignment_length','mismatch','gap_opens','q.start','q.end','s.start','s.end','evalue','bit_score')
#blast_result <- blast_result[which(blast_result$queryID!=blast_result$subjectID),]

orf_db <- read.fasta(paste0(taskID,"/database/orf_db.fasta"))

blast_result$querySeq <- orf_db$seq.text[match(blast_result$queryID,orf_db$seq.name)]
blast_result$subjectSeq <- orf_db$seq.text[match(blast_result$subjectID,orf_db$seq.name)]


blast_result_split <- split(blast_result,blast_result$queryID)

for (i in 1:length(blast_result_split)) {
  
  blast_result_i <- blast_result_split[[i]]
  blast_result_i <- blast_result_i[order(blast_result_i$evalue,decreasing = F),]
  ORF_name <- gsub(':.*','',unique(blast_result_i$queryID))
  write.table(blast_result_i,paste0(taskID,'/blastp/',ORF_name,'.blastp'),quote = F,row.names = F,col.names = T,sep = "\t")
}
