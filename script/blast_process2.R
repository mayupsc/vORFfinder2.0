suppressPackageStartupMessages(library(phylotools))
suppressPackageStartupMessages(library(jsonlite))

argv <- commandArgs(T)
usr_interest_blast <- argv[1]
taskID <- argv[2]

blast_json <- fromJSON(usr_interest_blast)

hits_number <- c()
for (i in 1:length(blast_json$BlastOutput2$report$program)) {
  hits_number<- c(hits_number,dim(blast_json$BlastOutput2$report$results$search$hits[[i]])[1])
}


orf_number <- length(blast_json$BlastOutput2$report$program)

for (i in 1:orf_number) {
  ORF_name <- gsub('\\:.*','',blast_json$BlastOutput2$report$results$search$query_title[i])
  orf_hits <- blast_json$BlastOutput2$report$results$search$hits[[i]] ###abbreviate code
  
  for (j in 1:dim(orf_hits)[1]){
    id <- orf_hits$description[[j]]$id
    #####add space or query showing in a wrong way
    n_gap_fill <- paste0(rep("-",orf_hits$hsps[[j]]$query_from-1),collapse = "")
    query <- paste0(n_gap_fill,orf_hits$hsps[[j]]$hseq)
    query_title <- paste0(">",id)
    
    system(paste0("echo '",query_title,"' >> ", paste0(taskID,'/blastp/',ORF_name,'.fa')))
    system(paste0("echo '",query,"' >> ",paste0(taskID,'/blastp/',ORF_name,'.fa')))
      
  }
}
