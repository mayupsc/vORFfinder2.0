library(phylotools)

argv <- commandArgs(T)

orf_file <- argv[1]
usr_species <- gsub('.orfs','',basename(orf_file))

orf_seq <- read.fasta(orf_file)

orf_seq$start <- as.numeric(apply(orf_seq,1,function(x){strsplit(x[1],":")[[1]][2]}))
orf_seq$end <- as.numeric(apply(orf_seq,1,function(x){strsplit(x[1],":")[[1]][3]}))
orf_seq$symbol <- gsub('_.*','',apply(orf_seq,1,function(x){strsplit(x[1],":")[[1]][1]}))
orf_seq$width <- abs(orf_seq$end - orf_seq$start)

orf_seq$strand <- '+'
orf_seq$strand[orf_seq$end<orf_seq$start] <- '-'

orf_seq$start2 <- orf_seq$start
orf_seq$start2[orf_seq$strand=="-"] <- orf_seq$end[orf_seq$strand=="-"]


## hits number
blast <- read.table(paste0("blastp/",usr_species,".blast"),sep = ",",stringsAsFactors = F)
colnames(blast) <- c('queryID','subjectID','identity','Alignment_length','mismatch','gap_opens','q.start','q.end','s.start','s.end','evalue','bit_score')

blast_hits_number <- aggregate(blast$subjectID,by=list(blast$queryID),length)
colnames(blast_hits_number) <- c('queryID','number')

orf_seq$hits <- blast_hits_number$number[match(orf_seq$seq.name,blast_hits_number$queryID)]


write.table(orf_seq,paste0('ORF/',usr_species,".orfviewer.txt"),sep = "\t",quote = F,row.names = F,col.names = T)
