##-- 1. receive user input file and parameters
taskID=$1
file=$2
orf_len=$3
start_codon=$4
genetic_code=$5
ignore_nest=$6

mkdir -p  $taskID/blastp $taskID/database $taskID/ORF $taskID/log
chmod 777 $taskID/blastp $taskID/database $taskID/ORF $taskID/log


printf "your taskID is : "$taskID"\n"
printf "ORFfinder parameters are setting to:\n"
printf "ORF length : "$orf_len"\n"
printf "start codon : "$start_codon"\n"
printf "genetic code : "$genetic_code"\n"
printf "ignore_nest : "$ignore_nest"\n"
printf "your file including:\n"

unzip $file -d $taskID/database/virus_seq
mv $taskID/database/virus_seq/*/* $taskID/database/virus_seq

##--  2. some prepare work

#format fasta files included in zip file
for id in $(ls $taskID/database/virus_seq/*fasta); do ID=`echo $id|sed 's/.*\\///g'|sed 's/.fasta//g'`;sed -i "s/>.*/>$ID/g" $id; done

cat $taskID/database/virus_seq/*fasta > $taskID/database/virus_db.fasta
Rscript script/virus_length.R $taskID/database/virus_db.fasta $taskID/database/virus_length.txt

/opt/conda/bin/mafft --quiet $taskID/database/virus_db.fasta > $taskID/database/virus_msa.fasta

##-- 3. construct phylogenetic tree : DECIPHER
Rscript script/phyloTree_DECIPHER.R $taskID/database/virus_msa.fasta $taskID/database/constree> $taskID/log/step3_phyloTree_DECIPHER.log 
sed -i 's/\"//g' $taskID/database/constree

##-- 4. make ORF database

for id in $(ls $taskID/database/virus_seq/*fasta); do { sample=`echo ${id}|sed 's/.*\///g'|sed 's/.fasta//g'`; software/ORFfinder -in $id -ml $orf_len -s $start_codon -g $genetic_code -n $ignore_nest -out $taskID/ORF/${sample}.orfs >> $taskID/log/step4_make_orf_db.log; }& done ;wait
sed -i 's/lcl|//g' $taskID/ORF/*.orfs
sed -i 's/ unnamed protein product.*//g' $taskID/ORF/*.orfs

cat $taskID/ORF/*.orfs > $taskID/database/orf_db.fasta
/opt/conda/bin/makeblastdb -in $taskID/database/orf_db.fasta -parse_seqids -dbtype prot >> $taskID/log/step4_make_orf_db.log
##-- 5. blastp
for id in $(ls $taskID/ORF/*orfs); do { sample=`echo ${id}|sed 's/.*\///g'|sed 's/.orfs//g'`; /opt/conda/bin/blastp -query $id -db $taskID/database/orf_db.fasta -outfmt 10 -evalue 0.05 > $taskID/blastp/${sample}.blast ; }& done ; wait
for id in $(ls $taskID/blastp/*blast); do { Rscript script/blast_process.R $id $taskID >> $taskID/log/step5_blastp.log ; }& done ; wait
for id in $(ls $taskID/blastp/*fa); do { /opt/conda/bin/mafft --quiet $id > ${id%%.fa}.msa.fa; }& done; wait
#for id in $(ls blastp/*msa.fa); do Rscript script/msa_html.R $id>> log/step5_blastp.log ; done

## files for orf viewer
for id in $(ls $taskID/ORF/*orfs); do { Rscript script/orf_process.R $id $taskID>> $taskID/log/step6_orfviewer.log; }& done

