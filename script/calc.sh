##-- 0. init
rm -r database blastp ORF
mkdir -p  blastp database/par  ORF
chmod 777 blastp database ORF

##-- 1. receive user input file and parameters
file=$1
orf_len=$2
start_codon=$3
genetic_code=$4
ignore_nest=$5

filename=`basename $file|sed 's/.zip//g'`

printf "ORFfinder parameters are setting to:\n"
printf "ORF length : "$orf_len"\n"
printf "start codon : "$start_codon"\n"
printf "genetic code : "$genetic_code"\n"
printf "ignore_nest : "$ignore_nest"\n"
printf "your file is "$file"\nfiles including:\n"

unzip $file -d database/virus_seq
mv database/virus_seq/*/* database/virus_seq

##--  2. some prepare work 

#format fasta files included in zip file
for id in $(ls database/virus_seq/*fasta); do ID=`echo $id|sed 's/.*\\///g'|sed 's/.fasta//g'`;sed -i "s/>.*/>$ID/g" $id; done

cat database/virus_seq/*fasta > database/virus_db.fasta
Rscript /srv/shiny-server/script/virus_length.R

/opt/conda/bin/mafft database/virus_db.fasta > database/virus_msa.fasta


##-- 3. construct phylogenetic tree : DECIPHER
Rscript /srv/shiny-server/script/phyloTree_DECIPHER.R
sed -i 's/\"//g' database/constree

##-- 4. make ORF database

for id in $(ls database/virus_seq/*fasta); do { sample=`echo ${id}|sed 's/.*\///g'|sed 's/.fasta//g'`; /srv/shiny-server/software/ORFfinder -in $id -ml $orf_len -s $start_codon -g $genetic_code -n $ignore_nest -out ORF/${sample}.orfs; }& done ;wait
sed -i 's/lcl|//g' ORF/*.orfs
sed -i 's/ unnamed protein product.*//g' ORF/*.orfs

cat ORF/*.orfs > database/orf_db.fasta
/opt/conda/bin/makeblastdb -in database/orf_db.fasta -parse_seqids -dbtype prot

##-- 5. /opt/conda/bin/blastp
for id in $(ls ORF/*orfs); do { sample=`echo ${id}|sed 's/.*\///g'|sed 's/.orfs//g'`; /opt/conda/bin/blastp -query $id -db database/orf_db.fasta -outfmt 10 -evalue 0.05 > blastp/${sample}.blast; }& done ; wait
for id in $(ls blastp/*blast); do { Rscript /srv/shiny-server/script/blast_process.R $id ; }& done ; wait
for id in $(ls blastp/*fa); do { /opt/conda/bin/mafft $id > ${id%%.fa}.msa.fa; }& done; wait
#for id in $(ls blastp/*msa.fa); do Rscript /srv/shiny-server/script/msa_html.R $id; done


## files for orf viewer
for id in $(ls ORF/*orfs); do { Rscript /srv/shiny-server/script/orf_process.R $id; }& done

