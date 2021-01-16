
##-- 0. init
rm -r database blastp ORF
mkdir -p  blastp database/par  ORF
chmod 777 blastp database ORF output

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

mafft database/virus_db.fasta > database/virus_msa.fasta


##-- 3. construct phylogenetic tree : phylip

Rscript /srv/shiny-server/script/fa2phy.R
awk '{printf("%-10s%s\n",$1,$2)}' database/virus.phy > database/virus.phy2
mv database/virus.phy2 database/virus.phy

echo -e "database/virus.phy\\nR\\n100\\nY\\n9" > database/par/seqboot.par
echo -e "database/seqboot.out\\nT\\n2.3628\\nM\\nD\\n100\\n2\\nY\\n" > database/par/dnadist.par
echo -e "database/dnadist.out\\nM\\n100\\n9\\nY" > database/par/neighbor.par
echo -e "database/nei.tree\\nY" > database/par/consense.par

phylip seqboot < database/par/seqboot.par 
mv outfile database/seqboot.out 
phylip dnadist < database/par/dnadist.par 
mv outfile database/dnadist.out 
phylip neighbor < database/par/neighbor.par 
mv outfile database/nei.out
mv outtree database/nei.tree 
phylip consense < database/par/consense.par 
mv outfile database/cons.out
mv outtree database/constree 


##-- 4. make ORF database

for id in $(ls database/virus_seq/*fasta); do { sample=`echo ${id}|sed 's/.*\///g'|sed 's/.fasta//g'`; /srv/shiny-server/software/ORFfinder -in $id -ml $orf_len -s $start_codon -g $genetic_code -n $ignore_nest -out ORF/${sample}.orfs; }& done ;wait
sed -i 's/lcl|//g' ORF/*.orfs
sed -i 's/ unnamed protein product.*//g' ORF/*.orfs

## files for orf viewer
for id in $(ls ORF/*orfs); do { Rscript /srv/shiny-server/script/orf_process.R $id; }& done ;wait

cat ORF/*.orfs > database/orf_db.fasta
makeblastdb -in database/orf_db.fasta -parse_seqids -dbtype prot

##-- 5. blastp
for id in $(ls ORF/*orfs); do { sample=`echo ${id}|sed 's/.*\///g'|sed 's/.orfs//g'`; blastp -query $id -db database/orf_db.fasta -outfmt 10 -evalue 0.05 > blastp/${sample}.blast; }& done ; wait
for id in $(ls blastp/*blast); do { Rscript /srv/shiny-server/script/blast_process.R $id ; }& done ; wait
for id in $(ls blastp/*fa); do { mafft $id > ${id%%.fa}.msa.fa; }& done; wait
#for id in $(ls blastp/*msa.fa); do Rscript /srv/shiny-server/script/msa_html.R $id; done


