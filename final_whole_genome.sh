#https://www.ncbi.nlm.nih.gov/sra/?term=SRX146853

##Data Doownloading

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501130/SRR501130_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501130/SRR501130_2.fastq.gz

#STEP1: Extract the data
gunzip *.gz

#Step2: quality check uisng fastqc

sudo apt-get update
sudo apt-get install fastqc

fastqc *.fastq
#There are three main parameters to decide the quality
#1. Per base sequence quality
#2. Overrepresented seq
#3. Adapter content

#create assebly using spades

#install

wget -c https://github.com/ablab/spades/archive/refs/tags/v4.2.0.tar.gz

tar -xzf v4.2.0.tar.gz
cd spades-4.2.0/
sudo apt-get install g++
sudo apt-get install zlib1g-dev
sudo apt-get install cmake
sudo apt install -y build-essential cmake libbz2-dev libz-dev libcurl4-openssl-dev libssl-dev

./spades_compile.sh
ls
cd bin
./spades.py


##spades cmds

/mnt/f/Dr_OmicsLab/NGS/Whole-genome-Assembly/spades-4.2.0/bin/spades.py --isolate -m 250 -1 SRR501130_1.fastq -2 SRR501130_2.fastq -o spades_result

 
###Genome annotation using prokka

#installation

sudo apt update
sudo apt install git bioperl perl-doc hmmer barrnap aragorn ncbi-blast+ parallel prodigal wget unzip -y
sudo apt install cpanminus -y
sudo cpanm Bio::SearchIO::hmmer3

git clone https://github.com/tseemann/prokka.git
cd prokka
./bin/prokka --setupdb
cd bin
./prokka

#annotation by prokka

./prokka contigs.fasta --outdir prokka_result

##BLAST against protein database

sudo apt-get install ncbi-blast+
mkdir blast
cd blast

wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz


#a) creation of ref db:

makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprodb   

#b) mapping of query with reference:

blastp -query longest_orfs.pep -db uniprodb -num_alignments 1 -num_threads 4 -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen" -out result.txt


#%Query_covg= (Exact_match/qlen)*100

#Exact_match= length-mismatch-gaps







