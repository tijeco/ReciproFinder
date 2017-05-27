
working_dir=$(date|sed 's/ /_/g'|sed 's/:/./g')
echo $working_dir
mkdir $working_dir

cd $working_dir
python ../makeseed.py
/media/BigRAID/JeffTemp/EvolvAGene4Package/EvolveAGene4-linux-x86-64 -f seed.fa -n 4 -o Phy -b 0.32
python3 ../splitTogene.py

for f in gene*

do
	cd $f
	/media/BigRAID/JeffTemp/EvolvAGene4Package/EvolveAGene4-linux-x86-64 -f $f.fa -n 4 -o Phy -b 0.32
	cd ..


done

for f in gene*/*pep_Unaligned.FASTA
do
	cat $f |sed "/^>/s/$/_$(dirname $f)/" >> Sim_genomes.fasta

done
makeblastdb -in Sim_genomes.fasta -out Sim_genomes.fasta.seq.db -dbtype prot
blastp -db Sim_genomes.fasta.seq.db -query Sim_genomes.fasta -outfmt 6 -out Sim_genomes.fasta.blastall -num_threads 13 -evalue 1E-5
cd ..
