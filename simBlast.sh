
working_dir=$(date|sed 's/ /_/g')
echo $working_dir
mkdir $working_dir

cd $working_dir
python ../makeseed.py
/media/BigRAID/JeffTemp/EvolvAGene4Package/EvolveAGene4-linux-x86-64 -f seed.fa -n 4 -o Phy -b 0.32
python3 ../splitTogene.py
for f in gene*fa
do

mkdir $f

donewnew
cd ..
