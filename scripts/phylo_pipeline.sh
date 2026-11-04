home=$
thread=$

script=$
References=$
refname=D39
ref=SpnSero02D39.fasta
ref1=SpnSero02D39.fna
Index1=`echo $References`"/SpnSero02D39/Index1"

# Build Index (speed up the alignment)
echo "Make OneLine Fasta for faidx indexing" # Used here
python $script\/MakeOneLineFasta.py $References\/SpnSero02D39/$ref $References\/SpnSero02D39/$ref1
mkdir $References\/SpnSero02D39/Index1/
bowtie2-build $References\/SpnSero02D39/$ref1 $References\/SpnSero02D39/Index1/index
samtools faidx $References\/SpnSero02D39/$ref1

# Define phage, repeat and variation regions (from Genbank and PHASTER)
python $script\/upload/extract_position_from_genbank.py --genbank $References\/SpnSero02D39/SpnSero02D39.gbk > $References\/SpnSero02D39/SpnSero02D39.phage_repeat_variation_all
cut -f 1 $References\/SpnSero02D39/SpnSero02D39.phage_repeat_variation_all | tr ',' '\t'  > $References\/SpnSero02D39/SpnSero02D39.phage_repeat_variation_from_genbank.txt	

## from PHASTER
cut -f 1 -d ' ' $References\/SpnSero02D39/SpnSero02D39.phage_regions_from_phaster_raw.txt | grep . | grep -v "#" | grep -v "CDS" | grep -v ">"| grep -v "-" | sed -e 's/#//g' | sed -e 's/complement(//g' | sed -e 's/)//g' | tr '..' '#' | sed -e 's/##/#/g'| tr '#' '\t' > $References\/SpnSero02D39/SpnSero02D39.phage_regions_from_phaster.txt

# need to be manually reviewed
## merge them all together
cat $References\/SpnSero02D39/SpnSero02D39.phage_regions_from_phaster.txt $References\/SpnSero02D39/SpnSero02D39.phage_repeat_variation_from_genbank.txt > $References\/SpnSero02D39/SpnSero02D39.phage_repeat_variation_final.txt

path=$


for dir in $path\/*
do
	name=$(basename $dir)
	echo $name
	mkdir -p $dir\/2_Trimmed
	fastp --in1 $dir\/1_RawData/$name*\_1.fastq.gz --in2 $dir\/1_RawData/$name*\_2.fastq.gz --out1 $dir\/2_Trimmed/$name\_R1_TrmP.fastq.gz --out2 $dir\/2_Trimmed/$name\_R2_TrmP.fastq.gz --average_qual 20 --length_required 35 --trim_poly_g --thread 16 --trim_poly_x
done


mkdir -p $home\/Z4_VCFphageD/

for dir in $path\/$head*
do
	name=`basename $dir`
	#name=`basename $dir| cut -f 1 -d "_"`
	echo $name
	mkdir -p $dir\/4_Mapping_$refname/
	time(

	echo "# Mapping the reads"
	bowtie2 -x $Index1\/index -1 $dir\/2_Trimmed/$name\_R1_TrmP.fastq.gz -2 $dir\/2_Trimmed/$name\_R2_TrmP.fastq.gz -S $dir\/4_Mapping_$refname/$name\.sam -p $thread
	#minimap2 -ax asm5 $References\/SpnSero02D39/SpnSero02D39.fasta  $dir\/1_RawData/$name\.fasta | samtools view -ubS --threads $thread > $dir\/4_Mapping_$refname/$name\.bam
	
	echo "# Convert .sam to .bam"
	samtools view -ubS --threads $thread $dir\/4_Mapping_$refname/$name\.sam > $dir\/4_Mapping_$refname/$name\.bam
	
	echo "# Sort the .bam"
	samtools sort --threads $thread $dir\/4_Mapping_$refname/$name\.bam -o $dir\/4_Mapping_$refname/$name\.sorted.bam
	
	echo "# Index sorted BAM alignments"
	samtools index -@ $thread $dir\/4_Mapping_$refname/$name\.sorted.bam
	
	echo "# Call SNPs and Indels"
	samtools mpileup -d 1000 -t DP -Q 20 -t SP -ugBf $References\/SpnSero02D39/$ref $dir\/4_Mapping_$refname/$name\.sorted.bam | bcftools call -cv --threads $thread | bcftools view -v snps --threads $thread > $dir\/4_Mapping_$refname/$name\.raw.vcf
	
	echo "# Filter for SNP Quality for_phylo (homogenous and phred 20)"
	bcftools view --threads $thread -g hom $dir\/4_Mapping_$refname/$name\.raw.vcf > $dir\/4_Mapping_$refname/$name\.hom.vcf
	python $script\/Filter_SNP_quality_for_SPN_phylo.py --qual 20 --vcf $dir\/4_Mapping_$refname/$name\.hom.vcf --output $dir\/4_Mapping_$refname/$name\.phrd.vcf
	
	echo "# Filter for biases"
	vcfutils.pl varFilter -Q 20 -1 0.001 -3 0.001 -4 0.001  $dir\/4_Mapping_$refname/$name\.phrd.vcf > $dir\/4_Mapping_$refname/$name\.biasF.vcf
	
	echo "# Filter out SNPs in phage_regions"
	python $script\/Filter_SNP_phage_repeat_regions.py --position $References\/SpnSero02D39/SpnSero02D39.phage_repeat_variation_final.txt --vcf $dir\/4_Mapping_$refname/$name\.biasF.vcf --output $dir\/4_Mapping_$refname/$name\.phageD.vcf
	
	echo "# Delete .bam file & Gzip .sam & sorted.bam file"
	pigz -9 -p $thread -f $dir\/4_Mapping_$refname/$name\.sam
	rm $dir\/4_Mapping_$refname/$name.bam

	echo '# Generate consensus sequence.'
	cp $dir\/4_Mapping_$refname/$name\.phageD.vcf $home\/Z4_VCFphageD/$name\.phageD.vcf
	cp $home\/Z4_VCFphageD/$name\.phageD.vcf $home\/Z4_VCFphageD/$name\.phageD_1.vcf
	bgzip -f --threads $thread $home\/Z4_VCFphageD/$name\.phageD_1.vcf
	bcftools index --threads $thread $home\/Z4_VCFphageD/$name\.phageD_1.vcf.gz
	bcftools consensus -f $References\/SpnSero02D39/$ref1 $home\/Z4_VCFphageD/$name\.phageD_1.vcf.gz > $home\/Z4_VCFphageD/$name\.consensus.fasta
	sed -i "1 s/.*/>$name/" $home\/Z4_VCFphageD/$name\.consensus.fasta
	perl -pe '/^>/ ? print "\n" : chomp' $home\/Z4_VCFphageD/$name\.consensus.fasta > $home\/Z4_VCFphageD/$name\.consensus1.tmp
	grep . $home\/Z4_VCFphageD/$name\.consensus1.tmp > $home\/Z4_VCFphageD/$name\.consensus.fasta
	rm $home\/Z4_VCFphageD/$name\.consensus1.tmp $home\/Z4_VCFphageD/$name\.phageD_1.vcf.gz*
	)
done

#***********************************************#
#					GUBBINS						#
#***********************************************#
## Gather all consensus sequences
echo "# Gather all consensus sequences."
mkdir -p $home\/Z5_RecombGubbins

thread=92

rm $home\/Z4_VCFphageD/GBS_0445.consensus.fasta

cat $home\/Z4_VCFphageD/*.consensus.fasta > $home\/Z5_RecombGubbins/Sero02_n170_consensus.fasta
grep -c ">" $home\/Z5_RecombGubbins/Sero02_n170_consensus.fasta && wc -l $home\/Z5_RecombGubbins/Sero02_n170_consensus.fasta
grep ">" $home\/Z5_RecombGubbins/Sero02_n170_consensus.fasta | sed -e 's/>//g' > $home\/Z5_RecombGubbins/Sero02_n170.id

# Run_Gubbins
echo "# Run Gubbins"
date +"%d-%b-%Y %r"
cd $home\/Z5_RecombGubbins/
time(run_gubbins.py --threads $thread --first-tree-builder raxml --converge-method recombination --model-fitter raxml --first-model-fitter raxml --model GTRGAMMA --verbose Sero02_n170_consensus.fasta)

date +"%d-%b-%Y %r"

cd && mkdir $home\/Z5_RecombGubbins/GubbinsFiles

mv $home\/Z5_RecombGubbins/Sero02_n170* $home\/Z5_RecombGubbins/GubbinsFiles/
cp $home\/Z5_RecombGubbins/GubbinsFiles/Sero02_n170*.gff $home\/Z5_RecombGubbins/
echo "# Zip the large files"
pigz -p 12 $home\/Z5_RecombGubbins/GubbinsFiles/*

#################################################################################
#		Remove Gubbins (Recombinant Sites) from VCF files		#
#################################################################################
echo "Make Parsed Location file."
wc -l $home\/Z5_RecombGubbins/Sero02_n170.recombination_predictions.gff
grep -v '#' $home\/Z5_RecombGubbins/Sero02_n170.recombination_predictions.gff | cut -f 4,5 > $home\/Z5_RecombGubbins/Sero02_n170.recombinant_predicted_regions.txt
wc -l $home\/Z5_RecombGubbins/Sero02_n170.recombinant_predicted_regions.txt

#### Final Filter for all isolates
mkdir $home\/Z6_VCFphylo
mkdir $home\/Z6_VCFphylopos

for filename in $home\/Z4_VCFphageD/*.phageD.vcf
do
	name=`basename $filename| cut -f 1 -d '.'`
	echo $name
	time(
	echo "# Filter SNP in recombinant (gubbins) sites"
	python $script\/Filter_SNP_phage_repeat_regions.py --position $home\/Z5_RecombGubbins/Sero02_n170.recombinant_predicted_regions.txt --vcf $filename --output $home\/Z6_VCFphylo/$name\.Phylo.vcf
	echo "# Extract position"
	grep  -v '#' $home\/Z6_VCFphylo/$name\.Phylo.vcf | cut -f 2 > $home\/Z6_VCFphylopos/$name\.Phylo.pos
	wc -l $home\/Z6_VCFphylopos/$name\.Phylo.pos
	)
done

# Check final number of SNPs (Remove OUTGROUP FILE FROM POSITON FOLDER BEFORE DOING IT).
## Outgroup: GBS_0445
rm $home\/Z6_VCFphylopos/GBS_0445.Phylo.pos

cat $home\/Z6_VCFphylopos/*.Phylo.pos | sort -ug > $home\/Sero02_n170_SNPall.pos
wc -l $home\/Sero02_n170_SNPall.pos
head -1 $home\/Sero02_n170_SNPall.pos
tail -1 $home\/Sero02_n170_SNPall.pos

#################################################################
#				Extract SNPS alleles & Run  RAxML				#
#################################################################
# Make SNPallele
mkdir $home\/Z7_SNPallele/

for filename in $home\/Z6_VCFphylo/*.Phylo.vcf
do
	name=`basename $filename| cut -f 1 -d '.'`
	echo $name
	time(
	echo '# Generate SNPalleles.'
	python $script\/extract_position_specific_SNP.py -r $References\/SpnSero02D39/$ref1 -s $home\/Sero02_n170_SNPall.pos -v $filename -i $name -o $home\/Z7_SNPallele/$name\.fasta
	)
done

## Gather all SNPallele sequences
echo "# Gather all global SNPallele sequences."
mkdir -p $home\/Z8_Phylogeny/

cat $home\/Z7_SNPallele/*.fasta > $home\/Z8_Phylogeny/Sero02_n170_SNPallele.fasta
sed -i -e 's/>/\n>/g' $home\/Z8_Phylogeny/Sero02_n170_SNPallele.fasta
grep -c ">" $home\/Z8_Phylogeny/Sero02_n170_SNPallele.fasta && wc -l $home\/Z8_Phylogeny/Sero02_n170_SNPallele.fasta
grep ">" $home\/Z8_Phylogeny/Sero02_n170_SNPallele.fasta | sed -e 's/>//g' > $home\/Z8_Phylogeny/Sero02_n170_SNPallele.idd


echo "Convert to phylip"
python $script\/fasta2phylip.py -i $home\/Z8_Phylogeny/Sero02_n170_SNPallele.fasta -o $home\/Z8_Phylogeny/Sero02_n170_SNPallele.phylip
pigz -p 12 $home\/Z8_Phylogeny/Sero02_n170_SNPallele.fasta

mkdir $home\/Z8_Phylogeny/raxml
cd $home\/Z8_Phylogeny/

# Compute the best topology ("RAxML_bestTree" file):
time(raxmlHPC-PTHREADS-AVX -p 1206 -s Sero02_n170_SNPallele.phylip -m GTRGAMMA -T $thread -o GBS_0445 -n Sero02_n170_SNPallele.bestTree)
# Compute 100 bootstraps to assess tree structure ( "RAxML_bootstrap" file):
time(raxmlHPC-PTHREADS-AVX -p 1206 -x 4221 -s Sero02_n170_SNPallele.phylip -N 100 -m GTRGAMMA -T $thread -o GBS_0445 -n Sero02_n170_SNPallele.bsSupport)
# Compute final tree: (with RAxML_bestTree & RAxML_bootstrap)
time(raxmlHPC-PTHREADS-AVX -p 1206 -m GTRGAMMA -f b -t RAxML_bestTree.Sero02_n170_SNPallele.bestTree -z RAxML_bootstrap.Sero02_n170_SNPallele.bsSupport -o GBS_0445 -n Sero02_n170_SNPallele.bestTree.bsSupport)

mv $home\/Z8_Phylogeny/*Sero02_n170_SNPallele* $home\/Z8_Phylogeny/raxml/
cp $home\/Z8_Phylogeny/raxml/RAxML_bipartitionsBranchLabels.* $home\/Z8_Phylogeny/
cp $home\/Z8_Phylogeny/RAxML_bipartitionsBranchLabels.Sero02_n170_SNPallele.bestTree.bsSupport $home\/Z8_Phylogeny/Sero02_n170_SNPallele.bestTree.bsSupport.corrected
pigz -p 12 $home\/Z8_Phylogeny/raxml/*



#################################################################
#																#
#						GPSC96 RAxML							#
#																#
#################################################################
Sero02_Home=$
home=$
thread=48

script=$
References=$
refname=D39
ref=SpnSero02D39.fasta
ref1=SpnSero02D39.fna

mkdir -p $home

## Copy all PhageD files from Sero02_Home to GPSC96 Home.
cp -r $Sero02_Home\/Z4_VCFphageD/ $home\/

## Delete all non-GPSC96 files
ls $home\/Z4_VCFphageD/*.phageD.vcf | grep -f $home\/Sero02_nonGPSC96_n4.txt | sed -e 's/^/rm /' | sed -e 's/phageD.vcf/*/g' > $home\/Sero02_nonGPSC96_n4_delete.sh
# manual edit for the script
bash $home\/Sero02_nonGPSC96_n4_delete.sh


#***********************************************#
#					GUBBINS						#
#***********************************************#
## Gather all consensus sequences
echo "# Gather all consensus sequences."
mkdir -p $home\/Z5_RecombGubbins

thread=48

rm $home\/Z4_VCFphageD/GBS_0445.consensus.fasta

cat $home\/Z4_VCFphageD/*.consensus.fasta > $home\/Z5_RecombGubbins/GPSC96_n166_consensus.fasta
grep -c ">" $home\/Z5_RecombGubbins/GPSC96_n166_consensus.fasta && wc -l $home\/Z5_RecombGubbins/GPSC96_n166_consensus.fasta
grep ">" $home\/Z5_RecombGubbins/GPSC96_n166_consensus.fasta | sed -e 's/>//g' > $home\/Z5_RecombGubbins/GPSC96_n166.id

# Run_Gubbins
echo "# Run Gubbins"
date +"%d-%b-%Y %r"
cd $home\/Z5_RecombGubbins/
time(run_gubbins.py --threads $thread --first-tree-builder raxml --converge-method recombination --model-fitter raxml --first-model-fitter raxml --model GTRGAMMA --verbose --prefix GPSC96_n166 GPSC96_n166_consensus.fasta)
date +"%d-%b-%Y %r"

cd && mkdir $home\/Z5_RecombGubbins/GubbinsFiles

mv $home\/Z5_RecombGubbins/GPSC96_n166* $home\/Z5_RecombGubbins/GubbinsFiles/
cp $home\/Z5_RecombGubbins/GubbinsFiles/GPSC96_n166*.gff $home\/Z5_RecombGubbins/
echo "# Zip the large files"
pigz -p 12 $home\/Z5_RecombGubbins/GubbinsFiles/*

#################################################################################
#		Remove Gubbins (Recombinant Sites) from VCF files		#
#################################################################################
echo "Make Parsed Location file."
wc -l $home\/Z5_RecombGubbins/GPSC96_n166.recombination_predictions.gff
grep -v '#' $home\/Z5_RecombGubbins/GPSC96_n166.recombination_predictions.gff | cut -f 4,5 > $home\/Z5_RecombGubbins/GPSC96_n166.recombinant_predicted_regions.txt
wc -l $home\/Z5_RecombGubbins/GPSC96_n166.recombinant_predicted_regions.txt

#### Final Filter for all isolates
mkdir $home\/Z6_VCFphylo
mkdir $home\/Z6_VCFphylopos

for filename in $home\/Z4_VCFphageD/*.phageD.vcf
do
	name=`basename $filename| cut -f 1 -d '.'`
	echo $name
	time(
	echo "# Filter SNP in recombinant (gubbins) sites"
	python $script\/Filter_SNP_phage_repeat_regions.py --position $home\/Z5_RecombGubbins/GPSC96_n166.recombinant_predicted_regions.txt --vcf $filename --output $home\/Z6_VCFphylo/$name\.Phylo.vcf
	echo "# Extract position"
	grep  -v '#' $home\/Z6_VCFphylo/$name\.Phylo.vcf | cut -f 2 > $home\/Z6_VCFphylopos/$name\.Phylo.pos
	wc -l $home\/Z6_VCFphylopos/$name\.Phylo.pos
	)
done

# Check final number of SNPs (Remove OUTGROUP FILE FROM POSITON FOLDER BEFORE DOING IT).
## Outgroup: GBS_0445, Serotype 35B
rm $home\/Z6_VCFphylopos/GBS_0445.Phylo.pos

cat $home\/Z6_VCFphylopos/*.Phylo.pos | sort -ug > $home\/GPSC96_n166_SNPall.pos
wc -l $home\/GPSC96_n166_SNPall.pos
head -1 $home\/GPSC96_n166_SNPall.pos
tail -1 $home\/GPSC96_n166_SNPall.pos

#################################################################
#				Extract SNPS alleles & Run  RAxML				#
#################################################################

# Make SNPallele
mkdir $home\/Z7_SNPallele/

for filename in $home\/Z6_VCFphylo/*.Phylo.vcf
do
	name=`basename $filename| cut -f 1 -d '.'`
	echo $name
	time(
	echo '# Generate SNPalleles.'
	python $script\/extract_position_specific_SNP.py -r $References\/SpnSero02D39/$ref1 -s $home\/GPSC96_n166_SNPall.pos -v $filename -i $name -o $home\/Z7_SNPallele/$name\.fasta
	)
done

## Gather all SNPallele sequences
echo "# Gather all global SNPallele sequences."
mkdir -p $home\/Z8_Phylogeny/

cat $home\/Z7_SNPallele/*.fasta > $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.fasta
sed -i -e 's/>/\n>/g' $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.fasta
grep -c ">" $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.fasta && wc -l $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.fasta
grep ">" $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.fasta | sed -e 's/>//g' > $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.idd

echo "Convert to phylip"
python $script\/fasta2phylip.py -i $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.fasta -o $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.phylip
pigz -p 12 $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.fasta

mkdir $home\/Z8_Phylogeny/raxml
cd $home\/Z8_Phylogeny/

# Compute the best topology ("RAxML_bestTree" file):
time(raxmlHPC-PTHREADS-AVX -p 1206 -s GPSC96_n166_SNPallele.phylip -m GTRGAMMA -T $thread -o GBS_0445 -n GPSC96_n166_SNPallele.bestTree)
# Compute 100 bootstraps to assess tree structure ( "RAxML_bootstrap" file):
time(raxmlHPC-PTHREADS-AVX -p 1206 -x 4221 -s GPSC96_n166_SNPallele.phylip -N 100 -m GTRGAMMA -T $thread -o GBS_0445 -n GPSC96_n166_SNPallele.bsSupport)
# Compute final tree: (with RAxML_bestTree & RAxML_bootstrap)
time(raxmlHPC-PTHREADS-AVX -p 1206 -m GTRGAMMA -f b -t RAxML_bestTree.GPSC96_n166_SNPallele.bestTree -z RAxML_bootstrap.GPSC96_n166_SNPallele.bsSupport -o GBS_0445 -n GPSC96_n166_SNPallele.bestTree.bsSupport)

mv $home\/Z8_Phylogeny/*GPSC96_n166_SNPallele* $home\/Z8_Phylogeny/raxml/
cp $home\/Z8_Phylogeny/raxml/RAxML_bipartitionsBranchLabels.* $home\/Z8_Phylogeny/
cp $home\/Z8_Phylogeny/RAxML_bipartitionsBranchLabels.GPSC96_n166_SNPallele.bestTree.bsSupport $home\/Z8_Phylogeny/GPSC96_n166_SNPallele.bestTree.bsSupport.corrected
pigz -p 12 $home\/Z8_Phylogeny/raxml/*