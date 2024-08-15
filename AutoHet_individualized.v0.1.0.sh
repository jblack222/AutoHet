#!/bin/bash

#### Prior to running this code, you should have individual sample.bam files from tsv2bam (if de novo) or samtools/similar (if reference-based)

#### THIS CODE SHOULD BE RUN FROM THE DIRECTORY CONTAINING YOUR sample.bam FILES
# samtools merge is annoying to run from a location outside of your samples

# >>>> List all your sample.bam files here, seperated by a tab with the suffix removed

samples="sample_1	sample_2	sample_3"

# >>>> SAMPLE LIST: sample_list.txt is a required list of all sample.bam files with 1 column, and each row should contain one file name
# eg
# sample_1.bam
# sample_2.bam
# sample_3.bam
# etc
# sample_list.txt should be placed in the bam_dir

sample_list="example.sample_list.txt"

# >>>> POPMAP: popmap.txt is a required list of all sample names and population assingments with 2 tab-seperated columns
# Each row should contain one file sample name and one group ID
# eg
# sample_1	Melbourne
# sample_2	Sydney
# sample_3	Darwin
# etc

popmap="example.popmap.txt"

# >>>> Set your directories

bam_dir="~/path/to/bams"
ref_dir="~/path/to/ref"
bcf_dir="~/path/to/bcfs"
vcf_dir="~/path/to/vcfs"
results_dir="~/path/to/results"
popmap_dir="~/path/to/popmap"

echo "Making a pseudo-reference file with individual files"
samtools merge $ref_dir/all_samples.merged.bam -b $bam_dir/$sample_list
samtools consensus -f fasta -o $ref_dir/all_samples.merged.fasta $ref_dir/all_samples.merged.bam

#### Output individual (observed) heterozygosity.
#### This produces VCF files for each individual, then filters these based on missing data, star alleles, and minimum and maximum depth, then calls each confidently genotyped site as HomRef, HomAlt, or Het. 
#### This uses a default min depth of 15 and a max of 71. See Nielsen R et al 2011 Nat Rev Genet, and Li H 2014 Bioinformatics for details on setting these cutoffs. 
#### The second half of this code is built from user Kevin Blighe's (https://www.biostars.org/u/41557/) comment at https://www.biostars.org/p/298361/

for K in $samples;
do
	echo "Making BCF and VCF of sample:"
	echo ${K}
	bcftools mpileup -a FORMAT/AD -a FORMAT/DP -I -O b -o $bcf_dir/${K}.calls.bcf.gz -f $ref_dir.merged.fasta $bam_dir/${K}.bam
	bcftools call -a GQ -m -O z -o $vcf_dir/${K}.calls.vcf.gz $bcf_dir/${K}.calls.bcf.gz

	echo "Filtering and checking heterozygosity for sample:"
	echo ${K}

	bcftools view -s ${K} $vcf_dir/${K}.calls.vcf.gz | \
		bcftools filter -e 'FORMAT/DP < 15' |  \
		bcftools filter -e 'FORMAT/DP > 71' |  \
		bcftools filter -e 'ALT="*"' | \
		bcftools filter -e 'GT="mis"'| \
		bcftools filter -e 'QUAL < 25' | \
		bcftools filter -e 'FORMAT/AD[*:0] < 3' | \
		bcftools filter -e 'FORMAT/AD[*:1] < 3' | \
	bcftools norm -a --atom-overlaps . -O z -o $vcf_dir/${K}.vcf.gz \

	paste <(bcftools view $vcf_dir/${K}.vcf.gz | \
		awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
		!/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
		\
	<(bcftools query -f '[\t%SAMPLE=%GT]\n' $vcf_dir/${K}.vcf.gz | \
		awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
		\
	<(bcftools query -f '[\t%SAMPLE=%GT]\n' $vcf_dir/${K}.vcf.gz | \
		awk 'BEGIN {print "nHomAlt"} {print gsub(/1\|1|1\/1/, "")}') \
		\
	<(bcftools query -f '[\t%SAMPLE=%GT]\n' $vcf_dir/${K}.vcf.gz | \
		awk 'BEGIN {print "nHomRef"} {print gsub(/0\|0|0\/0/, "")}') \
		\
	<(bcftools view $vcf_dir/${K}.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
		!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
		\
	<(bcftools view $vcf_dir/${K}.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
		!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\|1|1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
		\
	<(bcftools view $vcf_dir/${K}.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
		!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|0|0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
		\
	| sed 's/,\t/\t/g' | sed 's/,$//g' > $vcf_dir/${K}_het.txt \

	paste <( awk '{ sum += $6 } END { print sum/NR }' $vcf_dir/${K}_het.txt) > $vcf_dir/${K}_mean_het.txt  \

	paste <( wc -l $vcf_dir/${K}_het.txt) >  $vcf_dir/${K}_sites.txt  \

done

rm $results_dir/*.txt

mv $vcf_dir/*.txt $results_dir


cat $results_dir/*_mean_het.txt > $results_dir/all_samples.mean_het.no_labels.txt
paste $popmap_dir/$popmap $results_dir/all_samples.mean_het.no_labels.txt >  $results_dir/all_samples.mean_het.txt

cat $results_dir/*_sites.txt > $results_dir/all_samples.sites.no_labels.txt
paste  $results_dir/all_samples.mean_het.txt $results_dir/all_samples.sites.no_labels.txt > $results_dir/all_samples.mean_het.sites.txt
