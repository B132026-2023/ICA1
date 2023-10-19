#!usr/bin/bash

#CREATE DIRECTORY FOR ALL GENERATED RESULTS
rm -rf B132026_results
mkdir -p B132026_results

#PATH TO PAIRED END FASTQ SEQUENCES (CHANGE AS REQUIRED)
SEQPATH=/localdisk/data/BPSM/ICA1/fastq/

#PATH TO REFERENCE GENOME (CHANGE AS REQUIRED)
REFGEN=/localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz

#PATH TO THE GENE LOCI BED FILE (CHANGE AS REQUIRED)
REFBED=/localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed

#FASTQC EXTRACTION AND ANALYSIS
#fastqc ${SEQPATH}/*.fq.gz -o ./B132026_results/ -t 6

mkdir -p fastqc_results

fastqc --extract ${SEQPATH}/*.fq.gz -o ./fastqc_results -t 32

rm ./fastqc_results/*.html
rm ./fastqc_results/*.zip
rm fastq_summary.txt

touch fastq_summary.txt
echo -e "Sequence\tPass Count\tWarn Count\tFail Count" >> fastq_summary.txt
for dir in ./fastqc_results/*/; do

    fail_count=$(awk -F'\t' '$1 == "FAIL" { count++ } END { print (count > 0) ? count : 0 }' "$dir/summary.txt")
    warn_count=$(awk -F'\t' '$1 == "WARN" { count++ } END { print (count > 0) ? count : 0 }' "$dir/summary.txt")
    pass_count=$(awk -F'\t' '$1 == "PASS" { count++ } END { print (count > 0) ? count : 0 }' "$dir/summary.txt")

    file=$(echo $dir | awk -F'/' '{print $3}' | awk -F'_' '{print $1 "_" $2}')
    echo -e "$file\t$pass_count\t$warn_count\t$fail_count" >> fastq_summary.txt

done

rm -rf fastqc_results

#ALIGNMENT OF FASTQ READS WITH THE REFERENCE GENOME
bowtie2-build ${REFGEN} ./B132026_results/reference --threads 64

for sequence in ${SEQPATH}*_1.fq.gz; do

    sample_name=$(echo $sequence|cut -d"_" -f 1)

    bowtie2 --local \
            -x ./B132026_results/reference \
            -1 ${sample_name}_1.fq.gz -2 ${sample_name}_2.fq.gz \
            -S ./B132026_results/$(basename ${sample_name}).sam \
            --threads 64

    samtools sort ./B132026_results/$(basename ${sample_name}).sam \
            -o ./B132026_results/$(basename ${sample_name}).bam \
            -@ 64

    samtools index ./B132026_results/$(basename ${sample_name}).bam -@ 64

    rm ./B132026_results/*.sam

#GENERATE COUNT PER GENE FOR EACH ALIGNMENT
    bedtools coverage \
        -a ${REFBED} \
        -b ./B132026_results/$(basename ${sample_name}).bam \
	-counts > ./B132026_results/$(basename ${sample_name})_counts.txt

    rm ./B132026_results/*.bam*

done

rm ./B132026_results/reference*
