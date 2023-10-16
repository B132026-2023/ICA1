#!usr/bin/bash

#CREATE DIRECTORY FOR ALL GENERATED RESULTS
mkdir -p B132026_results

#PATH TO PAIRED END FASTQ SEQUENCES (CHANGE AS REQUIRED)
SEQPATH=/localdisk/data/BPSM/ICA1/fastq/

#PATH TO REFERENCE GENOME (CHANGE AS REQUIRED)
REFGEN=/localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz

#GENERATE FASTQC REPORTS
#fastqc ${SEQPATH}/*.fq.gz -o ./B132026_results/ -t 6

#ALIGNMENT OF FASTQ READS WITH THE REFERENCE GENOME 
bowtie2-build ${REFGEN} ./B132026_results/reference --threads 64

for sequence in ${SEQPATH}*_1.fq.gz; do

    sample_name=$(echo $sequence|cut -d"_" -f 1)

    bowtie2 -x ./B132026_results/reference \
            -1 ${sample_name}_1.fq.gz -2 ${sample_name}_2.fq.gz \
            -S ./B132026_results/$(basename ${sample_name}).sam \
            --threads 64

    samtools view -bS ./B132026_results/$(basename ${sample_name}).sam > \
            ./B132026_results/$(basename ${sample_name}).bam \
            -@ 64

    rm *.sam
done
