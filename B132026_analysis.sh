#!usr/bin/bash


#PATH VARIABLES (CHANGE AS REQUIRED)
SEQPATH=/localdisk/data/BPSM/ICA1/fastq/ #path to fastq sequences
REFGEN=/localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz
#path to reference genome^
REFBED=/localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed #path to bed file
INFOFILE=/localdisk/data/BPSM/ICA1/fastq/Tco2.fqfiles #path to information file


#empty arrays to store preferences for group definitions
#exclusion preferences for group 1
types_inputs1=()
times_inputs1=()
treatments_inputs1=()

#exclusion preferences for group 2
types_inputs2=()
times_inputs2=()
treatments_inputs2=()


#while getopts loop to enable flags for groupwise comaprisons of gene expression levels
while getopts :a:b:c:x:y:z: opt; do
    case $opt in
        a)
            IFS=',' read -r -a types_inputs1 <<< "$OPTARG"
            ;;
        b)
            IFS=',' read -r -a times_inputs1 <<< "$OPTARG"
            ;;
        c)
            IFS=',' read -r -a treatments_inputs1 <<< "$OPTARG"
            ;;
        x)
            IFS=',' read -r -a types_inputs2 <<< "$OPTARG"
            ;;
        y)
            IFS=',' read -r -a times_inputs2 <<< "$OPTARG"
            ;;
        z)
            IFS=',' read -r -a treatments_inputs2 <<< "$OPTARG"
            ;;
       \?)
            echo "Invalid option -$OPTARG"
            exit 1
            ;;
     esac
done


#function to subset fastq samples based in user defined preferences
function subset() {

    local types=(${!1}) #defining parameter arrays
    local times=(${!2})
    local treatments=(${!3})
    files=()

    #while loop to filter for files of interest
    while IFS=$'\t' read -r SampleName SampleType Replicate Time Treatment End1 End2; do
       if [[ ! ${types[@]} =~ ${SampleType} ]] &&
          [[ ! ${times[@]} =~ ${Time} ]] &&
          [[ ! ${treatments[@]} =~ ${Treatment} ]] &&
          [[ ! "SampleName" =~ ${SampleName}  ]]; then
               files+=(${SampleName}_counts.txt)
        fi
    done < ${INFOFILE}
}


#function to calculate the counts and means of counts for each gene for a group
function generate_means() {

    local files=${!1} #input array
    local output_file=${2} #input string

    mkdir -p 6th_col #temp file

    for file in ./B132026_results/*.txt; do

	#determines counts of specified group
        if [[ ${files[@]} =~ $(basename ${file//-/}) ]]; then
            cut -f6 $file > ./6th_col/$(basename $file).6th_col
        fi

    done

    #generate a tab deliminated file of counts for specified group
    paste ./6th_col/*.6th_col > ${output_file}.counts
    rm -rf 6th_col

    #Calcualte groupwise count mean
    awk 'BEGIN{FS="\t";}{sum=0; for(i=1; i<=NF; i++) sum+=$i; \
     mean=sum/NF; print mean}' ${output_file}.counts \
     > ${output_file}.means
}


#creates directory to store output files
rm -rf ./B132026_results
mkdir -p ./B132026_results
mkdir -p ./fastqc_results
mkdir -p ./fastqc_summary

#generates fastqc reports and extracts data
fastqc --extract ${SEQPATH}/*.fq.gz -o ./fastqc_results -t 32

#remove unwanted files
rm ./fastqc_results/*.html
rm ./fastqc_results/*.zip

#creates output file and adds header
echo -e "Sequence\tPass Count\tWarn Count\tFail Count" > fastqc_summary.txt

#iterates through all fastqc reports and extracts quility control data
for dir in ./fastqc_results/*/; do

    #counts number of declared PASS, WARN, and FAIL outputs
    fail_count=$(awk -F'\t' '$1 == "FAIL" { count++ } END { print (count > 0) ? count : 0 }' "$dir/summary.txt")
    warn_count=$(awk -F'\t' '$1 == "WARN" { count++ } END { print (count > 0) ? count : 0 }' "$dir/summary.txt")
    pass_count=$(awk -F'\t' '$1 == "PASS" { count++ } END { print (count > 0) ? count : 0 }' "$dir/summary.txt")

    file=$(echo $dir | awk -F'/' '{print $3}' | awk -F'_' '{print $1 "_" $2}')
    echo -e "$file\t$pass_count\t$warn_count\t$fail_count" >> fastqc_summary.txt #add fastqc counts

    mv ${dir}summary.txt ./fastqc_summary/${file}_summary.txt

done

#moves output file to fastqc_summary directory
mv ./fastqc_summary.txt ./fastqc_summary/


#removes tmp directory
rm -rf fastqc_results


#builds an index of the reference genome for subsequent alignment
bowtie2-build ${REFGEN} ./B132026_results/reference --threads 64

#adds header line
echo -e "Sample Name\tAlignment Rate" > alignment_log.txt.output

#iterate through fastq sequences
for sequence in ${SEQPATH}*_1.fq.gz; do

    sample_name=$(echo $sequence|cut -d"_" -f 1)

    # local alignment of paired-end reads to the reference genome
    (bowtie2 --local \
            -x ./B132026_results/reference \
            -1 ${sample_name}_1.fq.gz -2 ${sample_name}_2.fq.gz \
            -S ./B132026_results/$(basename ${sample_name}).sam \
            --threads 64) 2> alignment_log.tmp #recording the standard output to a tmp file

    #extracting the % alignment for each sample and recording in tab-delimited format
    seq_name=$(echo ${sample_name} |cut -d"/" -f 7)
    alignment_rate=$(grep "overall alignment rate" alignment_log.tmp | awk '{print $1}')
    echo -e "$seq_name\t$alignment_rate" >> alignment_log.txt.output
    echo -e "$seq_name\t$alignment_rate"

    #sorting the sam file and outputting a bam file
    samtools sort ./B132026_results/$(basename ${sample_name}).sam \
            -o ./B132026_results/$(basename ${sample_name}).bam \
            -@ 64

    #indexing the bam file
    samtools index ./B132026_results/$(basename ${sample_name}).bam -@ 64

    #removes sam files
    rm ./B132026_results/*.sam


    #generates count per gene for each alignment
    bedtools coverage \
        -a ${REFBED} \
        -b ./B132026_results/$(basename ${sample_name}).bam \
	-counts > ./B132026_results/$(basename ${sample_name})_counts.txt

    #removes bam and bam.bai files
    rm ./B132026_results/*.bam*

done

#removes the reference genome's index files
rm ./B132026_results/reference*


#calls subset, and counts and mean functions
#for group 1
subset types_inputs1[@] times_inputs1[@] treatments_inputs1[@]
generate_means files[@] "group_1"


#for group 2
subset types_inputs2[@] times_inputs2[@] treatments_inputs2[@]
generate_means files[@] "group_2"


#removes all alignment count files
rm ./B132026_results/*_counts.txt


#generates tab delimited files detailing the mean number of counts per gene with gene descriptions 
paste -d "\t" <(cut -f5 ${REFBED}) <(cut -f1 "group_1.means") > "group_1_means.output"
paste -d "\t" <(cut -f5 ${REFBED}) <(cut -f1 "group_2.means") > "group_2_means.output"
#cat "group_2.means" | cut -f1 "group_2.means"


#adds a small constant value to each mean value
awk -v constant=0.1 '{print $1 + constant}' group_1.means > group_1.constant.tmp
awk -v constant=0.1 '{print $1 + constant}' group_2.means > group_2.constant.tmp
paste -d "\t" "group_1.constant.tmp" "group_2.constant.tmp" > means.tmp


#generates a sorted fold change tab delimited file
awk '{printf "%.2f\n", $2/$1}' means.tmp > fold_change.tmp
paste -d "\t" <(cut -f5 ${REFBED}) <(cut -f1 "fold_change.tmp") | sort -t$'\t' -k2 -nr \
      --output=fold_change.sorted.output


#moves all output files to the B132026_results directory
mv *.output ./B132026_results/
mv ./fastqc_summary ./B132026_results/
mv ./group_1.counts ./B132026_results/group_1.counts.output
mv ./group_2.counts ./B132026_results/group_2.counts.output

#removes all temporary files
rm *.tmp
rm *.means
