#!usr/bin/bash


#PATH VARIABLES (CHANGE AS REQUIRED)
SEQPATH=/localdisk/data/BPSM/ICA1/fastq/ #path to fastq sequences
REFGEN=/localdisk/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz
#path to reference genome^
REFBED=/localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed #path to bed file
INFOFILE=/localdisk/data/BPSM/ICA1/fastq/Tco2.fqfiles #path to information file


#EMPTY ARRAYS TO STORE PREFERENCE ARGUMENTS FOR GROUPWISE COMPARISONS
#exclusion preferences for group 1
types_inputs1=()
times_inputs1=()
treatments_inputs1=()

#exclusion preferences for group 2
types_inputs2=()
times_inputs2=()
treatments_inputs2=()


#WHILE GETOPTS LOOP TO ENABLE GROUPWISE PREFERENCES VIA DEFINED FLAGS
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


#FUNCTION TO SUBSET GENOMIC SAMPLES BASED ON USER DEFINED GROUPWISE PREFERENCES
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


#FUNCTION TO CALCUALTE THE GROUPWISE COUNTS AND MEANS OF COUNTS FOR EACH GENE
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

    paste ./6th_col/*.6th_col > ${output_file}.counts #generate a tab deliminated file of counts for specified group
    rm -rf 6th_col

    #Calcualte groupwise count mean
    awk 'BEGIN{FS="\t";}{sum=0; for(i=1; i<=NF; i++) sum+=$i; \
     mean=sum/NF; print mean}' ${output_file}.counts \
     > ${output_file}.means
}


#FASTQC EXTRACTION AND ANALYSIS
#fastqc ${SEQPATH}/*.fq.gz -o ./B132026_results/ -t 6

#CREATE DIRECTORY FOR ALL GENERATED RESULTS
rm -rf B132026_results
mkdir -p B132026_results #reset and create tmp directory
mkdir -p fastqc_results
mkdir -p fastqc_summary
fastqc --extract ${SEQPATH}/*.fq.gz -o ./fastqc_results -t 32

rm ./fastqc_results/*.html #remove unwanted files
rm ./fastqc_results/*.zip
rm ./fastqc_summary.txt #reset file

touch fastqc_summary.txt
echo -e "Sequence\tPass Count\tWarn Count\tFail Count" >> fastqc_summary.txt # add headers
for dir in ./fastqc_results/*/; do

    fail_count=$(awk -F'\t' '$1 == "FAIL" { count++ } END { print (count > 0) ? count : 0 }' "$dir/summary.txt")
    warn_count=$(awk -F'\t' '$1 == "WARN" { count++ } END { print (count > 0) ? count : 0 }' "$dir/summary.txt")
    pass_count=$(awk -F'\t' '$1 == "PASS" { count++ } END { print (count > 0) ? count : 0 }' "$dir/summary.txt")

    file=$(echo $dir | awk -F'/' '{print $3}' | awk -F'_' '{print $1 "_" $2}')
    echo -e "$file\t$pass_count\t$warn_count\t$fail_count" >> fastqc_summary.txt #add fastqc counts

    mv ${dir}summary.txt ./fastqc_summary/${file}_summary.txt

done

mv ./fastqc_summary.txt ./fastqc_summary/
rm -rf fastqc_results


#ALIGNMENT OF FASTQ READS WITH THE REFERENCE GENOME
bowtie2-build ${REFGEN} ./B132026_results/reference --threads 64

#iterate through fastq sequences
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


#CALLING SUBSET, AND COUNTS AND MEAN FUNCTION
#for group 1
subset types_inputs1[@] times_inputs1[@] treatments_inputs1[@]
generate_means files[@] "group_1"

#for group 2
subset types_inputs2[@] times_inputs2[@] treatments_inputs2[@]
generate_means files[@] "group_2"

#PASTING GENERATING MEAN FILES WITH GENE DESCRIPTIONS SPECIFIED IN THE BED FILE FOR EACH GROUP
paste -d "\t" <(cut -f5 ${REFBED}) <(cut -f1 "group_1.means") > "group_1_means.output"
paste -d "\t" <(cut -f5 ${REFBED}) <(cut -f1 "group_2.means") > "group_2_means.output"
#cat "group_2.means" | cut -f1 "group_2.means"

#ADDING A SMALL CONSTANT VALUE (0.1) TO EACH MEAN VLAUE TO AVOID DIVIDING WITH 0 ERRORS
awk -v constant=0.1 '{print $1 + constant}' group_1.means > group_1.constant.means
awk -v constant=0.1 '{print $1 + constant}' group_2.means > group_2.constant.means
paste -d "\t" "group_1.constant.means" "group_2.constant.means" > means.tmp

#GENERATES A SORTED FOLD CHANGE FILE BASED ON THE MEAN FILES
awk '{printf "%.2f\n", $2/$1}' means.tmp > fold_change.tmp
paste -d "\t" <(cut -f5 ${REFBED}) <(cut -f1 "fold_change.tmp") | sort -t$'\t' -k2 -nr \
      --output=fold_change.sorted
#{ echo -e "Gene descriptions\tFold Change"; cat fold_change.sorted; } > fold_change.sorted.output
#rm fold_change.sorted

#REMOVES ALL TEMP FILES
rm *.means
rm *.tmp
