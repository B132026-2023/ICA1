
#PATH VARIABLES (CHANGE AS REQUIRED)
INFOFILE=/localdisk/data/BPSM/ICA1/fastq/Tco2.fqfiles
REFBED=/localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed

#EMPTY ARRAYS TO STORE PREFERENCE ARGUMENTS FOR GROUPWISE COMPARISONS
#preferences for group 1
types_inputs1=()
times_inputs1=()
treatments_inputs1=()

#preferences for group 2
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


#FUNCTION TO FILTER FOR GENOMIC SAMPLES BASED ON USER DEFINED GROUPWISE PREFERENCES 
function filter() {

    local types=(${!1})
    local times=(${!2})
    local treatments=(${!3})
    files=()


    while IFS=$'\t' read -r SampleName SampleType Replicate Time Treatment End1 End2; do
       if [[ ! ${types[@]} =~ ${SampleType} ]] &&
          [[ ! ${times[@]} =~ ${Time} ]] &&
          [[ ! ${treatments[@]} =~ ${Treatment} ]] &&
          [[ ! "SampleName" =~ ${SampleName}  ]]; then
               files+=(${SampleName}_counts.txt)
	fi
    done < ${INFOFILE}
}

#FUNCTION TO CALCUALTE THE COUNTS AND MEANS FOR EACH GENE FOR FILTERED GENOMIC FILES
function generate_means() {

    local files=${!1}
    local output_file=${2}

    mkdir -p 6th_col

    for file in ./B132026_results/*.txt; do

        if [[ ${files[@]} =~ $(basename ${file//-/}) ]]; then
            cut -f6 $file > ./6th_col/$(basename $file).6th_col
        fi

    done

    paste ./6th_col/*.6th_col > ${output_file}.counts
    rm -rf 6th_col

    awk 'BEGIN{FS="\t";}{sum=0; for(i=1; i<=NF; i++) sum+=$i; \
     mean=sum/NF; print mean}' ${output_file}.counts \
     > ${output_file}.means
}

#CALLING FILTER, AND COUNTS AND MEAN FUNCTION
#for group 1
filter types_inputs1[@] times_inputs1[@] treatments_inputs1[@]
generate_means files[@] "group_1"

#for group 2
filter types_inputs2[@] times_inputs2[@] treatments_inputs2[@]
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
      --output=fold_change.sorted.output

#REMOVES ALL TEMP FILES
rm fold_change.tmp
rm *.means
rm *.tmp
