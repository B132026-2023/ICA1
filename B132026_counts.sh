
INFOFILE=/localdisk/data/BPSM/ICA1/fastq/Tco2.fqfiles

types_inputs1=()
times_inputs1=()
treatments_inputs1=()

types_inputs2=()
times_inputs2=()
treatments_inputs2=()

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

function filter() {

    local types=(${!1})
    local times=(${!2})
    local treatments=(${!3})
    files=()

	#echo ${types[@]}
	#echo ${times[@]}
	#echo ${treatments[@]}

    while IFS=$'\t' read -r SampleName SampleType Replicate Time Treatment End1 End2; do
       if [[ ! ${types[@]} =~ ${SampleType} ]] &&
          [[ ! ${times[@]} =~ ${Time} ]] &&
          [[ ! ${treatments[@]} =~ ${Treatment} ]] &&
          [[ ! "SampleName" =~ ${SampleName}  ]]; then
               files+=(${SampleName}_counts.txt)
	      # echo $SampleName $SampleType $Time $Treatment
	fi

    done < ${INFOFILE}
}




function alignment_count() {

    local files=${!1}
    local output_name=${2}

    mkdir -p 6th_col

    for file in ./B132026_results/*.txt; do
        if [[ ${files[@]} =~ $(basename ${file//-/}) ]]; then
            cut -f6 $file > ./6th_col/$(basename $file).6th_col
        fi
    done

    paste ./6th_col/*.6th_col > ${output_name}.counts
    rm -rf 6th_col
}

filter types_inputs1[@] times_inputs1[@] treatments_inputs1[@]
alignment_count files[@] "combined_group_1"

filter types_inputs2[@] times_inputs2[@] treatments_inputs2[@]
alignment_count files[@] "combined_group_2"
