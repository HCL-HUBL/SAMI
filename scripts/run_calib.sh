#!/usr/bin/env bash

set -eo pipefail

usage="Usage: run_calib.sh -s <SAMPLE> -1 <R1> -2 <R2> -l <UMI_LENGTH> -r <READ_LENGTH> -e <CALIB_E_PARAMETER> -k <CALIB_K_PARAMETER> -m <CALIB_M_PARAMETER> -t <CALIB_T_PARAMETER> -c <NBR_CPU>"

### Read options
while getopts ":s:1:2:l:r:e:k:m:t:c:" opt
do
	case ${opt} in
		s) sample="$OPTARG"
		   ;;
		1) r1="$OPTARG"
		   ;;
		2) r2="$OPTARG"
		   ;;
		l) calib_umilength="$OPTARG"
		   ;;
		r) calib_readlength="$OPTARG"
		   ;;
		e) calib_e="$OPTARG"
		   ;;
		k) calib_k="$OPTARG"
		   ;;
		m) calib_m="$OPTARG"
		   ;;
		t) calib_t="$OPTARG"
		   ;;
		c) cpu="$OPTARG"
		   ;;
		*) printf "Invalid option -%s %s\n%s\n" "${opt}" "${OPTARG}" "${usage}" >&2
		   exit 1
		   ;;
	esac
done

if [[ -z ${sample} || -z ${r1} || -z ${r2} || -z ${calib_umilength} || -z ${calib_readlength} || -z ${calib_e} || -z ${calib_k} || -z ${calib_m} || -z ${calib_t} || -z ${cpu} ]]
then
	printf "One mandatory argument is empty\n%s\n" "${usage}" >&2
	exit 2
fi

### Split parameters per read bin
tab_calib_readlength=($(echo "${calib_readlength}" | sed 's/,/ /g'))
tab_calib_e=($(echo "${calib_e}" | sed 's/,/ /g'))
tab_calib_k=($(echo "${calib_k}" | sed 's/,/ /g'))
tab_calib_m=($(echo "${calib_m}" | sed 's/,/ /g'))
tab_calib_t=($(echo "${calib_t}" | sed 's/,/ /g'))

### Decompress FASTQ files
gunzip --stdout "${r1}" > "${r1%.gz}"
gunzip --stdout "${r2}" > "${r2%.gz}"

### Bin reads of different lengths depending on the parameters
BinReads "${r1%.gz}" "${r2%.gz}" "0,${calib_readlength}"

### Output file names
outR1="${sample}_R1.consensus.fastq"
outR2="${sample}_R2.consensus.fastq"
rm -f "${outR1}"
rm -f "${outR2}"

for((i=0; i < ${#tab_calib_readlength[*]}; i++))
do
	### Define the min length of read sequences
	min_size=${tab_calib_readlength[ ${i} ]}
	
	### Run the cluster (molecule) detection
	calib \
		--threads ${cpu} \
		-l ${calib_umilength} \
		-f "R1_gte${min_size}.fastq" \
		-r "R2_gte${min_size}.fastq" \
		-e ${tab_calib_e[ ${i} ]} \
		-k ${tab_calib_k[ ${i} ]} \
		-m ${tab_calib_m[ ${i} ]} \
		-t ${tab_calib_t[ ${i} ]} \
		-o "${sample}_gte${min_size}." 1> calib_${min_size}.stdout 2> calib_${min_size}.stderr
	
	### Get the consensus
	calib_cons \
		--threads ${cpu} \
		-c "${sample}_gte${min_size}.cluster" \
		-q "R1_gte${min_size}.fastq" \
		-q "R2_gte${min_size}.fastq" \
		-o "R1_gte${min_size}.consensus" \
		-o "R2_gte${min_size}.consensus" 1> calib-cons_${min_size}.stdout 2> calib-cons_${min_size}.stderr

	### Change the cluster read name to avoid issue
	### Remove UMIs
	### Concatenate
	awk -v umilength=${calib_umilength} -v min_size=${min_size} 'NR%4 == 1 { $1=$1"_gte"min_size } NR%2 == 0 { $0 = substr($0, umilength+1) } { print }' "R1_gte${min_size}.consensus.fastq" >> "${outR1}"
	awk -v umilength=${calib_umilength} -v min_size=${min_size} 'NR%4 == 1 { $1=$1"_gte"min_size } NR%2 == 0 { $0 = substr($0, umilength+1) } { print }' "R2_gte${min_size}.consensus.fastq" >> "${outR2}"
done

### Compress
gzip "${outR1}"
gzip "${outR2}"

### Cleanup
rm -f *.fastq

