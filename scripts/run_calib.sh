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

### Split the FASTQ into READ of different lenght depending on the parameters
tab_calib_readlength=($(echo "${calib_readlength}" | sed 's/,/ /g'))
tab_calib_e=($(echo "${calib_e}" | sed 's/,/ /g'))
tab_calib_k=($(echo "${calib_k}" | sed 's/,/ /g'))
tab_calib_m=($(echo "${calib_m}" | sed 's/,/ /g'))
tab_calib_t=($(echo "${calib_t}" | sed 's/,/ /g'))

for((ilength=0; ilength < ${#tab_calib_readlength[*]}; ilength++))
do
    ### Define the min and max lenght of read sequence
    istart=${tab_calib_readlength[ ${ilength} ]}
    if [[ $((${ilength}+1)) -eq ${#tab_calib_readlength[*]} ]]
    then
        iend=1000
    else
        iend=${tab_calib_readlength[ $((${ilength} + 1)) ]}
    fi

    ### Select only the reads of the good length
    zcat "${r1}" | \
        awk -v imin=${istart} -v imax=${iend} \
        'NR%4==1 {pline=$0}; NR%4==2 && length($0)>=imin && length($0)<imax {print pline; print $0; getline; print; getline; print}' > "${sample}_R1_nbr${istart}.fastq"
    zcat "${r2}" | \
        awk -v imin=${istart} -v imax=${iend} \
        'NR%4==1 {pline=$0}; NR%4==2 && length($0)>=imin && length($0)<imax {print pline; print $0; getline; print; getline; print}' > "${sample}_R2_nbr${istart}.fastq"

    ### Run the cluster (molecule) detection
    calib \
        --threads ${cpu} \
        -l ${calib_umilength} \
        -f "${sample}_R1_nbr${istart}.fastq" \
        -r "${sample}_R2_nbr${istart}.fastq" \
        -e ${tab_calib_e[ ${ilength} ]} \
        -k ${tab_calib_k[ ${ilength} ]} \
        -m ${tab_calib_m[ ${ilength} ]} \
        -t ${tab_calib_t[ ${ilength} ]} \
        -o "${sample}_nbr${istart}."

    ### Get the consensus
    calib_cons \
        --threads ${cpu} \
        -c "${sample}_nbr${istart}.cluster" \
        -q "${sample}_R1_nbr${istart}.fastq" \
        -q "${sample}_R2_nbr${istart}.fastq" \
        -o "${sample}_R1_nbr${istart}.consensus" \
        -o "${sample}_R2_nbr${istart}.consensus"
done

### Merge the consensus and remove the UMI sequence for the sequence and from the quality
outR1=$(basename "${r1}" .fastq.gz)".consensus.fastq.gz"
outR2=$(basename "${r2}" .fastq.gz)".consensus.fastq.gz"

awk -v umilength=${calib_umilength} '{ if(NR%4==2 || NR%4==0){ print substr($0,umilength+1) }else{ print $0 } }' "${sample}_R1_nbr"*".consensus.fastq" | gzip > "${outR1}"
awk -v umilength=${calib_umilength} '{ if(NR%4==2 || NR%4==0){ print substr($0,umilength+1) }else{ print $0 } }' "${sample}_R2_nbr"*".consensus.fastq" | gzip > "${outR2}"

### Clean the temporary files
rm -f "${sample}_R?_nbr*.fastq" "${sample}_R?_nbr*.consensus.fastq"
