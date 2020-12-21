#!/usr/bin/env bash

installdir=$(dirname $(readlink -f $0))

referenceSequences=${installdir}"/targets_SIRV109/sequences.fa"
bowtie2ref=${installdir}"/targets_SIRV109/sequences"
gtffile=${installdir}"/targets_SIRV109/annotation.gtf"

idemuxCPPbin=idemuxCPP
bowtie2bin=bowtie2
cutadaptbin=cutadapt
featurecountsbin=featureCounts

PRE_nparallel=4
PRE_samplesheet="sample-sheet.csv"

usage() { echo "Usage: $0 [-i|--inputFile <R1.fastq.gz>] [-o|--outputDirectory <path/to/outputdirectory>] (Optional: [-p|--parallel] <# of parallel jobs/threads> [-s|--samplesheet </path/to/sample-sheet.csv>])" 1>&2; exit 1; }

OPT=$(getopt -o i:o:p:s: --long inputFile:,outputDirectory:,parallel:,sampleSheet:, -- "$@")

eval set -- "$OPT"

while true; do
  case "$1" in
    -i | --inputFile ) read1="$2"; shift 2;;
    -o | --outputDirectory ) outdir="$2"; shift 2;;
    -p | --parallel ) nparallel="$2"; shift 2;;
    -s | --sampleSheet ) samplesheet="$2"; shift 2;;
    -- ) shift; break;;
    * ) echo "invalid state" ; usage;;
  esac
done

[ $# -gt 0 ] && { echo "no command line arguments outside options." 1>&2; usage; }

[ -z ${read1+x} ] && { echo "path to R1 has to be specified." 1>&2; usage; }
[ -z ${outdir+x} ] && { echo "path to output directory has to be specified." 1>&2; usage; }

[ -z ${nparallel+x} ] && { nparallel=${PRE_nparallel} ; }
[ -z ${samplesheet+x} ] && { samplesheet=${PRE_samplesheet} ; }

master_script=${outdir}"/start_processing.sh"
worker_script=${outdir}"/parallel_jobs.sh"
summaryfile=${outdir}"/counts_summary.tsv"

mkdir -p ${outdir}'/idemuxxed'

# ### master script
# demux with idemuxCPP
echo "${idemuxCPPbin} --i1-read=1 --i1-start=1 -1 ${read1} -s ${samplesheet} -w ${nparallel} -p ${nparallel} -o ${outdir}/idemuxxed/" > ${master_script}
# write a call of the worker script to the master script, which constitutes a call of the worker script with parallel.
echo "xargs -a ${worker_script} -d $'\n' -n1 -P${nparallel} bash -c" >> ${master_script}

# ### worker script
# all tool calls which can be processed in parallel are placed into the worker script
samples=$(awk -F ',' 'NR==1{tc=1 ; for(ii=1;ii<=NF;ii++){ if($ii=="sample_name"){ tc=ii } } } ; NR>1{print $tc}' ${samplesheet})
for sample in ${samples[@]}; do
    outprefix=${outdir}"/analysis/"${sample}/
    # trimming with cutadapt, alignment with bowtie2, counting with featurecounts; per sample
    echo "mkdir -p ${outprefix} && ${cutadaptbin} --quiet -m 20 -O 20 -a \"polyA=A{20}\" -a \"QUALITY=G{20}\" -n 2 ${outdir}/idemuxxed/${sample}.fastq.gz | ${cutadaptbin} --quiet -m 20 -O 3 --nextseq-trim=10 -a \"r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000\" - | ${cutadaptbin} --quiet -m 20 -O 3 -a \"r1polyA=A{18}\" - | ${cutadaptbin} --quiet -m 20 -O 20 -g \"r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20\" --discard-trimmed -o ${outprefix}/R1.fastq.gz - && ${bowtie2bin} -x ${bowtie2ref} -L 10 --mm -U ${outprefix}/R1.fastq.gz --quiet -S ${outprefix}/Aligned.out.sam && ${featurecountsbin} -a ${gtffile} -o ${outprefix}/src.byfeature -t feature -g target_id -M --fraction -p -B -d 30 -s 2 -C ${outprefix}/Aligned.out.sam 2> /dev/null && samtools view -f 4 ${outprefix}/Aligned.out.sam | wc -l > ${outprefix}/unmapped.count.txt && zcat ${outprefix}/R1.fastq.gz | grep -e \"^@\" | wc -l > ${outprefix}/input.count.txt"
done > ${worker_script}

chmod u+x ${master_script}
# all ready, call the master script to demultiplex, trim, align and count the samples.
${master_script}

# All samples which contained valid reads were aligned and counted. The following section loads the results of the alignment and counting into a summary table ${summaryfile}
tmp_results=$(mktemp ${outdir}/XXXXXXtmp)
rm -f ${summaryfile}
echo "id" > ${summaryfile}
echo "input reads" >> ${summaryfile}
echo "unmapped reads" >> ${summaryfile}
grep -P -e "^>" ${referenceSequences} | sed 's/^>//' >> ${summaryfile}
for sample in ${samples[@]}
do
    outprefix=${outdir}"/analysis/"${sample}
    echo $sample
    echo ${sample} > ${summaryfile}_tmp
    if [[ -f "${outprefix}/src.byfeature" ]]; then
        nir=$(cat ${outprefix}/input.count.txt)
        nonaligned=$(cat ${outprefix}/unmapped.count.txt)
        echo -e ${nir}"\n"${nonaligned} >> ${summaryfile}_tmp
        cut -f 7 ${outprefix}/src.byfeature | tail -n+3 >> ${summaryfile}_tmp
    fi
    paste ${summaryfile} ${summaryfile}_tmp > ${tmp_results}
    mv ${tmp_results} ${summaryfile}
done
awk -F '\t' 'BEGIN{nrmax=0 ; nfmax=0 } { for(ii=1;ii<=NF;ii++){content[NR][ii]=$ii} ; nrmax=NR ; if(NF>nfmax){nfmax=NF} } END{ for(jj=1;jj<=nfmax;jj++){ for(ii=1;ii<=nrmax;ii++){  printf "%s\t",content[ii][jj] } ; printf "\n" } }' ${summaryfile} > ${summaryfile}_tmp
mv ${summaryfile}_tmp ${summaryfile}

# assess SARS-CoV-2 status in samples
python3 ${installdir}/predict.py --input ${summaryfile} --output ${outdir}/prediction.tsv --sample-sheet ${samplesheet} 
