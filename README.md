# SARS-CoV-2 analysis script
**This analysis script is implemented as a bash routine and is intended to work in combination with Lexogen's SARS-CoV-2 test kit**

Lexogen developed a mass screening test to detect copies of the SARS-CoV-2 virus in human samples and is optimized for the application of thousands of samples in one sequencing run. To enable a quick and structured analysis of the samples after they were sequenced and demultiplexed, Lexogen provides this analysis script as a working solution and recommendation on the analysis procedure.
## Requirements
The script uses the following software tools and scripts in its routine:
```
    - idemuxCPP (https://github.com/Lexogen-Tools/idemuxcpp)
    - samtools (https://www.htslib.org/)
    - cutadapt (https://cutadapt.readthedocs.io/en/stable/)
    - bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    - featureCounts (http://subread.sourceforge.net/)
    - GNU awk (https://www.gnu.org/software/gawk/)
    - python 3.6 (https://www.python.org/downloads/)
```
**You can either install iDemuxCPP by using a pre-compiled package that matches your operating system or you can compile it yourself by following the instructions on the installation page of iDemuxCPP on github (https://github.com/Lexogen-Tools/idemuxcpp#1-installation).**

The other tools can either be installed and made available in the user's PATH or installed into a conda environment and used in this environment. If you use your own installations, please carefully read the content of the environment.yml file to check if the versions of your installed tools are compatible with the versions listed in this file. Also, please ensure that the tools are indeed available in your path.

You can use conda (anaconda, miniconda) to install the tools into a new conda environment and use this environment, 'lex_sars-cov-2', by calling 
```
$ conda env create -f /path/to/environment.yml && conda activate lex_sars-cov-2
```
in your command shell. The script further uses commands which come with bash: paste, cut, getopt, mv, rm, xargs. We assume that these are available to you.
## How to run this script
Provided you have the required tools installed, cloned the content of the git repository onto your machine and changed into this cloned directory, you can start the script with a call like
```
$ ./start_analysis.sh -i path/to/demultiplexed_data/R1.fastq.gz -s path/to/sample-sheet.csv -o path/to/output_directory/ -p [nr of parallel jobs]
```
**path/to/demultiplexed_data/R1.fastq.gz** is the path to a gzipped fastq file, which contains the reads of the sample(s) of this experiment.

**path/to/sample-sheet.csv** is the path to a sample sheet that lists the samples in the fastq file and their assigned barcodes to use for demultiplexing.

**path/to/output_directory** is the path to a directory to write intermediate analysis results, the quantification summary and the prediction table.

**[nr of parallel jobs]** is the number of parallel jobs to be processed with xargs and the number of threads available to idemux. Depending on the resources available on your machine, increasing this value might significantly shorten the duration of the analysis.

## The analysis procedure
The analysis script is a bash routine that uses minimal user input to process demultiplexed, gzipped fastq files. The script is written with the assumption that a relatively high number of samples are to be processed with this script, with moderate to low read depth per individual sample. The script uses xargs to submit tasks to be processed in parallel to sub processes, if possible. We recommend to use the script on a workstation or server that has a high number of available CPUs to shorten the effective duration of the calculations with a sensible use of the -p argument (There is little use in operating the script with 30 parallel jobs on a machine with 4 CPU cores).

The following sections outline the analysis steps that are implemented in the script.
### Demultiplexing with idemux
```
$ idemuxCPP --i1-read=1 --i1-start=1 -1 path/to/demultiplexed_data/R1.fastq.gz -s path/to/sample-sheet.csv -w [nr of threads] -p [nr of threads] -o path/to/outputdir/idemuxxed/
```
The library preparation includes the addition of inline barcodes to enhance the sample resolution in the RNA-Seq experiment. The script therefore starts with a call of iDemuxCPP (https://github.com/Lexogen-Tools/idemuxcpp) to demultiplex the pooled samples from the fastq file into separate fastq files per sample. The provided sample-sheet.csv file contains a set of 120 inline barcode sequences (i1) and generic sample names (sample_name). The names will be used for temporary files and directories. You can modify the names to make your samples easily identifiable. Please do not use characters that are not allowed for file and directory names in a *nix environment in the sample names.

If you are familiar with the usage of iDemuxCPP and know how to demultiplex your sequencing run into a single file, preserving the index information in your read headers, you can instead provide the result of this demultiplexing and list the all the samples of your run in the sample-sheet.csv file with their complete combination of i7, i5, and i1 barcode sequences instead. The analysis script will then process all of them in one call of the script. Please see the iDemuxCPP Readme (https://github.com/Lexogen-Tools/idemuxcpp) for more information.
The following is an example for the format of a sample-sheet
```
sample_name,i7,i5,i1
sample_001,,,CGGGAACCCGCA
sample_002,,,AAACGTTCATCC
sample_003,,,TTGTCCGATATG
sample_004,,,GTTTAAAGGCAG
sample_005,,,CCAAAGAGGGAT
sample_006,,,TCCTCTCTTCTA
sample_007,,,GAAGGGTAAAGC
sample_008,,,AGTCTCAGCAAA
sample_009,,,ATCGACTTGTGT
sample_010,,,AACCCTGGGAAG
sample_011,,,AGGTGGTTCTAC
sample_012,,,TACGCCCACGTG
...
```
### Trimming with cutadapt
```
$ cutadapt --quiet -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 path/to/outputdir/idemuxxed/${sample}.fastq.gz | \
> cutadapt --quiet -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | \
> cutadapt --quiet -m 20 -O 3 -a "r1polyA=A{18}" - | \
> cutadapt --quiet -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o path/to/outputdir/analysis/${sample}/R1.fastq.gz -
```
Prior to the alignment the reads are analyzed and potentially trimmed with cutadapt to increase the mappability of the sample. Read-outs of fragments which end up too short contain adapter sequences which do not match to the reference and would end up unmappable. Reads which result in a length <20 or consist entirely of adapter sequences are removed. The script calls cutadapt 4 times in succession, which are connected by piping the result of one call into the subsequent one. 
### Alignment with bowtie2
```
$ bowtie2 -x targets_SIRV109/sequences --mm -U path/to/outputdir/analysis/${sample}/R1.fastq.gz --quiet -S path/to/outputdir/analysis/${sample}/Aligned.out.sam
```
The library preparation generates reads of 4 locations on the SARS-CoV-2 virus and on one location of a spike-in control of Lexogen's RNA-Seq spike in controls (SIRVs, https://www.lexogen.com/sirvs/). The alignment is performed with bowtie2. The standard settings of bowtie already satisfy the requirements for this analysis.
### Counting with featureCounts
```
$ featureCounts -a targets/annotation.gtf -o path/to/outputdir/analysis/${sample}/src.byfeature -t feature -g target_id -M --fraction -p -B -d 30 -C path/to/outputdir/analysis/${sample}/Aligned.out.sam
```
The alignments are counted with featureCounts using the reference file annotation.gtf. This file contains features for each target. We expect the reads to align in anti-sense direction to the reference due to the library prepration. The quantification with featueCounts therefore counts the reads on the reverse strand (argument -s 2).
### Summarization and assessment
The summary of the quantification is written into a count matrix which lists the counts per target in the individual samples. This count matrix is used to assess the presence of the virus in the samples based on the read counts.
## Results
The results of the analysis pipeline are two tables in tab-separated text format
- A summary of the quantification per sample per target, **counts_summary.tsv**
```
id      input reads     unmapped reads  ORF_1a  5p_spike_cds    O-linked_glycan_residue_region  charite-similar_gene_E  SIRV109
sample_0001   24187   792     4448    6603    3133    7329    1882
sample_0002   42235   1673    15222   3321    5325    7865    8829
sample_0003   42433   812     13131   9829    3494    10878   4289
sample_0004   30549   20      11860   2111    5719    7195    3644
sample_0005   39564   2547    14323   9648    3479    7763    1804
sample_0006   42899   1547    17100   7554    8257    4441    4000
sample_0007   25251   18      6945    3928    1957    4780    7623
sample_0008   48188   1430    9306    9435    5273    4797    17947
sample_0009   0       0       0       0       0       0       0
sample_0010   0       0       0       0       0       0       0
sample_0011   0       0       0       0       0       0       0
sample_0012   0       0       0       0       0       0       0
...
```
- An assessment of the presence of the virus in the sample, **prediction.tsv**
```
id      status
sample_0001   True
sample_0002   True
sample_0003   True
sample_0004   True
sample_0005   True
sample_0006   True
sample_0007   True
sample_0008   True
sample_0009   False
sample_0010   False
sample_0011   False
sample_0012   False
...
```
