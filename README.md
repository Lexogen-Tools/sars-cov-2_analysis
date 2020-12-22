# SARS-CoV-2 analysis script

**This analysis script is implemented as a bash routine and is intended to work in combination with Lexogen's SARS-CoV-2 test kit**

The SARS-CoV-2 mass screening test developed by Lexogen currently analyzes up to 9216 samples in a single RNA-Seq run. Each sample in a run is identified via one of 96 unique dual i7/i5 indices and one of 96 Lexogen i1 sample indices (for further information please visit https://www.lexogen.com/quantseq-pool-sample-barcoded-3mrna-sequencing/). The script [start_analysis.sh](start_analysis.sh) in this repository allows fast and easy classification of each sample after sequencing and i7/i5 demultiplexing. The script can also be used to perform combined i7/i5 and i1 demultiplexing and therefore process a complete RNA-Seq run. The latter requires demultiplexing the run into a single fastq file, which retains the i7/i5 index information in the read header (see section "Demultiplexing with idemux").

## Installation
The script [start_analysis.sh](start_analysis.sh) uses the following software tools and scripts in its routine:
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

The other tools can either be installed and made available in the user's PATH or installed into a conda environment and used in this environment. If you use your own installations, please make sure that the versions of your installed tools correspond to the versions in [environment.yml](environment.yml). Also, please ensure that the tools are available in your path.

You can use conda (anaconda, miniconda) to install the tools into a new conda environment and use this environment, e.g. named 'lex_sars-cov-2', by calling 
```
$ conda env create -f /path/to/environment.yml && conda activate lex_sars-cov-2
```
The script further uses commands which come with bash: paste, cut, getopt, mv, rm, xargs. We assume that these are available on your system.
## How to run this script
Provided you have the required tools installed, cloned the content of the git repository onto your machine and changed into this cloned directory, you can start the script with a call like
```
$ ./start_analysis.sh -i path/to/demultiplexed_data/R1.fastq.gz -s path/to/sample-sheet.csv -o path/to/output_directory/ -p [nr of parallel jobs]
```
**path/to/demultiplexed_data/R1.fastq.gz** is the path to a gzipped fastq file, which contains the reads of the sample(s) of this experiment.

**path/to/sample-sheet.csv** is the path to a sample sheet that lists the samples in the fastq file and their associated Lexogen i1 sample barcodes (see section "Demultiplexing with idemux").

**path/to/output_directory** is the path to a directory to write intermediate analysis results, the quantification summary and the prediction table.

**[nr of parallel jobs]** is the number of parallel jobs to be processed with xargs and the number of threads available to idemux. Depending on the resources available on your machine, increasing this value might significantly shorten the duration of the analysis.

## The analysis procedure
The analysis script is a bash routine, which processes i7/i5 demultiplexed, gzipped fastq files. The script allows to analyze up to 96 samples within a single gzipped fastq file. We recommend a low to moderate read depth for each sample and to use a server with a sufficient number of cores for effective parallelization. The maximal number of parallel processes used by [start_analysis.sh](start_analysis.sh) can be specified with the -p argument.

The following sections outline the analysis steps implemented in the script.
### Demultiplexing with idemux
```
$ idemuxCPP --i1-read=1 --i1-start=1 -1 path/to/demultiplexed_data/R1.fastq.gz -s path/to/sample-sheet.csv -w [nr of threads] -p [nr of threads] -o path/to/outputdir/idemuxxed/
```
In addition to the use of i7/i5 indices, the library preparation introduces Lexogen's i1 sample barcodes (see https://www.lexogen.com/quantseq-pool-sample-barcoded-3mrna-sequencing/) to increase the number of simultaneously analyzed samples. The script starts after i7/i5 demultiplexing has been performed, e.g. with bcl2fastq or, more preferably, Lexogen's iDemuxCPP software (https://github.com/Lexogen-Tools/idemuxcpp). The first step in the script uses iDemuxCPP to demultiplex the gzipped fastq file according to a sample sheet which lists the names of the samples together with their associated Lexogen i1 sample barcode. A sample sheet for a set of 96 samples is given in [sample-sheet.csv](sample-sheet.csv), the first lines of which are shown below. Please note that sample names are used for the creation of temporary files and directories and must therefore conform with file name requirements in a *nix environment.

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
If you are familiar with the use of iDemuxCPP and know how to demultiplex your sequencing run into a single fastq file, preserving the index information in your read headers, you can instead 

* provide this single demultiplexed fastq file to the analysis script and 
* list all the samples of your run in the sample-sheet.csv file together with their complete combination of i7, i5, and i1 barcode sequences

The analysis script will then perform demultiplexing of i7/i5 and i1 indices together and analyze all corresponding samples in one go. For further information on how to use iDemuxCPP and demutliplex into a single fastq file preserving the index information in read headers, please refer to the iDemuxCPP Readme (https://github.com/Lexogen-Tools/idemuxcpp).

### Trimming with cutadapt

```
$ cutadapt --quiet -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 path/to/outputdir/idemuxxed/${sample}.fastq.gz | \
> cutadapt --quiet -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | \
> cutadapt --quiet -m 20 -O 3 -a "r1polyA=A{18}" - | \
> cutadapt --quiet -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o path/to/outputdir/analysis/${sample}/R1.fastq.gz -
```
Prior to alignment reads are trimmed with cutadapt. Reads which result in a length <20 or consist entirely of adapter sequences are removed.
### Alignment with bowtie2
```
$ bowtie2 -x targets_SIRV109/sequences --mm -L 10 -U path/to/outputdir/analysis/${sample}/R1.fastq.gz --quiet -S path/to/outputdir/analysis/${sample}/Aligned.out.sam
```
The library preparation generates reads at 4 locations on the SARS-CoV-2 virus genome and on one location of a spike-in control of Lexogen's RNA-Seq spike in controls (SIRVs, https://www.lexogen.com/sirvs/). The alignment is performed with bowtie2. 
### Counting with featureCounts
```
$ featureCounts -a targets/annotation.gtf -o path/to/outputdir/analysis/${sample}/src.byfeature -t feature -g target_id -M --fraction -p -B -d 30 -s 2 -C path/to/outputdir/analysis/${sample}/Aligned.out.sam
```
The alignments are counted with featureCounts using the reference file [annotation.gtf](annotation.gtf). This file contains features for each target. We expect the reads to align in anti-sense direction to the reference due to the library preparation. The quantification with featureCounts therefore counts the reads on the reverse strand (argument -s 2).

### Sample classification

The read counts per target are used to classify samples as either positive or negative. This is done with script [predict.py](predict.py) and command line

```
$ python3 ${installdir}/predict.py --input ${summaryfile} --output ${outdir}/prediction.tsv --sample-sheet ${samplesheet}
```

### Output

The main output of the script [start_analysis.sh](start_analysis.sh) is a table with the positive status for each sample, i.e. a status of TRUE corresponds to positive and a status of FALSE corresponds to negative. An example for such a table is given below

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

The script also produces a table with the target read counts per sample. This can be used for further analysis or for adjustment of the parameters in [predict.py](predict.py), which influence the sensitivity and specificity of the classification.

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
## Testing the pipeline

To verify that the pipeline has been properly installed, the gzipped fastq file in [Resources/testdata_R1.fastq.gz](Resources/testdata_R1.fastq.gz) can be used as a test case. The expected outcome of this test is that samples *sample_009-sample_016* are positive, while the remaining samples in *samples_001-samples_096* are negative. Below you can find the command line for performing this test and the first 20 lines of the test output. The complete test results can be found in [Resources/prediction.tsv](Resources/prediction.tsv).

```
$ mkdir -p tests_SARS-CoV-2 && ./start_analysis.sh -i Resources/testdata_R1.fastq.gz -s sample-sheet.csv -o tests_SARS-CoV-2 -p 30;
$ head -n 20 tests_SARS-CoV-2/prediction.tsv 
id      status
sample_001      False
sample_002      False
sample_003      False
sample_004      False
sample_005      False
sample_006      False
sample_007      False
sample_008      False
sample_009      True
sample_010      True
sample_011      True
sample_012      True
sample_013      True
sample_014      True
sample_015      True
sample_016      True
sample_017      False
sample_018      False
sample_019      False
```
