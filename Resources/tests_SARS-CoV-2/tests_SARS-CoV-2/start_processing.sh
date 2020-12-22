idemuxCPP --i1-read=1 --i1-start=1 -1 Resources/testdata_R1.fastq.gz -s sample-sheet.csv -w 30 -p 30 -o tests_SARS-CoV-2/idemuxxed/
xargs -a tests_SARS-CoV-2/parallel_jobs.sh -d $'\n' -n1 -P30 bash -c
