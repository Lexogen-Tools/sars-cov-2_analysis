idemuxCPP --i1-read=1 --i1-start=1 -1 Resources/testdata_R1.fastq.gz -s sample-sheet.csv -w 30 -p 30 -o /home/gwiktorin/tests/SARS-CoV-2//idemuxxed/
xargs -a /home/gwiktorin/tests/SARS-CoV-2//parallel_jobs.sh -d $'\n' -n1 -P30 bash -c
