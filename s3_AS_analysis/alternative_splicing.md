# salmon calculate isoform TPM (salmon version 1.8.0)

salmon index -t combined_f1_f.cds_add_ncbiRefSeq.fa -i combined_f1_f.cds_add_ncbiRefSeq.idx -p 50 --keepDuplicates
salmon quant -p 50 -i combined_f1_f.cds_add_ncbiRefSeq.idx --gcBias -l A -1 R1.fastq -2 R2.fastq --validateMappings -o ./salmon/sample1

# SUPPA2 identify splicing events

python /SUPPA-2.3/suppa.py generateEvents -b V -i combined_f1_f.cds_add_ncbiRefSeq.filter.gtf -o events -e SE SS MX RI FL -f ioe

## AF AL USE Variable 50
## A3 A5 MX RI SE use strict

cat ../suppa2/events_A3_strict.ioe ../suppa2/events_A5_strict.ioe ../suppa2/events_MX_strict.ioe ../suppa2/events_RI_strict.ioe ../suppa2/events_SE_strict.ioe events_AF_variable_50.ioe events_AL_variable_50.ioe > events_variable_50.ioe

# Calculate PSI for each splicing events

python /SUPPA-2.3/suppa.py psiPerEvent -i events_variable_50.ioe -e [tpm expression file from salmon] -o events

