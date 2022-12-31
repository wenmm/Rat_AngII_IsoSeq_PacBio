## Input files and format

| Format      | Description | Example     |
| :---        |    :----:   |          ---: |
| plain text      | mapping reads file       | mapping_depth.txt   |
| Bedgraph   | files in this format store the reads alignment information, which can be converted from BAM files by bedtools (e.g.:"bedtools genomecov -ibam *.bam -bga -split -trackline") | aligned_bg_files      |
| BED   | a bed file contains position of 3'UTR | 3utr_for_dapars.bed|
| plain text      | file contains position of tandem 3'UTR (Long 3UTRs correspond to short 3UTRs)        | 3utr_for_dapars.middle |

## Usage

```
for i in `less chrList.txt`
do
	grep -w $i 3utr_for_dapars.bed > $i\_3utr_for_dapars.bed
	python pacbio_pdui.py 3utr_for_dapars.middle mapping_depth.txt aligned_bg_files $i\_3utr_for_dapars.bed > $i\_dapars_result
doen
```