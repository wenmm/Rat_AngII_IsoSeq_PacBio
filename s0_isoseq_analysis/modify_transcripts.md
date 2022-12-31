### 1. SQANTI3 annotate and remove intrapriming isoforms
```
python sqanti3_qc.py --polyA_motif_list polya.txt -t 40 -o combined -d combined_sqanti3 -c ../2rd_analysis/total.rn7.SJ.out.tab --genename combined.collapsed.gff ncbiRefSeq.gtf /mnt/data1/refgenome/ucsc_mRatBN7_2/rn7.fa

python sqanti3_RulesFilter.py combined_classification.txt combined_corrected.fasta combined_corrected.gtf.cds.gff

python filter_gff_PB.py intrapriming_PB_list combined_corrected.gtf.cds.gff > combined_corrected.gtf.cds.s1.gff

```
### 2. mark isoforms with the same structure
```
gtfToGenePred combined_corrected.gtf.cds.s1.gff combined_corrected.gtf.cds.s1.gp -geneNameAsName2 -genePredExt

python self_compare.py combined_corrected.gtf.cds.s1.gp combined_corrected.gtf.cds.s1.gff > combined_corrected.gtf.cds.s1.selfcompare

grep -w exact combined_corrected.gtf.cds.s1.selfcompare > combined_corrected.gtf.cds.s1.selfcompare.exact

python net.py combined_corrected.gtf.cds.s1.selfcompare.exact > combined_corrected.gtf.cds.s1.selfcompare.exact.net
```

### 3. Select the most representative isoform (least distance)
```
python calculate_distance.py ref.pas.pos rn7_3READ_PolyASeq_summits.bed iso_seq_isoform_pas.bed > iso_seq_isoform_pas.distance
python calculate_distance.py ref.tss.pos rn7.cage_summits.bed iso_seq_isoform_tss.bed > iso_seq_isoform_tss.distance
```

### 5. program list
- SQANTI3-4.2 (https://github.com/ConesaLab/SQANTI3/tags)
- gtfToGenePred (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)
- minimap2 (2.17-r941) (https://github.com/lh3/minimap2)
- fmlrc2 (0.1.6) (https://github.com/HudsonAlpha/fmlrc2)