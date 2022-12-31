### 1. calculate isoform FLNC reads support
```
/tools/cDNA_Cupcake/sequence/sam_to_gff3.py -i control_corrected_final.fa -s rn7 control_corrected_final.sam
/tools/cDNA_Cupcake/sequence/sam_to_gff3.py -i angii_corrected_final.fa -s rn7 angii_corrected_final.sam

python self_compare50.py control.flnc.aligned.gname control.flnc.aligned.gtf combined_f1_f.cds_add_ncbiRefSeq.gtf > control.aligned.flnc
python self_compare50.py angii.flnc.aligned.gname angii.flnc.aligned.gtf combined_f1_f.cds_add_ncbiRefSeq.gtf > angii.aligned.flnc

grep -w exact control.aligned.flnc > control.aligned.flnc.exact
grep -w exact angii.aligned.flnc > angii.aligned.flnc.exact

python count.py control.aligned.flnc.exact > control.aligned.flnc.exact.count
python count.py angii.aligned.flnc.exact > angii.aligned.flnc.exact.count

python calculate_percent.py PacBio_angii_control_tr_count > PacBio_angii_control_tr_count_percent

Rscript isoform_add_tag.r
```