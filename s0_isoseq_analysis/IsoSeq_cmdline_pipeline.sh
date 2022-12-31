tag=$1

## Step 1 - Circular Consensus Sequence calling
ccs --numThreads 15 --noPolish --minLength=50 --minPasses=1 --maxPoaCoverage 10 --minPredictedAccuracy=0.7 --minZScore=-999 --maxDropFraction=0.8 --minPredictedAccuracy=0.8 --minSnr=3.75 ${tag}.subreads.bam ${tag}.ccs.bam

## Step 2 - Primer removal and demultiplexing
$lima ${tag}.ccs.bam primer.fasta ${tag}.fl.bam --isoseq --no-pbi

## Step 3 - Refine
isoseq3 refine ${tag}.fl.primer_5p--primer_3p.bam primer.fasta ${tag}.flnc.bam --require-polya

## Step 3b - Merge SMRT Cells
dataset create --type TranscriptSet merged_new.flnc.xml /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171016_100035/1_A01/m54069_171016_112355.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171016_100035/2_A02/m54069_171016_212317.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171016_100035/3_A03/m54069_171017_073307.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171016_100035/4_A04/m54069_171017_174259.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171018_103019/1_B01/m54069_171018_104524.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171018_103019/2_C01/m54069_171018_205255.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171018_103019/3_D01/m54069_171019_071229.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171018_103019/4_E01/m54069_171019_173231.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171025_131159/1_F01/m54069_171025_132032.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171025_131159/2_G01/m54069_171025_232753.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/control/r54069_20171025_131159/3_H01/m54069_171026_094732.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/angll/r54069_20171021_023051/1_B02/m54069_171021_023917.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/angll/r54069_20171021_023051/2_B03/m54069_171021_124654.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/angll/r54069_20171021_023051/3_B04/m54069_171021_230647.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/angll/r54069_20171021_023051/4_B05/m54069_171022_092642.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/angll/r54069_20171023_054927/1_C02/m54069_171023_055752.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/angll/r54069_20171023_054927/2_C03/m54069_171023_160527.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/angll/r54069_20171023_054927/3_C04/m54069_171024_022511.flnc.new.bam /mnt/data2/rat_Ang_3seq/sequel_raw_data/angll/r54069_20171023_054927/4_C05/m54069_171024_124507.flnc.new.bam

## Step 4 - Clustering
isoseq3 cluster -j 20 merged.flnc.xml unpolished.bam --split-bam 12

## Step 5 - Parallel Polishing
isoseq3 polish -j 30 unpolished.0.bam merged.subreadset.xml polished.0.bam
isoseq3 polish -j 30 unpolished.1.bam merged.subreadset.xml polished.1.bam
isoseq3 polish -j 30 unpolished.2.bam merged.subreadset.xml polished.2.bam
isoseq3 polish -j 30 unpolished.3.bam merged.subreadset.xml polished.3.bam
isoseq3 polish -j 30 unpolished.4.bam merged.subreadset.xml polished.4.bam
isoseq3 polish -j 30 unpolished.5.bam merged.subreadset.xml polished.5.bam
isoseq3 polish -j 30 unpolished.6.bam merged.subreadset.xml polished.6.bam
isoseq3 polish -j 30 unpolished.7.bam merged.subreadset.xml polished.7.bam
isoseq3 polish -j 30 unpolished.8.bam merged.subreadset.xml polished.8.bam
isoseq3 polish -j 30 unpolished.9.bam merged.subreadset.xml polished.9.bam
isoseq3 polish -j 30 unpolished.10.bam merged.subreadset.xml polished.10.bam
isoseq3 polish -j 30 unpolished.11.bam merged.subreadset.xml polished.11.bam

cat polished.*.hq.fasta > total_polished_hq.fasta

## Remove redundant
minimap2 -t 38 -ax splice -uf --secondary=no -C5 -O6,24 -B4 /mnt/data1/refgenome/ucsc_mRatBN7_2/rn7.fa total_polished_hq.fasta > total_polished_hq.sam

collapse_isoforms_by_sam.py --input total_polished_hq.fasta  -s total_polished_hq.sam  -o total_polished_hq_c0.99_i0.85_5_5000_3_diff50 --max_3_diff 50 --max_5_diff 5000 --gen_mol_count -c 0.99 -i 0.85

## we also mantain isoforms which got at least 5 FLNC reads support.

gunzip -c reads.fq.gz [reads2.fq.gz ...] | \
    awk 'NR % 4 == 2' | \
    sort | \
    tr NT TN | \
    ropebwt2 -LR | \
    tr NT TN | \
    fmlrc2-convert comp_msbwt.npy


fmlrc2 -p 20 -e 400 angii_comp_msbwt.npy angii_merged_flnc.fa angii_corrected_final.fa
fmlrc2 -p 20 -e 400 control_comp_msbwt.npy control_merged_flnc.fa control_corrected_final.fa

cat angii_corrected_final.fa control_corrected_final.fa >total_corrected_final.fa

minimap2 -t 38 -ax splice -uf --secondary=no -C5 -O6,24 -B4 /mnt/data1/refgenome/ucsc_mRatBN7_2/rn7.fa total_corrected_final.fa > total_corrected_final.sam

collapse_isoforms_by_sam.py --input total_corrected_final.fa  -s total_corrected_final.sam  -o total_flnc_c0.99_i0.85_5_5000_3_diff50 --max_3_diff 50 --max_5_diff 5000 --gen_mol_count -c 0.99 -i 0.85