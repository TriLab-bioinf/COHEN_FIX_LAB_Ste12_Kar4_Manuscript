#!/bin/bash

#bed=../pooled_groups_no_ctrl/Ste12_WT_AF.txt_win24/Ste12_WT_AF.txt_win24_Ste12_WT_AF.events.curated_notin_1_1.txt_dupsmerged_241bp.bed
bed=../1_1_Kar4dep_Ste12_AF.txt_win24/1_1_Kar4dep_Ste12_AF.txt_win24_1_1_Kar4dep_Ste12_AF.events.curated.txt_dupsmerged_241bp.bed

b1=../../dataset_668-670_MRose/36203_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.proper.q30_top.bw
b2=../../dataset_668-670_MRose/36203_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.proper.q30_bottom.bw
b3=../../dataset_745-MRose/39270_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.proper.q30_top.bw
b4=../../dataset_745-MRose/39270_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.proper.q30_bottom.bw
b5=../../dataset_668-670_MRose/36204_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_-_XO_FilteredBAM.proper.q30_top.bw
b6=../../dataset_668-670_MRose/36204_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_-_XO_FilteredBAM.proper.q30_bottom.bw
b7=../../dataset_745-MRose/39271_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_-_XO_FilteredBAM.proper.q30_top.bw
b8=../../dataset_745-MRose/39271_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_-_XO_FilteredBAM.proper.q30_bottom.bw

#out=Ste12_WT_AF.txt_win24_Ste12_WT_AF.events.curated_notin_1_1.txt_dupsmerged_241bp.bed
out=1_1_Kar4dep_Ste12_AF.txt_win24_1_1_Kar4dep_Ste12_AF.events.curated.txt_dupsmerged_241bp.bed

~/.local/bin/multiBigwigSummary BED-file --BED $bed -b $b1 $b2 $b3 $b4 $b5 $b6 $b7 $b8 -o ${out}_counts.npz --outRawCounts ${out}_counts.txt


