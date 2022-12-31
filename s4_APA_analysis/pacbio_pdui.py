import numpy as np
import os
import sys
import datetime
import threading
import scipy as sp
import scipy.stats
from multiprocessing import Pool
from bisect import bisect


import math
import time

import multiprocessing

check_point_dict = {}
with open(sys.argv[1],'r') as f:
    for line in f:
        raw = line.rstrip().split()
        check_point_dict.setdefault(raw[0],[]).append(int(raw[1]))

def load_sequencing_depth(depth_file):
    seq_depth_list = []
    for line in open(depth_file, 'r'):
        fields = line.strip('\n').split('\t')
        seq_depth_list.append(int(fields[-1]))

    return np.array(seq_depth_list)


All_samples_sequencing_depths = load_sequencing_depth(sys.argv[2])



All_sample_coverage_weights = All_samples_sequencing_depths/np.mean(All_samples_sequencing_depths)


def load_wig_files(wig_file):
    Aligned_Wig_files=[]
    for line in open(wig_file,'r'):
        line = line.rstrip().split()[0]
        Aligned_Wig_files.append(line)
    return Aligned_Wig_files


All_Wig_files = load_wig_files(sys.argv[3])


def Load_Target_Wig_files_Multiple_threads_shared_dict_sampleid_key(All_Wig_files,UTR_Annotation_file,curr_processing_chr):
    num_samples = len(All_Wig_files)
    UTR_events_dict = {}
    for line in open(UTR_Annotation_file, 'r'):
        fields = line.strip('\n').split('\t')
        curr_chr = fields[0]
        if curr_chr == curr_processing_chr:
            region_start = fields[1]
            region_end = fields[2]

            curr_strand = fields[-1]
            UTR_pos = "%s:%s-%s" %(curr_chr, region_start, region_end)
            end_shift = int(round(abs(int(region_start) - int(region_end)) * 0.2))
            if curr_strand == "+":
                region_end = str(int(region_end) - end_shift)
            else:
                region_start = str(int(region_start) + end_shift)
            region_start = int(region_start) + 1
            region_end = int(region_end) - 1
            if region_start + 50 < region_end:
                UTR_events_dict[fields[3]] = [fields[0],region_start,region_end,fields[-1],UTR_pos]

    All_samples_extracted_3UTR_coverage_dict = {} # create only 1 dict
    
    for i in range(num_samples):
        curr_wig_file = All_Wig_files[i]
        curr_sample_All_chroms_coverage_dict = {}
        with open(curr_wig_file, 'r') as fin:
            for line in fin:
                if line[0] != '#' and line[0] != 't':
                    fields = line.strip('\n').split('\t')
                    chrom_name = fields[0]
                    if chrom_name == curr_processing_chr:
                        region_start = int(fields[1])
                        region_end = int(fields[2])


                        if chrom_name not in curr_sample_All_chroms_coverage_dict:
                            curr_sample_All_chroms_coverage_dict[chrom_name] = [[0],[0]]
                        if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]:
                            curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start)
                            curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
                        curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)
                        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(int(float(fields[-1])))
                    else:
                        if len(curr_sample_All_chroms_coverage_dict)>0:
                            break
            fin.close()
        if curr_processing_chr not in curr_sample_All_chroms_coverage_dict:
            print('no wig: ' + curr_wig_file, file=sys.stderr)
        else:
            curr_sample_All_chroms_coverage_dict[curr_processing_chr][1].append(0)

        curr_sample_coverage_dict = {}

        for curr_3UTR_event_id in UTR_events_dict:
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr_local = curr_3UTR_structure[0]
            if curr_chr_local in curr_sample_All_chroms_coverage_dict:
                curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr_local]
                region_start = curr_3UTR_structure[1]
                region_end = curr_3UTR_structure[2]
                left_region_index = bisect(curr_chr_coverage[0],region_start)
                right_region_index = bisect(curr_chr_coverage[0],region_end)

                extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                extracted_3UTR_region.insert(0,region_start)
                extracted_3UTR_region.append(region_end)

                curr_event_info = [extracted_coverage,extracted_3UTR_region]
                All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id,i] = curr_event_info
                    
    return All_samples_extracted_3UTR_coverage_dict, UTR_events_dict
    

Annotated_3UTR_file = sys.argv[4]
curr_processing_chr = sys.argv[4].split('_')[0]


All_samples_Target_3UTR_coverages, UTR_events_dict = Load_Target_Wig_files_Multiple_threads_shared_dict_sampleid_key(All_Wig_files, Annotated_3UTR_file,curr_processing_chr)

All_events_ids = list(UTR_events_dict.keys())

def Convert_wig_into_bp_coverage(extracted_coverage,extracted_3UTR_region,strand_info):
    bp_coverage = np.zeros(extracted_3UTR_region[-1] - extracted_3UTR_region[0])
    relative_start = extracted_3UTR_region[0]
    for i in range(len(extracted_coverage)):
        curr_region_start = extracted_3UTR_region[i] - relative_start
        curr_region_end = extracted_3UTR_region[i+1] - relative_start
        bp_coverage[curr_region_start:curr_region_end] = extracted_coverage[i]
    if strand_info == '-':
        bp_coverage = bp_coverage[::-1]

    return bp_coverage


def Estimation_abundance(Region_Coverage, break_point):
    Long_UTR_abun = np.mean(Region_Coverage[break_point:])
    Short_UTR_abun = np.mean(Region_Coverage[0:break_point] - Long_UTR_abun)
    if Short_UTR_abun < 0:
        Short_UTR_abun = 0
    Coverage_diff = Region_Coverage[0:break_point] - Long_UTR_abun - Short_UTR_abun
    Coverage_diff= np.append(Coverage_diff, Region_Coverage[break_point:] - Long_UTR_abun)
    Mean_Squared_error = np.mean(Coverage_diff**2)

    return Mean_Squared_error, Long_UTR_abun, Short_UTR_abun


result_dict = {}
for curr_3UTR_id in All_events_ids:
    num_samples = 6
    for check_point in check_point_dict[curr_3UTR_id]:
        curr_3UTR_structure = UTR_events_dict[curr_3UTR_id]
        region_start = curr_3UTR_structure[1]
        region_end = curr_3UTR_structure[2]
        curr_strand = curr_3UTR_structure[-2]
        UTR_pos = curr_3UTR_structure[-1]
        curr_3UTR_all_samples_bp_coverage = []
        
        for i in range(num_samples):
            curr_sample_curr_3UTR_coverage_wig = All_samples_Target_3UTR_coverages[curr_3UTR_id, i]
            curr_3UTR_curr_sample_bp_coverage = Convert_wig_into_bp_coverage(curr_sample_curr_3UTR_coverage_wig[0], curr_sample_curr_3UTR_coverage_wig[1], curr_strand)
            curr_3UTR_all_samples_bp_coverage.append(curr_3UTR_curr_sample_bp_coverage)
            
        Region_Coverages = []
        Pass_threshold_index = []
        
        for i in range(num_samples):
            curr_Region_Coverage_raw = curr_3UTR_all_samples_bp_coverage[i]
            curr_Region_Coverage = curr_Region_Coverage_raw/All_sample_coverage_weights[i]
            
            curr_first_100_coverage = np.mean(curr_Region_Coverage_raw[0:99])
            if curr_first_100_coverage > 0.1:
                Pass_threshold_index.append(i)
                Region_Coverages.append(curr_Region_Coverage)
            
            Mean_squared_error_list = []
            Estimated_3UTR_abundance_list = []
            
            if curr_strand == '-':
                curr_search_point = region_end - check_point
            else:
                curr_search_point = check_point - region_start
            
            All_samples_result = [[],[],[]]

            for curr_sample_region_coverage in Region_Coverages:
                Mean_Squared_error, Long_UTR_abun, Short_UTR_abun = Estimation_abundance(curr_sample_region_coverage, curr_search_point)
                All_samples_result[0].append(Mean_Squared_error)
                All_samples_result[1].append(Long_UTR_abun)
                All_samples_result[2].append(Short_UTR_abun)
            Mean_Squared_error = np.mean(np.array(All_samples_result[0]))
            Mean_squared_error_list.append(Mean_Squared_error)
            Estimated_3UTR_abundance_list.append([All_samples_result[1],All_samples_result[2]])    
        
        min_ele_index = Mean_squared_error_list.index(min(Mean_squared_error_list))
            

        select_mean_squared_error = Mean_squared_error_list[min_ele_index]

        UTR_abundances = [['NA']*num_samples, ['NA']*num_samples]
        UTR_abundances_passed = Estimated_3UTR_abundance_list[min_ele_index]
        for k in range(len(Pass_threshold_index)):
            UTR_abundances[0][Pass_threshold_index[k]] = UTR_abundances_passed[0][k]
            UTR_abundances[1][Pass_threshold_index[k]] = UTR_abundances_passed[1][k]
            
        try:
            group1_coverages = np.zeros(2)
            group2_coverages = np.zeros(2)
    
            group1_coverages[0] =   UTR_abundances[0][0] + UTR_abundances[0][1]
            group1_coverages[1] =  UTR_abundances[1][0] + UTR_abundances[1][1]
            
            group2_coverages[0] =   UTR_abundances[0][2] + UTR_abundances[0][3]
            group2_coverages[1] =   UTR_abundances[1][2] + UTR_abundances[1][3]
        
            ratio_val,P_val = sp.stats.fisher_exact([group1_coverages/2,group2_coverages/2])
        except:
            ratio_val,P_val = 'Na','Na'
            
        if str(select_mean_squared_error) != "Na":
            num_non_zero = 1
            if num_non_zero > 0:
                All_long_inclusion_ratios = []

                for i in range(num_samples):
                    if UTR_abundances[0][i] != 'NA':
                        # long 3'UTR percentage
                        curr_sample_ratio = float(UTR_abundances[0][i])/(float(UTR_abundances[0][i]) + float(UTR_abundances[1][i]))
                        All_long_inclusion_ratios.append(curr_sample_ratio)
        try:
            a1,a2,a3,a4 = All_long_inclusion_ratios[0],All_long_inclusion_ratios[1],All_long_inclusion_ratios[2],All_long_inclusion_ratios[3]            
	       pdui_diff = np.mean([a3,a4]) - np.mean([a1,a2])
            result_dict.setdefault(curr_3UTR_id,[]).append([pdui_diff,check_point,All_long_inclusion_ratios[0],All_long_inclusion_ratios[1],All_long_inclusion_ratios[2],All_long_inclusion_ratios[3],ratio_val,P_val])
        except:
            pass
        
for k,v in result_dict.items():
    for i in v:
        print(k,'\t' ,i)