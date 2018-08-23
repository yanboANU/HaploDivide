#########################################################################
# File Name: step4.sh
# Author: ma6174
# mail: ma6174@163.com
# Created Time: Thu Aug 16 21:05:14 2018
#########################################################################
#!/bin/bash

#/home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first_10M_NIST/
nohup python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/calc_sensitive/calc_senstive_specifity_ref.py /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first_10M_NIST/SNP_record 1_snp_mutation_* 0 *_*_cov >snp_result &

######delete
nohup python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/calc_sensitive/calc_senstive_specifity_ref_delete.py ~/bio/Project/Guofei/reference/GRCh37_hg19/first_10M_NIST/Delete_record 1_delete_* 0 /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first*.fa *_cov $1 >delete_nofilter.log &


#####insert
nohup python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/calc_sensitive/calc_senstive_specifity_ref_insert.py ~/bio/Project/Guofei/reference/GRCh37_hg19/first_10M_NIST/Insert_record 1_insert_* 0 /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first*.fa *_cov $1 >insert_nofilter.log &


#/home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/chr1_multiple_pos
python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/indel/remove_delete_in_multiple_range.py ~/bio/Project/Guofei/reference/GRCh37_hg19/chr1_multiple_pos 1_delete_* 10 > remove.delete.log

python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/calc_sensitive/calc_senstive_specifity_ref_delete.py ~/bio/Project/Guofei/reference/GRCh37_hg19/first_10M_NIST/filter_Delete_record filter_in_multiple_range_delete 0 /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first*.fa *_cov $1 >delete_filter.log


python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/indel/remove_insert_in_multiple_range.py ~/bio/Project/Guofei/reference/GRCh37_hg19/chr1_multiple_pos 1_insert_* 10 > remove.insert.log

python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/calc_sensitive/calc_senstive_specifity_ref_insert.py ~/bio/Project/Guofei/reference/GRCh37_hg19/first_10M_NIST/filter_Insert_record filter_in_multiple_range_insert 0 /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first*.fa *_cov $1 >insert_filter.log

#get distribution
#python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/get_distribution/get_delete_rate_distribution.py /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first_10M_NIST/Delete_record 1_delete_0_10000000 10 /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first*.fa *_cov $1 >get_delete_distribution.log

#python /home/yulin/liyanbo/script/HaploDivide/HaploCheck/get_distribution/get_insert_rate_distribution.py /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first_10M_NIST/Insert_record 1_insert_0_10000000 10 /home/yanbo/bio/Project/Guofei/reference/GRCh37_hg19/first*.fa *_cov $1 >get_insert_distribution.log

#python3 /home/yulin/liyanbo/script/HaploDivide/HaploCheck/select_indel_slope.py TP_delete_rate_distribution_30 FP_delete_rate_distribution_30 >best_delete_30

#python3 /home/yulin/liyanbo/script/HaploDivide/HaploCheck/select_indel_slope.py TP_insert_rate_distribution_30 FP_insert_rate_distribution_30 >best_insert_30
