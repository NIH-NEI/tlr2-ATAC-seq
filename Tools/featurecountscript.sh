module load multiqc
module load bedtools
module load subread

## This script prepares feature counts for differential expression analysis using deseq2
## Tested using - sinteractive --cpus-per-task=64 --mem=64g --gres=lscratch:300 

## Concatenate (merge) genrich called peaks files
#cat *.narrowPeak > cat_merged.narrowPeak

## Merged peaks file sorted by chr column location start column and piped to bedtools merge
## awk converts the bedtools merge output to saf format for featureCounts use, generating the first column combining chr-start-stop-peakid
bedtools merge -i <(sort -k1,1V -k2,2n cat_merged.narrowPeak) | awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' > tko_cmst_pbs_wt_merged.peaks.saf

## featureCounts with paired end -p option and SAF format -F options, used -T 64 threads
featureCounts -p -a tko_cmst_pbs_wt_merged.peaks.saf -F SAF -o tko_cmt_pbs_wt_featureCounts.txt -T 64 ../Alignment/*sorted.bam

## multiqc of featurecounts summary
multiqc tko_cmt_pbs_wt_featureCounts.txt.summary -n tko_cmt_pbs_wt_featureCounts_multiqc


