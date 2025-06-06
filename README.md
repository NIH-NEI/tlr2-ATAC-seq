# ATAC-Seq Analysis of TLR2_Cmast_WT_PBS

## Data Analysis Tools/Methods

1. [atac-fastp.swarm](Tools/atac-fastp.swarm)  : swarm commands to run fastp to process fastq data, to remove adapters
2. [bowtie2.swarm](Tools/bowtie2.swarm) : swarm commands to run bowtie alignments
3. [chrMremovedBam.swarm](Tools/chrMremovedBam.swarm) : swarm commands to remove chrM data from alignments
4. [bamCoverageCPMNormalization.swarm](Tools/bamCoverageCPMNormalization.swarm) : swarm commands to normalize alignment data
5. [genrichPeakCalling.swarm](Tools/genrichPeakCalling.swarm) : swarm commands to call peaks
6. [genrich_stats.sh](Tools/genrich_stats.sh) : shell script that generates peak statistics
7. [featurecountscript.sh](Tools/featurecountscript.sh) : shell script to do feature counting (count matrix)
8. [deseq2.Rmd](Tools/deseq2.Rmd) : R script to run differential peak analysis
9. [extract-normalized-counts-with-annotation.R](Tools/extract-normalized-counts-with-annotation.R) : R script to generate normalized counts with annotation
10. [tko_cmst_peakannotation.Rmd](Tools/tko_cmst_peakannotation.Rmd) : R script to generate peak annotations
11. [sample_correlation.swarm](Tools/sample_correlation.swarm) : swarm command to generate sample correlation heatmap
12. [enrichment.Rmd](Tools/enrichment.Rmd) : R script to run enrichment analysis
