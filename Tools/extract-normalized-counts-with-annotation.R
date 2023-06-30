library(tidyverse)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rmarkdown)
library(ChIPseeker)
library(dplyr)

# Read raw feature counts data
cts <- read.csv("tko_cmt_pbs_wt_featureCounts_for_deseq2.txt",sep="\t")
head(cts)

tkopbscmst<-read.csv("TKO_PBS_vs_TKO_Cmst_sig_normalized_data.csv")
wt_pbs_cmst<-read.csv("WT_PBS_vs_WT_Cmst_sig_normalized_data.csv")
cmst_wt_tko<-read.csv("WT_Cmst_vs_TKO_Cmst_sig_normalized_data.csv")

head(tkopbscmst)

# Set the column header
names(tkopbscmst)[1] <- "peakName"
names(wt_pbs_cmst)[1] <- "peakName"
names(cmst_wt_tko)[1] <- "peakName"
names(cts)[1] <- "peakName"

head(tkopbscmst)
head(cts)

# Split first column to coordinate parts
tkopbscmst <- tkopbscmst %>% separate(peakName,c("Chr","Start","End"))
wt_pbs_cmst <- wt_pbs_cmst %>% separate(peakName,c("Chr","Start","End"))
cmst_wt_tko <- cmst_wt_tko %>% separate(peakName,c("Chr","Start","End"))
cts <- cts %>% separate(peakName,c("Chr","Start","End"))

head(tkopbscmst)
head(cts)

# Add peakname column
res1 = mutate(tkopbscmst, concated_column = paste(Chr, Start, End, sep = '.'))
res2 = mutate(wt_pbs_cmst, concated_column = paste(Chr, Start, End, sep = '_'))
res3 = mutate(cmst_wt_tko, concated_column = paste(Chr, Start, End, sep = '_'))
res4 = mutate(cts, concated_column = paste(Chr, Start, End, sep = '_'))

head(res4)

res1 = res1 %>% relocate(concated_column, .after = End)
res2 = res2 %>% relocate(concated_column, .after = End)
res3 = res3 %>% relocate(concated_column, .after = End)
res4 = res4 %>% relocate(concated_column, .after = End)

head(res4)

# Create granges object
myRanges1<- GRanges(seqnames=res1$Chr,ranges=IRanges(as.numeric(res1$Start),as.numeric(res1$End)),strand=NULL,res1$concated_column,res1$TKO_Cmt1, res1$TKO_Cmt2, res1$TKO_PBS1, res1$TKO_PBS2, res1$WT_Cmst1, res1$WT_Cmst2, res1$WT_PBS1, res1$WT_PBS2)
myRanges2<- GRanges(seqnames=res2$Chr,ranges=IRanges(as.numeric(res2$Start),as.numeric(res2$End)),strand=NULL,res2$concated_column,res2$TKO_Cmt1, res2$TKO_Cmt2, res2$TKO_PBS1, res2$TKO_PBS2, res2$WT_Cmst1, res2$WT_Cmst2, res2$WT_PBS1, res2$WT_PBS2)
myRanges3<- GRanges(seqnames=res3$Chr,ranges=IRanges(as.numeric(res3$Start),as.numeric(res3$End)),strand=NULL,res3$concated_column,res3$TKO_Cmt1, res3$TKO_Cmt2, res3$TKO_PBS1, res3$TKO_PBS2, res3$WT_Cmst1, res3$WT_Cmst2, res3$WT_PBS1, res3$WT_PBS2)
myRanges4<- GRanges(seqnames=res4$Chr,ranges=IRanges(as.numeric(res4$Start),as.numeric(res4$End)),strand=NULL,res4$concated_column,res4$TKO_Cmt1, res4$TKO_Cmt2, res4$TKO_PBS1, res4$TKO_PBS2, res4$WT_Cmst1, res4$WT_Cmst2, res4$WT_PBS1, res4$WT_PBS2)


# Annotate OCR's
tkopbscmst_annot<-annotatePeak(myRanges1, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
wt_pbs_cmst_annot<-annotatePeak(myRanges2, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
cmst_wt_tko_annot<-annotatePeak(myRanges3, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
cts_annot<-annotatePeak(myRanges4, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
summary(cts_annot)
cts_annot

# Convert to data frame
tkopbscmst_annot_df<-as.data.frame(tkopbscmst_annot)
wt_pbs_cmst_annot_df<-as.data.frame(wt_pbs_cmst_annot)
cmst_wt_tko_annot_df<-as.data.frame(cmst_wt_tko_annot)
cts_annot_df<-as.data.frame(cts_annot)

head(tkopbscmst_annot_df)
head(cts_annot_df)
View(head(cts_annot_df))

tkopbscmst_annot_dfo=tkopbscmst_annot_df %>% rename(concated_column=res1.concated_column,TKO_Cmt1=res1.TKO_Cmt1, TKO_Cmt2=res1.TKO_Cmt2, TKO_PBS1=res1.TKO_PBS1,TKO_PBS2=res1.TKO_PBS2,WT_Cmst1=res1.WT_Cmst1,WT_Cmst2=res1.WT_Cmst2,WT_PBS1=res1.WT_PBS1,WT_PBS2=res1.WT_PBS2)
wt_pbs_cmst_annot_dfo=wt_pbs_cmst_annot_df %>% rename(concated_column=res2.concated_column,TKO_Cmt1=res2.TKO_Cmt1, TKO_Cmt2=res2.TKO_Cmt2, TKO_PBS1=res2.TKO_PBS1,TKO_PBS2=res2.TKO_PBS2,WT_Cmst1=res2.WT_Cmst1,WT_Cmst2=res2.WT_Cmst2,WT_PBS1=res2.WT_PBS1,WT_PBS2=res2.WT_PBS2)
cmst_wt_tko_annot_dfo=cmst_wt_tko_annot_df %>% rename(concated_column=res3.concated_column,TKO_Cmt1=res3.TKO_Cmt1, TKO_Cmt2=res3.TKO_Cmt2, TKO_PBS1=res3.TKO_PBS1,TKO_PBS2=res3.TKO_PBS2,WT_Cmst1=res3.WT_Cmst1,WT_Cmst2=res3.WT_Cmst2,WT_PBS1=res3.WT_PBS1,WT_PBS2=res3.WT_PBS2)

# Covert symbol column to all upper case
cts_annot_df=cbind(data.frame(toupper(cts_annot_df$SYMBOL)),cts_annot_df)
colnames(cts_annot_df)

# Rename column headers
cts_annot_dfo=cts_annot_df %>% dplyr::rename(Symbol=toupper.cts_annot_df.SYMBOL., concated_column=res4.concated_column,TKO_Cmt1=res4.TKO_Cmt1, TKO_Cmt2=res4.TKO_Cmt2, TKO_PBS1=res4.TKO_PBS1,TKO_PBS2=res4.TKO_PBS2,WT_Cmst1=res4.WT_Cmst1,WT_Cmst2=res4.WT_Cmst2,WT_PBS1=res4.WT_PBS1,WT_PBS2=res4.WT_PBS2, SYMBOLlower=SYMBOL)
colnames(cts_annot_dfo)

# Select only gene symbol and data columns
cts_annot_selected= cts_annot_dfo %>% dplyr::select(Symbol,TKO_Cmt1,TKO_Cmt2,TKO_PBS1,TKO_PBS2,WT_Cmst1,WT_Cmst2,WT_PBS1,WT_PBS2)
colnames(cts_annot_selected)

# Remove rows without gene names
cts_annot_selected=cts_annot_selected %>% drop_na(Symbol)

# Group duplicate genes and sum column
head(cts_annot_selected) %>% group_by(Symbol) %>% summarise_at(c("TKO_Cmt1","TKO_Cmt2"),sum, na.rm = TRUE)
cts_duplicate_summed=cts_annot_selected %>% group_by(Symbol) %>% summarise_at(c("TKO_Cmt1","TKO_Cmt2","TKO_PBS1","TKO_PBS2","WT_Cmst1","WT_Cmst2","WT_PBS1","WT_PBS2"),sum)
View(cts_duplicate_summed)

# Write output to file
write.csv(wt_pbs_cmst_annot_dfo, file="WT_PBS_vs_WT_Cmst_sig_normalized_data_annotated.csv")
write.csv(cmst_wt_tko_annot_dfo, file="WT_Cmst_vs_TKO_Cmst_sig_normalized_data_annotated.csv")
write.table(cts_duplicate_summed, file="tko_cmt_pbs_wt_featureCounts_annotated_duplicate_summed.txt", quote = FALSE, row.names = FALSE, sep = "\t")

