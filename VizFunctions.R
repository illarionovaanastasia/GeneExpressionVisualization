library(readxl)
library(GenomicFeatures)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(readr)
library(data.table)
library(stringr)
library(plotly)
library(viridis)
library(hrbrthemes)
library(ggbreak) 
library(patchwork)
library(OUTRIDER)
library(forcats)
library(ggforce)
library(ggtranscript)
library(ComplexHeatmap)


# Required dataframes:
# cts - normilised counts from bRNAseq, scRNAseq, ATACseq, methylation profile assays
# vcf_gt - df with 0,1,2 genotypes and SVs (SV column with SV ids required)
# sample_map - lookup df if different sample ids were used for different assays
# homer - df with ATACseq/methylation regions annotations with homer

# _______ FUNCTIONS _____BLOCKS ____________________________________________________

# Function to plot gene expression from bRNAseq

plot_gt_expr = function(gts=vcf_gt, cts, smap=sample_map, SV, gene_id, gene_name){
        SV = as.character(SV)
        df = gts[gts$SV == SV,]
        print(paste0("bRNAseq: Plotting for SV ", SV))
        df = as.data.frame(t(df))
        df$CDI = rownames(df)
        df = merge(df, smap[,2:3], by.x = 2, by.y =2)
        colnames(df)[2] = "GT"
        #df$GT[df$GT == "0"] = "0/0"
        #df$GT[df$GT == "1"] = "+"
        #df$GT[df$GT == "2"] = "+"
        sr_df = cts[cts$...1 %like% gene_id,]
        if (dim(sr_df)[1] > 0){
                sr_df = as.data.frame(t(sr_df))
                sr_df$PPMI = rownames(sr_df)
                sr_df = sr_df[-1,]
                df = merge(df, sr_df, by.x = 3, by.y = 2)
                colnames(df)[4] = "norm_counts_SR"
                df$norm_counts_SR = as.numeric(df$norm_counts_SR)
                
                p = ggplot(df, mapping=aes(x=GT, y=norm_counts_SR)) +
                        geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
                        geom_boxplot(width=0.1) + 
                        theme_minimal() + scale_fill_brewer(palette="RdBu") + theme_minimal()+xlab(gene_name)
                png(paste0("brna/", SV, "_", gene_name,".png"))
                print(p)
                dev.off()
        } else {print(paste0(ensg, " was not detected"))}
        
}

# Function to plot gene expression from scRNAseq

plot_gt_expr_sc = function(gts=vcf_gt, cts, SV, gene_name){
        SV = as.character(SV)
        df = gts[gts$SV == SV,]
        print(paste0("scRNA: Plotting for SV ", SV))
        df = as.data.frame(t(df))
        df$CDI = rownames(df)
        colnames(df)[1] = "GT"
        #df$GT[df$GT == "0"] = "0/0"
        #df$GT[df$GT == "1"] = "+"
        #df$GT[df$GT == "2"] = "+"
        sr_df = cts[cts$gene == gene_name,]
        if (dim(sr_df)[1] > 0){
                sr_df = as.data.frame(t(sr_df))
                sr_df$PPMI = rownames(sr_df)
                sr_df = sr_df[-1,]
                df = merge(df, sr_df, by.x = 2, by.y = 2)
                colnames(df)[3] = "norm_counts"
                df$norm_counts = as.numeric(df$norm_counts)
                
                p = ggplot(df, mapping=aes(x=GT, y=norm_counts)) +
                        geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
                        geom_boxplot(width=0.1) + 
                        theme_minimal() + scale_fill_brewer(palette="RdBu") + theme_minimal()+xlab(gene_name)
                png(paste0("scrna/", SV, "_", gene_name,"_sc_DA.png"))
                print(p)
                dev.off()
        } else {print(paste0(gene_name, " was not detected"))}
        
}

# Function to plot ATACseq profile

plot_gt_atac_expr = function(gts=vcf_gt, cts, homer, SV, gene_id, gene_name, an = F_annot){
        print(paste0("ATACseq: Plotting for SV ", SV))
        mh_sb = subset(homer, homer$`Entrez ID` == gene_id)
        if (dim(mh_sb)[1] > 0){
                mh_sb = mh_sb[,c(1,8,3)]
                mcts_sb = subset(cts,cts$Region %in% mh_sb$PeackID)
                mcts_sbt = as.data.frame(t(mcts_sb))
                regions = mcts_sbt[1,]
                ppmi = rownames(mcts_sbt)[-1]
                mcts_sbt = as.data.frame(mcts_sbt[-1,])
                mcts_sbt$PPMI = ppmi 
                colnames(mcts_sbt)[1:(length(colnames(mcts_sbt))-1)] = regions
                mcts_sbt_long <- gather(mcts_sbt, Region, Norm_counts, colnames(mcts_sbt)[1]:colnames(mcts_sbt)[length(colnames(mcts_sbt))-1], factor_key=TRUE)
                mcts_sbt_long = merge(mcts_sbt_long, mh_sb, by.x = 2, by.y = 1)
                df = gts[gts$SV == SV,]
                df = as.data.frame(t(df))
                df$CDI = rownames(df)
                df = df[-1,]
                df = merge(df, an[,1:2], by.x = 2, by.y =1)
                colnames(df)[2] = "GT"
                mcts_sbt_long = merge(mcts_sbt_long, df, by.x = 2, by.y = 3)
                colnames(mcts_sbt_long)[5] = "Coordinate"
                mcts_sbt_long$Norm_counts = as.numeric(mcts_sbt_long$Norm_counts)
                mcts_sbt_long$Coordinate = as.factor(mcts_sbt_long$Coordinate)
                mcts_sbt_long$GT[mcts_sbt_long$GT == "0"] = "0/0"
                mcts_sbt_long$GT[mcts_sbt_long$GT != "0/0"] = "+"
                mcts_sbt_long$Annotation = as.factor(mcts_sbt_long$Annotation)
                p =ggplot(mcts_sbt_long, aes(x=Coordinate, y=Norm_counts, fill=GT)) +
                        geom_boxplot(position=position_dodge(1)) + 
                        theme_minimal()+xlab(paste0(gene_name))
                png(paste0("atac/", SV, "_", gene_name,"_atac_65.png"),width = 1016, height = 700)
                print(p)
                dev.off()  
        } else { (print(paste0(gene_name, " not in the ATACseq array")))}
        
        
}

# Function to plot methylation profile

plot_gt_meth_expr = function(gts=vcf_gt, cts, homer, SV, gene_id, gene_name, an = F_annot){
        print(paste0("METH: Plotting for SV ", SV))
        mh_sb = subset(homer, homer$`Entrez ID` == gene_id)
        if (dim(mh_sb)[1] > 0){
                mh_sb = mh_sb[,c(1,8,3)]
                mcts_sb = subset(cts,cts$ID %in% mh_sb$PeackID)
                mcts_sbt = as.data.frame(t(mcts_sb))
                regions = mcts_sbt[1,]
                ppmi = rownames(mcts_sbt)[-1]
                mcts_sbt = as.data.frame(mcts_sbt[-1,])
                mcts_sbt$PPMI = ppmi 
                colnames(mcts_sbt)[1:(length(colnames(mcts_sbt))-1)] = regions
                mcts_sbt_long <- gather(mcts_sbt, Region, Norm_counts, colnames(mcts_sbt)[1]:colnames(mcts_sbt)[length(colnames(mcts_sbt))-1], factor_key=TRUE)
                mcts_sbt_long = merge(mcts_sbt_long, mh_sb, by.x = 2, by.y = 1)
                df = gts[gts$SV == SV,]
                df = as.data.frame(t(df))
                df$CDI = rownames(df)
                df = df[-1,]
                df = merge(df, an[,1:2], by.x = 2, by.y =1)
                colnames(df)[2] = "GT"
                mcts_sbt_long = merge(mcts_sbt_long, df, by.x = 2, by.y = 3)
                colnames(mcts_sbt_long)[5] = "Coordinate"
                mcts_sbt_long$Norm_counts = as.numeric(mcts_sbt_long$Norm_counts)
                mcts_sbt_long$Coordinate = as.factor(mcts_sbt_long$Coordinate)
                mcts_sbt_long$GT[mcts_sbt_long$GT == "0"] = "0/0"
                mcts_sbt_long$GT[mcts_sbt_long$GT != "0/0"] = "+"
                mcts_sbt_long$Annotation = as.factor(mcts_sbt_long$Annotation)
                p =ggplot(mcts_sbt_long, aes(x=Coordinate, y=Norm_counts, fill=GT)) +
                        geom_boxplot(position=position_dodge(1)) + theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
                        xlab(paste0(gene_name))
                png(paste0("meth/", SV, "_", gene_name,"_meth_65.png"),width = 1016, height = 700)
                print(p)
                dev.off()   
        } else { (print(paste0(gene_name, " not in the methylation array")))}
        
}

# Function to plot transcript structure and expression

plot_diu = function(gts = vcf_gt, cts = cts_tr_long, SV, gene_name, gene_id, an = annot){
        print(paste0("DIU: Plotting for SV ", SV))
        gt_sv = gts[gts$SV == SV,]
        gt_sv = data.frame(Samples = colnames(gt_sv), GT = t(gt_sv[1,]))
        cts_gt_sv = merge(cts, gt_sv, by.x = 6, by.y = 1)
        cts_gt_sv$GT = as.factor(cts_gt_sv$GT)
        #cts_gt_sv = subset(cts_gt_sv, cts_gt_sv$norm_counts > 0)
        exons <- subset(an, an$type == "exon" & an$gene_name %like% gene_id)
        #exons = subset(exons, exons$transcript_name %in% cts_gt_sv$transcript)
        exp = ggplot(data=cts_gt_sv[cts_gt_sv$gene == gene_id,], aes(x=transcript, y=log10(norm_counts), fill = GT)) +
                coord_flip()+geom_boxplot() + theme_light() + xlab("")+theme(axis.title.y=element_blank(),
                                                                             axis.text.y=element_blank(),
                                                                             axis.ticks.y=element_blank())
        struc = ggplot(data = exons, aes(
                xstart = start,
                xend = end,
                y = transcript_name)) +
                geom_intron(
                        data = to_intron(exons, "transcript_name"),
                        aes(strand = strand)
                ) + xlab("bp")
        p = ggarrange(struc, exp, ncol = 2, nrow = 1, legend="right", align='v')
        
        png(paste0("diu/", SV, "_", gene_name,"_diu.png"),width=950, height=740, res=100)
        print(p)
        dev.off()
}

# _________________WRAP-UP  FUNCTION _________________________________________________


plot_assays = function(sv, geneid, genename){
        print(paste0("Started analysis for ", sv))
        plot_gt_expr(SV = sv, gene_name = genename, gene_id = geneid)
        plot_gt_expr_sc(SV = sv, gene_name = genename)
        plot_gt_atac_expr(SV = sv, gene_name = genename, gene_id = geneid)
        plot_gt_meth_expr(SV = sv, gene_name = genename, gene_id = geneid)
        plot_diu(SV = sv, gene_name = genename, gene_id = geneid)
}

# _________________  FUNCTION METADATA ____________________________________

get_metadata = function(SV, meta = sample_metadata, gts = vcf_gt){
        SV = as.character(SV)
        df = gts[gts$SV == SV,]
        df = as.data.frame(t(df))
        df$CDI = rownames(df)
        meta_df = merge(df, meta, by.x = 2, by.y = 1)
        return(meta_df)
}
