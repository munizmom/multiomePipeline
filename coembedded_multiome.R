###########################################################################################
###########################################################################################
##  Script for multi-ome embedded analysis: setting up the pipeline using Seurat    	 ##
###########################################################################################
###########################################################################################
# README
# Developped and maintained by Mar Muniz Moreno.Last update July 2022.
# This scrip have a serie of functions that will allow to perform multiome analysis
# in the case where a co embedded multiome is needed due to different processing times 
# of the  samples, for more info about the tools, tutorial and main packages used:
# https://satijalab.org/signac/articles/pbmc_multiomic.html
# https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Multiome.html
# https://biofam.github.io/MOFA2/tutorials.html
# main pipeline joint rnaseq and ATAC: https://satijalab.org/seurat/articles/atacseq_integration_vignette.html
# multiome analysis:
# https://satijalab.org/signac/articles/pbmc_multiomic.html
# https://satijalab.org/seurat/articles/multimodal_vignette.html
# https://satijalab.org/signac/articles/pbmc_vignette.html
# https://satijalab.org/signac/articles/integrate_atac.html
# https://satijalab.org/seurat/articles/atacseq_integration_vignette.html#annotate-scatac-seq-cells-via-label-transfer-1
###########################################################################################
# NOTES: data 
###########################################################################################
#  R version 4.1.2 (2021-11-01) 
############################################################################################
# preliminary step: got the input matrix files from cluster:
############################################################################################
##################################### Packages needed ######################################
############################################################################################
#load libraries
library("Signac");library("Seurat");library("EnsDb.Mmusculus.v79"); library("ggplot2");
library("biovizBase"); library("hdf5r");library("dplyr");library("tidyr");library("gridExtra");
library("patchwork"); library("grid");library("gridBase"); library("cowplot");
library("BSgenome.Mmusculus.UCSC.mm10"); library("readxl");library("matrixStats");
library("VennDiagram");



#export R_LIBS=/software/anaconda3/envs/R/lib/R/library/ #cluster


#Before begining:
#rm(list =ls()) ## erasing all the enviroment variables
set.seed(22); # need to set it to get always the same random results and plots
#sessionInfo()

#wd:
#####################################################
wd<- "/home/Documents/omicsData/multiome/"
setwd(wd);

#project name and folder
##################################
name <- "project"
f_name<- paste0(name, "/");
nameSub <- "seurat"
##################################

#folder structure
#################################
#set input directories
f_input <- "input/"; 
#set output directories
f_results <- "/results/";
f_Rdata <- "RData/";
f_qc <- "QC/";
f_violin <-"violinPLot/";
f_quantiles <-"quantiles/";
f_cellcycle<- "cellcycle/";
f_pca<- "pca/";
f_preIntegration <-"PreIntegration/";
f_after_enbedding <-"After_embedding/";
f_ridgeplot <- "RidgePlot/";
f_foundClusters <- "FoundClusters/";
f_tables <- "tables/";
f_plots <- "plots/";
f_markers <- "markers/";
f_tssEnrichment <- "tssEnrichment/"
f_nucleosomeSignal <- "nucleosomeSignal/"
f_violin <-"violinPlot/"
f_DEA <- "DEA/"
f_DAA <- "DAA/"
f_atacCoverage <- "ATAC_coverage/"
f_venn <- "vennDiagram/"
f_vdata <- "vennTables/"
f_vplots <- "vennPlots/"


#################################

#creation of the fodler structure
dir.create(file.path(wd, f_results), showWarnings = F);
dir.create(file.path(wd, f_input), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration,f_tssEnrichment), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration,f_nucleosomeSignal), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration,f_violin), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration,f_ridgeplot), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration,f_quantiles), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration,f_cellcycle), showWarnings = F);
#
dir.create(file.path(wd, f_results,f_qc ,f_after_enbedding), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_after_enbedding,f_tables), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration,f_pca), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_after_enbedding,f_pca), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata ,f_after_enbedding), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata ,f_preIntegration), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_tables), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_plots), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers), showWarnings = F);

dir.create(file.path(wd, f_results,f_markers), showWarnings = F);
dir.create(file.path(wd, f_results,f_markers,f_plots), showWarnings = F);
dir.create(file.path(wd, f_results,f_markers,f_tables), showWarnings = F);

dir.create(file.path(wd, f_results,f_DEA ), showWarnings = F);
dir.create(file.path(wd, f_results,f_DEA ,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_DAA ), showWarnings = F);
dir.create(file.path(wd, f_results,f_DAA ,f_Rdata), showWarnings = F);

dir.create(file.path(wd, f_results,f_atacCoverage), showWarnings = F);
dir.create(file.path(wd, f_results,f_atacCoverage ,f_plots), showWarnings = F);

dir.create(file.path(wd, f_results,f_venn), showWarnings = F);
##############################################################################


###################################################################################################
# 1) STEP1: seurat creating RNASeq and ATACseq objects
###################################################################################################
#multiple samples to integrate in multiome analysis following pipeline recommended:
#https://github.com/satijalab/seurat/issues/5346
#https://github.com/satijalab/seurat/issues/4064


metafile <- data.frame(Sample=c("1","2","3","4"),Genotype=c("Mut","Mut","Control","Control"),
	Replicate=c(1,2,1,2),Condition=c("Mut_1","Mut_2","Control_1","Control_2"));
date <- "072622";

#rnaseq samples
rnaObject_generation.function <- function(metafile,EnsDb_annotation,nameOutput,min.cells,min.features,date){
	dataList <- list();
	dfResList <- list();
	#note: setting up min.cells = 3, min.features = 200 as in the seurat package vignette

	for (i in 1:length(metafile$Sample)){
		print(i)
		sample <- Read10X_h5(paste0(wd,f_input,metafile$Sample[i],"_filtered_feature_bc_matrix.h5"));
		metadata <- read.csv(file=paste0(wd,f_input,metafile$Sample[i],"_per_barcode_metrics.csv"), header = TRUE, row.names = 1);
		metadata <- metadata[,c(1,3:20)]
		metadata$CB <- rownames(metadata)
		# get gene annotations for hg38
		annotation <- GetGRangesFromEnsDb(ensdb = EnsDb_annotation);
		seqlevelsStyle(annotation) <- "UCSC";
		# create a Seurat object containing the RNA adata
		samp.list <- CreateSeuratObject( counts = sample$`Gene Expression`,project = paste0("RNAseq_",metafile$Condition[i]), assay = "RNA", meta.data = metadata,min.cells = min.cells,
		min.features = min.features);

		## adding RNASeq QC parameters
		#####################################
		DefaultAssay(samp.list) <- "RNA";
		samp.list$percent.mt_RNA <-  PercentageFeatureSet(samp.list, pattern = "^mt-");
		samp.list$percent.Mut_RNA <-  PercentageFeatureSet(samp.list, pattern = "^Mut$");
		samp.list$percent.Gapdh_RNA <-  PercentageFeatureSet(samp.list, pattern = "^Gapdh$");
 
		#expression Mut
		#############################
	    resCounts <- as.matrix(GetAssayData(samp.list), slot = "data");
	    resCountsSel <-resCounts[which(rownames(resCounts) %in% "Mut" ),]; #row:genes, cols:cells
	    samp.list$Counts_Mut <- resCountsSel;
	    samp.list$Mut_expression <- ifelse(samp.list$Counts_Mut>0,"Mut_expressed","Mut_not_Expressed");
	    not0 <-resCountsSel[which(resCountsSel>0)];
	    GenesNb <-dim(resCounts)[1];
	    cellsNb<- length(resCountsSel);
	    cells_expressing <- length(not0);
	    res_Mut <-data.frame(GenesNb=GenesNb,cells=cellsNb,cells_expressing=cells_expressing,
	    	maxValue_RNAExpression=max(not0),minValue_RNAExpression=min(not0),
	 		maxPercValue=max(unique(samp.list$percent.Mut_RNA)), minPercValue=min(unique(samp.list$percent.Mut_RNA)));
		openxlsx::write.xlsx(res_Mut, file = paste0(wd,f_results,f_qc,f_preIntegration,metafile$Condition[i], "_Mut_stats",date,".xlsx"),rowNames =TRUE, colNames =TRUE);

		#expression Gapdh to consider
		#  only Gapdh + cells
		#############################
	    resCountsSel <-resCounts[which(rownames(resCounts) %in% "Gapdh" ),]; #row:genes, cols:cells
	    samp.list$Counts_Gapdh <- resCountsSel;
	    samp.list$Gapdh_expression <- ifelse(samp.list$Counts_Gapdh>0,"Gapdh_expressed","Gapdh_not_Expressed");
	    not0 <-resCountsSel[which(resCountsSel>0)];
	    GenesNb <-dim(resCounts)[1];
	    cellsNb<- length(resCountsSel);
	    cells_expressing <- length(not0);
	    res_Gapdh <-data.frame(GenesNb=GenesNb,cells=cellsNb,cells_expressing=cells_expressing,
  			maxValue_RNAExpression=max(not0),minValue_RNAExpression=min(not0),
	 		maxPercValue=max(unique(samp.list$percent.Gapdh_RNA)), minPercValue=min(unique(samp.list$percent.Gapdh_RNA)));
				openxlsx::write.xlsx(res_Gapdh, file = paste0(wd,f_results,f_qc,f_preIntegration,metafile$Condition[i], "_Gapdh_stats",date,".xlsx"),rowNames =TRUE, colNames =TRUE);


		#threshold absed on perc expression
		samp.list$high.Mut <- ifelse(samp.list$percent.Mut_RNA > 2, 'High', 'Low')

		#cell cycle dissociation genes
		#############################
		s.genes <- (cc.genes$s.genes)
		g2m.genes <- (cc.genes$g2m.genes)
		#s.genes <- s.genes[s.genes %in% rownames(seurat)] # genes in dataset if the dataset is merge of all samples. not this case we work with a list
		#g2m.genes <- g2m.genes[g2m.genes %in% rownames(seurat)] # genes in dataset if the dataset is merge of all samples. not this case we work with a list

		# Genes upregulated during dissociation of tissue into single cells.

		samp.list <- CellCycleScoring(samp.list, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE);

		dfRes<- as.data.frame(table(samp.list@meta.data$Phase))
		colnames(dfRes)[2] <- paste0("nb of cells");
		Idents(samp.list) <- "orig.ident"
		colnames(dfRes)[1] <- "cellCycle Stage";

		dfResList[[i]]<- dfRes;
		names(dfResList)[i] <- paste0(metafile$Condition[i]);
		dataList[[i]] <- samp.list;
		names(dataList)[i] <- paste0(metafile$Condition[i]);
		i <- i+1;
	};
	assign(paste0("rnaList_",nameOutput),dataList,.GlobalEnv);
	assign(paste0("dfRes_",nameOutput),dfResList,.GlobalEnv);

};

#EnsDb_annotation <- EnsDb.Hsapiens.v86
EnsDb_annotation <- EnsDb.Mmusculus.v79

rnaObject_generation.function(metafile,EnsDb_annotation,name, 3,200,date)

#atac samples
atacObject_generation.function <- function(metafile,EnsDb_annotation,nameOutput,min.cells,min.features,date){
	dataList <- list();

	for (i in 1:length(metafile$Condition)){
		print(i);
		counts <-  Read10X_h5(paste0(wd,f_input,metafile$Sample[i],"_filtered_feature_bc_matrix.h5"));
		fragPath <- paste0(wd,f_input,metafile$Sample[i],"_atac_fragments.tsv.gz");
		metadata <- read.csv(file=paste0(wd,f_input,metafile$Sample[i],"_per_barcode_metrics.csv"), header = TRUE, row.names = 1);
		metadata <- metadata[,c(1:3,21:30)];
		metadata$CB <- rownames(metadata);
		total_fragments <- Signac::CountFragments(fragPath);
		rownames(total_fragments) <-total_fragments$CB;
		# Note: total_Fragments column frequency_counts is the same that the cellRanger atac_fragments,so combining
		#total_fragments[which(total_fragments$CB=="AAACAGCCATGGTTAT-1"),]
	 	#head(samp.list@metadata)
		colnames(total_fragments)[c(2)]<- c("atac_fragments");
		colnames(total_fragments)[c(3:5)] <- paste0("atac_", colnames(total_fragments)[c(3:5)])

		# get gene annotations for the genome
		annotation <- GetGRangesFromEnsDb(ensdb = EnsDb_annotation);
		seqlevelsStyle(annotation) <- "UCSC";

		# create ATAC assay and the seuratObject for atac
		atacAssay <- CreateChromatinAssay(
		  counts = counts$Peaks,
		  sep = c(":", "-"),
		  fragments = fragPath,
		  annotation = annotation,min.cells = min.cells,
	  	  min.features = min.features);

		samp.list_atac <- CreateSeuratObject(
	     counts = atacAssay,
	     assay = "ATAC",
	     project = paste0("ATACseq_",metafile$Condition[i]),
	     meta.data = metadata,
	     annotation = annotation);

	    # call peaks using MACS2
		peaks <- CallPeaks(samp.list_atac, macs2.path = "/software/anaconda3/envs/r4-base/bin/macs3")


		# remove peaks on nonstandard chromosomes and in genomic blacklist regions
		peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
		peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE) #blacklist_hg38_unified

		# quantify counts in each peak
		macs3_counts <- FeatureMatrix(
		  fragments = Fragments(samp.list_atac),
		  features = peaks,
		  cells = colnames(samp.list_atac));

		# create a new assay using the MACS2 peak set and add it to the Seurat object
		samp.list_atac[["peaks"]] <- CreateChromatinAssay(
		  counts = macs3_counts,
		  fragments = fragPath,
		  annotation = annotation);

		#needed line for scATACSeq:  when we merge the atacs later on, in orig.ident is not stored the sample name as happends in rnaseq is "seuratObject"
		# to avoiud losing the info of from whcih sample the peaks come from adding this lin. After merging the orig.ident should be set to sample.identity.
		#samp.list_atac$sample.identity <- metafile$Condition[i]; 


		## adding ATAC QC parameters
		#####################################
		# Nucleosome banding pattern: The histogram of DNA fragment sizes (determined from the paired-end sequencing reads)
		#   should exhibit a strong nucleosome banding pattern corresponding to the length of DNA wrapped around a single 
		#   nucleosome. We calculate this per single cell, and quantify the approximate ratio of mononucleosomal to
		#   nucleosome-free fragments (stored as nucleosome_signal)
		#
		# Transcriptional start site (TSS) enrichment score. The ENCODE project has defined an ATAC-seq targeting score based 
		#   on the ratio of fragments centered at the TSS to fragments in TSS-flanking regions
		#   (see https://www.encodeproject.org/data-standards/terms/). Poor ATAC-seq experiments typically will have a
		#   low TSS enrichment score. We can compute this metric for each cell with the TSSEnrichment() function, and the
		#    results are stored in metadata under the column name TSS.enrichment.
		#
		# Total number of fragments in peaks: A measure of cellular sequencing depth / complexity. Cells with very few reads
		#    may need to be excluded due to low sequencing depth. Cells with extremely high levels may represent doublets,
		#    nuclei clumps, or other artefacts.
		#
		# Fraction of fragments in peaks: Represents the fraction of all fragments that fall within ATAC-seq peaks. Cells with
		#    low values (i.e. <15-20%) often represent low-quality cells or technical artifacts that should be removed. Note
		#    that this value can be sensitive to the set of peaks used.
		#
		#Ratio reads in genomic blacklist regions The ENCODE project has provided a list of blacklist regions, representing 
		#   reads which are often associated with artefactual signal. Cells with a high proportion of reads mapping to these 
		#   areas (compared to reads mapping to peaks) often represent technical artifacts and should be removed. 
		#   ENCODE blacklist regions for human (hg19 and GRCh38), mouse (mm10), Drosophila (dm3), and C. elegans (ce10) are 
		#   included in the Signac package.

		#DefaultAssay(samp.list_atac) <- "ATAC";
		#head(samp.list_atac@meta.data)
	 	samp.list_atac$atac_mononucleosomal <- total_fragments[colnames(samp.list_atac), "atac_mononucleosomal"];
	 	samp.list_atac$atac_nucleosome_free <- total_fragments[colnames(samp.list_atac), "atac_nucleosome_free"];
	 	samp.list_atac$atac_reads_count <- total_fragments[colnames(samp.list_atac), "atac_reads_count"];
		
		# FRiP() will count the fraction of reads in peaks for each cell
		samp.list_atac <- FRiP(object = samp.list_atac, assay = 'ATAC',total.fragments = 'atac_fragments');
		#	
		samp.list_atac$pct_reads_in_peaks <- samp.list_atac$atac_peak_region_fragments / samp.list_atac$atac_fragments  * 100;
		samp.list_atac$percent.mt_ATAC <- samp.list_atac$atac_mitochondrial_reads / samp.list_atac$atac_fragments  * 100;
		
		#
		#samp.list_atac$blacklist_ratio <- samp.list_atac$blacklist_region_fragments / samp.list_atac$peak_region_fragments 
		##NOTE: is in the pipeline of analysis for atac alone, as we dont have this info in the metadata generated thus using other function to create it
		## The ratio of reads in genomic blacklist regions, that are known to artifactually accumulate reads in 
		## genome sequencing assays, can be diagnostic of low-quality cells. https://www.nature.com/articles/s41598-019-45839-z
		
		samp.list_atac$blacklist_fraction <- FractionCountsInRegion(object = samp.list_atac, assay = 'ATAC',regions = blacklist_mm10);
		
		samp.list_atac <- NucleosomeSignal(samp.list_atac);
		samp.list_atac <- TSSEnrichment(samp.list_atac,fast = FALSE);

		samp.list_atac$high.tss <- ifelse(samp.list_atac$TSS.enrichment > 2, 'High', 'Low');
		samp.list_atac$nucleosome_group <- ifelse(samp.list_atac$nucleosome_signal > 2, 'NS > 2', 'NS < 2');


		dataList[[i]] <- samp.list_atac;
		names(dataList)[i] <- paste0(metafile$Condition[i]);
		i <- i+1;
	};
	assign(paste0("atactList_",nameOutput),dataList,.GlobalEnv);

};

atacObject_generation.function(metafile,EnsDb_annotation,name,3,200,date)
	
save(atactList,rnaList,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_preintegration_prefiltering_Step1b.RData"))




###################################################################################################
# STEP2: qc ANALYSIS OF counts data before prefiltering
###################################################################################################

'''
##################
#Quality control
##################

bamBased_atacseqQC.function <- function(bamPath,nameOutput){
#based in bam file
	library(ATACseqQC)
	## input the bamFile from the ATACseqQC package 
	bamfile <- system.file("extdata", "GL1.bam", 
	                        package="ATACseqQC", mustWork=TRUE)
	bamfile.labels <- gsub(".bam", "", basename(bamfile))
	## generate fragement size distribution
	fragSize <- fragSizeDist(bamfile, bamfile.labels)
	source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"))
	# 1) Estimate the library complexity
	#bamQC(bamfile, outPath=NULL)
	estimateLibComplexity(readsDupFreq(bamfile))
	
	#2) generate fragement size distribution
	fragSize <- fragSizeDist(bamfile, bamfile.labels)
};


#min(samp.list@meta.data$pct_reads_in_peaks) #should be higher than 15-20 percent! 

'''


qc_multiome.function <- function(atacList,rnaseqList,nameOutput,date){
	options(scipen=999);

	for (i in 1:length(atacList)){
		print(i)

		pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_tssEnrichment,names(atacList)[i],"TSS_enrichments_groups",date,".pdf"), width = 5, height = 3,onefile=TRUE);
		plot(TSSPlot(atacList[[i]], group.by = 'high.tss') + NoLegend() +  ggtitle("TSS enrichment per signal group")+ geom_hline(yintercept = 1,linetype="dashed"));
		dev.off();
		pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_tssEnrichment,names(atacList)[i],"TSS_enrichments",date,".pdf"), width = 4, height = 3,onefile=TRUE);
		plot(TSSPlot(atacList[[i]], assay = "ATAC")  +  ggtitle("TSS enrichment all") + geom_hline(yintercept = 1,linetype="dashed"));
		dev.off();


		pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,names(rnaseqList)[i],"qc_vlnPlots_RNAbased",date,".pdf"), width = 10, height = 4)

		plot(VlnPlot(object = rnaseqList[[i]],
			features = c("nCount_RNA","nFeature_RNA" ,"percent.mt_RNA","blacklist_fraction","genes_dissoc1"),
		 	ncol = 5,pt.size = 0));
		dev.off();

		pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,names(atacList)[i],"qc_vlnPlots_ATACbased",date,".pdf"), width = 14, height = 4.5)
		plot(VlnPlot(object = atacList[[i]],
			features = c("nCount_ATAC","nFeature_ATAC" ,"percent.mt_ATAC","TSS.enrichment","nucleosome_signal","pct_reads_in_peaks"),
		 	ncol = 6,pt.size = 0));
		dev.off();


		pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,names(atacList)[i],"qc_expression_GenesInterest",date,".pdf"), width = 14, height = 4.5)
		plot(VlnPlot(object = rnaseqList[[i]],
			features = c("Counts_Mut","Counts_Gapdh"),
		 	ncol = 6,pt.size = 0));
		dev.off();


		p1a <- VlnPlot(atacList[[i]], features = "nCount_ATAC",y.max=quantile(atacList[[i]]$nCount_ATAC, 0.98)) +
		 geom_hline(yintercept = 1000,colour = "#000080", size = 1.2,linetype="dashed") + geom_hline(yintercept = 100000,colour = "#000080",size=1,linetype="dashed")+
		stat_summary(fun = mean, geom='point', size = 25, colour =  "#1034A6", shape = 95,show.legend = FALSE)


		p1b <- VlnPlot(rnaseqList[[i]], features = "nCount_RNA",y.max=quantile(rnaseqList[[i]]$nCount_RNA, 0.98)) + 
		geom_hline(yintercept = 1000,colour = "#000080", size = 1.2,linetype="dashed") + geom_hline(yintercept = 25000,colour = "#000080",size=1,linetype="dashed")+
		stat_summary(fun = mean, geom='point', size = 25, colour =  "#1034A6", shape = 95,show.legend = FALSE)



		p2a <- VlnPlot(atacList[[i]], features = "nFeature_ATAC") & geom_hline(yintercept = 25000,colour = "#000080", size = 1.2,linetype="dashed");
		p2b <- VlnPlot(rnaseqList[[i]], features = "nFeature_RNA") + stat_summary(fun = mean, geom='point', size = 25, colour =  "#1034A6", shape = 95,show.legend = FALSE);
		p3a <- VlnPlot(atacList[[i]], features = "percent.mt_ATAC",y.max=quantile(atacList[[i]]$percent.mt_ATAC, 0.95));
		p3b <- VlnPlot(rnaseqList[[i]], features = "percent.mt_RNA",y.max=quantile(rnaseqList[[i]]$percent.mt_RNA, 0.95)) & 
		geom_hline(yintercept = 10, colour = "#000080", size = 1.2,linetype="dashed");
		p3c <- VlnPlot(rnaseqList[[i]], features = "percent.Mut_RNA",y.max=quantile(rnaseqList[[i]]$percent.Mut_RNA, 0.95)) & 
		geom_hline(yintercept = 10, colour = "#000080", size = 1.2,linetype="dashed");
		p3d <- VlnPlot(rnaseqList[[i]], features = "percent.Gapdh_RNA",y.max=quantile(rnaseqList[[i]]$percent.Gapdh_RNA, 0.95)) & 
		geom_hline(yintercept = 10, colour = "#000080", size = 1.2,linetype="dashed");
		p3e <- VlnPlot(rnaseqList[[i]], features = "Counts_Mut",y.max=quantile(rnaseqList[[i]]$Counts_Mut, 0.95)) & 
		geom_hline(yintercept = 10, colour = "#000080", size = 1.2,linetype="dashed");
		p3f <- VlnPlot(rnaseqList[[i]], features = "Counts_Gapdh",y.max=quantile(rnaseqList[[i]]$Counts_Gapdh, 0.95)) & 
		geom_hline(yintercept = 10, colour = "#000080", size = 1.2,linetype="dashed");

		p4 <- VlnPlot(atacList[[i]], features = "TSS.enrichment") & geom_hline(yintercept = 1, colour = "#000080", size = 1.2,linetype="dashed");
		p5 <- VlnPlot(atacList[[i]], features = "nucleosome_signal") & geom_hline(yintercept = 2, colour = "#000080", size = 1.2,linetype="dashed");
		p6 <- VlnPlot(atacList[[i]], features = "pct_reads_in_peaks") + geom_hline(yintercept = 20, colour = "#000080", size = 1.2,linetype="dashed"); 
		p7 <- VlnPlot(atacList[[i]], features = "atac_peak_region_fragments",y.max=quantile(atacList[[i]]$atac_peak_region_fragments, 0.98)) + 
		geom_hline(yintercept = 3000, colour = "#000080", size = 1.2,linetype="dashed") + geom_hline(yintercept = 20000, colour = "#000080", size = 1.2,linetype="dashed");
		
		p8 <- VlnPlot(atacList[[i]], features = "blacklist_fraction",y.max=quantile(atacList[[i]]$blacklist_fraction, 0.98)) + geom_hline(yintercept = 0.05, colour = "#000080", size = 1.2,linetype="dashed") 
		plts <- wrap_plots(p1a, p1b,p2a,p2b, p3a,p3b,p4,p5,p6,p7,p8, ncol = 3);

		pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_individual",date,".pdf"), width = 4, height = 4.5,onefile=TRUE);
		plot(p1a);
		plot(p1b);
		plot(p2a);
		plot(p2b);
		plot(p3a);
		plot(p3b);
		plot(p3c);
		plot(p3d);
		plot(p3e);
		plot(p3f);
		plot(p4);
		plot(p5);
		plot(p6);
		plot(p7);
		plot(p8);
		dev.off();
		pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,names(atacList)[i],"_together",date,".pdf"), width = 20, height = 25.5,onefile=TRUE)
		plot(plts)
		dev.off();


		#tables quantiles
		##################3

		output.nCount.RNA <- as.data.frame(quantile(rnaseqList[[i]]$nCount_RNA, probs=seq(0,1,0.05)));
		output.nFeatures.RNA <- as.data.frame(quantile(rnaseqList[[i]]$nFeature_RNA, probs=seq(0,1,0.05)));
		output.percent.mt.RNA <- as.data.frame(quantile(rnaseqList[[i]]$percent.mt_RNA, probs=seq(0,1,0.05)));
		output.percent.Mut.RNA <- as.data.frame(quantile(rnaseqList[[i]]$percent.Mut_RNA, probs=seq(0,1,0.05)));
		output.percent.Gapdh.RNA <- as.data.frame(quantile(rnaseqList[[i]]$percent.Gapdh_RNA, probs=seq(0,1,0.05)));

		output.nCount.ATAC <- as.data.frame(quantile(atacList[[i]]$nCount_ATAC, probs=seq(0,1,0.05)));
		output.nFeatures.ATAC <- as.data.frame(quantile(atacList[[i]]$nFeature_ATAC, probs=seq(0,1,0.05)));
		output.percent.mt.ATAC <- as.data.frame(quantile(atacList[[i]]$percent.mt_ATAC, probs=seq(0,1,0.05)));

		output.nucleosome_signal.ATAC <- as.data.frame(quantile(atacList[[i]]$nucleosome_signal, probs=seq(0,1,0.05)));
		output.pct_reads_in_peaks.ATAC <- as.data.frame(quantile(atacList[[i]]$pct_reads_in_peaks, probs=seq(0,1,0.05)));
		output.blacklist_fraction.RNA <- as.data.frame(quantile(atacList[[i]]$blacklist_fraction, probs=seq(0,1,0.05)));
		output.TSS.enrichment.ATAC <- as.data.frame(quantile(atacList[[i]]$TSS.enrichment, probs=seq(0,1,0.05)));

		res_quantiles<- list(nCountRNA=output.nCount.RNA,nFeatures.RNA=output.nFeatures.RNA,percent.mt.RNA=output.percent.mt.RNA,
		nCountATAC=output.nCount.ATAC,nFeatures.ATAC=output.nFeatures.ATAC,percent.mt.ATAC=output.percent.mt.ATAC,
		nucleosome_signal.ATAC=output.nucleosome_signal.ATAC,pct_reads_in_peaks.ATAC=output.pct_reads_in_peaks.ATAC,
		blacklist_fraction.RNA=output.blacklist_fraction.RNA,TSS.enrichment.ATAC=output.TSS.enrichment.ATAC,
		percent.Mut.RNA=output.percent.Mut.RNA,
		percent.Gapdh.RNA=output.percent.Gapdh.RNA); #,percent.Mut.ATAC=output.percent.Mut.ATAC,percent.Gapdh.ATAC=output.percent.Gapdh.ATAC  );

		openxlsx::write.xlsx(res_quantiles, file = paste0(wd,f_results,f_qc,f_preIntegration,f_quantiles,names(atacList)[i], "_QC_quantiles",date,".xlsx"),rowNames =TRUE, colNames =TRUE);
		i <- i+1;
	};
};

qc_multiome.function(atactList,rnaList,name,date)





QCPlots_RNAseq.function <- function(samp.list,nameOutput,date, width,height,ncol){
	#library(grid)
	#library(gridBase)
	#library(gridExtra)
	options(scipen=5)

	percentMut.List<- list();
	percentMut.List.ymax<- list();
	percentGapdh.List<- list();
	percentGapdh.List.ymax<- list();
	ncount.List <-list();
	nFeature.List<-list();
	percentmt.List<- list();
	ncount.List.ymax <-list();
	nFeature.List.ymax<-list();
	percentmt.List.ymax<- list();
	perSample.List <-list();
	perSampleRidge.List <-list();
	perSampleScatter1.List <-list();
	perSampleScatter2.List <-list();
	perSampleScatter3.List <-list();
	perSampleScatter4.List <-list();
	perSampleScatter5.List <-list();
	perSampleScatter6.List <-list();
	perSampleScatter7.List <-list();
	perSampleScatter8.List <-list();
	perSampleScatter9.List <-list();
	cellsPerGenePlot.List <- list();

	for ( i in 1:length(samp.list)) {
		if (i==1) {
			ymax_nCount_vector <-	 max(samp.list[[1]]$nCount_RNA);	
			ymax_nFeature_vector <-	 max(samp.list[[1]]$nFeature_RNA);	
			ymax_percent.mt_vector <-	 max(samp.list[[1]]$percent.mt_RNA);	
		} else{
			ymax_nCount_vector.tmp <-	 max(samp.list[[i]]$nCount_RNA);	
			ymax_nFeature_vector.tmp <-	 max(samp.list[[i]]$nFeature_RNA);	
			ymax_percent.mt_vector.tmp <-	 max(samp.list[[i]]$percent.mt_RNA);	

			ymax_nCount_vector <- append(ymax_nCount_vector,ymax_nCount_vector.tmp);
			ymax_nFeature_vector <- append(ymax_nFeature_vector,ymax_nFeature_vector.tmp);
			ymax_percent.mt_vector <- append(ymax_percent.mt_vector,ymax_percent.mt_vector.tmp);
		};
		i <- i+1;
	};

	ymax_nCount <- max(ymax_nCount_vector);
	ymax_nFeature <- max(ymax_nFeature_vector);
	ymax_percent.mt <- max(ymax_percent.mt_vector);

	rm(i);

	for ( i in 1:length(samp.list)) {

		ncount.List.ymax[[i]] <- VlnPlot( samp.list[[i]], features = c("nCount_RNA"), ncol = 1, y.max=ymax_nCount);
		ncount.List[[i]]  <-VlnPlot( samp.list[[i]], features = c("nCount_RNA"), ncol = 1);

		nFeature.List.ymax[[i]] <- VlnPlot( samp.list[[i]], features = c("nFeature_RNA"), ncol = 1, y.max=ymax_nFeature);
		nFeature.List[[i]]  <-VlnPlot( samp.list[[i]], features = c("nFeature_RNA"), ncol = 1);

		percentMut.List.ymax[[i]] <- VlnPlot( samp.list[[i]], features = c("percent.Mut_RNA"), ncol = 1, y.max=ymax_percent.mt);
		percentMut.List[[i]]  <-VlnPlot( samp.list[[i]], features = c("percent.Mut_RNA"), ncol = 1);

		percentGapdh.List.ymax[[i]] <- VlnPlot( samp.list[[i]], features = c("percent.Gapdh_RNA"), ncol = 1, y.max=ymax_percent.mt);
		percentGapdh.List[[i]]  <-VlnPlot( samp.list[[i]], features = c("percent.Gapdh_RNA"), ncol = 1);

		percentmt.List.ymax[[i]] <- VlnPlot( samp.list[[i]], features = c("percent.mt_RNA"), ncol = 1, y.max=ymax_percent.mt);
		percentmt.List[[i]]  <-VlnPlot( samp.list[[i]], features = c("percent.mt_RNA"), ncol = 1);


		perSample.List[[i]]  <-VlnPlot( samp.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt_RNA"), ncol = 3);
		perSampleRidge.List[[i]]  <-RidgePlot(object =samp.list[[i]], features=c("nFeature_RNA","nCount_RNA", "percent.mt_RNA"), ncol = 1);
		perSampleScatter1.List[[i]]  <- FeatureScatter(object = samp.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA");
		perSampleScatter2.List[[i]]  <- FeatureScatter(object = samp.list[[i]], feature1 = "percent.mt_RNA", feature2 = "nFeature_RNA");
		perSampleScatter3.List[[i]]  <- FeatureScatter(object = samp.list[[i]], feature1 = "percent.mt_RNA", feature2 = "nCount_RNA");
		cellsPerGenePlot.List[[i]] <- sort(Matrix::rowSums(GetAssayData(samp.list[[i]]) >= 3)) #Plot the distribution of number of cells each gene is represented by

		perSampleScatter4.List[[i]]  <- FeatureScatter(object = samp.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.Mut_RNA");
		perSampleScatter5.List[[i]]  <- FeatureScatter(object = samp.list[[i]], feature1 = "nFeature_RNA", feature2 = "percent.Mut_RNA");
		perSampleScatter6.List[[i]]  <- FeatureScatter(object = samp.list[[i]], feature1 = "percent.mt_RNA", feature2 = "percent.Mut_RNA");

		perSampleScatter7.List[[i]]  <- FeatureScatter(object = samp.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.Gapdh_RNA");
		perSampleScatter8.List[[i]]  <- FeatureScatter(object = samp.list[[i]], feature1 = "nFeature_RNA", feature2 = "percent.Gapdh_RNA");
		perSampleScatter9.List[[i]]  <- FeatureScatter(object = samp.list[[i]], feature1 = "percent.mt_RNA", feature2 = "percent.Gapdh_RNA");



		names(ncount.List)[i] <- names(samp.list)[i];
		names(nFeature.List)[i] <- names(samp.list)[i];
		names(percentmt.List)[i] <- names(samp.list)[i];
		names(percentMut.List)[i] <- names(samp.list)[i];
		names(percentGapdh.List)[i] <- names(samp.list)[i];
		names(ncount.List.ymax)[i] <- names(samp.list)[i];
		names(nFeature.List.ymax)[i] <- names(samp.list)[i];
		names(percentmt.List.ymax)[i] <- names(samp.list)[i];
		names(perSample.List)[i] <- names(samp.list)[i];
		names(perSampleRidge.List)[i] <- names(samp.list)[i];
		names(perSampleScatter2.List)[i] <- names(samp.list)[i];
		names(perSampleScatter3.List)[i] <- names(samp.list)[i];
		names(perSampleScatter1.List)[i] <- names(samp.list)[i];
		names(perSampleScatter4.List)[i] <- names(samp.list)[i];
		names(perSampleScatter5.List)[i] <- names(samp.list)[i];
		names(perSampleScatter6.List)[i] <- names(samp.list)[i];
		names(perSampleScatter7.List)[i] <- names(samp.list)[i];
		names(perSampleScatter8.List)[i] <- names(samp.list)[i];
		names(perSampleScatter9.List)[i] <- names(samp.list)[i];
		i <- i+1;
	};

	#plotting

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_plots_perSample_",date,".pdf"), width = 8, height = 5,onefile=TRUE)
		for (i in 1:length(perSample.List)) {
			print(perSample.List[[i]]);
		i<-i+1;
	}

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_ridgeplot,nameOutput,"_qc_RidgePlots_perSample_",date,".pdf"), width = 8, height = 12,onefile=TRUE)
		for (i in 1:length(perSampleRidge.List)) {
			print(perSampleRidge.List[[i]]);
		i<-i+1;
	}

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_counts_feature_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(perSampleScatter1.List)) {
			print(perSampleScatter1.List[[i]]);
		i<-i+1;
	}

	dev.off();


	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percentMt_feature_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(perSampleScatter2.List)) {
			print(perSampleScatter2.List[[i]]);
		i<-i+1;
	}

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percentMt_counts_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(perSampleScatter3.List)) {
			print(perSampleScatter3.List[[i]]);
		i<-i+1;
	}

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_cells_per_GenePlot_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(cellsPerGenePlot.List)) {
			print(plot(cellsPerGenePlot.List[[i]] , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")) #Plot the distribution of number of cells each gene is represented by
		i<-i+1;
	}

	dev.off();

	#
	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percMut_counts_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(perSampleScatter4.List)) {
			print(perSampleScatter4.List[[i]]);
		i<-i+1;
	}

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percMut_features_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(perSampleScatter5.List)) {
			print(perSampleScatter5.List[[i]]);
		i<-i+1;
	}

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percMut_percentMt_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(perSampleScatter6.List)) {
			print(perSampleScatter6.List[[i]]);
		i<-i+1;
	}

	dev.off();
	#

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percGapdh_counts_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(perSampleScatter7.List)) {
			print(perSampleScatter7.List[[i]]);
		i<-i+1;
	}

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percGapdh_features_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(perSampleScatter8.List)) {
			print(perSampleScatter8.List[[i]]);
		i<-i+1;
	}

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percGapdh_percentMt_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
		for (i in 1:length(perSampleScatter9.List)) {
			print(perSampleScatter9.List[[i]]);
		i<-i+1;
	}

	dev.off();


	
	#all samples ins ame pdf one for ncount, one for nfeature and one for percent.mt
	nrow <- round(length(ncount.List)/ncol,digits=0)
	grid.arrange(grobs = ncount.List, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_nCounts.pdf"), 
	   plot = marrangeGrob(ncount.List, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);

	grid.arrange(grobs = nFeature.List, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_nFeature.pdf"), 
	   plot = marrangeGrob(nFeature.List, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);


	grid.arrange(grobs = percentmt.List, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_percent_mt.pdf"), 
	   plot = marrangeGrob(percentmt.List, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);
	#same max
	grid.arrange(grobs = ncount.List.ymax, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_nCounts_same_ymax.pdf"), 
	   plot = marrangeGrob(ncount.List.ymax, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);

	grid.arrange(grobs = nFeature.List.ymax, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_nFeature_same_ymax.pdf"), 
	   plot = marrangeGrob(nFeature.List.ymax, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);


	grid.arrange(grobs = percentmt.List.ymax, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_percent_mt_same_ymax.pdf"), 
	   plot = marrangeGrob(percentmt.List.ymax, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);

};



ncol <-4 #for the plot of all samples together 
QCPlots_RNAseq.function(rnaList,name,date,15,20,ncol)


combined_ridgePlot.function <- function(samp.list_RNA,samp.list_ATAC, nameOutput,date){

	combinedQC_RNA <- merge(x = samp.list_RNA[[1]], y=c(samp.list_RNA[[2]],samp.list_RNA[[3]], samp.list_RNA[[4]]))
	combinedQC_ATAC <- merge(x = samp.list_ATAC[[1]], y=c(samp.list_ATAC[[2]],samp.list_ATAC[[3]], samp.list_ATAC[[4]]))



	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_ridgeplot,nameOutput,"_qc_RidgePlots_allSamplesMerged_",date,".pdf"), width = 20, height = 9,onefile=TRUE);
		RidgePlot(object =combinedQC_RNA, features=c("nFeature_RNA","nCount_RNA", "percent.mt_RNA"), ncol = 3);
		RidgePlot(object =combinedQC_ATAC, features=c("nFeature_ATAC","nCount_ATAC", "percent.mt_ATAC"), ncol = 3);
	dev.off();

};

combined_ridgePlot.function(rnaList,name,atactList,date)
	combinedQC_ATAC <- merge(x = atactList[[1]], y=c(atactList[[2]],atactList[[3]], atactList[[4]]))


########################################################################################################
# STEP3 :combine datasets and filter out low quality cells
########################################################################################################

combinedQC_ATAC <- merge(x = atactList[[1]], y=c(atactList[[2]],atactList[[3]], atactList[[4]]))
combinedQC_RNA <- merge(x = rnaList[[1]], y=c(rnaList[[2]],rnaList[[3]], rnaList[[4]]))

save.image(file=paste0(wd,f_results, f_Rdata,"session_multiome_step3.RData"))

#
rnaList_subset <-  subset(
  x = combinedQC_RNA,
  subset = nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt_RNA < 10 );




atactList_subset <- subset(
  x = combinedQC_ATAC,
  subset = nCount_ATAC < 100000 &
    nCount_ATAC > 1000 &
    nucleosome_signal < 2 &
    blacklist_fraction < 0.05 &
    TSS.enrichment > 1);

save(rnaList_subset, atactList_subset,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_preintegration_filtered_step3.RData"))

########################################################################################################
# STEP4 :Redo the QC after filtering
########################################################################################################

QCPlots_RNAseq_merged.function <- function(samp.list,nameOutput,date, width,height,ncol){
	#library(grid)
	#library(gridBase)
	#library(gridExtra)
	options(scipen=5)
	
	ymax_nCount_vector <-	 max(samp.list$nCount_RNA);	
	ymax_nFeature_vector <-	 max(samp.list$nFeature_RNA);	
	ymax_percent.mt_vector <-	 max(samp.list$percent.mt_RNA);	
	ymax_nCount <- max(ymax_nCount_vector);
	ymax_nFeature <- max(ymax_nFeature_vector);
	ymax_percent.mt <- max(ymax_percent.mt_vector);

		ncount.List.ymax <- VlnPlot( samp.list, features = c("nCount_RNA"), ncol = 1, y.max=ymax_nCount);
		ncount.List  <-VlnPlot( samp.list, features = c("nCount_RNA"), ncol = 1);

		nFeature.List.ymax <- VlnPlot( samp.list, features = c("nFeature_RNA"), ncol = 1, y.max=ymax_nFeature);
		nFeature.List  <-VlnPlot( samp.list, features = c("nFeature_RNA"), ncol = 1);

		percentMut.List.ymax <- VlnPlot( samp.list, features = c("percent.Mut_RNA"), ncol = 1, y.max=ymax_percent.mt);
		percentMut.List  <-VlnPlot( samp.list, features = c("percent.Mut_RNA"), ncol = 1);

		percentGapdh.List.ymax <- VlnPlot( samp.list, features = c("percent.Gapdh_RNA"), ncol = 1, y.max=ymax_percent.mt);
		percentGapdh.List  <-VlnPlot( samp.list, features = c("percent.Gapdh_RNA"), ncol = 1);

		percentmt.List.ymax <- VlnPlot( samp.list, features = c("percent.mt_RNA"), ncol = 1, y.max=ymax_percent.mt);
		percentmt.List  <-VlnPlot( samp.list, features = c("percent.mt_RNA"), ncol = 1);


		perSample.List  <-VlnPlot( samp.list, features = c("nFeature_RNA", "nCount_RNA", "percent.mt_RNA"), ncol = 3);
		perSampleRidge.List  <-RidgePlot(object =samp.list, features=c("nFeature_RNA","nCount_RNA", "percent.mt_RNA"), ncol = 1);
		perSampleScatter1.List  <- FeatureScatter(object = samp.list, feature1 = "nCount_RNA", feature2 = "nFeature_RNA");
		perSampleScatter2.List  <- FeatureScatter(object = samp.list, feature1 = "percent.mt_RNA", feature2 = "nFeature_RNA");
		perSampleScatter3.List  <- FeatureScatter(object = samp.list, feature1 = "percent.mt_RNA", feature2 = "nCount_RNA");
		cellsPerGenePlot.List <- sort(Matrix::rowSums(GetAssayData(samp.list) >= 3)) #Plot the distribution of number of cells each gene is represented by

		perSampleScatter4.List  <- FeatureScatter(object = samp.list, feature1 = "nCount_RNA", feature2 = "percent.Mut_RNA");
		perSampleScatter5.List  <- FeatureScatter(object = samp.list, feature1 = "nFeature_RNA", feature2 = "percent.Mut_RNA");
		perSampleScatter6.List  <- FeatureScatter(object = samp.list, feature1 = "percent.mt_RNA", feature2 = "percent.Mut_RNA");

		perSampleScatter7.List  <- FeatureScatter(object = samp.list, feature1 = "nCount_RNA", feature2 = "percent.Gapdh_RNA");
		perSampleScatter8.List  <- FeatureScatter(object = samp.list, feature1 = "nFeature_RNA", feature2 = "percent.Gapdh_RNA");
		perSampleScatter9.List  <- FeatureScatter(object = samp.list, feature1 = "percent.mt_RNA", feature2 = "percent.Gapdh_RNA");

	#plotting

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_plots_perSample_",date,".pdf"), width = 8, height = 5,onefile=TRUE)
			print(perSample.List);

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_ridgeplot,nameOutput,"_qc_RidgePlots_perSample_",date,".pdf"), width = 8, height = 12,onefile=TRUE)
			print(perSampleRidge.List);


	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_counts_feature_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(perSampleScatter1.List);

	dev.off();


	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percentMt_feature_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(perSampleScatter2.List);


	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percentMt_counts_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(perSampleScatter3.List);

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_cells_per_GenePlot_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(plot(cellsPerGenePlot.List , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")) #Plot the distribution of number of cells each gene is represented by

	dev.off();

	#
	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percMut_counts_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(perSampleScatter4.List);


	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percMut_features_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(perSampleScatter5.List);


	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percMut_percentMt_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(perSampleScatter6.List);

	dev.off();
	#

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percGapdh_counts_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(perSampleScatter7.List);

	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percGapdh_features_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(perSampleScatter8.List);
	dev.off();

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_violin,nameOutput,"_qc_Scatter_percGapdh_percentMt_",date,".pdf"), width = 8, height = 8,onefile=TRUE)
			print(perSampleScatter9.List);

	dev.off();


	
	#all samples ins ame pdf one for ncount, one for nfeature and one for percent.mt
	nrow <- round(length(ncount.List)/ncol,digits=0)
	grid.arrange(grobs = ncount.List, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_nCounts.pdf"), 
	   plot = marrangeGrob(ncount.List, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);

	grid.arrange(grobs = nFeature.List, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_nFeature.pdf"), 
	   plot = marrangeGrob(nFeature.List, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);


	grid.arrange(grobs = percentmt.List, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_percent_mt.pdf"), 
	   plot = marrangeGrob(percentmt.List, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);
	#same max
	grid.arrange(grobs = ncount.List.ymax, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_nCounts_same_ymax.pdf"), 
	   plot = marrangeGrob(ncount.List.ymax, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);

	grid.arrange(grobs = nFeature.List.ymax, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_nFeature_same_ymax.pdf"), 
	   plot = marrangeGrob(nFeature.List.ymax, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);


	grid.arrange(grobs = percentmt.List.ymax, nrow = nrow) 
	ggsave(filename = paste0(wd,f_results,f_qc,f_preIntegration,f_violin, nameOutput,"_qc_percent_mt_same_ymax.pdf"), 
	   plot = marrangeGrob(percentmt.List.ymax, ncol=ncol,nrow=nrow,as.table=FALSE), 
  		width = width, height = height);
};

qc_multiome_merged.function <- function(atacList,rnaseqList,nameOutput,date){
	options(scipen=999);


	#tables quantiles
	##################
	output.nCount.RNA <- as.data.frame(quantile(rnaseqList$nCount_RNA, probs=seq(0,1,0.05)));
	output.nFeatures.RNA <- as.data.frame(quantile(rnaseqList$nFeature_RNA, probs=seq(0,1,0.05)));
	output.percent.mt.RNA <- as.data.frame(quantile(rnaseqList$percent.mt_RNA, probs=seq(0,1,0.05)));
	output.percent.Mut.RNA <- as.data.frame(quantile(rnaseqList$percent.Mut_RNA, probs=seq(0,1,0.05)));
	output.percent.Gapdh.RNA <- as.data.frame(quantile(rnaseqList$percent.Gapdh_RNA, probs=seq(0,1,0.05)));

	output.nCount.ATAC <- as.data.frame(quantile(atacList$nCount_ATAC, probs=seq(0,1,0.05)));
	output.nFeatures.ATAC <- as.data.frame(quantile(atacList$nFeature_ATAC, probs=seq(0,1,0.05)));
	output.percent.mt.ATAC <- as.data.frame(quantile(atacList$percent.mt_ATAC, probs=seq(0,1,0.05)));
	#output.percent.Mut.ATAC <- as.data.frame(quantile(atacList$percent.Mut_ATAC, probs=seq(0,1,0.05)));
	#output.percent.Gapdh.ATAC <- as.data.frame(quantile(atacList$percent.Gapdh_ATAC, probs=seq(0,1,0.05)));

	output.nucleosome_signal.ATAC <- as.data.frame(quantile(atacList$nucleosome_signal, probs=seq(0,1,0.05)));
	output.pct_reads_in_peaks.ATAC <- as.data.frame(quantile(atacList$pct_reads_in_peaks, probs=seq(0,1,0.05)));
	output.blacklist_fraction.RNA <- as.data.frame(quantile(atacList$blacklist_fraction, probs=seq(0,1,0.05)));
	output.TSS.enrichment.ATAC <- as.data.frame(quantile(atacList$TSS.enrichment, probs=seq(0,1,0.05)));

	res_quantiles<- list(nCountRNA=output.nCount.RNA,nFeatures.RNA=output.nFeatures.RNA,percent.mt.RNA=output.percent.mt.RNA,
	nCountATAC=output.nCount.ATAC,nFeatures.ATAC=output.nFeatures.ATAC,percent.mt.ATAC=output.percent.mt.ATAC,
	nucleosome_signal.ATAC=output.nucleosome_signal.ATAC,pct_reads_in_peaks.ATAC=output.pct_reads_in_peaks.ATAC,
	blacklist_fraction.RNA=output.blacklist_fraction.RNA,TSS.enrichment.ATAC=output.TSS.enrichment.ATAC,
	percent.Mut.RNA=output.percent.Mut.RNA,
	percent.Gapdh.RNA=output.percent.Gapdh.RNA); #,percent.Mut.ATAC=output.percent.Mut.ATAC,percent.Gapdh.ATAC=output.percent.Gapdh.ATAC  );

	openxlsx::write.xlsx(res_quantiles, file = paste0(wd,f_results,f_qc,f_preIntegration,f_quantiles,nameOutput, "_QC_quantiles",date,".xlsx"),rowNames =TRUE, colNames =TRUE);

};
qc_multiome_merged.function(atactList_subset,rnaList_subset,"selFinal",date);

rnaObject_quantfyGenes.function <- function(samp.list,nameOutput,date){

	#expression Mut
	#############################
    resCounts <- as.matrix(GetAssayData(samp.list), slot = "data");
    colnames(resCounts) <-samp.list$orig.ident
    resCountsSel <-resCounts[which(rownames(resCounts) %in% "Mut" ),]; #row:genes, cols:cells
    samp.list$Counts_Mut <- resCountsSel;
    samp.list$Mut_expression <- ifelse(samp.list$Counts_Mut>0,"Mut_expressed","Mut_not_Expressed");
    not0 <-resCountsSel[which(resCountsSel>0)];
  
    cond1 <- data.frame(counts_per_cell=not0[grep("Control_1",names(not0))],condition="Control_1");
    cond2 <- data.frame(counts_per_cell=not0[grep("Control_2",names(not0))],condition="Control_2");
    cond3 <- data.frame(counts_per_cell=not0[grep("Mut_1",names(not0))],condition="Mut_1");
    cond4 <- data.frame(counts_per_cell=not0[grep("Mut_2",names(not0))],condition="Mut_2");
    df.tmp <- bind_rows(cond1,cond2,cond3,cond4);
    cells_expressing_Mut_per_condition <- as.data.frame(df.tmp  %>% dplyr::group_by(condition ) %>% dplyr::summarise(cells_Expressing = n()));
    cells_expressing_Mut_per_condition$max_value_expression <- as.data.frame(df.tmp  %>% dplyr::group_by(condition ) %>% dplyr::summarise( max(counts_per_cell)))[,2];
    cells_expressing_Mut_per_condition$min_value_expression <- as.data.frame(df.tmp  %>% dplyr::group_by(condition ) %>% dplyr::summarise( min(counts_per_cell)))[,2];
    cells_expressing_Mut_per_condition$GenesNb <- dim(resCounts)[1];
    cells_expressing_Mut_per_condition$cellsNb <- length(resCountsSel);
    cells_expressing_Mut_per_condition$cells_expressing_in_all_conditions_together <- length(not0);



	#expression Gapdh
	#############################
    resCountsSel <-resCounts[which(rownames(resCounts) %in% "Gapdh" ),]; #row:genes, cols:cells
    samp.list$Counts_Gapdh <- resCountsSel;
    samp.list$Gapdh_expression <- ifelse(samp.list$Counts_Gapdh>0,"Gapdh_expressed","Gapdh_not_Expressed");   
    not0 <-resCountsSel[which(resCountsSel>0)];
    
    cond1 <- data.frame(counts_per_cell=not0[grep("Control_1",names(not0))],condition="Control_1");
    cond2 <- data.frame(counts_per_cell=not0[grep("Control_2",names(not0))],condition="Control_2");
    cond3 <- data.frame(counts_per_cell=not0[grep("Mut_1",names(not0))],condition="Mut_1");
    cond4 <- data.frame(counts_per_cell=not0[grep("Mut_2",names(not0))],condition="Mut_2");
    df.tmp <- bind_rows(cond1,cond2,cond3,cond4)
    cells_expressing_Gapdh_per_condition <- as.data.frame(df.tmp  %>% dplyr::group_by(condition ) %>% dplyr::summarise(cells_Expressing = n()));
    cells_expressing_Gapdh_per_condition$max_value_expression <- as.data.frame(df.tmp  %>% dplyr::group_by(condition ) %>% dplyr::summarise(max(counts_per_cell)))[,2];
    cells_expressing_Gapdh_per_condition$min_value_expression <- as.data.frame(df.tmp  %>% dplyr::group_by(condition ) %>% dplyr::summarise(min(counts_per_cell)))[,2];
    cells_expressing_Gapdh_per_condition$GenesNb <- dim(resCounts)[1];
    cells_expressing_Gapdh_per_condition$cellsNb <- length(resCountsSel);
    cells_expressing_Gapdh_per_condition$cells_expressing_in_all_conditions_together <- length(not0);

    res1 <- list(Mut=cells_expressing_Mut_per_condition,gapdh=cells_expressing_Gapdh_per_condition);

	openxlsx::write.xlsx(res1, file = paste0(wd,f_results,f_qc,f_preIntegration,nameOutput, "_Mut_Gapdh_stats",date,".xlsx"),rowNames =FALSE, colNames =TRUE);

};
#rnaObject_quantfyGenes.function(rnaList_subset,"selected","06282022");
rnaObject_quantfyGenes.function(rnaList_subset,"selFinal",date);


combined_ridgePlot_merged.function <- function(samp.list_RNA,samp.list_ATAC, nameOutput,date){

	pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,f_ridgeplot,nameOutput,"_qc_RidgePlots_allSamplesMerged_",date,".pdf"), width = 20, height = 9,onefile=TRUE);
		RidgePlot(object =samp.list_RNA, features=c("nFeature_RNA","nCount_RNA", "percent.mt_RNA"), ncol = 3);
		RidgePlot(object =samp.list_ATAC, features=c("nFeature_ATAC","nCount_ATAC", "percent.mt_ATAC"), ncol = 3);
	dev.off();

};
combined_ridgePlot_merged.function(rnaList_subset,"selFinal",atactList_subset,date);


##########################################################################
#stats about nb of counts, mt percentage etc based on quantiles
##########################################################################

#here we calculate que quantiles, after ordering the data in ascending orther, 0% in the min value. 
# 100% max value. the others are calculated based on the ranking
output_preFilteringRNAobject.nCountRNA <- as.data.frame(lapply(rnaList, function(x){quantile(x$nCount_RNA, probs=seq(0,1,0.05))}));
output_preFilteringRNAobject.nFeatures <- as.data.frame(lapply(rnaList, function(x){quantile(x$nFeature_RNA, probs=seq(0,1,0.05))}));
output_preFilteringRNAobject.percent.mt <- as.data.frame(lapply(rnaList, function(x){quantile(x$percent.mt_RNA, probs=seq(0,1,0.05))}));

#for the subseted objects
output_RNA_subsetObject.nCountRNA <- as.data.frame(quantile(rnaList_subset$nCount_RNA, probs=seq(0,1,0.05)));
output_RNA_subsetObject.nFeature_RNA <- as.data.frame(quantile(rnaList_subset$nFeature_RNA, probs=seq(0,1,0.05)));
output_RNA_subsetObject.percent.mt_RNA <- as.data.frame(quantile(rnaList_subset$percent.mt_RNA, probs=seq(0,1,0.05)));

#FOR ATACSeq
####################
output_preFilteringATACobject.nCountATAC <- as.data.frame(lapply(atactList, function(x){quantile(x$nCount_ATAC, probs=seq(0,1,0.05))}));
output_preFilteringATACobject.nFeatures <- as.data.frame(lapply(atactList, function(x){quantile(x$nFeature_ATAC, probs=seq(0,1,0.05))}));
output_preFilteringATACobject.percent.mt <- as.data.frame(lapply(atactList, function(x){quantile(x$percent.mt_ATAC, probs=seq(0,1,0.05))}));

#for the subseted objects
output_ATAC_subsetObject.nCountATAC <- as.data.frame(quantile(atactList_subset$nCount_ATAC, probs=seq(0,1,0.05)));
output_ATAC_subsetObject.nFeature_ATAC <- as.data.frame(quantile(atactList_subset$nFeature_ATAC, probs=seq(0,1,0.05)));
output_ATAC_subsetObject.percent.mt_ATAC <- as.data.frame(quantile(atactList_subset$percent.mt_ATAC, probs=seq(0,1,0.05)));

res_Basic_QC_stats_pre_post_filtering.List <- list(RNApreFiltering.nCount=output_preFilteringRNAobject.nCountRNA,
RNApreFiltering.nFeature=output_preFilteringRNAobject.nFeatures,
RNApreFiltering.mtPerc=output_preFilteringRNAobject.percent.mt,
RNApostFiltering.nCount=output_RNA_subsetObject.nCountRNA,
RNApostFiltering.nFeature=output_RNA_subsetObject.nFeature_RNA,
RNApostFiltering.mtPerc=output_RNA_subsetObject.percent.mt_RNA,
ATACpreFiltering.nCount=output_preFilteringATACobject.nCountATAC,
ATACpreFiltering.nFeature=output_preFilteringATACobject.nFeatures,
ATACpreFiltering.mtPerc=output_preFilteringATACobject.percent.mt,
ATACpostFiltering.nCount=output_ATAC_subsetObject.nCountATAC,
ATACpostFiltering.nFeature=output_ATAC_subsetObject.nFeature_ATAC,
ATACpostFiltering.mtPerc=output_ATAC_subsetObject.percent.mt_ATAC);
openxlsx::write.xlsx(res_Basic_QC_stats_pre_post_filtering.List, file = paste0(wd, f_results,f_qc ,f_after_enbedding,f_tables,name, "_QCquantiles_pre_post_filtering.xlsx"),rowNames =TRUE, colNames =TRUE);



###################################################################################################
##STEP5 snRNASeq proccessing steps before integration with snATAC-Seq
###################################################################################################
# Standard analysis of RNASeq
rnaList_subset <- NormalizeData(rnaList_subset)

#rnaList_subset <- FindVariableFeatures(rnaList_subset)
# Instead using find var genes using all common genes expressed in more than 2 counts in at least 10 cells:

find_commonGenes.function <- function(min.value,min.cells,samp.list_RNA,nameOutput){
	num.cells <- Matrix::rowSums(GetAssayData(samp.list_RNA, slot = "count") > min.value)
	genes.use <- names(num.cells[which(num.cells >= min.cells)])
	print(paste0("Nb of genes in common: ", length(genes.use))); #3783 8612
	VariableFeatures(samp.list_RNA) <- genes.use
	assign("nameOutput",samp.list_RNA,.GlobalEnv);
	assign("nameOutput",samp.list_RNA,.GlobalEnv);
};
find_commonGenes.function(2,10,rnaList_subset,"rnaList_subset");


rnaList_subset <- ScaleData(rnaList_subset);
rnaList_subset <- RunPCA(rnaList_subset, npcs = 100, verbose = FALSE);
rnaList_subset <- RunUMAP(rnaList_subset, reduction = "pca", dims = 1:100)
rnaList_subset <- JackStraw(rnaList_subset, num.replicate=100,dims=50);

elbowPlot<- ElbowPlot(rnaList_subset,ndims=50);

pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,"elbowPlot_nb_dimensions_rna.pdf"), width = 8, height = 8,onefile=TRUE)
print(elbowPlot);
dev.off();

#rerun UMAP to remove the noise in the clustering that using higher dimensions that needed will create
rnaList_subset <- RunPCA(rnaList_subset, npcs = 30, verbose = FALSE)
rnaList_subset <- RunUMAP(rnaList_subset, reduction = "pca", dims = 1:30)
rnaList_subset <- FindNeighbors(object = rnaList_subset,reduction = 'pca',dims =c(1:30));

save(rnaList_subset,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"rnaSeq_multiome_processing_preintegration_step5.RData"))
save.image(file=paste0(wd,f_results, f_Rdata,"session_multiome_step5.RData"))





###################################################################################################
##STEP6 snATAC-Seq processing steps before integration witn snRNA-Seq
###################################################################################################
## ATAC analysis add gene annotation information
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotations) <- "UCSC"

# We exclude the first dimension as this is typically correlated with sequencing depth
atactList_subset <- RunTFIDF(atactList_subset)
atactList_subset <- FindTopFeatures(atactList_subset, min.cutoff = "q10"); #same threshold I sused for the rnaseq
atactList_subset <- RunSVD(atactList_subset)

#identify the lsi components correlated with sequencing. Removing those component with abs(r)>0.5
corR_atac <- DepthCor(atactList_subset,n=30) #is a very strong correlation between the first and third LSI component and the total number of counts for the cell, so we will perform downstream steps without those components
pdf(file=paste0(wd,f_results,f_qc,f_preIntegration,"corrPlot_nb_dimensions_atac.pdf"), width = 8, height = 8,onefile=TRUE)
print(corR_atac);
dev.off();

# UMAP using dims=30, considering the res parameter optimal for rnaseq and to treat the atac in a similar way
atactList_subset <- RunUMAP(atactList_subset, reduction = "lsi", dims = c(2,4:30), reduction.name = "umap.atac", reduction.key = "atacUMAP_") 
atactList_subset <- FindNeighbors(object = atactList_subset,reduction = 'lsi',dims =c(2,4:30));


# quantify gene activity
gene.activities <- GeneActivity(atactList_subset, features = VariableFeatures(rnaList_subset)) 
# NOTE : as I used as variable features all the genes quantified gives error in the quantifications we need to use selected top 2000 max

# add gene activities as a new assay
atactList_subset[['geneActivity']] <- CreateAssayObject(counts = gene.activities);

save(atactList_subset,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"atacSeq_multiome_processing_preintegration_step6.RData"))

###################################################################################################
##STEP7 a) Find clusters based on different resolution parameters
##      b) Identify cell types markers in snRNASeq and transfer that info to the snATACSeq
###################################################################################################

######### a) Find clusters based on different resolution parameters
###################################################################################################
findClusters.function <- function(combined.samp.list,dataType,nameSampleList, resolutionVector,nameOutput,height,width,date,reductionName){
	findClustersList <- list();
	clusterPlotsList <- list();
	comb.markersList <- list();
	cellcyclePlotsList <- list();

	##### check the nb of clusters  depending on the resolution ####
	## NOTE these results change depending in the function runUMASP how many clusters you allowed.!
	combined.samp.list <- FindClusters(combined.samp.list, resolution = resolutionVector);
	plot_perCluster <- DimPlot(object = combined.samp.list, group.by=grep("_res",colnames(combined.samp.list@meta.data),value = TRUE), ncol=2 , pt.size=3.0, reduction = reductionName, label = T);

	pdf(file=paste0("Clustering_perResolution_",nameOutput,"_UMAP_",date,".pdf"), height = height, width = width);
	print(plot_perCluster);
	dev.off();

	#####################another way
	for (i in 1:length(resolutionVector)){
		findClusters<- FindClusters(combined.samp.list, resolution = resolutionVector[i]); #can eb done in only one step
		findClustersList[[i]] <- findClusters;
		names(findClustersList)[i] <- gsub("\\.","",resolutionVector[i]);
		
		clusterPlotsList[[i]] <- DimPlot(findClusters, reduction = reductionName, label = TRUE, raster = FALSE) +
		labs(title = paste0("Clustering with ",resolutionVector[i], " Resolution"));

		if( dataType=="RNA") {
			cellcyclePlotsList[[i]] <- 	DimPlot(object = findClusters, pt.size=0.5, group.by = "Phase", reduction = reductionName );	
		};

		i <- i+1;
	};

	pdf(file=paste0(wd,f_results,f_foundClusters,f_plots, nameOutput,"_","Umap_diff_res_",date,".pdf"),onefile=TRUE,width = width, height = height)
		for (i in 1:length(clusterPlotsList)) {
			print(clusterPlotsList[[i]]);
			i<-i+1;
		};
	dev.off();

	if( dataType=="RNA") {
		pdf(file=paste0(wd,f_results,f_foundClusters,f_plots, nameOutput,"_","cellCycle_umap_diff_res_",date,".pdf"),onefile=TRUE,width = 9, height = 6)
			for (i in 1:length(cellcyclePlotsList)) {
				print(cellcyclePlotsList[[i]]);
				i<-i+1;
			};
		dev.off();
	};

	#assign(paste0("clusterPlotsList_",nameOutput),clusterPlotsList,.GlobalEnv);
	assign(paste0("findClustersList_",nameOutput),findClustersList,.GlobalEnv);
	assign(paste0(nameSampleList),combined.samp.list,.GlobalEnv);
};
DefaultAssay(atactList_subset) <- "ATAC"
resolutionVector <- c(0.01,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.4,1.5,1.6,1.8,2,2.4,3,3.5);
findClusters.function(atactList_subset,"ATAC","atactList_subset", resolutionVector,paste0(name,"_atac"),8,12,date,"umap.atac");
resolutionVector <- c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.4,1.5,1.6,1.8,2,2.4,3,3.5); 
findClusters.function(rnaList_subset,"RNA","rnaList_subset", resolutionVector,paste0(name,"_rna"),8,12,date,"umap");
	


# How many cluster each resolution produces?

sapply(grep("_res",colnames(atactList_subset@meta.data),value = TRUE),
       function(x) length(unique(atactList_subset@meta.data[,x])))

#ATAC_snn_res.0.05  ATAC_snn_res.0.1 
#                9                15
sapply(grep("_res",colnames(rnaList_subset@meta.data),value = TRUE),
       function(x) length(unique(rnaList_subset@meta.data[,x])))

#RNA_snn_res.0.05  RNA_snn_res.0.1  
#               8               10               

#redo the cell cycle plot with smaller size
cellcyclePlots <- 	DimPlot(object = rnaList_subset, pt.size=0.5, group.by = "Phase", reduction = "umap" );	
pdf(file=paste0(wd,f_results,f_foundClusters,f_plots, name,"_","cellCycle_umap_diff_res_",date,"_vsmall.pdf"),onefile=TRUE,width = 6, height = 4)
plot(cellcyclePlots)
dev.off();
pdf(file=paste0(wd,f_results,f_foundClusters,f_plots, name,"_","cellCycle_umap_diff_res_",date,"_vBig.pdf"),onefile=TRUE,width = 8, height = 6)
plot(cellcyclePlots)
dev.off();


save(rnaList_subset, atactList_subset,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_preintegration_filtered_step7.RData"))
save.image(file=paste0(wd,f_results, f_Rdata,"session_multiome_step7.RData"))

######### b)Identify cell types markers in snRNASeq and transfer that info to the snATACSeq
###################################################################################################

##########################
#functions for plotting
##########################


FeaturePlot_markers.function <- function(findClustersList,assay,width,height,nameOutput,date, metadataCondition,featureName,folderName,nameMarker,factor){
	#note: the factor argument is to be set as "yes" when more than 4 genes will be plotted that way all the plots are scaled 
	# in a similar way doesnt matter if you want to plot 15 or 30, the plot becomes bigger to maintain the plottting ratio

	dir.create(file.path(wd, f_results,f_markers,f_plots,folderName), showWarnings = F);
	DefaultAssay(findClustersList) <- assay;

	if (factor=="yes"){
		factorPlot <-  round(length(featureName)/4,digits=0);
	} else{
		factorPlot<-1
	};
	#
	if (metadataCondition=="Sample") {
		plot <- FeaturePlot(findClustersList, features = featureName, raster = FALSE,split.by = metadataCondition);
	
	pdf(file=paste0(wd,f_results,f_markers,f_plots,folderName,nameMarker,"_umap_",nameOutput,
		"_",metadataCondition,"_",date,".pdf"),onefile=TRUE,width = width*factorPlot, height = height*factorPlot)
	print(plot);
	dev.off();


	} else if (metadataCondition=="NULL"){
		plot <- FeaturePlot(findClustersList, features = featureName, raster = FALSE);

	pdf(file=paste0(wd,f_results,f_markers,f_plots,folderName,nameMarker,"_umap_",nameOutput,
		"_",metadataCondition,"_",date,".pdf"),onefile=TRUE,width = width*factorPlot, height = height*factorPlot)
	print(plot);
	dev.off();


	} else{
		plot <- FeaturePlot(findClustersList, features = featureName, raster = FALSE,split.by = metadataCondition);
	pdf(file=paste0(wd,f_results,f_markers,f_plots,folderName,nameMarker,"_umap_",nameOutput,
		"_",metadataCondition,"_",date,".pdf"),onefile=TRUE,width = width*factorPlot, height = height*factorPlot)
	print(plot);
	dev.off();

	};
	#
};


FeaturePlot_FinalClustering.function <- function(findClustersList,width,height,nameOutput,date, metadataCondition,featureName,folderName){
  library(RColorBrewer);

  dir.create(file.path(wd, f_results,f_markers,f_plots,folderName), showWarnings = F);
  DefaultAssay(findClustersList) <- "RNA";
  #
  if (metadataCondition=="Sample") {
    plot <- FeaturePlot(findClustersList, features = featureName, raster = FALSE,split.by = metadataCondition) & scale_color_viridis_c() 
  
  pdf(file=paste0(wd,f_results,f_markers,f_plots,folderName,featureName,"_umap_",nameOutput,
    "_",metadataCondition,"_",date,".pdf"),onefile=TRUE,width = width, height = height)
  print(plot);
  dev.off();


  } else if (metadataCondition=="NULL"){
    plot <- FeaturePlot(findClustersList, features = featureName, raster = FALSE) & scale_color_viridis_c() & 
      theme(legend.text=element_text(size=15),legend.spacing.y = unit(2, 'cm'),legend.key.height = unit(2.3, "cm"),
        axis.title = element_text(size=22,face="bold"),axis.text.x = element_text(size=18, color="black"),axis.text.y = element_text(size=18, color="black"), 
        axis.line = element_line(size=1.5))


  png(file=paste0(wd,f_results,f_markers,f_plots,folderName,featureName,"_umap_",nameOutput,
    "_",metadataCondition,"_",date,".png"), height = height, width = width, units = "in", res = 1200);
  print(plot);
  dev.off();


  } else{
    plot <- FeaturePlot(findClustersList, features = featureName, raster = FALSE,split.by = metadataCondition) && scale_color_viridis_c() & 
    theme(legend.text=element_text(size=15),legend.spacing.y = unit(2, 'cm'),legend.key.height = unit(2.3, "cm"),
            axis.title = element_text(size=22,face="bold"),axis.text.x = element_text(size=18, color="black"),axis.text.y = element_text(size=18, color="black"), 
            axis.line = element_line(size=1.5))
;
  png(file=paste0(wd,f_results,f_markers,f_plots,folderName,featureName,"_umap_",nameOutput,
    "_",metadataCondition,"_",date,".png"), height = height, width = width, units = "in", res = 1200);
  print(plot);
  dev.off();

  };
  #
};

#########################################################
#defining marker genes to test in different 
# vector lists and run the FeaturePlot_markers.function
#########################################################
#Markers

bcells.markers <- c("Ins1", "Ins2","Hadh","G6pc2"); 
acells.markers <- c("Gcg", "Ttr", "Gc");
delta.markers <- c("Sst");


#other markers possible from azymuth
#https://azimuth.hubmapconsortium.org/references/#Human%20-%20Pancreas
#bcells.markers <- c("Ins1", "Ins2") #IAPP, INS, DLK1, INS-IGF2, G6PC2, HADH, ADCYAP1, GSN, NPTX2, C12orf75
#acells.markers <- c(GCG, TTR, PPP1R1A, CRYBA2, TM4SF4, MAFB, GC, GPX3, PCSK2, PEMT "MPxi1","Ins","ck19","Amylase")
#gamma.markers <-c(PPY, AQP3, MEIS2, ID2, GPC5-AS1, CARTPT, PRSS23, ETV1, PPY2, TUBB2A)
#delta.markers <- c(SST, RBP4, SERPINA1, RGS2, PCSK1, SEC11C, HHEX, LEPR, MDK, LY6H)
#acinar.markers <- c(REG1A, PRSS1, CTRB2, CTRB1, REG1B, CELA3A, PRSS2, REG3A, CPA1, CLPS)
#ductal.markers <- c(SPP1, MMP7, IGFBP7, KRT7, ANXA4, SERPINA1, LCN2, CFTR, KRT19, SERPING1)
#epsilon.markers <- c(BHMT, VSTM2L, PHGR1, TM4SF5, ANXA13, ASGR1, DEFB1, GHRL, COL22A1, OLFML3)
#swann.markers <- c(NGFR, CDH19, UCN2, SOX10, S100A1, PLP1, TSPAN11, WNT16, SOX2, TFAP2A)
#active_stellate.markers <- C(COL1A1, COL1A2, COL6A3, COL3A1, TIMP3, TIMP1, CTHRC1, SFRP2, BGN, LUM)
#immune.markers.markers <- c(ACP5, APOE, HLA-DRA, TYROBP, LAPTM5, SDS, FCER1G, C1QC, C1QB, SRGN)
#cycling.markers <- c(UBE2C, TOP2A, CDK1, BIRC5, PBK, CDKN3, MKI67, CDC20, CCNB2, CDCA3)
#quiescent_stellate.markers <-c(RGS5, C11orf96, FABP4, CSRP2, IL24, ADIRF, NDUFA4L2, GPX3, IGFBP4, ESAM)
#endotelial.markers <-c(PLVAP, RGCC, ENG, PECAM1, ESM1, SERPINE1, CLDN5, STC1, MMP1, GNG11)

#

###################################################################################################
##STEP8 Identify anchors
###################################################################################################

# normalize gene activities
DefaultAssay(atactList_subset) <- "geneActivity"
atactList_subset <- NormalizeData(atactList_subset,
	assay = 'geneActivity',
	normalization.method = 'LogNormalize',
	scale.factor = median(atactList_subset$nCount_RNA));

atactList_subset <- ScaleData(atactList_subset, features = rownames(atactList_subset))

save(rnaList_subset, atactList_subset,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_preintegration_filtered_step8.RData"))

transfer.anchors <- FindTransferAnchors(reference = rnaList_subset, query = atactList_subset, features = VariableFeatures(object = rnaList_subset),
    reference.assay = "RNA", query.assay = "geneActivity", reduction = "cca")
save.image(file=paste0(wd,f_results, f_Rdata,"session_multiome_step8.RData"))

save(transfer.anchors,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"transfer_anchors_step8.RData"))

###################################################################################################
##STEP9  coembedding of scRNASeq with scATACSeq 
###################################################################################################
#from bigmem QUEU,  step SRONGLY consuming resources: script multiome_coembbeding_step9.R #merge
#
# NOTE: The imputation can be done using 1) only the variable genes from scRNA-seq, or the full transriptome
# I am doing the full transcriptome.
# Refdata contains scRNA-seq expression matrix for the scRNA-seq cells.  
# imputation  will contain an imputed scRNA-seq matrix for each of the ATAC cells


genes.use <- VariableFeatures(rnaList_subset)
refdata <- GetAssayData(rnaList_subset, assay = "RNA", slot = "data");
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
	weight.reduction = atactList_subset[["lsi"]],
    dims = c(2,4:30));
atactList_subset[["RNA"]] <- imputation; #Transfering 23179 features onto reference data

coembed_multiome <- merge(x = rnaList_subset, y = atactList_subset);
coembed_multiome <- ScaleData(coembed_multiome, features = rownames(coembed_multiome), do.scale = FALSE);
coembed_multiome <- RunPCA(coembed_multiome, features = genes.use, verbose = FALSE,dims = 1:30);
coembed_multiome <- RunUMAP(coembed_multiome, dims = 1:30, reduction.name = "umap.multiome.embed", reduction.key = "multiomeUMAP_");
coembed_multiome <- FindNeighbors(object = coembed_multiome,reduction = 'pca',dims =c(1:30));

#save objects
save(rnaList_subset, atactList_subset,
	file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_preintegration_filtered_step9.RData"));
save(coembed_multiome, file=paste0(wd, f_results, f_Rdata,f_after_enbedding,"coembedded_multiome_step9.RData"));
save.image(file=paste0(wd,f_results, f_Rdata,"session_multiome_step9.RData"));

#DimPlot(coembed_multiome, group.by = c("orig.ident", "seurat_annotations"))

##########################################################
#plot those markers and MT percentage
##########################################################

FeaturePlot_FinalClustering.function(rnaList_subset,8,6,"MtDysfunction",date, "NULL","percent.mt_RNA","MTpercentage_RNA/");
FeaturePlot_FinalClustering.function(atactList_subset,8,6,"MtDysfunction",date, "NULL","percent.mt_ATAC","MTpercentage_RNA/");

FeaturePlot_markers.function(rnaList_subset,"RNA",9,7,paste0(name,"_rna"),date, "NULL",bcells.markers,"Bcells/","bcells","no");
FeaturePlot_markers.function(atactList_subset,"geneActivity",9,7,paste0(name,"_atac"),date, "NULL",bcells.markers,"Bcells/","bcells","no");

FeaturePlot_markers.function(rnaList_subset,"RNA",9,7,paste0(name,"_rna"),date, "NULL",acells.markers,"Acells/","acells","no");
FeaturePlot_markers.function(atactList_subset,"geneActivity",9,7,paste0(name,"_atac"),date, "NULL",acells.markers,"Acells/","acells","no");

FeaturePlot_markers.function(rnaList_subset,"RNA",9,7,paste0(name,"_rna"),date, "NULL",delta.markers,"Dcells/","dcells","no");
FeaturePlot_markers.function(atactList_subset,"geneActivity",9,7,paste0(name,"_atac"),date, "NULL",delta.markers,"Dcells/","dcells","no");



FeaturePlot_markers.function(atactList_subset,"RNA",9,7,paste0(name,"_atac"),date, "NULL",bcells.markers,"Bcells/","bcells_atac_RNAimput","no");
FeaturePlot_markers.function(atactList_subset,"RNA",9,7,paste0(name,"_atac"),date, "NULL",acells.markers,"Acells/","acells_atac_RNAimput","no");
FeaturePlot_markers.function(atactList_subset,"RNA",9,7,paste0(name,"_atac"),date, "NULL",delta.markers,"Dcells/","dcells_atac_RNAimput","no");


###################################################################################################
## reassesing the ATACSeq object to decrease level of comlexity show in umap to 
#       reflect better the populations we want to highlight
###################################################################################################

#DEA_thresholds adj p-value< 0.05 -1<log2(FC)> 058  

DefaultAssay(atactList_subset) <- "ATAC"

atactList_subset_lessdim <- RunUMAP(atactList_subset, reduction = "lsi", dims = c(2,4:15), reduction.name = "umap.atac", reduction.key = "atacUMAP_") 
atactList_subset_lessdim <- FindNeighbors(object = atactList_subset_lessdim,reduction = 'lsi',dims =c(2,4:15));


# quantify gene activity
gene.activities <- GeneActivity(atactList_subset_lessdim, features = VariableFeatures(rnaList_subset))

# add gene activities as a new assay
atactList_subset_lessdim[['geneActivity']] <- CreateAssayObject(counts = gene.activities);
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
	weight.reduction = atactList_subset_lessdim[["lsi"]],
    dims = c(2,4:15));
atactList_subset_lessdim[["RNA"]] <- imputation; #Transfering 23179 features onto reference data

save(atactList_subset_lessdim,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"atacSeq_multiome_processing_preintegration_less_dims_step9_BIS.RData"))
#load(file=paste0(wd, f_results, f_Rdata,f_preIntegration,"atacSeq_multiome_processing_preintegration_less_dims_step9_BIS.RData")) #atactList_subset_lessdim,

#plot the umap
resolutionVector <- c(0.02,0.03,0.04,0.05); 
findClusters.function(atactList_subset_lessdim,"ATAC","atactList_subset_lessdim", resolutionVector,paste0(name,"_atac_lessDims"),8,12,date,"umap.atac");
sapply(grep("_res",colnames(atactList_subset_lessdim@meta.data),value = TRUE),
       function(x) length(unique(atactList_subset_lessdim@meta.data[,x])))

FeaturePlot_FinalClustering.function(atactList_subset_lessdim,8,6,"MtDysfunction",date, "NULL","percent.mt_ATAC","MTpercentage_RNA/");
FeaturePlot_markers.function(atactList_subset_lessdim,"RNA",12,8,paste0(name,"_atac_lessdim"),date, "NULL",bcells.markers,"Bcells/","bcells","no");
FeaturePlot_markers.function(atactList_subset_lessdim,"RNA",12,8,paste0(name,"_atac_lessdim"),date, "NULL",acells.markers,"Acells/","acells","no");
FeaturePlot_markers.function(atactList_subset_lessdim,"RNA",6,6,paste0(name,"_atac_lessdim"),date, "NULL",delta.markers,"Dcells/","dcells","no");
#NOTE: 
######################################################################################################
# conclusion : althought the clustering is a bit better per cell type there is one community split 
# that looks off in the centre. So i will still prefer to use the UMAAP with more dims and all the
# communities split than this one where two diff communities are so closely together. 
#Thus continuing downstream analyses with atactList_subset not with the lower dims.
###################################################################################################


#plots
################
resolutionVector <- c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.4,1.5,1.6,1.8,2,2.4,3,3.5); 
#findClusters.function(coembed_multiome,"ATAC","coembed_multiome", resolutionVector,paste0(name,"_coembed_multiome"),7,10,"071522","umap.multiome.embed");
findClusters.function(coembed_multiome,"RNA","multiome_RNA", resolutionVector,paste0("multiome"),7,10,"071522","umap.multiome.embed");

# How many cluster each resolution produces?

sapply(grep("_res",colnames(coembed_multiome@meta.data),value = TRUE),
       function(x) length(unique(coembed_multiome@meta.data[,x])))

# RNA_snn_res.0.05   RNA_snn_res.0.1     
#                9                11                                
#ATAC_snn_res.0.05  ATAC_snn_res.0.1  
#               10                16                



FeaturePlot_FinalClustering.function(coembed_multiome,8,6,"MtDysfunction","071522", "NULL","percent.mt_RNA","MTpercentage_RNA/");
FeaturePlot_markers.function(coembed_multiome,"RNA",12,8,paste0(name,"_multiomeEmbed"),"071522", "NULL",bcells.markers,"Bcells/","bcells","no");
FeaturePlot_markers.function(coembed_multiome,"RNA",12,8,paste0(name,"_multiomeEmbed"),date, "NULL",acells.markers,"Acells/","acells","no");
FeaturePlot_markers.function(coembed_multiome,"RNA",6,6,paste0(name,"_multiomeEmbed"),date, "NULL",delta.markers,"Dcells/","dcells","no");
bcells.markers <- c("Ins1", "Ins2","Hadh","G6pc2");
acells.markers <- c("Gcg", "Ttr", "Gc");
delta.markers <- c("Sst");


###################################################################################################
##STEP10A Identification of clusters
###################################################################################################
#using res at 0.05 for multiome, crnaseq to 0.1 res, and atac to res=0.05, storing these resolutions 
# in metadata as finalClustering


rnaList_subset$final_cluster <- rnaList_subset$RNA_snn_res.0.1; 
#coembed_multiome$final_cluster <- coembed_multiome$ 

#to label with cellType info each cluster
###############################################
Idents(rnaList_subset) <- "final_cluster";
#Idents(coembed_multiome) <- "final_cluster"


rnaList_subset  <- RenameIdents(
  object = rnaList_subset,
  '0' = 'Beta_c0',
  '1' = 'Unkn_c1',
  '2' = 'Delta_C2',
  '3' = 'Alpha_c3',
  '4' = 'Unkn_c4',
  '5' = 'Alpha_c5',
  '6' = 'Beta_c6',
  '7' = 'Alpha_c7',
  '8' = 'Alpha_c8',
  '9' = 'Unkn_c9');
rnaList_subset$cellType <- Idents(rnaList_subset);

#adding info of samples metadata to perform the differential analysis between conditions
rnaList_subset$condition <- gsub("_.*","",gsub("RNAseq_","",rnaList_subset$orig.ident));
unique(rnaList_subset$condition); # "Mut"  "Control"
rnaList_subset$cluster_condition<- paste0(rnaList_subset$cellType, "__",rnaList_subset$condition);

Idents(rnaList_subset) <- "final_cluster"

#######################
###################################################################################################
##STEP10b  transfer identity of clusters from scRNASeq to scATACSeq to work  with scATAC not embedded
###################################################################################################
#to be performed depending on the markers and res cluster chosen
load(file=paste0(wd, f_results, f_Rdata,f_preIntegration,"transfer_anchors_step8.RData")) #transfer.anchors,

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rnaList_subset$cellType,
    weight.reduction = atactList_subset[["lsi"]], dims = c(2,4:30));

atactList_subset <- AddMetaData(atactList_subset, metadata = celltype.predictions)


atactList_subset$condition <- gsub("_.*","",gsub("ATACseq_","",atactList_subset$orig.ident));
atactList_subset$cluster_condition<- paste0(atactList_subset$predicted.id, "__",atactList_subset$condition);
atactList_subset$cluster_condition2<- paste0(atactList_subset$final_cluster, "__",atactList_subset$condition);

save(rnaList_subset, atactList_subset,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_step10.RData"))


#check the clusters predicted using rnaseq data criceicri
##################################################################
resolutionVector <- c(0.01,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.4,1.5,1.6,1.8,2,2.4,3,3.5); 

plot_perCluster <- DimPlot(object = atactList_subset, group.by=grep("_res",colnames(combined.samp.list@meta.data),value = TRUE), ncol=2 , pt.size=3.0, reduction = reductionName, label = T);
findClusters<- FindClusters(atactList_subset, resolution = 0.05); #can eb done in only one step
plot<- DimPlot(findClusters, reduction = reductionName, label = TRUE, raster = FALSE) +
labs(title = paste0("Clustering with ",resolutionVector[i], " Resolution"));

pdf(file=paste0("Clustering_perResolution_",nameOutput,"_UMAP_",date,".pdf"), height = height, width = width);
print(plot_perCluster);
dev.off();


findClusters.function(atactList_subset,"ATAC","atactList_subset", resolutionVector,paste0(name,"_atac_"),8,12,date,"umap.atac");
sapply(grep("_res",colnames(atactList_subset@meta.data),value = TRUE),
       function(x) length(unique(atactList_subset@meta.data[,x])))

FeaturePlot_FinalClustering.function(atactList_subset,8,6,"MtDysfunction",date, "NULL","percent.mt_ATAC","MTpercentage_RNA/");
FeaturePlot_markers.function(atactList_subset,"RNA",12,8,paste0(name,"_atac_lessdim"),date, "NULL",bcells.markers,"Bcells/","bcells","no");
FeaturePlot_markers.function(atactList_subset,"RNA",12,8,paste0(name,"_atac_lessdim"),date, "NULL",acells.markers,"Acells/","acells","no");
FeaturePlot_markers.function(atactList_subset,"RNA",6,6,paste0(name,"_atac_lessdim"),date, "NULL",delta.markers,"Dcells/","dcells","no");





#violin plot:visualize the co-embedding of both datasets
coembedViolin<- DimPlot(coembed_multiome, group.by = c("orig.ident", "seurat_annotations"))
png(file=paste0(wd,f_results,f_markers,f_plots,folderName,featureName,"_umap_",nameOutput,".png"), height = height, width = width, units = "in", res = 1200);
print(coembedViolin);
dev.off();


######################################################################################################################
# STEP10C: checking if the cell type label predictions in atac based in rnaseq converge 
#          with the labels we would have use if we didnt count with rnaseq
######################################################################################################################
#NOTE: prediction.score.max  quantifies the uncertainty associated with the predicted annotations. 
# cells  correctly annotated are  associated with high prediction scores (>90%), 
# cells incorrectly annotated have with  lower scores (<50%). 
#Incorrect annotation  tend to reflect closely cell types (i.e. Intermediate vs. Naive B cells).

#to really be able to do this plots the labels in predicted vs real throuth should be the same, in this case is not.
# because we didnt wanted to annotate the atac by itself
######################################################################################################################
predictions <- table(atactList_subset$final_cluster, atactList_subset$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

correct <- length(which(atactList_subset$final_cluster == atactList_subset$predicted.id))
incorrect <- length(which(atactList_subset$final_cluster != atactList_subset$predicted.id))
atactList_subset$annotation_correct <- atactList_subset$predicted.id == atactList_subset$final_cluster
pred1 <- DimPlot(atactList_subset, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
pred2 <- DimPlot(atactList_subset, group.by = "final_cluster", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
data <- FetchData(atactList_subset, vars = c("prediction.score.max", "annotation_correct"))


#plots
#A) umap comparison
png(file=paste0(wd,f_results,f_foundClusters,f_plots,"prediction_based_rnaseq_vs_atac_clusters_umap.png"), height = 8, width = 12, units = "in", res = 1200);
print(pred1 | pred2);
dev.off();

#B) heatmap how close is the prediction based on rnaseq compared to the clustering of just atacseq by percentage of cells
# library(cowplot);
pheatmap.prediction <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdensityPlot <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
    geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")

png(file=paste0(wd,f_results,f_foundClusters,f_plots,"prediction_based_rnaseq_vs_atac_clusters_heatmap.png"), height = 10, width = 16, units = "in", res = 1200);
print(pheatmap.prediction + pdensityPlot);
dev.off();
######################################################################################################################

Idents(atactList_subset) <- "final_cluster"
save(rnaList_subset, atactList_subset,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_step10.RData"))
#save.image(file=paste0(wd,f_results, f_Rdata,"session_multiome_step10.RData"))



###################################################################################################
##STEP11 DEA/DAA: differential expression/accessibility analyses
###################################################################################################

##########
# DEA -->
##########

library("Seurat");library("dplyr");library("cowplot");library("ggplot2");library("MAST");
library("openxlsx");library("grid");library("gridBase");library("gridExtra");library("patchwork");



#Before begining:
#rm(list =ls()) ## erasing all the enviroment variables
set.seed(22); # need to set it to get always the same random results and plots
#sessionInfo()

#wd:
#####################################################
wd<- "/home/omicData/multiome/";
setwd(wd);

#set output directories
f_results <- "/results/";
f_Rdata <- "RData/";
f_DEA <- "DEA/"
f_preIntegration <-"PreIntegration/";
#############################

dir.create(file.path(wd, f_results), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_DEA ), showWarnings = F);
dir.create(file.path(wd, f_results,f_DEA ,f_Rdata), showWarnings = F);
####################################################################
####################################################################


load(file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_step10.RData")); #rnaList_subset, atactList_subset
date <- date


####################################################################
####################################################################


vector_cluster.condition <- unique(rnaList_subset$cluster_condition);
vector_cluster.condition <- vector_cluster.condition[order(vector_cluster.condition )];
vector_cluster.condition_U <- vector_cluster.condition[grep("__Mut$",vector_cluster.condition)];
vector_cluster.condition_V <- vector_cluster.condition[grep("__Control$",vector_cluster.condition)];


DEA_MAST_perCondition.function <- function(experiment.clustersF,vector_cluster_ident1,vector_cluster_ident2,nameOutput,IdentName,date){

    markers.resList <- list();
    markers.EGsList <- list();
	Idents(experiment.clustersF) <- IdentName;

    for(i in 1:length(vector_cluster_ident1)) {

            overall.LA.EGs<- FindMarkers(experiment.clustersF, ident.1=vector_cluster_ident1[i] , ident.2 = vector_cluster_ident2[i] , test.use = "MAST",logfc.threshold =0);
            overall.LA.markers<- FindMarkers(experiment.clustersF, ident.1=vector_cluster_ident1[i] , ident.2 = vector_cluster_ident2[i] , test.use = "MAST");
            openxlsx::write.xlsx(overall.LA.EGs, file = paste0(wd,f_results,f_DEA,nameOutput,"_markers_",vector_cluster_ident1[i],"_",date, ".xlsx"),rowNames =TRUE, colNames =TRUE);

            markers.EGsList[[i]] <-overall.LA.EGs;
            markers.resList[[i]] <-overall.LA.markers;
            names(markers.EGsList)[i] <- vector_cluster_ident1[i];
            names(markers.resList)[i] <- vector_cluster_ident1[i];
            i<-i+1;
    };
    save(markers.EGsList, file=paste0(wd,f_results,f_DEA,f_Rdata,"markers.EGsList_",nameOutput,"_",date,".RData"));
    save(markers.resList, file=paste0(wd,f_results,f_DEA,f_Rdata,"markers.resList_",nameOutput,"_",date,".RData"));

};
DEA_MAST_perCondition.function(rnaList_subset,vector_cluster.condition_U,vector_cluster.condition_V,"Mut_vs_Control","cluster_condition",date);


##########
#<-- DEA 
##########



##########
# DAA  -->
##########

library("Seurat");library("dplyr");library("cowplot");library("ggplot2");library("MAST");
library("openxlsx");library("grid");library("gridBase");library("gridExtra");library("patchwork");



#Before begining:
#rm(list =ls()) ## erasing all the enviroment variables
set.seed(22); # need to set it to get always the same random results and plots
#sessionInfo()

#wd:
#####################################################
wd<- "/home/omicData/multiome/";
setwd(wd);


#set output directories
f_results <- "/results/";
f_Rdata <- "RData/";
f_DAA <- "DAA/"
f_preIntegration <-"PreIntegration/";
#############################

dir.create(file.path(wd, f_results), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_DAA ), showWarnings = F);
dir.create(file.path(wd, f_results,f_DAA ,f_Rdata), showWarnings = F);
####################################################################
####################################################################

load(file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_step10.RData")); #rnaList_subset, atactList_subset
date <- "072622"


####################################################################
####################################################################



DAA_MAST_perCondition.function <- function(experiment.clustersF,vector_cluster_ident1,vector_cluster_ident2,nameOutput,IdentsName,date){

    accesible.DAList <- list();
    accesible.AList <- list();

    for(i in 1:length(vector_cluster_ident1)) {

		DefaultAssay(experiment.clustersF) <- 'peaks'
		Idents(experiment.clustersF) <- IdentsName

		da_peaks <- FindMarkers(
		  object = experiment.clustersF,
		  ident.1 = vector_cluster_ident1[i],
		  ident.2 = vector_cluster_ident2[i],
		  min.pct = 0.05,
		  test.use = 'LR',
		  latent.vars = 'atac_peak_region_fragments');
		
		a_peaks <- FindMarkers(
		  object = experiment.clustersF,
		  ident.1 = vector_cluster_ident1[i],
		  ident.2 = vector_cluster_ident2[i],
		  min.pct = 0.05,
		  test.use = 'LR',
		  latent.vars = 'atac_peak_region_fragments',logfc.threshold =0);

        openxlsx::write.xlsx(a_peaks, file = paste0(wd,f_results,f_DAA,nameOutput,"_peaks_",vector_cluster_ident1[i],"_",date, ".xlsx"),rowNames =TRUE, colNames =TRUE);

            accesible.AList[[i]] <-a_peaks;
            accesible.DAList[[i]] <-da_peaks;
            names(accesible.AList)[i] <- vector_cluster_ident1[i];
            names(accesible.DAList)[i] <- vector_cluster_ident1[i];
            i<-i+1;
    };
    save(accesible.AList, file=paste0(wd,f_results,f_DAA,f_Rdata,"accesible.AList_",nameOutput,"_",date,".RData"));
    save(accesible.DAList, file=paste0(wd,f_results,f_DAA,f_Rdata,"accesible.DAList_",nameOutput,"_",date,".RData"));

};

#1) DAA based on cell types predicted from scRNASeq
#############################################################3
vector_cluster.condition <- unique(atactList_subset$cluster_condition);
vector_cluster.condition <- vector_cluster.condition[order(vector_cluster.condition )];
vector_cluster.condition_U <- vector_cluster.condition[grep("__Mut$",vector_cluster.condition)];
vector_cluster.condition_V <- vector_cluster.condition[grep("__Control$",vector_cluster.condition)];


DAA_MAST_perCondition.function(atactList_subset,vector_cluster.condition_U,vector_cluster.condition_V,"Mut_vs_Control_per_clltype","cluster_condition",date);

#2) DAA based on res cluster 0.05
#############################################################3
vector_cluster.condition <- unique(atactList_subset$cluster_condition2);
vector_cluster.condition <- vector_cluster.condition[order(vector_cluster.condition )];
vector_cluster.condition_U <- vector_cluster.condition[grep("__Mut$",vector_cluster.condition)];
vector_cluster.condition_V <- vector_cluster.condition[grep("__Control$",vector_cluster.condition)];

DAA_MAST_perCondition.function(atactList_subset,vector_cluster.condition_U,vector_cluster.condition_V,"Mut_vs_Control_perclstr","cluster_condition2",date);


##########
#<-- DAA 
##########



###################################################################################################
##STEP11b calculate log2FC per pair of mut and control samples
###################################################################################################
load(file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_step10.RData")); #rnaList_subset, atactList_subset
date <- "072622"

####################################################################
######## calculate log2FC per pair of mut and control samples ######
############ plus a std of all the log2FC calculations #############
####################################################################

pair1 <- subset(rnaList_subset, orig.ident ==c("RNAseq_Mut_1","RNAseq_Control_1"));
pair2 <- subset(rnaList_subset, orig.ident ==c("RNAseq_Mut_2","RNAseq_Control_2"));
pair3 <- subset(rnaList_subset, orig.ident ==c("RNAseq_Mut_2","RNAseq_Control_1"));
pair4 <- subset(rnaList_subset, orig.ident ==c("RNAseq_Mut_1","RNAseq_Control_2"));

vector_cluster.condition <- unique(rnaList_subset$cluster_condition);
vector_cluster.condition <- vector_cluster.condition[order(vector_cluster.condition )];
vector_cluster.condition_U <- vector_cluster.condition[grep("__Mut$",vector_cluster.condition)];
vector_cluster.condition_V <- vector_cluster.condition[grep("__Control$",vector_cluster.condition)];


DEA_MAST_perCondition_perPair.function <- function(experiment.clustersF,vector_cluster_ident1,vector_cluster_ident2,nameOutput,IdentName,date,namePair){

    markers.EGsList <- list();
	Idents(experiment.clustersF) <- IdentName;

    for(i in 1:length(vector_cluster_ident1)) {

            overall.LA.EGs<- FindMarkers(experiment.clustersF, ident.1=vector_cluster_ident1[i] , ident.2 = vector_cluster_ident2[i] , test.use = "MAST",logfc.threshold =0);

            markers.EGsList[[i]] <-overall.LA.EGs;
            names(markers.EGsList)[i] <- vector_cluster_ident1[i];
            i<-i+1;
    };

    for (ii in 1:(length(markers.EGsList))){
        if(ii==1){
          df <- markers.EGsList[[ii]];
      	  df$cluster <- gsub("__.*","",names(markers.EGsList)[ii]);
          EGs=df;
        } else{
          df.tmp <- markers.EGsList[[ii]];
      	  df.tmp$cluster <- gsub("__.*","",names(markers.EGsList)[ii]);
          EGs<- bind_rows(EGs,df.tmp);
        };
        ii <- ii+1;
    };
    colnames(EGs) <-  paste0(namePair,"_",  colnames(EGs));
    EGs$Gene <- gsub("\\...*","",rownames(EGs));
    rownames(EGs)<- NULL;
    colnames(EGs)[6]<- "cluster"
    assign(paste0("EGs_",namePair),EGs,.GlobalEnv);

};
DEA_MAST_perCondition_perPair.function(pair1,vector_cluster.condition_U,vector_cluster.condition_V,"pair1","cluster_condition",date,"pair1");
DEA_MAST_perCondition_perPair.function(pair2,vector_cluster.condition_U,vector_cluster.condition_V,"pair2","cluster_condition",date,"pair2");
DEA_MAST_perCondition_perPair.function(pair3,vector_cluster.condition_U,vector_cluster.condition_V,"pair3","cluster_condition",date,"pair3");
DEA_MAST_perCondition_perPair.function(pair4,vector_cluster.condition_U,vector_cluster.condition_V,"pair4","cluster_condition",date,"pair4");

allEGs <-mget(ls(pattern="EGs_"));
allEGs_df<- full_join(allEGs[[1]],allEGs[[2]],by=c("Gene","cluster"))
allEGs_df<- full_join(allEGs_df,allEGs[[3]],by=c("Gene","cluster"))
allEGs_df<- full_join(allEGs_df,allEGs[[4]],by=c("Gene","cluster"))

save(allEGs, file=paste0(wd, f_results,f_DEA, f_Rdata,"EGs_perPair_log2fc_calculation.RData"));
#library(matrixStats)
allEGs_df$log2FC_std = rowSds(as.matrix(allEGs_df[,c("pair1_avg_log2FC","pair2_avg_log2FC","pair3_avg_log2FC","pair4_avg_log2FC")]))
openxlsx::write.xlsx(allEGs_df, file = paste0(wd,f_results,f_DEA,"EGs_perPair_log2fc_calculation.xlsx"),rowNames =FALSE, colNames =TRUE);

###################################################################################################
##STEP12 Peaks annotation with the most proximal gene
###################################################################################################


library("Seurat");library("dplyr");library("cowplot");library("ggplot2");library("MAST");library("Signac");
library("openxlsx");library("grid");library("gridBase");library("gridExtra");library("patchwork");



#Before begining:
#rm(list =ls()) ## erasing all the enviroment variables
set.seed(22); # need to set it to get always the same random results and plots
#sessionInfo()

#wd:
#####################################################
wd<- "/home/omicData/multiome/";
setwd(wd);


#set output directories
f_results <- "/results/";
f_Rdata <- "RData/";
f_DAA <- "DAA/"
f_preIntegration <-"PreIntegration/";

load(file=paste0(wd,f_results,f_DAA,f_Rdata,"accesible.AList_Mut_vs_Control_per_clltype_",date,".RData")); #accesible.AList
load(file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_step10.RData")); #rnaList_subset, atactList_subset
date <- "072622"

for (ii in 1:(length(accesible.AList))){
    if(ii==1){
      df <- accesible.AList[[ii]];
      df$cluster <- gsub("__.*","",names(accesible.AList)[ii]);
      accesPeaks=df;
    } else{
      df.tmp <- accesible.AList[[ii]];
      df.tmp$cluster <- gsub("__.*","",names(accesible.AList)[ii]);
      accesPeaks<- bind_rows(accesPeaks,df.tmp);
    };
    ii <- ii+1;
};

    accesPeaks$PeakID <- gsub("\\...*","",rownames(accesPeaks));
    rownames(accesPeaks)<- NULL;
    accesPeaks$PeakID2 <- accesPeaks$PeakID; 
    accesPeaks<-tidyr::separate(accesPeaks, PeakID2, c("Peak_chr" ,"Peak_start","Peak_end"  ),sep="-");
    accesPeaks$Peak_size_bps <- as.numeric(accesPeaks$Peak_end) - as.numeric(accesPeaks$Peak_start);

closest_genes <- ClosestFeature(atactList_subset, regions = accesPeaks$PeakID)
colnames(closest_genes)[7] <-"PeakID";
accesPeaksAnnot <- left_join(accesPeaks,closest_genes,by="PeakID")

openxlsx::write.xlsx(accesPeaksAnnot, file = paste0(wd,f_results,f_DAA, "Mut_vs_Control_peaks_annotated_all_clusters_per_clltype.xlsx"),rowNames =FALSE, colNames =TRUE);

###################################################################################################
##STEP13 subsetting atac,rnaseq and multiome objects to see the cell types we are interested on 
###################################################################################################
# A) Annotate and subset ATACseq: selecting alpha, beta & delta cells
#########################################################################
atactList_subset$cluster_orig.ident <- paste0(gsub("__.*","",atactList_subset$cluster_condition),"__", gsub("ATACseq_","",atactList_subset$orig.ident));
Idents(atactList_subset) <- "predicted.id"
atac_selected_CellTypes <-subset(x = atactList_subset, idents = c("Beta_c0","Delta_C2", "Alpha_c3","Alpha_c5","Beta_c6","Alpha_c7","Alpha_c8"), invert = FALSE)

# B) Annotate and subset RNAseq: selecting alpha, beta & delta cells
#########################################################################
rnaList_subset$cluster_orig.ident <- paste0(gsub("__.*","",rnaList_subset$cluster_condition),"__", gsub("RNAseq_","",rnaList_subset$orig.ident));
Idents(rnaList_subset) <- "cellType"
rna_selected_CellTypes <-subset(x = rnaList_subset, idents = c("Beta_c0","Delta_C2", "Alpha_c3","Alpha_c5","Beta_c6","Alpha_c7","Alpha_c8"), invert = FALSE)

save(rnaList_subset, atactList_subset,
	file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_step13.RData")); 

# C) Annotate and subset coembed multiome: selecting alpha, beta & delta cells
#########################################################################
#multiome #as we generated the multiome before annotating rnaseq now annotating

load( file=paste0(wd, f_results, f_Rdata,f_after_enbedding,"coembedded_multiome_step9.RData")); #coembed_multiome

coembed_multiome$final_cluster <- coembed_multiome$RNA_snn_res.0.1; 
Idents(coembed_multiome) <- "final_cluster";

coembed_multiome  <- RenameIdents(
  object = coembed_multiome,
  '0' = 'Beta_c0',
  '1' = 'Unkn_c1',
  '2' = 'Delta_C2',
  '3' = 'Alpha_c3',
  '4' = 'Unkn_c4',
  '5' = 'Alpha_c5',
  '6' = 'Beta_c6',
  '7' = 'Alpha_c7',
  '8' = 'Alpha_c8',
  '9' = 'Unkn_c9');
coembed_multiome$cellType <- Idents(coembed_multiome);

#adding info of samples metadata to perform the differential analysis between conditions
coembed_multiome$condition <- gsub("_.*","",gsub("RNAseq_","",gsub("ATACseq_","",coembed_multiome$orig.ident)));
unique(coembed_multiome$condition); # "Mut"  "Control"
coembed_multiome$cluster_condition<- paste0(coembed_multiome$cellType, "__",coembed_multiome$condition);
unique(coembed_multiome$cluster_condition);

coembed_multiome$cluster_orig.ident <- paste0(gsub("__.*","",coembed_multiome$cluster_condition),"__", gsub("RNAseq_","",gsub("ATACseq_","",coembed_multiome$orig.ident)));
Idents(coembed_multiome) <- "cellType"
coembed_multiome_selected_CellTypes <-subset(x = coembed_multiome, idents = c("Beta_c0","Delta_C2", "Alpha_c3","Alpha_c5","Beta_c6","Alpha_c7","Alpha_c8"), invert = FALSE)

save(atac_selected_CellTypes, rna_selected_CellTypes,file=paste0(wd, f_results, f_Rdata,f_preIntegration,"multiome_processing_step13_rnaseq_atacSeq_subset.RData"));
save(coembed_multiome_selected_CellTypes,file=paste0(wd, f_results, f_Rdata,f_after_enbedding,"multiome_processing_step13_coembed_multiome_subset.RData"));
save(coembed_multiome, file=paste0(wd, f_results, f_Rdata,f_after_enbedding,"coembedded_multiome_step13.RData")); #coembed_multiome

# D) redoing UMAPS plots with the selected cell clusters
#########################################################################



bcells.markers <- c("Ins1", "Ins2","Hadh","G6pc2"); 
acells.markers <- c("Gcg", "Ttr", "Gc");
delta.markers <- c("Sst");

FeaturePlot_FinalClustering.function(rna_selected_CellTypes,8,6,"MtDysfunction_rna_selCellTypes",date, "NULL","percent.mt_RNA","MTpercentage_RNA/");
FeaturePlot_markers.function(rna_selected_CellTypes,"RNA",9,7,paste0(name,"_rna_selCellTypes"),date, "NULL",bcells.markers,"Bcells/","bcells","no");
FeaturePlot_markers.function(rna_selected_CellTypes,"RNA",9,7,paste0(name,"_rna_selCellTypes"),date, "NULL",acells.markers,"Acells/","acells","no");
FeaturePlot_markers.function(rna_selected_CellTypes,"RNA",9,7,paste0(name,"_rna_selCellTypes"),date, "NULL",delta.markers,"Dcells/","dcells","no");

FeaturePlot_FinalClustering.function(atac_selected_CellTypes,8,6,"MtDysfunction_atac_selCellType",date, "NULL","percent.mt_ATAC","MTpercentage_RNA/");
FeaturePlot_markers.function(atac_selected_CellTypes,"RNA",9,7,paste0(name,"_atac_selCellTypes"),date, "NULL",bcells.markers,"Bcells/","bcells_atac_RNAimput","no");
FeaturePlot_markers.function(atac_selected_CellTypes,"RNA",9,7,paste0(name,"_atac_selCellTypes"),date, "NULL",acells.markers,"Acells/","acells_atac_RNAimput","no");
FeaturePlot_markers.function(atac_selected_CellTypes,"RNA",9,7,paste0(name,"_atac_selCellTypes"),date, "NULL",delta.markers,"Dcells/","dcells_atac_RNAimput","no");


plot_perCluster <- DimPlot(object = rna_selected_CellTypes, group.by="cellType", ncol=2 , pt.size=3.0, reduction = "umap", label = F);
pdf(file=paste0(wd,f_results,f_foundClusters,f_plots,"Clustering_perResolution_",paste0(name,"_rna_selCellTypes"),"_UMAP_",date,".pdf"), height = 8, width = 12);
print(plot_perCluster);
dev.off();


plot_perClusteratac <- DimPlot(object = atac_selected_CellTypes, group.by="ATAC_snn_res.0.05", ncol=2 , pt.size=3.0, reduction = "umap.atac", label = F);
pdf(file=paste0(wd,f_results,f_foundClusters,f_plots,"Clustering_perResolution_",paste0(name,"_atac_selCellClusters"),"_UMAP_",date,".pdf"), height = 8, width = 12);
print(plot_perClusteratac);
dev.off();


plot_perClusteratac2 <- DimPlot(object = atac_selected_CellTypes, group.by="predicted.id", ncol=2 , pt.size=3.0, reduction = "umap.atac", label = F);
pdf(file=paste0(wd,f_results,f_foundClusters,f_plots,"Clustering_perResolution_",paste0(name,"_atac_selCellTypes"),"_UMAP_",date,".pdf"), height = 8, width = 12);
print(plot_perClusteratac2);
dev.off();



plot_perCluster_multiome <- DimPlot(object = coembed_multiome_selected_CellTypes, group.by="cellType", ncol=2 , pt.size=3.0, reduction = "umap.multiome.embed", label = F);
pdf(file=paste0(wd,f_results,f_foundClusters,f_plots,"Clustering_perResolution_",paste0(name,"_multiome_selCellTypes"),"_UMAP_",date,".pdf"), height = 8, width = 12);
print(plot_perCluster_multiome);
dev.off();


FeaturePlot_FinalClustering.function(coembed_multiome_selected_CellTypes,8,6,"MtDysfunction_multiome_selCellTypes",date, "NULL","percent.mt_RNA","MTpercentage_RNA/");
FeaturePlot_markers.function(coembed_multiome_selected_CellTypes,"RNA",12,8,paste0(name,"_multiomeEmbed_selCellTypes"),date, "NULL",bcells.markers,"Bcells/","bcells","no");
FeaturePlot_markers.function(coembed_multiome_selected_CellTypes,"RNA",12,8,paste0(name,"_multiomeEmbed_selCellTypes"),date, "NULL",acells.markers,"Acells/","acells","no");
FeaturePlot_markers.function(coembed_multiome_selected_CellTypes,"RNA",6,6,paste0(name,"_multiomeEmbed_selCellTypes"),date, "NULL",delta.markers,"Dcells/","dcells","no");


# E) plot RNA and geneActivity data side by side
#########################################################################
load(file=paste0(wd, f_results, f_Rdata,f_after_enbedding,"multiome_processing_step13_coembed_multiome_subset.RData")); #coembed_multiome_selected_CellTypes,
load(file=paste0(wd, f_results, f_Rdata,f_after_enbedding,"coembedded_multiome_step13.RData")); #coembed_multiome_selected_CellTypes,

coembed_multiome.TMP<-coembed_multiome_selected_CellTypes
coembed_multiome.TMP<-coembed_multiome
rnaList_subset2<- rnaList_subset
rnaList_subset2 <- FindVariableFeatures(rnaList_subset2, selection.method = "vst", nfeatures = 2000)
DefaultAssay(coembed_multiome.TMP) <- "peaks"
atactList_subset2<-atactList_subset
gene.activities <- GeneActivity(atactList_subset2, features = VariableFeatures(rnaList_subset2)) 

atactList_subset2[['geneActivity']] <- CreateAssayObject(counts = gene.activities);


DefaultAssay(coembed_multiome.TMP) <- "RNA"
coembed_multiome.TMP <- NormalizeData(coembed_multiome.TMP, normalization.method = "LogNormalize", margin = 2)
DefaultAssay(coembed_multiome.TMP) <- "geneActivity"

# Note that the following command is an alternative but returns the same result
coembed_multiome.TMP <- NormalizeData(coembed_multiome.TMP, normalization.method = "LogNormalize", margin = 2, assay = "RNA")
# Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(coembed_multiome.TMP) <- "RNA"
p1 <- FeaturePlot(coembed_multiome.TMP, "Ins2", cols = c("lightgrey", "red")) + ggtitle("Ins2 RNA")
DefaultAssay(coembed_multiome.TMP) <- "geneActivity"
p2 <- FeaturePlot(coembed_multiome.TMP, "Ins2") + ggtitle("Ins2 geneActivity")

png(file=paste0(wd,f_results,f_markers,f_plots,"RNa_vs_GeneActivity_Ins2_subsetCells.png"), res = 600, height = 7, width = 7,units = "in", ); #
print(p1 | p2);
dev.off();
#########################################################################



###################################################################################################
##STEP14 Identify and plot coverage of genes and regions of interest fro ATAC data 
###################################################################################################
#We can visualize these links using the CoveragePlot() function, 
# or alternatively we could use the CoverageBrowser() function in an interactive analysis:


# select peak close to a gene of interest: Example
Peaks <- unique(accesPeaksAnnot[which(accesPeaksAnnot$gene_name=="Mut"),"PeakID"]); #3  
#select a peak close to a gene of interest by genomic position
Peaks <- unique(accesPeaksAnnot[grep("chr19-296",accesPeaksAnnot$PeakID),]); ##position of Mut gene Chr13:30724295-30729299 
# set plotting order
orderPlot <- c( paste0(c("Beta_c0","Beta_c6","Alpha_c3", "Alpha_c5" ,"Alpha_c7", "Alpha_c8", "Delta_C2"),"__Mut"),
				paste0(c("Beta_c0","Beta_c6","Alpha_c3", "Alpha_c5" ,"Alpha_c7", "Alpha_c8", "Delta_C2"),"__Control"))


orderPlot_perSample <- c( paste0(c("Beta_c0","Beta_c6","Alpha_c3", "Alpha_c5" ,"Alpha_c7", "Alpha_c8", "Delta_C2"),"__Mut_1"),
				paste0(c("Beta_c0","Beta_c6","Alpha_c3", "Alpha_c5" ,"Alpha_c7", "Alpha_c8", "Delta_C2"),"__Mut_2"),
				paste0(c("Beta_c0","Beta_c6","Alpha_c3", "Alpha_c5" ,"Alpha_c7", "Alpha_c8", "Delta_C2"),"__Control_1"),
				paste0(c("Beta_c0","Beta_c6","Alpha_c3", "Alpha_c5" ,"Alpha_c7", "Alpha_c8", "Delta_C2"),"__Control_2"))

coveragePlot.function <- function(atacObject, IdentsName,orderPlot,extend.upstream,extend.downstream,nameOutput,PeaksID,width,height){

	DefaultAssay(atacObject) <- 'peaks';
	Idents(atacObject) <- IdentsName;

	# first compute the GC content for each peak
	###############################################
	#atacObject <- RegionStats(atacObject, genome = BSgenome.Mmusculus.UCSC.mm10); #this annotation is not up to date with the one using by UCSC so gives error
	
	# link peaks to genes #if not one before
	###############################################

	#atacObject <- LinkPeaks(
	#object = atacObject,
	#peak.assay = "peaks",
	#expression.assay = "RNA",
	#genes.use = PeaksID$gene_name);


	levels(atacObject) <- orderPlot;
	coverPlotList <- list();
	
	for ( i in 1:length(PeaksID)){
		coverPlotList[[i]] <- CoveragePlot(
	    object = atacObject,
	    region = PeaksID[i],
	    extend.upstream = extend.upstream,
	    extend.downstream = extend.downstream);
		i <- 1+1;
	};

	pdf(file=paste0(wd,f_results,f_atacCoverage,f_plots,nameOutput,"_cov_",extend.upstream,"_",extend.downstream,".pdf"), width = width, height = height,onefile=TRUE)
	print(coverPlotList);
	dev.off();
    assign(paste0("coverPlotList_",nameOutput),coverPlotList,.GlobalEnv);
    #p_ensemble<- patchwork::wrap_plots(coverPlotList[[1]], coverPlotList[[2]],coverPlotList[[3]], ncol = 1)

};

coveragePlot.function(atac_selected_CellTypes,"cluster_condition",orderPlot,15000,15000,"Pksnear_Mut_pcond",Peaks,8,12);
coveragePlot.function(atac_selected_CellTypes,"cluster_condition",orderPlot,5000,5000,"Pksnear_Mut_pcond",Peaks,8,12);
coveragePlot.function(atac_selected_CellTypes,"cluster_orig.ident",orderPlot_perSample,20000,20000,"Pksnear_Mut_psample",Peaks,8,24);
coveragePlot.function(atac_selected_CellTypes,"cluster_orig.ident",orderPlot_perSample,15000,15000,"Pksnear_Mut_psample",Peaks,8,24);
coveragePlot.function(atac_selected_CellTypes,"cluster_orig.ident",orderPlot_perSample,5000,5000,"Pksnear_Mut_psample",Peaks,8,24);

###################################################################################################
##STEP15 Intersection DEGs and DAGs  
###################################################################################################
#####################################
 library(readxl)
 accesPeaksAnnot<- as.data.frame(read_excel(paste0(wd,f_results,f_DAA, "Mut_vs_Control_peaks_annotated_all_clusters_per_clltype.xlsx"),sheet =1, col_names =T)); 
 load( file=paste0(wd,f_results,f_DEA,f_Rdata,"markers.EGsList_Mut_vs_ControldateRData"));#markers.EGsList

############################################
#1) list per cluster pval thressholds
############################################
DAGslist <-list();
DEGslist<- list();
DAGsReglist <-list();
DEGsReglist<- list();
for ( i in 1:length(unique(accesPeaksAnnot$cluster))){
	accesPeaksAnnot$Regulation <- ifelse(accesPeaksAnnot$avg_log2FC>0, "UpRegulated","DownRegulated")
	accesPeaksAnnot<- accesPeaksAnnot[which(accesPeaksAnnot$p_val<0.05),];
	DAGsReglist[[i]] <-unique(accesPeaksAnnot[which(accesPeaksAnnot$cluster==unique(accesPeaksAnnot$cluster)[i]),c("gene_name","Regulation","cluster")]);
	DAGslist[[i]] <-unique(accesPeaksAnnot[which(accesPeaksAnnot$cluster==unique(accesPeaksAnnot$cluster)[i]),c("gene_name","cluster")]);
	names(DAGslist)[i] <-paste0("ATACSeq_",unique(accesPeaksAnnot$cluster)[i]);
	names(DAGsReglist)[i] <-paste0("ATACSeq_",unique(accesPeaksAnnot$cluster)[i]);
}


for ( i in 1:length(markers.EGsList)){
	EGs<- markers.EGsList[[i]];
	EGs$gene_name <- rownames(EGs)
	EGs$cluster <- names(markers.EGsList)[i];
	EGs$Regulation <- ifelse(EGs$avg_log2FC>0, "UpRegulated","DownRegulated")
	EGs<- EGs[which(EGs$p_val<0.05),];
	DEGsReglist[[i]] <-unique(EGs[,c("gene_name","Regulation","cluster")]);
	DEGslist[[i]] <-unique(EGs[,c("gene_name","cluster")]);
	names(DEGslist)[i] <-paste0("RNASeq_",gsub("__Mut","",names(markers.EGsList)[i]));
	names(DEGsReglist)[i] <-paste0("RNASeq_",gsub("__Mut","",names(markers.EGsList)[i]));
}

#combine both
degDagList <- c(DEGslist, DAGslist)   # Merge two lists
degDagRegList <- c(DEGsReglist, DAGsReglist)   # Merge two lists

names(degDagList)
names(degDagRegList)


############################################
#2) list per cluster pval thressholds plus log2fc in atac
############################################
DAGslist2 <-list();
DAGsReglist2 <-list();
for ( i in 1:length(unique(accesPeaksAnnot$cluster))){
	accesPeaksAnnot$Regulation <- ifelse(accesPeaksAnnot$avg_log2FC>0, "UpRegulated","DownRegulated")
	accesPeaksAnnot<- accesPeaksAnnot[which(accesPeaksAnnot$p_val<0.05),];
	accesPeaksAnnot1<- accesPeaksAnnot[which(accesPeaksAnnot$avg_log2FC< -1),];
	accesPeaksAnnot2<- accesPeaksAnnot[which(accesPeaksAnnot$avg_log2FC> 1),];
	accesPeaksAnnotF <- bind_rows(accesPeaksAnnot1,accesPeaksAnnot2);

	DAGsReglist2[[i]] <-unique(accesPeaksAnnotF[which(accesPeaksAnnotF$cluster==unique(accesPeaksAnnotF$cluster)[i]),c("gene_name","Regulation","cluster")]);
	DAGslist2[[i]] <-unique(accesPeaksAnnotF[which(accesPeaksAnnotF$cluster==unique(accesPeaksAnnotF$cluster)[i]),c("gene_name","cluster")]);
	names(DAGslist2)[i] <-paste0("ATACSeq_",unique(accesPeaksAnnot$cluster)[i]);
	names(DAGsReglist2)[i] <-paste0("ATACSeq_",unique(accesPeaksAnnot$cluster)[i]);
}


#combine both
degDagList2 <- c(DEGslist, DAGslist2)   # Merge two lists
degDagRegList2 <- c(DEGsReglist, DAGsReglist2)   # Merge two lists


############################################
#2) list per cluster pval adj thressholds 
############################################

DAGslist3 <-list();
DEGslist3<- list();
DAGsReglist3 <-list();
DEGsReglist3<- list();
for ( i in 1:length(unique(accesPeaksAnnot$cluster))){
	accesPeaksAnnot$Regulation <- ifelse(accesPeaksAnnot$avg_log2FC>0, "UpRegulated","DownRegulated")
	accesPeaksAnnot<- accesPeaksAnnot[which(accesPeaksAnnot$p_val_adj<0.05),];
	DAGsReglist3[[i]] <-unique(accesPeaksAnnot[which(accesPeaksAnnot$cluster==unique(accesPeaksAnnot$cluster)[i]),c("gene_name","Regulation","cluster")]);
	DAGslist3[[i]] <-unique(accesPeaksAnnot[which(accesPeaksAnnot$cluster==unique(accesPeaksAnnot$cluster)[i]),c("gene_name","cluster")]);
	names(DAGslist3)[i] <-paste0("ATACSeq_",unique(accesPeaksAnnot$cluster)[i]);
	names(DAGsReglist3)[i] <-paste0("ATACSeq_",unique(accesPeaksAnnot$cluster)[i]);
}


for ( i in 1:length(markers.EGsList)){
	EGs<- markers.EGsList[[i]];
	EGs$gene_name <- rownames(EGs)
	EGs$cluster <- names(markers.EGsList)[i];
	EGs$Regulation <- ifelse(EGs$avg_log2FC>0, "UpRegulated","DownRegulated")
	EGs<- EGs[which(EGs$p_val_adj<0.05),];
	DEGsReglist3[[i]] <-unique(EGs[,c("gene_name","Regulation","cluster")]);
	DEGslist3[[i]] <-unique(EGs[,c("gene_name","cluster")]);
	names(DEGslist3)[i] <-paste0("RNASeq_",gsub("__Mut","",names(markers.EGsList)[i]));
	names(DEGsReglist3)[i] <-paste0("RNASeq_",gsub("__Mut","",names(markers.EGsList)[i]));
}

#combine both
degDagList3 <- c(DEGslist3, DAGslist3)   # Merge two lists
degDagRegList3 <- c(DEGsReglist3, DAGsReglist3)   # Merge two lists


save(accesPeaksAnnot,markers.EGsList,degDagList,degDagRegList,degDagList2,degDagRegList2,degDagList3,degDagRegList3,
 file=paste0(wd, f_results, f_Rdata,f_after_enbedding,"DEGs_DAGs_lists_for_venn_step15.RData")); 
load( file=paste0(wd, f_results, f_Rdata,f_after_enbedding,"DEGs_DAGs_lists_for_venn_step15.RData")); 



#vennDiagram per cluster
###################################

model_colors <- data.frame(matrix(nrow=length(names(degDagList)), ncol=4))
colnames(model_colors) <- c("model","fill","cat.col","labels")
model_colors$model <- names(degDagList)


#color for the filling of the shapes "fill"
model_colors[grep("RNASeq_Beta_c0",model_colors$model), 'fill'] <- "darksalmon"
model_colors[grep("ATACSeq_Beta_c0",model_colors$model), 'fill'] <- "darkorchid2"
model_colors[grep("RNASeq_Beta_c6",model_colors$model), 'fill'] <- "azure3"
model_colors[grep("ATACSeq_Beta_c6",model_colors$model), 'fill'] <- "deeppink"
model_colors[grep("RNASeq_Alpha_c3",model_colors$model), 'fill'] <- "blue2"
model_colors[grep("ATACSeq_Alpha_c3",model_colors$model), 'fill'] <- "lightgoldenrod2"
model_colors[grep("RNASeq_Alpha_c5",model_colors$model), 'fill'] <- "aquamarine2"
model_colors[grep("ATACSeq_Alpha_c5",model_colors$model), 'fill'] <- "thistle1"
model_colors[grep("RNASeq_Alpha_c7",model_colors$model), 'fill'] <- "indianred1"
model_colors[grep("ATACSeq_Alpha_c7",model_colors$model), 'fill'] <- "navajowhite2"
model_colors[grep("RNASeq_Alpha_c8",model_colors$model), 'fill'] <- "dodgerblue1"
model_colors[grep("ATACSeq_Alpha_c8",model_colors$model), 'fill'] <- "darkolivegreen3"

#####col names: cat.col

model_colors[grep("RNASeq_Beta_c0",model_colors$model), 'cat.col'] <- "darksalmon"
model_colors[grep("ATACSeq_Beta_c0",model_colors$model), 'cat.col'] <- "darkorchid3"
model_colors[grep("RNASeq_Beta_c6",model_colors$model), 'cat.col'] <- "grey28"
model_colors[grep("ATACSeq_Beta_c6",model_colors$model), 'cat.col'] <- "deeppink"
model_colors[grep("RNASeq_Alpha_c3",model_colors$model), 'cat.col'] <- "blue2"
model_colors[grep("ATACSeq_Alpha_c3",model_colors$model), 'cat.col'] <- "lightgoldenrod4"
model_colors[grep("RNASeq_Alpha_c5",model_colors$model), 'cat.col'] <- "aquamarine4"
model_colors[grep("ATACSeq_Alpha_c5",model_colors$model), 'cat.col'] <- "thistle4"
model_colors[grep("RNASeq_Alpha_c7",model_colors$model), 'cat.col'] <- "indianred4"
model_colors[grep("ATACSeq_Alpha_c7",model_colors$model), 'cat.col'] <- "navajowhite4"
model_colors[grep("RNASeq_Alpha_c8",model_colors$model), 'cat.col'] <- "dodgerblue4"
model_colors[grep("ATACSeq_Alpha_c8",model_colors$model), 'cat.col'] <- "coral4"
model_colors$labels <- gsub("ATACSeq_","",gsub("RNASeq_"," ",names(degDagList)))



vennD2_homologs.function<- function(model_colors,model1,model2,label1,label2,DEGsList,DEGsListReg,categoryNames,modelRef){	
    #library("VennDiagram");
	category <-c(label1,label2);          
	f_category <- paste0(modelRef,"_",model1,"_",model2,"/");
	name_category <- paste0(modelRef,"_",model1,"_",model2);
	#creating the folders needed
	dir.create(file.path(getwd (), f_results, f_venn), showWarnings = F);
	dir.create(file.path(getwd (), f_results, f_venn,f_category), showWarnings = F);
	dir.create(file.path(getwd (), f_results, f_venn,f_category,f_vdata), showWarnings = F);
	dir.create(file.path(getwd (), f_results, f_venn,f_category,f_vplots), showWarnings = F);

	#overlaps

	file1 <- unique(DEGsList[[grep(model1,names(DEGsList))]][,"gene_name"]);
	file2 <- unique(DEGsList[[grep(model2,names(DEGsList))]][,"gene_name"]);
	area1 <- length(unique(unlist(file1)));
	area2 <- length(unique(unlist(file2)));
	n12 <- length(Reduce(intersect, list(t(file1),t(file2))));
	
	#saving data
	data <- data.frame(area1,area2,n12);
	colnames(data) <- c("area1", "area2", "Intercept_12");
	write.table(data, file = paste0(wd, f_results, f_venn,f_category,f_vdata,name_category, "_venn_sumup.txt"), quote = F, row.names = F, col.names = T, sep = "\t");
	#saving each intercept list of DEGs
	write.table(data.frame(area_1=unique(unlist(file1))), file = paste0(wd, f_results, f_venn,f_category,f_vdata,name_category, "_venn_data_intercept_list_area1.txt"), quote = F, row.names = F, col.names = T, sep = "\t");
	write.table(data.frame(area_2=unique(unlist(file2))), file = paste0(wd, f_results, f_venn,f_category,f_vdata,name_category, "_venn_data_intercept_list_area2.txt"), quote = F, row.names = F, col.names = T, sep = "\t");
	write.table(data.frame(n12=Reduce(intersect, list(t(file1),t(file2)))), file = paste0(wd, f_results, f_venn,f_category,f_vdata,name_category, "_venn_data_intercept_list_n12.txt"), quote = F, row.names = F, col.names = T, sep = "\t");
	
	#ploting the venn diagram
	fill<- 	c(model_colors[grep(model1,model_colors$model),2],model_colors[grep(model2,model_colors$model),2]);
	cat.col <-c(model_colors[grep(model1,model_colors$model),3],model_colors[grep(model2,model_colors$model),3]);
	
	grDevices::cairo_pdf(file= paste0(wd,f_results, f_venn,f_category,f_vplots, name_category,".pdf"),height =3, width = 3) 
	venn.plot1 <- draw.pairwise.venn(
	  area1 = area1,
	  area2 = area2,
	  cross.area = n12,
	  category = categoryNames,
	  fill = fill,
	  cat.col = cat.col,
	  margin = 0.3,
	  cex = 0.7,
	  alpha = 0.7,
	  cat.cex = 0.9,
	  cat.default.pos = "outer",
	  cat.just=list(c(0.2,-2.2), c(0.9,10.2)), 
	  ext.pos = 2.9,
	  ext.dist = -0.05,
	  ext.length = 0.85,
	  euler.d=FALSE,
	  scaled=FALSE,
	  ind=TRUE
	);
	assign(paste0("venn.plot1"),venn.plot1,.GlobalEnv);
	assign(paste0("n12",model1,model2),data.frame(n12=Reduce(intersect, list(t(file1),t(file2)))),.GlobalEnv);
	grid.draw(venn.plot1);
	grid.draw(venn.plot1);
	dev.off();

	grid.draw(venn.plot1);
	grid.draw(venn.plot1);
	dev.off();
	save(venn.plot1,file=paste0(wd,f_results, f_venn,f_category,f_vplots, name_category,".Rdata")); 
	
	#stats shared region of DEGs    
	n12Genes<- data.frame(n12=Reduce(intersect, list(t(file1),t(file2))));
	
	cond1 <-DEGsListReg[[grep(model1,names(DEGsListReg))]][,c(1,2)];
	cond1 <-as.data.frame(cond1[which(cond1[,1] %in% n12Genes[,1]),]);
	cond1.Up <- cond1[which(cond1[,2]=="UpRegulated"),1];  
	cond1.Down <- cond1[which(cond1[,2]=="DownRegulated"),1]; 
	
	cond2 <- DEGsListReg[[grep(model2,names(DEGsListReg))]][,c(1,2)];
	cond2 <-as.data.frame(cond2[which(cond2[,1] %in% n12Genes[,1]),]);
	cond2.Up <- cond2[which(cond2[,2]=="UpRegulated"),1];  
	cond2.Down <- cond2[which(cond2[,2]=="DownRegulated"),1];  
	

	#total number hits
	############################################
	total.Shared <- unique(cond1[,1][cond1[,1] %in% cond2[,1]]);
	ltotal.Shared <- length(unique(cond1[,1][cond1[,1] %in% cond2[,1]]));

	# shared, up in both
	############################################
	lshared.Up_Up <-length(unique(cond1.Up[cond1.Up %in% cond2.Up]));
	shared.Up_Up <-unique(cond1.Up[cond1.Up %in% cond2.Up]);

	#shared, down in one cond, up in the other
	############################################
	lshared.1Up_2Down<- length(unique(cond1.Up[cond1.Up %in% cond2.Down]));
	shared.1Up_2Down<- unique(cond1.Up[cond1.Up %in% cond2.Down]);
	
	lshared.2Up_1Down<- length(unique(cond2.Up[cond2.Up %in% cond1.Down]));
	shared.2Up_1Down<- unique(cond2.Up[cond2.Up %in% cond1.Down]);

	# shared, down in both
	############################################
	lshared.Down_Down<-length(unique(cond1.Down[cond1.Down %in% cond2.Down])); 
	shared.Down_Down<-unique(cond1.Down[cond1.Down %in% cond2.Down]); 

	assign(paste0(f_category),f_category,.GlobalEnv);
	cond1$model <- model1;
	colnames(cond1)[2] <- paste0(colnames(cond1)[2],".",model1);
	
	cond2$model <- model2
	colnames(cond2)[2] <- paste0(colnames(cond2)[2],".",model2);
	colnames(cond2)[1] <- "GeneName";
	colnames(cond1)[1] <- "GeneName";
	a <- full_join(cond1,cond2, by=colnames(cond1)[1]);
	a$fmodel <- paste0(a[,3],": ",a[,5]);
	a <- a[,-c(3,5,7)];
	assign(paste0(modelRef,"_","sharedDev2models.",model1,"_",model2),a,.GlobalEnv);
	sharedList <- list(stats=data.frame(ltotal.Shared=ltotal.Shared,lshared.Up_Up=lshared.Up_Up,lshared.1Up_2Down=lshared.1Up_2Down,lshared.2Up_1Down=lshared.2Up_1Down,lshared.Down_Down=lshared.Down_Down),sharedDEGs=a,total.Shared=total.Shared,shared.Up_Up=shared.Up_Up,shared.1Up_2Down=shared.1Up_2Down,shared.2Up_1Down=shared.2Up_1Down,shared.Down_Down=shared.Down_Down);
	assign(paste0(modelRef,"_","sharedList.",model1,model2),sharedList,.GlobalEnv);
	openxlsx::write.xlsx(sharedList, file = paste0(wd,f_results, f_venn,f_category,f_vdata, "dataVenn_sharedList.xlsx"));

};

vennD2_homologs.function(model_colors,"RNASeq_Beta_c0","ATACSeq_Beta_c0","DEGs","DAGs",degDagList,degDagRegList,c("Bcells c0 DEGs","Bcells c0 DAGs"),"Beta_c0")
vennD2_homologs.function(model_colors,"RNASeq_Beta_c6","ATACSeq_Beta_c6","DEGs","DAGs",degDagList,degDagRegList,c("Bcells c6 DEGs","Bcells c6 DAGs"),"Beta_c6")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c3","ATACSeq_Alpha_c3","DEGs","DAGs",degDagList,degDagRegList,c("Acells c3 DEGs","Acells c3 DAGs"),"Alpha_c3")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c5","ATACSeq_Alpha_c5","DEGs","DAGs",degDagList,degDagRegList,c("Acells c5 DEGs","Acells c5 DAGs"),"Alpha_c5")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c7","ATACSeq_Alpha_c7","DEGs","DAGs",degDagList,degDagRegList,c("Acells c7 DEGs","Acells c7 DAGs"),"Alpha_c7")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c8","ATACSeq_Alpha_c8","DEGs","DAGs",degDagList,degDagRegList,c("Acells c8 DEGs","Acells c8 DAGs"),"Alpha_c8")

vennD2_homologs.function(model_colors,"RNASeq_Beta_c0","ATACSeq_Beta_c0","DEGs","DAGs",degDagList2,degDagRegList2,c("Bcells c0 DEGs","Bcells c0 DAGs"),"Beta_c0")
vennD2_homologs.function(model_colors,"RNASeq_Beta_c6","ATACSeq_Beta_c6","DEGs","DAGs",degDagList2,degDagRegList2,c("Bcells c6 DEGs","Bcells c6 DAGs"),"Beta_c6")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c3","ATACSeq_Alpha_c3","DEGs","DAGs",degDagList2,degDagRegList2,c("Acells c3 DEGs","Acells c3 DAGs"),"Alpha_c3")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c5","ATACSeq_Alpha_c5","DEGs","DAGs",degDagList2,degDagRegList2,c("Acells c5 DEGs","Acells c5 DAGs"),"Alpha_c5")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c7","ATACSeq_Alpha_c7","DEGs","DAGs",degDagList2,degDagRegList2,c("Acells c7 DEGs","Acells c7 DAGs"),"Alpha_c7")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c8","ATACSeq_Alpha_c8","DEGs","DAGs",degDagList2,degDagRegList2,c("Acells c8 DEGs","Acells c8 DAGs"),"Alpha_c8")


vennD2_homologs.function(model_colors,"RNASeq_Beta_c0","ATACSeq_Beta_c0","DEGs","DAGs",degDagList3,degDagRegList3,c("Bcells c0 DEGs","Bcells c0 DAGs"),"Beta_c0")
vennD2_homologs.function(model_colors,"RNASeq_Beta_c6","ATACSeq_Beta_c6","DEGs","DAGs",degDagList3,degDagRegList3,c("Bcells c6 DEGs","Bcells c6 DAGs"),"Beta_c6")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c3","ATACSeq_Alpha_c3","DEGs","DAGs",degDagList3,degDagRegList3,c("Acells c3 DEGs","Acells c3 DAGs"),"Alpha_c3")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c5","ATACSeq_Alpha_c5","DEGs","DAGs",degDagList3,degDagRegList3,c("Acells c5 DEGs","Acells c5 DAGs"),"Alpha_c5")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c7","ATACSeq_Alpha_c7","DEGs","DAGs",degDagList3,degDagRegList3,c("Acells c7 DEGs","Acells c7 DAGs"),"Alpha_c7")
vennD2_homologs.function(model_colors,"RNASeq_Alpha_c8","ATACSeq_Alpha_c8","DEGs","DAGs",degDagList3,degDagRegList3,c("Acells c8 DEGs","Acells c8 DAGs"),"Alpha_c8")

###########################################################################################
###########################################################################################
##  End script for multi-ome embedded analysis: setting up the pipeline using Seurat     ##
###########################################################################################
## Finalised on 07/27/22 by Mar Muniz Moreno.
## Any Qs, or feedback please contact me.
############################################################################################
