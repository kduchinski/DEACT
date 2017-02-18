## Introduction
Analyzing transcriptome changes between healthy and abnormal cells is key in understanding disease. Researchers who want to analyze parallel or independent RNA-Seq experiments are interested in **subsets of transcripts that stand out comparatively across experiments**, but identifying these can be problematic for those **not computationally trained**. We developed a user-friendly web application designed to **visualize**, **compare**, and **subset** differentially expressed transcripts in two complementary experiments. The **D**ataset **E**xploration **A**nd **C**uration **T**ool (**DEACT**) enables researchers to compare two RNA-Seq studies using a simple interface for **managing**, **viewing**, and **plotting** data. With DEACT, researchers with little or no programming experience can:
- analyze any two complementary studies (biological replicates or contrasting experimental conditions)
- instantly identify and select any set of genes with a level of engagement that neither scripts nor traditional plots offer
- quickly answer preliminary questions about new RNA-seq data to prompt downstream analyses
- and retain a flow of discussion in a collaborative setting.

![alt text](http://i.imgur.com/RusSLEn.png)

## A: File upload
-	Maximum file size: 5MB
-	Tab-delimited text file (.txt/.csv)
-	Columns:
  *	Gene IDs
  *	Gene symbols
  *	Change in expression for each condition
  *	Log2 fold change in expression for each condition
  *	Probability values of each condition
-	See Section E for a sample R script for file creation from Cufflinks data.
- An example dataset is also available. These data represent transcripts from a knockdown and an overexpression of the FLI1 transcription factor.

## B: Plot options
- This section is not available until settings have been saved.
-	The user may choose to graph points by fpkm difference or log2 fold change.
-	Gene expression data are commonly graphed by log¬2 fold change for the sake of normalization. However, if the control value is 0 fpkm, then the fold change constitutes a division by 0. These data points will have a fold change value of positive or negative infinity and cannot be plotted. If log2¬ fold change is selected, then the user can choose whether or not to include points with infinite fold changes.
-	By default, only genes significantly affected by both conditions are displayed on the graph. Choose the "Full dataset" option to graph every point.
- If options are changed after a selection has been made, the current selection will be erased.

## C: Regulation table
- This section is not available until settings have been saved.
-	This table shows how many genes are significantly affected by each condition.
-	The last four rows are divisions of the genes affected by both conditions. In order, these rows represent quadrants I, III, IV, and II of the scatterplot respectively.
-	The buttons to the right of these will add all genes in a category to the subset.

## D: Activity panel
1. Manage data
  - When a dataset is uploaded or the example dataset has been selected, the columns of the dataset will appear in drop-down menus.
  - First, use the drop-down menus to specify the parameters as labeled. Users may also name the conditions, which default to Condition 1 and Condition 2, respectively. Condition names may be edited at any time by returning to this panel.
  - When finished, saving settings unlocks all other functionality.
2. Plot data

  ![alt text](http://imgur.com/xxgxjE5.png)

  -	Plots Condition 1 vs. control on the x-axis and Condition 2 vs. control on the y-axis in the units chosen (Section B)
  -	Click and drag to select genes, or click one of the icons above the plot to zoom, pan, or export a still version of the plot.
  -	Hover over a point to show its coordinates (expression changes) and gene symbol. 
  -	If plotting the full dataset, the hover tag will also show which conditions significantly affected a gene. 
3. View data
  ![alt text](http://imgur.com/H4rvKD4.png)
  - A table will appear here once data have been selected from the regulation table (Section C) or from the scatterplot (Section D2).
  - Users may preview the curated subset prior to downloading. The table offers sorting and searching functionality.
4. Download data
  ![alt text](http://imgur.com/tou96kA.png)
  - The user-curated subset may be downloaded as a tab-delimited text file. Users must first choose the columns of the dataset to download, as some downstream analytic tools require specific formats. Multiple versions of the file may be downloaded if desired.

## 8. Input file script from Cufflinks data
```R
cuff <- readCufflinks(dir="/jobs/tmp/Fli1//cuffdiff",rebuild=F)
mySigGeneIdsLoss<-getSig(cuff,'X5324','shc1',alpha=0.05,level="genes")
mySigGeneIdsGain<-getSig(cuff,'AdFli1','AdGFP',alpha=0.05,level="genes")
mySigGeneIdsInter = intersect(mySigGeneIdsGain,mySigGeneIdsLoss)
mySigGeneIds = union(mySigGeneIdsGain, mySigGeneIdsLoss)

#Get expression values in fpkm of all genes
fpkm.avg = fpkm(genes(cuff))

#Calculate change in expression of significant genes for each condition
fpkmdiff.loss = c()
fpkmdiff.gain = c()
for (i in 1:length(mySigGeneIds)){
  id = mySigGeneIds[i]
  fpkm.loss = fpkm.avg[fpkm.avg$gene_id==id & fpkm.avg$sample_name == 'X5324',]$fpkm
  fpkm.lctrl = fpkm.avg[fpkm.avg$gene_id==id & fpkm.avg$sample_name == 'shc1',]$fpkm
  fpkm.gain = fpkm.avg[fpkm.avg$gene_id==id & fpkm.avg$sample_name == 'AdFli1',]$fpkm
  fpkm.gctrl = fpkm.avg[fpkm.avg$gene_id==id & fpkm.avg$sample_name == 'AdGFP',]$fpkm
  fpkmdiff.gain[length(fpkmdiff.gain)+1] = fpkm.gain - fpkm.gctrl
  fpkmdiff.loss[length(fpkmdiff.loss)+1] = fpkm.loss - fpkm.lctrl
  if (fpkm.loss > fpkm.lctrl & fpkm.gain > fpkm.gctrl) {
    posLossPosGain[length(posLossPosGain)+1] = id
  } else if (fpkm.loss < fpkm.lctrl & fpkm.gain > fpkm.gctrl) {
    negLossPosGain[length(negLossPosGain)+1] = id
  } else if (fpkm.loss > fpkm.lctrl & fpkm.gain < fpkm.gctrl){
    posLossNegGain[length(posLossNegGain)+1] = id
  } else {
    negLossNegGain[length(negLossNegGain)+1] = id
  }
}
#Sort by id
df.temp = data.frame("mySigGeneIds" = mySigGeneIds,"fpkmdiff.loss" = fpkmdiff.loss, "fpkmdiff.gain" = fpkmdiff.gain)
df.temp = df.temp[1:length(mySigGeneIds),]
df.temp = df.temp[with(df.temp, order(mySigGeneIds, fpkmdiff.loss)),]
fpkmdiff.loss = df.temp$fpkmdiff.loss
fpkmdiff.gain = df.temp$fpkmdiff.gain
mySigGeneIds = df.temp$mySigGeneIds

#Retrieve annotation information for significant genes
annotation.table = annotation(genes(cuff))
annotation.table = annotation.table[which(annotation.table$gene_id %in% mySigGeneIds),]

#Retrieve CuffDiff data for significant genes
gene.ids<-featureNames(genes(cuff))
myGene<-getFeatures(cuff,gene.ids,level='genes') 
diffTable = diffData(myGene)
diffTable = diffTable[which(diffTable$gene_id %in% mySigGeneIds),]

#Store probability values and fold change values
q.loss = c()
q.gain = c()
foldChange.loss = c()
foldChange.gain = c()
for (i in 1:(nrow(diffTable))){
  if (diffTable[i,2]=='X5324' & diffTable[i,3]=='shc1'){
    q.loss[length(q.loss)+1] = diffTable[i,10]
    foldChange.loss[length(foldChange.loss)+1] = diffTable[i,7]
  } else if (diffTable[i,2]=='AdFli1' & diffTable[i,3] == 'AdGFP'){
    q.gain[length(q.gain)+1] = diffTable[i,10]
    foldChange.gain[length(foldChange.gain)+1] = diffTable[i,7]
  }
}

#Create a table for use in DEACT
fpkm.df = data.frame(“gene_short_name” = annotation.table$gene_short_name, “gene_id” = annotation.table$gene_id, “fpkmdiff.loss” = fpkmdiff.loss, “fpkmdiff.gain” = fpkmdiff.gain, “foldChange.loss” = foldChange.loss, “foldChange.gain” = foldChange.gain, “qValueLoss” = q.loss, “qValueGain” = q.gain)
write.table(fpkm.df, file = "~/fpkm", sep = "\t", row.names = F, quote = F)
```
##Acknowledgements and Contributers
Katherine Duchinski<sup>1</sup> and [Margaret Antonio](https://github.com/antmarge)<sup>2</sup>

Dr. Paul Anderson<sup>1</sup>

Drs. Dennis Watson<sup>3</sup>, Patricia Watson<sup>3</sup>, and Robert Wilson<sup>3</sup> provided primary transcriptomic data which was used to develop this application.

We would like to thank the College of Charleston for hosting the NSF Omics REU which is supervised by the National Science Foundation DBI Award 1359301. We also acknowledge support from the Genomics Shared Resource, Hollings Cancer Center, Medical University of South Carolina. This shared resource is supported in part by the Hollings Cancer Center, Medical University of South Carolina Support Grant (P30 CA 138313).

Page maintained by Katherine Duchinski.

<sup>1</sup>College of Charleston
<sup>2</sup>Boston College
<sup>3</sup>Medical University of South Carolina
