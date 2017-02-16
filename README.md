## Introduction
Analyzing transcriptome changes between healthy and abnormal cells is key in understanding disease. Researchers who want to analyze parallel or independent RNA-Seq experiments are interested in **subsets of transcripts that stand out comparatively across experiments**, but identifying these can be problematic for those **not computationally trained**. We developed a user-friendly web application designed to **visualize**, **compare**, and **subset** differentially expressed transcripts in two complementary experiments. The **D**ataset **E**xploration **A**nd **C**uration **T**ool (**DEACT**) enables researchers to compare two RNA-Seq studies using a simple interface for **managing**, **viewing**, and **plotting** data. With DEACT, researchers with little or no programming experience can:
- analyze any two complementary studies (biological replicates or contrasting experimental conditions)
- instantly identify and select any set of genes with a level of engagement that neither scripts nor traditional plots offer
- quickly answer preliminary questions about new RNA-seq data to prompt downstream analyses
- and retain a flow of discussion in a collaborative setting.

## 1. File upload
-	Maximum file size: 5MB
-	Tab-delimited text file (.txt/.csv)
-	Columns:
  *	Gene IDs
  *	Gene symbols
  *	Change in expression for each condition
  *	Log2 fold change in expression for each condition
  *	Probability values of each condition
-	Code for file creation in Section 8 
## 2. Column selection
-	Use the drop-down arrows to select the column corresponding to the specified information or click on the box and type in the column name
-	The table in Section 3 will be updated automatically as information is filled in
## 3. Regulation table
-	This table shows how many genes are significantly affected by each condition
-	The last four rows are divisions of the genes affected by both conditions. In order, these rows represent quadrants I, III, IV, and II of the scatterplot respectively
-	The buttons to the right of these will add all genes in a category to the subset

## 4. Plot units toggle
-	Choose to graph points by fpkm difference or log2 fold change before selecting from the plot. Switching between them will erase the current selection.
-	Gene expression data are commonly graphed by log¬2 fold change. However, if the control value is 0 fpkm, then the fold change constitutes a division by 0. These data points will have a fold change value of positive or negative infinity and cannot be plotted. If log2¬ fold change is selected, then the user can choose whether or not to include points with infinite fold changes.

## 5. Scatterplot
-	When the parameters have been set through the drop-down menus in Section 2, set the plot toggle to “Make plot”
-	A scatterplot will appear in place of the column selection menus
-	Plots condition 1 vs. control on the x-axis and condition 2 vs. control on the y-axis in the units chosen (Section 3)
-	Regulation information for individual points may be derived from the plot. Points to the right of the y-axis were up-regulated in condition 1 and points above the x-axis 

## 6. Page navigation and “Full Dataset” tab
-	On the first page, only genes significantly affected by both conditions are displayed on the graph. 
-	Click on the “Full Dataset” tab to see the graph of every point. Points from this graph can be added to the subset as well. A table of these points will be shown below the scatterplot on this page. 
-	The button in the sidebar on the “Full Dataset” page will add these points to the table on the first page.

## 7. Download button
-	Download the selected data as a tab-delimited text file containing the log2 fold change, probability value, and gene symbol columns.

## 8. Input file script
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
Katherine Duchinski<sup>1</sup> and [Margaret Antonio](https://github.com/mmlantonio)<sup>2</sup>
Dr. Paul Anderson<sup>1</sup>

Drs. Dennis Watson<sup>3</sup>, Patricia Watson<sup>3</sup>, and Robert Wilson<sup>3</sup> provided primary transcriptomic data which was used to develop this application.

We would like to thank the College of Charleston for hosting the NSF Omics REU which is supervised by the National Science Foundation DBI Award 1359301. We also acknowledge support from the Genomics Shared Resource, Hollings Cancer Center, Medical University of South Carolina. This shared resource is supported in part by the Hollings Cancer Center, Medical University of South Carolina Support Grant (P30 CA 138313).

Page maintained by Katherine Duchinski.

<sup>1</sup>College of Charleston
<sup>2</sup>Boston College
<sup>3</sup>Medical University of South Carolina
