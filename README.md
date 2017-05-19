---
title: "IntramiRExploreR_Vignettes_ver05"
author: "Surajit Bhattacharyal and Daniel N. Cox"
date: "Friday, December 09, 2016"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{IntramiRExploreR}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

##Introduction

Micro RNAs (**miRNAs**) are a group of small non-coding RNAs (21-25 nucleotide long),which have been associated with post-transcriptional Gene Silencing, since its discovery in ***C.elegans*** as a regulator of larval  development (R. C. Lee, Feinbaum, & Ambros, 1993; Wightman, Ha, & Ruvkun, 1993). miRNAs play a major roles in various developmental and disease conditions, however one significant challenge in characterizing the mechanism via which miRNAs exert their post-transcriptional effect is the identification of biologically relevant target mRNAs, given that miRNAs exhibit a one-to-many relationship with their putative target mRNAs. Most microRNA Target prediction tools look at biophysical properties like seed sequence matching and Gibbs Free Energy, but that does not take into consideration the expression of the miRNA of the tissue of interest. 

This is an important parameter to take into account, as the expression of the miRNA, along with the physical properties discussed earlier, would determine effectively whether a particular miRNA plays a functional role in the process of interest. Although there have been tools available which predict targets for a given microRNA based on statistical correlation between miRNA expression and mRNA expression values. But this approach has 2 caveats:
a) The miRNA expression values both microarray and RNA-seq expression values are limited compared to those mRNAs.
b) Other than human and mouse, the number of miRNA expression data is quite limited for other model organisms like D.melanogaster, and predicting miRNA functionality across the whole genome is quite difficult

One method to bypass this impediment and use expression profiles to identify miRNA targets, in a model organism like Drosophila, is to focus on intragenic miRNAs, which are located within host protein coding genes. Intragenic miRNAs constitute approximately 60% of all miRNAs in ***Drosophila melanogaster*** (Fruit Flies), making these miRNAs an important component in post-transcriptional regulation of gene expression. Reports have confirmed that the expression of intragenic miRNAs is highly correlated with the expression of the host gene mRNA (Baskerville & Bartel 2005; Karali et al. 2007; Kim & Kim 2007). Based upon this correlation, it is possible to use the host gene expression values as a proxy for the expression of the intragenic miRNA (Tsang et al. 2007). Target prediction of intragenic miRNAs using the host gene expression has been successfully implemented in Humans, with the HOCTAR algorithm (Gennarino et al. 2008; Gennarino et al. 2011). Given that there are much larger available datasets for mRNA expression profiles, than miRNA expression profiles, by using the host gene as a proxy for the intragenic miRNA, one can significantly extend bioinformatic analyses and statistical power in predicting miRNA:mRNA target interactions, that are rooted not only in target prediction algorithms, but also in biologically relevant, and inversely correlated, patterns of expression between a miRNA and potential target mRNAs.

The IntramiRExploreR tool using 2 distinct correlation methods Distance and Pearson correlation finds targets for miRNA in Fruit flies using the availbale Affymetric microarray data available in the Gene Expression omnibus database. Other than the targets the tool also integrates Gene Ontology functionalities using FGNet(Bioconductor), Data from NCBI and visualisation tool using igraph.


###Installing the package

**IntramiEExploreR** is currently available from the github repository. Installation method would be as the following:

```{r eval=FALSE}
library("devtools")
devtools::install_github("sbhattacharya3/IntramiRExploreR")
```
```{r eval=TRUE}
library("IntramiRExploreR")
```
IntramiRExploreR has dependency on  R version (>= 3.1.2). To use the DAVID functionality for Gene Ontology functional classification (called from **GetGOS_ALL** function), user has to install the **RDAVIDWebService** package using the link below:
http://stackoverflow.com/questions/31480579/r-david-webservice-sudden-transport-error-301-error-moved-permanently.

###Target Prediction using expression data

For building up the intragenic miRNA target data base, we have used Affymetrix platform 1 & 2 microarray datasets for ***D.melanogaster***, from GEO database(Barrett et al., 2013). For the significance of the statistical analysis, experiments with greater than or equal to 5 assays were considered. The experiments were normalized using the Robust Multichip Average (RMA) (Bolstad, Irizarry, Astrand, & Speed, 2003) from the affy package (Gautier, Cope, Bolstad, & Irizarry, 2004)(Gautier et al., 2004) from the Bioconductor suite (Gentleman et al., 2004) in R. The statistical functions are then used to find the correlation between the host genes and each of the other genes in an experiment. The correlation methods used are Pearson Correlation (Lee Rodgers & Nicewander, 1988; Pearson, 1895) and distance correlation(Szekely & Rizzo, 2009). 

After the correlation analysis has been performed, a false discovery rate calculation, Benjamini Hochberg (BH) False Dicovery Rate (FDR) Calculation (Benjamini & Hochberg, 1995) is done on the p values obtained for each miRNA-mRNA pair for a particular experiment, using the p.adjust function in R. To identify statistically significant, anti-correlated mRNA targets (p<0) for a particular miRNA, all mRNAs with a q-value (FDR threshold) of less than 0.01 are selected across all experiments.  From these analyses, the top 25% most frequently occurring mRNAs are then compared with the targets predicted for a given miRNA in a variety of target databases (TargetScan, PITA, and Miranda).  A target gene which is found in the output list of both the statistical tests and also found in the target database can be called as a putative target for a given miRNA. To get the most important putative targets a scoring system has also been designed. The scoring system is a summation of  3 parameters:

1)  Probability of sequence conservation of both the targets and the miRNA, across the different species considered in TargetScan and Miranda databases.
2)  Number of complementary sites in a target for a given miRNA obtained from Pita Target Database.
3)  Probability of occurrence of the target:miRNA pair across the different experiments
 
These Statistically predicted targets for a given miRNA of interest can be obtained using **miRTargets_Stat** function, but can be visualized by the user using the **Visualisation** function.

These Statistically predicted targets for a given miRNA of interest can be obtained using **miRTargets_Stat** function.

```{r eval=TRUE}
miR="dme-miR-12"
a<-miRTargets_Stat(miR,method=c("Both"),Platform=c("Affy1"),Text=FALSE)
a[1:4,1:5]
```
The input to the function are single or multiple miRNAs, the Statitical method which predicts the target, and the platform. The method chosen here is "Both" which is an union of both the **Pearson** and the **Distance** correlation method. The platform is chosen as **Affy1** (Affymetrix platform1).  The  output from the function is targets that are statistically significant, the score associated to each target, the GEO accession IDS where the miRNA and the Targets are correlated and the function of the target genes from the flybase.

Similarly, **genes_Stat** is used to obtain statistically relevant miRNAs that target a gene of interest.

```{r eval=TRUE}
gene ="Ank2"
a<-genes_Stat(gene,geneIDType="GeneSymbol", method=c("Both"),Platform=c("Affy1"))
a[1:4,1:5]
```
**genes_Stat** has similar output format as **miRTargets_Stat**, the only difference is that it outs the miRNA function from flybase, instead of the genes.

**Visualisation** function has three output formats:
a)**text**: Output miRNA targets result obtained from **miRTargets_Stat**, in text format.
b)**Cytoscqape**: Output in the format of cytoscape input files.
c)**igraphs**: Output  miRNA:Target gene results in the form of network.
d)If no output format is chosen, a datframe containing the result returned to the user.

```{r eval=TRUE}
miR=c("dme-miR-12","dme-miR-283")
a<-Visualisation(miR,mRNA_type=c("GeneSymbol"),method=c("Both"),platform=c("Affy1"),
                visualisation=c("console"),thresh=10)
a[1:10,1:5]
```
The input to the function are single or multiple miRNAs, the Statitical method which predicts the target, and the platform. The method chosen here is "Both" which is an union of both the **Pearson** and the **Distance** correlation method. The platform is chosen as **Affy1** (Affymetrix platform1).  The  output from the function is targets that are statistically significant, the score associated to each target, the GEO accession IDS where the miRNA and the Targets are correlated and the function of the target genes from the flybase.

The output can be visualised using **igraph**.


Similarly, **Genes_Visualisation** is used to obtain statistically relevant miRNAs that target a gene of interest, as an output from **genes_Stat** function.

```{r eval=TRUE}
mRNA="Syb"
a<-Gene_Visualisation(mRNA,mRNA_type=c("GeneSymbol"),method=c("Pearson"),
            platform=c("Affy1"),visualisation= "console")
a[1:10,1:5]
```

The output can be visualised using **igraph**, similar to visualisation function.

####Gene Ontology

**GetGOS_ALL** function outputs functional network clusters, using FGNet. topGO and DAVID are the 2 available GO methods. 
```{r, eval=FALSE}
miR="dme-miR-12"
a<-Visualisation(miR,mRNA_type=c("GeneSymbol"),method=c("Both"),platform=c("Affy1"),thresh=100,
            visualisation="console")
genes<-a$Target_GeneSymbol
GetGOS_ALL(genes,GO=c("topGO"),term=c("GO_BP"),path="C://",filename="test")
```

####References

1. Baskerville, S., & Bartel, D. P. (2005). **Microarray profiling of microRNAs reveals frequent coexpression with neighboring miRNAs and host genes.** RNA (New York, N.Y.), 11(3), 241–7. http://doi.org/10.1261/rna.7240905
2. Gennarino, V. A., Sardiello, M., Avellino, R., Meola, N., Maselli, V., Anand, S., … Banfi, S. (2008). **MicroRNA target prediction by expression analysis of host genes.** Genome Research, 19(3), 481–490. http://doi.org/10.1101/gr.084129.108
3. Gennarino, V. A., Sardiello, M., Mutarelli, M., Dharmalingam, G., Maselli, V., Lago, G., & Banfi, S. (2011). **HOCTAR database: A unique resource for microRNA target prediction.** Gene, 480(1–2), 51–58. http://doi.org/10.1016/j.gene.2011.03.005
4. Karali, M., Peluso, I., Marigo, V., & Banfi, S. (2007). **Identification and characterization of microRNAs expressed in the mouse eye.** Investigative Ophthalmology & Visual Science, 48(2), 509–15. http://doi.org/10.1167/iovs.06-0866
5. Kim, Y.-K., & Kim, V. N. (2007). **Processing of intronic microRNAs.** The EMBO Journal, 26(3), 775–83. http://doi.org/10.1038/sj.emboj.7601512
6. Lee, R. C., Feinbaum, R. L., & Ambros, V. (1993). **The C. elegans heterochronic gene lin-4 encodes small RNAs with antisense complementarity to lin-14.** Cell, 75(5), 843–54. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/8252621
7. Tsang, J., Zhu, J., & van Oudenaarden, A. (2007). **MicroRNA-mediated feedback and feedforward loops are recurrent network motifs in mammals.** Molecular Cell, 26(5), 753–67. http://doi.org/10.1016/j.molcel.2007.05.018
8. Wightman, B., Ha, I., & Ruvkun, G. (1993). **Posttranscriptional regulation of the heterochronic gene lin-14 by lin-4 mediates temporal pattern formation in C. elegans.** Cell, 75(5), 855–62. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/8252622


