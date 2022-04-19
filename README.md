# Microarray_data_analysis_2022

To be updated 06/04/2022
---
title: "Full Project Workflow"
author: "Greta Paulauskaite"
output: word_document: default
mainfont: Times New Roman
---

To begin analysis, check the version of R and update if needed. In this workflow, I used the RStudio integrated development environment and the R version can be updating by checking the 'Help' tab. Once the version is up to date, install packages if needed.

```
Instalation can be done either by install.packages(__) function or by BiocManager::install(__). 
```

Load libraries that will be used in this analysis. 

```{r message = FALSE, warning=FALSE}
library(GEOquery)
library(Biobase)
library(AnnotationDbi)
library(BiocGenerics)
library(dplyr)
library(hgu133a.db)
library(readr)
library(limma)
library(sparklyr)
library(knitr)
library(reticulate)
library(pvca)
library(stringr)
library(topGO)
library(ReactomePA)
library(ggnewscale)
library(clusterProfiler)
library(genefilter)
```


#### connection to Apache Spark

To address the memory problems with handling big datasets on local computer, and to optimise the time it takes for R to solve complicated commands, Apache Spark extension is used. Connection with Spark will be done later in the pipeline. To use the application, you need to install Java. Check Windows configuration, in this case 64-bit, and install the application from java.com. 

```
install.packages("spkarlyr", repos = "http://cran.us.r-project.org")
library(sparklyr)
spark_available_versions() ## check for the latest versions available at the time.
spark_install("3.2")
Sys.setenv(JAVA_HOME="")
system("java -version") ## check the java version installed.

```

#### Connection to Python

Installation of package 'reticulate' and loading the library is necessary. To connect the Python packages scikit-learn, pandas, numpy and matplotlib, first the download of Anaconda to local computer needs to be done. During installation, make sure to tick the 'Add Anaconda3 to my PATH environment variable' as it is critical to obtain the environment in RStudio. After the installation is complete, in Terminal type the following:
```
conda create -n py4.10 python=4.10 scikit-learn pandas numpy matplotlib
```
Always check the current and newest version of conda. Connect the environment and check if the connection is active with `conda list env` in the Terminal. 
Back in Console, after loading the package, check the connectivity.
```
reticulate::conda_list()
```
Next, set your environment. The connection of python syntax will be indicated by '>>>' instead of '>' in the console. 
```
use_condaenv("py3.8", required = TRUE)
py_config() ## for double checking the connection
```
#################################################

#### 1. Data Exploration and Filtering

Get GSE dataset via GEO query. Setting to attaching GPL soft file as an annotation, however additional annotation steps will done at a later stage.

``` {r results='hide', warning= FALSE, message=FALSE}
Sys.setenv("VROOM_CONNECTION_SIZE" = 27983872 * 2) ## to solve VROOM_CONNECTION_SIZE error.
gset <- getGEO("GSE68468", GSEMatrix = TRUE, getGPL = TRUE)
length(gset)
gse <- gset[[1]]
```
Explore the dataset.Expression set is a complex matrix, and the Bioconductor projects have built-in functions to extract the data that was submitted to the GEO database.

```{r results='hide', warning=FALSE, message=FALSE}
pData(gse) ## phenotypic data.
featureNames(gse) ## check the gene ID's attached to the experiment.
varLabels(gse) ## retrieve variable labels.
anyMissing(gse) ## checks for missing values in an object.
sampleNames(gse) ## retrieve the sample ID's.
featureData(gse) ## check the object class.
fData(gse) ## access the feature data.
experimentData(gse) ## access the details of the experiment, including abstract from the original published paper.
validObject(gse) ## quality check the data. All assayData components have the same number of features and samples and the number, and names of phenoData rows match the number and names of assayData columns. 
exprs(gse) ## access the expression readings. Expression set produces G x N, where G is the number of genes on a chip and N is the number of tissues analyzed.
annotation(gse) ## access the platform that was used to carry out experiment, in this case GPL96.
```
```{r results='hide', warning=FALSE, message=TRUE}
ncol(gse) 
nrow(gse)
```


Log transformation will be needed throughout every step, thus this will occur first.

```{r results='hide', warning=FALSE, message=FALSE}
exprs(gse) <- log2(exprs(gse))
```

To increase reproducibility and quality control data if issues arise, each change made to the dataset will be assigned to a new variable. Separation of control probes will take place.

Remove Affymetrix control probes. The number of control probes will be available by looking at platform specifications. GPL96 has 68 probes all starting with AFFX suffix. Whether probes will stay or not in analysis will be based in the type analysis will occur. For now, control probes will stay in the main 'gse' expression set, but a value of ID numbers will be separated for later use.

```{r results='hide', warning=FALSE, message=FALSE}
controlProbes <- grep("AFFX", featureNames(gse))
controlProbes ## quality control that all 68 have been selected.

```

#####################################################
#### Annotation unfinnished

Mapping gene ID's to the affymetrix probe ID's is a step that will aid to recognize gene by their consensus names.
The package quality control information can be accessed by calling capture output function. The following result will illuminate how many of the keys are mapped in different hgu133a.db package objects. Since in this study we are looking at the gene expression and are aiming to map as many genes as possible, the hgu113aACCNUM has the most IDs matched.
Looking at the documentation of hgu133a.db package, the ACCNUM is a map of manufacturer identifiers to accession numbers. The identifiers can be also downloaded from NCBI GEO databses under the platform GPL96. 

```{r include=TRUE, message=TRUE, comment="", prompt=TRUE, echo=FALSE}
qcdata = capture.output(hgu133a())
kable(head(qcdata, 20), format = "markdown", caption = "Summary of the annotation package identifiers", align ="lr")

```


```{r results='hide', warning=FALSE, message=FALSE}
columns(hgu133a.db)
keytypes(hgu133a.db)
keys <- featureNames(gse)
anno <- AnnotationDbi::select(hgu133a.db, 
                       keys=keys, 
                        columns=c("SYMBOL"), 
                         keytype = "PROBEID")
annoGSE <- subset(anno, !is.na(SYMBOL))
anyMissing(annoGSE$SYMBOL)
anno_grouped <- group_by(annoGSE, PROBEID)
anno_summarized <- dplyr::summarize(anno_grouped, no_of_mathces = n_distinct(SYMBOL))
no_of_matches <- n_distinct(anno_grouped$SYMBOL)
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
```

To understand better, let's run some test to see how well the annotation worked. 

```{r results='hide', warning=FALSE, message=TRUE}
multi <- toggleProbes(hgu133aSYMBOL, "all")
dim(multi)
multiOnly <- toggleProbes(multi, "multiple")
dim(multiOnly)
singleOnly <- toggleProbes(multiOnly, "single")
dim(singleOnly)
hasMultiProbes(multiOnly)
hasSingleProbes(multiOnly)
hasMultiProbes(singleOnly)
hasSingleProbes(singleOnly)
multiTable <- toTable(multiOnly)
singleTable <- toTable(singleOnly)
count.mappedkeys(singleOnly)
count.mappedkeys(multiOnly)
table(multiTable$symbol)
multiCount <- as.data.frame(table(multiTable$symbol)) %>% filter(Freq == 1)
```
#####################################################


#### Evaluating the data 

Before we make any changes to the data, we must inspect if the data is up to standards to analyse in the first place. There are a couple techniques to evaluate if any transformations are needed to be done to the data. 
First step is to evaluate the whole dateset as a whole.
One of the techniques to do that is to get a detail report using Array Quality Metrics.

```{r results='hide', warning=FALSE, message=FALSE}
arrayQualityMetrics(expressionset = gse,
                 outdie = tempdir(),
                 force = TRUE, 
                 do.logtransform = FALSE,
                 intgroup = c("characteristics_ch1.4"))
```

Since the experiment has 12 different groups attached totalling in 390 samples the overall effect of the Array Quality Metrics is less clear. To better visualise and inspect if data needs smoothing filter data by histology of the samples.


Inspect the most suitable variables for filtering. The aim is to make sure that filtering does not miss out a single variable. In this experiment, the best filtering method can be done by sorting out the "characteristics_ch1.4" column, which specifies the sample histology. Columns like disease_state:ch1, characteristics_ch1 and others fails to summarize groups exclusively by the type of disease it comes from. 
To quality control the filtering, check where all 390 samples come from and how many representative samples of a group we have.

Paste the phenotypic data into a new variable.

```{r results='hide', warning=FALSE, message=TRUE}
pheno <- pData(gse)
head(pheno) ## check the columns of interest again.
if ((rownames(pheno)) == colnames(exprs(gse)))  
{
  print("TRUE")
} 
```

```{r include=TRUE, message=TRUE, comment="", prompt=TRUE, echo=FALSE}
kable(table(pheno$characteristics_ch1.4), format = "markdown", caption = "Summary of samples in the experiment", col.names=c("Histology", "Sample Count"), align="lr")
```

Based on the table, we select the groups of interest for analysis. The more samples of a histology group, the more meaningful and powerful are the statistical tools. First, separate the each histology group into a new ExpressionSet. 

```{r results='hide', warning=FALSE, message=TRUE}
filter1 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon cancer"]

if (filter1 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon cancer"])
{
   print("TRUE")
}

coloncancer <- gse[, filter1]

table(coloncancer@phenoData$characteristics_ch1.4)

filter2 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal colon mucosa"]

if (filter2 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal colon mucosa"])
{
  print("TRUE")
}

normalcolon <- gse[, filter2]

table(normalcolon@phenoData$characteristics_ch1.4)

filter3 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon polyp"]

if (filter3 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon polyp"])
{
  print("TRUE")
}

colonpolyps <- gse[, filter3]

table(colonpolyps@phenoData$characteristics_ch1.4)

filter4 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the liver"]

if (filter4 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the liver"])
{ 
 print("TRUE")  
}

livermet <- gse[, filter4]

table(livermet@phenoData$characteristics_ch1.4)

filter5 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal liver"]

if (filter5 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal liver"])
{
  print("TRUE")
}

normalliver <- gse[, filter5]

table(normalliver@phenoData$characteristics_ch1.4)

filter6 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the lung"]

if (filter6 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the lung"])
{
  print("TRUE")
}

lungmet <- gse[, filter6]

table(lungmet@phenoData$characteristics_ch1.4)

filter7 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal lung"]

if (filter7 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal lung"])
{
  print("TRUE")
}

normallung <- gse[, filter7]

table(normallung@phenoData$characteristics_ch1.4)
```

Next, for the limma analysis it's uself to have the groups of interest in comparison in the dataset. Thus, prepare the groups for the differential expression study.

```{r results='hide', warning=FALSE, message=TRUE}

vs1 <- colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa"]

if (vs1 == colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa"])
{
  print("TRUE")
}

cancerVSnormal <- gse[, vs1]

table(cancerVSnormal@phenoData$characteristics_ch1.4)

vs2 <- colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon polyp" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa"]

if (vs2 == colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon polyp" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa"])
{
  print("TRUE")
}

polypVSnormal <- gse[, vs2]

table(polypVSnormal@phenoData$characteristics_ch1.4)

vs3 <- colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon polyp" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa" | gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer"]

if (vs3 == colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon polyp" | gse@phenoData@data$"characteristics_ch1.4"=="histology: normal colon mucosa" | gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer"])
{
  print("TRUE")
}

normalVSpolypVScancer <- gse[, vs3]

table(normalVSpolypVScancer@phenoData$characteristics_ch1.4)
 

vs4 <- colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon polyp" | gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer"]

if (vs4 == colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon polyp" | gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer"])
{
  print("TRUE")
}

cancerVSpolyp <- gse[, vs4]

table(cancerVSpolyp@phenoData$characteristics_ch1.4)

vs5 <- colnames(gse)[ gse@phenoData$"characteristics_ch1.4"=="histology: normal liver" | gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the liver"]

if (vs5 == colnames(gse)[ gse@phenoData$"characteristics_ch1.4"=="histology: normal liver" | gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the liver"])
{
  print("TRUE")
}

normalliverVSlivermet <- gse[, vs5]

table(normalliverVSlivermet@phenoData$characteristics_ch1.4)

vs6 <- colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the liver"]

if (vs6 == colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the liver"])
{
  print("TRUE")
}

cancervslivermet <- gse[, vs6]

table(cancervslivermet@phenoData$characteristics_ch1.4)

vs7 <- colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal lung" | gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the lung"]

if (vs7 == colnames(gse)[gse@phenoData$"characteristics_ch1.4"=="histology: normal lung" | gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the lung"])
{
  print("TRUE")
}

normallungVSlungmet <- gse[, vs7]

table(normallungVSlungmet@phenoData$characteristics_ch1.4)

vs8 <- colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the lung"]

if (vs8 == colnames(gse)[gse@phenoData@data$"characteristics_ch1.4"=="histology: colon cancer" | gse@phenoData$"characteristics_ch1.4"=="histology: colon carcinoma metastatic to the lung"])
{
  print("TRUE")
}

cancerVSlungmet <- gse[, vs8]

table(cancerVSlungmet@phenoData$characteristics_ch1.4)
```

Let's use limma's function of Plot Densities to quality control before applying fitting. 

```{r results='hide', warning=FALSE, message=FALSE}
plotDensities(cancerVSnormal, legend = FALSE, main = "Colon cancer versus normal colon density plot")
plotDensities(polypVSnormal, legend = FALSE, main = "Colon polyp versus normal colon density plot")
plotDensities(normalVSpolypVScancer, legend = FALSE, main = "Colon polyp, colon cancer and normal colon mucosa density plot")
plotDensities(cancerVSpolyp, legend = FALSE, main = "Colon cancer versus colon polyp density plot")
plotDensities(normalliverVSlivermet, legend = FALSE, main = "Normal liver versus liver metastasis density plot")
plotDensities(cancervslivermet, legend = FALSE, main = "Colon cancer versus liver metastasis density plot")
plotDensities(normallungVSlungmet, legend = FALSE, main = "Normal lung versus lung metastasis density plot")
plotDensities(cancerVSlungmet, legend = FALSE, main = "Colon cancer versus lung metastasis density plot")
```


```{r results='hide', warning=FALSE, message=FALSE}
arrayQualityMetrics(cancerVSnormal, outdit = tempdir(), do.logtransform = FALSE)
expressioncancerVSnormal <- exprs(cancerVSnormal)

```

For limma pipeline, we will use the datasets created that have group contrasts that we are interested in,
To visualise why the quantile normalisation will be applied to the limma pipeline, the effect can be easily seen by plotting the boxplot of the data. For large datasets, the effect is harder to see due to the number of samples present. In R, you can double check that none of the data was quantile normalized in the original study by plotting it. For report purposes, I will include the smaller dataset that will be used to visualize the effect. 

```{r results='hide', warning=FALSE, message=FALSE}
oligo::boxplot(cancerVSnormal,  main = "Colon cancer versus normal colon boxplot")
oligo::boxplot(polypVSnormal,  main = "Colon polyp versus normal colon boxplot")
oligo::boxplot(normalVSpolypVScancer,  main = "Colon polyp, colon cancer and normal colon mucosa boxplot")
oligo::boxplot(cancerVSpolyp, main = "Colon cancer versus colon polyp boxplot")
oligo::boxplot(normalliverVSlivermet, main = "Normal liver versus liver metastasis boxplot")
oligo::boxplot(cancervslivermet,  main = "Colon cancer versus liver metastasis boxplot")
oligo::boxplot(normallungVSlungmet, main = "Normal lung versus lung metastasis boxplot")
oligo::boxplot(cancerVSlungmet, main = "Colon caner versus lung metastasis boxplot")

```


The observed effect will skew downstream analysis process. Thus, we will perform quantile normalisation between arrays to even out the values.

```{r results='hide', warning=FALSE, message=FALSE}
QnormallungVSlungmet <- normallungVSlungmet
exprs(QnormallungVSlungmet) <- normalizeBetweenArrays(exprs(QnormallungVSlungmet))
oligo::boxplot(QnormallungVSlungmet, main ="Normal lung versus lung metastasis normalised")
```

```{r fig.align='left`}
par(mfrow=c(2,1))
oligo::boxplot(normallungVSlungmet, main = "Before normalisation")
oligo::boxplot(QnormallungVSlungmet, main = "After normalisation")
```

```{r results='hide', warning=FALSE, message=FALSE}

QnormalliverVSlivermet <- normalliverVSlivermet
exprs(QnormalliverVSlivermet) <- normalizeBetweenArrays(exprs(QnormalliverVSlivermet))
oligo::boxplot(QnormalliverVSlivermet)

QcancerVSlungmet <- cancerVSlungmet
exprs(QcancerVSlungmet) <- normalizeBetweenArrays(exprs(QcancerVSlungmet))
oligo::boxplot(QcancerVSlungmet)

Qcancervslivermet <- cancervslivermet
exprs(Qcancervslivermet) <- normalizeBetweenArrays(exprs(Qcancervslivermet))
oligo::boxplot(Qcancervslivermet)

QcancerVSpolyp <- cancerVSpolyp
exprs(QcancerVSpolyp) <- normalizeBetweenArrays(exprs(QcancerVSpolyp))
oligo::boxplot(QcancerVSpolyp)

QnormalVSpolypVScancer <- normalVSpolypVScancer
exprs(QnormalVSpolypVScancer) <- normalizeBetweenArrays(exprs(QnormalVSpolypVScancer))
oligo::boxplot(QnormalVSpolypVScancer)

QpolypVSnormal <- polypVSnormal
exprs(QpolypVSnormal) <- normalizeBetweenArrays(exprs(QpolypVSnormal))
oligo::boxplot(QpolypVSnormal)

QcancerVSnormal <- cancerVSnormal
exprs(QcancerVSnormal) <- normalizeBetweenArrays(exprs(QcancerVSnormal))
oligo::boxplot(QcancerVSnormal)
```

Another check how efficient is quantile normalisation is more clearly visible when re-plotting density plots. For example, easiest to visualise is the smaller datasest, so again, let's choose normal lung and secondary lung metastasis dataset.

```{r fig.align='left'}
par(mfrow=c(1,2))
plotDensities(normallungVSlungmet, legend=FALSE, "Density plot before normalisation", main = "Before normalisation")
plotDensities(QnormallungVSlungmet, legend=FALSE, "Density plot after normalisation", main = "After normalisation")
```

Checking the other datasets, and observing if the normalisation had the same effect on all the datasets.Density plots occupy a similar density pattern, thus we can safely conclude that the data points from our dataset are now coming from the same distribution.

```{r results='hide', warning=FALSE, message=FALSE}
plotDensities(QpolypVSnormal, legend = FALSE)
plotDensities(QnormalVSpolypVScancer, legend = FALSE)
plotDensities(QnormalliverVSlivermet, legend = FALSE)
plotDensities(QcancerVSpolyp, legend = FALSE)
plotDensities(QcancerVSnormal, legend = FALSE)
plotDensities(QcancerVSlungmet, legend = FALSE)
plotDensities(Qcancervslivermet, legend = FALSE)
```

Visualing plots with PCA.

```{r fig.align='left'}
EX <- Biobase::exprs(QcancerVSnormal)
PCA1 <- prcomp(t(EX), scale = FALSE)
percentVar <- round(100*PCA1$sdev^2/sum(PCA1$sdev^2),1)
sd_ratio <- sqrt(percentVar[2]/percentVar[1])
dataGG <- data.frame(PC1 = PCA1$x[,1], PC2 = PCA1$x[,2],
            Histology = Biobase::pData(QcancerVSnormal)$characteristics_ch1.4)
ggplot(dataGG, aes(PC1, PC2)) + 
          geom_point(aes(colour = Histology)) + 
          ggtitle("PCA plot of normalised data") +
          xlab(paste0("PC1, VarExp:", percentVar[1], "%")) +
          ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          coord_fixed(ratio = sd_ratio) +
          scale_color_manual(values = c("darkorange2", "dodgerblue4"))

```


Viualising batch effects. 

```{r results='hide', warning=FALSE, message=FALSE}
pct_threshold <- 0.6
colnames(pData(Qcancervslivermet))
batch.factors <-c("source_name_ch1", "characteristics_ch1.4")
pvca0bj <- pvcaBatchAssess(Qcancervslivermet, batch.factors, pct_threshold)
```

Visualing most variables genes.

```{r fig.align='left'}
sds <- apply(exprs(QcancerVSnormal), 1, sd)
sds0 <- sort(sds)
plot(1:length(sds0), sds0, main = "Distribution of variability for all genes", 
      sub="Vertical lines represent 90% and 95% percentiles",
      xlab="Gene index (from least to most variable)",
      ylab="Standard deviation")
 abline(v=length(sds)*c(0.9,0.95))
```
Filter duplicates and lowly expressed genes.

```{r}
AQcancervsnormal <- QcancerVSnormal
annotation(AQcancervsnormal) <- "hgu133a.db"
filtered <- nsFilter(AQcancervsnormal, 
                      require.entrez = TRUE, 
                      remove.dupEntrez = TRUE,
                      var.filter = TRUE,
                      var.func = IQR,
                      var.cutoff = 0.95,
                      filterByQuantile = TRUE,
                      feature.exclude = "^AFFX")
print(filtered$filter.log)
filtered_cancervsnorm <- filtered$eset
nrow(filtered_cancervsnorm)
head(filtered_cancervsnorm@featureData$"Gene Symbol")
 
```

Design the matrix.

```{r results='hide', warning=FALSE, message=FALSE}
metadata <- pData(filtered_cancervsnorm)
design <- model.matrix(~0+metadata$characteristics_ch1.4)
colnames(design) <- c("Cancer", "Normal")
fit <- lmFit(filtered_cancervsnorm, design)
contrasts <- makeContrasts(Cancer - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
class(fit2)
results <- decideTests(fit2)
summary(results)
topTab <- topTable(fit2, number = nrow(fit2), adjust = "fdr")
cancerVSnormTopTable <- topTab[, ! names(topTab) %in% c("ID", "GB_ACC", "SPOT_ID", "Species.Scientific.Name", "Annotation.Date", "Sequence.Source", "Target.Description", "Representative.Public.ID", "Refeq.Transcript.ID", "Gene.Ontology.Biological.Process", "Gene.Ontology.Cellular.Component", "Gene.Ontology.Molecular.Function")]
cancerVSnormTopTable <- cancerVSnormTopTable[, ! names(cancerVSnormTopTable) %in% c("Sequence.Type", "Gene.Title", "RefSeq.Transcript.ID")]

results1 <- decideTests(fit2, method = "separate", adjust.method = "fdr", p.value = 0.1, lfc = 1)

sum.res.rows <- apply(abs(results1), 1, sum)

res.selected <- results1[sum.res.rows!=0,]

```

```{r include=TRUE, message=TRUE, comment="", prompt=TRUE, echo=FALSE}
summary(results)
```


```{r include=TRUE, message=TRUE, comment="", prompt=TRUE, echo=FALSE}
print(summary(results1))
```

```{r results='hide', warning=FALSE, message=FALSE}
colnames(cancerVSnormTopTable)
names(cancerVSnormTopTable)[names(cancerVSnormTopTable) == "ENTREZ_GENE_ID"] <- "ENTREZID"

listOfTables <- list(ColonCancerVSNormalColon = cancerVSnormTopTable)
listofSelected <- list()

for (i in 1:length(listOfTables)) 
  {
  topTable <- listOfTables[[i]]
  whichGenes <- topTable["adj.P.Val"] < 0.15
  selectedIDs <- rownames(topTable)[whichGenes]
  EntrezIDs <- AnnotationDbi::select(hgu133a.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listofSelected[[i]] <- EntrezIDs
   names(listofSelected)[i] <- names(listOfTables)[i]
}
sapply(listofSelected, length)
listOfSelected <- listofSelected
```

Mapping

```{r results='hide', warning=FALSE, message=FALSE}

mapped_genes2GO <- mappedkeys(hgu133aGO)
mapped_genes2KEGG <- mappedkeys(hgu133aPATH)
mapped_genes <- union(mapped_genes2GO, mapped_genes2KEGG)
library(ReactomePA)

dir.create("results")
listOfData <- listOfSelected[1]
comparisonsNames <- names(listOfData)
universe <- mapped_genes

for (i in 1:length(listOfData)){
genesIn <- listOfData[[i]]
comparison <- comparisonsNames[i]
enrich.result <- enrichPathway(gene = genesIn,
                                 pvalueCutoff = 0.05,
                                qvalueCutoff = 0.9,
                                 readable = T,
                                 pAdjustMethod = "BH",
                                 organism = "human",
                                 universe = universe)
  
  cat("##################################")
  cat("\nComparison: ", comparison,"\n")
  print(head(enrich.result))
 
   if (length(rownames(enrich.result@result)) != 0) {
   write.csv(as.data.frame(enrich.result), 
              file =paste0("./results/","ReactomePA.Results.",comparison,".csv"), 
             row.names = FALSE)
   
   pdf(file=paste0("./results/","ReactomePABarplot.",comparison,".pdf"))
     print(barplot(enrich.result, showCategory = 15, font.size = 4, 
            title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
  dev.off()
  
  pdf(file = paste0("./results/","ReactomePAcnetplot.",comparison,".pdf"))
     print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
          vertex.label.cex = 0.75))
   dev.off()
   }
 }

entrez1 <- mapIds(hgu133a.db, 
                 keys = rownames(results),
                 keytype = "PROBEID",
                 colum = "ENTREZID")

reactome <- enrichPathway(gene = entrez1[back], 
                          universe = entrez1[c(back, 
                                              backs)],
                          organism = "human",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.9, 
                          readable = TRUE)
reactome @result$Description <- paste0(str_sub(
    reactome @result$Description, 1, 20),
    "...")

```

Alternative mapping 

```{r results='hide', warning=FALSE, message=FALSE}
##KEGG

entrez <- fit2$genes[, "ENTREZ_GENE_ID"]
enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")
topKEGG(enrich_kegg)

##GO
enrich_go <- goana(fit2, geneid = entrez, species = "Hs")
topGO(enrich_go, ontology = "BP")

## Filtering for genes involved in the biological processes
bpGO <- cancerVSnormTopTable
x <- org.Hs.egGO2ALLEGS
Rkeys(x) <- "GO:0008150"
ID <- mappedLkeys(x)
i <- bpGO$ENTREZID %in% ID
bpGO <- bpGO[i,] 

## Filtering for genes involved in xenobiotic metabolic process
xmpGO <- cancerVSnormTopTable
x1 <- org.Hs.egGO2ALLEGS
Rkeys(x1) <- "GO:0006805"
ID1 <- mappedLkeys(x1)
i1 <- xmpGO$ENTREZID %in% ID1
i1_d <- as.numeric(i1)
xmpGO <- xmpGO[i1,]

## can reapeat the process for each pathway of interest, which will in turn pull the table of all the genes symbols and their values. 



```

Alternative mapping 2.

```{r results='hide', warning=FALSE, message=FALSE}
library(clusterProfiler)
## CC - cellular  component, MF - molecular function, BP - biological process.

ggoCC <- groupGO(gene = entrez,
               OrgDb = org.Hs.eg.db,
               ont = "CC",
               level = 3,
               readable = TRUE)

ggoMF <- groupGO(gene = entrez,
               OrgDb = org.Hs.eg.db,
               ont = "MF",
               level = 3,
               readable = TRUE)

ggoBP <- groupGO(gene = entrez,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               level = 3,
               readable = TRUE)

## Over-representation analysis

egoCC <- enrichGO(gene = entrez,
                  OrgDb = org.Hs.eg.db,
                  universe      = names(geneList), ##
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

egoMF <- enrichGO(gene = entrez,
                  OrgDb = org.Hs.eg.db,
                  universe      = names(geneList), ##
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

egoBP <- enrichGO(gene = entrez,
                  OrgDb = org.Hs.eg.db,
                  universe      = names(geneList), ##
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

gene.df <- bitr(entrez, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

egoCC2 <- enrichGO(gene = gene.df$ENSEMBL,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'ENSEMBL',
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

egoMF2 <- enrichGO(gene = gene.df$ENSEMBL,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'ENSEMBL',
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

egoBP2 <- enrichGO(gene = gene.df$ENSEMBL,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'ENSEMBL',
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

geneList = cancerVSnormTopTable[,4]
names(geneList) = as.character(cancerVSnormTopTable[,3])
geneList = sort(geneList, decreasing = TRUE)

GOaCC <- gseGO(geneList = geneList,
               OrgDb = org.Hs.eg.db,
               ont = "CC",
               minGSSize = 100,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               verbose = FALSE)

```
Alternative mapping 3.
```{r results='hide', warning=FALSE, message=FALSE}


back <- subset(cancerVSnormTopTable, adj.P.Val < 0.1)$PROBEID
backidx <- genefilter::genefinder(filtered_cancervsnorm , as.character(back), method = "manhattan", scale = "none")
backidx <- sapply(backidx, function(x)x$indices)
backs <- featureNames(filtered_cancervsnorm )[backidx]
backs <- setdiff(backs, back)
intersect(backs, back)

IDs <- rownames(cancerVSnormTopTable)
un <- IDs %in% c(backs, back)
sl <- IDs %in% back

genes <- sl[un]
genes <- factor(as.integer(sl[un]))
names(genes) <- IDs[un]

top_Go_data <- new("topGOdata", ontology = "BP", allGenes = genes, nodeSize = 10, annot = annFUN.db, affyLib = "hgu133a.db")

topGoElim <- runTest(top_Go_data, algorithm = "elim", statistic = "Fisher")

topGoClassic <- runTest(top_Go_data, algorithm = "classic", statistic = "Fisher")


resTopGo <- GenTable(top_Go_data, Fisher.elim = topGoElim, Fisher.classic = topGoClassic,
                     orderBy = "Fisher.elim", topNodes = 100)

genesTopGo <- printGenes(top_Go_data, whichTerms = resTopGo$GO.ID, chip = "hgu133a.db", geneCutOff = 1000)

resTopGo$sig_genes <- sapply(genesTopGo, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"], ";"),
        collapse = "")
})


```



```{r}
goplot(egoCC)
cnetplot(egoCC, showCategory = 10)
```

### Linear Discriminant Analysis

```{r results='hide', warning=FALSE, message=FALSE}

packages <- c("knitr")
for (pkg in packages) {
  if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

pkg <- "qvalue"
if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
  source("http://bioconductor.org/biocLite.R")
  biocLite();
  biocLite(pkg)
}
library(pkg, character.only = TRUE)

```

Prepare the data with multiple testing. I will use already cleaned up version of the data that I prepared in the firsts steps of the project. Only thing I will change is the variability cutoff, in order to eliminate more genes. 

```{r results='hide', warning=FALSE, message=FALSE}

LDRfilter <- nsFilter(AQcancervsnormal, 
                      require.entrez = TRUE, 
                      remove.dupEntrez = TRUE,
                      var.filter = TRUE,
                      var.func = IQR,
                      var.cutoff = 0.75,
                      filterByQuantile = TRUE,
                      feature.exclude = "^AFFX")
print(LDRfilter$filter.log)
LDR_canVSnorm <- LDRfilter$eset
nrow(LDR_canVSnorm)
head(LDR_canVSnorm@featureData$"Gene Symbol")


table(LDR_canVSnorm$characteristics_ch1.4)
      
canVSnormal <- exprs(LDR_canVSnorm)

canVSnorm <- as.data.frame(canVSnormal)

vv <- LDR_canVSnorm$characteristics_ch1.4

colnames(canVSnorm) <- paste(colnames(canVSnorm), vv, sep = "_")

```

Single Student Test. Select random gene to test null hyptohesis.

```{r fig.align='right'}

n.Cancer <- 186
n.Normal <- 55

g <- sample(1:nrow(canVSnorm), 1)
g.profile <-as.vector(as.matrix(canVSnorm[g,]))
plot.col <- c('histology: colon cancer' = '#4444BB', 'histology: normal colon mucosa' = '#FFFF88')
par(mar = c(5.1, 4.1, 4.1, 11))
barplot(g.profile, col=plot.col[vv], main = "One random gene expression in the data", ylab = "Relative expression", xlab = "Individual samples n=65" , xlim = c(1, 65))
legend('topright', c("histology: colon cancer", "histology: normal colon mucosa"), col=plot.col[c("histology: colon cancer", "histology: normal colon mucosa")], pch=15, bty="o", bg='white', cex = 0.65)

sample.Cancer <- g.profile[vv == "histology: colon cancer"]
sample.Normal <- g.profile[vv == "histology: normal colon mucosa"]

mean.est.Cancer <- mean(sample.Cancer)
mean.est.Normal <- mean(sample.Normal)

sample.sd.Cancer <- sd(sample.Cancer) * sqrt((n.Cancer-1)/n.Cancer)
sample.sd.Normal <- sd(sample.Normal) * sqrt((n.Normal-1)/n.Normal)

sd.est.Cancer <- sd(sample.Cancer)
sd.est.Normal <- sd(sample.Normal)


sd.err.est.Cancer <- sd(sample.Cancer) / sqrt(n.Cancer)
sd.err.est.Normal <- sd(sample.Normal) / sqrt(n.Normal)

diff.sd.est <- sqrt((n.Cancer * sample.sd.Cancer^2 + n.Normal * sample.sd.Normal^2) * (1/n.Cancer + 1/n.Normal) / (n.Cancer+n.Normal-2))

d <- abs(mean.est.Cancer - mean.est.Normal)
t.obs.Student <- d/diff.sd.est
P.val.Student <- 2*pt(q=t.obs.Student,
                      df=n.Cancer+n.Normal-2,
                      lower.tail = F)
t.student <- t.test(sample.Cancer, sample.Normal, var.equal = TRUE)

print(t.student)

t.welch <- t.test(sample.Cancer, sample.Normal, var.equal = FALSE)
 
print(t.welch)

t.statistics <- vector()
P.values <- vector()

for (g in 1:nrow(canVSnorm)) {
  print(paste("Random gene", g))
  g.profile <- as.vector(canVSnorm[g,])
  sample.Cancer <- g.profile[vv == "histology: colon cancer"]
  sample.Normal <- g.profile[vv == "histology: normal colon mucosa"]
  t <- t.test(sample.Cancer, sample.Normal)
  t.statistics <- append(t.statistics, t$statistic)
  P.values <- append(P.values, t$p.value)
}
print(P.values)

canvsnorm.t.result <- t.test.multi(canVSnorm, vv)
dim(canvsnorm.t.result)
names(canvsnorm.t.result)
sum(canvsnorm.t.result$E.value <=1)
canvsnorm.E <- canVSnorm[canvsnorm.t.result$E.value <=1,]

write.table(canvsnorm.E, file = file.path(dir.results, "Evalue_Table"))

```

```{r results='hide', warning=FALSE, message=FALSE}

dir.base <- 'http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics'
dir.base.ASG1 <- 'http://pedagogix-tagc.univ-mrs.fr/courses/ASG1'

source(file.path(dir.base, 'R-files', 'config.R'))

setwd(dir.results)
print(paste("Result directory", dir.results))

packages <- c("knitr")
for (pkg in packages) {
  if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

pkg <- "qvalue"
if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
  source("http://bioconductor.org/biocLite.R")
  biocLite();
  biocLite(pkg)
}
library(pkg, character.only = TRUE)


url.stats4bioinfo <- "http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics"
source(file.path(url.stats4bioinfo, 'R-files/config.R'))

source(file.path(url.stats4bioinfo, 'R-files/util/util_student_test_multi.R'))

CVSNMpheno <- phenoData(LDR_canVSnorm)
colnames(CVSNMpheno)
sample.labels <- as.vector(CVSNMpheno@data$characteristics_ch1.4)
sample.colors <- plot.col
group.descriptions <- rbind(sample.labels, sample.colors[vv])
group.descriptions <- t(group.descriptions)
group.descriptions <- as.data.frame(group.descriptions)
group.labels <- group.descriptions$sample.labels
names(group.labels) <- row.names(group.descriptions)
group.colors <- group.descriptions$V2
names(group.colors) <- row.names(group.descriptions)
group.legend <- paste(sep ="", group.labels, row.names(group.descriptions))

## Variance per gene

var.per.gene <- apply(canVSnorm, 1, var)
sd.per.gene <- apply(canVSnorm, 1, sd)

```

```{r}

par(mar = c(2, 2, 2, 2))
par(mfrow=c(1,2))

hist(var.per.gene, breaks=100, col="#BBFFDD", main="Variance per gene", xlab="Variance", ylab="Number of genes")
hist(sd.per.gene, breaks=100, col="#BBFFDD", main="Standard dev. per gene", xlab="Standard deviation", ylab="Number of genes")

par(mar = c(2, 2, 2, 2))
par(mfrow=c(1,2))
```

```{r results='hide', warning=FALSE, message=FALSE}

## Select 20 top-ranking gens for training

genes.by.decr.var <- sort(var.per.gene,decreasing=TRUE)
top.nb <- 20
genes.selected.by.var <- names(genes.by.decr.var[1:top.nb])

## Rank by cross-sample variance

gene.ranks <- data.frame(var=var.per.gene)
gene.ranks$var.rank <- rank(-gene.ranks$var, ties.method='random')
head(gene.ranks, n=10)

```

```{r}
kable(gene.ranks[names(genes.by.decr.var[1:5]),], caption = "Five genes with the  highest variance.")

kable(gene.ranks[names(tail(genes.by.decr.var)),], caption = "Five genes with the  lowest variance.")

```

```{r results='hide', warning=FALSE, message=FALSE}
g1 <- 236
g2 <- 1213

maxvar.g1 <- names(genes.by.decr.var[1])
maxvar.g2 <- names(genes.by.decr.var[2])

x <- as.vector(as.matrix(canVSnorm[maxvar.g1,]))
y <- as.vector(as.matrix(canVSnorm[maxvar.g2,]))

plot(x,y,
      col=sample.colors,
      type='n',
      panel.first=grid(col='black'), 
      main="2 genes with the highest variance", 
      xlab=paste('gene', maxvar.g1), 
      ylab=paste('gene', maxvar.g2))

text(x, y,labels=sample.labels,col=sample.colors,pch=0.5)
legend('topright',col=group.colors, 
         legend=group.legend,pch=1,cex=0.6,bg='white',bty='o')

var.per.gene <- apply(canVSnorm, 1, var)
sd.per.gene <- apply(canVSnorm, 1, sd)


par(mfrow=c(1,2))
boxplot(x ~ sample.labels, las=1,
        horizontal=TRUE,
        main=maxvar.g1, 
        col="#BBBBBB")
data.frame(group.colors)

boxplot(y ~ sample.labels, las=1,
        horizontal=TRUE,
        main=maxvar.g2, col="#BBBBBB")
par(mfrow=c(1,1))


## Load Welch's

source(file.path(dir.util, "util_student_test_multi.R"))
group.of.interest <- "histology: colon cancer"
one.vs.others <- sample.labels
## one.vs.others[sample.labels != group.of.interest] <- "histology: normal colon mucosa"
print(table(one.vs.others))


welch.one.vs.others <- t.test.multi(canVSnorm, 
                                    one.vs.others,
                                    volcano.plot = FALSE)

kable(head(welch.one.vs.others), caption = "Head of the Welch result table. Each row corrresponds to one probeset (gene), each column to one statistics used for the Welch test.")

test.name <- paste(group.of.interest, '.vs.others.sig', sep='')
gene.ranks[,test.name] <- welch.one.vs.others$sig
gene.ranks[,paste(test.name, ".rank", sep="")] <- 
    rank(-welch.one.vs.others$sig, ties.method='random')

group.of.interest1 <- as.data.frame(group.of.interest)
ample.labels1 <- as.data.frame(sample.labels)
one.vs.others2 <- as.vector(one.vs.others)

## Apply the Welch test for the 3 other majority groups
## Not necessary for the cancer vs normal dataset, as there are only two groups that we are looking at.
for (group.of.interest1 in c("histology: colon cancer", "histology: normal colon mucosa")) {
    print(paste("Selecting differentially expressed genes for", group.of.interest1, "versus others"))
    one.vs.others2 <- sample.labels1
    one.vs.others2[sample.labels1 != group.of.interest1] <- "histology: normal colon mucosa"
    
    welch.one.vs.others <- t.test.multi(canVSnorm,                           one.vs.others2,
                         volcano.plot = FALSE)

    test.name <- paste(group.of.interest1, '.vs.others.sig', sep='')
    gene.ranks[,test.name] <- welch.one.vs.others$sig
    gene.ranks[,paste(test.name, ".rank", sep="")] <- 
      rank(-welch.one.vs.others$sig, ties.method='random')
}


head(gene.ranks)

write.table(gene.ranks, file=file.path(dir.results, 'gene_ranks.tab'), sep='\t', quote=F, col.names=NA)
 

par(mfrow=c(1,1))


## ANOVA ordering

g <- 1234 ## random gene

g.expr <- unlist(canVSnorm[g,])

g.for.anova <- data.frame("expr"=g.expr, "group"=sample.labels)

g.aov.result <- aov(formula = expr ~ group, data = g.for.anova)

print(g.aov.result)

g.anova.result <- anova(lm(formula = expr ~ group, data = g.for.anova))

print(g.anova.result)

attributes(g.anova.result)

pval <- as.numeric(unlist(g.anova.result)["Pr(>F)1"])
print(pval)

eval <- pval * nrow(canVSnorm)

print(eval)

g.anova.summary <- data.frame("g"=g, 
                     "name"=row.names(canVSnorm[g,]),
                     "pval"=pval,
                     "eval"=eval,
                     "sig"=-log(eval, base=10))
kable(g.anova.summary, caption = "Anova result for an aribtrary gene. ")


## PCA

expr.prcomp <- prcomp(t(canVSnorm))

attributes(expr.prcomp) 
names(expr.prcomp) 

plot(expr.prcomp, xlab='Component', col="#BBDDFF")

sd.per.pc <- expr.prcomp$sdev

var.per.pc <- sd.per.pc^2

sd.per.pc.percent <- sd.per.pc/sum(sd.per.pc)

var.per.pc.percent <- var.per.pc/sum(var.per.pc)

barplot(var.per.pc.percent[1:10], main='Can vs norm, Percent of variance  per component', xlab='Component', ylab='Percent variance', col='#BBDDFF')


## LINEAR DISCRIMINANT ANALYSIS LDA

library(MASS)

one.vs.others.lda.allvars <- lda(t(canVSnorm), one.vs.others, CV=FALSE)

## Incorrect usage - as we used a training set to train on itself.

predict.lda.allvars <- predict(object = one.vs.others.lda.allvars, newdata = t(canVSnorm))
hits <- sum(one.vs.others == predict.lda.allvars$class)
errors <- sum(one.vs.others != predict.lda.allvars$class)
total <- hits + errors
(hit.rate.internal <- hits / total)
(error.rate.internal <- errors / total)

## Leave-one-out

one.vs.others.lda.allvars.loo <- lda(t(canVSnorm),one.vs.others,CV=TRUE)

(hit.rate.loo <- sum(one.vs.others == one.vs.others.lda.allvars.loo$class) / total)
(error.rate.loo <- 1 - hit.rate.loo) ## classification rate falls with extra testing

## Random expectation for the hit rate

random.hit.rates <- vector()
for (rep in 1:10000) {
  random.hit.rates <- append(random.hit.rates, sum(one.vs.others == sample(one.vs.others)) / total)
}
(random.hit.rates.mean <- mean(random.hit.rates))

prior <- as.vector(table(one.vs.others))/length(one.vs.others)
(hit.rate.expect <- sum(prior^2))

hist(random.hit.rates, breaks=(0:total)/total, col="lightgrey", 
     freq=TRUE,
     main="Hit rate analysis",
     xlab="Hit rate",
     ylab="Frequency")
arrows(x0=hit.rate.loo, y0 = 1000, x1=hit.rate.loo, y1=100, 
       col="darkgreen", lwd=2, code=2, , length=0.2, angle=20)
arrows(x0=hit.rate.internal, y0 = 1000, x1=hit.rate.internal, y1=100, 
       col="red", lwd=2, code=2, , length=0.2, angle=20)
arrows(x0=hit.rate.expect, y0 = 1000, x1=hit.rate.expect, y1=100, 
       col="darkblue", lwd=2, code=2, , length=0.2, angle=20)

legend("topleft", legend=c("random", "random expectation", "LOO", "internal"), 
       col=c("grey", "darkblue", "darkgreen", "red"), 
       lwd=c(4,2,2,2))


## Selecting a subset of the variables (feature selection)

welch.Bo.vs.others <- t.test.multi(canVSnorm, 
                                   one.vs.others,
                                   volcano.plot = FALSE)
print(rownames(welch.Bo.vs.others[welch.Bo.vs.others$rank <= 20,]))

welch.Bo.vs.others.sorted <- welch.Bo.vs.others[order(welch.Bo.vs.others$sig, decreasing=TRUE),]

sorted.names <- rownames(welch.Bo.vs.others.sorted)
print(sorted.names[1:20])

welch.Bo.vs.others[sorted.names[0:20], c("E.value","sig")]

g1 <- sorted.names[1]
g2 <- sorted.names[2]
x <- as.vector(as.matrix(canVSnorm[g1,]))
y <- as.vector(as.matrix(canVSnorm[g2,]))

top.variables <- 20
selected.genes <- sorted.names[1:top.variables]

one.vs.others.lda.classifier <- lda(t(canVSnorm[selected.genes,]),one.vs.others,CV=FALSE) 

one.vs.others.lda.loo <- lda(t(canVSnorm[selected.genes,]),one.vs.others,CV=TRUE) 

one.vs.others.loo.predicted.class <- as.vector(one.vs.others.lda.loo$class)

(one.vs.others.lda.loo.xtab <- table(one.vs.others, one.vs.others.loo.predicted.class))

(one.vs.others.lda.loo.xtab2 <- table(sample.labels, one.vs.others.loo.predicted.class))


library(lattice)

levelplot(one.vs.others.lda.loo.xtab)

hits <- one.vs.others == one.vs.others.loo.predicted.class
errors <- one.vs.others  != one.vs.others.loo.predicted.class

(nb.hits <- sum(na.omit(hits)))

(nb.pred <- length(na.omit(hits)))

(hit.rate <- nb.hits / nb.pred )


## Multi-group classification

multigroup.lda.loo <- lda(t(canVSnorm[selected.genes,]),sample.labels,CV=TRUE) 

multigroup.loo.predicted.class <- as.vector(multigroup.lda.loo$class)

table(multigroup.loo.predicted.class)

(multigroup.lda.loo.xtab <- table(sample.labels, multigroup.loo.predicted.class))

library(lattice)
levelplot(multigroup.lda.loo.xtab)

hits <- sample.labels == multigroup.loo.predicted.class

errors <- sample.labels != multigroup.loo.predicted.class

(nb.hits <- sum(na.omit(hits)))

(nb.pred <- length(na.omit(hits)))

(hit.rate <- nb.hits / nb.pred )


## Random expectation of a hit rate

n <- length(sample.labels)
n.goi<- sum(sample.labels == group.of.interest1)
n.others <- sum(sample.labels != group.of.interest1)
prior <- c("goi"=n.goi/n,
           "others" = n.others/n)

exp.hits <- c("goi"=(n.goi*prior["goi"]),
              "other"=(n.others*prior["others"]))

print(exp.hits)

exp.hit.rate <- sum(prior^2)

print(exp.hit.rate)

multi.samples.per.class <- unlist(table(sample.labels))

multi.prior <- (multi.samples.per.class)/sum(multi.samples.per.class)
(multi.expect.hit.rate <- sum(multi.prior^2))

## Training a classifier with permuted labels

sample.labels.perm <- as.vector(sample(sample.labels))

table(sample.labels, sample.labels.perm)

lda.loo.labels.perm <- lda(t(canVSnorm[selected.genes,]),sample.labels.perm,CV=TRUE) 

loo.predicted.class.labels.perm <- as.vector(lda.loo.labels.perm$class)
lda.loo.labels.perm.xtab <- table(sample.labels.perm, loo.predicted.class.labels.perm)

hits.label.perm <- sample.labels.perm == loo.predicted.class.labels.perm
(nb.hits.label.perm <- sum(na.omit(hits.label.perm)))

(nb.pred.label.perm <- length(na.omit(hits.label.perm)))

(hit.rate.label.perm <- nb.hits.label.perm / nb.pred.label.perm )


## Label permutation test for two-group classification

sample.labels.perm.2gr <- as.vector(sample(one.vs.others))
table(one.vs.others, sample.labels.perm.2gr)

permuted.equal <- sum(diag(table(one.vs.others, sample.labels.perm.2gr)))
(permuted.equal.rate <- permuted.equal/length(one.vs.others))

lda.loo.labels.perm.2gr <- lda(t(canVSnorm[selected.genes,]),sample.labels.perm.2gr,CV=TRUE) 

loo.predicted.class.labels.perm.2gr <- as.vector(lda.loo.labels.perm.2gr$class)
lda.loo.labels.perm.2gr.xtab <- table(sample.labels.perm.2gr, loo.predicted.class.labels.perm.2gr)
print(lda.loo.labels.perm.2gr.xtab)


(hit.rate.label.perm.2gr <- sum(diag(lda.loo.labels.perm.2gr.xtab)) / sum(lda.loo.labels.perm.2gr.xtab))


## Quadratic discriminant analys QDA

(n.variables <- length(selected.genes))

one.vs.others.qda.loo <- qda(t(canVSnorm[selected.genes,]),one.vs.others,CV=TRUE)

table(one.vs.others, one.vs.others.qda.loo$class)

(one.vs.others.qda.loo.hit.rate <- sum(one.vs.others == one.vs.others.qda.loo$class)/n)

label.perm.one.vs.others.qda.loo <- qda(t(canVSnorm[selected.genes,]),
                                        sample(one.vs.others, replace=FALSE),  CV=TRUE) 
  

table(one.vs.others, label.perm.one.vs.others.qda.loo$class)

(label.perm.one.vs.others.qda.loo.hit.rate <- 
   sum(one.vs.others == label.perm.one.vs.others.qda.loo$class)/n)

one.vs.others.lda.loo <- lda(t(canVSnorm[selected.genes,]),one.vs.others,CV=TRUE) 

(one.vs.others.lda.loo.hit.rate <- sum(one.vs.others == one.vs.others.lda.loo$class)/n)


write.csv(canVSnorm, "C:/Users/greta/Documents/LDR/data.csv", row.names = TRUE)
types <- as.data.frame(vv)
write.csv(types, "C:/Users/greta/Documents/LDR/types.csv", row.names = TRUE)

```


#### Connection to Python

Installation of package 'reticulate' and loading the library is necessary. To connect the Python packages scikit-learn, pandas, numpy and matplotlib, first the download of Anaconda to local computer needs to be done. During installation, make sure to tick the 'Add Anaconda3 to my PATH environment variable' as it is critical to obtain the environment in RStudio. After the installation is complete, in Terminal type the following:
```
library(reticulate)
use_condaenv("py3.8", required = TRUE)

py_run_string("import os as os")
py_run_string("os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = 'C:/Users/2015524/Anaconda3/envs/py3.8/Library/plugins/platforms'")

## ANACONDA TERMINAL


conda create -n py3.8 python=3.8 scikit-learn pandas numpy matplotlib
```
Always check the current and newest version of conda. Connect the environment and check if the connection is active with `conda env list` in the Terminal. 
Back in Console, after loading the package, check the connectivity.
```
reticulate::conda_list()
```
Next, set your environment. The connection of python syntax will be indicated by '>>>' instead of '>' in the console.

In the console, run

```
conda activate py3.8
```
In Rstudio console:conda create -n py3.8 python=3.8 scikit-learn pandas numpy matplotlib
```
Always check the current and newest version of conda. Connect the environment and check if the connection is active with `conda env list` in the Terminal. 
Back in Console, after loading the package, check the connectivity.
```
reticulate::conda_list()
```
Next, set your environment. The connection of python syntax will be indicated by '>>>' instead of '>' in the console.

In the console, run

```
conda activate py3.8
```
In Rstudio console:
```
use_condaenv
use_condaenv("C:/Users/greta/anaconda3/envs/py3.8", required = TRUE)

py_config() ## for double checking the connection

repl_python()


```
use_condaenv("py3.8", required = TRUE)
py_config() ## for double checking the connection

```

```{r}
library(caret)

```



```{python}
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn
import GEOparse

gse = GEOparse.get_GEO(geo="GSE68468", destdir="./")

for gsm_name, gsm in gse.gsms.items():
    print("Name: ", gsm_name)
    print("Metadata:",)
    for key, value in gsm.metadata.items():
        print(" - %s : %s" % (key, ", ".join(value)))
    print ("Table data:",)
    print (gsm.table.head())
    break
  
for gpl_name, gpl in gse.gpls.items():
    print("Name: ", gpl_name)
    print("Metadata:",)
    for key, value in gpl.metadata.items():
        print(" - %s : %s" % (key, ", ".join(value)))
    print("Table data:",)
    print(gpl.table.head())
    break
  
  

df = pd.read_excel(r'C:/Users/greta/Documents/LDR/data.xlsx')
vv = pd.read_excel(r'C:/Users/greta/Documents/LDR/types.xlsx') 
vv = vv.iloc[:,1:]
vv.head()
 
## SVM classifier

array = df.values
array1 = vv.values

X = array[:,1:241]
Y = array[:,241]

from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression

model = LogisticRegression()
rfe = RFE(model, 3)
fit = rfe.fit(X, Y)


```


