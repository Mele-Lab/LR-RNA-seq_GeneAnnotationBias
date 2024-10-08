## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/07_differential_expressions"
##
##    load setup sources
mn5source <- "/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/Z_resources/myutils.R"
localsource <-"/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/Z_resources/myutils.R"
ifelse(file.exists(mn5source), source(mn5source), source(localsource))
rm(mn5source, localsource)
##    sets wd, load tidyverse, data.table, sets DT threads
setup_script(relative_path, 3, 48)
##
##    catch arguments from standard input 
##      catch_args(number_of_args, "obj_name1", "obj_name2")
catch_args(0)
##
## 0----------------------------END OF HEADER----------------------------------0

# load data
metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==FALSE]

annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcripidv_geneidv_match.tsv")

counts <- fread("../../novelannotations/v47_kallisto_quant/matrix.abundance.tsv")
counts <- annot[counts, on=c("transcriptid.v"="transcript_id")]
counts <- counts[, lapply(.SD, sum), by = geneid.v, .SDcols = patterns("_")]
counts <- column_to_rownames(counts, var="geneid.v")
colnames(counts) <- gsub("_.*", "", colnames(counts))
samplenames <- metadata$sample
names(samplenames) <- metadata$cell_line_id
colnames(counts) <- samplenames[colnames(counts)]
metadata <- metadata[, .(sample, sex, ooa)]
setcolorder(counts, sort(colnames(counts)))
setorder(metadata, sample)

# library(tweeDEseqCountData)
# data(pickrell1)
# counts<-exprs(pickrell1.eset)


# data(genderGenes)
# data(annotEnsembl63)
# annot <- annotEnsembl63[,c("Symbol","Chr")]
# rm(annotEnsembl63)
# We now have the counts, gender of each sample and annotation (gene symbol and chromosome) for each Ensemble gene. We can form a DGElist object using the edgeR package.
library(missMethyl)
library(edgeR)
y <- DGEList(counts=counts) # potser em demanar que li doni els gens també

# drop lowly expressed genes by keeping genes with at least 1 count per million reads in at least 20 samples. Finally we perform scaling normalisation.
isexpr <- rowSums(cpm(y)>1) >= 5
y <- y[isexpr,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

#We set up the design matrix and test for differential variability.

design <- model.matrix(~sex+ooa, metadata) # does it accept mixed models?
fitvar <- varFit(y, design = design)
## Converting counts to log counts-per-million using voom.


  summary(decideTests(fitvar.hapmap))
##        (Intercept) gendermale
## Down             0          2
## NotSig           0      17308
## Up           17310          0
topDV <- topVar(fitvar,coef=3)

##                       Symbol Chr  SampleVar LogVarRatio DiffLevene         t
## ENSG00000213318 RP11-331F4.1  16 5.69839463   -2.562939 -0.9859943 -8.031243
## ENSG00000129824       RPS4Y1   Y 2.32497726   -2.087025 -0.4585620 -4.957005
## ENSG00000233864       TTTY15   Y 6.79004140   -2.245369 -0.6085233 -4.612934
## ENSG00000176171        BNIP3  10 0.41317384    1.199292  0.3632133  4.219404
## ENSG00000197358      BNIP3P1  14 0.39969125    1.149754  0.3353288  4.058147
## ENSG00000025039        RRAGD   6 0.91837213    1.091229  0.4926839  3.977022
## ENSG00000103671        TRIP4  15 0.07456448   -1.457139 -0.1520583 -3.911300
## ENSG00000171100         MTM1   X 0.44049558   -1.133295 -0.3334619 -3.896490
## ENSG00000149476          DAK  11 0.13289523   -1.470460 -0.1919880 -3.785893
## ENSG00000064886       CHI3L2   1 2.70234584    1.468059  0.8449434  3.782010
##                      P.Value  Adj.P.Value
## ENSG00000213318 7.238039e-12 1.252905e-07
## ENSG00000129824 3.960855e-06 3.428120e-02
## ENSG00000233864 1.496237e-05 8.633290e-02
## ENSG00000176171 6.441668e-05 2.787632e-01
## ENSG00000197358 1.147886e-04 3.973982e-01
## ENSG00000025039 1.527695e-04 4.375736e-01
## ENSG00000103671 1.921104e-04 4.375736e-01
## ENSG00000171100 2.022293e-04 4.375736e-01
## ENSG00000149476 2.956364e-04 4.425050e-01
## ENSG00000064886 2.995692e-04 4.425050e-01

genesDV <- rownames(topDV.hapmap)
# jitter expresison
