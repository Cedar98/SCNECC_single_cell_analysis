install.packages("BiocManager")
BiocManager::install("GSVA")
library(GSVA)

p <- 10000 ## number of genes
n <- 30    ## number of samples
## simulate expression values from a standard Gaussian distribution
X <- matrix(rnorm(p*n), nrow=p,
            dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
X[1:5, 1:5]

## sample gene set sizes
gs <- as.list(sample(10:100, size=100, replace=TRUE))
## sample gene sets
gs <- lapply(gs, function(n, p)
  paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))

gsvaPar <- gsvaParam(X, gs)
gsvaPar

gsva.es <- gsva(gsvaPar, verbose=FALSE)
dim(gsva.es)

gsva.es[1:5, 1:5]

class(gs)

length(gs)

head(lapply(gs, head))

library(org.Hs.eg.db)

goannot <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="GO")
head(goannot)

genesbygo <- split(goannot$ENTREZID, goannot$GO)
length(genesbygo)

head(genesbygo)

library(GSEABase)
library(GSVAdata)

data(c2BroadSets)
class(c2BroadSets)

attr(,"package")

c2BroadSets
GeneSetCollection

library(Biobase)

data(commonPickrellHuang)

stopifnot(identical(featureNames(huangArrayRMAnoBatchCommon_eset),
                    featureNames(pickrellCountsArgonneCQNcommon_eset)))
stopifnot(identical(sampleNames(huangArrayRMAnoBatchCommon_eset),
                    sampleNames(pickrellCountsArgonneCQNcommon_eset)))

canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      grep("^REACTOME", names(c2BroadSets)),
                                      grep("^BIOCARTA", names(c2BroadSets)))]
canonicalC2BroadSets
GeneSetCollection

data(genderGenesEntrez)

MSY <- GeneSet(msYgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"),
               setName="MSY")
MSY
setName: MSY 
geneIdType: EntrezId
collectionType: Broad
bcCategory: c2 (Curated)
bcSubCategory: NA
XiE <- GeneSet(XiEgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"),
               setName="XiE")
XiE
setName: XiE 
geneIdType: EntrezId
collectionType: Broad
bcCategory: c2 (Curated)
bcSubCategory: NA

canonicalC2BroadSets <- GeneSetCollection(c(canonicalC2BroadSets, MSY, XiE))
canonicalC2BroadSets
GeneSetCollection

huangPar <- gsvaParam(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets,
                      minSize=5, maxSize=500)
esmicro <- gsva(huangPar)
pickrellPar <- gsvaParam(pickrellCountsArgonneCQNcommon_eset,
                         canonicalC2BroadSets, minSize=5, maxSize=500,
                         kcdf="Poisson")
esrnaseq <- gsva(pickrellPar)

library(edgeR)

lcpms <- cpm(exprs(pickrellCountsArgonneCQNcommon_eset), log=TRUE)

genecorrs <- sapply(1:nrow(lcpms),
                    function(i, expmicro, exprnaseq)
                      cor(expmicro[i, ], exprnaseq[i, ], method="spearman"),
                    exprs(huangArrayRMAnoBatchCommon_eset), lcpms)
names(genecorrs) <- rownames(lcpms)

pwycorrs <- sapply(1:nrow(esmicro),
                   function(i, esmicro, esrnaseq)
                     cor(esmicro[i, ], esrnaseq[i, ], method="spearman"),
                   exprs(esmicro), exprs(esrnaseq))
names(pwycorrs) <- rownames(esmicro)

par(mfrow=c(1, 2), mar=c(4, 5, 3, 2))
hist(genecorrs, xlab="Spearman correlation",
     main="Gene level\n(RNA-seq log-CPMs vs microarray RMA)",
     xlim=c(-1, 1), col="grey", las=1)
hist(pwycorrs, xlab="Spearman correlation",
     main="Pathway level\n(GSVA enrichment scores)",
     xlim=c(-1, 1), col="grey", las=1)

par(mfrow=c(1, 2))
rmsy <- cor(exprs(esrnaseq)["MSY", ], exprs(esmicro)["MSY", ])
plot(exprs(esrnaseq)["MSY", ], exprs(esmicro)["MSY", ],
     xlab="GSVA scores RNA-seq", ylab="GSVA scores microarray",
     main=sprintf("MSY R=%.2f", rmsy), las=1, type="n")
fit <- lm(exprs(esmicro)["MSY", ] ~ exprs(esrnaseq)["MSY", ])
abline(fit, lwd=2, lty=2, col="grey")
maskPickrellFemale <- pickrellCountsArgonneCQNcommon_eset$Gender == "Female"
maskHuangFemale <- huangArrayRMAnoBatchCommon_eset$Gender == "Female"
points(exprs(esrnaseq["MSY", maskPickrellFemale]),
       exprs(esmicro)["MSY", maskHuangFemale],
       col="red", pch=21, bg="red", cex=1)
maskPickrellMale <- pickrellCountsArgonneCQNcommon_eset$Gender == "Male"
maskHuangMale <- huangArrayRMAnoBatchCommon_eset$Gender == "Male"
points(exprs(esrnaseq)["MSY", maskPickrellMale],
       exprs(esmicro)["MSY", maskHuangMale],
       col="blue", pch=21, bg="blue", cex=1)
legend("topleft", c("female", "male"), pch=21, col=c("red", "blue"),
       pt.bg=c("red", "blue"), inset=0.01)
rxie <- cor(exprs(esrnaseq)["XiE", ], exprs(esmicro)["XiE", ])
plot(exprs(esrnaseq)["XiE", ], exprs(esmicro)["XiE", ],
     xlab="GSVA scores RNA-seq", ylab="GSVA scores microarray",
     main=sprintf("XiE R=%.2f", rxie), las=1, type="n")
fit <- lm(exprs(esmicro)["XiE", ] ~ exprs(esrnaseq)["XiE", ])
abline(fit, lwd=2, lty=2, col="grey")
points(exprs(esrnaseq["XiE", maskPickrellFemale]),
       exprs(esmicro)["XiE", maskHuangFemale],
       col="red", pch=21, bg="red", cex=1)
points(exprs(esrnaseq)["XiE", maskPickrellMale],
       exprs(esmicro)["XiE", maskHuangMale],
       col="blue", pch=21, bg="blue", cex=1)
legend("topleft", c("female", "male"), pch=21, col=c("red", "blue"),
       pt.bg=c("red", "blue"), inset=0.01)

data(gbm_VerhaakEtAl)
gbm_eset
ExpressionSet (storageMode: lockedEnvironment)

protocolData: none
phenoData

varLabels: subtype
featureData: none

Annotation:  
  head(featureNames(gbm_eset))

table(gbm_eset$subtype)

data(brainTxDbSets)
lengths(brainTxDbSets)

lapply(brainTxDbSets, head)

gbmPar <- gsvaParam(gbm_eset, brainTxDbSets, maxDiff=FALSE)
gbm_es <- gsva(gbmPar)

library(RColorBrewer)
subtypeOrder <- c("Proneural", "Neural", "Classical", "Mesenchymal")
sampleOrderBySubtype <- sort(match(gbm_es$subtype, subtypeOrder),
                             index.return=TRUE)$ix
subtypeXtable <- table(gbm_es$subtype)
subtypeColorLegend <- c(Proneural="red", Neural="green",
                        Classical="blue", Mesenchymal="orange")
geneSetOrder <- c("astroglia_up", "astrocytic_up", "neuronal_up",
                  "oligodendrocytic_up")
geneSetLabels <- gsub("_", " ", geneSetOrder)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol <- hmcol[length(hmcol):1]

heatmap(exprs(gbm_es)[geneSetOrder, sampleOrderBySubtype], Rowv=NA,
        Colv=NA, scale="row", margins=c(3,5), col=hmcol,
        ColSideColors=rep(subtypeColorLegend[subtypeOrder],
                          times=subtypeXtable[subtypeOrder]),
        labCol="", gbm_es$subtype[sampleOrderBySubtype],
        labRow=paste(toupper(substring(geneSetLabels, 1,1)),
                     substring(geneSetLabels, 2), sep=""),
        cexRow=2, main=" \n ")
par(xpd=TRUE)
text(0.23,1.21, "Proneural", col="red", cex=1.2)
text(0.36,1.21, "Neural", col="green", cex=1.2)
text(0.47,1.21, "Classical", col="blue", cex=1.2)
text(0.62,1.21, "Mesenchymal", col="orange", cex=1.2)
mtext("Gene sets", side=4, line=0, cex=1.5)
mtext("Samples          ", side=1, line=4, cex=1.5)

data(leukemia)
leukemia_eset
ExpressionSet (storageMode: lockedEnvironment)

protocolData

varLabels: ScanDate
varMetadata: labelDescription
phenoData

varLabels: subtype

featureData: none

Annotation: hgu95a 

cgpC2BroadSets <- c2BroadSets[c(grep("_UP$", names(c2BroadSets)),
                                grep("_DN$", names(c2BroadSets)))]
cgpC2BroadSets
leukPar <- gsvaParam(leukemia_eset, cgpC2BroadSets,
                     minSize=10, maxSize=500)
leukemia_es <- gsva(leukPar)

class(leukemia_es)

attr(,"package")

leukemia_es
ExpressionSet (storageMode: lockedEnvironment)


protocolData: none
phenoData


varLabels: subtype

featureData: none

Annotation:  
  head(featureNames(leukemia_es))

library(limma)

mod <- model.matrix(~ factor(leukemia_es$subtype))
colnames(mod) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.01)
summary(res)

gssizes <- geneSetSizes(leukemia_es)
plot(sqrt(gssizes), sqrt(fit$sigma), xlab="Sqrt(gene sets sizes)",
     ylab="Sqrt(standard deviation)", las=1, pch=".", cex=3)

fit <- eBayes(fit, trend=gssizes)
res <- decideTests(fit, p.value=0.01)
summary(res)

tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]
plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75),
     main="", xlab="GSVA enrichment score difference", las=1,
     ylab=expression(-log[10]~~Raw~P-value))
abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])),
       col=grey(0.5), lwd=1, lty=2)
points(tt$logFC[match(DEpwys, rownames(tt))],
       -log10(tt$P.Value[match(DEpwys, rownames(tt))]),
       pch=".", cex=5, col="darkred")
text(max(tt$logFC)*0.85, -log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])),
     "1% FDR", pos=3)

DEpwys_es <- exprs(leukemia_es[DEpwys, ])
colorLegend <- c("darkred", "darkblue")
names(colorLegend) <- c("ALL", "MLL")
sample.color.map <- colorLegend[pData(leukemia_es)[, "subtype"]]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                            method="complete")
heatmap(DEpwys_es, ColSideColors=sample.color.map, xlab="samples",
        ylab="Pathways", margins=c(2, 20),
        labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "",
                                          rownames(DEpwys_es))), 1, 35),
        labCol="", scale="row", Colv=as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering))
legend("topleft", names(colorLegend), fill=colorLegend, inset=0.01, bg="white")

res <- igsva()


