# -----------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
   install.packages("BiocManager")
BiocManager::install("stageR")

## -----------------------------------------------------------------------------
library(stageR)

## ----echo=TRUE,warning=FALSE--------------------------------------------------
library(edgeR) ; library(Biobase) ; library(limma) ; library(utils) 

## -----------------------------------------------------------------------------
data(hammer.eset, package="stageR")
eset <- hammer.eset ; rm(hammer.eset)

## -----------------------------------------------------------------------------
pData(eset)$Time #typo. Will do it ourself
time <- factor(rep(c("mo2","w2"),each=4),levels=c("w2","mo2"))
pData(eset)$protocol
treat <- factor(c("control","control","SNL","SNL","control","control","SNL","SNL"),levels=c("control","SNL"))
design <- model.matrix(~time*treat)
rownames(design) = paste0(time,treat,rep(1:2,4))
colnames(design)[4] = "timeMo2xTreatSNL"
design

## -----------------------------------------------------------------------------
cpmOffset <- 2
keep <- rowSums(cpm(exprs(eset))>cpmOffset)>=2 #2cpm in 2 samples
dge <- DGEList(exprs(eset)[keep,])
colnames(dge) = rownames(design)
dge <- calcNormFactors(dge)

## -----------------------------------------------------------------------------
## regular analysis
voomObj <- voom(dge,design,plot=TRUE)
fit <- lmFit(voomObj,design)
contrast.matrix <- makeContrasts(treatSNL, treatSNL+timeMo2xTreatSNL, timeMo2xTreatSNL, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res <- decideTests(fit2)
summary.TestResults(res) #nr of significant up-/downregulated genes
colSums(summary.TestResults(res)[c(1,3),]) #total nr of significant genes

topTable(fit2)
topTableF(fit2)

## -----------------------------------------------------------------------------
uniqueGenesRegular <- which(res[,1]!=0 | res[,2]!=0 | res[,3]!=0)
length(uniqueGenesRegular) #total nr of significant genes


################### STAGER analys ##########################################
## -----------------------------------------------------------------------------
alpha <- 0.05
nGenes <- nrow(dge)
tableF <- topTable(fit2, number=nGenes, sort.by="none") #screening hypothesis
pScreen <- tableF$P.Value
names(pScreen) = rownames(tableF)

## -----------------------------------------------------------------------------
pConfirmation <- sapply(1:3,function(i) topTable(fit2, coef=i, number=nGenes, sort.by="none")$P.Value)
dimnames(pConfirmation) <- list(rownames(fit2),c("t1","t2","t1t2"))
stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=FALSE)

## -----------------------------------------------------------------------------
stageRObj1 <- stageWiseAdjustment(object=stageRObj, method="none", alpha=0.05)
stageRObj2 <- stageWiseAdjustment(object=stageRObj, method="holm", alpha=0.05)

## -----------------------------------------------------------------------------
head(getAdjustedPValues(stageRObj1, onlySignificantGenes=FALSE, order=FALSE))
head(getAdjustedPValues(stageRObj1, onlySignificantGenes=TRUE, order=TRUE))

head(getAdjustedPValues(stageRObj2, onlySignificantGenes=FALSE, order=FALSE))
head(getAdjustedPValues(stageRObj2, onlySignificantGenes=TRUE, order=TRUE))
## -----------------------------------------------------------------------------
res <- getResults(stageRObj1)
head(res)
colSums(res) #stage-wise analysis results

res <- getResults(stageRObj2)
head(res)
colSums(res) #stage-wise analysis results

## -----------------------------------------------------------------------------
DE_info

DE_info$celltype_de$de_output$table$`M-S`$L_NK_CD56._CD16. %>% str()
DE_info$celltype_de$de_output$data
DE_info$celltype_de$de_output$fit$L_NK_CD56._CD16.


glmQLFTest_output = glmQLFTest(DE_info$celltype_de$de_output$fit$L_NK_CD56._CD16., contrast = DE_info$celltype_de$de_output$args$contrast, poisson.bound=TRUE)
screen_tbl = topTags(glmQLFTest_output, n = glmQLFTest_output$coefficients %>% nrow(), adjust.method = "BH", sort.by = "PValue", p.value = 1)
pScreen = screen_tbl$table$PValue
names(pScreen) = rownames(screen_tbl$table)

glmQLFTest_output_confirmation_1 = glmQLFTest(DE_info$celltype_de$de_output$fit$L_NK_CD56._CD16., contrast = DE_info$celltype_de$de_output$args$contrast[,1], poisson.bound=TRUE) 
confirmation_tbl1 = topTags(glmQLFTest_output_confirmation_1, n = glmQLFTest_output_confirmation_1$coefficients %>% nrow(), adjust.method = "BH", sort.by = "PValue", p.value = 1) 
pConfirmation1 = confirmation_tbl1$table$PValue
names(pConfirmation1) = rownames(confirmation_tbl1$table)

glmQLFTest_output_confirmation_2 = glmQLFTest(DE_info$celltype_de$de_output$fit$L_NK_CD56._CD16., contrast = DE_info$celltype_de$de_output$args$contrast[,2], poisson.bound=TRUE) 
confirmation_tbl2 = topTags(glmQLFTest_output_confirmation_2, n = glmQLFTest_output_confirmation_2$coefficients %>% nrow(), adjust.method = "BH", sort.by = "PValue", p.value = 2) 
pConfirmation2 = confirmation_tbl2$table$PValue
names(pConfirmation2) = rownames(confirmation_tbl2$table)

glmQLFTest_output_confirmation_3 = glmQLFTest(DE_info$celltype_de$de_output$fit$L_NK_CD56._CD16., contrast = DE_info$celltype_de$de_output$args$contrast[,3], poisson.bound=TRUE) 
confirmation_tbl3 = topTags(glmQLFTest_output_confirmation_3, n = glmQLFTest_output_confirmation_3$coefficients %>% nrow(), adjust.method = "BH", sort.by = "PValue", p.value = 3) 
pConfirmation3 = confirmation_tbl3$table$PValue
names(pConfirmation3) = rownames(confirmation_tbl3$table)

pConfirmation = matrix(c(pConfirmation1[names(pScreen)], pConfirmation2[names(pScreen)], pConfirmation3[names(pScreen)]), ncol = 3)
rownames(pConfirmation) = names(pScreen)
colnames(pConfirmation) = c("M","S", "A")

stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=FALSE)
stageRObj <- stageWiseAdjustment(object=stageRObj, method="holm", alpha=0.05)

head(getAdjustedPValues(stageRObj, onlySignificantGenes=TRUE, order=TRUE))
DE_info$celltype_de$de_output_tidy %>% filter(contrast == "A-(S+M)/2" & cluster_id == "L_NK_CD56._CD16.") %>% arrange(p_adj)
DE_info$celltype_de$de_output_tidy %>% filter(contrast == "A-(S+M)/2" & cluster_id == "L_NK_CD56._CD16.") %>% arrange(p_adj) %>% filter(p_adj <= 0.05) %>% nrow()
getAdjustedPValues(stageRObj, onlySignificantGenes=TRUE, order=TRUE) %>% data.frame() %>% as_tibble() %>% filter(A <= 0.05) %>% nrow()

