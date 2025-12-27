#Package check & install
packages <- c("BiocManager", "edgeR", "limma", "affy")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "BiocManager") {
      install.packages("BiocManager")
    } else {
      BiocManager::install(pkg, ask = FALSE)
    }
  }
}

#Load libraries
library(edgeR)
library(limma)
library(affy)


#Load count matrix
counts <- read.csv("samples.csv", row.names=1)
head(counts)

#Create DGE object
# Stores count data and sample information
data1 <- DGEList(counts)
data1 <- calcNormFactors(data1)

#Filter low-expression genes
#Removes genes with very low counts across samples
cutoff <- 5
drop <- which(apply(cpm(data1) , 1, max) < cutoff)
d <- data1[-drop,]
dim(d)


snames <- colnames(counts)
group <- interaction(substr(snames , 1 , nchar(snames)-2))
plotMDS(d, col = as.numeric(group))
mm <- model.matrix(~0 + group)
y <- voom(d , mm , plot = T)
fit <- lmFit(y , mm)
mm
contr <- makeContrasts(groupLum.A. - groupTNBC. , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit , contr)
tmp <- eBayes(tmp)

#Differential expression results
#Extract top genes ranked by p-value
top.table <- topTable(tmp , sort.by = "P" , n = Inf)
head(top.table , 20)

contr <- makeContrasts(groupLum.A. - groupLum.B. , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit , contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp , sort.by = "P" , n = Inf)
head(top.table , 20)