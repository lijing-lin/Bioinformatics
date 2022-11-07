##### Part 1: Estimate the proportion of true null hypotheses from a vector of p-values.
# Use function propTrueNull from limma package

if(!require('limma')) {
  install.packages('limma')
  library('limma')
}

#  construct gene expression data from samples from 2 groups
n.genes <- 4000
n.samples <- 20

group <- factor(rep(c('A','B'), each=10))

data <- matrix(rnorm(n.genes*n.samples), ncol = n.samples) # gene X sample

data[c(1:200), c(1:10)]  = data[c(1:200),c(1:10)] + 2

# apply Wilcoxon-Mann-Whitney test on each row and abtain a list of 
p = rep(0, n.genes)
for (i in seq(1,n.genes)){
  r <- wilcox.test(data[i,] ~ group)
  p[i] <-r$p.value
}

1-200/4000

propTrueNull(p, method="lfdr")




##### Part 2: to identify not differentially expressed genes with equivalence test
# The following code were adapted from the link below 
# (Note: the old dataset is no longer available and therefore replaced)
# https://github.com/crazyhottommy/RNA-seq-analysis/blob/master/Not_differentially_expressed_genes_DESeq2.md

if(!require('DESeq2')) {
  install.packages('DESeq2')
  library('DESeq2')
}

library( DESeq2 )

# Read count data
countsName <- "https://bioconnector.github.io/workshops/data/airway_scaledcounts.csv"
download.file(countsName, destfile = "airway_scaledcounts.csv", method = "auto")

countData <- read.csv('airway_scaledcounts.csv', header = TRUE, sep = ",")
head(countData)

# Read experiment design data 
metaDataName <- "https://bioconnector.github.io/workshops/data/airway_metadata.csv"
download.file(metaDataName, destfile = "airway_metadata.csv", method = "auto")

metaData <- read.csv('airway_metadata.csv', header = TRUE, sep = ",")

metaData

 
# create a DESeqDataSet object from count data
dse <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ dex, tidy = TRUE)

# Perform a standard DESeq2 analysis
dse <- DESeq(dse)

# The log2 fold changes are found here
beta <- results(dse)$log2FoldChange
 
# The log fold change estimates all come with standard error (SE) information
# which we could find in both the 'rowData' and 'results')

a =  mcols(rowData(dse))
a$description

betaSE <- rowData(dse)$SE_dex_treated_vs_control

betaSE2 <- results(dse)$lfcSE  

# check if they are identical (should return TRUE)
all( betaSE == betaSE2, na.rm=TRUE )

# Internally, the Wald test is implemented as a simple two-sided
# z test of beta/betaSE. Two demonstrate this, we do the test
# manually and compare with the standard DESeq2 analysis output
 
pvalDE <- 2 * pnorm( abs( beta ), sd = betaSE, lower.tail=FALSE )

all( abs( pvalDE - results(dse)$pvalue ) < 1e-15, na.rm=TRUE )

# (returns TRUE)

# This was the test for DE, of course, i.e., small pvalDE means that
# the gene's expression change (the true value of beta) is not zero

# What we want is the opposite, namely find gene, for which abs(beta)
# is smaller than some threshold, theta
thr <- 0.4
# LL: to get a estimate of thr above, we could have a look at the percentage of
# True Null hypothesis using propTrueNull function introduced in PART 1
propTrueNull(pvalDE, method="lfdr")
h = hist(beta,breaks = 100)
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE)
# it looks 0.4 is a good estimate of threshold

# So, we do our two one-sided tests. For a one-sided z test, we
# simply use tail probabilities from the normal distribution.

# First, the test of H0_A: true_beta > thr
pA <- pnorm( beta, mean = thr, sd = betaSE, lower.tail=TRUE )

# Next, the test of H0_B: true_beta < -thr
pB <- pnorm( beta, mean = -thr, sd = betaSE, lower.tail=FALSE )

# The overall p value is the maximum, because we want to reject H0_A
# and H0_B simultaneously
pvalTOST <- pmax( pA, pB )


# Let's adjust our two p values with BH:
sigDE <- p.adjust( pvalDE, "BH" ) < .1
sigSmall <- p.adjust( pvalTOST, "BH" ) < .5


# And make an MA plot:
# DE genes plotted in red (col = 2 when sigDE==TRUE, sigSmall==FALSE) 
# non-DE genes plotted in black  (col = 1 when sigSmall==TRUE, sigDE==FALSE)
# genes neither DE nor non-DE plotted in green  (col = 0 when sigSmall==TRUE, sigDE==TRUE)


plot(
  rowMeans( counts(dse,normalized=TRUE) ), beta,
  log="x", pch= 0, cex= (1 + sigDE + 2*sigSmall)/2,
  col = 3 - sigDE - 2*sigSmall ,
  xlab = 'gene expression (normalised)', ylab = 'logFC')

# Plot is attached.

table(sigDE,sigSmall)


