\name{perm.path}
\alias{perm.path}

\title{
 Perform Permutation Based Pathway Analysis
}
\description{
 This is the package main function.
}
\usage{
perm.path(expr, y, local.test, global.test="wilcoxon", B, gset, min.num=2, max.num, 
imputeval=NULL, transfun=function(x){x}, sort="pval", anno=NULL)
}

\arguments{
  \item{expr}{
    An \eqn{K \times n} matrix of gene expression data, where \eqn{K} is the number of genes 
    and \eqn{n} is the number of samples.  
}
  \item{y}{
    An outcome vector of length \eqn{n}. 
}
  \item{local.test}{
    Local test statistic of each gene. Current possible choices are \eqn{t-test}, \eqn{Wilcoxon}
    test, \eqn{Pearson}, \eqn{Spearman} and \eqn{JT} test.
}
  \item{global.test}{
    Global test statictic, used to compute the score. Current possible choices are \eqn{mean},
    \eqn{meanabs} (mean of absolute values) and \eqn{maxmean}.
}
  \item{B}{
    specifies the number of random permutations to be performed.
}
  \item{gset}{
    A list of pathways. Each element is a vector of gene names. 
The list element names are the pathway names.
}
  \item{min.num}{
    Specifies the minimum number of genes that a pathway should have. 
Pathways with smaller number of genes will be excluded.
}
  \item{max.num}{
    Specifies the maximum number of genes that a pathway should have. 
Pathways with larger number of genes will be excluded.
}
  \item{imputeval}{
    The gene expression value to be imputed in case of missing values. 
The default choice is \eqn{NULL} in which case no imputation is done.
}
  \item{transfun}{
    Specifies transformation of the gene expression data. 
The default option is untransformed gene expression data.
}
  \item{sort}{
    Specifies sorting of the results. 
If \eqn{sort="pval"} sorting is done in order of increasing \eqn{p-values}. 
    If \eqn{sort="score"} sorting is done in order of decreasing \eqn{scores}.
}
  \item{anno}{
    If \eqn{TRUE} the output contains annotation of each pathway.
}
}

\value{
  This function returns a list consisting of the following elements:
  \item{res}{Data frame consisting of the pathway names (Pathway), the genes 
involved in each pathway (Genes), the number of genes in each pathway (Size), 
the score for each pathway (Score), the permutation raw p-value (pval), 
the FWER-adjusted permutation p-value (pfwer), the FDR-adjusted permutation p-value, 
the Bonferroni-adjusted permutation p-value (bonferroni)}
  \item{stats}{The individual test statistic for each gene}
  \item{scores}{A matrix of scores. The matrix is of dimension \eqn{(B+1)\times K}, 
where K is the number of pathways. The first column contains the unpermuted scores, 
the remaining \eqn{B} columns contain the scores computed after each permutation.}
}
\references{
   B. Efron, R. Tibshirani (2007) 
On Testing the Significance of Sets of Genes. \emph{The Annals of Applied Statistics}. 
\bold{Vol. 1}, No 1, 107--129.

   A. Subramanian, P. Tamayo, V. K. Mootha, S. Mukherjee, B. L. Ebert, M. A. Gillette, 
   A. Paulovich, S. L. Pomeroy, T. R. Golub, E. S. Lander and J. P. Mesirov (2005), 
   Gene Set Enrichment Analysis: A knowledge-based Approach for Interpreting 
   Genome-Wide Expression Profiles. \emph{Proc. Natl. Acad. Sci. USA}. \bold{Vol. 102}, No 43, 15545--15550.
}

\examples{
set.seed(1234)

## Generate toy phenotype and gene expression data sets
## This example consists of 40 genes grouped into 5 pathways and 100 patients
## grp is a binary trait (e.g., case vs control)
## bp is a continuous trait (e.g., blood pressure)

n = 100
K = 40
grp = rep(1:0,each=n/2)
bp = rnorm(n)

pdat = data.frame(grp, bp)
rm(grp, bp)
expdat = matrix(rnorm(K*n),K,n)

## Assign marker names g1,...,gK to the expression data set and
## patient ids id1,...,idn to the expression and phenotype data
gnames = paste("g",1:K,sep="")
rownames(expdat) = gnames
patid = paste("id",1:n,sep="")
rownames(pdat) = patid
colnames(expdat) = patid

#Group the K genes into M pathways of sizes n1,...,nM
M = 5
p = runif(M)
p = p/sum(p)
nM = rmultinom(1, size=K, prob=p)
gset = lapply(nM, function(x){gnames[sample(x)]})
names(gset) = paste("pathway",1:M,sep="")

## Carry out permutation analysis with grp as the outcome
## using the two-sample Wilcoxon with B=100 random permutations
perm.path(expdat, y=pdat[["grp"]], local.test="wilcoxon", global.test="maxmean", B=100, 
gset=gset, min.num=2, max.num=50, sort="score")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
