\name{permPATH2HTML}
\alias{permPATH2HTML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 This is a function for creating an HTML file
}
\description{
  The function creates an HTML file.
}
\usage{
permPATH2HTML(dat, dir, fname, title=NULL, bgcolor="#BBBBEE")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    A data frame.
  }
  \item{dir}{
    Directory in which to store the file.
  }
  \item{fname}{
    File name.
  }
  \item{title}{
    The title of the html file.
  }
  \item{bgcolor}{
    Color for the html background.
  }
}

\examples{
## Generate toy phenotype and gene expression data sets
## This example consists of 40 genes grouped into 5 pathways and 100 patients
## grp is a binary trait (e.g., case vs control)
## bp is a continuous trait (e.g., blood pressure)
set.seed(1234)
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
res = perm.path(expdat, y=pdat[["grp"]], local.test="wilcoxon", global.test="maxmean", 
B=100, gset=gset, min.num=2, max.num=50, sort="score")

# create an html file
#epermPATH2HTML(rstab, dir="/dir/", fname="tophits")
}
