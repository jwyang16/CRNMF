# CRNMF: Constrained Robust Non-negative Matrix Factorization

CRNMF is a software tool designed for dropout imputation and dimension reduction of scRNA-seq data. CRNMF imputes scRNA-seq datasets using a constrained non-negative matrix factorization model.


## CRNMF Installation

CRNMF can be installed via Github.
To install the latest version of CRNMF package via Github, run following commands in R:

	if (!require("devtools"))
	   install.packages("devtools")
	devtools::install_github("jwyang16/CRNMF")

## Getting Started

Users can load CRNMF as follows:
	
	library(CRNMF)

## Example
Run CRNMF to impute and cluster scRNA-seq datasets

	data(deng)
	data <- deng$deng.data
	label <- deng$deng.label

	mat <- CRNMF(data,maxitr=30,lambda=1.5,r=4)
	W <- mat$W
	H <- mat$H
	S <- mat$S

	H2 <- H/rep(colSums(H),each=nrow(H))
	set.seed(2019)
	I <- kmeans(t(H2),6,iter.max=250,nstart=150)
	evalt(label,I$cluster)

