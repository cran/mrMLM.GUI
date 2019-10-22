\name{ISIS}
\alias{ISIS}
\title{To perform GWAS with ISIS EM-BLASSO method}
\description{
Iterative Sure Independence Screening EM-Bayesian LASSO}
\usage{
ISIS(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,CLO)
}
\arguments{
	\item{gen}{genotype matrix.}
	\item{phe}{phenotype matrix.}
	\item{outATCG}{genotype for code 1.}
	\item{genRaw}{raw genotype.}
	\item{kk}{kinship matrix.}
	\item{psmatrix}{population structure matrix.}
	\item{svpal}{Critical P-value for selecting variable.}
	\item{svmlod}{Critical LOD score for significant QTN.}
	\item{Genformat}{Format for genotypic codes.}
      \item{CLO}{number of CPU.}
      	}
\author{
Zhang Ya-Wen, Li Pei, and Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(mrgen)
P1=data(mrphe)
G2=data(mrgenraw)
result=ISIS(mrgen,mrphe,outATCG=NULL,mrgenraw,kk=NULL,psmatrix=NULL,
0.01,3,1,CLO=1)
}