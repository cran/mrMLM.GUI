\name{pKWmEB}
\alias{pKWmEB}
\title{To perform GWAS with pKWmEB method}
\description{
Kruskal-Wallis test with empirical Bayes under polygenic background control}
\usage{
pKWmEB(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,CLO)
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
result=pKWmEB(mrgen,mrphe,outATCG=NULL,mrgenraw,kk=NULL,psmatrix=NULL,
0.05,3,1,CLO=1)
}