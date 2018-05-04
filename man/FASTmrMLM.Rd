\name{FASTmrMLM}
\alias{FASTmrMLM}
\title{To perform GWAS with FASTmrMLM method}
\description{
FAST multi-locus random-SNP-effect Mixed Linear Model}
\usage{
FASTmrMLM(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svrad,svmlod,Genformat,CLO)
}
\arguments{
	\item{gen}{genotype matrix.}
	\item{phe}{phenotype matrix.}
	\item{outATCG}{genotype for code 1.}
	\item{genRaw}{raw genotype.}
	\item{kk}{kinship matrix.}
	\item{psmatrix}{population structure matrix.}
	\item{svpal}{Critical P-value for selecting variable.}
      \item{svrad}{Search Radius in search of potentially associated QTN.}
      \item{svmlod}{Critical LOD score for significant QTN.}
	\item{Genformat}{Format for genotypic codes.}
      \item{CLO}{number of CPU.}
	}
\author{
Zhang Ya-Wen, Li Pei, Ren Wen-Long, Ni Yuan-Li, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(mrgen)
P1=data(mrphe)
G2=data(mrgenraw)
FASTmrMLM(mrgen,mrphe,outATCG=NULL,mrgenraw,kk=NULL,psmatrix=NULL,
0.01,20,3,1,CLO=1)
}


