\name{mrMLM.GUI-package}
\alias{mrMLM.GUI-package}
\alias{mrMLM.GUI}
\docType{package}
\title{
Multi-Locus Random-SNP-Effect Mixed Linear Model Tools for Genome-Wide Association Study with Graphical User Interface
}
\description{
Description: Conduct multi-locus genome-wide association study under the framework of multi-locus random-SNP-effect mixed linear model (mrMLM). First, each marker on the genome is scanned. Bonferroni correction is replaced by a less stringent selection criterion for significant test. Then, all the markers that are potentially associated with the trait are included in a multi-locus genetic model, their effects are estimated by empirical Bayes and all the nonzero effects were further identified by likelihood ratio test for ture QTL. 
}
\details{
\tabular{ll}{
Package: \tab mrMLM.GUI\cr
Type: \tab Package\cr
Version: \tab 4.0.1\cr
Date: \tab 2020-9-29\cr
Depends: \tab shiny,lars\cr
Imports: \tab methods,foreach,ncvreg,coin,sampling,data.table,doParallel,shinyjs,bigmemory,mrMLM\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
Users can use library(mrMLM.GUI) to start the GUI and use 'mrMLM.GUI()' to restart the programe.
}
\author{
Zhang Ya-Wen, Li Pei, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
1. Zhang YM, Mao Y, Xie C, Smith H, Luo L, Xu S*. Genetics 2005,169:2267-2275.
2. Wang SB, Feng JY, Ren WL, Huang B, Zhou L, Wen YJ, et al. Sci Rep 2016,6:19444.
3. Tamba CL, Ni YL, Zhang YM*. PLoS Comput Biol 2017,13(1):e1005357.
4. Zhang J, Feng JY, Ni YL, Wen YJ, Niu Y, Tamba CL, et al. Heredity 2018,118(6):517-524.
5. Ren WL, Wen YJ, Dunwell JM, Zhang YM*. Heredity 2018,120(3): 208-218.
6. Wen YJ, Zhang H, Ni YL, Huang B, Zhang J, Feng JY, et al. Brief Bioinform 2018,19(4):700-712.
7. Tamba CL, Zhang YM. bioRxiv,preprint first posted online Jun. 7, 2018, doi:https://doi.org/10.1101/341784.
8. Zhang YW, Tamba CL, Wen YJ, Li P, Ren WL, Ni YL, et al. Genomics, Proteomics & Bioinformatics, Accept.
}
\examples{
\dontrun{ mrMLM.GUI() }
}