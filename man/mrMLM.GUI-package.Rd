\name{mrMLM.GUI-package}
\alias{mrMLM.GUI-package}
\alias{mrMLM.GUI}
\docType{package}
\title{
Multi-Locus Random-SNP-Effect Mixed Linear Model for Multi-Locus GWAS and Multi-QTL Mapping
}
\description{
Conduct multi-locus genome-wide association study under the framework of random-SNP-effect mixed linear model (mrMLM). First, each marker on the genome is scanned. Bonferroni correction is replaced by a less stringent selection criterion for significant test. Then, all the markers that are potentially associated with the trait are included in a multi-locus model, their effects are estimated by empirical Bayes and true QTN are identified by likelihood ratio test.
}
\details{

\tabular{ll}{
Package: \tab mrMLM.GUI\cr
Type: \tab Package\cr
Version: \tab 3.1\cr
Date: \tab 2018-4-25\cr
Depends: \tab shiny,MASS,data.table,doParallel,foreach\cr
Imports: \tab methods,openxlsx,stringr,qqman,ggplot2,lars,ncvreg,coin\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
Users can use library(mrMLM.GUI) to start the GUI and use 'mrMLM.GUI()' to restart the programe.
}
\author{
Zhang Ya-Wen, Li Pei, Ren Wen-Long, Ni Yuan-Li, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
Improving power and accuracy of genome-wide association studies via a multi-locus mixed linear model methodology.
Wang Shi-Bo,Feng Jian-Ying,Ren Wen-Long,Huang Bo,Zhou Ling,Wen Yang-Jun,Zhang Jin,Jim M. Dunwell,Xu Shizhong*,Zhang Yuan-Ming*
}

\keyword{
 R
mrMLM.GUI
package
}

\examples{
\dontrun{ mrMLM.GUI() }
}
