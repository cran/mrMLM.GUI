\name{mrMLM.GUI-package}
\alias{mrMLM.GUI-package}
\alias{mrMLM.GUI}
\docType{package}
\title{
Multi-Locus Random-SNP-Effect Mixed Linear Model for Multi-Locus GWAS and Multi-QTL Mapping
}
\description{
Conduct multi-locus genome-wide association study under the framework of random-SNP-effect mixed linear model (mrMLM). First, each marker on the genome is scanned. Bonferroni correction is replaced by a less stringent selection criterion for significant test. Then, all the markers that are potentially associated with the trait are included in a multi-locus model, their effects are estimated by empirical Bayes and true QTNs are identified by likelihood ratio test.
}
\details{

\tabular{ll}{
Package: \tab mrMLM.GUI\cr
Type: \tab Package\cr
Version: \tab 3.2\cr
Date: \tab 2018-8-4\cr
Depends: \tab shiny,MASS,data.table,doParallel,foreach\cr
Imports: \tab methods,openxlsx,stringr,qqman,ggplot2,lars,ncvreg,coin\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
Users can use library(mrMLM.GUI) to start the GUI and use 'mrMLM.GUI()' to restart the programe.
}
\author{
Zhang Ya-Wen, Li Pei, Ren Wen-Long, Ni Yuan-Li, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
Methodological implementation of mixed linear models in multi-locus genome-wide association studies.
Wen Yang-Jun, Zhang Hanwen, Ni Yuan-Li, Huang Bo, Zhang Jin, Feng Jian-Ying, Wang Shi-Bo, Dunwell Jim M., Zhang Yuan-Ming, Wu Rongling.
}

\keyword{
 R
mrMLM.GUI
package
}

\examples{
\donttest{ mrMLM.GUI() }
}
