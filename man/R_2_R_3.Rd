% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_way_crossed_ANOVA.R
\name{R_2_R_3}
\alias{R_2_R_3}
\title{Function to evaluate the observed, critical and P-values of the test statistics R_2 and R_3}
\usage{
R_2_R_3(error, a, b, cellsize, cellmean, cellvariance)
}
\arguments{
\item{error}{Numeric, the level of significance/type-I error}

\item{a}{Numeric, the number of levels of the factor A}

\item{b}{Numeric, the number of levels of the factor B}

\item{cellsize}{matrix having a rows and b columns, and (i,j)-th element is the number of observations in the (i,j)-th cell}

\item{cellmean}{matrix having a rows and b columns, and (i,j)-th element is the sample mean of the observations in the (i,j)-th cell}

\item{cellvariance}{matrix having a rows and b columns, and (i,j)-th element is the sample variance of the observations in the (i,j)-th cell}
}
\value{
numeric vector consisting of test statistic values, critical values and P-values of the test R_2 and R_3. Also, the decisions of the tests R_2 and R_3.
}
\description{
Function to evaluate the observed, critical and P-values of the test statistics R_2 and R_3
}
\examples{
R_2_R_3(0.05,3,2, rbind(c(35,35),c(56,56),c(51,51)), rbind(c(52.7429,99.9143),c(58.8750,126.7679),c(95.8039,261.5490)), rbind(c(67.2555,274.0807),c(222.8023,955.0542),c(815.3608,8597.5325)))
}
