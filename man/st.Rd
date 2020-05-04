\name{Student}
\alias{stu}
\alias{man.whi}
\alias{pettit}
\title{Mann-Whitney-Wilcoxon Test, Student t-test and Pettit Test }
\description{
  Compute Rolling Mann-Whitney-Wilcoxon Test, Rolling Student t-test and Pettit test for homogeneity, NAs allow.
}
\usage{
  stu(serie,n_period=10)

  man.whi(serie,n_period=10)
  
  pettit(serie,n_period=10)
}
\arguments{
  \item{serie}{numeric vector where the breakpoint is looked for}

  \item{n_period}{an integer specifying the minimal length of a complete period to consider}

}



\value{
  \code{pettit}, \code{man.whi} and \code{stu} returns a list with the breakpoint index and its p value
  \describe{
  \item{breaks}{index where the breakpoint is found}
  \item{p.value}{p value of the test}
  }
}

\details{
\code{man.whi} compute Mann-Whitney-Wilcoxon Test and \code{stu} the Student t-test for all possible breaks in the serie and get the most significant break.In both test NA values are allow.

\code{pettit} Compute de Pettit Test for Homogeneity. NA values are allow. Pettitt 1979 in Applied statistics.

}
\examples{
x <- c(rnorm(100,3,5),rnorm(200,5,5))
break1 <- stu(serie=x)
break2 <- man.whi(serie=x)
plot(x)
abline(v=break1$breaks, col='blue') # to plot a line where the breakpoint is found by the Student test
abline(v=break2$breaks, col='red') # to plot a line where the breakpoint is found by the Mann-Whitney-Wilcoxon test
}
