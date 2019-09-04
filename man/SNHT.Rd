\name{SNHT}
\alias{SNHT}
\alias{Buishand_R}
\title{Buishand Range Test and Standard Normal Homogeneity Test}
\description{
  Compute Buishand Range Test or Standard Normal Homogeneity Test for a serie, NAs allow in both Test
}
\usage{
  SNHT(serie,n_period=10,dstr='norm',simulations = 1000)

  Buishand_R(serie,n_period=10,dstr='norm',simulations = 1000)

}
\arguments{
  \item{serie}{numeric vector where the breakpoint is looked for}

  \item{n_period}{an integer specifying the minimal length of a complete period to consider}

  \item{dstr}{character specifying which distribution should be used for test simulations, 'norm' (default; normal distribution), 'gamma', and 'self' (will compute bootstrap)}

  \item{simulations}{an integer specifying how many Monte Carlo simulations to perform, default is 1000.}

}



\value{
  \code{SNHT} and \code{Buishand_R} returns a list with the breakpoint index and it's p value
  \describe{
    \item{breaks}{index where the breakpoint is found}
    \item{p.value}{p value of the test}
  }
}

\details{
  \code{SNHT} compute Standard Normal Homogeneity Test where NA values are allow.

  \code{Buishand_R} Compute Buishand Range Test for Homogeneity where NA values are allow. Publish by Buishand in 1982 in Journal of Hydrology.

}
\examples{
  x <- c(runif(100,3,5),runif(200,5,5))
  break1 <- SNHT(serie=x)
  break2 <- Buishand_R(serie=x)
  plot(x)
  abline(v=break1$breaks, col='blue') # to plot a line where the breakpoint is found by the SNHT test
  abline(v=break2$breaks, col='red') # to plot a line where the breakpoint is found by the Buishand Range test
}
