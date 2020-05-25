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

  \code{Buishand_R} Compute Buishand Range Test for Homogeneity where NA values are allow.

}

\references{
- Alexandersson, H., jan 1986. A homogeneity test applied to precipitation data. Journal of Climatology
6 (6), 661–675.
URL http://doi.wiley.com/10.1002/joc.3370060607

- Buishand, T., aug 1982. Some methods for testing the homogeneity of rainfall records. Journal of Hydrology
58 (1-2), 11–27.
URL https://www.sciencedirect.com/science/article/pii/002216948290066Xhttps://linkinghub.elsevier.
com/retrieve/pii/002216948290066X
}
