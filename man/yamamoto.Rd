\name{yamamoto}
\alias{yamamoto}

\title{Yamamoto}
\description{
  Compute Yamamoto Test for a serie.
}
\usage{
  yamamoto(serie, alpha = 0.1, n_period = 10)

}
\arguments{
  \item{serie}{numeric vector where the breakpoint is looked for}

  \item{n_period}{an integer specifying the minimal length of a complete period to consider}

  \item{alpha}{numeric, p value to use}

}



\value{
  \code{yamamoto} returns a list with the breakpoints indexes and the amount
  \describe{
    \item{breaks}{vector of indexes where the breakpoint is found}
    \item{n}{Amount of breakpoints}
  }
}
\details{
  \code{yamamoto} compute the Yamamoto Test.
}
