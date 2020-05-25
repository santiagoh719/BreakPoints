\name{yamamoto}
\alias{yamamoto}
\alias{sum_yam}

\title{Yamamoto}
\description{
  Compute Yamamoto Test for climate jumps in a serie.
}
\usage{
  yamamoto(serie, alpha = 0.1, n_period = 10)

}
\arguments{
  \item{serie}{numeric vector where the breakpoint is looked for}

  \item{n_period}{an integer specifying the length of the window to use, can not bet odd}

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
\references{
  Yamamoto, R., Iwashima, T., Kazadi, S. N., & Hoshiai, M. (1985). Climatic jump: a hypothesis in climate diagnosis. Journal of the Meteorological Society of Japan. Ser. II, 63(6), 1157-1160.
}
