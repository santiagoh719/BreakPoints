\name{N_break_point}
\alias{N_break_point}
\title{N_break_point}
\description{
  Look for several breakpoints
}
\usage{
  N_break_point(serie, n_max = 1, n_period=10, seed='T',method='student',dstr='norm', auto_select = F, alpha = NULL)
}
\arguments{
  \item{serie}{numeric vector where the breakpoint is looked for}

  \item{n_max}{integer up to how many breakpoint should be test}

  \item{n_period}{an integer specifying the minimal length of a complete period to consider}

  \item{seed}{logical should the first breakpoint be creating using the inhomogeinity method (T, default) or equaly separate point should be use (FALSE) }

  \item{method}{which method should be use for breakpoint detection, supported: 'student' (default), 'mann-whitney','SNHT','buishand' and 'pettit'}

  \item{dstr}{character specifying wich distribution should be used for test simulations,if method is SNHT or buishand}

  \item{auto_select}{ lofical, should an automatic selection of how many breakpoints are be made, default FALSE}

  \item{alpha}{numeric, critical value to use for auto_select}

}



\value{
  \code{N_break_point} returns a list with the breakpoints index, it's p value and how many breakpoints are. If auto_select = F, a list with one list as specify for each n breakpoint tried
  \describe{
  \item{breaks}{index where the breakpoints are found}
  \item{p.value}{p value of the test}
  \item{n}{how many breakpoints are looked for}
  }
}

\details{
Compute homogeneity test for all possible breaks in the serie considering several breakpoints. NA values are allow.

}
\examples{

}
