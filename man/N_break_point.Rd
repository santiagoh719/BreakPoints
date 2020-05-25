\name{N_break_point}
\alias{N_break_point}
\title{N_break_point}
\description{
  Look for several breakpoints
}
\usage{
  N_break_point(serie, n_max = 1, n_period=10,
                          seed=FALSE, auto_select = FALSE,
                          alpha = NULL,method='SNHT',dstr='norm')
}
\arguments{
  \item{serie}{numeric vector where the breakpoint is looked for}

  \item{n_max}{integer up to how many breakpoint should be test}

  \item{n_period}{an integer specifying the minimal length of a complete period to consider}

  \item{seed}{ if FALSE (default) the first breakpoints to inicialize the iteration are calculated by splitting the serie in equal separated parts. If seed is given must be a list() of length n_max where it is specified the first breakpoints to take for each iteration. }

  \item{method}{which method should be used for breakpoint detection, supported: 'student', 'mann-whitney','SNHT' (default),'buishand' and 'pettit'}

  \item{dstr}{character specifying which distribution should be used for test simulations, only used if method is SNHT or buishand. Posible distributions 'norm' (default, normal dist), 'gamma', and 'self' (bootstrap)}

  \item{auto_select}{ logical, should an automatic selection of how many breakpoints are be made, default FALSE}

  \item{alpha}{numeric, critical p value to use for auto_select}

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

\references{
Hurtado, S. I., Zaninelli, P. G., & Agosta, E. A. (2020). A multi-breakpoint methodology to detect changes in climatic time series. An application to wet season precipitation in subtropical Argentina. Atmospheric Research, 104955.
}
