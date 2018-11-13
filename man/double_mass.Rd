\name{double_mass}
\alias{double_mass}

\title{Double Mass Curve }
\description{
    Compute Double Mass Curve comparing your serie with simulated serie.
}
\usage{
  double_mass(serie,ploting=T,date_axis=NULL,simulations = 10000, alpha = 0.5)

}
\arguments{
  \item{serie}{numeric vector where the breakpoint is looked for}

  \item{ploting}{ logical, should the curve be plot (default TRUE), if FALSE the output is a data.frame with the values for the curves}

  \item{date_axis}{ vector of type Date or POSIX with the date of the serie in same order. If given will be used for x-axis in the plot }

  \item{simulations}{an integer specifying how many montecarlo simulations to perform}

  \item{alpha}{numeric, critical value to select upper and lower boundaries in simulations}
}



\value{
  \code{double_mass} returns either a ggplot objet or a list with the commulative sum and the boundaries of comulative sum of simulated data
  \describe{
    \item{real}{ vector of the commulative sum of the serie}
    \item{lower_norm}{ vector of lower boundary of normal distribution simulations comulative sum}
    \item{upper_norm}{ vector of upper boundary of normal distribution simulations comulative sum}
    \item{lower_gamma}{ vector of lower boundary of gamma distribution simulations comulative sum}
    \item{upper_gamma}{ vector of upper boundary of gamma distribution simulations comulative sum}
  }
}

\details{
  

}
\examples{

}
