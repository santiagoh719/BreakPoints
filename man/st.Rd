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

\code{pettit} Compute the Pettit Test for Homogeneity. NA values are allow.

}
\examples{
x <- c(rnorm(100,3,5),rnorm(200,5,5))
break1 <- stu(serie=x)
break2 <- man.whi(serie=x)
plot(x)
# plot a line where the breakpoint is found by the Student test
abline(v=break1$breaks, col='blue')
# plot a line where the breakpoint is found by the Mann-Whitney-Wilcoxon test
abline(v=break2$breaks, col='red')
}
\references{
- Ruxton, G. D., jul 2006. The unequal variance t-test is an underused alternative to Student’s t-test and
the Mann–Whitney U test. Behavioral Ecology 17 (4), 688–690.
URL: http://academic.oup.com/beheco/article/17/4/688/215960/The-unequal-variance-ttest-is-an-
underused

- Mann, H. B., Whitney, D. R., mar 1947. On a Test of Whether one of Two Random Variables is
Stochastically Larger than the Other. The Annals of Mathematical Statistics 18 (1), 50–60.
URL http://projecteuclid.org/euclid.aoms/1177730491

- Pettitt, A. N., 1979. A Non-Parametric Approach to the Change-Point Problem. Applied Statistics 28 (2),
126.
URL https://www.jstor.org/stable/10.2307/2346729?origin=crossref
}
