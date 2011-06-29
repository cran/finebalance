\name{function-internal}
\Rdversion{1.1}
\alias{minchisq}
\alias{mintotal}
\alias{minorder}
\alias{bound}
\alias{minchisq.force}
\alias{mintotal.force}
\alias{minorder.force}
\alias{bound.force}
\alias{smahal}
\title{
Internal finebalance functions .
}
\description{
Internal functions called by finebalance.
}
\usage{
minchisq(n.treated, n.control, abs.diff = NULL, perc.diff = NULL)
mintotal(n.treated, n.control, abs.diff = NULL, perc.diff = NULL)
minorder(n.treated, n.control, abs.diff = NULL, perc.diff = NULL)
}
\arguments{
  \item{n.treated}{
a vector of the number of treated subjects in each stratum
}
  \item{n.control}{
a vector of the number of control subjects in each stratum
}
  \item{abs.diff}{
maximum absolute difference between the number of treated and number of matched controls in each stratum
}
  \item{perc.diff}{
The maximum ratio of the difference between number of matched controls and number of treated units to the number of treated units in a stratum allowed
}
}
\details{
 These are not to be called directly by the users.
}
\value{
  \item{n.matched.up}{upper bound for the number of matched controls in each stratum}
  \item{n.matched.low }{lower bound for the number of matched controls in each stratum}
}
\references{
Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes
}
\author{
Dan Yang
}
\keyword{internal}
