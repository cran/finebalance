\name{finebalance}
\alias{finebalance}
\title{
Approximate fine balance when exact fine balance is not achievable
}
\description{
Uses the algorithm of minimum cost network flow.
}
\usage{
finebalance(distance.matrix, z, f, type = c("total", "order", "chisq"), abs.diff = NULL, perc.diff = NULL, force = NULL)
}
\arguments{
  \item{distance.matrix}{
nt by nc matrix, the (i,j)the entry is the distance of the ith treated and jth control
}
  \item{z}{
a vector of length n=nt+nc, 0,1 valued, 1 being treated, 0 control
}
  \item{f}{
a vector of length nt+nc, factor to be balanced
}
  \item{type}{
Specifies the criterion of the imbalance; see the paper for details
}
  \item{abs.diff}{
The maximum difference between number of matched controls and number of treated units in a stratum allowed 
}
  \item{perc.diff}{
The maximum ratio of the difference between number of matched controls and number of treated units to the number of treated units in a stratum allowed 
}
  \item{force}{
vector of length nt+nc, 0, 1 valued, 1 means the individual must be matched, 0 means the individual can be either matched or not.  
If force is provided, abs.diff and perc.diff will be ignored
}
}

\details{
"total" minimizes the total distance while minimizing the sum of the absolute difference between number of matched controls and number of treated units.
"order" minimizes the maximum deviation from number of treated units.
"chisq" minimizes (n.matched.control-n.treated)^2/n.treated
}
\value{
\item{matchpair}{The matchpair[i]th person is matched to the ith treated}
\item{matchpairnames}{Make use of the colnames and rownames of the distance matrix, and matchpairnames[i] is matched to the ith treated}
\item{total.distance }{Total distance after match}
\item{matchtable}{Returns a table description of the number of treated, number of control, and number of matched control in each stratum}
}
\references{
Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes
}
\author{
Dan Yang
}
\examples{
treatment=c(rep(1,21),rep(0,26));
age=c(38,44,39,33,35,39,27,43,39,43,41,36,35,37,39,34,35,53,38,37,38,48,63,44,40,50,52,56,47,38,34,42,36,41,41,31,56,51,36,44,35,34,39,45,42,30,35);
race=c(0,0,0,1,rep(0,5),1,rep(0,11),1,rep(0,12),rep(1,4),rep(0,9));
smoker=c(rep(0,2),rep(1,4),0,1,1,0,1,0,0,0,1,0,1,0,1,0,1,0,0,1,rep(0,8),1,0,1,1,1,0,1,0,0,1,1,0,0,0,1);
Xmat=cbind(age,race,smoker);
dpc=c(1.77,1.02,1.44,.65,2.08,.61,2.86,4.19,4.88,1.08,2.03,2.81,.94,1.43,1.25,2.97,1.01,2.07,1.15,1.07,1.63,1.08,1.09,1.1,1.1,.93,1.11,.98,2.2,.88,1.55,.55,1.04,1.66,1.49,1.36,1.02,.99,.65,.42,2.33,.97,.62,1.02,1.78,.95,1.59);
agesq=age^2;
age.race=age*race;
age.smoker=age*smoker;
race.smoker=race*smoker;
model3=glm(treatment~age+agesq+race+smoker+age.race+age.smoker+race.smoker,family=binomial);
propscore.model3=predict(model3,type="response")
distance.matrix=smahal(treatment,Xmat);
race.and.smoker=rep(1,length(race))
race.and.smoker[race==0&smoker==0]=1
race.and.smoker[race==1&smoker==0]=2
race.and.smoker[race==0&smoker==1]=3
race.and.smoker[race==1&smoker==1]=4
finebalance(distance.matrix, treatment, as.factor(race.and.smoker), type="total")
}
\keyword{nonparametric}
\keyword{network optimization}
\keyword{multivariate matching}
