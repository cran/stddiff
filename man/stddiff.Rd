\name{stddiff}
\alias{stddiff}
\alias{stddiff.numeric}
\alias{stddiff.binary}
\alias{stddiff.category}
\title{Calculate the Standardized Difference for Numeric, Binary and Category Variables}
\description{Contains three main functions including stddiff.numeric(), stddiff.binary() and stddiff.category(). These are used to calculate the standardized difference between two groups. It is especially used to evaluate the balance between two groups before and after propensity score matching.
}
\usage{
stddiff.numeric(data,gcol,vcol)
stddiff.binary(data,gcol,vcol)
stddiff.category(data,gcol,vcol)
}
\arguments{
  \item{data}{a dataframe}
  \item{gcol}{a column number of group variable in data, 0 for control group, 1 for treatment group}
  \item{vcol}{one or more column numbers of different types variables in data}
}
\details{
stddiff.numeric() is used for numeric variables.
stddiff.binary() is used for binomial variables.
stddiff.category() is used for categorical variables.
stddiff should be less than 0.2.
}
\value{
for stddiff.numeric function
\item{mean.c}{the mean of control group}
\item{sd.c}{the standard deviation of control group}
\item{mean.t}{the mean of treatment group}
\item{sd.t}{the standard deviation of treatment group}
\item{missing.c}{the counts of missing value of control group}
\item{missing.t}{the counts of missing value of treatment group stddiff: the standardized difference between two groups}
\item{stddiff.l}{the lower limit of the 95 percentage confidence interval of standardized difference}
\item{stddiff.u}{the upper limit of the 95 percentage confidence interval of standardized difference}

for stddiff.binary function:
\item{p.c}{the proportion of last level in the control group}
\item{p.t}{the proportion of last level in the treatment group}
\item{missing.c}{the counts of missing value of control group}
\item{missing.t}{the counts of missing value of treatment group}
\item{stddiff}{the standardized difference between two groups}
\item{stddiff.l}{the lower limit of the 95 percentage confidence interval of standardized difference}
\item{stddiff.u}{the upper limit of the 95 percentage confidence interval of standardized difference}

for stddiff.category function:
\item{p.c}{the proportion of each level in the control group}
\item{p.t}{the proportion of each level in the treatment group}
\item{missing.c}{the counts of missing value of control group}
\item{missing.t}{the counts of missing value of treatment group}
\item{stddiff}{the standardized difference between two groups}
\item{stddiff.l}{the lower limit of the 95 percentage confidence interval of standardized difference}
\item{stddiff.u}{the upper limit of the 95 percentage confidence interval of standardized difference}
}
\references{
Yang DS, Dalton JE. A Unified Approach to Measuring the Effect Size Between Two Groups Using SAS. SAS Global Forum 2012. paper 335
}
\author{
Zhicheng Du<dgdzc@hotmail.com>, Yuantao Hao<haoyt@mail.sysu.edu.cn>
}
\note{
nothing
}
\seealso{
nothing
}
\examples{
#set.seed(2016)
#treat=round(abs(rnorm(100)+1)*10,0) %% 2
#numeric=round(abs(rnorm(100)+1)*10,0)
#binary=round(abs(rnorm(100)+1)*10,0) %% 2
#category=round(abs(rnorm(100)+1)*10,0) %% 3
#data=data.frame(treat,numeric,binary,category)
#stddiff.numeric(data=data,gcol=1,vcol=c(2,2))
#stddiff.binary(data=data,gcol=1,vcol=c(3,3))
#stddiff.category(data=data,gcol=1,vcol=c(4,4))
}
\keyword{standardized difference,propensity score matching}
