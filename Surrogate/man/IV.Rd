
\name{IV}
\alias{IV}
\title{IV Analysis Using the Transformation Probit Model in Case-Cohort Studie}
\usage{
IV(ppt, data, m0=0, m1=1, PX=NULL, tol=0.005, SS=1, mu=0, bin=0.5, setting=NULL, level=NULL)
}
\arguments{
\item{data}{data.frame(X1,X2,D,Z)}

\item{int_theta}{initial value of theta for iteration, nonnegative values vector of length 2}

\item{tol}{maximum tolerance of change during the iteration}

\item{step}{maximum  number of the iteration}

\item{interpolation}{time can be vector or scalar}
}
\description{
Estimating the direct and indirect of the frailty model
}
\examples{
ee=runif(30000,0,0.5);data[,1]=data[,1]+ee# separating the ties
colnames(data)[1:5] = c("Ts","D","xi","Z","M")
ppt=seq(0,max(data$Ts[data$D==1]),by=0.01) 

ans_IV=IV(ppt,data,m0=0.5,m1=2.5,SS=2.25^0.5,tol=0.001,PX=c(12,13,14),setting="median")
 plot(ppt,ans_med$WIE,type="l)
 points(ppt,ans_med$WIE+1.96*ans_med$WIE.v)
 points(ppt,ans_med$WIE-1.96*ans_med$WIE.v)
 abline(h=0,col=2)
}
