
\name{IV}
\alias{IV}
\title{IV Analysis Using the Transformation Probit Model in Case-Cohort Studie}
\usage{
IV(ppt, data, m0=0, m1=1, PX=NULL, tol=0.005, SS=1, mu=0, bin=0.5, setting=NULL, level=NULL)
}
\arguments{
\item{data}{data.frame(Ts,D,xi,Z,M), with remaining columns as covariates}

\item{ppt}{A vector of time points of interest for the x-axis}

\item{setting}{Maximum tolerance of change during the iteration}

\item{steting}{Optional: "median", "mean", or "level"}

\item{level}{When "setting" is set to "level", specify a vector of the same length as the covariates}
}
\description{
Estimating the direct and indirect of the frailty model
}
\examples{
data(expdata)
data=expdata
ee=runif(30000,0,0.5);data[,1]=data[,1]+ee# separating the ties
colnames(data)[1:5] = c("Ts","D","xi","Z","M")
ppt=seq(0,max(data$Ts[data$D==1]),by=0.01) 

ans_IV=IV(ppt,data,m0=0.5,m1=2.5,SS=2.25^0.5,tol=0.001,PX=c(12,13,14),setting="median")
 plot(ppt,ans_IV$WIE,type="l")
 points(ppt,ans_IV$WIE+1.96*ans_IV$WIE.v)
 points(ppt,ans_IV$WIE-1.96*ans_IV$WIE.v)
 abline(h=0,col=2)
}
