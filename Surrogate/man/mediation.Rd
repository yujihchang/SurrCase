
\name{mediation}
\alias{mediation}
\title{mediation Analysis Using the Transformation Probit Model in Case-Cohort Studie}
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
data(expdata)
data=expdata
ee=runif(30000,0,0.5);data[,1]=data[,1]+ee# separating the ties
colnames(data)[1:5] = c("Ts","D","xi","Z","M")
ppt=seq(0,max(data$Ts[data$D==1]),by=0.01) 

ans_med=mediation(ppt,data,int="FALSE",PX=c(12,13,14),setting="median")

 plot(ppt,ans_med$omega.IE,col=2,type="l",ylim=c(-0.03,0.03))
 points(ppt,ans_med$omega.IE+1.96*ans_med$omega.IE.v,col=2,type="l",ylim=c(-0.03,0.03))
 points(ppt,ans_med$omega.IE-1.96*ans_med$omega.IE.v,col=2,type="l",ylim=c(-0.03,0.03))
 points(ppt,ans_med$omega.DE-1.96*ans_med$omega.DE.v,col=1,type="l",ylim=c(-0.03,0.03))
 points(ppt,ans_med$omega.DE+1.96*ans_med$omega.DE.v,col=1,type="l",ylim=c(-0.03,0.03))
 points(ppt,ans_med$omega.DE,col=1,type="l",ylim=c(-0.03,0.03))
 abline(h=0,col=3)


}
