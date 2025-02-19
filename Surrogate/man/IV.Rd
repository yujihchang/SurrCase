
\name{IV}
\alias{IV}
\title{IV Analysis Using the Transformation Probit Model in Case-Cohort Studie}
\usage{
IV(ppt, data, m0=0, m1=1, PX=NULL, tol=0.005, SS=1, mu=0, bin=0.5, setting=NULL, level=NULL)
}
\arguments{
\item{data}{A data frame containing Ts, D, xi, Z, and M, with the remaining columns representing covariates}  
\item{m0 m1}{The mediation level of interest}  
\item{PX}{A vector specifying the selected covariate columns to be included in the analysis (columns must be after the 5th column)}  
\item{bin}{Variance correction tolerance. If greater than bin, the correction follows the formula on page 7}  
\item{tol}{Convergence tolerance for the iterative algorithm}  
\item{SS}{The probit model parameter. Specifically, different SS correspond to different assumed models in the probit-transformation model class}  
\item{ppt}{A vector of time points of interest for the x-axis}  
\item{setting}{Maximum tolerance for change during the iteration}  
\item{steting}{Optional: "median", "mean", or "level"}  
\item{level}{When "setting" is set to "level", specify a vector of the same length as the covariates}  
}
\description{
Estimating the direct and indirect of the frailty model
}
\examples{
data(expdata)

ee=runif(30000,0,0.5);data[,1]=data[,1]+ee# separating the ties
colnames(data)[1:5] = c("Ts","D","xi","Z","M")
ppt=seq(0,max(data$Ts[data$D==1]),by=0.01) 

ans_IV=IV(ppt,data,m0=0.5,m1=2.5,SS=2.25^0.5,tol=0.001,PX=c(12,13,14),setting="median")
 plot(ppt,ans_IV$WIE,type="l")
 points(ppt,ans_IV$WIE+1.96*ans_IV$WIE.v)
 points(ppt,ans_IV$WIE-1.96*ans_IV$WIE.v)
 abline(h=0,col=2)
}
