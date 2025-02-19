
\name{mediation}
\alias{mediation}
\title{mediation Analysis Using the Transformation Probit Model in Case-Cohort Studie}
\usage{
mediation(ppt,data,int="TRUE",PX=NULL,tol=0.005,setting,level=lv)
}
\arguments{
\item{data}{A data frame containing Ts, D, xi, Z, and M, with the remaining columns representing covariates}  
\item{PX}{A vector specifying the selected covariate columns to be included in the analysis (columns must be after the 5th column)}  
\item{tol}{Convergence tolerance for the iterative algorithm}  
\item{ppt}{A vector of time points of interest for the x-axis}  
\item{int}{TRUE indicates that the interaction between mediation and exposure is considered}  
\item{steting}{Optional: "median", "mean", or "level"}  
\item{level}{When "setting" is set to "level", specify a vector of the same length as the covariates}  
}
\description{
Surrogate Marker Assessment Using Mediation and Instrumental Variable Analyses in a Case-cohort Design (Annals of Applied Statistics, 2023).
}
\examples{
data(expdata)

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
