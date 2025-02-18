### OUTPUT###
#WIE: the point-wise estimated effect obtained from IV analyses through the W_IE() function (page 7)
#WIE.v: the point-wise estimated sd of WIE
#delta.IV: estimated bz/az (page 7)
#delta.IV.v: estimated sd. of delta.IV

#rm(list=ls(all=TRUE))

#source("probitW.r")
#data=IV_data(n=2000,a0=0.5,az=1,au=1,sm=2^0.5,b0=0,bm=1,bu=1,st=1,pxi=0.2,cr=0.325)

lv=c( -0.0118, -0.0100,  0.0068,  0.5771,  0.6140,  0.6129,  1, 47,  1)

IV=function(ppt,data,m0=0,m1=1,PX=NULL,tol=0.005,SS=1,mu=0,bin=0.5,setting,level=lv){

colnames(data)[1:5] = c("Ts","D","xi","Z","M")
Data=data
Data$Ts=data$Ts; D=data$D; xi=data$xi; X=data[,-(1:5)]; Z=data$Z; M=data$M

pxi=mean(xi)

if(length(PX)>0){
WW=as.matrix(data[,PX])

Wcorf=coef(glm(xi~as.matrix(WW),data,family=binomial))
n=NROW(WW)
WW_I=cbind(rep(1,n),WW)
if(length(PX)>2)  log_pxi=apply(t(WW_I)*Wcorf,2,sum)
if(length(PX)==1) log_pxi=WW_I*Wcorf
pxi=exp(log_pxi)/(exp(log_pxi)+1)
}

xdim=NCOL(data)-5
if(xdim>0)X=as.matrix(data[,-(1:5)]) else X=NULL

pi=D+(1-D)*xi/pxi
sub_cohort_i=which(pi>0)
SCH=Data[sub_cohort_i,]
X=X[sub_cohort_i,]
pi=pi[sub_cohort_i]
##################################################################################
sid=sort(SCH$Ts,index.return=TRUE)$ix
SCH=SCH[sid,]
if(xdim>1)X=X[sid,] 
if(xdim==1)X=X[sid]
pi=pi[sid]

dd=min(which(SCH$D==1))

if(dd>1){
SCH=SCH[-(1:(dd-1)),]
if(xdim>1) X=as.matrix(X[-(1:(dd-1)),])
if(xdim==1) X=as.matrix(X[-(1:(dd-1))])
pi=pi[-(1:(dd-1))]
}

##################################################################################
if(xdim>0) regM=lm(M~X+Z,SCH,weights=pi) else regM=lm(M~Z,SCH,weights=pi)
#dim(SCH)
nsch=dim(SCH)[1]
XZ=cbind(rep(1,nsch),X,SCH$Z)



if(xdim>0) XX=cbind(SCH$Z,X) else XX=as.matrix(SCH$Z)
pb=tran.npmle_W(Y=SCH$Ts,D=SCH$D,X=XX,r=-1,p.sigma=SS,MU=mu,weight=pi,TOL=tol)

#plot(pb$jump,pb$bas,ylim=c(0,max(SCH$Ts)))


if(xdim>0){
px=xdim

ax=regM$coef[1:(px+1)]
az=regM$coef[-(1:(px+1))]
bz =pb$coef[1]
bx =pb$coef[-1]
}

if(xdim==0){
ax=regM$coef[1]
az=regM$coef[-1]
bz =pb$coef[1]
bx =0
}

#sum((pi*regM$residuals)^2)
CX=1
if(setting=="median"){
if(xdim>0){
if(px>1) CX=c(1,apply(X,2,median))
if(px==1) CX=c(1,mean(X))
}
}
if(setting=="mean"){
if(xdim>0){
if(px>1) CX=c(1,apply(X,2,mean))
if(px==1) CX=c(1,mean(X))
}
}
if(setting=="level"){
if(xdim>0){
if(px>1) CX=c(1,level)
if(px==1) CX=c(1,level)
}
}

XZ_BB=t(XZ)*c(ax,az)

sm=(sum(pi*( ( SCH$M -apply(XZ_BB,2,sum) ))^2)/sum(pi))

#sch=SCH[SCH$D>0,]
#xx=XX[SCH$D>0,]
#schix=sort(sch$Ts,index.return=TRUE)$ix
#xx=xx[schix,]
#st=sd(log(LA)-apply(xx*(-1)*pb$coef,1,sum))
# 這邊是為了估計 st, 但他是已知所以不用估計



####################################################################
# For variance
Uai=
cbind( 
XZ*matrix( (SCH$M -apply(XZ_BB,2,sum)  ) ,dim(XZ)[1],dim(XZ)[2])
)/sqrt(nsch)

DMa=(t(XZ)%*%diag(pi)%*%XZ)/nsch
##------------------------------------------------------
Umsi=(pi*sm-pi*(SCH$M -apply(XZ_BB,2,sum) )^2)/sqrt(nsch)

DMsi=mean(pi)

Jpb=(pb$JM)

SS2=SS^2

regp=NCOL(DMa)
Hp=NCOL(Jpb)
allp=1+regp+NCOL(Jpb)
AD=matrix(0,allp,allp)
AD[1,1]=DMsi
AD[2:(1+regp),2:(1+regp)]=DMa
AD[(2+regp):allp,(2+regp):allp]=Jpb

AUU=cbind(Umsi,Uai,t(pb$Uis))

# all variance and covariance matrix
#
VIV=solve(AD)%*%t(AUU)%*%AUU%*%solve(AD)/nsch

################################################################
# Effect
LA=pb$haz
Wm=function(m){
A=log(LA)+sum(bx*CX[-1])+(bz/az)*(-sum(ax*CX) + m)
Bin=SS2-(bz/az)^2*sm
if(Bin>bin) B=sqrt(Bin) else B=bin^0.5
pnorm(A/B,lower.tail=FALSE)
}
################################################################
#D Effect
#SS=1
THL=NCOL(VIV)
LAL=length(LA)
ppx=length(CX)

Bin=SS2-(bz/az)^2*sm

if(Bin>bin){
dWm=function(m){

#m=1
A=log(pb$haz)+sum(bx * CX[-1] )+(bz/az)*(-sum(ax*CX) + m)
B=sqrt(SS2-(bz/az)^2*sm)
B2=1/2*1/B
C=A/B


p.sm=-dnorm(C)*A/B^2*B2*(bz/az)^2
if(length(CX)>1) p.ax=t(matrix(dnorm(C)*1/B*(bz/az),length(C),(1+xdim)))*CX else p.ax=dnorm(C)*1/B*(bz/az)*CX
p.az=  dnorm(C)*(1/B* bz/(az^2)* (-sum(ax*CX) + m)  +2* A/(B^2)*B2*(bz^2/az^3)*sm    )
p.bz= -dnorm(C)*(1/B*  1/(az  )* (-sum(ax*CX) + m)  +2* A/(B^2)*B2*(  bz/az^2)*sm    )

if(length(CX)>1)   p.bx=t(matrix(-dnorm(C)*1/B,length(C),xdim))*CX[-1] else p.bx=-dnorm(C)*1/B*CX[-1]

length(CX[-1])->p.bx.L
p.rr=-dnorm(C)*1/B*1/LA


vv=matrix(NA,THL,LAL)
for(ss in 1:LAL){
vv[1,ss]=p.sm[ss]
if((xdim+1)>1) vv[2:(1+ppx),ss]=p.ax[,ss] else vv[2:(1+ppx),ss]=p.ax[ss]
vv[(1+ppx)+1,ss]=p.az[ss]
vv[(1+ppx)+2,ss]=p.bz[ss]
if(p.bx.L>0){
vv[( (1+ppx)+3 ) : ( (1+ppx)+2+p.bx.L ),ss]=p.bx[,ss]
vv[ ( (1+ppx)+3+p.bx.L ):THL,ss] = p.rr[ss]*c(rep(1,ss),rep(0,LAL-ss))
}
if(p.bx.L==0){
vv[ ( (1+ppx)+3):THL,ss]=p.rr[ss]*c(rep(1,ss),rep(0,LAL-ss))
}

}
vv
}
}


if(Bin<=bin){

dWm=function(m){

#m=1
A=log(pb$haz)+sum(bx*CX[-1])+(bz/az)*(-sum(ax*CX) + m)
B=bin
B2=0
C=A/B


p.sm=-dnorm(C)*A/B^2*B2*(bz/az)^2
if(length(CX)>1) p.ax=t(matrix(dnorm(C)*1/B*(bz/az),length(C),(1+xdim)))*CX else p.ax=dnorm(C)*1/B*(bz/az)*CX
p.az=  dnorm(C)*(1/B* bz/(az^2)* (-sum(ax*CX) + m)  +2* A/(B^2)*B2*(bz^2/az^3)*sm    )
p.bz= -dnorm(C)*(1/B*  1/(az  )* (-sum(ax*CX) + m)  +2* A/(B^2)*B2*(  bz/az^2)*sm    )

if(length(CX)>1)   p.bx=t(matrix(-dnorm(C)*1/B,length(C),xdim))*CX[-1] else p.bx=-dnorm(C)*1/B*CX[-1]

length(CX[-1])->p.bx.L
p.rr=-dnorm(C)*1/B*1/LA


vv=matrix(NA,THL,LAL)
for(ss in 1:LAL){
vv[1,ss]=p.sm[ss]
if((xdim+1)>1) vv[2:(1+ppx),ss]=p.ax[,ss] else vv[2:(1+ppx),ss]=p.ax[ss]
vv[(1+ppx)+1,ss]=p.az[ss]
vv[(1+ppx)+2,ss]=p.bz[ss]
if(p.bx.L>0){
vv[( (1+ppx)+3 ) : ( (1+ppx)+2+p.bx.L ),ss]=p.bx[,ss]
vv[ ( (1+ppx)+3+p.bx.L ):THL,ss] = p.rr[ss]*c(rep(1,ss),rep(0,LAL-ss))
}
if(p.bx.L==0){
vv[ ( (1+ppx)+3):THL,ss]=p.rr[ss]*c(rep(1,ss),rep(0,LAL-ss))
}

}
vv
}
}


WIE=Wm(m1)-Wm(m0)
dWIE=dWm(m1)-dWm(m0)

WIE.v=rep(NA,LAL)
for(ss in 1:LAL){
WIE.v[ss]=t(dWIE[,ss])%*%VIV%*%dWIE[,ss]
}
################################################################# IV delta~
delta.IV=bz/az
d_delta.IV=rep(0,NCOL(VIV))
d_delta.IV[1+ppx+1]=-bz/(az^2)
d_delta.IV[1+ppx+2]=1/az

delta.IV.v=t(d_delta.IV)%*%VIV%*%d_delta.IV

################################################################
list(
WIE=predit.haz(ppt,WIE,pb$jump.time),
WIE.v=predit.haz(ppt,WIE.v,pb$jump.time)^0.5,
delta.IV=delta.IV,
delta.IV.v=delta.IV.v^0.5
#Bin=SS2-(bz/az)^2*sm
)

}


















