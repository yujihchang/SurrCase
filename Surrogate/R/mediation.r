### OUTPUT###
#omega.DE: the point-wise estimated direct effect obtained from Omega function (page 4)
#omega.IE: the point-wise estimated indirect effect obtained from Omega function (page 4)
#omega.TE: the point-wise estimated total effect obtained from Omega function (page 4)
#Delta.DE: the point-wise estimated direct effect obtained from Delta (page 6)
#Delta.IE: the point-wise estimated indirect effect obtained from Delta (page 6)
#Delta.TE:  the point-wise estimated total effect obtained from Delta (page 6)
#omega.DE.v: the point-wise estimated sd of Omega.DE
#omega.IE.v: the point-wise estimated sd of Omega.IE
#omega.TE.v: : the point-wise estimated sd of Omega.TE
#Delta.DE.v: the estimated sd of Delta.DE
#Delta.IE.v: the estimated sd of Delta.IE
#Delta.TE.v: : the estimated sd of Delta.IE
#log.pm: the point-wise estimated log proportion of mediation obtained by function pi_omega (page 4) 
#log.pm.sd: the point-wise estimated sd of log.pm
#log.pm.delta: the point-wise estimated log proportion of mediation obtained by Delta (page 4)
#log.pm.delta.sd: the estimated sd of og.pm.delta
#regM: the estimated regression of model (2) in page 5
#coefP: the estimated regression of model (1) in page 5


source("probit.r")


lv=c( -0.0118, -0.0100,  0.0068,  0.5771,  0.6140,  0.6129,  1, 47,  1)
mediation=function(ppt,data,int="TRUE",PX=NULL,tol=0.005,setting,level=lv){
colnames(data)[1:5] = c("Ts","D","xi","Z","M")

Data=data

D=Data$D; xi=Data$xi; X=Data[,-(1:5)]; 
Z=Data$Z; M=Data$M
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

xdim=NCOL(Data)-5
if(xdim>0) X=as.matrix(Data[,-(1:5)]) else X=NULL


pi=round(D+(1-D)*xi/pxi,5)
sub_cohort_i=which(pi>0)
SCH=Data[sub_cohort_i,]
X=X[sub_cohort_i,]
pi=pi[sub_cohort_i]

if(length(X)>0){
d.na=which( is.na(apply(cbind(SCH,X),1,sum)) ==1)
if(length(d.na)>0){
SCH=SCH[-d.na,]
X=X[-d.na,]
pi=pi[-d.na]
}
}
if(length(X)==0){
d.na=which( is.na(apply(cbind(SCH),1,sum)) ==1)
if(length(d.na)>0){
SCH=SCH[-d.na,]
pi=pi[-d.na]
}
}


if(xdim>0) regM=lm(SCH$M~X+SCH$Z,weights=pi) else regM=lm(SCH$M~SCH$Z,weights=pi)
#dim(SCH)
nsch=dim(SCH)[1]
XZ=cbind(rep(1,nsch),X,SCH$Z)

#solve(t(XZ)%*%diag(pi)%*%XZ)%*%t(XZ)%*%(SCH$M*pi)


if(int=="TRUE"){
if(xdim>0)XX=cbind(SCH$Z,SCH$M,SCH$Z*SCH$M,X) else XX=cbind(SCH$Z,SCH$M,SCH$Z*SCH$M)
pb=tran.npmle(Y=SCH$Ts,D=SCH$D,X=XX,r=-1,p.sigma=1,weight=pi,TOL=tol)
}

if(int!="TRUE"){
if(xdim>0)XX=cbind(SCH$Z,SCH$M,X) else XX=cbind(SCH$Z,SCH$M)
pb=tran.npmle(Y=SCH$Ts,D=SCH$D,X=XX,r=-1,p.sigma=1,weight=pi,TOL=0.001,iter=500)
}


#plot(pb$jump,pb$haz)

if(int=="TRUE"){
if(xdim>0){
px=NCOL(X)

ax=regM$coef[1:(px+1)]
az=regM$coef[-(1:(px+1))]
bz   =pb$coef[1]
bm   =pb$coef[2]
bzm  =pb$coef[3]
bx   =pb$coef[-(1:3)]
}

if(xdim==0){
ax=regM$coef[1]
az=regM$coef[-1]
bz    =pb$coef[1]
bm   =pb$coef[2]
bzm =pb$coef[3]
bx =0
}
}


if(int!="TRUE"){
if(xdim>0){
px=NCOL(X)

ax=regM$coef[1:(px+1)]
az=regM$coef[-(1:(px+1))]
bz   =pb$coef[1]
bm  =pb$coef[2]
bx   =pb$coef[-(1:2)]
}

if(xdim==0){
ax=regM$coef[1]
az=regM$coef[-1]
bz    =pb$coef[1]
bm   =pb$coef[2]
bx =0
}
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

bzm=if(int=="TRUE") bzm=bzm else bzm=0

LA=pb$haz
pLA=predit.haz(ppt,LA,pb$jump.time)
omega=function(za,zb,R){
pnorm((
      log(R)+sum(bx*CX[-1])+bz*za+(bm+bzm*za)*(sum(ax*CX)+ az*zb )
      ) / sqrt(1+(bm+bzm*za)^2*sm),lower.tail=FALSE)
}



effects=function(za,zb,R){
#za=0;zb=1;R=pLA
DDE=-(zb-za)*(bz+bzm*( sum(ax*CX)+az*za ) )
DIE =-(zb-za)*az*(bm+bzm*zb)
ODE=omega(za,za,R)-omega(zb,za,R)
OIE=omega(zb,za,R)-omega(zb,zb,R)
list(ODE=ODE,OIE=OIE,DDE=DDE,DIE=DIE)
}

EFF=effects(0,1,pLA)

####################################################################

Uai=
cbind( 
XZ*matrix( (SCH$M -apply(XZ_BB,2,sum)  ) ,dim(XZ)[1],dim(XZ)[2])
)

CMa_1=t(Uai)%*%diag(pi)
CMa=CMa_1%*%t(CMa_1)/nsch
DMa=(t(XZ)%*%diag(pi)%*%XZ)/nsch
#
#DDD=diag( solve(DMa)%*%CMa%*%solve(DMa)/nsch)^0.5
Vreg=solve(DMa)%*%CMa%*%solve(DMa)/nsch
#####################################################################

Umsi=pi*sm-pi*(SCH$M -apply(XZ_BB,2,sum) )^2
#(mean(pi)^(-1)* mean(Umsi^2)* mean(pi)^{-1}/nsch)^0.5
Vreg_11=(mean(pi)^(-1)* mean(Umsi^2)* mean(pi)^{-1}/nsch)


Vpb=pb$VarM


regp=NCOL(Vreg)
Hp=NCOL(Vpb)
allp=1+NCOL(Vreg)+NCOL(Vpb)
ALLV=matrix(0,allp,allp)
ALLV[1,1]=Vreg_11
ALLV[2:(1+regp),2:(1+regp)]=Vreg
ALLV[(2+regp):allp,(2+regp):allp]=Vpb

#za=zb=1


################################################################
THL=NCOL(ALLV)
LAL=length(LA)
ppx=length(CX)
if(int=="TRUE") L.bzm=1 else L.bzm=0
omega_v=function(za,zb){
#za=0;zb=1
A= log(LA)+sum(bx*CX[-1])+bz*za+(bm+bzm*za)*(sum(ax*CX)+az*zb)
B= sqrt(1+(bm+bzm*za)^2*sm)
B2=(1/2)*1/sqrt(1+(bm+bzm*za)^2*sm)

C=A/B

p.sm=dnorm(C)*A/B^2*B2*(bm+bzm*za)^2
if(length(CX)>1) p.ax=t(matrix(-dnorm(C)*1/B*(bm+bzm*za),length(C),length(CX)))*CX else p.ax=-dnorm(C)*1/B*(bm+bzm*za)*CX
p.az=-dnorm(C)*1/B*(bm+bzm*za)*zb
p.bz=-dnorm(C)*1/B*za
p.bm= -dnorm(C)*(1/B*   ( sum(ax*CX)+az*zb)  -A/(B^2)*2*(bm+bzm*za)*B2*sm    )
p.bzm=NULL
if(int=="TRUE"){
p.bzm=-dnorm(C)*(1/B*za*( sum(ax*CX)+az*zb)  -A/(B^2)*2*(bm+bzm*za)*B2*za*sm )
}
if(length(CX)>1){
p.bx=t(matrix(-dnorm(C)*1/B,length(C),(length(CX)-1) ))* CX[-1];length(CX[-1])->p.bx.L
}
if(length(CX)==1){
p.bx= 0;length(CX[-1])->p.bx.L
}

p.rr=-dnorm(C)*1/B*1/LA

if(length(CX)>1) p.ax=p.ax else p.ax=t(as.matrix(p.ax))

vv=matrix(NA,THL,LAL)
for(ss in 1:LAL){
vv[1,ss]=p.sm[ss]
vv[2:(1+ppx),ss]=p.ax[,ss]
vv[(1+ppx)+1,ss]=p.az[ss]
vv[(1+ppx)+2,ss]=p.bz[ss]
vv[(1+ppx)+3,ss]=p.bm[ss]
if(int=="TRUE"){vv[(1+ppx)+4,ss]=p.bzm[ss]}

if(p.bx.L>0){
vv[( (1+ppx)+3+L.bzm+1 ) : ( (1+ppx)+3+L.bzm+p.bx.L ),ss]=p.bx[,ss]
vv[ ( (1+ppx)+3+L.bzm+p.bx.L+1 ):THL,ss]=p.rr[ss]*c(rep(1,ss),rep(0,LAL-ss))
}
if(p.bx.L==0){
vv[ ( (1+ppx)+4+L.bzm):THL,ss]=p.rr[ss]*c(rep(1,ss),rep(0,LAL-ss))
}

}
vv
}
##################################################################


dDE=omega_v(1,0)-omega_v(0,0)
dIE=omega_v(1,1)-omega_v(1,0)
dTE=omega_v(1,1)-omega_v(0,0)


omegaDE.v=omegaIE.v=omegaTE.v=rep(NA,LAL)
for(ss in 1:LAL){
omegaDE.v[ss]=t(dDE[,ss])%*%ALLV%*%dDE[,ss]
omegaIE.v[ss]=t(dIE[,ss])%*%ALLV%*%dIE[,ss]
omegaTE.v[ss]=t(dTE[,ss])%*%ALLV%*%dTE[,ss]

}




################################################################
#
Delta_DE=function(za,zb){
dde.sm=0
dde.ax=-(zb-za)*bzm*CX
dde.az=-(zb-za)*bzm*za
dde.bz=-(zb-za)
dde.bm=0
dde.bzm=NULL
if(int=="TRUE") dde.bzm=-(zb-za)*(sum(ax*CX)+az*za) 
if(length(CX[-1])==0)  dde.bx=NULL else dde.bx=0*CX[-1]
dde.rr=rep(0,LAL)
c(
dde.sm,
dde.ax,
dde.az,
dde.bz,
dde.bm,
dde.bzm,
dde.bx,
dde.rr
)
}

DDE=Delta_DE(0,1)

Delta_DE_v=t(DDE)%*%ALLV%*%DDE



##################################################################
################################################################
#
Delta_IE=function(za,zb){
#za=0;zb=1
die.sm=0
die.ax=0*CX
die.az=-(zb-za)*(bm+bzm*zb)
die.bz=0
die.bm=-(zb-za)*az
die.bzm=NULL
if(int=="TRUE")die.bzm=-(zb-za)*az*zb
if(length(CX[-1])==0)  die.bx=NULL else die.bx=0*CX[-1]
die.rr=rep(0,LAL)
c(
die.sm,
die.ax,
die.az,
die.bz,
die.bm,
die.bzm,
die.bx,
die.rr
)
}

DIE=Delta_IE(0,1)

Delta_IE_v=t(DIE)%*%ALLV%*%DIE
###############################################################
#variance of total effect 

Delta_TE_v=t(DDE+DIE)%*%ALLV%*%(DDE+DIE)

################################################################
#proportional mediation
#estimation
Tiv=diff(c(0,pb$jump.time))
L.AB=length(pb$jump.time)

pmA=(cumsum( (omega(1,0,LA)-omega(1,1,LA))*Tiv ))
pmB=(cumsum( (omega(0,0,LA)-omega(1,1,LA))*Tiv ))
pmD=(cumsum( (omega(0,0,LA)-omega(1,0,LA))*Tiv  ))


log.pm=log( abs( pmA)/abs(pmB ) )
P.log.pm=ifelse(pmB!=0 , log.pm, 0) 
P.log.pm=predit.haz(ppt, P.log.pm, pb$jump.time) 
P.log.pm=ifelse(P.log.pm!=0, P.log.pm, -Inf) # pmD can be 0 if the integral region is {0}

#direction_omega=ifelse(pmA*pmD>0,"Y","N")

#variance

dTE=omega_v(1,1)-omega_v(0,0)

dpmA=dpmB=NULL

for(iiii in 1:L.AB){
TTT=rep(0,L.AB)
TTT[1:iiii]=Tiv[1:iiii]
temp.dIE_TTT=apply(t(dIE)*TTT,2,sum)/pmA[iiii]
temp.dTE_TTT=apply(t(dTE)*TTT,2,sum)/pmB[iiii]

dIE_TTT=ifelse(abs(temp.dIE_TTT)=="Inf"|temp.dIE_TTT=="NaN" ,0,temp.dIE_TTT)
dTE_TTT=ifelse(abs(temp.dTE_TTT)=="Inf"|temp.dTE_TTT=="NaN" ,0,temp.dTE_TTT)

dpmA=cbind(dpmA,dIE_TTT)
dpmB=cbind(dpmB,dTE_TTT)
}


dpmAB=dpmA-dpmB
dpmAB=ifelse(dpmAB=="NaN",0,dpmAB)

vpm=rep(0,L.AB)
for(vvv in 1:L.AB){
temp.dpmAB=dpmAB[,vvv]
vpm[vvv]=t(temp.dpmAB)%*%ALLV%*%temp.dpmAB
}

P.log.pm_vpm=predit.haz(ppt, vpm, pb$jump.time) 

###############################################################
#proportional mediation for Delta version
#estimation

log.pm.delta=log(abs(EFF$DIE/(EFF$DDE+EFF$DIE)))

direction_delta=ifelse(EFF$DIE*EFF$DDE>0,"Y","N")

dpm.delta_IE=DIE/EFF$DIE
dpm.delta_TE=(DIE+DDE)/(EFF$DDE+EFF$DIE)
d_log.pm.delta=dpm.delta_IE-dpm.delta_TE

vpm.delta=t(d_log.pm.delta)%*%ALLV%*%d_log.pm.delta





################################################################
list(
omega.DE =  EFF$ODE,
omega.IE =  EFF$OIE,
omega.TE = EFF$OIE+EFF$ODE,
Delta.DE = EFF$DDE,
Delta.IE = EFF$DIE,
Delta.TE = EFF$DIE+EFF$DDE,
omega.DE.v=predit.haz(ppt,omegaDE.v,pb$jump.time)^0.5,
omega.IE.v=predit.haz(ppt,omegaIE.v,pb$jump.time)^0.5,
omega.TE.v=predit.haz(ppt,omegaTE.v,pb$jump.time)^0.5,
Delta.DE.v=Delta_DE_v^0.5,
Delta.IE.v=Delta_IE_v^0.5,
Delta.TE.v=Delta_TE_v^0.5,
log.pm=P.log.pm,
log.pm.sd=P.log.pm_vpm^0.5,
#direction.omega=direction_omega,
log.pm.delta=log.pm.delta,
log.pm.delta.sd=vpm.delta^0.5,
regM=regM$coef,
coefP=pb$coef
)
}












