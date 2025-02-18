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










########################################################################################
#This is the package for calculating the npmle without considering the                    
#left truncation under the transformation model.
#In the package the sd. of the R is calculated by its approximated fish information matrix.
#Remark:
#1. the general odds transformation is adopted here.
#2. Y is the observed event times.
#3. D is the censored index.
#4. r is the setting parameter for determining the transformation. 
#Output:
#coef is the regression coefficient.
#bas.haz is the baseline cumulative hazard function.
#iter is the total times of literation.
#JM is the Fisher information matrix for the regression coefficient.
#sd.R is the s.e. of hazard function at all observed event times. 
#
######################################################

predit.haz=function(x,R,Y){
Ys=sort(Y)
prefun=function(t){
cid=sum(t>=Ys)
if(cid<1) ans=0
if(1<=cid& cid<length(Ys)) ans=R[cid]
if(length(Ys)<=cid) ans=R[length(Ys)]
ans
}
apply(matrix(x), 1 ,prefun)
}

########################################################################################



tran.npmle_W=function(Y,D,X,r=-1,p.sigma=1,MU=0,TOL=0.05,weight,iter=50){
#r=-1;p.sigma=p.sigma;liter=25;Bint=NULL;Rint=NULL
#r=-1;p.sigma=1;liter=25;Bint=NULL;Rint=NULL
#r=-1;p.sigma=SS^0.5;TOL=0.0001;Bint=NULL;Rint=NULL;Y=SCH$Ts;D=SCH$D;X=XX;weight=pi;speed=10;TOL=0.005

Bint=NULL;Rint=NULL

######################################################

B=Bint
R=Rint
n=length(Y)
p=length(X)/length(Y)

Y=sort(Y,index.return=TRUE)
D=D[Y$ix]
weight=weight[Y$ix]
if(p==1){X=matrix(X[Y$ix],n,p)}                    
if(p>1){ X=matrix(X[Y$ix,],n,p)}
Y=Y$x
dd=min(which(D==1))

if(dd>1){
Y=Y[-(1:(dd-1))]
X=as.matrix(X[-(1:(dd-1)),])
D=D[-(1:(dd-1))]
weight=weight[-(1:(dd-1))]
n=length(Y)
}
#########################################################################################


if(length(Bint)==0) B=coxph(Surv(Y,D)~X,weights=weight)$coef
if(length(Rint)==0) {
BB=coxph(Surv((Y),D)~X,ties="breslow")
R=basehaz(BB, centered=FALSE)
R=predit.haz(Y,R[,1],R[,2])
}           
dR=diff(c(0,R))
can.jup=which(dR>0)
Yp=Y[can.jup]
dR=dR[dR>0]
R=cumsum(dR)


###########################  transformation function G ##################################
##
##如果要增加不同的 transformations 可在這塊統修改。
##
#--------------------------------------------------------------------------------------- 
# LA is cumulative based-line
# la is dLA
#---------------------------------------------------------------------------------------
G=function(t){
if(r==-1){AA=-log(plnorm(t,meanlog = MU,sdlog =p.sigma,lower.tail=FALSE))}
if(r==0){AA=t}
if(r>0){AA=log(1+r*t)/r }
AA}

g=function(t){
if(r==-1){AA=dlnorm(t,sdlog =p.sigma)/(plnorm(t,meanlog = MU,sdlog =p.sigma,lower.tail=FALSE))}
if(r==0){AA=1}
if(r>0){AA=1/(1+r*t) }
AA}

dg=function(t){
temp=t
mean = MU
    sig2=p.sigma^2
    sig = sqrt(sig2)
    Phi = pnorm(log(temp),mean,sig)
    temp2 = dnorm(log(temp),mean,sig)
    phi = temp2/temp
    dphi = -temp2/(temp^2)-((log(temp)-mean)/(temp*sig2))*temp2/temp
    d2phi = 2*(temp^-3)*temp2+2*(temp^-2)*temp2*((log(temp)-mean)/(temp*sig2))+(temp^-1)*temp2*((log(temp)-mean)/(temp*sig2))^2-(temp^-1)*temp2*(1/(sig2*(temp^2)))+(temp^-1)*temp2*(log(temp)-mean)/((temp^2)*sig2)
    d2G = dphi/(1-Phi)+(phi/(1-Phi))^2
  d2G}


d2g=function(t){
temp=t
    mean = MU
    sig2=p.sigma^2
    sig = sqrt(sig2)
    Phi = pnorm(log(temp),mean,sig)
    temp2 = dnorm(log(temp),mean,sig)
    phi = temp2/temp
    dphi = -temp2/(temp^2)-((log(temp)-mean)/(temp*sig2))*temp2/temp
    d2phi = 2*(temp^-3)*temp2+2*(temp^-2)*temp2*((log(temp)-mean)/(temp*sig2))+(temp^-1)*temp2*((log(temp)-mean)/(temp*sig2))^2-(temp^-1)*temp2*(1/(sig2*(temp^2)))+(temp^-1)*temp2*(log(temp)-mean)/((temp^2)*sig2)
    d3G = d2phi/(1-Phi)+(dphi*phi)/((1-Phi)^2)+2*phi*dphi/((1-Phi)^2)+(2*phi^3)/((1-Phi)^3)
d3G  }





kpa=function(t){ dg(t)/g(t) }
#tt=2;g(tt)-kpa(tt)
####################################################

estB=function(B,R){
Rtise=predit.haz(Y,R,Yp)
EBX=exp(apply(t(X)*B,2,sum))
EBR=EBX*Rtise
part1=(D*kpa(EBR)-g(EBR))*EBR+D
apply(t(X)%*%diag(part1*weight),1,mean)
}

####################################################

d1i=function(B,R){
Rtise=predit.haz(Y,R,Yp)
EBR=exp(apply(t(X)*B,2,sum))*Rtise
part1=d2g(EBR)/g(EBR)-(dg(EBR)/g(EBR))^2 
(part1*D-dg(EBR))*exp(2*apply(t(X)*B,2,sum))
}
#d1i(B,R)
d0i=function(B,R){
Rtise=predit.haz(Y,R,Yp)
EBR=exp(apply(t(X)*B,2,sum))*Rtise
part1=D*dg(EBR)/g(EBR)-g(EBR) 
(part1)*exp(apply(t(X)*B,2,sum))
}
#d0i(B,R)

Ubb=function(B,R){
Rtise=predit.haz(Y,R,Yp)
t(X)%*%
diag(
n^-1*
weight*(d1i(B,R)*Rtise^2
+d0i(B,R)*Rtise
)
)%*%X
}
#Ubb(B,R)


###############################################################################
njump=table(Y[D==1])
njump_V=D
njump_V[D>0]=njump

estR=function(B,R){
#B=BBnew;R=RRold
Rtise=predit.haz(Y,R,Yp)
dRk=diff(c(0,Rtise))

jps=sum(dRk!=0)
EBX=exp(apply(t(X)*B,2,sum))
EBR=EBX*Rtise
Mi=matrix(1:n,jps,n,byrow=TRUE)
Mj=matrix((1:n)[dRk!=0],jps,n)

gg=(D*kpa(EBR)-g(EBR))*EBX*weight

part1=diag(njump_V/dRk)[dRk!=0,]

part2=(  matrix( gg,jps,n,byrow=TRUE) * (Mj<=Mi) )


apply(part2+part1,1,mean)
}
#estR(B,R)
#length(equR(B,R))

Urr=function(B,R){
Rtise=predit.haz(Y,R,Yp)
dRk=diff(c(0,Rtise))
njp=sum(dRk!=0)
Mi=matrix(1:n,n,n)
Mj=t(Mi)
part1=(matrix(d1i(B,R)*weight,n,n,byrow=TRUE))*(Mj>=Mi)
part1=part1[,n:1]
part2=t(apply(part1,1,cumsum))
part2=part2[,n:1]
#part2[1:5,1:5]
part3=part2-diag(njump_V/dRk^2)
(part3)[dRk!=0,dRk!=0]/n
}
#Urr(B,R)
#########################################################
Ubr=function(B,R){
Rtise=predit.haz(Y,R,Yp)
dRk=diff(c(0,Rtise))
part1=(d1i(B,R)*Rtise+d0i(B,R))*weight
part2=t(X)*matrix(part1,p,n,byrow=TRUE)
#----
if(NCOL(X)>1){
part3=apply(part2[,n:1],1,cumsum)
part3=t(part3[n:1,])
aans=part3[,dRk!=0]/n
               }
#----
if(NCOL(X)==1){
part3=cumsum(part2[n:1])
part3=t(as.matrix(part3[n:1]))
aans=t(as.matrix(part3[,dRk!=0]/n))
               }
aans
}

JM=function(B,R){
JML=rbind(
Ubb(B,R),
t(Ubr(B,R))
)
JMR=rbind(
Ubr(B,R),
Urr(B,R)
)
JM=cbind(JML,JMR)
}


##############################################################

equ=function(B,R){
c(estB(B,R),estR(B,R))
}

goal=equ(B,R)
goal_B=estB(B,R)
BBold=B*0
BBnew=B
kk=0
RRold=R+0.1
RRnew=R
Rtise=predit.haz(Y,RRold,Yp)
kk=1
##############################################################
##############################################################
##############################################################
while( max(abs(goal_B))>0.05 | kk<(iter*2)){
BBold=BBnew
kk=kk+1
#print(kk)
BBnew=BBold-c(solve(Ubb(BBold,RRnew))%*%estB(BBold,RRnew))
#goal_B=estB(BBnew,RRnew)

#while(abs(max(goal_R))>TOL){

RRold=RRnew
Rtise=predit.haz(Y,RRold,Yp)
dRk_old=diff(c(0,Rtise))
dRk_old=dRk_old[dRk_old!=0]
if(kk<10) {dRk_new=dRk_old-c(solve(Urr(BBold,RRold))%*%estR(BBold,RRold))*0.001}
if(kk>=10) {dRk_new=dRk_old-c(solve(Urr(BBold,RRold))%*%estR(BBold,RRold))*0.1}


#dRk_m=min(dRk_new[dRk_new>0])
dRk_new=ifelse(dRk_new>0,dRk_new,dRk_old)
RRnew=cumsum(dRk_new)

goal_B=estB(BBnew,RRnew)
#goal_R=estR(BBnew,RRnew)
#print(max(abs(goal_R)))
#print(max(abs(goal_B)))
#goal=max(abs(goal_R),abs(goal_B))
#print(c(max(abs(goal_R)),max(abs(goal_B))))
#print(abs(max(c(BBnew,RRnew)-c(BBold,RRold))))
#print(max(abs(goal_B)))
#plot(Y,Rtise);abline(a=0,b=1)

}
#BBnew=BBold;RRnew=RRold
#B=BBold;R=RRold

##############################################################
B=BBold;R=RRold


goal=equ(B,R)
goal_B=estB(B,R)
BBold=B*0
BBnew=B
kk=0
RRold=R+0.1
RRnew=R
Rtise=predit.haz(Y,RRold,Yp)
kk=1
##############################################################
while(max(abs(c(BBnew,RRnew)-c(BBold,RRold)))>TOL | max(abs(goal_B))>TOL&kk<iter){
BBold=BBnew
kk=kk+1
BBnew=BBold-c(solve(Ubb(BBold,RRnew))%*%estB(BBold,RRnew))
#goal_B=estB(BBnew,RRnew)

#while(abs(max(goal_R))>TOL){

RRold=RRnew
Rtise=predit.haz(Y,RRold,Yp)
dRk_old=diff(c(0,Rtise))
dRk_old=dRk_old[dRk_old!=0]
dRk_new=dRk_old-c(solve(Urr(BBold,RRold))%*%estR(BBold,RRold))

#dRk_m=min(dRk_new[dRk_new>0])
dRk_new=ifelse(dRk_new>0,dRk_new,dRk_old)
RRnew=cumsum(dRk_new)

goal_B=estB(BBnew,RRnew)
#goal_R=estR(BBnew,RRnew)
#print(max(abs(goal_R)))
#print(max(abs(goal_B)))
#goal=max(abs(goal_R),abs(goal_B))
#print(c(max(abs(goal_R)),max(abs(goal_B))))
#print(abs(max(c(BBnew,RRnew)-c(BBold,RRold))))
#print(max(abs(goal_B)))
#plot(Y,Rtise);abline(a=0,b=1)
}
#BBnew=BBold;RRnew=RRold
#B=BBold;R=RRold

##############################################################

jm=JM(BBnew,RRnew)

UU=function(B,R){

Rtise=predit.haz(Y,R,Yp)
EBX=exp(apply(t(X)*B,2,sum))
EBR=EBX*Rtise
part11=(D*kpa(EBR)-g(EBR))*EBR+D
#t(X)%*%diag(part11*weight)

#################################################
#B=BBnew;R=RRold
Rtise=predit.haz(Y,R,Yp)
dRk=diff(c(0,Rtise))
jps=sum(dRk!=0)
EBX=exp(apply(t(X)*B,2,sum))
EBR=EBX*Rtise
Mi=matrix(1:n,jps,n,byrow=TRUE)
Mj=matrix((1:n)[dRk!=0],jps,n)
gg=(D*kpa(EBR)-g(EBR))*EBX*weight
part1=diag(njump_V/dRk)[dRk!=0,]
part2=(  matrix( gg,jps,n,byrow=TRUE) * (Mj<=Mi) )
##################################################
ans=rbind(
t(X)%*%diag(part11*weight),
part1+part2
)
ans/sqrt(n)
}

Uis=UU(BBnew,RRnew)
Vall= solve(jm)%*%Uis%*%t(Uis)%*%solve(jm)/n
BB_var=(diag((Vall))^0.5)[1:p]


arrow=c( rep(0,p) , rep(1,length(njump) ) )
LA_var=rep(NA,length(njump))
for(jj in (p+1):length(njump)){
LA_var[jj]= t(arrow[1:jj]) %*% Vall[1:jj,1:jj] %*% arrow[1:jj] 
}

collectNA=NULL
for(cc in 1:length(njump) ){
if(is.na(LA_var[cc])==1){ collectNA=c(collectNA,cc)}
if(is.na(LA_var[cc])==0) {LA_var[collectNA]=LA_var[cc];RRnew[collectNA]=RRnew[cc];collectNA=NULL }
}



report=list(coef=BBnew,
            coef.sd=BB_var,
            haz=RRnew,
            haz.sd=LA_var^0.5,
            jump.time=Yp,
            VarM=Vall,
            JM=jm,
            Uis=Uis,kk=kk)
###############################################################
}

tran.npmle=cmpfun(tran.npmle)











