library(compiler)
library(survival)
library(MASS)
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


