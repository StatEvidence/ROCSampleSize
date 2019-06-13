
roc.ss1=function(p=0.85,k=1,len=0.1,alpha=0.05,p1=0.9,beta=0.2,set.low=4,use.low="F",reader=1,rho=0.3){
######################################
## Author: Jeffrey D. Blume
## Date:   July 2007
## Revised for R January 2012
## Version 3.0 
##
## Name: roc.ss1
## Description: Projects sample size with bounds
## 	for estimating and testing a single AUC
#######################################

#######################################	
## Function Inputs:
## p 	    is AUC (and AUC under null hyp.)
## k 	    is ratio of sample sizes (k=m/n); can be any positive number
## 		m=normals; n=abnormal ; AUC=P(normal<abnormal)
## 		typically k>1 "more normals than abnormals"
## len 	    is the length of the CI Len=2*MOE
## alpha    CI coefficient (and Type one error - 2 sided)
## p1 	    is AUC under alternative hypothesis
## beta     is the Type two error
## reader   is the number of readers 
## rho 	    is the correlation between readers (r.drsm)
## set.low  is the bottom of the search interval for 'unitroot' command
##		This is the smallest sample size to search, set.low >= 2
## use.low  When =="T", uses set.low instead of automatic lower limit 
##		may be needed if k is large, for example
##
## To double check: set b=0.5 and len=2*(p1-p) (K,alpha do not matter)
##	EST and HYPT sample sizes will be the same
##	why? b/c p1 = p0+1.96SEs always has 50% power to reject
#######################################

##################
## Set functions
##################

lb.1=function(n,m,p){ min(n,m) - ((min(n,m)-1)^2)/(12*p*(1-p)*(max(m,n)-1))}
assign("lb.1",lb.1)

lb.2=function(n,m,p){1-(m+n-2)*min(p,1-p)/max(p,1-p)+4*sqrt(2*min(p,1-p)*(m-1)*(n-1))/(3*max(p,1-p))}
assign("lb.2",lb.2)

cond=function(n,m,p){ifelse((min(n,m)-1)/(max(n,m)-1)<=2*min(p,1-p),1,0)}
assign("cond",cond)

lb=function(n,m,p){lb.1(n,m,p)*cond(n,m,p)+lb.2(n,m,p)*(1-cond(n,m,p))}
assign("lb",lb)

v.exp=function(n,m,p){(p*(1-p) + (n-1)*((p/(2-p))-p^2)+(m-1)*((2*(p^2)/(1+p))-p^2))/(n*m)}
assign("v.exp",v.exp)

v.nor=function(n,m,p){
(0.0099*exp(-0.5*(1.414*qnorm(p))^2))*((5*(1.414*qnorm(p))^2+8)/n+((1.414*qnorm(p))^2+8)/m)}
assign("v.nor",v.nor)

## Inflation factor for multiple readers
infl=(1/reader+rho*(reader-1)/reader)

##################
## Bounds
##################

n.max=p*(1-p)*((2*qnorm(1-alpha/2)/len)^2)
## n.max is the minimum need in the smaller of the two groups

n.max=ceiling(max(n.max,n.max/k))
## inflates if more abnormals than normal are availible.

thing=function(n,p,length,a,kr){n*n*kr/lb(n,m=n*kr,p=p)-p*(1-p)*((2*qnorm(1-a/2)/length)^2)}

n.min=uniroot(thing,c(0,n.max),p=p,length=len,a=alpha,kr=k)$root

##################
## Calculate for Bi-exponential
##################
## n indexes abnormals -- m indexes normals

thing2=function(n,p,length,a,kr){4*(qnorm(1-a/2)^2)*v.exp(n,m=n*kr,p=p)-length^2}

n.min2=uniroot(thing2,c(1,n.max),p=p,length=len,a=alpha,kr=k)$root

##################
## Calculate for Bi-normal model (Taylor series approx; equal variances)
##################

thing3=function(n,p,length,a,kr){4*(qnorm(1-a/2)^2)*v.nor(n,m=n*kr,p=p)-length^2}

n.min3=uniroot(thing3,c(1,n.max),p=p,length=len,a=alpha,kr=k)$root

##################
## Power
##################
## Bounds

n.pow=(qnorm(beta)*sqrt(p1*(1-p1))+qnorm(alpha/2)*sqrt(p*(1-p)))^2/(p1-p)^2

n.pow=ceiling(max(n.pow,n.pow/k))
## inflates if more abnormals than normal are availible.

## Set search limits for uniroot (sometimes causes a problem)
lower=max(k,1/k,2)
if (use.low=="T") {lower4=set.low} else lower4=lower

## Power computations below use only one tail
## thingX is upper tail; thingXa is lower tail
## use upper when p1>p0 and lower when p1<p0

thing4=function(n,pp0,pp1,b,a,kr)
{(qnorm(1-a/2)*sqrt(pp0*(1-pp0)*lb(n,m=n*kr,p=pp0)/(n*n*kr))+(pp0-pp1))/sqrt(pp1*(1-pp1)*lb(n,m=n*kr,p=pp1)/(n*n*kr))-qnorm(b)}

thing4a=function(n,pp0,pp1,b,a,kr)
{(-qnorm(1-a/2)*sqrt(pp0*(1-pp0)*lb(n,m=n*kr,p=pp0)/(n*n*kr))+(pp0-pp1))/sqrt(pp1*(1-pp1)*lb(n,m=n*kr,p=pp1)/(n*n*kr))-qnorm(1-b)}

if (p1<=p) {n.min4=uniroot(thing4a,c(lower4,n.pow),pp0=p,pp1=p1,b=beta,a=alpha,kr=k)$root}
else {n.min4=uniroot(thing4,c(lower4,n.pow),pp0=p,pp1=p1,b=beta,a=alpha,kr=k)$root}

## Bi-exponential

thing5=function(n,pp0,pp1,b,a,kr)
{(qnorm(1-a/2)*sqrt(v.exp(n,m=n*kr,p=pp0))+(pp0-pp1))/sqrt(v.exp(n,m=n*kr,p=pp1))-qnorm(b)}

thing5a=function(n,pp0,pp1,b,a,kr)
{(-qnorm(1-a/2)*sqrt(v.exp(n,m=n*kr,p=pp0))+(pp0-pp1))/sqrt(v.exp(n,m=n*kr,p=pp1))-qnorm(1-b)}

if (p1<=p) {n.min5=uniroot(thing5a,c(lower,n.pow),pp0=p,pp1=p1,b=beta,a=alpha,kr=k)$root}
else {n.min5=uniroot(thing5,c(lower,n.pow),pp0=p,pp1=p1,b=beta,a=alpha,kr=k)$root}

## Bi-normal (Taylor series approx; equal variances)

thing6=function(n,pp0,pp1,b,a,kr)
{(qnorm(1-a/2)*sqrt(v.nor(n,m=n*kr,p=pp0))+(pp0-pp1))/sqrt(v.nor(n,m=n*kr,p=pp1))-qnorm(b)}

thing6a=function(n,pp0,pp1,b,a,kr)
{(-qnorm(1-a/2)*sqrt(v.nor(n,m=n*kr,p=pp0))+(pp0-pp1))/sqrt(v.nor(n,m=n*kr,p=pp1))-qnorm(1-b)}

if (p1<=p) {n.min6=uniroot(thing6a,c(lower,n.pow),pp0=p,pp1=p1,b=beta,a=alpha,kr=k)$root}
else {n.min6=uniroot(thing6,c(lower,n.pow),pp0=p,pp1=p1,b=beta,a=alpha,kr=k)$root}

##################
## Output
##################

cat("#########################################################################\n")
cat(paste("##  Bounds on sample size projections: Single AUC = p(normal < abnormal)\n"))
cat("##  Ratio of normals to abnormals, k , is",k,"(usually k>=1) \n")
cat("##  Number of Readers =",reader,"; Reader corr. =",rho,"; Inflation factor =",round(infl,4),"\n")
cat("##  Table indicates number of abnormals required \n")
cat("##\n")
cat(paste("##  ESTIMATION  : AUC =",p,"; CI Length =",len,"; alpha =",alpha,"\n"))
cat("## ---------------------------------------------------------------------\n")
cat("##  Bounds      : Max     = ",ceiling(n.max*infl),";  Min    = ",ceiling(n.min*infl),"\n")
cat("##  Projections : Bi-Nor  = ",ceiling(n.min3*infl),";  Bi-Exp = ",ceiling(n.min2*infl),"\n")
cat("##\n")
cat("##  HYP-TESTING : AUC0 =",p,"; AUC1 =",p1,"; beta =", beta,"; alpha =",alpha,"\n")
cat("## ---------------------------------------------------------------------\n")
cat("##  Bounds      : Max     = ",ceiling(n.pow*infl),";  Min    = ",ceiling(n.min4*infl),"\n")
cat("##  Projections : Bi-Nor  = ",ceiling(n.min6*infl),";  Bi-Exp = ",ceiling(n.min5*infl),"\n")
cat("##\n")
cat("#########################################################################\n")
}

##
#