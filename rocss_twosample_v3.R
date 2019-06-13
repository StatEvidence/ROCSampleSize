
roc.ss2=function(pa=0.7,pb=0.8,pa1=0.7,pb1=0.85,rho=0.3,k=1,len=0.15,
alpha=0.05,beta=0.2,set.low=4,use.low="F",use.alt="F",reader=1,rho2=0.3){
######################################
## Author: Jeffrey D. Blume
## Date:   July 2007
## Revised for R January 2012
## Version 3.0 
##
## Name: roc.ss2
## Description: Projects sample size with bounds
## 	for estimating and testing the difference b/w two AUCs
#######################################

#######################################	
## Function Inputs:
## pa       is AUC for modality A
## pa1      is AUC for modality A under alternative
## pb       is AUC for modality B
## pb1      is AUC for modality B under alternative
## rho      is correlation between tests (r.srdm) "same reader different modality"
##		i.e. rho is the test correlation
## rho2     is correlation between the differences in AUC tests amoung different
##		readers (r.difffdr)
##		i.e. rho is the reader correlation
## k        is ratio of sample sizes(k=m/n) ; can be any positive number
## 		m=normals; n=abnormal ; AUC=P(normal<abnormal)
## 		typically k>1 "more normals than abnormals"
## len      is the length of the CI Len=2*MOE
## alpha    CI coefficient	(and Type one error - 2 sided)
## beta     is the Type two error
## reader   is the number of readers
## set.low  is the bottom of the search interval for 'unitroot' command
##		This is the smallest sample size to search, set.low >= 2
## use.low  If tRue uses set.low instead of automatic lower limit 
##		may be needed if k is large, for example
## use. alt if TRUE uses the alternative proportion for CI projections
##
## To double check: 
##	set b=0.5, pb=pa, len=2*(pb1-pa1) 
##	use.alt must be FALSE ; K,alpha do not matter
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

vd.nor=function(n,m,pa,pb,rho)
{v.nor(n=n,m=m,p=pa)+v.nor(n=n,m=m,p=pb)-2*rho*sqrt(v.nor(n=n,m=m,p=pa)*v.nor(n=n,m=m,p=pb))}
assign("vd.nor",vd.nor)

vd.exp=function(n,m,pa,pb,rho)
{v.exp(n=n,m=m,p=pa)+v.exp(n=n,m=m,p=pb)-2*rho*sqrt(v.exp(n=n,m=m,p=pa)*v.exp(n=n,m=m,p=pb))}
assign("vd.exp",vd.exp)

vd.top=function(n,m,pa,pb,rho)
{(pa*(1-pa)+pb*(1-pb))/min(m,n)-2*rho*sqrt(pa*(1-pa)*lb(n=n,m=m,p=pa)*pb*(1-pb)*lb(n=n,m=m,p=pb))/(m*n)}
assign("vd.top",vd.top)

vd.bot=function(n,m,pa,pb,rho)
{(pa*(1-pa)*lb(n=n,m=m,p=pa)+pb*(1-pb)*lb(n=n,m=m,p=pb))/(m*n)-2*rho*sqrt(pa*(1-pa)*pb*(1-pb))/min(n,m)}
assign("vd.bot",vd.bot)

vd.hi=function(n,m,pa,pb,rho)
{(pa*(1-pa)+pb*(1-pb)-2*rho*sqrt(pa*(1-pa)*pb*(1-pb)))/min(n,m)}
assign("vd.hi",vd.hi)

vd.low=function(n,m,pa,pb,rho)
{max(lb(n=n,m=m,p=pa),lb(n=n,m=m,p=pb))*(pa*(1-pa)+pb*(1-pb)-2*rho*sqrt(pa*(1-pa)*pb*(1-pb)))/(m*n)}
assign("vd.low",vd.low)

## Inflation factor for multiple readers
infl=(1/reader+rho2*(reader-1)/reader)

##################
## Bounds
##################

## Set probabilities to alternative if desired
if (use.alt=="T") {
pa.save=pa
pb.save=pb
pa=pa1
pb=pb1
} 

## Bounds on CIs
n.hi=((2*qnorm(1-alpha/2)/len)^2)*(pa*(1-pa)+pb*(1-pb)-2*rho*sqrt(pa*(1-pa)*pb*(1-pb)))

n.hi=ceiling(max(n.hi,n.hi/k))
if (use.low=="T") {low=set.low} else low=2

## inflates if more abnormals than normal are availible.
## n.max is the minimum need in the smaller of the two groups
## also possible to use uniroot, as follows
## thing0=function(n,pa,pb,length,a,kr,rho){4*(qnorm(1-a/2)^2)*vd.hi(n,m=kr*n,pa=pa,pb=pb,rho=rho)-length^2}
## n.max=uniroot(thing0,c(0,1000),pa=pa,pb=pb,length=len,a=alpha,kr=k,rho=rho)$root

thing.low=function(n,pa,pb,length,a,kr,rho){4*(qnorm(1-a/2)^2)*vd.low(n,m=kr*n,pa=pa,pb=pb,rho=rho)-length^2}

n.low=uniroot(thing.low,c(low,n.hi),pa=pa,pb=pb,length=len,a=alpha,kr=k,rho=rho)$root

thing.bot=function(n,pa,pb,length,a,kr,rho){4*(qnorm(1-a/2)^2)*vd.bot(n,m=kr*n,pa=pa,pb=pb,rho=rho)-length^2}

n.bot=uniroot(thing.bot,c(low,n.hi),pa=pa,pb=pb,length=len,a=alpha,kr=k,rho=rho)$root

thing.top=function(n,pa,pb,length,a,kr,rho){4*(qnorm(1-a/2)^2)*vd.top(n,m=kr*n,pa=pa,pb=pb,rho=rho)-length^2}

n.top=uniroot(thing.top,c(low,n.hi*10),pa=pa,pb=pb,length=len,a=alpha,kr=k,rho=rho)$root

##################
## Calculate for Bi-exponential
##################
## n indexes abnormals -- m indexes normals

thing.exp=function(n,pa,pb,length,a,kr,rho){4*(qnorm(1-a/2)^2)*vd.exp(n,m=kr*n,pa=pa,pb=pb,rho=rho)-length^2}

n.exp=uniroot(thing.exp,c(1,n.hi),pa=pa,pb=pb,length=len,a=alpha,kr=k,rho=rho)$root

##################
## Calculate for Bi-normal model (Taylor series approx; equal variances)
##################

thing.nor=function(n,pa,pb,length,a,kr,rho){4*(qnorm(1-a/2)^2)*vd.nor(n,m=kr*n,pa=pa,pb=pb,rho=rho)-length^2}

n.nor=uniroot(thing.nor,c(1,n.hi),pa=pa,pb=pb,length=len,a=alpha,kr=k,rho=rho)$root

##################
## Power
##################

## Set probabilities back if changed
if (use.alt=="T") {
pa=pa.save
pb=pb.save
} 

## Bounds
n.pow=((qnorm(beta)*sqrt((pa1*(1-pa1)+pb1*(1-pb1)-2*rho*sqrt(pa1*(1-pa1)*pb1*(1-pb1))))+qnorm(alpha/2)*sqrt((pa*(1-pa)+pb*(1-pb)-2*rho*sqrt(pa*(1-pa)*pb*(1-pb)))))/((pb1-pa1)-(pb-pa)))^2

## Set search limits for uniroot (sometimes causes a problem)
upper=ceiling(max(n.pow*3,n.pow*3/k))
lower=max(k,1/k,2)
if (use.low=="T") {lower4=set.low} else lower=2

## Power computations below use only one tail
## thingX is upper tail; thingXa is lower tail
## use upper when p1>p0 and lower when p1<p0

## Top bound
thing1=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(qnorm(1-a/2)*sqrt(vd.top(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.top(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(b)}

thing1a=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(-qnorm(1-a/2)*sqrt(vd.top(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.top(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(1-b)}

if ((pb1-pa1)<=(pb-pa)) {np.top=uniroot(thing1a,c(lower,upper),pa0=pa, pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}
else {np.top=uniroot(thing1,c(lower,upper),pa0=pa,pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}

## Hi bound
thing2=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(qnorm(1-a/2)*sqrt(vd.hi(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.hi(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(b)}

thing2a=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(-qnorm(1-a/2)*sqrt(vd.hi(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.hi(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(1-b)}

if ((pb1-pa1)<=(pb-pa)) {np.hi=uniroot(thing2a,c(lower,upper),pa0=pa, pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}
else {np.hi=uniroot(thing2,c(lower,upper),pa0=pa,pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}

## Low bound
thing3=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(qnorm(1-a/2)*sqrt(vd.low(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.low(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(b)}

thing3a=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(-qnorm(1-a/2)*sqrt(vd.low(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.low(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(1-b)}

if ((pb1-pa1)<=(pb-pa)) {np.low=uniroot(thing3a,c(lower,upper),pa0=pa, pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}
else {np.low=uniroot(thing3,c(lower,upper),pa0=pa,pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}

## bottom bound
thing4=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(qnorm(1-a/2)*sqrt(vd.bot(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.bot(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(b)}

thing4a=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(-qnorm(1-a/2)*sqrt(vd.bot(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.bot(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(1-b)}

if ((pb1-pa1)<=(pb-pa)) {np.bot=uniroot(thing4a,c(lower,upper),pa0=pa, pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}
else {np.bot=uniroot(thing4,c(lower,upper),pa0=pa,pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}

## Bi-exponential
thing5=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(qnorm(1-a/2)*sqrt(vd.exp(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.exp(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(b)}

thing5a=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(-qnorm(1-a/2)*sqrt(vd.exp(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.exp(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(1-b)}

if ((pb1-pa1)<=(pb-pa)) {np.exp=uniroot(thing5a,c(lower,upper),pa0=pa, pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}
else {np.exp=uniroot(thing5,c(lower,upper),pa0=pa,pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}

## Bi-normal (Taylor series approx; equal variances)
thing6=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(qnorm(1-a/2)*sqrt(vd.nor(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.nor(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(b)}

thing6a=function(n,pa0,pb0,pa1,pb1,b,a,kr,rho)
{(-qnorm(1-a/2)*sqrt(vd.nor(n,m=n*kr,pa=pa0,pb=pb0,rho=rho))+((pb0-pa0)-(pb1-pa1)))/sqrt(vd.nor(n,m=n*kr,pa=pa1,pb=pb1,rho=rho))-qnorm(1-b)}

if ((pb1-pa1)<=(pb-pa)) {np.nor=uniroot(thing6a,c(lower,upper),pa0=pa, pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}
else {np.nor=uniroot(thing6,c(lower,upper),pa0=pa,pb0=pb,pa1=pa1,pb1=pb1,b=beta,a=alpha,kr=k,rho=rho)$root}

##################
## Output
##################

cat("#########################################################################\n")
cat("##  Bounds on sample size projections: Comparing Two AUCs (A & B) \n")
cat("##  Ratio of normals to abnormals is k =",k," for both tests \n")
cat("##  Number of Readers =",reader,"; Reader corr. =",rho2,"; Inflation factor =",round(infl,4),"\n")
cat("##  Table indicates number of abnormals required \n")
cat("##\n")
cat("##  ESTIMATION  : CI Length =",len,"; alpha =",alpha,"; rho =",rho,"\n")
cat("##              : AUC(A) =",ifelse(use.alt=="T",pa1,pa),"; AUC(B) =",ifelse(use.alt=="T",pb1,pb),"; change =",round(pb-pa,4),"\n")
cat("## ---------------------------------------------------------------------\n")
cat("##  Bounds      : Max     = ",ceiling(n.top*infl),";  Min    = ",ceiling(n.bot*infl),"\n")
cat("##  Bounds**    : Max     = ",ceiling(n.hi*infl) ,";  Min    = ",ceiling(n.low*infl),"\n")
cat("##  Projections : Bi-Nor  = ",ceiling(n.nor*infl),";  Bi-Exp = ",ceiling(n.exp*infl),"\n")
cat("##\n")
cat("##  HYP-TESTING : beta =", beta,"; alpha =",alpha,"; rho =",rho,"\n")
cat("##              : AUC(A)0 =",pa,"; AUC(B)0 =",pb,"; Null change =",round(pb-pa,4),"\n")
cat("##              : AUC(A)1 =",pa1,"; AUC(B)1 =",pb1,"; Alt  change =",round(pb1-pa1,4),"\n")
cat("## ---------------------------------------------------------------------\n")
cat("##  Bounds      : Max     = ",ceiling(np.top*infl),";  Min    = ",ceiling(np.bot*infl),"\n")
cat("##  Bounds**    : Max     = ",ceiling(np.hi*infl) ,";  Min    = ",ceiling(np.low*infl),"\n")
cat("##  Projections : Bi-Nor  = ",ceiling(np.nor*infl),";  Bi-Exp = ",ceiling(np.exp*infl),"\n")
cat("##\n")
cat("#########################################################################\n")
}

##
#