rm(list=ls())     

load_file = function(file_name){
     load(file_name)
     obj_names = ls()
     obj_names = obj_names[obj_names != 'file_name']
     if(length(obj_names) > 1){
          #return a list of objects
          dat = sapply(obj_names, function(x)get(x), simplify=FALSE, USE.NAMES=TRUE)
     }else{
          #return a single object
          dat = get(obj_names)     
     }
     return(dat)
}


fptcdf=function(z,x0max,chi,driftrate,sddrift) {
     if (x0max<1e-10) return(pnorm(chi/z,mean=driftrate,sd=sddrift,lower.tail=F))
     zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu ; xx=chiminuszu-x0max
     chizu=chiminuszu/zs ; chizumax=xx/zs
     tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
     tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
     1+(tmp1+tmp2)/x0max
}

fptpdf=function(z,x0max,chi,driftrate,sddrift) {
     if (x0max<1e-10) return( (chi/z^2)*dnorm(chi/z,mean=driftrate,sd=sddrift))
     zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu
     chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
     (driftrate*(pnorm(chizu)-pnorm(chizumax)) + 
               sddrift*(dnorm(chizumax)-dnorm(chizu)))/x0max
}


n1PDFfixedt0=function(t,x0max,chi,drift,sdI,truncdrifts=TRUE) {
     # Generates defective PDF for responses on node #1.
     # "truncdrifts" sets whether that part of the multi-variate
     # normal distribution on drift rates which would otherwise
     # lead to non-terminating trials is truncated.
     G=1-fptcdf(z=t,x0max=x0max[2],chi=chi[2],driftrate=drift[,2],sddrift=sdI[2])
     out=G*fptpdf(z=t,x0max=x0max[1],chi=chi[1],driftrate=drift[,1],sddrift=sdI[1])
     if (truncdrifts) {
          out=out/(1-prod(pnorm(-drift/sdI)))
          out[t<=0]=0
          return(out)
     } else {
          return(out)
     }
}


get.dens.2choice=function(rt,crct,b,A,v,s,t0){
     out=numeric(length(rt))
     out[crct] = n1PDFfixedt0(rt[crct]-t0,x0max=c(A,A),chi=c(b,b),drift=cbind(v[crct],1-v[crct]),sdI=s)
     out[!crct] = n1PDFfixedt0(rt[!crct]-t0,x0max=c(A,A),chi=c(b,b),drift=cbind(1-v[!crct],v[!crct]),sdI=s[2:1])
     pmax(out,1e-10)
}

log.dens.prior=function(x,hyper){
     out=0
     for (p in names(x)) out =
               out+dtnorm(x[p],hyper[paste(p,"mu",sep=".")],hyper[paste(p,"sigma",sep=".")],0,Inf,log=TRUE)
     out
}


log.dens.like=function(x,data,par.names){
     names(x)=par.names
     if (!doesStartPointVary) {
          A=x["A"]
     }
     if (!doesThresholdVary & !doesStartPointVary) {
          b=x["b"] + x["A"]
     }
     if (!doesT0Vary) {
          t0=x["t0"]
     }
     if (!doesDriftRateCorrectVary & !doesDriftRateErrorVary) {
          vs=c(x["vc"],x["ve"])
     }
     if (!doesStDevErrorVary) {
          s=c(1,x["sve"])
     }
     
     out=0
     
     for (cond in conds) { 						# Change to your number of conditions (within)
          
          if (doesStartPointVary) {
               A=x[paste("A",cond,sep=".")]
          }
          if (doesThresholdVary & doesStartPointVary) {
               b=x[paste("b",cond,sep=".")] + x[paste("A",cond,sep=".")]
          } else if (doesThresholdVary & !doesStartPointVary) {
               b=x[paste("b",cond,sep=".")] + x["A"]
          } else if (!doesThresholdVary & doesStartPointVary) {
               b=x["b"] + x[paste("A",cond,sep=".")]
          }
          if (doesT0Vary) {
               t0=x[paste("t0",cond,sep=".")]
          }
          if (doesDriftRateCorrectVary & doesDriftRateErrorVary) {
               vs=c(x[paste("vc",cond,sep=".")],x[paste("ve",cond,sep=".")])
          } else if (!doesDriftRateCorrectVary & doesDriftRateErrorVary) {
               vs=c(x["vc"],x[paste("ve",cond,sep=".")])
          } else if (doesDriftRateCorrectVary & !doesDriftRateErrorVary) {
               vs=c(x[paste("vc",cond,sep=".")],x["ve"])
          }
          if (doesStDevErrorVary) {
               s=c(1,x[paste("sve",cond,sep=".")])
          }
          
          s = c(x[paste('svc')],x[paste('sve')])
          c = x['c']
          beta = x['beta']
          sim = exp(-c*data$dst)
          p = sim / (beta + sim)
          tmp=get.dens.2choice(rt=data$rt,crct=data$response=='s',b=b,A=A,v=p,s=s,t0=t0)
          if(any(is.nan(p))){
               out = -Inf
          }else{
               out=out+sum(log(tmp))
          }
     }
     out
     
}

log.dens.hyper=function(theta,phi,prior){
     sum((dtnorm(theta,phi[1],phi[2],0,Inf,log=TRUE))) + 
          (dtnorm(phi[1],prior$mu[1],prior$mu[2],0,Inf,log=TRUE)) + 
          (dtnorm(phi[2],prior$sigma[1],prior$sigma[2],0,Inf,log=TRUE))
}

crossover=function(i,pars,use.theta,use.like,data,hyper,par.names){
     #require(msm)
     use.weight=use.like[i] + log.dens.prior(use.theta[i,],hyper[i,])
     gamma = 2.38/sqrt(2*length(pars))
     index=sample(c(1:n.chains)[-i],2,replace=F)
     theta=use.theta[i,]						
     theta[pars]=use.theta[i,pars] + gamma*(use.theta[index[1],pars]-use.theta[index[2],pars]) + runif(1,-b,b)
     #  theta=matrix(theta,1,length(theta))
     like=log.dens.like(theta,data,par.names=par.names)
     weight=like + log.dens.prior(theta,hyper[i,])
     if(is.na(weight))weight=-Inf
     if(runif(1) < exp(weight-use.weight)) {							
          use.theta[i,]=theta
          use.like[i]=like
     }
     c(use.like[i],use.theta[i,])
}

crossover_hyper=function(i,pars,use.theta,use.phi,prior){
     #require(msm)
     use.weight=log.dens.hyper(use.theta[i,],use.phi[i,pars],prior)
     gamma = 2.38/sqrt(2*length(pars))
     index=sample(c(1:n.chains)[-i],2,replace=F)
     phi=use.phi[i,]
     phi[pars]=use.phi[i,pars] + gamma*(use.phi[index[1],pars]-use.phi[index[2],pars]) + runif(1,-b,b)
     #  phi=matrix(phi,1,length(phi))
     weight=log.dens.hyper(use.theta[i,],phi[pars],prior)
     if(is.na(weight))weight=-Inf
     if(runif(1) < exp(weight-use.weight)) use.phi[i,]=phi
     use.phi[i,]
}

migration.crossover=function(pars,use.theta,use.like,data,hyper,par.names){
     # migration step
     n.migration.chains=ceiling(runif(1,0,n.chains))
     use.chains=sample(1:n.chains,n.migration.chains)
     migration.use.weight=rep(NA,n.migration.chains)
     migration.weight=rep(NA,n.migration.chains)
     for (mi in 1:n.migration.chains) {
          migration.use.weight[mi]=use.like[use.chains[mi]] + log.dens.prior(use.theta[use.chains[mi],],hyper[use.chains[mi],])
          newChain = mi - 1
          if (newChain == 0) newChain = n.migration.chains
          migration.weight[mi]=use.like[use.chains[newChain]] + log.dens.prior(use.theta[use.chains[newChain],],hyper[use.chains[mi],])
          if(runif(1) < exp(migration.weight[mi]-migration.use.weight[mi])) {            	
               use.theta[use.chains[mi],]=use.theta[use.chains[newChain],]
               use.like[use.chains[mi]]=use.like[use.chains[newChain]]
          }
     }
     cbind(use.like,use.theta)
}

migration.crossover_hyper=function(pars,use.theta,use.phi,prior){
     # migration step
     n.migration.chains=ceiling(runif(1,0,n.chains))
     use.chains=sample(1:n.chains,n.migration.chains)
     migration.use.weight=rep(NA,n.migration.chains)
     migration.weight=rep(NA,n.migration.chains)
     for (mi in 1:n.migration.chains) {
          migration.use.weight[mi]=log.dens.hyper(use.theta[use.chains[mi],],use.phi[use.chains[mi],pars],prior)
          newChain = mi - 1
          if (newChain == 0) newChain = n.migration.chains
          migration.weight[mi]=log.dens.hyper(use.theta[use.chains[mi],],use.phi[use.chains[newChain],pars],prior)
          if(runif(1) < exp(migration.weight[mi]-migration.use.weight[mi])) {          		
               use.phi[use.chains[mi],pars]=use.phi[use.chains[newChain],pars]
          }
     }
     use.phi
}



library(msm)
library(optparse)

#setwd('C:/Users/Jeffrey/Dropbox/backup/Research Projects/mds_3M_v2/no_viewpoint/dst_cat_sep_model/')

option_list = list(
     make_option(c('-sp','--start_point'), type="numeric", default=0, 
                 help="vary startpoint", metavar="numeric"),
     make_option(c('-b','--threshold'), type="numeric", default=0, 
                 help="vary threshold", metavar="numeric"),
     make_option(c('-c','--drift_correct'), type="numeric", default=0, 
                 help="vary drift correct", metavar="numeric"),
     make_option(c('-e','--drift_error'), type="numeric", default=0, 
                 help="vary drift error", metavar="numeric"),
     make_option(c('-std','--st_dev'), type="numeric", default=0, 
                 help="vary st_dev", metavar="numeric"),
     make_option(c('-t','--t0'), type="numeric", default=0, 
                 help="vary t0", metavar="numeric"),
     make_option(c('-n','--nmc'), type="numeric", default=3000, 
                 help="number of iterations", metavar="numeric"),
     make_option(c('-nc','--n_cores'), type="numeric", default=1, 
                 help="number of cores", metavar="numeric"),
     make_option(c('-o','--file_out'), type="character", default='', 
                 help="output file", metavar="character"),
     make_option(c('-i','--file_in'), type="character", default='', 
                 help="input file", metavar="character"),
     make_option(c('-cat','--category'),type='numeric',default=1,help='',metavar='character'),
     make_option(c('-v','--viewpoint'), default=FALSE, action = 'store_true', #false unless -v is on cmd line
                 help="viewpoint condition")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dat = load_file(opt$file_in)
model_name = 'dst_cat_sep_model'
net_name = dat$net
color_mode = dat$color_mode
dst_euc_mat = dat$dst
dat = dat$dat

#make directory if none exists
if (opt$viewpoint) {
     viewpoint = 'viewpoint'
}else{
     viewpoint = 'no_viewpoint'    
}
opt$file_out =  file.path(dirname(opt$file_out),net_name,model_name,'samples',color_mode,viewpoint,basename(opt$file_out))
if (!dir.exists(dirname(opt$file_out))) {
     dir.create(dirname(opt$file_out),recursive = TRUE)
}

if(!opt$viewpoint){
     dat = dat[dat$viewpoint == 'same' | (dat$viewpoint == 'diff' & dat$object == 'diff'),] #same objects must have same viewpoint
}

dat = dat[dat$trial > 180, ]
dat = dat[dat$response != 'timeout',]

cat_list = list(c(1,100),c(101,200),c(201,300),c(301,400),c(401,500))

cat = opt$cat
dat = dat[dat$ID_1 >= cat_list[[cat]][1]  & dat$ID_1 <= cat_list[[cat]][2],]

library(plyr)

data = llply(unique(dat$subject), 
             function(x){
                  data.frame(subject = dat$subject[dat$subject==x], 
                             rt=dat$rt[dat$subject==x] / 1000, 
                             response = dat$response[dat$subject==x],
                             dst = dst_euc_mat[cbind(dat$ID_1[dat$subject==x], dat$ID_2[dat$subject==x])]
                  )
             })



doesStartPointVary = opt$start_point     # Does start point vary over conditions for this model?
doesThresholdVary = opt$threshold    # Does threshold vary over conditions for this model? (NOTE: This MUST also vary if start point varies)
doesDriftRateCorrectVary = opt$drift_correct     # Does the mean drift rate for the correct accumulator vary over conditions for this model?
doesDriftRateErrorVary = opt$drift_error     # Does the mean drift rate for the error accumulator vary over conditions for this model?
doesStDevErrorVary = opt$st_dev     # Does the standard deviation of the drift rate of the error accumulator vary over conditions for this model?
doesT0Vary = opt$t0     # Does t0 vary over conditions for this model?

thin_steps = 10
n.pars= 7
n.hpars=n.pars*2
conds=1
n.cond=1
S=length(data)

n.chains=2*n.pars
nmc=opt$nmc

migration.freq = 20
migration.start = 400
migration.end = 800

b=.001

A.start=1
B.start=1
t0.start=.2
vc.start=3
ve.start=2
sve.start=1

start.points=c(rep(A.start,ifelse(doesStartPointVary,n.cond,1)),
               rep(B.start,ifelse(doesThresholdVary,n.cond,1)),
               rep(t0.start,ifelse(doesT0Vary,n.cond,1)),
               rep(vc.start,ifelse(doesDriftRateCorrectVary,n.cond,1)),
               rep(ve.start,ifelse(doesDriftRateErrorVary,n.cond,1)),
               rep(sve.start,ifelse(doesStDevErrorVary,n.cond,1)))

theta=array(NA,c(n.chains,n.pars,S,nmc))
phi=array(NA,c(n.chains,n.hpars,nmc))
weight=array(-Inf,c(nmc,n.chains,S))

theta.names=NULL

if (doesStartPointVary) {theta.names=c(theta.names,paste("A",conds,sep="."))} else {theta.names=c(theta.names,"A")}
if (doesThresholdVary) {theta.names=c(theta.names,paste("b",conds,sep="."))} else {theta.names=c(theta.names,"b")}
if (doesT0Vary) {theta.names=c(theta.names,paste("t0",conds,sep="."))} else {theta.names=c(theta.names,"t0")}
if (doesStDevErrorVary) {theta.names=c(theta.names,paste("sve",conds,sep="."))} else {theta.names=c(theta.names,"sve")}

theta.names = c(theta.names,'beta')
theta.names = c(theta.names,'c')
theta.names = c(theta.names,'svc')
colnames(theta) = theta.names

phi.names=paste(rep(theta.names,each=2),c("mu","sigma"),sep=".") # Change accordingly.
colnames(phi) <- phi.names


for(i in 1:n.chains){
     tmp=grep("A",phi.names)
     phi[i,tmp,1]=rtnorm(n=length(tmp),mean=2,sd=1,0,Inf)
     tmp=grep("b",phi.names)
     phi[i,tmp,1]=rtnorm(n=length(tmp),mean=1,sd=.5,0,Inf)
     tmp=grep("t0",phi.names)
     phi[i,tmp,1]=rtnorm(n=length(tmp),mean=1,sd=.5,0,Inf)
     tmp=grep("sve",phi.names)
     phi[i,tmp,1]=rtnorm(n=length(tmp),mean=.8,sd=.4,0,Inf)
     tmp=grep("svc",phi.names)
     phi[i,tmp,1]=rtnorm(n=length(tmp),mean=.8,sd=.4,0,Inf)
     tmp=grep("c",phi.names)
     phi[i,tmp,1]=rtnorm(n=length(tmp),mean=1,sd=1,0,Inf)
     tmp=grep("beta",phi.names)
     phi[i,tmp,1]=rtnorm(n=length(tmp),mean=1,sd=1,0,Inf)
}
prior=list()

tmp=grep("A",theta.names,value=TRUE)
for (n in 1:length(tmp)) {
     tmp2=tmp[n]
     prior[[tmp2]]=list(mu=c(1,1),sigma=c(1,1))
}

tmp=grep("b",theta.names,value=TRUE)
for (n in 1:length(tmp)) {
     tmp2=tmp[n]
     prior[[tmp2]]=list(mu=c(0.4,0.4),sigma=c(0.4,0.4))
}


tmp=grep("t0",theta.names,value=TRUE)
for (n in 1:length(tmp)) {
     tmp2=tmp[n]
     prior[[tmp2]]=list(mu=c(.3,.3),sigma=c(.3,.3))
}

tmp=grep("c",theta.names,value=TRUE)
for (n in 1:length(tmp)) {
     tmp2=tmp[n]
     prior[[tmp2]]=list(mu=c(1,1),sigma=c(1,1))
}


tmp=grep("beta",theta.names,value=TRUE)
for (n in 1:length(tmp)) {
     tmp2=tmp[n]
     prior[[tmp2]]=list(mu=c(1,1),sigma=c(1,1))
}


tmp=grep("sve",theta.names,value=TRUE)
for (n in 1:length(tmp)) {
     tmp2=tmp[n]
     prior[[tmp2]]=list(mu=c(1,1),sigma=c(1,1))
}

tmp=grep("svc",theta.names,value=TRUE)
for (n in 1:length(tmp)) {
     tmp2=tmp[n]
     prior[[tmp2]]=list(mu=c(1,1),sigma=c(1,1))
}


begin = date()  # useful for checking how long your fit runs

for(i in 1:n.chains){
     for(j in 1:S){
          while (weight[1,i,j]==-Inf) {
               # Change accordingly, and make sure the mean/sd vectors are n.pars long.
               theta[i,,j,1]=rtnorm(n=n.pars,mean=start.points,sd=start.points/5,0,Inf)
               weight[1,i,j]=log.dens.like(theta[i,,j,1],data=data[[j]],par.names=theta.names)
          }
     }
     
}

for(i in 2:nmc){
     if(i%%10==0)cat("\n ",i,'/',nmc,' ')
     phi[,,i]=phi[,,i-1]
     rand.samp=sample(1:n.chains,n.chains)
     for (p in theta.names) {
          which.theta=match(x=p,table=theta.names)
          which.phi=match(x=paste(p,c("mu","sigma"),sep="."),table=phi.names)
          if (i %% migration.freq == 0 & i > migration.start & i < migration.end) {
               phi[,,i]=migration.crossover_hyper(pars=which.phi,use.theta=theta[rand.samp,which.theta,,i-1],use.phi=phi[,,i],prior=prior[[p]])
          } else {
               phi[,,i]=t(sapply(1:n.chains,crossover_hyper,pars=which.phi,use.theta=theta[rand.samp,which.theta,,i-1],use.phi=phi[,,i],prior=prior[[p]]))
          }
     }
     rand.samp=sample(1:n.chains,n.chains)
     hyper=phi[rand.samp,,i]
     for(j in 1:S){
          if (i %% migration.freq == 0 & i > migration.start & i < migration.end) {
               temp = migration.crossover(pars=1:n.pars,use.theta=theta[,,j,i-1],use.like=weight[i-1,,j],data=data[[j]],hyper=hyper,par.names=theta.names)
          } else {
               temp = t(sapply(1:n.chains,crossover,pars=1:n.pars,use.theta=theta[,,j,i-1],use.like=weight[i-1,,j],data=data[[j]],
                               hyper=hyper,par.names=theta.names))
          }
          
          weight[i,,j]=temp[,1]
          theta[,,j,i]=temp[,2:(n.pars+1)]
          
     }
     
}


end = date()
begin
end			
nmc_thin = seq(1,nmc,by = thin_steps)
theta = theta[,,,nmc_thin]
phi = phi[,,nmc_thin]

samples = list(theta=theta,phi=phi,begin=begin,end=end,opt=opt)
save(samples,file=opt$file_out)


