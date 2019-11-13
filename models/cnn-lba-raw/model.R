rm(list=ls())     

setwd('~/lba-cnn')
source('models/cnn-lba-raw/helper_functions.R')

library(msm)

##### Parameters #####

network = 'inception' #options are: resnet, inception, vgg
num_samples = 3000
category = 1 #options are: (1) greebles 1, (2) greebles 2, (3) zig 1, (4) zig 2, (5) sheingbugs

######################



opt = list()
opt$nmc = num_samples
opt$file_out = file.path('output',network,'raw',paste0('samples_cat',category,'.RData'))
opt$file_in = file.path('data',network,'raw','dat.RData')
opt$category = category

cat('\n')
cat('Samples will be save to ', opt$file_out)


dat = load_file(opt$file_in)
dst_euc_mat = dat$dst
dat = dat$dat

dat = dat[dat$trial > 180, ]
dat = dat[dat$response != 'timeout',]

cat_list = list(c(1,100),c(101,200),c(201,300),c(301,400),c(401,500))

cat = opt$cat
dat = dat[dat$ID_1 >= cat_list[[cat]][1]  & dat$ID_1 <= cat_list[[cat]][2],]

data = lapply(unique(dat$subject), 
             function(x){
                  data.frame(subject = dat$subject[dat$subject==x], 
                             rt=dat$rt[dat$subject==x] / 1000, 
                             response = dat$response[dat$subject==x],
                             dst = dst_euc_mat[cbind(dat$ID_1[dat$subject==x], dat$ID_2[dat$subject==x])]
                  )
             })


doesStartPointVary = 0    # Does start point vary over conditions for this model?
doesThresholdVary = 0    # Does threshold vary over conditions for this model? (NOTE: This MUST also vary if start point varies)
doesDriftRateCorrectVary = 0     # Does the mean drift rate for the correct accumulator vary over conditions for this model?
doesDriftRateErrorVary = 0   # Does the mean drift rate for the error accumulator vary over conditions for this model?
doesStDevErrorVary = 0    # Does the standard deviation of the drift rate of the error accumulator vary over conditions for this model?
doesT0Vary = 0    # Does t0 vary over conditions for this model?

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

if (!dir.exists(dirname(opt$file_out)))
{
        dir.create(dirname(opt$file_out),recursive = TRUE)
}

save(samples,file=opt$file_out)


