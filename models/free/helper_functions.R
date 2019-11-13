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
     N=length(drift) # Number of responses.
     if (N>2) {
          tmp=array(dim=c(length(t),N-1))
          for (i in 2:N) tmp[,i-1]=fptcdf(z=t,x0max=x0max[i],chi=chi[i],
                                          driftrate=drift[i],sddrift=sdI[i])
          G=apply(1-tmp,1,prod)
     } else {
          G=1-fptcdf(z=t,x0max=x0max[2],chi=chi[2],driftrate=drift[2],sddrift=sdI[2])
     }
     out=G*fptpdf(z=t,x0max=x0max[1],chi=chi[1],driftrate=drift[1],sddrift=sdI[1])
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
     out[crct] =n1PDFfixedt0(rt[crct]-t0,x0max=c(A,A),chi=c(b,b),drift=v,sdI=s)
     out[!crct]=n1PDFfixedt0(rt[!crct]-t0,x0max=c(A,A),chi=c(b,b),drift=v[2:1],sdI=s[2:1])
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
               vs=c(x[paste("vc",cond,sep=".")],1-x[paste("vc",cond,sep=".")])
          }
          if (doesStDevErrorVary) {
               s=c(1,x[paste("sve",cond,sep=".")])
          }
          
          s = c(x[paste('svc')],x[paste('sve')])
          tmp=data$cond==cond
          tmp=get.dens.2choice(rt=data$rt[tmp],crct=data$acc[tmp]==1,b=b,A=A,v=vs,s=s,t0=t0)
          
          out=out+sum(log(tmp))
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
