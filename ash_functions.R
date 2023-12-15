### The Functions Used in MATH5472 Project

## Choose component variance from the fixed grid
get_sigmas=function(est_beta,est_s){
# est_beta: estimated beta
# est_s: estimated standard error
  sigma_min=min(est_s)/10   
  sigma_max=2*sqrt(max(est_beta^2-est_s^2))   
  sigmas=sigma_min
  while(sqrt(2)*max(sigmas)<sigma_max){
    sigmas=c(sigmas,sqrt(2)*max(sigmas))
  }
  return(sigmas)
}

# EM algorithm for parameter estimation
EM=function(J,K,est_beta,est_s){
# J: dimension
# K: number of components
# est_beta: estimated beta
# est_s: estimated standard error
  # initialization
  iter_max=100   # max iterations
  tol=10^(-6)    # tolerance
  pi=c(1-K/J,rep(1/J,K))
  lambda=c(10,rep(1,K))  # penalty coefficient
  l=matrix(0,K+1,J)
  w=matrix(0,K+1,J)
  for(j in 1:J)
    for(k in 0:K){
      if(k==0)
        l[k+1,j]=pnorm(est_beta[j],0,est_s[j])
      else
        l[k+1,j]=pnorm(est_beta[j],0,est_s[j]+sigmas[k])
    }
  for(iter in 1:iter_max){
    temp=matrix(rep(pi,J),K+1,J)*l
    w=temp/matrix(rep(colSums(temp),K+1),K+1,J,byrow=T)
    # E-step
    n=rowSums(w)+lambda-1
    # M-step
    new_pi=n/sum(n)
    if(sum((new_pi-pi)^2)<tol)
      break
    else
      pi=new_pi # update
  }
  return(list(pi=pi,w=w))
}

# Posterior distribution of effects
poster=function(beta,j,est_beta,est_s,est_w){
# beta: the beta value
# j: j-th dimension
# est_beta: estimated beta
# est_s: estimated standard error
# est_w: estimated weight matrix
  mix_comp=rep(0,K+1)
  mix_comp[1]=1*(beta>0)
  for(k in 1:K){
    mix_comp[k+1]=pnorm(est_beta[k],est_s[k])
  }
  poster=sum(est_w[,j]*mix_comp)
  return(poster)
}

## Samples for plot
get_samples=function(est_beta,est_s){
# est_beta: estimated beta
# est_s: estimated standard error
  z=est_beta/est_s                # z-scores
  p=2*pmin(pnorm(z),1-pnorm(z))   # p-values
  pp=runif(J/2)
  zz=sample(z,J/2)
  return(list(p=p,pp=pp,z=z,zz=zz))
}

# Local false discovery rate
get_lfdr=function(j,est_beta,est_s,est_w){
# j: j-th dimension
# est_beta: estimated beta
# est_s: estimated standard error
# est_w: estimated weight matrix
  pr1=poster(0.001,j,est_beta,est_s,est_w)
  pr2=poster(-0.001,j,est_beta,est_s,est_w)
  return(pr1-pr2)
}

# Local false sign rate
get_lfsr=function(j,est_beta,est_s,est_w){
# j: j-th dimension
# est_beta: estimated beta
# est_s: estimated standard error
# est_w: estimated weight matrix
  pr=poster(0,j,est_beta,est_s,est_w)
  return(max(pr,1-pr))
}

## ash results
ash_results=function(J,est_beta,est_s){
# J: dimension
# est_beta: estimated beta
# est_s: estimated standard error
  lfsr=rep(0,J)                      # local false sign rate
  sigmas=get_sigmas(est_beta,est_s)  # sigmas from grid
  K=length(sigmas)                   # number of components
  EM_result=EM(J,K,est_beta,est_s)   # parameter estimation by EM algorithm
  est_pi=EM_result$pi
  est_w=EM_result$w
  for(j in 1:J)
    lfsr[j]=get_lfsr(j,est_beta,est_s,est_w)
  return(lfsr)
}
