library(nimble)
library(lattice)
library(mvtnorm)
library(matrixStats)
library(fields)
library(boot)
library(maps)
library(mapproj)

##### DAGAR precision #####
dagarcov <- nimbleFunction(run =function(rho=double(0),Minc=double(2),ni=double(1),maxn=double(0),
  neimat=double(2), nijvec=double(1), intersectmat=double(2)) {returnType(double(2))
    n=length(ni)
    u=rho^2
    sumfuncvec=nimArray(0,maxn)
    for(i in 1:maxn)  sumfuncvec[i]=i/(1-u+i*u)
    cumsumfuncvec=nimArray(0,maxn)
    cumsumfuncvec[1]=sumfuncvec[1]
    for(i in 2:maxn) cumsumfuncvec[i]=cumsumfuncvec[i-1]+sumfuncvec[i]
    Qd=matrix(0,n,n)
    for(i in 1:n){
        s=0
        if(ni[i]>0) for(ck in 1:ni[i]){
            k=neimat[i,ck]
            s=s+ u*cumsumfuncvec[ni[k]]/(ni[k]*(ni[k]+1))
            Qd[i,k]=Qd[i,k]-rho
        }
        Qd[i,i]=(1-u)+u*ni[i]/2+s
        if(i < n){
            for(j in (i+1):n){
                t=0
                jc=0
                counter=(i-1)*n+j
                nij=nijvec[counter]
                if(nij>0){
                    jointn=intersectmat[counter,1:nij]
                    for(ck in 1:length(jointn))
                    {
                        k=jointn[ck]
                        t=t+1/(2*(ni[k]+1))+(ni[k]-cumsumfuncvec[ni[k]])/(ni[k]*(ni[k]+1)*(ni[k]-1))
                    }
                }
                Qd[i,j]=Qd[i,j]+t
                Qd[j,i]=Qd[j,i]+t
            }
        }
    }
    Qd=Qd/(1-rho^2)
    return(Qd)})

##### as.vector ####
nimbleVector <- nimbleFunction(run=function(vM=double(2),zero=double(1)){
    returnType(double(1)) 
    v=zero
    for(i in 1:length(v)) v[i]=vM[i,1]
    return(vM)
})

#### model blocks #####
### order free model ##
dagarCode <- nimbleCode({
    G[,] <- dagarcov(rho,Minc[,],ni[], maxn, neimat[,], nijvec[], intersectmat[,])
    Vw[,] <- G[,]*invsigmasq
    w[] ~ dmnorm(muw[], prec=Vw[,])
    Vy[,] <- id[,]*invtausq
    mu[] <- X[,]%*%beta[]+w[]
    y[] ~ dmnorm(mu[], prec=Vy[,])
    beta[] ~ dmnorm(mubeta[],cholesky=Cbeta[,], prec_param=0)
    invsigmasq ~ dgamma(2, rate=1)
    invtausq ~ dgamma(2, rate=0.2)
    rho ~dunif(0,1)
})

### ordered model ####
onedagarCode <- nimbleCode({
    w[1] ~ dnorm(0, invsigmasq)
    y[1] ~ dnorm(X[1,1]*beta[1]+X[1,2]*beta[2] + w[1], invtausq)
    for (i in 2:n){
        weight <- rho/(1+(dni[i-1]-1)*rho^2)
        
        mw[i-1] <- get_mean(udnei=udnei[],d=dni[i-1],c=cni[i-1],wv=w[1:(i-1)])
        w[i] ~ dnorm(mw[i-1]*weight, invsigmasq*(1+(dni[i-1]-1)*rho^2)/(1-rho^2))
        
        y[i] ~ dnorm(X[i,1]*beta[1]+X[i,2]*beta[2] + w[i], invtausq)
    }
    beta[1] ~ dnorm(0,sd=100)
    beta[2] ~ dnorm(0,sd=100)
    invsigmasq ~ dgamma(2, rate=1)
    invtausq ~ dgamma(2, rate=0.2)
    rho ~dunif(0,1)
})

get_mean <- nimbleFunction(
    run = function(udnei = double(1), d = double(0), c=double(0), wv = double(1)) {
        mw <- 0
        for(i in 1:d){
            ind <- udnei[c-d+i]
            mw <- mw + wv[ind]
        }
        
        returnType(double(0))
        return(mw)
    }
)


### creates the adjacency matrix for grids ###
inclattice=function(m){
    n=m^2
    Minc=matrix(0,n,n)
    for(i in 1:(m-1))	for(j in 1:(m-1)) Minc[(i-1)*m+j,(i-1)*m+j+1]=Minc[(i-1)*m+j,i*m+j]=1
    for(i in 1:(m-1)) Minc[(i-1)*m+m,i*m+m]=1
    for(j in 1:(m-1)) Minc[(m-1)*m+j,(m-1)*m+j+1]=1
    Minc+t(Minc)
}

######## Basic settings ########
Nmcmc=10000
Nburn=5001

rho=0.5
r=0.1
seed=1

m=10
Minc=inclattice(m)


### defining some constants for the DAGAR model
n=nrow(Minc)
ni=rowSums(Minc)
maxn=max(ni)
neimat=matrix(0,n,maxn)
neighbors=lapply(1:n,function(x) which(Minc[x,]==1))
dneighbors=sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
dni=sapply(dneighbors,length)
nmax=max(dni)
cni=cumsum(dni)
dneimat=sapply(dneighbors, function(nei,nmax,n) c(nei,rep(n+1,nmax+1-length(nei))),nmax,n)
udnei=unlist(dneighbors)
for(i in 1:n) neimat[i,1:ni[i]]=neighbors[[i]]
intersectmat=matrix(0,n^2,maxn)
nijvec=rep(0,n^2)
for(i in 1:n) for(j in 1:n)
{ neighij=intersect(neighbors[[i]],neighbors[[j]])
nij=length(neighij)
nijvec[n*(i-1)+j]=nij
if(nij >0) intersectmat[n*(i-1)+j,1:nij]=neighij
}


##### data generation #####
muw <- rep(0,n)
id <- diag(n)

sigs=4
set.seed(seed)

s=cbind(rep(1:m,m),kronecker(1:m,rep(1,m)))
dmat=as.matrix(dist(s))
G=rho^dmat
w=as.vector(rmvnorm(1,rep(0,n),sigs*G))

X=matrix(rnorm(2*n),ncol=2)
beta=c(1,5)
Y=as.vector(X%*%beta)+w+sqrt(r*sigs)*rnorm(n)

##### ordered dagar model #####
set.seed(1)
onedagarModel <- nimbleModel(onedagarCode, constants = list(n = n, dmat=dmat, muw=muw, id = id, ni=ni, Minc=Minc, maxn=maxn, neimat=neimat, mubeta=c(0,0), 
    Cbeta=100*identityMatrix(2),nijvec=nijvec,intersectmat=intersectmat,zero=rep(0,n),dneighbors=dneighbors,udnei=udnei,dni=dni,p=ncol(X),cni=cni,dneimat=dneimat,
    nmax=nmax),dimensions=list(y=n, X=c(n,2), mu=n, Vw=c(n, n), w=n, Vy=c(n, n), G = c(n, n), beta=2,mubeta=2, Cbeta=c(2,2), M=c(n+1,n)), check=FALSE)

onedagarModel$setData(list(y = Y, X=X))

onedagarModel$invsigmasq <- 1
onedagarModel$invtausq <- 0.1
onedagarModel$rho <- 0.5
onedagarModel$w <- rep(0,n)

ConedagarModel <- compileNimble(onedagarModel)

onedagarconf <- configureMCMC(ConedagarModel, print=TRUE)
onedagarconf$printSamplers()
onedagarconf$addMonitors(c('w'))	

onedagarMCMC <- buildMCMC(onedagarconf)

ConedagarMCMC <- compileNimble(onedagarMCMC, project = ConedagarModel)

tod1=Sys.time()
ConedagarMCMC$run(Nmcmc)
tod2=Sys.time()

### MCMC output 
onedagarMCMCsamples <- as.matrix(ConedagarMCMC$mvSamples)

### summary statistics 
rhohatod=mean(onedagarMCMCsamples[Nburn:Nmcmc,'rho'])
rhatod=mean(onedagarMCMCsamples[Nburn:Nmcmc,'invsigmasq']/onedagarMCMCsamples[Nburn:Nmcmc,'invtausq'])
b1hatod=mean(onedagarMCMCsamples[Nburn:Nmcmc,'beta[1]'])
b2hatod=mean(onedagarMCMCsamples[Nburn:Nmcmc,'beta[2]'])

wpostod=onedagarMCMCsamples[,paste0("w[",1:n,"]")]
whatod=apply(wpostod[Nburn:Nmcmc,],2,mean)
msevecod=mean((w-whatod)^2)


##### dagar model #####
set.seed(1)
dagarModel <- nimbleModel(dagarCode, constants = list(n = n, dmat=dmat, muw=muw, id = id, ni=ni, Minc=Minc, maxn=maxn, neimat=neimat, mubeta=c(0,0), 
    Cbeta=100*identityMatrix(2),nijvec=nijvec,intersectmat=intersectmat,zero=rep(0,n),p=ncol(X)), dimensions=list(y=n, X=c(n,2), mu=n, Vw=c(n, n), w=n, 
    Vy=c(n, n), G = c(n, n), beta=2,mubeta=2, Cbeta=c(2,2)), check=FALSE)

dagarModel$setData(list(y = Y,X=X))

dagarModel$invsigmasq <- 1
dagarModel$invtausq <- 0.1
dagarModel$rho <- 0.5
dagarModel$w <- rep(0,n)

CdagarModel <- compileNimble(dagarModel)

dagarconf <- configureMCMC(CdagarModel, print=TRUE)
dagarconf$printSamplers()
dagarconf$addMonitors(c('w'))	

dagarMCMC <- buildMCMC(dagarconf)

CdagarMCMC <- compileNimble(dagarMCMC, project = CdagarModel)

td1=Sys.time()
CdagarMCMC$run(Nmcmc)
td2=Sys.time()

### MCMC output 
dagarMCMCsamples <- as.matrix(CdagarMCMC$mvSamples)

### summary statistics
rhohatd=mean(dagarMCMCsamples[Nburn:Nmcmc,'rho'])
rhatd=mean(dagarMCMCsamples[Nburn:Nmcmc,'invsigmasq']/dagarMCMCsamples[Nburn:Nmcmc,'invtausq'])
b1hatd=mean(dagarMCMCsamples[Nburn:Nmcmc,'beta[1]'])
b2hatd=mean(dagarMCMCsamples[Nburn:Nmcmc,'beta[2]'])

wpostd=dagarMCMCsamples[,paste0("w[",1:n,"]")]
whatd=apply(wpostd[Nburn:Nmcmc,],2,mean)
msevecd=mean((w-whatd)^2)
