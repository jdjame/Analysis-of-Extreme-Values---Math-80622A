#############################################
# Jean-Yves Djamen                          #
# Math 80622A Final project code            #
#############################################

# Sourcing data for analysis 

library(evir)
data(nidd.thresh)
nidd.data<-nidd.thresh
detach(package:evir)

library(texmex)
pharma.data<-liver[c('TBL.B','TBL.M','dose')]
detach(package:texmex)

library(mev)
library(ggplot2)
library(quantreg)
library(gridExtra)


########################################################  Simulation Study########################################################  

run.mini.exp<-function(exp.data,us,n,i,gp.rmse.T1,gp.rmse.T2,egp1.rmse.T1,egp1.rmse.T2,egp2.rmse.T1,egp2.rmse.T2,egp1.kap, egp2.kap){
  for (u in 1:length(us)){
    #for each u calculate return level
    temp.data<-exp.data[exp.data>us[u]]
    
    print('---------------gp---------------')
    #quantile based on gp
    
    gp.out<-extRemes::fevd(x=temp.data, threshold = us[u], type='GP',method='MLE', verbose=F)
    gp.rmse.T1[i,u]<-extRemes::return.level(gp.out, return.period=n*1.5)
    gp.rmse.T2[i,u]<-extRemes::return.level(gp.out, return.period=n*5)
    #gpd.mle(temp.data-us[u],args='quant', p=1/(n*1.5))
    #gpd.mle(temp.data-us[u],args='quant', p=1/(n*5))
    
    
    #quantile based on egp1 
    error<- F
    tryCatch(
      {
        print('---------------egp1---------------')
        egp1.out<-fit.egp(temp.data, us[u], model='egp1',show=F)#,init=c(1,as.vector(gp.out$estimate)))
        
        print('---------------egp2---------------')
        #quantile based on egp2
        egp2.out<-fit.egp(temp.data, us[u], model='egp2',show=F)#,init=c(1,as.vector(gp.out$estimate)))
        
      },
      error=function(e){
        #if there is an error terminate
        error<<-T
      }
    )
    
    if(! error ){
      egp1.kap[i,u]<-egp1.out$estimate[1]
      egp2.kap[i,u]<-egp2.out$estimate[1]
      egp1.rmse.T1[i,u]<-egp.retlev(temp.data,us[u],as.vector(egp1.out$estimate),model='egp1',p=1/(n*1.5),plot=FALSE)[1]
      egp1.rmse.T2[i,u]<-egp.retlev(temp.data,us[u],as.vector(egp1.out$estimate),model='egp1',p=1/(n*5),plot=FALSE)[1]
      
      egp2.rmse.T1[i,u]<-egp.retlev(temp.data,us[u],as.vector(egp2.out$estimate),model='egp2',p=1/(n*1.5),plot=FALSE)[1]
      egp2.rmse.T2[i,u]<-egp.retlev(temp.data,us[u],as.vector(egp2.out$estimate),model='egp2',p=1/(n*5),plot=FALSE)[1]
    }
    else{
      return(list('passed'=F))
    }
  }
  
  return(list('gp.rmse.T1'=gp.rmse.T1,'gp.rmse.T2'=gp.rmse.T2,
              'egp1.rmse.T1'=egp1.rmse.T1,'egp1.rmse.T2'=egp1.rmse.T2,
              'egp2.rmse.T1'=egp2.rmse.T1,'egp2.rmse.T2'=egp2.rmse.T2,'passed'=T,
              'egp1.kap'=egp1.kap, 'egp2.kap'=egp2.kap))
}


run.exp<-function(n, num.exp=4, gen.method=1,opt.seed=80622,num.pts=20,range.flag=F,range){
  
  #############################################
  #n: number of samples                       #
  #num.exp: number of experiments to run      #
  #gen.method: method used to generate data   #
  #               1->Normal(0,1)              #
  #               2->Student t2               #
  #               1->Beta(1.5,3)              #
  #############################################
  
  set.seed(opt.seed)
  
  #sequence of quantiles to use as thresholds
  quants<-c(1/n, 1-30/n)
  if(gen.method==1){
    #thresholds
    us<-qnorm(quants)
    
    #true values of z 
    z.t1= qnorm(1-1/(n*1.5))
    z.t2= qnorm(1-1/(n*5))
  }
  else if(gen.method==2){
    #thresholds
    us<-qt(quants,df=2)
    
    #true values of z 
    z.t1= qt(1-1/(n*1.5),df=2)
    z.t2= qt(1-1/(n*5),df=2)
  }
  else{
    #thresholds
    us<-qbeta(quants,shape1=1.5, shape2=3)
    
    #true values of z 
    z.t1= qbeta(1-1/(n*1.5),shape1=1.5, shape2=3)
    z.t2= qbeta(1-1/(n*5),shape1=1.5, shape2=3)
  }
  
  us<-seq(us[1],us[2],(us[2]-us[1])/num.pts)
  
  if(range.flag){
    us=seq(range[1], range[2], (range[2]-range[1])/num.pts)
  }
  
  #initialize matrices to hold estimates for return levels based on threshold
  gp.rmse.T1= matrix(nrow=num.exp, ncol=length(us))
  egp1.rmse.T1= matrix(nrow=num.exp, ncol=length(us))
  egp2.rmse.T1= matrix(nrow=num.exp, ncol=length(us))
  
  gp.rmse.T2= matrix(nrow=num.exp, ncol=length(us))
  egp1.rmse.T2= matrix(nrow=num.exp, ncol=length(us))
  egp2.rmse.T2= matrix(nrow=num.exp, ncol=length(us))
  
  #initialize matrix of kappa parameters
  egp1.kap<-matrix(nrow=num.exp, ncol=length(us))
  egp2.kap<-matrix(nrow=num.exp, ncol=length(us))
  
  runs<-1
  #for(i in 1:num.exp){
  while(runs<=num.exp){
    #for each experiment, generate new data
    if(gen.method==1){
      exp.data<-rnorm(n)
    }
    else if(gen.method==2){
      exp.data<-rt(n,df=2)
    }
    else{
      exp.data<-rbeta(n,shape1=1.5, shape2=3)
    }
    exp.data<-as.vector(exp.data)
    
    temp.out=run.mini.exp(exp.data,us,n,runs,
                          gp.rmse.T1,gp.rmse.T2,
                          egp1.rmse.T1,egp1.rmse.T2,
                          egp2.rmse.T1,egp2.rmse.T2,egp1.kap, egp2.kap)
    
    if(temp.out$passed){
      gp.rmse.T1=temp.out$gp.rmse.T1
      gp.rmse.T2=temp.out$gp.rmse.T2
      egp1.rmse.T1=temp.out$egp1.rmse.T1
      egp1.rmse.T2=temp.out$egp1.rmse.T2
      egp2.rmse.T1=temp.out$egp2.rmse.T1
      egp2.rmse.T2=temp.out$egp2.rmse.T2
      egp1.kap=temp.out$egp1.kap
      egp2.kap=temp.out$egp2.kap
      runs=runs+1
      print(runs)
    }
    else{
      print("failed, trying again")
    }
  }
  
  return(list('gp.rmse.T1'= gp.rmse.T1,'gp.rmse.T2'= gp.rmse.T2,
              'egp1.rmse.T1'= egp1.rmse.T1,'egp1.rmse.T2'= egp1.rmse.T2,
              'egp2.rmse.T1'= egp2.rmse.T1,'egp2.rmse.T2'= egp2.rmse.T2,
              'z.t1'=z.t1,'z.t2'=z.t2,'us'=us,'egp1.kap'=egp1.kap,'egp2.kap'=egp2.kap))
  
}





####### Normal Data #######
norm.100<-run.exp(100,num.exp=250,opt.seed=8,num.pts=35,range.flag=T,range=c(-2.5,0.5))
norm.1000<-run.exp(1000,num.exp=250,opt.seed=0,num.pts=50,range.flag=T,range=c(-3,2))
norm.10000<-run.exp(10000,num.exp=250,opt.seed=6,num.pts=65,range.flag=T,range=c(-3.5,3))

####### Student T data ####### 
st.100<-run.exp(100,num.exp=250,gen.method = 2,opt.seed=22,num.pts=20,range.flag=T,range=c(-2,0.5))
st.1000<-run.exp(1000,num.exp=250,gen.method = 2,opt.seed=80,num.pts=6,range.flag=T,range=c(-0.5,4))
#st.10000<-run.exp(10000,num.exp=100,gen.method = 2,opt.seed=100)

####### Beta data ####### 
beta.100<-run.exp(n=100,num.exp=250,gen.method = 3,opt.seed=62,num.pts=10,range.flag=T,range=c(0,0.5))
beta.1000<-run.exp(n=1000,num.exp=250,gen.method = 3,opt.seed=2,num.pts=20,range.flag=T,range=c(0,0.7))
#beta.10000<-run.exp(10000,num.exp=100,gen.method = 3,opt.seed=100)


get.exp.plots<-function(exp.output, ret.lev=1,y.lim=c(0,0.5),x.lab="",y.lab="",titl=""){
  
  if(ret.lev==1){
    gp.rmse<-colMeans((sqrt((exp.output$gp.rmse.T1-exp.output$z.t1)^2)))
    egp1.rmse<-colMeans((sqrt((exp.output$egp1.rmse.T1-exp.output$z.t1)^2)))
    egp2.rmse<-colMeans((sqrt((exp.output$egp2.rmse.T1-exp.output$z.t1)^2)))
  }
  else{
    gp.rmse<-colMeans((sqrt((exp.output$gp.rmse.T2-exp.output$z.t2)^2)))
    egp1.rmse<-colMeans((sqrt((exp.output$egp1.rmse.T2-exp.output$z.t2)^2)))
    egp2.rmse<-colMeans((sqrt((exp.output$egp2.rmse.T2-exp.output$z.t2)^2)))

  }
  #print(egp1.rmse)
  plot(exp.output$us,gp.rmse, lty=1,type='l',col='grey' ,ylim=y.lim,xlab=x.lab, ylab=y.lab,main=titl,lwd=2)
  lines(exp.output$us,egp1.rmse,  lty=1,lwd=2)
  lines(exp.output$us, egp2.rmse, lty=2,lwd=2)
}

par(mfrow=c(2,3))
#locator(2)
get.exp.plots(norm.100,titl='100', y.lab=expression('T'[1]),y.lim=c(0,1))
get.exp.plots(norm.1000,titl='1000',y.lim=c(0,1))
get.exp.plots(norm.10000,titl='10000',y.lim=c(0,1))
get.exp.plots(norm.100,2, y.lab=expression('T'[2]),y.lim=c(0,1))
get.exp.plots(norm.1000,2,x.lab='Threshold',y.lim=c(0,1))
get.exp.plots(norm.10000,2,y.lim=c(0,1))


par(mfrow=c(2,4))
get.exp.plots(beta.100,y.lim=c(0.03,0.086),titl=expression('n=100,T'[1])) 
get.exp.plots(beta.100,2,y.lim=c(0,0.11),titl=expression('n=100,T'[2]))
get.exp.plots(beta.1000,y.lim=c(0,0.05),titl=expression('n=1000,T'[1]))
get.exp.plots(beta.1000,2,y.lim=c(0,0.042),titl=expression('n=1000,T'[2]))

get.exp.plots(st.100,y.lim=c(0,15),titl=expression('n=100,T'[1]))
get.exp.plots(st.100,2,y.lim=c(7,18),titl=expression('n=100,T'[2]))
get.exp.plots(st.1000,y.lim=c(5,30),titl=expression('n=1000,T'[1]))
get.exp.plots(st.1000,2,y.lim=c(0,200),titl=expression('n=1000,T'[2]))


library(matrixStats)

get.kap.plots<-function(inpt,num.pts,range,colstop){
  
  in.data<-inpt$egp1.kap[,1:colstop]
  ks<-colMedians(in.data)
  
  #errs<-apply(inpt$egp1.kap, 2, sd)
  
  lb<-c()#ks-1.96*errs
  ub<-c()#quantile(inpt$egp1.kap, 0.975)#ks+1.96*errs
  
  for(i in 1:dim(in.data)[2]){
    lb[i]<-quantile(in.data[,i], 0.025)
    ub[i]<-quantile(in.data[,i], 0.975)
  }
  
  us<-seq(range[1],range[2],(range[2]-range[1])/num.pts)
  
  plot1<-ggplot()+
    geom_point(mapping = aes(x = us,y=ks))+
    geom_ribbon(aes(x=us,ymin= lb, ymax=ub),fill='grey50',alpha=0.5,linetype=2, colour='red')+
    geom_hline(yintercept=1)
    
  
  return(plot1)
}
kap1<-get.kap.plots(norm.100,num.pts=35,range=c(-2.5,0.5),colstop=36)+xlab("")+ylab(expression(hat(kappa)))+ggtitle('100')+
  theme(plot.title = element_text(hjust = 0.5))
kap2<-get.kap.plots(norm.1000,num.pts=49,range=c(-3,2),colstop=50)+xlab("Threshold")+ylab(expression(hat(kappa)))+ggtitle('1000')+#+ylim(c(0,5.5))
  theme(plot.title = element_text(hjust = 0.5))
kap3<-get.kap.plots(norm.10000,num.pts=63,range=c(-3.5,2),colstop=64)+xlab("")+ylab(expression(hat(kappa)))+ggtitle('10000')+#ylim(c(0,8))+
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(kap1,kap2,kap3, nrow=1, ncol=3)


########################################################  River Nidd Data########################################################  

############# Parammeter Plots #############
#list of thresholds to use

us<-seq(65.08,88.61,(88.61-(65.08))/39)
#initialize list of items to plot
gp.shape<- gp.scale<- c()
gp.shape.ci<-gp.scale.ci<- matrix(NA, nrow=length(us), ncol = 2)

egp.scale<-egp.shape<-egp.kap<-c()
egp.shape.ci<-egp.scale.ci<-egp.kap.ci<- matrix(NA, nrow=length(us), ncol = 2)

#for each threshold, fit gp and EGP1, the find parameter estimates and their CI (asymptotic normality)
for (i in 1:length(us)){
  
  # Standard GP params 
  gp.out<-fit.gpd(nidd.data, us[i])
  gp.scale<-c(gp.scale,gp.out$estimate[[1]] )
  gp.shape<-c(gp.shape,gp.out$estimate[[2]] )
  
  gp.ci<-confint(gp.out)
  gp.scale.ci[i,]<- as.vector(gp.ci[1,])#c(gp.out$estimate[[1]]-1.96*gp.out$std.err[[1]],gp.out$estimate[[1]]+1.96*gp.out$std.err[[1]])
  gp.shape.ci[i,]<- as.vector(gp.ci[2,])#c(gp.out$estimate[[2]]-1.96*gp.out$std.err[[2]],gp.out$estimate[[2]]+1.96*gp.out$std.err[[2]])
  
  # EGP1 params
  temp.data<-nidd.data[nidd.data>us[i]]
  egp.out<-fit.egp(temp.data, us[i], model='egp1')
  egp.kap<-c(egp.kap,egp.out$estimate[[1]] )
  egp.scale<-c(egp.scale,egp.out$estimate[[2]] )
  egp.shape<-c(egp.shape,egp.out$estimate[[3]] )
  egp.ci<-confint(egp.out,level=0.95)
  egp.kap.ci[i,]<- as.vector(egp.ci[1,])#c(egp.out$estimate[[1]]-1.96*egp.out$std.err[[1]],egp.out$estimate[[1]]+1.96*egp.out$std.err[[1]])
  egp.scale.ci[i,]<- as.vector(egp.ci[2,])#c(egp.out$estimate[[2]]-1.96*egp.out$std.err[[2]],egp.out$estimate[[2]]+1.96*egp.out$std.err[[2]])
  egp.shape.ci[i,]<- as.vector(egp.ci[3,])#c(egp.out$estimate[[3]]-1.96*egp.out$std.err[[3]],egp.out$estimate[[3]]+1.96*egp.out$std.err[[3]])
}

# modify the scale and its cis
egp.scale=egp.scale -egp.shape * as.numeric(us)
egp.scale.ci=egp.scale.ci - egp.shape * as.numeric(us)
gp.scale=gp.scale -gp.shape * as.numeric(us)
gp.scale.ci=gp.scale.ci - gp.shape * as.numeric(us)

# plot all graphs
par(mfrow=c(3,2))

plot.graph.exp2<-function(estimates, estimates.ci){
  plot1<-ggplot()+
    geom_point(mapping = aes(x = us,y=estimates))+
    geom_ribbon(aes(x=us,ymin= estimates.ci[,1], ymax=estimates.ci[,2]),fill='grey50',alpha=0.5,linetype=2, colour='red')
    
  return(plot1)
}

#plot egp1 shape
a<-ggplot()+ geom_point(mapping = aes(x = us,y=egp.shape))+
  geom_ribbon(aes(x=us,ymin= egp.shape.ci[,1], ymax=egp.shape.ci[,2]),fill='grey50',alpha=0.5,linetype=2, colour='red')+
  xlab('u')+ylab(expression(hat(xi)['EGP1']))
  
#plot egp1 scale
b<-ggplot()+ geom_point(mapping = aes(x = us,y=egp.scale))+ylim(-150,150)+
  geom_ribbon(aes(x=us,ymin= egp.scale.ci[,1], ymax=egp.scale.ci[,2]),fill='grey50',alpha=0.5,linetype=2, colour='red')+
  xlab('u')+ylab(expression(hat(sigma)['EGP1']))

c<-plot.graph.exp2(egp.kap,egp.kap.ci)+geom_hline(yintercept=1.)+
  xlab('u')+ylab(expression(hat(kappa)['EGP1']))


d<-plot.graph.exp2(gp.shape,gp.shape.ci)+ylim(-0.5,1.50)+
  xlab('u')+ylab(expression(hat(xi)['GP']))
e<-plot.graph.exp2(gp.scale,gp.scale.ci)+ylim(-100,70)+
  xlab('u')+ylab(expression(hat(sigma)['GP']))


# histogram with threshold excedences
egp1.out<-fit.egp(nidd.data, 65.08, model='egp1')
points.y<-c()
points.x<-seq(65,340, 0.25)
for(i in 1:length(points.x)){
  points.y<-c(points.y,egp.ll(points.x[i], 65.08,as.vector(egp1.out$estimate), model='egp1' ))
}


density.plt<-ggplot()+
  geom_line(data=data.frame('x'=points.x,'y'=exp(points.y)) ,aes(x=x,y=y))+
  geom_histogram(data=as.data.frame(as.vector(nidd.data)),aes(x=nidd.data,y=..density..),color="black",alpha=0.2,bins=57)+
  xlim(c(65,340))+
  ylim(c(0,0.04))

grid.arrange(a,d,b,e,c,density.plt, ncol=2)

############# #diagnostic of proposed egp and proposed gp

#qq plots for egp 
obs.data<-nidd.data[nidd.data>65.08]
n<-length(obs.data)
ps<-seq(1/(n+1),n/(n+1), 1/(n+1) )
egpd.out<-fit.egp(nidd.data,65.08 , model='egp1',show=T)
egpd.ci<-confint(egpd.out)
egp.estim<-as.vector(egpd.out$estimate)
fitted.data<-c()
for(i in 1:n){
  fitted.data[i]<-egp.retlev(obs.data,65.08,egp.estim,model='egp1',p=1-ps[i],plot=F)
}
plot(fitted.data, sort(obs.data), xlab='Fitted EGP1', ylab='Observed')
abline( a=0, b=1)

#qq plots for gp
n<-length(obs.data)
ps<-seq(1/(n+1),n/(n+1), 1/(n+1) )
get.quants<-function(a, u,p){
  u + (a[1] * (p^(-a[2]) - 1))/a[2]
}
gp.out<-fit.gpd(nidd.data, 75)
plot(get.quants(c(gp.out$estimate[[1]],gp.out$estimate[[2]]),65.08, 1-ps ), sort(obs.data),xlab='Fitted GP', ylab='Observed')
abline( a=0, b=1)



############# Return level plots #############
us= c(seq(65.08,75.3,(75.3- 65.08)/4), seq(75.3, 82.6,(82.6-75.3)/2))
#us= c(seq(65.08,80.41,(75.3- 65.08)/4))
return.levels<- c(10,50,100,200,500)
gp.out.mat<-matrix(ncol=length(return.levels), nrow=length(us))
egp1.out.mat<-matrix(ncol=length(return.levels), nrow=length(us))
egp2.out.mat<-matrix(ncol=length(return.levels), nrow=length(us))

for (u in 1:length(us)){
  temp.data<-nidd.data[nidd.data>us[u]]
  
  #fit egps 
  egp1.out<-fit.egp(temp.data, us[u], model='egp1',show=F)
  egp2.out<-fit.egp(temp.data, us[u], model='egp2',show=F)
  
  #calculate each of the return levels
  for(retlev in 1:length(return.levels)){
    egp1.out.mat[u,retlev]<-egp.retlev(temp.data,us[u],as.vector(egp1.out$estimate),model='egp1',p=1/return.levels[retlev],plot=F)
    egp2.out.mat[u,retlev]<-egp.retlev(temp.data,us[u],as.vector(egp2.out$estimate),model='egp2',p=1/return.levels[retlev],plot=FALSE)
    gp.out.mat[u,retlev]<-gpd.mle(nidd.data-us[u],args='quant',p=1/return.levels[retlev])
  }
}


plot.retlev<-function(retlev.mat, ret.levs=return.levels){
  init.vec=retlev.mat[1,]
  plot(ret.levs, init.vec, ylim=c(0,1300), pch='1',type="b",xlim=c(10,500),log='x',
       ylab='Estimated Return Level',
       xlab='Return Period (years)')
  
  for (i in 2:(dim(retlev.mat)[1])){
    lines(ret.levs,retlev.mat[i,],pch=as.character(i),type="b")
  }
}

par(mfrow=c(1,3))
plot.retlev(egp1.out.mat)
plot.retlev(gp.out.mat)
plot.retlev(egp2.out.mat)



########################################################  Pharmaceutical data########################################################  
pharma.data.A<-pharma.data[pharma.data$dose=='A', ]
pharma.data.B<-pharma.data[pharma.data$dose=='B', ]
pharma.data.C<-pharma.data[pharma.data$dose=='C', ]
pharma.data.D<-pharma.data[pharma.data$dose=='D', ]

#getting indexis
a.index=which(pharma.data$dose=='A')
a.index= c(a.index[1], tail(a.index,n=1))
b.index=which(pharma.data$dose=='B')
b.index= c(b.index[1], tail(b.index,n=1))
c.index=which(pharma.data$dose=='C')
c.index= c(c.index[1], tail(c.index,n=1))
d.index=which(pharma.data$dose=='D')
d.index= c(d.index[1], tail(d.index,n=1))

get.resids<-function(data.slice,y.lab){
  out<-rq(log(data.slice$TBL.M)~log(data.slice$TBL.B),tau = 0.5)
  res<-residuals(out)
  plot(res,ylab=y.lab, xlab='Patient', ylim = c(-0.67,0.67), pch=19)
  return(res)
}

#residual data to work on
par(mfrow=c(2,2))
pharma.data.A=get.resids(pharma.data.A, y.lab='Bilirubin A')
pharma.data.B=get.resids(pharma.data.B, y.lab='Bilirubin B')
pharma.data.C=get.resids(pharma.data.C, y.lab='Bilirubin C')
pharma.data.D=get.resids(pharma.data.D, y.lab='Bilirubin D')
full.pharma.data<-c(pharma.data.A,pharma.data.B,pharma.data.C,pharma.data.D)

#modified fitting function where tkappa and/or xi remail the same 
modified.egp.fit<-function (xdat,indices, thresh, model = c("egp1", "egp2", "egp3"), init, 
                            show = FALSE, commonality=0) {
  # if commonality == 0 -> same shape (kappa parameter)
  # if commonality == 1 -> same shape (kappa parameter) and tail index (xi parameter)
  if (!(model %in% c("egp1", "egp2", "egp3")) || length(model) != 
      1) {
    stop("Invalid model  argument: must be one of `egp1', `egp2' or `egp3'.")
  }
  if (length(thresh) > 1) {
    warning("Length of threshold vector greater than one. Selecting first component.")
    thresh <- thresh[1]
  }
  changinit <- missing(init)
  if (!changinit) {
    if (any(init[1:2] < 0)) {
      changinit <- TRUE
    }
  }
  if (changinit) {
    init <- c(kappa = 1.01, suppressWarnings(fit.gpd(xdat, 
                                                     threshold = thresh[1], show = FALSE)$est))
  }
  xdata = xdat[xdat > thresh]
  xmax <- max(xdata)
  #init<-c(1.5,1,-0.2)
  mle <- alabama::auglag(par = c(init,init,init,init), fn = function(par, xdat, 
                                                                     thresh, model) {
    -egp.ll(par = par[1:3], xdat = na.omit(xdat[indices[[1]][1]:indices[[1]][2]]), thresh = thresh, model = model)-
      egp.ll(par = par[4:6], xdat = na.omit(xdat[indices[[2]][1]:indices[[2]][2]]), thresh = thresh, model = model)-
      egp.ll(par = par[7:9], xdat = na.omit(xdat[indices[[3]][1]:indices[[3]][2]]), thresh = thresh, model = model)-
      egp.ll(par = par[10:12], xdat = na.omit(xdat[indices[[4]][1]:indices[[4]][2]]), thresh = thresh, model = model)
  }, hin = function(par, ...) {
    
    xi.restriction<-c(par[3]-par[6]- 1e-10,par[6]-par[3]- 1e-10,
                      par[6]-par[9]- 1e-10,par[9]-par[6] - 1e-10,
                      par[12]-par[9]- 1e-10,par[9]-par[12] - 1e-10)
    kappa.restriction<-c(par[1]-par[4]- 1e-10,par[4]-par[1]- 1e-10,
                         par[4]-par[7]- 1e-10,par[7]-par[4] - 1e-10,
                         par[10]-par[7]- 1e-10,par[7]-par[10] - 1e-10)
    restriction<-c()
    if(commonality==0){
      restriction<-kappa.restriction
    }
    else if(commonality ==1){
      restriction<-c(kappa.restriction,xi.restriction)
    }
    c(restriction,par[1] - 1e-10, par[2] - 1e-10, par[3] + 1, ifelse(par[3] < 
                                                                       0, thresh - par[2]/par[3] - xmax, 1))
  }, xdat = xdata, thresh = thresh, model = model, control.outer = list(trace = FALSE, 
                                                                        method = "BFGS"), control.optim = list(maxit = 500, reltol = 1e-10))
  
  fitted <- list()
  fitted$estimate <- fitted$param <- mle$par
  fitted$deviance <- 2 * mle$value
  fitted$nllh <- mle$value
  # if (mle$convergence == 0) {
  #   fitted$convergence <- "successful"
  #   fitted$vcov <- try(solve(mle$hessian))
  #   fitted$std.err <- try(sqrt(diag(fitted$vcov)))
  #   if (is.character(mle$se) || mle$par[3] < -0.5) {
  #     fitted$vcov <- NULL
  #     fitted$se <- rep(NA, 3)
  #   }
  # }
  # else {
  #   fitted$convergence <- mle$convergence
  #   warning("Maximization routine may have failed; check output and try providing better starting values")
  # }
  print(fitted)
  # names(fitted$estimate) <- names(fitted$std.err) <- c("kappa", 
  #                                                      "scale", "shape","kappa", 
  #                                                      "scale", "shape","kappa", 
  #                                                      "scale", "shape","kappa", 
  #                                                      "scale", "shape")
  # if (!is.null(fitted$vcov)) {
  #   colnames(fitted$vcov) <- rownames(fitted$vcov) <- c("kappa", 
  #                                                       "scale", "shape","kappa", 
  #                                                       "scale", "shape","kappa", 
  #                                                       "scale", "shape","kappa", 
  #                                                       "scale", "shape")
  # }
  fitted$counts <- mle$counts
  fitted$threshold <- thresh
  fitted$nat <- length(xdata)
  fitted$pat <- length(xdata)/length(xdat)
  fitted$exceedances <- xdata
  #fitted$hessian <- mle$hessian
  fitted$method <- "copt"
  fitted$model <- model
  class(fitted) <- c("mev_egp")
  if (show) {
    print(fitted)
  }
  return(invisible(fitted))
}



us<-seq(-0.65, .15,(0.15+0.65)/19 )

data.A.kappa<-data.B.kappa<-data.C.kappa<-data.D.kappa<-full.data.kappa<-c()
data.A.xi<-data.B.xi<-data.C.xi<-data.D.xi<-full.data.xi<-c()

kap.lb.ci<-kap.ub.ci<-c()
xi.lb.ci<-xi.ub.ci<-c()

for (u in 1:length(us)){
  
  #fit egp1 on each of the datasets
  a.out<-fit.egp(pharma.data.A, us[u] , model='egp1')
  b.out<-fit.egp(pharma.data.B, us[u] , model='egp1')
  c.out<-fit.egp(pharma.data.C, us[u] , model='egp1')
  d.out<-fit.egp(pharma.data.D, us[u] , model='egp1')
  full.out.same.kappa<-modified.egp.fit(full.pharma.data, list(a.index, b.index,c.index,d.index),
                             us[u] , model='egp1',commonality = 0)
  full.out.same.kappa.xi<-modified.egp.fit(full.pharma.data, list(a.index, b.index,c.index,d.index),
                                        us[u] , model='egp1',commonality = 1)
  
  #set confidence intervals
  a.ci=confint(a.out,level=0.95)
  b.ci=confint(b.out,level=0.95)
  c.ci=confint(c.out,level=0.95)
  d.ci=confint(d.out,level=0.95)
  kap.lb.ci[u]= min(a.ci[1,1],b.ci[1,1],c.ci[1,1],d.ci[1,1] ,na.rm = T)
  kap.ub.ci[u]= max(a.ci[1,2],b.ci[1,2],c.ci[1,2],d.ci[1,2] ,na.rm = T)
  xi.lb.ci[u]= min(a.ci[3,1],b.ci[3,1],c.ci[3,1],d.ci[3,1] ,na.rm = T)
  xi.ub.ci[u]= max(a.ci[3,2],b.ci[3,2],c.ci[3,2],d.ci[3,2] ,na.rm = T)

  data.A.kappa[u]=a.out$estimate[1]
  data.B.kappa[u]=b.out$estimate[1]
  data.C.kappa[u]=c.out$estimate[1]
  data.D.kappa[u]=d.out$estimate[1]
  full.data.kappa[u]=full.out.same.kappa$estimate[1]
  
  data.A.xi[u]=a.out$estimate[3]
  data.B.xi[u]=b.out$estimate[3]
  data.C.xi[u]=c.out$estimate[3]
  data.D.xi[u]=d.out$estimate[3]
  full.data.xi[u]=full.out.same.kappa.xi$estimate[3]
  
  
}

par(mfrow=c(1,2))
# plot kappa
plot(us, data.A.kappa, pch='A',type="b",cex=0.5, ylim=c(0,6), xlab='u', ylab=expression(hat(kappa)))
lines(us,data.B.kappa,  pch='B',type="b",cex=0.5)
lines(us,data.C.kappa,  pch='C',type="b",cex=0.5)
lines(us,data.D.kappa,  pch='D',type="b",cex=0.5)
lines(us,full.data.kappa, type='l', col='grey',lwd=3)
polygon(c(us,rev(us)),c(kap.lb.ci,rev(kap.ub.ci)),col=rgb(0,0,0,0.1), border = FALSE)
lines(us, kap.lb.ci, lty = 'dashed', col = 'red')
lines(us, kap.ub.ci, lty = 'dashed', col = 'red')

# plot xi
plot(us, data.A.xi, pch='A',type="b",cex=0.5, ylim=c(-0.85,0.4), xlab='u', ylab=expression(hat(xi)))
lines(us,data.B.xi,  pch='B',type="b",cex=0.5)
lines(us,data.C.xi,  pch='C',type="b",cex=0.5)
lines(us,data.D.xi,  pch='D',type="b",cex=0.5)
lines(us,full.data.xi, type='l', col='grey',lwd=3)
polygon(c(us,rev(us)),c(xi.lb.ci,rev(xi.ub.ci)), border = NA,col=rgb(0,0,0,0.1))
lines(us, xi.lb.ci, lty = 'dashed', col = 'red')
lines(us, xi.ub.ci, lty = 'dashed', col = 'red')


#qq plots for gp
obs.data<-full.pharma.data[full.pharma.data>0.1047831]
n<-length(obs.data)
ps<-seq(1/(n+1),n/(n+1), 1/(n+1) )
get.quants<-function(a, u,p){
  u + (a[1] * (p^(-a[2]) - 1))/a[2]
}

plot(get.quants(c(0.21,-0.27),0.1047831, 1-ps ), sort(obs.data),ylim=c(-0.15,0.75), xlim=c(-0.15,0.75),xlab='Fitted GP', ylab='Observed')
abline( a=0, b=1)

#qq plots for egp 
obs.data<-full.pharma.data[full.pharma.data>-0.1326027]
n<-length(obs.data)
ps<-seq(1/(n+1),n/(n+1), 1/(n+1) )
egpd.out<-fit.egp(full.pharma.data,-0.1326027 , model='egp1',show=T)
egpd.ci<-confint(egpd.out)
egp.estim<-as.vector(egpd.out$estimate)
fitted.data<-c()
for(i in 1:n){
  fitted.data[i]<-egp.retlev(obs.data,-0.1326027,egp.estim,model='egp1',p=1-ps[i],plot=F)
}
plot(fitted.data, sort(obs.data), ylim=c(-0.15,0.75), xlim=c(-0.15,0.75), xlab='Fitted EGP1', ylab='Observed')
abline( a=0, b=1)

########################################################  Precipitation data - Rain########################################################  
library(readr)

daily_rain_2015 <- read_csv("Desktop/McGill U/Fall 2020/Math 80622A/project/daily rain 2015.csv")
daily_rain_2016 <- read_csv("Desktop/McGill U/Fall 2020/Math 80622A/project/daily rain 2016.csv")
daily_rain_2017 <- read_csv("Desktop/McGill U/Fall 2020/Math 80622A/project/daily rain 2017.csv")
daily_rain_2018 <- read_csv("Desktop/McGill U/Fall 2020/Math 80622A/project/daily rain 2018.csv")
daily_rain_2019 <- read_csv("Desktop/McGill U/Fall 2020/Math 80622A/project/daily rain 2019.csv")
daily_rain_2020 <- read_csv("Desktop/McGill U/Fall 2020/Math 80622A/project/daily rain 2020.csv")

get.seasonal.data<-function(daily.rain, season){
  daily.rain$date <- as.Date(daily.rain$`Date/Time`, format= "%Y-%m-%d")
  return(subset(daily.rain, date>= season[1] & date <= season[2]))
}

summer.2015=get.seasonal.data(daily_rain_2015, c('2015-05-01','2015-08-31'))
summer.2016=get.seasonal.data(daily_rain_2016, c('2016-05-01','2016-08-31'))
summer.2017=get.seasonal.data(daily_rain_2017, c('2017-05-01','2017-08-31'))
summer.2018=get.seasonal.data(daily_rain_2018, c('2018-05-01','2018-08-31'))
summer.2019=get.seasonal.data(daily_rain_2019, c('2019-05-01','2019-08-31'))
summer.2020=get.seasonal.data(daily_rain_2020, c('2020-05-01','2020-08-31'))


all.summer<-rbind(summer.2015,summer.2016,summer.2017,summer.2018,summer.2019,summer.2020)

rain.data<-all.summer$`Total Rain (mm)`
rain.data<-rain.data[rain.data>0]
rain.data<-rain.data[!is.na(rain.data)]
plot(rain.data, pch=19, ylab='Total rain (mm)', main='Daily Rain (Summer 2015 - Summer 2020)')

us<-quantile(rain.data, c(0,.80))
us<- seq(us[1], us[2], (us[2]-us[1])/29)

#### analyzing parameter stability over grid of thresholds #### 

#initialize list of items to plot
gp.shape<- gp.scale<- c()
gp.shape.ci<-gp.scale.ci<- matrix(NA, nrow=length(us), ncol = 2)

egp.scale<-egp.shape<-egp.kap<-c()
egp.shape.ci<-egp.scale.ci<-egp.kap.ci<- matrix(NA, nrow=length(us), ncol = 2)

#for each threshold, fit gp and EGP1, the find parameter estimates and their CI (asymptotic normality)
for (i in 1:length(us)){
  
  # Standard GP params 
  gp.out<-fit.gpd(rain.data, us[i])
  gp.scale<-c(gp.scale,gp.out$estimate[[1]] )
  gp.shape<-c(gp.shape,gp.out$estimate[[2]] )
  
  gp.scale.ci[i,]<- c(gp.out$estimate[[1]]-1.96*gp.out$std.err[[1]],gp.out$estimate[[1]]+1.96*gp.out$std.err[[1]])
  gp.shape.ci[i,]<- c(gp.out$estimate[[2]]-1.96*gp.out$std.err[[2]],gp.out$estimate[[2]]+1.96*gp.out$std.err[[2]])
  
  # EGP1 params
  temp.data<-rain.data[rain.data>us[i]]
  egp.out<-fit.egp(temp.data, us[i], model='egp1')
  egp.kap<-c(egp.kap,egp.out$estimate[[1]] )
  egp.scale<-c(egp.scale,egp.out$estimate[[2]] )
  egp.shape<-c(egp.shape,egp.out$estimate[[3]] )
  egp.ci<-confint(egp.out,level=0.95)
  egp.kap.ci[i,]<- c(egp.out$estimate[[1]]-1.96*egp.out$std.err[[1]],egp.out$estimate[[1]]+1.96*egp.out$std.err[[1]])
  egp.scale.ci[i,]<- c(egp.out$estimate[[2]]-1.96*egp.out$std.err[[2]],egp.out$estimate[[2]]+1.96*egp.out$std.err[[2]])
  egp.shape.ci[i,]<- c(egp.out$estimate[[3]]-1.96*egp.out$std.err[[3]],egp.out$estimate[[3]]+1.96*egp.out$std.err[[3]])
}

#plot egp1 shape
e.shape=ggplot()+ geom_point(mapping = aes(x = us,y=egp.shape))+geom_point()+
  geom_ribbon(aes(x=us,ymin= egp.shape.ci[,1], ymax=egp.shape.ci[,2]),fill='grey50',alpha=0.5,linetype=2, colour='red')+
  xlab('u')+ylab(expression(hat(xi)['EGP1']))

#plot egp1 scale
e.scale=ggplot()+ geom_point(mapping = aes(x = us,y=egp.scale))+geom_point()+
  geom_ribbon(aes(x=us,ymin= egp.scale.ci[,1], ymax=egp.scale.ci[,2]),na.rm = T,fill='grey50',alpha=0.5,linetype=2, colour='red')+ylim(0,20)+
  xlab('u')+ylab(expression(hat(sigma)['EGP1']))#+xlim(c(0,16))

e.kap=plot.graph.exp2(egp.kap,egp.kap.ci)+geom_hline(yintercept=1.)+geom_point()+
  xlab('u')+ylab(expression(hat(kappa)['EGP1']))


g.shape=plot.graph.exp2(gp.shape,gp.shape.ci)+ylim(-0.5,1)+geom_point()+
  xlab('u')+ylab(expression(hat(xi)['GP']))
g.scale=plot.graph.exp2(gp.scale,gp.scale.ci)+ylim(0,20)+geom_point()+
  xlab('u')+ylab(expression(hat(sigma)['GP']))

grid.arrange(e.shape,g.shape,e.scale,g.scale, e.kap, nrow=3, ncol=2)

# histogram with threshold excedences
chosen.threshold=us[1]
egp1.out<-fit.egp(rain.data,chosen.threshold , model='egp1')
points.y<-c()
points.x<-seq(chosen.threshold,max(rain.data), 0.25)
for(i in 1:length(points.x)){
  points.y<-c(points.y,egp.ll(points.x[i], chosen.threshold,as.vector(egp1.out$estimate), model='egp1' ))
}
  
density.plt<-ggplot()+
  geom_line(data=data.frame('x'=points.x,'y'=exp(points.y)) ,aes(x=x,y=y))+
  geom_histogram(data=as.data.frame(rain.data),aes(x=rain.data,y=..density..),color="black",alpha=0.2,bins=60)

grid.arrange(e.shape,g.shape,e.scale,g.scale, e.kap,density.plt, nrow=3, ncol=2)


############# Return level plots #############

s<-quantile(rain.data, c(0,.80))
us<- seq(us[1], us[2], (us[2]-us[1])/8)


return.levels<- c(10,50,100,200,500)
gp.out.mat<-matrix(ncol=length(return.levels), nrow=length(us))
egp1.out.mat<-matrix(ncol=length(return.levels), nrow=length(us))
egp2.out.mat<-matrix(ncol=length(return.levels), nrow=length(us))

for (u in 1:length(us)){
  temp.data<-rain.data[rain.data>us[u]]
  
  #fit egps 
  egp1.out<-fit.egp(rain.data, us[u], model='egp1',show=F)
  egp2.out<-fit.egp(rain.data, us[u], model='egp2',show=F)
  
  #calculate each of the return levels
  for(retlev in 1:length(return.levels)){
    egp1.out.mat[u,retlev]<-egp.retlev(rain.data,us[u],as.vector(egp1.out$estimate),model='egp1',p=1/return.levels[retlev],plot=F)
    egp2.out.mat[u,retlev]<-egp.retlev(rain.data,us[u],as.vector(egp2.out$estimate),model='egp2',p=1/return.levels[retlev],plot=FALSE)
    gp.out.mat[u,retlev]<-gpd.mle(rain.data-us[u],args='quant',p=1/return.levels[retlev])
  }
}


plot.retlev<-function(retlev.mat, ret.levs=return.levels,titl=''){
  init.vec=retlev.mat[1,]
  plot(ret.levs, init.vec, ylim=c(0,100), pch='1',type="b",xlim=c(10,500),log='x',
       ylab='Estimated Return Level',
       xlab='Return Period (years)', main=titl)
  
  for (i in 2:(dim(retlev.mat)[1])){
    lines(ret.levs,retlev.mat[i,],pch=as.character(i),type="b")
  }
}

par(mfrow=c(1,3))
plot.retlev(egp1.out.mat,titl= 'EGP1')
plot.retlev(gp.out.mat,titl= 'GP')
plot.retlev(egp2.out.mat,titl= 'EGP2')
########################################################  Precipitation data - Snow########################################################  


























same.shape<-modified.egp.fit(full.pharma.data, list(a.index, b.index,c.index,d.index),
                             us[6] , model='egp1',commonality = 0)


