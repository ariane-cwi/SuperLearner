#############
# LIBRARIES #
#############

# library(devtools)
# install_github(repo="ariane-cwi/conformal", subdir="conformalInference")

setwd("~/GitHub/SLsurv")
libraries_sources = function(){
  library(survival)
  library(FastPseudo)
  library(pseudo)
  library(randomForest)
  library(randomForestSRC)
  library(tidyverse)
  library(gridExtra)
  library(haven)
  library(devtools)
  library(knitr)
  library(Hmisc)
  library(parallel)
  library(conformalInference)
  library(SuperLearner)
  source("~/GitHub/SLsurv/pseudomean2.R")
  source("~/GitHub/SuperLearner/R/SuperLearnerPseudo.R")
  source("~/GitHub/SuperLearner/R/internals.R")
}
libraries_sources()


########################
# SIMULATION FUNCTIONS #
########################

simulate.A = function(n,seed=NULL,p=2,a=3,beta_tilde0 = c(5.5,2.5,2.5), 
                      beta0 = c(5.5,2.097,2.097,3.16), lambda=0.07,
                      C.indep.Z=T, k=12,nu=6,b=c(2,1),tau=8.8){
  if(!is.null(seed)){set.seed(seed)}
  x = matrix(rbinom(n*p,1,0.5), nrow=n)
  t = cbind(rep(1,n), x) %*% beta_tilde0 + runif(n,-a,a)
  if(C.indep.Z){c = rexp(n,rate=lambda)}
  else{c = k*(-log(1-runif(n))*exp(-as.double(x%*%b)))**(1/nu)}
  tobs = pmin(t,c)
  delta = (t <= c)*1
  mu_tau = cbind(rep(1,n),x[,1]*(1-x[,2]),x[,2]*(1-x[,1]),x[,1]*x[,2])%*% beta0
  return(list(t=t,c=c,tobs=tobs,delta=delta,x=data.frame(x),mu_tau=mu_tau,tau=tau))
}


simulate.B = function(n,seed=NULL,a=5,k=2,nu=6,p=3,b=c(2,1,0),lambda=0.3,tau=3.6){
  if(!is.null(seed)){set.seed(seed)}
  x = matrix(runif(n*p,-a,a),n,p)
  u = runif(n)
  t = k*(-log(1-u)*exp(-as.double(x%*%b)))**(1/nu)
  c = rexp(n,lambda)
  tobs = pmin(t,c)
  delta = (t <= c)*1
  return(list(t=t,c=c,tobs=tobs,delta=delta,x=data.frame(x),tau=tau))
}


simulate.C = function(n,seed=NULL,a=5,k=2,nu=6,p.bin=6,p.cont=9,lambda=0.3,tau=2.8){
  if(!is.null(seed)){set.seed(seed)}
  xbin = matrix(rbinom(n*p.bin,1,prob=0.4),n,p.bin)
  xcont = matrix(runif(n*p.cont,0,1),n,p.cont)
  x = cbind(xcont[,1],xbin[,1],xcont[,2],xbin[,2],xcont[,3],xbin[,3],xcont[,c(4,5)],xbin[,4],xcont[,6],xbin[,-seq(1,4)],xcont[,-seq(1,6)])
  y = x[,3] - 3*x[,5] + 2*x[,1]*x[,10] + 4*x[,2]*x[,7] + 3*x[,4]*x[,5] - 5*x[,6]*x[,10] + 3*x[,8]*x[,9] + x[,1]*x[,4] - 2*x[,6]*x[,9] - 4*x[,3]*x[,4] - x[,7]*x[,8]
  u = runif(n)
  t = k*(-log(1-u)*exp(-as.double(y)))**(1/nu)
  c = rexp(n,lambda)
  tobs = pmin(t,c)
  delta = (t <= c)*1
  return(list(t=t,c=c,tobs=tobs,delta=delta,x=data.frame(x),tau=tau))
}


###################
# PARALLELISATION #
###################

parallelisation_init = function(){
  cl <<- makeCluster(getOption("cl.cores", detectCores()))
  invisible(clusterEvalQ(cl, {
    library(survival)
    library(FastPseudo)
    library(pseudo)
    library(randomForest)
    library(randomForestSRC)
    library(tidyverse)
    library(gridExtra)
    library(haven)
    library(devtools)
    library(knitr)
    library(Hmisc)
    library(parallel)
    library(conformalInference)
    library(SuperLearner)
    source("~/GitHub/SLsurv/pseudomean2.R")
    source("~/GitHub/SuperLearner/R/SuperLearnerPseudo.R")
    source("~/GitHub/SuperLearner/R/internals.R")
  }))
  clusterExport(cl=cl,c("simulate.A","simulate.B","simulate.C"))
}


############################################
# 1/ STANDARD VS SPLIT PSEUDO-OBSERVATIONS #
############################################

# On sélectionne un pourcentage rho de données sur lesquelles calculer les
# pseudo-observations split (D2), le reste est utilisé pour les calculer (D1).

simulate = simulate.A
N = c(20,50,100)
rho = c(0.2,0.5,0.8)
nb.repeat = 50
set.seed(0)

standard_VS_split = purrr::cross_df(list(n=N,p=rho,m=1:nb.repeat)) %>% 
  mutate(df = map(n, ~ simulate(.)),
         idx.pobs = map2(n,p, ~ sample(1:.x,floor(.x*.y))),
         pobs = map2(df,idx.pobs, ~ pseudomean(.x$tobs,.x$delta,.x$tau)[.y]),
         pobs2 = map2(df,idx.pobs, ~ pseudomean2(.x$tobs,.x$delta,.x$tau,.y)$pseudo)) %>%
  select(-df,-idx.pobs) %>%
  unnest(c(pobs, pobs2)) 

pdf("fig/stand_vs_split.pdf",h=3.5,w=5)
standard_VS_split %>% 
  # Comment for not reducing the number of points in the scatter plots
  group_by(n, p) %>% 
  slice_sample(n = min(min(N)*min(rho)*nb.repeat,300), replace = TRUE) %>%
  ##
  ggplot(aes(x = pobs, y = pobs2)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 0.4) + 
  labs(x = 'Pseudo-Observations (Standard)', y = 'Pseudo-Observations (Split)') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(n ~ p, labeller = label_both)
dev.off()


###############################################################
# 2/ STANDARD VS SPLIT PSEUDO-OBSERVATIONS FOR SUPER LEARNING #
###############################################################


SLcomparison = function(n,v,n0=N0,sl.library=SL.library){
  # Training set and pseudo-observations
  df = simulate(n)
  pobs = pseudomean(df$tobs,df$delta,df$tau)
  # Test set
  df0 = simulate(n0)
  df0$t.min.tau = pmin(df0$t,df0$tau)
  MSE = function(pred){mean((df0$t.min.tau-pred)**2)} 
  # Super-Learner training, on true restricted times / pobs / split pobs
  cvControl = list(V=v,stratifyCV = F,shuffle = T,validRows = NULL)
  out.SL = SuperLearner(pmin(df$t,df$tau),data.frame(df$x),newX = data.frame(df0$x),
                        cvControl = cvControl, SL.library = sl.library)
  out.SLpobs = SuperLearner(pobs,data.frame(df$x),newX = data.frame(df0$x),
                            cvControl = cvControl, SL.library = sl.library)
  out.SLsplit = SuperLearnerPseudo(df$tobs,df$delta,data.frame(df$x),tmax=df$tau,newX = data.frame(df0$x),
                                   cvControl = cvControl, SL.library = sl.library)
  
  m = length(sl.library)
  res <- data.frame(mse = c(sapply(1:m,function(s) MSE(out.SL$library.predict[,s])),
                            MSE(out.SL$SL.predict)),
                    beta = c(as.double(out.SL$coef),NA),
                    mseP = c(sapply(1:m,function(s) MSE(out.SLpobs$library.predict[,s])),
                             MSE(out.SLpobs$SL.predict)),
                    betaP = c(as.double(out.SLpobs$coef),NA),
                    mseS = c(sapply(1:m,function(s) MSE(out.SLsplit$library.predict[,s])),
                             MSE(out.SLsplit$SL.predict)),
                    betaS = c(as.double(out.SLsplit$coef),NA)
  )
  return(res)
}

SLcomparison.pll = function(x){
  set.seed(x)
  return(purrr::cross_df(list(n=N,v=V)) %>% 
           mutate(res = map2(n,v, ~ SLcomparison(.x,.y))) )
}

##################
simulate = simulate.C
SLtypes = c("True restricted event times","Pseudo-Observations","Split Pseudo-Observations")
SL.library = c("SL.lm", "SL.glmnet", "SL.gam", "SL.ranger", "SL.nnet")
models.names = c("LM","Lasso","GAM","RF","NN","SL")
color_mapping = c(LM = "#F8766D", Lasso = "#B79F00", GAM = "#00BA38", 
                  RF = "#00B6EB", NN = "#A58AFF", SL = "#F564E3")
N = c(100,200,300,400,500)
N0 = 1000 
V = 6
nb.repeat = 80
##################

parallelisation_init()
clusterExport(cl,c("simulate","SL.library","N","N0","V","SLcomparison"))
resSLcomp = parLapply(cl, X = 1:nb.repeat, fun = SLcomparison.pll)
stopCluster(cl)


mseSLcomp = purrr::cross_df(list(model = models.names, n=N, v=V, rep = 1:nb.repeat, SLtype=SLtypes)) %>%
  mutate(mse = as.double(unlist(data.frame(data.table::rbindlist(data.table::rbindlist(resSLcomp)$res))[,c(1,3,5)])))
mseSLcomp$SLtype = factor(mseSLcomp$SLtype, levels = SLtypes)

# pdf(file="fig/mse_SL.pdf",w=6,h=3)
# mseSLcomp %>% filter(model=="SL") %>%
#   ggplot(aes(x = as.factor(n), y = mse)) +
#   geom_boxplot(lwd=0.3,outlier.size = 0.3) +
#   labs(x = 'Training size', y = 'MSE') +
#   coord_cartesian(ylim = c(0.05, 0.25)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   facet_grid(~ SLtype)
# graphics.off()

mseSLcomp.tab = mseSLcomp %>% group_by(v,n,model,SLtype) %>% 
  dplyr::summarise(median = round(median(mse,na.rm = T),3),
                   Q1 = round(quantile(mse,0.25,na.rm = T),3),
                   Q3 = round(quantile(mse,0.75,na.rm = T),3))


pdf(file="fig/mse_SL_graph.pdf",w=7,h=4)
mseSLcomp.tab %>% ggplot(aes(x = n, y = median, color = model, group = model)) +
  geom_line(size = 0.4) +
  geom_ribbon(aes(ymin = Q1, ymax = Q3, fill = model), alpha = 0.1, size = 0.4, linetype = "dashed") +
  geom_point(shape = 4, size = 2) +
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  labs(x = "n", y = "MSE") +
  facet_grid(~SLtype) +
  # coord_cartesian(ylim = c(0.175, 0.33)) + # scheme C
  theme_minimal()
graphics.off()


betaSLcomp = purrr::cross_df(list(model = models.names, n=N, v=V, rep = 1:nb.repeat, SLtype=SLtypes)) %>%
  mutate(beta = as.double(unlist(data.frame(data.table::rbindlist(data.table::rbindlist(resSLcomp)$res))[,c(2,4,6)])))
betaSLcomp$SLtype = factor(betaSLcomp$SLtype, levels = SLtypes)

betaSLcomp.tab = betaSLcomp %>% filter(model != "SL") %>% 
  group_by(v,n,model,SLtype) %>% 
  dplyr::summarise(median = round(median(beta,na.rm = T),3),
                   Q1 = round(quantile(beta,0.25,na.rm = T),3),
                   Q3 = round(quantile(beta,0.75,na.rm = T),3))


pdf(file="fig/beta_SL_graph.pdf",w=7,h=4)
betaSLcomp.tab %>% ggplot(aes(x = n, y = median, color = model, group = model)) +
  geom_line(size = 0.4) +
  geom_ribbon(aes(ymin = Q1, ymax = Q3, fill = model), alpha = 0.1, size = 0.4, linetype = "dashed") +
  geom_point(shape = 4, size = 2) +
  scale_color_manual(values = color_mapping[-length(models.names)]) +
  scale_fill_manual(values = color_mapping[-length(models.names)]) +
  labs(x = "n", y = "Weight") +
  facet_grid(~SLtype) +
  theme_minimal()
graphics.off()


#############################################
# 3/ COMPARISON WITH OTHER SURVIVAL METHODS #
#############################################

### TRAIN AND PREDICT FUNCTIONS ###

integral = function(f,min=0,max){
  return(integrate(f,min,max,subdivisions=3000,rel.tol = 1e-4)$value)
}

train.fun = function(x,t,d,tau,pred.models=c("all"),v=10,
                     sl.library = c("SL.glmnet","SL.ranger","SL.nnet")) {
  bincols = colMeans((x==1|x==0), na.rm=T) == 1
  for(i in 1:length(bincols)){
    # Turn the binary columns into factors
    if(bincols[i]){x[[i]] = as.factor(x[[i]])}
  }
  datafr = data.frame(x)
  
  pred.list = list()

  # Cox
  if("cox" %in% pred.models || "all" %in% pred.models){
    pred.list$out.cox = coxph(Surv(time=t,event=d)~., data = datafr)
  }
  # Random Survival Forests
  if("rsf" %in% pred.models || "all" %in% pred.models){
    pred.list$out.rsf = rfsrc(Surv(time=t,event=d)~.,data=data.frame(t=t,d=d,x),
                              importance=T, splitrule="bs.gradient", mtry=ncol(datafr))
  }
  # Pseudo-observations + Super Learner
  if("sl" %in% pred.models || "all" %in% pred.models){
    pobs = pseudoKM(t,as.logical(d),tau)$pseudoval
    cvControl = list(V=v,stratifyCV = F,shuffle = T,validRows = NULL)
    pred.list$out.SLpobs = SuperLearner::SuperLearner(pobs,datafr,SL.library = sl.library,
                                                      cvControl = cvControl)
  }
  return(pred.list)
}

predict.fun = function(out,x0,tau,pred.models=c("all")) {
  predictions = numeric(0)
  
  # Cox
  if("cox" %in% pred.models || "all" %in% pred.models){
    surv.cox = survfit(out$out.cox,newdata=data.frame(x0))
    pred.cox = numeric(nrow(x0))
    for(i in 1:nrow(x0)){
      cox.curve = stepfun(surv.cox$time,c(1,surv.cox$surv[,i]))
      pred.cox[i] = integral(cox.curve,0,tau)
    }
    predictions = c(predictions,pred.cox)
  }
  # Random Survival Forests
  if("rsf" %in% pred.models || "all" %in% pred.models){
    obj.pred.rsf = predict(out$out.rsf,data.frame(x0))
    vals_list = split(obj.pred.rsf$survival, row(obj.pred.rsf$survival))
    surv.rsf = purrr::map(vals_list, ~ stepfun(obj.pred.rsf$time.interest,c(1,.),f=1))
    pred.rsf = purrr::map_dbl(surv.rsf, ~integral(.,0,tau))
    predictions = c(predictions,as.double(pred.rsf))
  }
  # Pseudo-observations + Super Learner
  if("sl" %in% pred.models || "all" %in% pred.models){
    pred.SLpobs = predict(out$out.SLpobs,data.frame(x0))$pred
    predictions = c(predictions,pred.SLpobs)
  }
  
  return(predictions)
}

### PARALLELISATION ###

SurvComparison = function(n,n0=N0){
  # Data sets
  df0 = simulate(n0)
  df0$t.min.tau = pmin(df0$t,df0$tau)
  df = simulate(n)
  
  # Training 
  out = train.fun(df$x,df$tobs,df$delta,df$tau,pred.models=models,v=V,SL.library)
  # Test
  x0 = df0$x
  bincols = colMeans((x0==1|x0==0), na.rm=T) == 1
  for(i in 1:length(bincols)){
    # Turn the binary columns into factors
    if(bincols[i]){x0[[i]] = as.factor(x0[[i]])}
  }
  preds = matrix(predict.fun(out,x0,df$tau,pred.models=models),nrow = n0)
  mse = apply((df0$t.min.tau - preds)**2,2,mean)
  
  return(mse)
}

SurvComparison.pll = function(x){
  set.seed(x)
  return(purrr::cross_df(list(n=N)) %>% 
           mutate(mse = map(n, ~ SurvComparison(.))) )
}

### APPLICATION ###

##################
simulate = simulate.C
models = c("cox","rsf","sl")
models.names = c("Cox","RSF","P.obs.+SL")
SL.library = c("SL.lm", "SL.glmnet", "SL.gam", "SL.ranger", "SL.nnet")
N = c(200,500,1000,1500,2000)
N0 = 1000 
V = 6
nb.repeat = 80
##################

#save(mseSurvComp,file="fig/mseSurvComp.Rdata")

parallelisation_init()
clusterExport(cl,c("simulate","N","N0","V","models","SL.library","train.fun",
                   "predict.fun","integral","SurvComparison"))
mseSurvComp = parLapply(cl, X = 1:nb.repeat, fun = SurvComparison.pll)
stopCluster(cl)

mseSurvComp = purrr::cross_df(list(model = models.names, n=N, rep = 1:nb.repeat)) %>%
  mutate(mse = unlist(data.table::rbindlist(mseSurvComp)$mse))
mseSurvComp$model = factor(mseSurvComp$model, levels = models.names)

mseSurvComp.tab = mseSurvComp %>% group_by(n,model) %>% 
  dplyr::summarise(mean = round(mean(mse),3),
                   median = round(median(mse),3),
                   Q1 = round(quantile(mse,0.25),3),
                   Q3 = round(quantile(mse,0.75),3)) 

pdf("fig/mse_surv_tab.pdf",h=nrow(mseSurvComp.tab)/3.5,w=6)
grid.table(mseSurvComp.tab)
dev.off()

pdf(file="fig/mse_surv_graph.pdf",w=3.5,h=4)
ggplot(mseSurvComp.tab, aes(x = n, y = median, color = model, group = model)) +
  geom_line() +
  geom_ribbon(aes(ymin = Q1, ymax = Q3, fill = model), alpha = 0.1, linetype = "dashed") +
  geom_point(shape = 4, size = 2) +
  labs(x = "n", y = "MSE") +
  theme_minimal()
graphics.off()


#############################
# 4/ REAL DATA APPLICATIONS #
#############################

############
# BRCANCER #
############


# The outcome of interest is the time to relapse, not the time to death
brcancer0 = read_dta(url("https://www.stata-press.com/data/r16/brcancer.dta"))
brcancer = as.data.frame(as.matrix(brcancer0))
brcancer = brcancer %>% select(-id, -x4, -x4b, -x5e)
brcancer$x2 = brcancer$x2 - 1

# Events
tobs.brcancer = brcancer$rectime
delta.brcancer = brcancer$censrec
tau.brcancer = quantile(tobs.brcancer,0.9,names=F)
x.brcancer = brcancer %>% select(-rectime,-censrec)
p = dim(x.brcancer)[2]


#####
models = c("cox","rsf","sl")
models.names = c("Cox","RSF","P.obs.+SL")
SL.library = c("SL.lm", "SL.glmnet", "SL.gam", "SL.ranger", "SL.nnet")
V = 6
m = length(models.names)
train.fun.br = function(x,t,d,tau){train.fun(x,t,d,tau,pred.models=models,v=V,sl.library=SL.library)}
predict.fun.br = function(out,x0,tau){predict.fun(out,x0,tau,pred.models=models)}
######

out.brcancer = 
  rmst.pred(x.brcancer, tobs.brcancer, delta.brcancer, tau.brcancer, 
            train.fun.br, predict.fun.br, n.folds=10, cens.model="km", 
            alpha=0.1, varsL=0, varsG=0, verbose=TRUE, seed=20,
            error=T,roo=T,vimpL=T,vimpG=T,rho=0.5)

### MSE ###

pdf(file="fig/mse_brcancer.pdf",w=3.5,h=3.5)
par(mfrow=c(1,1))
plot(out.brcancer, model.names = models.names, elements = c("mse"))
graphics.off()


### Local variable importance ###

pdf(file="fig/vimpL_brcancer.pdf",w=9,h=6)
par(mfrow=c(m,5))
plot(out.brcancer, model.names = models.names, elements = c("vimpL"),
     varsL=c(2,4,5,6,7))
graphics.off()

### Global variable importance ###

## Test ##

sink(file = "fig/test_brcancer.txt", type = "output")
print(out.brcancer, model.names = models.names, var.names = colnames(x.brcancer))
sink()

## Display ##

pdf(file="fig/vimpG_brcancer.pdf",w=12,h=3)
par(mfrow=c(1,m))
plot(out.brcancer, model.names = models.names, elements = c("vimpG"))
graphics.off()


### Multiple splitting ###

multiple.splits.br.pll = function(X){
  set.seed(X)
  out = loco.surv(x.brcancer, tobs.brcancer, delta.brcancer, tau.brcancer, 
                  train.fun.br, predict.fun.br, vars=0,rho=0.5)
  return(sapply(1:length(out$inf.sign), 
                FUN = function(l){out$inf.sign[[l]][,1]}))
}

N = 40

parallelisation_init()
clusterExport(cl,c("simulate","models","V","SL.library","integral",
                   "train.fun","predict.fun","train.fun.br","predict.fun.br",
                   "x.brcancer","tobs.brcancer","delta.brcancer","tau.brcancer"))

multiple.brcancer = parLapply(cl, X = 1:N, fun = multiple.splits.br.pll)

multiple.brcancer = array(unlist(multiple.brcancer), dim = c(p,m,N))

stopCluster(cl)



pval.twice.med.br = round(apply(multiple.brcancer,MARGIN=2,
                                function(x){apply(x,1,function(y){pmin(2*median(y),1)})}),3)

signif.code = function(v) {
  code = rep("",length(v))
  code[v < 0.1] = "."
  code[v < 0.05] = "*"
  code[v < 0.01] = "**"
  code[v < 0.001] = "***"
  return(code)
}

tab.pval.twice.med.br = colnames(x.brcancer)
for(l in 1:m){
  tab = round(pval.twice.med.br[,l],digits=3)
  code = signif.code(pval.twice.med.br[,l])
  tab = cbind(tab,code)
  colnames(tab) = c("P-value",models.names[l])
  tab.pval.twice.med.br = cbind(tab.pval.twice.med.br,tab)
}
colnames(tab.pval.twice.med.br)[1] = "Var"

sink(file = "fig/multi_test_brcancer.txt", type = "output")
print(tab.pval.twice.med.br,quote=F)
print("Significance codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",quote=F)
sink()


########
# PEAK #
########

# Events
data(peakVO2, package = "randomForestSRC")
tobs.peak = peakVO2$ttodead
delta.peak = peakVO2$died
tau.peak = quantile(tobs.peak,0.9,names=F)
x.peak = peakVO2 %>% select(-ttodead,-died)
p = dim(x.peak)[2]


#####
models = c("cox","rsf","sl")
models.names = c("Cox","RSF","P.obs.+SL")
SL.library = c("SL.lm", "SL.glmnet", "SL.gam", "SL.ranger", "SL.nnet")
V = 6
m = length(models.names)
train.fun.peak = function(x,t,d,tau){train.fun(x,t,d,tau,pred.models=models,v=V,sl.library=SL.library)}
predict.fun.peak = function(out,x0,tau){predict.fun(out,x0,tau,pred.models=models)}
######

out.peak = 
  rmst.pred(x.peak, tobs.peak, delta.peak, tau.peak, 
            train.fun.peak, predict.fun.peak, n.folds=10, cens.model="km", 
            alpha=0.1, varsL=0, varsG=0, verbose=TRUE, seed=20,
            error=T,roo=T,vimpL=T,vimpG=T,rho=0.5)

### MSE ###

pdf(file="fig/mse_peak.pdf",w=3.5,h=3.5)
par(mfrow=c(1,1))
plot(out.peak, model.names = models.names, elements = c("mse"))
graphics.off()


### Local variable importance ###

pdf(file="fig/vimpL_peak.pdf",w=25,h=6)
par(mfrow=c(m,13))
plot(out.peak, model.names = models.names, elements = c("vimpL"),
     varsL=c(1,22,23,26,28,29,30,31,33,34,35,36,39))
graphics.off()

### Global variable importance ###

## Test ##

sink(file = "fig/test_peak.txt", type = "output")
print(out.peak, model.names = models.names, var.names = colnames(x.peak))
sink()

## Display ##

pdf(file="fig/vimpG_peak.pdf",w=12,h=3)
par(mfrow=c(1,m))
plot(out.peak, model.names = models.names, elements = c("vimpG"))
graphics.off()


### Multiple splitting ###

multiple.splits.peak.pll = function(X){
  set.seed(X)
  out = loco.surv(x.peak, tobs.peak, delta.peak, tau.peak, 
                  train.fun.peak, predict.fun.peak, vars=0,rho=0.5)
  return(sapply(1:length(out$inf.sign), 
                FUN = function(l){out$inf.sign[[l]][,1]}))
}

N = 40

parallelisation_init()
clusterExport(cl,c("simulate","models","V","SL.library","integral",
                   "train.fun","predict.fun","train.fun.peak","predict.fun.peak",
                   "x.peak","tobs.peak","delta.peak","tau.peak"))

multiple.peak = parLapply(cl, X = 1:N, fun = multiple.splits.br.pll)

multiple.peak = array(unlist(multiple.peak), dim = c(p,m,N))

stopCluster(cl)



pval.twice.med.peak = round(apply(multiple.peak,MARGIN=2,
                                function(x){apply(x,1,function(y){pmin(2*median(y),1)})}),3)


tab.pval.twice.med.peak = colnames(x.peak)
for(l in 1:m){
  tab = round(pval.twice.med.peak[,l],digits=3)
  code = signif.code(pval.twice.med.peak[,l])
  tab = cbind(tab,code)
  colnames(tab) = c("P-value",models.names[l])
  tab.pval.twice.med.peak = cbind(tab.pval.twice.med.peak,tab)
}
colnames(tab.pval.twice.med.peak)[1] = "Var"

sink(file = "fig/multi_test_peak.txt", type = "output")
print(tab.pval.twice.med.peak,quote=F)
print("Significance codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",quote=F)
sink()


