library(SamplingStrata)
data("swissmunicipalities")
swissmun <- swissmunicipalities[,c("REG","COM","Nom","HApoly",
                                   "Surfacesbois","Surfacescult",
                                   "Airbat","POPTOT")]
swissmun$HApoly.cat <- var.bin(swissmun$HApoly,15)
table(swissmun$HApoly.cat)
swissmun$POPTOT.cat <- var.bin(swissmun$POPTOT,15,iter.max = 500)
table(swissmun$POPTOT.cat)
ndom <- length(unique(swissmun$REG))
frame <- buildFrameDF(df = swissmun,
                      id = "COM",

                      X = c("POPTOT.cat","HApoly.cat"),
                      Y = c("Airbat","Surfacesbois"),
                      domainvalue = "REG")
frame$domainvalue <-as.numeric(factor(frame$domainvalue, levels = unique(frame$domainvalue)))
head(frame<-(as.data.frame(frame)))
cv <- as.data.frame(list(DOM=rep("DOM1",ndom),
                         CV1=rep(0.10,ndom),
                         CV2=rep(0.10,ndom),
                         domainvalue=c(1:ndom) ))

cv
set.seed(1234)
strata <- buildStrataDF(frame, progress=F)
set.seed(1234)
kmean <- KmeansSolution(strata,errors=cv,maxclusters = 10)
nstrata <- tapply(kmean$suggestions,
                  kmean$domainvalue,
                  FUN=function(x) length(unique(x)))
nstrata

library(mlrMBO)
cv_experiment<-cv
HyperparametersFunction<-function(x){
  x1 <- data.frame(sapply(x, function(x) as.numeric(as.character(x))))
  jsize <-round(x1[1,],0)
  length_of_markov_chain  <-round(x1[2,],0)
  Temp<-x1[3,]
  decrement_constant<-x1[4,]
  kmax_percent <- x1[5,]
  ProbNewStratum <-x1[6,]
  cat("jsize", jsize,"\n")
  cat("length_of_markov_chain", length_of_markov_chain,"\n")
  cat("Temp", Temp,"\n")
  cat("decrement_constant", decrement_constant,"\n")
  cat("kmax_percent", kmax_percent,"\n")
  cat("ProbNewStratum", ProbNewStratum,"\n")

  solution <-  parallelSAA(strata=strata,errors=cv_experiment,
                                             sugg=kmean ,
                                             Temp=Temp,initialStrata=nstrata, decrement_constant=decrement_constant,
                                             end_time =Inf,
                                             showSettings = FALSE, jsize,length_of_markov_chain,
                                             verbose = FALSE,dominio=NULL,minnumstrat=2,
                                             kmax=kmax_percent,ProbNewStratum=ProbNewStratum,
                                             strcens=FALSE,writeFiles=FALSE, showPlot=TRUE,
                                             minTemp = 0.000000000001, realAllocation=TRUE)
  cat("sample size", sum(solution$best), "\n")
  return(sum(solution$best))
}
nvar <- length(grep("CV", names(cv)))
my_data <- strata[,3:(2+nvar)]
objfun2 = makeSingleObjectiveFunction(
  name = "mixed_example",
  fn = HyperparametersFunction,
  par.set = makeParamSet(
    makeDiscreteParam("jsize", values = seq(10,50,10)),
    makeDiscreteParam("length_of_markov_chain", values =seq(100,300,100)),
    makeNumericParam("Temp", lower = 0,upper = 0.001),
    makeNumericParam("decrement_constant", lower = 0.5,upper = 1),
    makeNumericParam("kmax_percent", lower = 0.0001,upper = 0.025),
    makeNumericParam("ProbNewStratum", lower = 0,upper = 0.1)
  ),
  noisy=TRUE,
  has.simple.signature = FALSE,
  minimize = TRUE
)



#objfun2(design2[1,])
surr.rf = makeLearner("regr.randomForest", predict.type = "se")

control2 = makeMBOControl()
control2 = setMBOControlInfill(
  control = control2,
  crit = makeMBOInfillCritCB(cb.lambda = 5),
  opt.focussearch.points =1000
)

control2 = setMBOControlTermination(
  control = control2,
  iters = 10
)

design2 = generateDesign(n = 10, par.set = getParamSet(objfun2))
#apply(as.matrix(design2),1,objfun2)
mlr::configureMlr(show.info = FALSE, show.learner.output = FALSE, on.learner.warning = "quiet")
ptm <- proc.time()
run2 = mbo(objfun2, design = design2, learner = surr.rf, control = control2, show.info = TRUE)
proc.time() - ptm
library(readr)


direnew<-getwd()
for(i in 1:length(run2)){
  write.csv(data.frame(run2[[i]]), file = paste0(direnew,'/', names(run2)[i], '.csv'))
}
fileConn<-file("final_opt_state.txt")
writeLines(run2$final.opt.state, fileConn)
close(fileConn)
zz <- file("final_opt_state.Rout", open="wt")
sink(zz)
sink(zz, type="message")
run2$final.opt.state
## back to the console
sink(type="message")
sink()
run2$final.opt.state
sink()
sink()

opt_path <- read_csv("opt.path.csv")
SATuningTime<- sum(opt_path$exec.time)

zz <- file("control.Rout", open="wt")
sink(zz)
sink(zz, type="message")
run2$control
## back to the console
sink(type="message")
sink()
run2$control
sink()

zz <- file("models.Rout", open="wt")
sink(zz)
sink(zz, type="message")
run2$models
## back to the console
sink(type="message")
sink()
run2$models
sink()

zz <- file("final_State.Rout", open="wt")
sink(zz)
sink(zz, type="message")
run2$final.state
## back to the console
sink(type="message")
sink()
run2$final.state
sink()

write(SATuningTime,"SATuningTime.csv")

