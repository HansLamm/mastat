#!/usr/bin/env Rscript
#single run HowBig simpact parameter search with ABC (parallel version)
args <- commandArgs(TRUE)

#parameters 1: amount cores, 2: population size, 3: mean agegap, 4: sd agegap, 5: eyecap amount, 6: seed
ncluster <- as.double(args[1])
population_size <- as.double(args[2])
meanagegap <- as.double(args[3])
agegapsd <- as.double(args[4])
eyecap_amount <- as.double(args[5])
seed <- as.double(args[6])
  
ABC_LogDir <- Sys.getenv("VSC_DATA")

library(EasyABC)
library(data.table)
library(nlme)
library(ggplot2)
library(RColorBrewer)
library(shape)
library(reshape2)
library(MASS)
library(grid)
library(stringr)
library(plyr)
library(pscl)
library(dplyr)
library(magrittr)
library(RSimpactCyan)
library(readcsvcolumns)

#simpact wrapper
simpact4ABC.wrapper <- function(inputvector){  
    library(data.table)
    library(nlme)
    library(ggplot2)
    library(RColorBrewer)
    library(shape)
    library(reshape2)
    library(MASS)
    library(grid)
    library(stringr)
    library(plyr)
    library(pscl)
    library(dplyr)
    library(magrittr)
    library(RSimpactCyan)
    library(readcsvcolumns)
    
    #code Wim Delva    
     readthedata <- function(modeloutput){
      path <- as.character(modeloutput["outputfile"])
      outputID <- as.character(modeloutput["id"])
      DestDir <- sub(pattern = paste0(outputID, "output.txt"), replacement = "", x = path, fixed=T)
      personlogfilename <- paste0(DestDir, outputID, "personlog.csv")
      relationlogfilename <- paste0(DestDir, outputID, "relationlog.csv")
      eventlogfilename <- paste0(DestDir, outputID, "eventlog.csv")
      treatmentlogfilename <- paste0(DestDir, outputID, "treatmentlog.csv")
      inputparamlogfilename <- paste0(DestDir, outputID, "settingslog.csv")
      periodiclogfilename <- paste0(DestDir, outputID, "periodiclog.csv")
      
      ptable <- data.table::fread(personlogfilename, sep = ",", skip = 0)
      rtable <- data.table::fread(relationlogfilename, sep = ",", skip = 0)
      readetable <- readcsvcolumns::read.csv.columns(eventlogfilename, has.header = FALSE, column.types = "rssiirsiir")
      etable <- data.table::setDT(readetable)
      etable.colnames <- colnames(etable)
      if (ncol(etable) == 10){
        data.table::setnames(etable, etable.colnames,
                             c("eventtime", "eventname", "p1name", "p1ID", "p1gender", "p1age", "p2name", "p2ID", "p2gender", "p2age"))
      } else {
        data.table::setnames(etable, etable.colnames,
                             c("eventtime", "eventname", "p1name", "p1ID", "p1gender", "p1age", "p2name", "p2ID", "p2gender", "p2age", "extradescript", "value"))
        
      }
      ttable <- data.table::fread(treatmentlogfilename, sep  = ",", skip = 0)
      itable <- data.table::fread(inputparamlogfilename, sep  = ",", skip = 0)
      
      if (file.exists(periodiclogfilename)){
        ltable <- data.table::fread(periodiclogfilename, sep = ",", skip = 0)
        outputtables <- list(ptable = ptable, rtable = rtable, etable = etable, ttable = ttable, itable = itable, ltable = ltable)
      } else {
        outputtables <- list(ptable = ptable, rtable = rtable, etable = etable, ttable = ttable, itable = itable)
      }
      return(outputtables)
     }
    
    #code Wim Delva   
    #Calculate population growth rate.
    pop.growth.calculator <- function(datalist = datalist,
                                      timewindow = c(0, 20)){
      end.popsize <- datalist$ltable %>% subset(Time==timewindow[2]) %>% select(PopSize) %>% as.numeric()
      start.popsize <- ifelse(timewindow[1]==0,
                              yes = datalist$itable$population.nummen + datalist$itable$population.numwomen,
                              no = datalist$ltable %>% subset(Time==timewindow[1]) %>% select(PopSize) %>% as.numeric())
      growth.rate <- log(end.popsize/start.popsize) / diff(timewindow)
      return(growth.rate)
   }
    
    
    ########
    # ABC parameter inference for Simpact model
    
    #data input
    ABC_ScratchDir <- Sys.getenv("VSC_SCRATCH_NODE")
    source(paste0(ABC_ScratchDir, "/simpact_inputparam.R"))
    
    #location simulation data
    ABC_DestDir <- Sys.getenv("VSC_SCRATCH_NODE")
    
    ##cfg settings
    cfgHowBig <- list()
    
    #simulation and logging settings
    targetwindow.start <- 10
    targetwindow.end <- 20 
    targetwindow.duration <- targetwindow.end - targetwindow.start
    targetwindow.midpoint <- targetwindow.start + targetwindow.duration/2
    cfgHowBig["population.simtime"] <- targetwindow.end #calibrate based on the period 10-20 years 
    cfgHowBig["periodiclogging.interval"] <- 1
    
    #population size
    cfgHowBig["population.numwomen"] <- population_size
    cfgHowBig["population.nummen"] <- population_size
    
    #eyecap: based on amount (300, 1000)
    eyecap_fraction<-eyecap_amount/(population_size*2)
    if (eyecap_fraction > 1) eyecap_fraction = 1
    cfgHowBig["population.eyecap.fraction"] <- eyecap_fraction
    
    ##formation
    cfgHowBig["person.agegap.man.dist.type"] <- "fixed"
    cfgHowBig["person.agegap.woman.dist.type"] <- "fixed"
    cfgHowBig["person.agegap.man.dist.fixed.value"] <- meanagegap
    cfgHowBig["person.agegap.woman.dist.fixed.value"] <- meanagegap
    cfgHowBig["formation.hazard.type"] <- "agegap"
    
    ##no HIV
    cfgHowBig["hivseed.time"] <- 100
    
    ##intervention -> none for ABC parameter estimation
    # Nobody will start ART during the simulations
    cfgHowBig["diagnosis.baseline"] <- -100 # This will result in timing of HIV diagnosis way beyond the simulation period.
    
    #parameters delivered through the ABC inputvector
    seedid <- inputvector[1]
    cfgHowBig["formation.hazard.agegap.gap_factor_man"] <- inputvector[2]
    cfgHowBig["formation.hazard.agegap.gap_factor_woman"] <- inputvector[2]
    cfgHowBig["formation.hazard.agegap.baseline"] <- inputvector[3]
    cfgHowBig["dissolution.alpha_0"] <- inputvector[4]  #baseline
    cfgHowBig["dissolution.alpha_4"] <- inputvector[5]  #weight age
    
    results <- simpact.run(cfgHowBig, ABC_DestDir, seed=seedid)
    datalist <- readthedata(results)
    
    #calculation summary statistics
    #sd agegap
    agegapsd <- sd(datalist$rtable$AgeGap)
    
    #Number of relationships formed per year
      targetwindow.rels <- subset(datalist$rtable, FormTime >=targetwindow.start & FormTime < targetwindow.end)
      nrel <- nrow(targetwindow.rels)
      halfnpeople <- nrow(subset(datalist$ptable, TOB <= targetwindow.midpoint & TOD > targetwindow.midpoint))/2
      relsperpersonperyear <- nrel / halfnpeople / targetwindow.duration
    
    #ongoing relations
      personsalive <- datalist$ptable[is.infinite(TOD), .(ID, Gender, TOB, agegroup = (2 + 5*floor((targetwindow.end - TOB)/5)))]
    #persons alive (for no ongoing-branch)
      ongoingrelations <-datalist$rtable[DisTime==Inf]
      if (!empty(ongoingrelations)) { 
        personsinrelation <- rbind(ongoingrelations[, .(ID = IDm)], ongoingrelations[, .(ID = IDw)])  #simpact 19.6 (HPC)
        #personsinrelation <- rbind(ongoingrelations[, .(ID = ID1)], ongoingrelations[, .(ID = ID2)])   #simpact 20.0
        setkey(personsinrelation, ID)
        personsinrelation <- unique(personsinrelation, by="ID") #remove doubles
        personsinrelation[, ':=' (relation = 1L)] #add relation status
        inrelation <- merge(personsalive, personsinrelation, by = "ID", all = TRUE)
        inrelation[is.na(relation), relation :=0L]
        inrelation1 <- inrelation[relation==1L, .N, .(agegroup)]
        inrelation11 <- inrelation1[, .(agegroup, yes=N, no=0L)]
        inrelation0 <- inrelation[relation==0L, .N, .(agegroup)]
        inrelation00 <- inrelation0[, .(agegroup, yes=0L, no=N)]
        inrelation_grouped <- rbind(inrelation11, inrelation00)
      } else {
        age_grouped <- personsalive[, .N, .(agegroup)]
        inrelation_grouped <- age_grouped[, .(agegroup, yes=N, no=0L)]
      }
      inrelationfit <- glm(cbind(yes, no) ~ agegroup, family = binomial, data=inrelation_grouped)
      ongoingrelations.logitintercept <- summary(inrelationfit)$coefficients["(Intercept)","Estimate"]
      ongoingrelations.logitslope <- summary(inrelationfit)$coefficients["agegroup","Estimate"]
    
    #population growth - double target
    pop.growth1 <- pop.growth.calculator(datalist = datalist, timewindow = c(targetwindow.start, targetwindow.end))
    pop.growth2 <- pop.growth.calculator(datalist = datalist, timewindow = c(targetwindow.start, targetwindow.midpoint))
    
    outputvector <- c(agegapsd, relsperpersonperyear, ongoingrelations.logitintercept, ongoingrelations.logitslope, pop.growth1, pop.growth2)
    return(outputvector)
 }

##input-parameters for nodes
ABC_ScratchDir <- Sys.getenv("VSC_SCRATCH_NODE")
filename_inputparam <- paste0(ABC_ScratchDir, "/simpact_inputparam.R")
writeLines(c(paste0("population_size <- ", population_size), paste0("meanagegap <- ", meanagegap), paste0("eyecap_amount <- ", eyecap_amount)), filename_inputparam)

##parameter estimation with ABC

#target summary statistics
#ongoing relations: logistic regression model (intercept = -2.05 and slope = 0.076)
av.age = c(17, 22, 27, 32, 37, 42, 47, 52, 57, 62)
DT = data.table(av.age.group=av.age, yes=c(73, 183, 252, 201, 165, 102, 59, 13, 6, 0), no=c(171, 308, 202, 118, 64, 34, 14, 7, 2, 2))
lrfit <- glm(cbind(yes, no) ~ av.age.group, family = binomial, data=DT)
ongoingrelations.logitintercept.target <- lrfit$coeff["(Intercept)"]
ongoingrelations.logitslope.target <- lrfit$coeff["av.age.group"]
#amount partners
newrelationsall.average <- 0.4907407
#population growth
pop.growth.target <- 0.016
#sd agegap -> varieert [1,2,3,5]
agegapsd.target<-agegapsd
sum_stat_obs <- c(agegapsd.target, newrelationsall.average, ongoingrelations.logitintercept.target, ongoingrelations.logitslope.target, pop.growth.target, pop.growth.target)

##prior (parameters to estimate)
#1: formation gap_factor, 2: formation baseline 3: dissolution baseline, 4: dissolution age_factor
simpact_prior <- list(c("unif", -3, 0), c("unif", 0, 8), c("unif", 0, 4), c("unif", -6, 0))

#ABC settings
n_init <- 100
alpha <- 0.25 # This is the proportion of particles kept at each step
pacc <- 0.1 # This is the stopping criterion of the algorithm. The smaller, the more strict the criterion.

#ABC call to wrapper
ABC_LenormandResult <- ABC_sequential(method="Lenormand",
                                      model=simpact4ABC.wrapper,
                                      prior=simpact_prior,
                                      nb_simul=n_init,
                                      summary_stat_target=sum_stat_obs,
                                      alpha=alpha,
                                      p_acc_min=pacc,
                                      use_seed=TRUE,
                                      seed_count=seed,
                                      n_cluster=ncluster)

##results
filename_results <- paste0(ABC_LogDir, "/ABC_LenormandResult-", population_size, "-", meanagegap, "-", agegapsd.target, "-", eyecap_amount, ".rds")
saveRDS(ABC_LenormandResult, filename_results)

filename_resultlog <- paste0(ABC_LogDir, "/ABC_resultlog-", population_size, "-", meanagegap, "-", agegapsd.target, "-", eyecap_amount, ".txt")
txt_nsim <- paste0("nsim: ", ABC_LenormandResult$nsim)
txt_computime <- paste0("computime: ", ABC_LenormandResult$computime)

#bestfit met ks
library(ks)
eval1 <- seq(min(ABC_LenormandResult$param[,1]), max(ABC_LenormandResult$param[,1]), length.out=25)
eval2 <- seq(min(ABC_LenormandResult$param[,2]), max(ABC_LenormandResult$param[,2]), length.out=25)
eval3 <- seq(min(ABC_LenormandResult$param[,3]), max(ABC_LenormandResult$param[,3]), length.out=25)
eval4 <- seq(min(ABC_LenormandResult$param[,4]), max(ABC_LenormandResult$param[,4]), length.out=25)
evalps <- data.frame(eval1, eval2, eval3, eval4)
fhat.ks <- kde(ABC_LenormandResult$param, eval.points=evalps)
fhat.ks$estimate
fhat.ks.max <- which.max(fhat.ks$estimate)
bestfit.ks <- ABC_LenormandResult$param[fhat.ks.max,]
txt_bestfit <-  paste0("bestfit: ", bestfit.ks)
writeLines(c(txt_nsim, txt_computime, txt_bestfit), filename_resultlog)
