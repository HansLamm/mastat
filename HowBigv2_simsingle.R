#!/usr/bin/env Rscript
#model: ongoing+amount+pop.growth
#v2: no intervention, but 2 versions of transmission parameter (in hazard with viral load -> different rates of HIV transmission)
#introduction HIV at year 10, continue till t=intro + 25y -> sim till year 35 
#100 runs / store incidence, prevalence, checks (population growth)
#
args <- commandArgs(TRUE)
#parameters 1:population, 2:agegap.mean, 3:agegap.sd, 4:eyecap.amount, 5:nsim, 6:computime, 7:formation.hazard.agegap.gap_factor, 8:formation.hazard.agegap.baseline, 9:dissolution.alpha_0, 10:dissolution.alpha_4
library(RSimpactCyan)
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
library(exactci)
library(magrittr)
library(readcsvcolumns)


#code Wim Delva
#included from https://github.com/wdelva/RSimpactHelp/tree/master/R:

#' Read the Simpact output files.
#'
#' Read the .csv files and combine them into one list.

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

#' Calculate population growth rate.
#'
#' Calculate the growth rate of the entire population, averaged over a particular time window.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound) that should be retained, e.g. c(20, 30)
#' @return population growth rate estimate, for the specified time window
#' @examples
#' growth.rate <- pop.growth.calculator(datalist = datalist, timewindow = c(0, 20))

pop.growth.calculator <- function(datalist = datalist,
                                  timewindow = c(0, 20)){
  end.popsize <- datalist$ltable %>% subset(Time==timewindow[2]) %>% select(PopSize) %>% as.numeric()
  start.popsize <- ifelse(timewindow[1]==0,
                          yes = datalist$itable$population.nummen + datalist$itable$population.numwomen,
                          no = datalist$ltable %>% subset(Time==timewindow[1]) %>% select(PopSize) %>% as.numeric())
  growth.rate <- log(end.popsize/start.popsize) / diff(timewindow)
  return(growth.rate)
}

#' Subset the datalist ptable to people alive at a point in time and add HIV status.
#'
#' Subset the datalist ptable to include only people who were alive at a point in time and add their HIV status at that point in time.
#'
#' @param DT The person table (ptable) that is produced by \code{\link{readthedata()}}
#' @param timepoint Point in time at which the subset should be created and HIV status should be evaluated.
#' @return a data.table that only includes people who were alive at the timepoint and that records their HIV status.
#' @examples
#' alive.twenty.dt <- alive.infected(DT = datalist$ptable, timepoint = 20)

alive.infected <- function(DT = datalist$ptable, timepoint = 20){ # arguments are the personlog data.table and a point in time
  DTalive <- subset(DT, TOB <= timepoint & TOD > timepoint)
  DTalive$Infected <- timepoint >= DTalive$InfectTime # Now we allocate infection status to all people in our table of living people
  return(DTalive)
}


#' Calculate HIV incidence, overall and stratified.
#'
#' Calculate the HIV incidence in a time window and for specific age groups and gender strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timewindow Boundaries of the time window (lower bound < time <= upper bound) that should be retained, e.g. c(20, 30)
#' @return a dataframe with cases, exposure time, incidence estimate and surrounding confidence bounds,
#' for the specified time window and age group, overall, and stratified by gender
#' @examples
#' incidence.df <- incidence.calculator(datalist = datalist, agegroup = c(15, 30), timewindow = c(20, 30))

incidence.calculator <- function(datalist = datalist,
                                 agegroup = c(15, 30),
                                 timewindow = c(20, 30)){
  time.of.lowerbound.agegroup <- datalist$ptable$TOB + agegroup[1]
  time.of.lowerbound.timewind <- timewindow[1]
  exposure.start <- pmax(time.of.lowerbound.agegroup, time.of.lowerbound.timewind)
  time.of.upperbound.agegroup <- datalist$ptable$TOB + agegroup[2]
  time.of.upperbound.timewind <- timewindow[2]
  time.of.HIV.infection <- datalist$ptable$InfectTime
  exposure.end <- pmin(time.of.HIV.infection, pmin(time.of.upperbound.agegroup, time.of.upperbound.timewind))
  exposure.time <- exposure.end - exposure.start # This is the naive exposure time, before tidying up
  real.exposure.time <- exposure.time > 0 # We create a vector to see who REALLY had exposure time
  exposure.time[real.exposure.time == FALSE] <- 0
  
  # Now we check who of the people with the real exposure time had the event
  # Their InfectTime must be after their exposure.time started, and before or at exposure.end
  infection.after.exposure.start <- datalist$ptable$InfectTime > exposure.start
  infection.before.or.at.exposure.end <- datalist$ptable$InfectTime <= exposure.end
  infection.in.timewindow <- infection.after.exposure.start & infection.before.or.at.exposure.end
  
  
  datalist$ptable$incident.case <- infection.in.timewindow
  datalist$ptable$exposure.times <- exposure.time
  
  raw.df <- data.frame(datalist$ptable)
  
  # Now we apply some dplyr function to get the sum of cases and sum of exposure.time per gender.
  incidence.df <- dplyr::summarise(group_by(raw.df, Gender),
                                   sum.exposure.time = sum(exposure.times),
                                   sum.incident.cases = sum(incident.case),
                                   incidence = sum(incident.case) / sum(exposure.times),
                                   incidence.95.ll = as.numeric(poisson.exact(x = sum(incident.case), T = sum(exposure.times))$conf.int)[1],
                                   incidence.95.ul = as.numeric(poisson.exact(x = sum(incident.case), T = sum(exposure.times))$conf.int)[2]
  )
  
  # Now we add the overall incidence to this dataframe
  incidence.all.df <- cbind(Gender = NA,
                            dplyr::summarise(raw.df,
                                             sum.exposure.time = sum(exposure.times),
                                             sum.incident.cases = sum(incident.case),
                                             incidence = sum(incident.case) / sum(exposure.times),
                                             incidence.95.ll = as.numeric(poisson.exact(x = sum(incident.case), T = sum(exposure.times))$conf.int)[1],
                                             incidence.95.ul = as.numeric(poisson.exact(x = sum(incident.case), T = sum(exposure.times))$conf.int)[2]
                            ))
  incidence.df <- rbind(incidence.df, incidence.all.df)
  return(incidence.df)
}

# Incidence is calculated in 3 steps.
# 1. Calculate PY of exposure per person
# 2. Calculate whether the person had the event or not
# 3. Divide events by sum of PY.

# 1. Calculation of PY of exposure per person. Exposure time starts at the max
# of the lower bound of the timewindow and the time at which the person reaches
# the lower bound of the age group. Exposure time ends at the min of the upper
# bound of the timewindow, the time at which the person reaches the upper bound
# of the age group, and the time at which the person gets infected with HIV.
# Any negative exposure time will be reset to zero.

#########
#' Calculate HIV prevalence, overall and stratified.
#'
#' Calculate the HIV prevalence at a point in time, for specific age groups and gender strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param timepoint Point in time at which the HIV prevalence should be calculated.
#' @return a dataframe with prevalence estimate and surrounding confidence bounds,
#' for the specified time point and age group, overall, and stratified by gender
#' @examples
#' prevalence.df <- prevalence.calculator(datalist = datalist, agegroup = c(15, 30), timepoint = 30)

prevalence.calculator <- function(datalist = datalist,
                                  agegroup = c(15, 30),
                                  timepoint = 30){
  DTP <- datalist$ptable
  DTalive.infected <- alive.infected(DT = DTP, timepoint = timepoint) # First we only take the data of people who were alive at the timepoint
  DTalive.infected.agegroup <- subset(DTalive.infected, TOB <= timepoint - agegroup[1] & TOB > timepoint - agegroup[2])
  raw.df <- data.frame(DTalive.infected.agegroup)
  # Now we apply some dplyr function to get the sum of cases and sum of exposure.time per gender.
  prevalence.df <- dplyr::summarise(group_by(raw.df, Gender),
                                    popsize = length(Gender),
                                    sum.cases = sum(Infected),
                                    pointprevalence = sum(Infected) / length(Gender),
                                    pointprevalence.95.ll = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[1],
                                    pointprevalence.95.ul = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[2]
  )
  prevalence.all.df <- cbind(Gender = NA,
                             dplyr::summarise(raw.df,
                                              popsize = length(Gender),
                                              sum.cases = sum(Infected),
                                              pointprevalence = sum(Infected) / length(Gender),
                                              pointprevalence.95.ll = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[1],
                                              pointprevalence.95.ul = as.numeric(binom.test(x = sum(Infected), n = length(Gender))$conf.int)[2]
                             )
  )
  prevalence.df <- rbind(prevalence.df, prevalence.all.df)
  return(prevalence.df)
}

##code Wil Delva till here



###########
# simulation

# single run on HPC
ABC_DestDir <- Sys.getenv("VSC_SCRATCH_NODE")
ABC_LogDir <- Sys.getenv("VSC_DATA")

#single run parameters - from data.csv-file
#ABC-result: $popsize $agegapmean $agegapsd $eyecap $nsim $computime $agegap_factor $agegap_baseline $dissalpha0 $dissalpha4
#popsize, agegapmean, agegapsd, eyecap / transmission (0 of 1)
population_size <- as.double(args[1])
meanagegap <- as.double(args[2])
agegap_sd <- as.double(args[3])
eyecap_amount<-as.double(args[4])
transmission_type<-as.double(args[11])

modelname <-paste0(population_size, "-", meanagegap, "-", agegap_sd, "-", eyecap_amount, "-t", transmission_type)
DestDir <- paste0(ABC_DestDir, "/", modelname)
LogDir <- ABC_LogDir

##intervention scenarios
#1. no intervention: CD4 threshold=0 -> initial situation remains

  modelmetrics <- list()
  modelreps<-100 #100 simulations of each parameter combination
  modelmetrics[["inc1"]] <- vector(mode = "list", length = modelreps)
  modelmetrics[["inc5"]] <- vector(mode = "list", length = modelreps)
  modelmetrics[["prev"]] <- vector(mode = "list", length = modelreps)
  modelmetrics[["checks"]] <- vector(mode = "list", length = modelreps)
  modelrepvect <- 1:modelreps
  
  ##cfg settings
  cfg <- list()
  cfgHowBig <- simpact.getconfig(cfg)
  cfgHowBig["population.numwomen"] <- population_size
  cfgHowBig["population.nummen"] <- population_size
  cfgHowBig["person.agegap.man.dist.type"] <- "fixed"
  cfgHowBig["person.agegap.woman.dist.type"] <- "fixed"
  cfgHowBig["person.agegap.man.dist.fixed.value"] <- meanagegap
  cfgHowBig["person.agegap.woman.dist.fixed.value"] <- meanagegap
  cfgHowBig["dissolution.Dp"] <- meanagegap
  
  #estimated parameters (ABC best fit):
  cfgHowBig["formation.hazard.agegap.gap_factor_man"] <- as.double(args[7])
  cfgHowBig["formation.hazard.agegap.gap_factor_woman"] <- as.double(args[7])
  cfgHowBig["formation.hazard.agegap.baseline"] <- as.double(args[8])
  cfgHowBig["dissolution.alpha_0"] <- as.double(args[9])  #baseline
  cfgHowBig["dissolution.alpha_4"] <- as.double(args[10])  #weight age
  
  #eyecap: based on eyecap_amount (cfr Dunbar's number 300 / 1000 for larger value)
  eyecap_fraction<-eyecap_amount/(population_size*2)
  if (eyecap_fraction > 1) eyecap_fraction = 1 
  cfgHowBig["population.eyecap.fraction"] <- eyecap_fraction
  
  #simtime
  simtime<-35 #simulation time = intro HIV + 25y
  cfgHowBig["population.simtime"] <- simtime
  cfgHowBig["periodiclogging.interval"] <- 1

  #introduction HIV
  cfgHowBig["hivseed.time"] <- 10
  cfgHowBig["hivseed.type"] <- "fraction"
  cfgHowBig["hivseed.fraction"] <- 0.05 #5% of agegroup 20-30 years -> appr 1% sexually active population
  cfgHowBig["hivseed.age.min"] <- 20
  cfgHowBig["hivseed.age.max"] <- 30
  #to make sure 1% HIV is present in small populations
  if (population_size <= 100){
    cfgHowBig$hivseed.fraction <- NULL
    cfgHowBig["hivseed.type"] <- "amount"
    ifelse(population_size < 100, cfgHowBig["hivseed.amount"]<-1, cfgHowBig["hivseed.amount"]<-2)
  }
  cfgHowBig["monitoring.cd4.threshold"] <-0 #initial situation

  ##transmission param
  #transmission_type=0 -> default values
  #transmission_type=1 -> new values
  if (transmission_type==1){
    cfgHowBig["transmission.param.a"] <- -1
    cfgHowBig["transmission.param.b"] <- -90
    cfgHowBig["transmission.param.c"] <- 0.5
  }
  
  
  for (modelrep in modelrepvect){
    simID <- paste0(modelname, "-", modelrep)
    identifier <- paste0("%T-%y-%m-%d-%H-%M-%S-", simID, "-")
    seedid <- 1000 + modelrep
    
    fullsim <- simpact.run(cfgHowBig,
                           DestDir,
                           identifierFormat = identifier,
                           seed = seedid)
       
    datalist <- readthedata(modeloutput = fullsim)
    
    ##toevoegen incidence + prevalence data

    #incidence -> rate at end simulation (over last 5 or 10 years)
    inc10 <- incidence.calculator(datalist=datalist, agegroup = c(15, 70), timewindow = c(simtime - 10, simtime))
    modelmetrics[["inc10"]][[modelrep]] <- inc10

    #evolution incidence rate 1: from year 15 till end sim over last 5 years
    inc5.list <- vector("list", simtime-14)
    for (t in seq(15, simtime)) {
      incidence <- incidence.calculator(datalist=datalist, agegroup = c(15, 70), timewindow = c(t - 5, t))
      inc5.list[t-14] <- incidence$incidence[3] #overall incidence
    }
    modelmetrics[["inc5"]][[modelrep]] <- inc5.list
    
    #evolution incidence rate 2: from year 11 till end sim over last year
    inc1.list <- vector("list", simtime-10)
    for (t in seq(11, simtime)) {
      incidence <- incidence.calculator(datalist=datalist, agegroup = c(15, 70), timewindow = c(t - 1, t))
      inc1.list[t-10] <- incidence$incidence[3] #overall incidence
    }
    modelmetrics[["inc1"]][[modelrep]] <- inc1.list

    #evolution prevalence rate: from year 10 till end sim
    prev.list <- vector("list", simtime-9)
    for (t in seq(10, simtime)) {
      prev <- prevalence.calculator(datalist=datalist, agegroup = c(15, 70), timepoint = t)
      prev.list[t-9] <- prev$pointprevalence[3] #overall prevalence at time t
    }
    modelmetrics[["prev"]][[modelrep]] <- prev.list
    
    #checks (population growth, agegap mean + sd, relation targets)
    #agegap
    model.agegap.mean <- mean(datalist$rtable$AgeGap)
    model.agegap.sd <- sd(datalist$rtable$AgeGap)
    # Number of relationships formed in the target window: similar to parameter search, as the epidemy can cause deviation at later stages
    targetwindow.start <- 10
    targetwindow.end <- 20
    targetwindow.duration <- targetwindow.end - targetwindow.start
    targetwindow.midpoint <- targetwindow.start + targetwindow.duration/2
    targetwindow.rels <- subset(datalist$rtable, FormTime >=targetwindow.start & FormTime < targetwindow.end)
    nrel <- nrow(targetwindow.rels)
    halfnpeople <- nrow(subset(datalist$ptable, TOB <= targetwindow.midpoint & TOD > targetwindow.midpoint))/2
    model.newrelationsall.average <- nrel / halfnpeople / targetwindow.duration
    #population growth -> at targetwindow to avoid strong influence of epidemy
    model.popgrowth1 <- pop.growth.calculator(datalist = datalist, timewindow = c(targetwindow.start, targetwindow.end))
    model.popgrowth2 <- pop.growth.calculator(datalist = datalist, timewindow = c(targetwindow.start, targetwindow.midpoint))
    #ongoing relations -> at end simulation as behaviour needs to remain constant
    personsalive <- datalist$ptable[is.infinite(TOD), .(ID, Gender, TOB, agegroup = (2 + 5*floor((simtime - TOB)/5)))]
    ongoingrelations <-datalist$rtable[DisTime==Inf]
    if (!empty(ongoingrelations)) { 
      personsinrelation <- rbind(ongoingrelations[, .(ID = IDm)], ongoingrelations[, .(ID = IDw)])
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
    model.ongoingrelations.logitintercept <- summary(inrelationfit)$coefficients["(Intercept)","Estimate"]
    model.ongoingrelations.logitslope <- summary(inrelationfit)$coefficients["agegroup","Estimate"]
    modelmetrics[["checks"]][[modelrep]] <- list(model.popgrowth1,
                                                 model.popgrowth2,
                                                 model.agegap.mean,
                                                 model.agegap.sd,
                                                 model.newrelationsall.average,
                                                 model.ongoingrelations.logitslope,
                                                 model.ongoingrelations.logitintercept)
    
  }
  
  
  #storage results
  filename_results <- paste0(LogDir, "/SimResult-", population_size, "-", meanagegap, "-", agegap_sd, "-", eyecap_amount, "-t", transmission_type, ".rds")
  saveRDS(modelmetrics, filename_results)
  
  filename_inc1 <- paste0(LogDir, "/inc1-", population_size, "-", meanagegap, "-", agegap_sd, "-", eyecap_amount, "-t", transmission_type, ".csv")
  filename_inc5 <- paste0(LogDir, "/inc5-", population_size, "-", meanagegap, "-", agegap_sd, "-", eyecap_amount, "-t", transmission_type, ".csv")
  filename_inc10 <- paste0(LogDir, "/inc10-", population_size, "-", meanagegap, "-", agegap_sd, "-", eyecap_amount, "-t", transmission_type, ".csv")
  filename_prev <- paste0(LogDir, "/prev-", population_size, "-", meanagegap, "-", agegap_sd, "-", eyecap_amount, "-t", transmission_type, ".csv")
  filename_checks <- paste0(LogDir, "/checks-", population_size, "-", meanagegap, "-", agegap_sd, "-", eyecap_amount, "-t", transmission_type, ".csv")
  
  
  inc1.frame <- data.frame()
  inc5.frame <- data.frame()
  inc10.frame <- data.frame()
  prev.frame <- data.frame()
  checks.frame <- data.frame()
  
  for (modelrep in modelrepvect){
    inc5.rep <- data.frame(modelmetrics[["inc5"]][[modelrep]][1:(simtime-14)], row.names = NULL)
    colnames(inc5.rep)<-seq(15, simtime)
    inc5.frame <- rbind(inc5.frame, inc5.rep)

    inc1.rep <- data.frame(modelmetrics[["inc1"]][[modelrep]][1:(simtime-10)], row.names = NULL)
    colnames(inc1.rep)<-seq(11, simtime)
    inc1.frame <- rbind(inc1.frame, inc1.rep)
    
    rep.frame <- as.data.frame(rep(modelrep, 3))
    colnames(rep.frame)<-"modelrep"
    inc10.rep<-cbind(rep.frame, modelmetrics[["inc10"]][[modelrep]])
    inc10.frame <- rbind(inc10.frame, inc10.rep)
    
    prev.rep <- data.frame(modelmetrics[["prev"]][[modelrep]][1:(simtime-9)], row.names = NULL)
    colnames(prev.rep)<-seq(10, simtime)
    prev.frame <- rbind(prev.frame, prev.rep)
    
    checks.rep <- data.frame(modelmetrics[["checks"]][[modelrep]][1:7], row.names = NULL)
    colnames(checks.rep)<-c("popgrowth10y", "popgrowth5y", "agegap.mean", "agegap.sd", "newrelationsall.average", "ongoingrelations.logitslope", "ongoingrelations.logitintercept")
    checks.frame <- rbind(checks.frame, checks.rep)
  }
  
  write.csv(inc1.frame, filename_inc1, sep=",", row.names=FALSE)
  write.csv(inc5.frame, filename_inc5, sep=",", row.names=FALSE)
  write.csv(inc10.frame, filename_inc10, sep=",", row.names=FALSE)
  write.csv(prev.frame, filename_prev, sep=",", row.names=FALSE)
  write.csv(checks.frame, filename_checks, sep=",", row.names=FALSE)


