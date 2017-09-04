##changepoint+graphs sim analysis
library(data.table)
library(ggplot2)
library(reshape2)
library(grid)
library(stringr)
library(multcomp)
library(bcp)

sourcedir <- "/home/hans/HowBigHPC/simresultv2/"
prev.filelist <- list.files(path=sourcedir, pattern = "prev")
inc10.filelist <- list.files(path=sourcedir, pattern = "inc10")
inc5.filelist <- list.files(path=sourcedir, pattern = "inc5")
inc1.filelist <- list.files(path=sourcedir, pattern = "inc1-")
checks.filelist <- list.files(path=sourcedir, pattern = "checks")

graphdir <- "/home/hans/HowBigHPC/graphv2/"

##### reading simulation data for changepoint analysis - no graphs
simcompvaluesextra<-function(filename) {
  result_table <- fread(paste0(sourcedir, filename), header=TRUE, stringsAsFactors=FALSE)
  
  #parameter values
  namesplit<-str_split_fixed(filename, "-", 6)
  type<-namesplit[,1]
  population<-as.numeric(namesplit[,2])
  meanagegap<-as.numeric(namesplit[,3])
  sdagegap<-as.numeric(namesplit[,4])
  eyecap<-as.numeric(namesplit[,5])
  transmission<-as.numeric(substring(namesplit[,6], 2, 2))
  
  #simulation values 
  result35_mean <- mean(result_table$"35")
  result35_low <- sort(result_table$"35")[6]
  result35_high <- sort(result_table$"35")[95]
  result35_var <- var(result_table$"35")
  result35_reducedvar <- var(sort(result_table$"35")[6:95])  #for weighted means
  result35_IQR <- IQR(result_table$"35")   #for weighted means
  resultmaxvalue <- mean(apply(result_table,1,max))
  resultmaxvalue_var <- var(apply(result_table,1,max))
  resultYmax <- mean(apply(result_table,1,which.max))
  resultYmax_var <- var(apply(result_table,1,which.max))

  return(c(population, meanagegap, sdagegap, eyecap, transmission, result35_mean, result35_low, result35_high, result35_var, result35_reducedvar, result35_IQR, resultmaxvalue, resultmaxvalue_var, resultYmax, resultYmax_var))
}

for (type in c("inc1-", "inc5-", "prev-")) {
  filelist <- list.files(path=sourcedir, pattern = type)
  targetfilename <- paste0(graphdir, "simcomp_ext_", substring(type, 1, 4), ".csv")
  simcompdataextra <- t(sapply(filelist, simcompvaluesextra))
  colnames(simcompdataextra)<-c("population","meanagegap","sdagegap","eyecap","transmission","Y35.mean","Y35.low","Y35.high", "Y35.var", "Y35.reducedvar", "Y35.IQR", "max.value", "max.value.var", "max.Y", "max.Y.var")
  write.csv(simcompdataextra, targetfilename, row.names = FALSE)
}

####changepoint analysis

type<-"inc5"
#type<-"prev"
simcompfilename <- paste0(graphdir, "simcomp_ext_", type, ".csv")
simcomp_extra <- fread(paste0(simcompfilename), header=TRUE, stringsAsFactors=FALSE)
simdata<-simcomp_extra[transmission==1 & population>=1000,]
simtotal<-simcomp_extra[transmission==1 & population>50,]
simtotal.order<-simtotal[order(population,eyecap,meanagegap,sdagegap)]
  
##changepoint with segmented
weights<-1/simtotal$Y35.reducedvar
#weights<-1/simtotal$Y35.IQR
simmodel.tot<-lm(Y35.mean ~ population + meanagegap + sdagegap + eyecap, simtotal, weights=weights)
seg.tot<-segmented(simmodel.tot, seg.Z=~population, psi=1000)
summary(seg.tot)
#conclusion: changepoint with segmented strongly impacted by choice correction variance

##changepoint with bcp
bcp.sim<-bcp(simtotal.order$Y35.mean, simtotal.order$population)
#bcp.sim<-bcp(simtotal.order$Y35.mean, as.matrix(simtotal.order[,population:eyecap]))
plot(bcp.sim, main="Incidence")
#plot(bcp.sim, main="Prevalence")


##analysis above changepoint
simmodel<-lm(Y35.mean ~ population + meanagegap + sdagegap + eyecap, simdata)
summary(simmodel)
#check variance for weighted means
simdata<-simcomp_extra[transmission==1 & population>1000,]
summary(simdata)
#weights with reducedvar
weights<-1/simdata$Y35.reducedvar
simmodel<-lm(Y35.mean ~ population + meanagegap + sdagegap + eyecap, simdata, weights=weights)
#weights with IQR
weights<-1/simdata$Y35.IQR
simmodel<-lm(Y35.mean ~ population + meanagegap + sdagegap + eyecap, simdata, weights=weights)

##adjustment for multiple comparison
#glht->factor needed in model
simdata2<-simdata
simdata2$population<-factor(simdata2$population)
weights<-1/simdata$Y35.reducedvar
simmodel<-lm(Y35.mean ~ population + meanagegap + sdagegap + eyecap, simdata2, weights=weights)

sim.mc<-glht(simmodel, linfct = mcp(population = "Tukey")) 
summary(sim.mc)
#no significant difference between population factors

#contrasts
simdata3<-simdata
simdata4<-simdata
simdata5<-simdata
simdata3$meanagegap<-factor(simdata3$meanagegap)
simdata4$sdagegap<-factor(simdata4$sdagegap)
simdata5$eyecap<-factor(simdata5$eyecap)
weights<-1/simdata$Y35.reducedvar
simmodel3<-lm(Y35.mean ~ population + meanagegap + sdagegap + eyecap, simdata3, weights=weights)
simmodel4<-lm(Y35.mean ~ population + meanagegap + sdagegap + eyecap, simdata4, weights=weights)
simmodel5<-lm(Y35.mean ~ population + meanagegap + sdagegap + eyecap, simdata5, weights=weights)






#### graph evolution incidence/prevalence (each parameter combination + intervention scenario)
evolutiongraph <- function(filename) {
  result_table <- fread(paste0(sourcedir, filename))
  plotdata <- melt(result_table)
  plotdata$rowid <- 1:100
  plotdata$model <- 1

  ptext1 <- strsplit(filename, "[/.]")[[1]]
  ptext2 <- ptext1[length(ptext1)-1]
  ptext3 <- strsplit(ptext2, "-")[[1]]
  paramtext <- paste0("pop.size: ", ptext3[2], "\nmean age gap: ", ptext3[3], "\nsd age gap: ", ptext3[4], "\neyecap: ", ptext3[5], "\nintervention: ", ptext3[6])
  if (grepl("prev", ptext3[1])) {ylabel<-"prevalence"}
  else if (grepl("inc1", ptext3[1])) {ylabel<-"incidence1"}
  else if (grepl("inc5", ptext3[1])) {ylabel<-"incidence5"}
  ggplot(plotdata, aes(variable, value, group=factor(rowid))) + 
    geom_line(color='red') +  
    geom_smooth(aes(group = model), color='black') +
    labs(
      x = "year",
      y = ylabel) +
      annotate("text",  x=Inf, y = Inf, label = paramtext, vjust=1, hjust=1, size=3) 
  ggsave(paste0(graphdir, substr(filename, 1, nchar(filename)-4), ".eps"))
}

sapply(inc1.filelist, evolutiongraph)
sapply(inc5.filelist, evolutiongraph)
sapply(prev.filelist, evolutiongraph)


##graph overview inc/prev

meansim<-function(filename) {
  result_table <- fread(paste0(sourcedir, filename), header=TRUE, stringsAsFactors=FALSE)
  result_mean <- apply(result_table,2,mean)
  return(result_mean)
}

meangraph <- function(data, dest) {
  plotdata <- melt(data)
  plotdata$Var2<-as.character(plotdata$Var2)
  
  population<-as.numeric(str_split_fixed(plotdata$Var2, "-",3)[,2])
  meanagegap<-as.numeric(str_split_fixed(plotdata$Var2, "-",6)[,3])
  sdagegap<-as.numeric(str_split_fixed(plotdata$Var2, "-",6)[,4])
  eyecap<-as.numeric(str_split_fixed(plotdata$Var2, "-",6)[,5])
  transmission<-str_split_fixed(plotdata$Var2, "-",6)[,6]
  transmission<-str_split_fixed(transmission, "\\.",2)[,1]

  if (grepl("prev", plotdata$Var2[1])) {ylabel<-"prevalence"} 
  else if (grepl("inc1", plotdata$Var2[1])) {ylabel<-"incidence1"}
  else if (grepl("inc5", plotdata$Var2[1])) {ylabel<-"incidence5"}
  xlabel<-paste0("year","\n\nmean age gap: ", meanagegap[[1]], "  sd age gap: ", sdagegap[[1]], "  eyecap: ", eyecap[[1]], "  transmission: ", transmission[[1]])
  
  plotdata$Var2<-population
  names(plotdata)[names(plotdata) == 'Var2'] <- 'population'

  ggplot(plotdata, aes(Var1, value, group=factor(population))) + 
    geom_line(aes(color=factor(population))) +
    scale_colour_brewer(palette="Spectral") +
    labs(
      x = xlabel, 
      y = ylabel) 
  ggsave(paste0(graphdir, dest), device="eps")
}

#select for each parameter combination
for (type in c("inc1-", "inc5", "prev")) {
  for (transmission in c("t0", "t1")){
    for (eyecap in c(300, 1000)) {
      for (agegapsd in c(1,2,3,5)){
        for (meanagegap in c(0, 1, 2, 5, 10)){
          filelist <- list.files(path=sourcedir, pattern = type)
          searchstring<-paste0("-", meanagegap,"-", agegapsd,"-", eyecap,"-", transmission, ".csv")
          filename<-paste0("overview2-", substring(type, 1, 4), substr(searchstring, 1, nchar(searchstring)-4), ".eps")
          
          filelist2 <- filelist[grep(searchstring, filelist)]
          if (length(filelist2) != 0) {
            mean_overview <- sapply(filelist2, meansim)
            meangraph(mean_overview, filename)
          }
        }
      }
    }
  }
}



###age mixing comparisons
#for each parameter combination: mean(inc5/prev at y25) + 90% interval (sorted: value 6-95)

simcompvalues<-function(filename) {
  result_table <- fread(paste0(sourcedir, filename), header=TRUE, stringsAsFactors=FALSE)
  
  #parameter values
  namesplit<-str_split_fixed(filename, "-", 6)
  type<-namesplit[,1]
  population<-as.numeric(namesplit[,2])
  meanagegap<-as.numeric(namesplit[,3])
  sdagegap<-as.numeric(namesplit[,4])
  eyecap<-as.numeric(namesplit[,5])
  transmission<-as.numeric(substring(namesplit[,6], 2, 2))

  #simulation values 
  result35_mean <- mean(result_table$"35")
  result35_low <- sort(result_table$"35")[6]
  result35_high <- sort(result_table$"35")[95]
  return(c(population, meanagegap, sdagegap, eyecap, transmission, result35_mean, result35_low, result35_high))
}


for (type in c("inc1-", "inc5-", "prev-")) {
  filelist <- list.files(path=sourcedir, pattern = type)
  targetfilename <- paste0(graphdir, "simcomp_", substring(type, 1, 4), ".csv")
  simcompdata <- t(sapply(filelist, simcompvalues))
  colnames(simcompdata)<-c("population","meanagegap","sdagegap","eyecap","transmission","Y35.mean","Y35.low","Y35.high")
  write.csv(simcompdata, targetfilename, row.names = FALSE)
  
  transmission.s<-1
    for (eyecap.s in c(300, 1000)) { 
      #comparison meanagegap
      for (sdagegap.s in c(1,2,3,5)){
        plotdata <- subset(as.data.frame(simcompdata), eyecap==eyecap.s & transmission==transmission.s & sdagegap==sdagegap.s, select=population:Y35.high)
        
        if (type=="prev-") {ylabel<-"Prevalence at Y35"} 
        else if (type=="inc1-") {ylabel<-"incidence1 at Y35"}
        else if (type=="inc5-") {ylabel<-"incidence5 at Y35"}
        xlabel<-paste0("Population size","\n\n sd age gap: ", sdagegap.s, "  eyecap: ", eyecap.s, "  transmission: ", transmission.s)
        graphfilename<-paste0("compmeanagegap-", type, sdagegap.s,"-", eyecap.s,"-t", transmission.s, ".eps")
        
        pd <- position_dodge(0.2)
        ggplot(plotdata, aes(x=as.factor(population), y=Y35.mean, group=factor(meanagegap), color=factor(meanagegap))) + 
          geom_errorbar(aes(ymin=Y35.low, ymax=Y35.high), width=.2, position=pd) +
          geom_line(position=pd) +
          geom_point(position=pd) +
          labs(
            x = xlabel,
            y = ylabel)
        ggsave(paste0(graphdir, graphfilename), device="eps")
      }
      #comparison sd
      for (meanagegap.s in c(0, 1, 2, 5, 10)){
        plotdata <- subset(as.data.frame(simcompdata), eyecap==eyecap.s & transmission==transmission.s & meanagegap==meanagegap.s, select=population:Y35.high)
        if (type=="prev-") {ylabel<-"Prevalence at Y35"} 
        else if (type=="inc1-") {ylabel<-"incidence1 at Y35"}
        else if (type=="inc5-") {ylabel<-"incidence5 at Y35"}
        xlabel<-paste0("Population size","\n\nmean age gap: ", meanagegap.s, "  eyecap: ", eyecap.s, "  transmission: ", transmission.s)
        graphfilename<-paste0("compsdagegap-", type, meanagegap.s,"-", eyecap.s,"-t", transmission.s, ".eps")
        
        pd <- position_dodge(0.2)
        ggplot(plotdata, aes(x=as.factor(population), y=Y35.mean, group=factor(sdagegap), color=factor(sdagegap))) + 
          geom_errorbar(aes(ymin=Y35.low, ymax=Y35.high), width=.2, position=pd) +
          geom_line(position=pd) +
          geom_point(position=pd) +
          labs(
            x = xlabel,
            y = ylabel)
        ggsave(paste0(graphdir, graphfilename), device="eps")
      }
    }
  #comparison eyecap
  for (meanagegap.s in c(0, 1, 2, 5, 10)){
    for (sdagegap.s in c(1,2,3,5)){
      plotdata <- subset(as.data.frame(simcompdata), transmission==transmission.s & meanagegap==meanagegap.s & sdagegap==sdagegap.s, select=population:Y35.high)
      if (type=="prev-") {ylabel<-"Prevalence at Y35"} 
      else if (type=="inc1-") {ylabel<-"incidence1 at Y35"}
      else if (type=="inc5-") {ylabel<-"incidence5 at Y35"}
      xlabel<-paste0("Population size","\n\nmean age gap: ", meanagegap.s, "  sd age gap: ", sdagegap.s, "  transmission: ", transmission.s)
      graphfilename<-paste0("compeyecap-", type, meanagegap.s,"-", sdagegap.s,"-t", transmission.s, ".eps")
      
      pd <- position_dodge(0.2)
      ggplot(plotdata, aes(x=as.factor(population), y=Y35.mean, group=factor(eyecap), color=factor(eyecap))) + 
        geom_errorbar(aes(ymin=Y35.low, ymax=Y35.high), width=.2, position=pd) +
        geom_line(position=pd) +
        geom_point(position=pd) +
        labs(
          x = xlabel,
          y = ylabel)
      ggsave(paste0(graphdir, graphfilename), device="eps")
    }
  }
}







######extra analysis
# 1. qqplot -> normality
# 2. variance~popsize (compmeanagegap/compsdagegap)
# 3. Ymax~popsize (compmeanagegap/compsdagegap)  #analysis speed epidemic - population size
###### 

###qqplot
for (type in c("inc1-", "inc5-", "prev-")) {
  filelist <- list.files(path=sourcedir, pattern = type)
  for (file in filelist) {
    targetfilename <- paste0(graphdir, "extra/qqnorm-", str_split(file, "\\.")[[1]][1], ".svg")
    result_table <- fread(paste0(sourcedir, file), header=TRUE, stringsAsFactors=FALSE)
    svg(targetfilename)
    qqnorm(result_table$"35")
    dev.off()
  }
}

###variance graphs
for (type in c("inc1", "inc5", "prev")) {
  #inlezen simcompdataextra
  simcompfilename <- paste0(graphdir, "simcomp_ext_", type, ".csv")
  simcomp_extra <- fread(paste0(simcompfilename), header=TRUE, stringsAsFactors=FALSE)
  transmission.s<-1
  
  for (eyecap.s in c(300, 1000)) { 
    #comparison meanagegap
    for (sdagegap.s in c(1,2,3,5)){
      plotdata <- subset(as.data.frame(simcomp_extra), eyecap==eyecap.s & transmission==transmission.s & sdagegap==sdagegap.s, select=population:max.Y.var)
      if (type=="prev") {ylabel<-"Variance prevalence at Y35"} 
      else if (type=="inc1") {ylabel<-"Variance incidence1 at Y35"}
      else if (type=="inc5") {ylabel<-"Variance incidence5 at Y35"}
      xlabel<-paste0("Population size","\n\n sd age gap: ", sdagegap.s, "  eyecap: ", eyecap.s, "  transmission: ", transmission.s)
      graphfilename<-paste0("compvar-meanagegap-", type,"-" , sdagegap.s,"-", eyecap.s,"-t", transmission.s, ".eps")
    
      pd <- position_dodge(0.2)
      ggplot(plotdata, aes(x=as.factor(population), y=Y35.var, group=factor(meanagegap), color=factor(meanagegap))) + 
        geom_line(position=pd) +
        geom_point(position=pd) +
        labs(
          x = xlabel,
          y = ylabel)
      ggsave(paste0(graphdir, "extra/", graphfilename), device="eps")
    }
  #comparison sd
    for (meanagegap.s in c(0, 1, 2, 5, 10)){
      plotdata <- subset(as.data.frame(simcomp_extra), eyecap==eyecap.s & transmission==transmission.s & meanagegap==meanagegap.s, select=population:max.Y.var)
      if (type=="prev") {ylabel<-"Variance prevalence at Y35"} 
      else if (type=="inc1") {ylabel<-"Variance incidence1 at Y35"}
      else if (type=="inc5") {ylabel<-"Variance incidence5 at Y35"}
      xlabel<-paste0("Population size","\n\nmean age gap: ", meanagegap.s, "  eyecap: ", eyecap.s, "  transmission: ", transmission.s)
      graphfilename<-paste0("compvar-sdagegap-", type, "-", meanagegap.s,"-", eyecap.s,"-t", transmission.s, ".eps")
    
      pd <- position_dodge(0.2)
      ggplot(plotdata, aes(x=as.factor(population), y=Y35.var, group=factor(sdagegap), color=factor(sdagegap))) + 
        geom_line(position=pd) +
        geom_point(position=pd) +
        labs(
          x = xlabel,
          y = ylabel)
      ggsave(paste0(graphdir, "extra/", graphfilename), device="eps")
    }
  }
}

### Ymax~popsize (compmeanagegap/compsdagegap)  #analysis speed epidemic - population size

for (type in c("inc1", "inc5", "prev")) {
  simcompfilename <- paste0(graphdir, "simcomp_ext_", type, ".csv")
  simcomp_extra <- fread(paste0(simcompfilename), header=TRUE, stringsAsFactors=FALSE)
  transmission.s<-1
  
  for (eyecap.s in c(300, 1000)) { 
    #comparison meanagegap
    for (sdagegap.s in c(1,2,3,5)){
      plotdata <- subset(as.data.frame(simcomp_extra), eyecap==eyecap.s & transmission==transmission.s & sdagegap==sdagegap.s, select=population:max.Y.var)
      if (type=="prev") {ylabel<-"Y max prevalence"} 
      else if (type=="inc1") {ylabel<-"Y max incidence1"}
      else if (type=="inc5") {ylabel<-"Y max incidence5"}
      xlabel<-paste0("Population size","\n\n sd age gap: ", sdagegap.s, "  eyecap: ", eyecap.s, "  transmission: ", transmission.s)
      graphfilename<-paste0("compYmax-meanagegap-", type,"-" , sdagegap.s,"-", eyecap.s,"-t", transmission.s, ".eps")
      
      pd <- position_dodge(0.2)
      ggplot(plotdata, aes(x=as.factor(population), y=max.Y, group=factor(meanagegap), color=factor(meanagegap))) + 
        geom_line(position=pd) +
        geom_point(position=pd) +
        labs(
          x = xlabel,
          y = ylabel)
      ggsave(paste0(graphdir, "extra/", graphfilename), device="eps")
    }
    #comparison sd
    for (meanagegap.s in c(0, 1, 2, 5, 10)){
      plotdata <- subset(as.data.frame(simcomp_extra), eyecap==eyecap.s & transmission==transmission.s & meanagegap==meanagegap.s, select=population:max.Y.var)
      if (type=="prev") {ylabel<-"Y max prevalence"} 
      else if (type=="inc1") {ylabel<-"Y max incidence1"}
      else if (type=="inc5") {ylabel<-"Y max incidence5"}
      xlabel<-paste0("Population size","\n\nmean age gap: ", meanagegap.s, "  eyecap: ", eyecap.s, "  transmission: ", transmission.s)
      graphfilename<-paste0("compYmax-sdagegap-", type, "-", meanagegap.s,"-", eyecap.s,"-t", transmission.s, ".eps")
      
      pd <- position_dodge(0.2)
      ggplot(plotdata, aes(x=as.factor(population), y=max.Y, group=factor(sdagegap), color=factor(sdagegap))) + 
        geom_line(position=pd) +
        geom_point(position=pd) +
        labs(
          x = xlabel,
          y = ylabel)
      ggsave(paste0(graphdir, "extra/", graphfilename), device="eps")
    }
  }
}



