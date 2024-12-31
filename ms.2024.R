#==== Dependencies
library("tidyverse")

reddtimingoffset <- 7 # days to shift a carcass timing back to spawning/redd deposition day
allyears <- 2004:2023

#==== FUNCTIONS ====

zAdjust <- function(obs=rawcarcass2021,cc=ccmatrix){
  x1 <- sum(obs[,2])
  x2 <- sum(obs[,3])
  x3 <- sum(obs[,4])
  x4 <- sum(obs[,5])
  A1 <- x1/(1-cc[2,2]- cc[3,2] - cc[4,2])
  A2 <- (x2 - A1*cc[2,2]) / (1 - cc[3,3] - cc[4,3])
  A3 <- (x3 - A2*cc[3,3]-A1*cc[3,2]) / (1 - cc[4,4])
  A4 <- x4 - A3*cc[4,4] - A2*cc[4,3] - A1*cc[4,2]
  return(c(A1,A2,A3,A4))
}

zExpand <- function(raw=c(0,10,10,0),scal = 1,method="round"){
  if(scal <= 0){out <- rep(0,length(raw));return(out)}
  out <- raw
  tot <- sum(raw)*scal
  # step through each day, add up extra carcasses (proportion on date/ total) until summed to one, then remove from total
  tally <- 0
  for(j in 1:length(raw)){
    tally <- tally + raw[j]*scal
    # ceiling method
    if(method == "round"){
      if(tally > 0){
        out[j] <- round(tally)
      }
      tally <- tally-round(tally);
    } else {
      if(tally > -1){
        out[j] <- ceiling(tally) ; tally <- tally-ceiling(tally);
      }
    }
  }
  return(out)
}

#==== Graphing functions ====
zridges <- function(usescale=40,maxyr = 2024,minyr=2004,type="",use1=aerial2,t1="",
                    hist=NULL,colo =NULL,xlimm=NULL,localaxis=NULL){
  par(mar=c(4.5,3,3,1),mgp=c(2.5,1,0))
                 cexx <- par()$cex*1.0
  if(type=="png")cexx <- par()$cex*1.2
  if(type=="pdf")cexx <- par()$cex*1.2
  cexaxis <- cexx*0.9
  yrz <- minyr:maxyr
  totalz <- rep(0,length(yrz)) ; for(y in yrz){totalz[match(y,yrz)] <- sum(use1$N[use1$year==y])}
  # 108 is earliest survey DOY in carcass or aerial
  
  if(is.null(xlimm))xlimm <- c(89,250)
  if(is.null(colo))colo <- "darkgreen"
  plot(0,0,xlim=xlimm,ylim=c(minyr,maxyr+1),axes=FALSE,xlab="",ylab="",cex=cexx,xaxs="i", yaxs="i");
  abline(v=c(90,121,152,182,213,244),lwd=1,lty=2,col="grey90")
  axis(1,at=c(90,121,152,182,213,244),labels=c("Apr 1","May 1","Jun 1","Jul 1","Aug 1","Sep 1"),cex.axis=cexaxis*1)
  j <- 0 + minyr-2004
  for(y in yrz){
    if(usescale=="local"){ drawscale <- ceiling(max(10,max(use1$N[use1$year==y])*1.4)/10)*10 } 
    if(usescale=="limited"){ drawscale <- ifelse(max(use1$N[use1$year==y]) > 40,300,50) }
    if(is.numeric(usescale)){drawscale=usescale}
    # print(drawscale) ; print(paste(y,"drawscale",drawscale,max(use1$N[use1$year==y])))
    zdraw(use2=use1[use1$year==y,],scaleit=drawscale,y=y,hist=hist,colo=colo,xlimm=xlimm,localaxis=TRUE,cexaxis=cexaxis)
    abline(h=y)
  }
  mtext(t1,3,cex=1.5,line=3)
  text(rep(95,length(yrz)),yrz+0.5,yrz,cex=cexx,adj=0) ; # paste0(yrz,"\nN =",totalz)
}


zdraw <- function(use2,scaleit=30,y=year,hist=NULL,colo=NULL,xlimm=NULL,localaxis=NULL,cexaxis=1){
  if(is.null(xlimm)) I <- 118:243 else I <- min(xlimm):max(xlimm)
  j <- j + 1
  yy <- use2$N[I]
  yy <- y +yy/scaleit
  a <- min(50+j*10,255); 
  b <- 100 
  cc <-  max(250-j*10,10)
  if(is.null(colo))colo <- rgb(a,b,cc,150,NULL,255)
  if(is.null(hist)){
    polygon(c(I,rev(I)),c(rep(y,length(I)),rev(yy)),col=colo,border=colo )
    lines(I,yy,col=coll)
  } else{
    for(i in 1:length(I)){  
      arrows(I[i],y,I[i],yy[i],length=0,col=colo,lwd=2,lend=2) 
    }
  }
  if(!is.null(localaxis)){# draw the y axis
    if(scaleit== 50) axlabs <- c(0,10,20,30,40) else axlabs=c(0,50,100,150,200,250)
     ats <- y + axlabs/scaleit
    # cat(y," ");cat(max(use2$N[I])," "); cat(max(yy)," ") ; cat(max(axlabs)," ") ; cat(max(ats),par()$usr[4],"\n")
    axis(2,at=ats,labels=axlabs,las=2,cex.axis=cexaxis) # needs to be only on the year range specific areas
  }
}

#==== Read files ====
ccmatrix <- read.table("carcass.conversion.matrix.csv",sep=",",row.names=NULL,as.is=TRUE,header=TRUE)
x <- read.csv("carcasses.2004.2023.csv")  # reboot.carcasses.2024.csv saved from Killam workbook WR SPAWNER spatial-temporal distrbution 9-14-21.xlsx
x <- x %>% mutate(year=year(as.Date(mdy(Date))))
print(head(x))


#==== Process the raw carcass files to get adjusted numbers ====
annualtotalscarcass <- NULL
for(year in allyears){
  annualtotalscarcass <- rbind.data.frame(annualtotalscarcass,cbind.data.frame(year,raw=sum(x$SumOfCount[x$year==year])))
}
print("Annualtotals from carcasses");print(annualtotalscarcass)

carcassraw <- list()
carcassrawsections <- NULL
carcassAdjustsections <- NULL
for(year in allyears){
  w <- x[x$year==year,]
  w1 <- w %>% mutate(Date2 = as.Date(mdy(Date))) %>% select(-Date) %>% rename(Date=Date2)
  xx <- w1 %>% group_by(year,Date,Section) %>% summarise(SumOfCount= sum(SumOfCount))%>% pivot_wider(names_from = Section,
        values_from=SumOfCount,values_fill=0) %>% 
        mutate(Day=yday(Date)) %>% ungroup() %>% select(-c(year,Date))
  # Order the data for consistency 
  if(match("4",names(xx),nomatch=0) == 0){xx <- bind_cols(xx,`4`=0)}
  xx <- xx %>% relocate(`1`,.after = `Day`)
  xx <- xx %>% relocate(`2`,.after = `1`)
  xx <- xx %>% relocate(`3`,.after = `2`)
  xx <- xx %>% relocate(`4`,.after = `3`)
  carcassraw[[as.character(year)]] <- xx
  x2 <- cbind.data.frame("year"=year,"one"=sum(xx$`1`),"two"=sum(xx$`2`),"three"=sum(xx$'3'),"four"=ifelse('4' %in% names(xx),sum(xx$`4`),0))
  carcassrawsections <- rbind.data.frame(carcassrawsections,x2) 
  j <- round(zAdjust(obs=x2,cc=ccmatrix))
  carcassAdjustsections <- rbind.data.frame(carcassAdjustsections,cbind.data.frame("year"=year,"one"=j[1],"two"=j[2],"three"=j[3],"four"=j[4]))
}
carcassrawsections$Total= rowSums(carcassrawsections[,-1])
carcassAdjustsections$Total= rowSums(carcassAdjustsections[,-1])
print(names(carcassraw))
print(carcassrawsections)
print(carcassAdjustsections)

#==== Generate carcassexpand list ====

carcassexpand <- list() 
carcassexpand.type478 <- list() 
carcassraw.fakeexpand <- list()
expand.raw.diffs <- NULL
outDirectory <- "generated/"
for(year in allyears){ 
   # use<- carcassraw[[as.character(year)]] %>% rename(RKM483='1') %>% rename(RKM478='2' ) %>% rename(RKM474='3') %>% rename(RKM454='4')
  use<- carcassraw[[as.character(year)]] %>% rename(RKM483='1') %>% rename(RKM479='2' ) %>% rename(RKM474='3') %>% rename(RKM454='4')
  
  holdDat <- cbind.data.frame(Year=year,X1=NA,X2=NA,X3=NA,X4=NA,Totals=NA)
    new <- zAdjust(use,ccmatrix)
    holdDat[1,2:5] <-new
    holdDat[1,6] <-sum(new)
    newuse <- use
    newuseceiling <- use
    for(colu in 2:ncol(use)){
      surveyraw <- sum(use[,colu])
      expp <- holdDat[holdDat$Year == year,colu]
      scalar <- 0
      if(surveyraw > 0) scalar <- expp/surveyraw
      # print(paste(sect,colu,expp,surveyraw,scalar))
      if(scalar > 0){
        newuse[[colu]] <- zExpand(use[[colu]] ,scalar,method="round") 
        newuseceiling[[colu]] <- zExpand(use[[colu]] ,scalar,method="ceiling") 
      } else {
        newuse[[colu]]        <- rep(0,length=length(newuse[[colu]]))
        newuseceiling[[colu]] <- rep(0,length=length(newuseceiling[[colu]])) 
      }
    }
    # find blank rows
    newuse2 <- newuse[!apply(newuse[,-1] == 0, 1, all),]  # works
    newuse2[[1]] <- newuse2[[1]]-reddtimingoffset
    carcassraw.fakeexpand[[as.character(year)]] <- use
    # carcassexpand[[as.character(year)]] <- newuse2
    # carcassexpand.type478[[as.character(year)]] <- newuse2
    write_csv(use,paste0(outDirectory,"carcassraw.fakeexpand.",year,".csv"))
    #write_csv(newuse2,paste0(outDirectory,"carcassExpand.",year,".csv"))
    # write_csv(newuse2,paste0(outDirectory,"carcassExpand.type478.",year,".csv"))
    
    if(year==min(allyears))print(paste("year","raw","expand"))
    print(paste(year,sum(use[,-1]),sum(newuse2[,-1]),sum(use[,-1]) - sum(newuse2[,-1])))
    expand.raw.diffs <- rbind.data.frame(expand.raw.diffs,cbind.data.frame(year,raw=sum(use[,-1]),expand=sum(newuse2[,-1]),diff=sum(use[,-1]) - sum(newuse2[,-1])))
}
print(summary(expand.raw.diffs))

#==== WEB call to CBR database for aerial redds ====
# if(0) will prevent web calls and load from static file locally generated
if(0){
  usereddsaerial <- list()
  i <- 0
  for(year in allyears){
    print(year)
    raw <- 5
    string <- "http://cbr.washington.edu/sac-bin/fishmodel/getandplottemp.pl?dirUseId=test47&temponly=yes&tempsource=dbtemp&hatchmechanism=atu"  
    string2 <- paste(string,"&tempyear=",year,"&raw=",raw,"&redds=dbredds&reddyear=",year,sep="")
    usereddsaerial[[as.character(year)]] <- read.csv(string2)
  }
  save(usereddsaerial,file="usereddsaerial.Rdata")
} else {
  load("usereddsaerial.Rdata",verbose = TRUE)
}


#=== Tally aerial by in/out of carcass area (RKM444 or above) ====
# RKM names in aerial redd  data are locations used WITHIN the reaches. Compare to table 1 in manuscript
allabove <- allbelow <- 0
for(year in allyears){
  junk <- usereddsaerial[[as.character(year)]]
  junk <- junk[,-1]
  above<- sum(junk[,match(c("RKM483","RKM479","RKM470","RKM450"),names(junk),nomatch=0)])
  below<- sum(junk[,match(c("RKM440","RKM430","RKM415") ,names(junk),nomatch=0)])
  print(paste(year,above,below,round(below/(above+below),3)))
  allabove <- allabove + above
  allbelow<- allbelow + below
}
print(paste(allabove,allbelow,allbelow/(allabove+allbelow)))

#=== Process aerial ====
allspawn <- gaps <- annualtotalsaerial <-  aerialsections <- NULL
for(y in allyears){ 
  rr <- usereddsaerial[[as.character(y)]]
  rrbyday <- as_tibble(rr) %>%  group_by(Day) %>% transmute(N = sum(across(starts_with("RKM"))))
  newgaps <- rrbyday$Day[2:length(rrbyday$Day)] - rrbyday$Day[1:(length(rrbyday$Day)-1)]
  allspawn[[as.character(y)]] <- rrbyday
  gaps <- c(gaps,newgaps)
  annualtotalsaerial <- rbind.data.frame(annualtotalsaerial,
                              cbind.data.frame(y,raw=sum(rrbyday$N)))
  aerialsections <- rbind.data.frame(aerialsections,cbind.data.frame( year=y,"R1"=sum(rr$RKM483),"R2"=sum(rr$RKM479),"R3"=sum(rr$RKM470),"R4"=sum(rr$RKM450)))
}
print("Annual Totals Aerial");print(annualtotalsaerial)
print(aerialsections)

#==== Aerial process to spread out across an entire year, than back fill =====
aerial1 <- aerial2 <- NULL
for(year in allyears){
  aerial1 <- rbind.data.frame(aerial1,cbind.data.frame(year=year,allspawn[[as.character(year)]]))
  aerial2 <- rbind.data.frame(aerial2,cbind.data.frame(year=rep(year,365),doy=1:365,N=0))
}
for(i in 1:nrow(aerial1)){
  I <- aerial2$year == aerial1$year[i] & aerial2$doy==aerial1$Day[i]; 
  aerial2$N[I] <- aerial1$N[i];
}



for(y in allyears){
  aerial2$cumN[aerial2$year == y] <- cumsum(aerial2$N[aerial2$year == y]) 
}

LCL <- UCL <- LCL95 <- UCL95 <- MedianAer <- MedianAerDate <- rep(NA,length(allyears))
for(year in allyears){
  i <- year-min(allyears) +1
  mx <- max(aerial2$cumN[aerial2$year==year])
  LCL[i] <- min(aerial2$doy[aerial2$year==year & aerial2$cumN >= 0.1*mx])
  UCL[i] <- max(aerial2$doy[aerial2$year==year & aerial2$cumN <= 0.9*mx])
  # median is a little tricky. want the mean day if the median balue is on multiple days
  # Find maximum day below the median and minimum day above the median and average the gap
  low <- max(aerial2$doy[aerial2$year==year & aerial2$cumN <= 0.5*mx])
  hi <- min(aerial2$doy[aerial2$year==year & aerial2$cumN >= 0.5*mx])
  MedianAer[i] <- mean(low,hi)
  junk <- as.Date(paste0(year-1,"-12-31")) + days(MedianAer[i])
  MedianAerDate[i] <- paste(month(junk,label=TRUE,abbr=TRUE),day(junk)) 
  LCL95[i] <- min(aerial2$doy[aerial2$year==year & aerial2$cumN >= 0.025*mx])
  UCL95[i] <- max(aerial2$doy[aerial2$year==year & aerial2$cumN <= 0.975*mx])
}

#==== process carcass ====
# Don't need gaps, because we can smooth 3 days twice
# Also, normalize the spawning distributions to visualize inter-annual normalcy
allspawnCar <- allnormalized <- allmeans <- allstds <- alljunk <- NULL
for(y in allyears){ 
  rr <- carcassexpand[[as.character(y)]]
  rrbyday <- as_tibble(rr) %>%  group_by(Day) %>% transmute(N = sum(across(starts_with("RKM"))))
  allspawnCar[[as.character(y)]] <- rrbyday
  junk <- rep(rrbyday$Day,rrbyday$N)
  thismean <- mean(junk)
  thisstd <- sqrt(var(junk))
  junk2 <- (junk-thismean)/thisstd
  allnormalized <- rbind.data.frame(allnormalized,cbind.data.frame("year"=y,"Nday"=rrbyday$Day - thismean,"N"=rrbyday$N))
  allmeans <- c(allmeans,thismean)
  allstds <- c(allstds,thisstd)
  alljunk <- c(alljunk,junk)
}
allGstd <- sqrt(var(allnormalized$Nday))
allGmean <- mean(allmeans)
hist(rep(allnormalized$Nday,allnormalized$N),breaks=50,freq=FALSE,main="")
x<- allGmean+seq(-80,80,by=5)
lines(x-allGmean,dnorm(x,mean=allGmean,sd=allGstd))
curve(dnorm(x,mean=allGmean,sd=allGstd),add=TRUE)

carcass1 <- carcass2 <- NULL
for(year in allyears){
  carcass1 <- rbind.data.frame(carcass1,cbind.data.frame(year=year,allspawnCar[[as.character(year)]]))
  carcass2 <- rbind.data.frame(carcass2,cbind.data.frame(year=rep(year,365),doy=1:365,N=0))
}
for(i in 1:nrow(carcass1)){
  I <- carcass2$year == carcass1$year[i] & carcass2$doy==carcass1$Day[i]; 
  carcass2$N[I] <- carcass1$N[i];
}
carcass2$cumN <- NA
# carcass3.1 formerly used here. change back to carcass2
MedianCar <- MedianCarDate <- LCLc <- UCLc <- LCLc95 <- UCLc95 <- inferredCarRedds <- rep(NA,length(allyears))
for(year in allyears){
  carcass2$cumN[carcass2$year==year] <- cumsum(carcass2$N[carcass2$year==year])
  i <- year-min(allyears)+1
  mx <- max(carcass2$cumN[carcass2$year==year])
  inferredCarRedds[i] <- mx
  LCLc[i] <- min(carcass2$doy[carcass2$year==year & carcass2$cumN >= 0.1*mx])
  UCLc[i] <- max(carcass2$doy[carcass2$year==year & carcass2$cumN <= 0.9*mx])
  # median is a little tricky. want the mean day if the median value is on multiple days
  # Find maximum day below the median and minimum day above the median and average the gap
  low <- max(carcass2$doy[carcass2$year==year & carcass2$cumN <= 0.5*mx])
  hi <- min(carcass2$doy[carcass2$year==year & carcass2$cumN >= 0.5*mx])
  MedianCar[i] <- mean(low,hi)
  junk <- as_date(paste0(year-1,"-12-31")) + days(MedianCar[i])
  MedianCarDate[i] <- paste(month(junk,label=TRUE,abbr=TRUE),day(junk)) 
  LCLc95[i] <- min(carcass2$doy[carcass2$year==year & carcass2$cumN >= 0.025*mx])
  UCLc95[i] <- max(carcass2$doy[carcass2$year==year & carcass2$cumN <= 0.975*mx])
}

#==== Additional diagnostics ====
results <-cbind.data.frame("year"=allyears,"MedianAerial"=MedianAer,"MedianAerDate"=MedianAerDate,"MedianCarcass"=MedianCar,"MedianCarcassDate"=MedianCarDate,
                           "DiffsMedian"=MedianAer-MedianCar, "Aerial95"=UCL95-LCL95,"Carcass95"=UCLc95-LCLc95,
                           "Diffs95"= (UCL95-LCL95) - (UCLc95-LCLc95),aerialcount=annualtotalsaerial[,2],carcasscount=annualtotalscarcass$raw,
                       adjust=inferredCarRedds,delta=inferredCarRedds-annualtotalscarcass$raw,
                       car.aer.ratio=annualtotalscarcass$raw/annualtotalsaerial$raw)

print(results)
summary(results)
summary(abs(results$DiffsMedian))
summary(abs(results$Diffs95))

# SD inferredCarRedds
print(sqrt(var(results$adjust)))
print(sqrt(var(results$carcasscount)))
print(sqrt(var(results$aerialcount)))
print(sqrt(var(results$MedianAerial)))
print(sqrt(var(results$MedianCarcass)))
print(sqrt(var(results$Aerial95)))
print(sqrt(var(results$Carcass95)))

write.csv(results,file="results.csv",row.names=F)


#==== Draw fig with  counts in timeseries ====
for(fig in c("pdf")){   # """png" ){   #  
  if(fig == "pdf") pdf(file="counttimeseries.pdf",width=7,height=7)
  if(fig == "png") png(file="counttimeseries.png",width=700,height=700)
  colos <- c("darkgreen","orange","maroon","darkblue")
  par(mfrow=c(2,1))
  par(mar=c(3,4,1,1))
  x <- aerialsections
  x[x<0] = 0
  spacespace <- 0.2
  barplot(t(x[,c(5:2)]),border=NA,col=colos,ylim=c(0,2000),space=spacespace,ylab="Counts")
  legend("topright",bty="n",fill=rev(colos),legend=c("Reach 1","Reach 2","Reach 3","Reach 4"))
  axis(1,labels=allyears,at=(allyears-2003.5)*(1+ spacespace))
  abline(h=0)
  text(4,1800,"Aerial Survey Distributions by Reach",adj=0)
  x <- carcassAdjustsections
  x[x<0] = 0
  barplot(t(x[,c(5:2)]),border=NA,col=colos,ylim=c(0,3000),space=spacespace,ylab="Counts")
  legend("topright",bty="n",fill=rev(colos),legend=c("Reach 1","Reach 2","Reach 3","Reach 4"))
  axis(1,labels=allyears,at=(allyears-2003.5)*(1  + spacespace))
  abline(h=0)
  text(4,2800,"Adjusted Carcass Survey Distributions by Reach",adj=0)
  if(fig=="pdf" | fig == "png") dev.off()
}

#==== Draw figure for  ms. Four Panel Figure====
for(fig in c("" )){  #   c("","png","pdf")){  #   
  if(fig == "pdf")  pdf(file="bothdistrib.pdf",width=14,height=18)
  if(fig == "png")  png(file="bothdistrib.png",width=1000,height=1200)
  par(mfrow=c(1,4))
  par(mar=c(4.5,3,2,1),mgp=c(2.5,1,0))
  par(cex=1)
  cexaxis <- 1 *par()$cex ; cexoutside <- par()$cex ;  cexpanel <-  1.2*par()$cex
  if(fig == "pdf"){
    cexaxis <- 1 *par()$cex ; cexoutside <- 1.5*par()$cex ;  cexpanel <-  1.5*par()$cex
  }

  # For shading
  usecol <- rgb(10,10,10,3,NULL,25)
  # usescale=65
  usescale="limited"
  coll1 <- coll2 <- NULL
  coll1 <- rgb(5,10,5,20,NULL,25)
  coll2 <- rgb(10,5,10,20,NULL,25)
  usehist="hist" ; # else usehist=NULL

  want <- 2014:2023 ; I <- aerial2$year >= min(want) & aerial2$year <= max(want)
  zridges(usescale=usescale,maxyr=max(want),minyr=min(want),type=fig,use1 =aerial2[I,] ,t1 = "", colo=coll1,hist=usehist)
  #polygon(c(LCL[match(want,allyears)],rev(UCL[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  polygon(c(LCL95[match(want,allyears)],rev(UCL95[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  lines(MedianAer[match(want,allyears)],want)
      
  mtext("(a) Aerial Surveys",line=0.5,adj=1,cex=cexpanel)

  want <- 2004:2013 ; I <- aerial2$year >= min(want) & aerial2$year <= max(want)
  zridges(usescale=usescale,maxyr=max(want),minyr=min(want),type=fig,use1 =aerial2[I,] ,t1 = "",colo=coll1,hist=usehist)
  #polygon(c(LCL[match(want,allyears)],rev(UCL[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  polygon(c(LCL95[match(want,allyears)],rev(UCL95[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  lines(MedianAer[match(want,allyears)],want)
  
  want <- 2014:2023 ; I <- carcass2$year >= min(want) & carcass2$year <= max(want)
  zridges(usescale=usescale,maxyr=max(want),minyr=min(want),type=fig,use1 = carcass2[I,] ,t1 ="", colo=coll2,hist=usehist)
  #polygon(c(LCLc[match(want,allyears)],rev(UCLc[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  polygon(c(LCLc95[match(want,allyears)],rev(UCLc95[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  lines(MedianCar[match(want,allyears)],want)
  
  mtext("(b) Carcass Surveys",line=0.5,adj=1,cex=cexpanel)
  
  want <- 2004:2013 ; I <- carcass2$year >= min(want) & carcass2$year <= max(want)
  zridges(usescale=usescale,maxyr=max(want),minyr=min(want),type=fig,use1 = carcass2[I,] ,t1 ="", colo=coll2,hist=usehist)
  #polygon(c(LCLc[match(want,allyears)],rev(UCLc[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  polygon(c(LCLc95[match(want,allyears)],rev(UCLc95[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  lines(MedianCar[match(want,allyears)],want)

  mtext("Date",side=1,outer=TRUE,line=-2,cex=cexoutside)
  if(fig != "") dev.off()
}


