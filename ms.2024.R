#==== Dependencies
library("tidyverse")

reddtimingoffset <- 7 # days to shift a carcass timing back to spawning/redd deposition day
allyears <- 2004:2023

#==== FUNCTIONS ====
zrunavg <- function(x, stretch = 1){
  if(stretch < 1 || (floor(stretch) != stretch)) {
    return("Stretch must be integer: 1 or greater")
  }
  y <- x
  for(i in 1:stretch){y[i] <- mean(x[1:i], na.rm = TRUE)}
  for(i in stretch:length(x)) {
    y[i] <- mean(x[(i - stretch + 1):i], na.rm = TRUE)
  }
  y
}
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
zExpand <- function(raw=c(0,10,10,0),scal = 1){
  if(scal <= 0){out <- rep(0,length(raw));return(out)}
  out <- raw
  tot <- sum(raw)*scal
  # step through each day, add up extra carcasses (proportion on date/ total) until summed to one, then remove from total
  tally <- 0
  tallytally <- 0
  for(j in 1:length(raw)){
    tally <- tally + raw[j]*scal
    if(tally > -1){
      out[j] <- ceiling(tally) ; tally <- tally-ceiling(tally);
    } else {print(paste("Tally",tally))}
  }
  return(out)
}

#==== Graphing functions ====
zridges <- function(scale=40,maxyr = 2024,minyr=2004,type="",use1=aerial3,t1="",colo =NULL,xlimm=NULL){
  par(mar=c(4,2,1,0),mgp=c(1.8,0.4,0))
  cexx <- 1
  if(type=="png")cexx <- par()$cex*1.8
  if(type=="pdf")cexx <- par()$cex*1
  cexaxis <- cexx
  cexaxis <- 1.3 *cexx
  yrz <- minyr:maxyr
  totalz <- rep(0,length(yrz)) ; for(y in yrz){totalz[match(y,yrz)] <- sum(use1$N[use1$year==y])}
   # 108 is earliest survey DOY in carcass or aerial
  
  if(is.null(xlimm))xlimm <- c(108,249)
  if(is.null(colo))colo <- "darkgreen"
  plot(0,0,xlim=xlimm,ylim=c(minyr,maxyr+0.8),axes=FALSE,xlab="",ylab="",cex=cexx,xaxs="i", yaxs="i");
  
  text(rep(125,length(yrz)),yrz+0.5,yrz,cex=cexaxis,adj=0) ; # paste0(yrz,"\nN =",totalz)
  abline(v=c(121,152,182,213,244),lwd=1,lty=2,col="grey90")
  axis(1,at=c(121,152,182,213,244),labels=c("May 1 ","Jun 1 ","Jul 1 ","Aug 1 ","Sep 1  "),cex.axis=cexaxis)
   j <- 0 + minyr-2004
  for(y in yrz){
    zdraw(use2=use1[use1$year==y,],scaleit=scale,y=y,colo=colo,xlimm=xlimm)
  }
  text(260,maxyr+0.7,t1,cex=1.5,adj=0)
}

zdraw <- function(use2,scaleit=30,y=year,colo=NULL,hist=NULL,xlimm=NULL){
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
      # segments(I[i],y,I[i],yy[i],col=colo,lwd=4) 
      arrows(I[i],y,I[i],yy[i],length=0,col=colo,lwd=2,lend=2) 
      
    }
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
  # moot for allyears 2004 - 2023 bu retained for available pre2004 processing
  w$Section[w$Section==2.1] <- 2
  w$Section[w$Section==1.2] <- 2
  w$Section[w$Section==3.2] <- 3
  w$Section[w$Section==2.3] <- 3
  w$Section[w$Section==4.3] <- 4
  w$Section[w$Section==3.4] <- 4
  w1 <- w %>% mutate(Date2 = as.Date(mdy(Date))) %>% select(-Date) %>% rename(Date=Date2)
  xx <- w1 %>% group_by(year,Date,Section) %>% summarise(SumOfCount= sum(SumOfCount))%>% pivot_wider(names_from = Section,
                                                                                                     values_from=SumOfCount,values_fill=0) %>% mutate(Day=yday(Date)) %>% ungroup() %>% select(-c(year,Date))
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
print(carcassrawsections)[-c(1:4),]
print(carcassAdjustsections)[-c(1:4),]

# Use the carcassraw list to generate a carcassexpand list

carcassexpand <- list()
expand.raw.diffs <- NULL
outDirectory <- "generated/"
for(year in allyears){ 
    use<- carcassraw[[as.character(year)]] %>% rename(RKM483='1') %>% rename(RKM479='2' ) %>% rename(RKM474='3') %>% rename(RKM454='4')
    holdDat <- cbind.data.frame(Year=year,X1=NA,X2=NA,X3=NA,X4=NA,Totals=NA)
    new <- zAdjust(use,ccmatrix)
    holdDat[1,2:5] <-new
    holdDat[1,6] <-sum(new)
    newuse <- use
    newuse[[1]] <- newuse[[1]]-reddtimingoffset
    for(colu in 2:(dim(use)[2])){
      sect <- colu - 1
      surveyraw <- sum(use[,colu])
      expp <- holdDat[holdDat$Year == year,colu]
      scalar <- 0
      if(surveyraw > 0) scalar <- expp/surveyraw
      # print(paste(sect,colu,expp,surveyraw,scalar))
      if(scalar > 0){
          newuse[[colu]] <- zExpand(use[[colu]] ,scalar) # OK because x is a tibble, right?
      } else {
        newuse[[colu]] <- rep(0,length=length(newuse[[colu]]))
      }
    }
    # find blank rows
    newuse2 <- newuse[!apply(newuse[,-1] == 0, 1, all),]
    carcassexpand[[as.character(year)]] <- newuse2
    write_csv(newuse2,paste0(outDirectory,"carcassExpand.",year,".csv"))
    if(year==min(allyears))print(paste("year","raw","expand"))
    print(paste(year,sum(use[,-1]),sum(newuse2[,-1]),sum(use[,-1]) - sum(newuse2[,-1])))
    expand.raw.diffs <- c(expand.raw.diffs,sum(use[,-1]) - sum(newuse2[,-1]))
}
print(summary(expand.raw.diffs))

# WEB call to CBR database for aerial redds
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


#=== comb through usereddsaerial list for a tally by in or out of carcass area (RKM444 or above) ====
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

#=== process aerial ====
allspawn <- gaps <- annualtotalsaerial <-  NULL
for(y in allyears){ 
  rr <- usereddsaerial[[as.character(y)]]
  rrbyday <- as_tibble(rr) %>%  group_by(Day) %>% transmute(N = sum(across(starts_with("RKM"))))
  newgaps <- rrbyday$Day[2:length(rrbyday$Day)] - rrbyday$Day[1:(length(rrbyday$Day)-1)]
  allspawn[[as.character(y)]] <- rrbyday
  gaps <- c(gaps,newgaps)
  annualtotalsaerial <- rbind.data.frame(annualtotalsaerial,
                              cbind.data.frame(y,raw=sum(rrbyday$N)))
  
}
print("Annual Totals Aerial");print(annualtotalsaerial)
#==== Take survey observations, spread out across an entire year, than back fill =====
aerial1 <- aerial2 <- NULL
for(year in allyears){
  aerial1 <- rbind.data.frame(aerial1,cbind.data.frame(year=year,allspawn[[as.character(year)]]))
  aerial2 <- rbind.data.frame(aerial2,cbind.data.frame(year=rep(year,365),doy=1:365,N=0))
}
for(i in 1:nrow(aerial1)){
  I <- aerial2$year == aerial1$year[i] & aerial2$doy==aerial1$Day[i]; 
  aerial2$N[I] <- aerial1$N[i];
}

aerial3 <- aerial2
# work backwards from the end
i <- nrow(aerial2)
while(i > 1){
  thisN <- aerial2$N[i]
  thisD <- aerial2$doy[i]
  if(thisN > 0){
    # divide these redds backwards in time
    # find next index
    j <- i-1
    while(j > 1){
      if(aerial2$N[j] > 0){
        break
      } else {j <- j-1 }
      # limit to 7 days
      if(i-j > 7) break
    }
    gap <- i - j
    k <- i
    remains <- thisN
    while(k > j){
      if(remains < 1){
        break;
      } else {
        aerial3$N[k] <- round(remains/(k-j) + 0.5)
        remains <- remains - round(remains/(k-j) + 0.5) 
        k <- k - 1
      }
    }
    i <- k
  } else {i <- i-1}
  
}
aerial3$year  <- as.numeric(aerial3$year)
aerial3$cumN <- NULL

for(y in allyears){
  aerial3$cumN[aerial3$year == y] <- cumsum(aerial3$N[aerial3$year == y]) 
}

LCL <- UCL <- LCL95 <- UCL95 <- MedianAer <- rep(NA,length(allyears))
for(year in allyears){
  i <- year-min(allyears) +1
  mx <- max(aerial3$cumN[aerial3$year==year])
  LCL[i] <- min(aerial3$doy[aerial3$year==year & aerial3$cumN >= 0.1*mx])
  UCL[i] <- max(aerial3$doy[aerial3$year==year & aerial3$cumN <= 0.9*mx])
  # median is a little tricky. want the mean day if the median balue is on multiple days
  # Find maximum day below the median and minimum day above the edian and average the gap
  low <- max(aerial3$doy[aerial3$year==year & aerial3$cumN <= 0.5*mx])
  hi <- min(aerial3$doy[aerial3$year==year & aerial3$cumN >= 0.5*mx])
  MedianAer[i] <- mean(low,hi)
  LCL95[i] <- min(aerial3$doy[aerial3$year==year & aerial3$cumN >= 0.025*mx])
  UCL95[i] <- max(aerial3$doy[aerial3$year==year & aerial3$cumN <= 0.975*mx])
}

#==== process carcass ====
# Don't need gaps, because we smooth 3 days twice
allspawnCar <- NULL
for(y in allyears){ 
  rr <- carcassexpand[[as.character(y)]]
  rrbyday <- as_tibble(rr) %>%  group_by(Day) %>% transmute(N = sum(across(starts_with("RKM"))))
  allspawnCar[[as.character(y)]] <- rrbyday
}

carcass1 <- carcass2 <- NULL
for(year in allyears){
  carcass1 <- rbind.data.frame(carcass1,cbind.data.frame(year=year,allspawnCar[[as.character(year)]]))
  carcass2 <- rbind.data.frame(carcass2,cbind.data.frame(year=rep(year,365),doy=1:365,N=0))
}
for(i in 1:nrow(carcass1)){
  I <- carcass2$year == carcass1$year[i] & carcass2$doy==carcass1$Day[i]; 
  carcass2$N[I] <- carcass1$N[i];
}
carcass3.1 <- carcass3 <- carcass2

for(year in allyears){
  # Smooth by 3 days twice
  I <- carcass2$year == year
  carcass3$N[I][-length(carcass2$N[I])] <- zrunavg(carcass2$N[I],3)[-1]
  carcass3.1$N[I][-length(carcass3$N[I])] <- zrunavg(carcass3$N[I],3)[-1]
  carcass3.1$cumN[I] <- cumsum(carcass3.1$N[I]) 
}

MedianCar <- LCLc <- UCLc <- LCLc95 <- UCLc95 <- inferredCarRedds <- rep(NA,length(allyears))
for(year in allyears){
  i <- year-min(allyears)+1
  mx <- max(carcass3.1$cumN[carcass3.1$year==year])
  inferredCarRedds[i] <- mx
  LCLc[i] <- min(carcass3.1$doy[carcass3.1$year==year & carcass3.1$cumN >= 0.1*mx])
  UCLc[i] <- max(carcass3.1$doy[carcass3.1$year==year & carcass3.1$cumN <= 0.9*mx])
  # median is a little tricky. want the mean day if the median value is on multiple days
  # Find maximum day below the median and minimum day above the median and average the gap
  low <- max(carcass3.1$doy[carcass3.1$year==year & carcass3.1$cumN <= 0.5*mx])
  hi <- min(carcass3.1$doy[carcass3.1$year==year & carcass3.1$cumN >= 0.5*mx])
  MedianCar[i] <- mean(low,hi)
  LCLc95[i] <- min(carcass3$doy[carcass3.1$year==year & carcass3.1$cumN >= 0.025*mx])
  UCLc95[i] <- max(carcass3$doy[carcass3.1$year==year & carcass3.1$cumN <= 0.975*mx])
}

#==== additional diagnostics ====
results <-cbind.data.frame("year"=allyears,"inferredCarRedds"=inferredCarRedds, "MedianAerial"=MedianAer,"MedianCarcass"=MedianCar,
                           "DiffsMedian"=MedianAer-MedianCar, "Aerial95"=UCL95-LCL95,"Carcass95"=UCLc95-LCLc95,
                           "Diffs95"= (UCL95-LCL95) - (UCLc95-LCLc95))
print(results)
summary(abs(results$DiffsMedian))
summary(abs(results$Diffs95))

print(cbind.data.frame(annualtotalsaerial,annualtotalscarcass$raw,results$inferredCarRedds,annualtotalscarcass$raw/annualtotalsaerial$raw))
summary(cbind.data.frame(annualtotalsaerial,annualtotalscarcass$raw,results$inferredCarRedds,annualtotalscarcass$raw/annualtotalsaerial$raw))
# SD inferredCarRedds
print(sqrt(var(results$inferredCarRedds)))
print(sqrt(var(annualtotalsaerial$raw)))

print(cbind(allyears,LCL95,UCL95,LCLc95,UCLc95,MedianAer,MedianCar))
write.csv(cbind(allyears,RangeAerial=results$Aerial95,RangeCarcass=results$Carcass95,MedianAer,MedianCar),file="results.csv",row.names=F)

#==== Draw figure for  ms. Four Panel Figure====

for(fig in  "pdf"){  #  c("pdf","png")){   # ""){    #  
     cexaxis <- 1 ; cexoutside <- 0.8 ;  cexpanel <-  0.8
   if(fig == "pdf"){
     pdf(file="bothdistrib.pdf",width=7,height=9)
     cexaxis <- 1 ; cexoutside <- 0.8 ;  cexpanel <-  0.8
  }
  if(fig == "png") png(file="bothdistrib.png",width=700,height=900)
  par(mfrow=c(1,4))
  # For shading
  usecol <- rgb(10,10,10,4,NULL,25)
  usescale=60
  coll1 <- coll2 <- NULL
  coll1 <- rgb(5,10,5,15,NULL,25)
  coll2 <- rgb(10,5,10,15,NULL,25)
  {
  want <- 2014:2023 ; I <- aerial3$year >= min(want) & aerial3$year <= max(want)
  zridges(scale=usescale,maxyr=max(want),minyr=min(want),type=fig,use1 =aerial3[I,] ,t1 = "", colo=coll1)
  polygon(c(LCL[match(want,allyears)],rev(UCL[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  polygon(c(LCL95[match(want,allyears)],rev(UCL95[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  lines(MedianAer[match(want,allyears)],want)
      # y axis a little more complex. repeat scaling across all year only on first of layout
      for(k in min(want):max(want)){
        scaleaxis <- pretty(0:usescale); scaleaxis <- scaleaxis[-length(scaleaxis)]
        yats <- k+scaleaxis/usescale
        axis(2,at=yats,label=scaleaxis,las=2,lwd=0.5,lwd.ticks=1,col.ticks="black",col.axis="black",mgp=c(3,0,-1),cex.axis=cexaxis)
      }
  mtext("(a) Aerial Surveys",line=-0.5,adj=1,cex=cexpanel)

  want <- 2004:2013 ; I <- aerial3$year >= min(want) & aerial3$year <= max(want)
  zridges(scale=usescale,maxyr=max(want),minyr=min(want),type=fig,use1 =aerial3[I,] ,t1 = "",colo=coll1)
  polygon(c(LCL[match(want,allyears)],rev(UCL[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  polygon(c(LCL95[match(want,allyears)],rev(UCL95[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  lines(MedianAer[match(want,allyears)],want)
  
  want <- 2014:2023 ; I <- carcass3$year >= min(want) & carcass3$year <= max(want)
  zridges(scale=usescale,maxyr=max(want),minyr=min(want),type=fig,use1 = carcass3[I,] ,t1 ="", colo=coll2)
  polygon(c(LCLc[match(want,allyears)],rev(UCLc[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  polygon(c(LCLc95[match(want,allyears)],rev(UCLc95[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  lines(MedianCar[match(want,allyears)],want)
  
  mtext("(b) Carcass Surveys",line=-0.5,adj=1,cex=cexpanel)
  
  want <- 2004:2013 ; I <- carcass3$year >= min(want) & carcass3$year <= max(want)
  zridges(scale=usescale,maxyr=max(want),minyr=min(want),type=fig,use1 = carcass3[I,] ,t1 ="", colo=coll2)
  polygon(c(LCLc[match(want,allyears)],rev(UCLc[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  polygon(c(LCLc95[match(want,allyears)],rev(UCLc95[match(want,allyears)])),c(want,rev(want)),col=usecol,border=usecol)
  lines(MedianCar[match(want,allyears)],want)
  }
  mtext("Estimated Spawning Date",side=1,outer=TRUE,line=-2,cex=cexoutside)
  if(fig != "") dev.off()
}

#==== Supplement figure with each year separated ====
for(fig in c("pdf")){   #  ""){    #  
  if(fig == "pdf") pdf(file="yeardetails.pdf",width=7,height=7)
  if(fig == "png") png(file="yeardetails.png",width=800,height=800)
  
for(year in 2004:2023) { # 2023){  # 2019){   #  
  par(mfrow=c(2,2))  
  for(j in 1:4){
    # cycle through the smooth, raww aerial and caracss combos nad plot
    par(mar=c(4,3,0,0),mgp=c(2.2,0.7,0))
    if(j==1){use1 <- aerial3;t1 <- paste("(a)",year,"smoothed\naerial survey")}
    if(j==2){use1 <- aerial2;t1 <- paste("(b)",year,"raw\naerial survey")}
    if(j==3){use1 <- carcass3;t1 <- paste("(c)",year,"smoothed\ncarcass survey")}
    if(j==4){use1 <- carcass2;t1 <- paste("(d)",year,"raw\ncarcass survey")}
  if(j ==3 | j==4){
    par(mar=c(6,3,0,0))
  }
  cexx <- 1.2
  if(fig=="png")cexx <- par()$cex*1.8
  if(fig=="pdf")cexx <- par()$cex*1
  cexaxis <- 1.3 *cexx
  
  xlimm <- c(108,249)
  histcolor <- rgb(.20,.20,.60,1)
  
  plot(0,0,xlim=xlimm,ylim=c(year,year+1),axes=FALSE,xlab="",ylab="",cex=cexx,xaxs="i", yaxs="i");
  abline(v=c(121,152,182,213),lwd=1,lty=2,col="grey90")
  use2 <- use1[use1$year==year,]
  scalar <- max(c(aerial2$N[aerial2$year==year],carcass2$N[carcass2$year==year]),na.rm=TRUE)
  
    vals <-pretty(0:scalar)
    # vals <- vals[vals > 0]
    vals <- vals[-length(vals)]
    ats <- year + vals/scalar
    abline(h=ats,lwd=1,lty=2,col="grey90")
    if(j==1 | j==3){  
      axis(2,at=ats,labels=vals, las=2,lwd=1,lwd.ticks=1,cex.axis=cexaxis) ;  
      axis(2,at=ats,labels=vals, las=2,lwd=1,lwd.ticks=1,col="white",col.ticks="black",col.axis="white",line=3) ; 
    } # draw Y axis. Needs some scaling to the year.
    # box()
  zdraw(use2,scaleit=scalar,hist=TRUE,colo=histcolor,xlimm=xlimm)
  text(130,year+.9,t1,cex=cexaxis,adj=0)
  # if(j==1)text(120,year+1,year,cex=cexx)

  if(j==3 | j ==4){ 
    axis(1,at=c(121,152,182,213,244),labels=c("May 1 ","Jun 1 ","Jul 1 ","Aug 1 ","Sep 1  "),
       outer=FALSE,cex.axis=cexaxis,lwd=1,lwd.ticks=1)
  }
  } 
  I <- 118:243 
 # par(mfrow=c(1,1))
 # plot(carcass2$doy[carcass2$year==year][I],cumsum(carcass2$N[carcass2$year==year][I]),type="l",lwd=8,col="grey80")
 # lines(carcass3$doy[carcass2$year==year][I],cumsum(carcass3$N[carcass2$year==year][I]),type="l",lwd=4,col="darkgreen")
 # lines(carcass3.1$doy[carcass2$year==year][I],cumsum(carcass3.1$N[carcass2$year==year][I]),type="l",lwd=2,col="pink")
 # 
  mtext("Estimated Spawning Date",side=1,outer=TRUE,line=-2,cex=cexaxis)
  
}

if(fig != "")dev.off()
}

