#######################################################################################################################################
######################### 	STORM PETREL DENSITY ESTIMATION ON FILFLA ISLAND 	################################################
#######################################################################################################################################

## written by steffen.oppel@rspb.org.uk on 13 July 2013
## modified for Malta on 3 January 2014


################################################################################################################
###################### LOAD PACKAGES AND CUSTOM SCRIPTS   ######################################################
################################################################################################################


library(RODBC)
library(secr)
library(foreign)
library(sp)
library(maptools)
library(rgdal)
library(SPACECAP)
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Statistics\\SPACECAP.r")
source("C:\\STEFFEN\\RSPB\\Statistics\\SPACECAP.r")
source("C:\\STEFFEN\\RSPB\\Statistics\\UTM_conversion_old.r")
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Statistics\\UTM_conversion_old.r")




################################################################################################################
###################### LOAD AND MANIPULATE DATA 	  ######################################################
################################################################################################################

####################### LOADING RINGING DATA           #########################################################
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Malta\\Raw_data")
setwd("C:\\STEFFEN\\RSPB\\Malta\\Raw_data")

SP<-odbcConnectAccess2007('ESP_Malta.accdb')
ESP <- sqlQuery(SP, "SELECT * FROM SECR_captures")  ## changed query to only have data from 2000
trap <- sqlQuery(SP, "SELECT * FROM SECR_input_detectors")
odbcClose(SP)


### CONVERT COORDINATES INTO UTM COORDINATES
trap[,20:21]<-UTM_conversion(trap$LAT, trap$LONG,33)



####################### LOADING DATA FOR CREATING MASKS #########################################################
##### these are polygons created specifically for SECR in GIS
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Malta\\Analysis\\ESP_abundance_survival")
setwd("C:\\STEFFEN\\RSPB\\Malta\\Analysis\\ESP_abundance_survival")

data2<-readOGR("Filfla.kml", layer="Filfla.kml")
plot(data2)
proj4string(data2)<-CRS("+proj=longlat +ellps=WGS84")
filfla<-spTransform(data2, CRS=CRS("+proj=utm +zone=33 +ellps=WGS84 +units=m"))
#writeOGR(obj=filfla,dsn="Filfla.shp",driver="ESRI Shapefile",layer="Filfla")
Filflapoly <- readShapePoly("Filfla.shp") 



#######################################################################################
############### Data Preparation for 'secr'   #########################################
#######################################################################################
# NOTE: for 'secr' the input file must be in this order: island, Ring_nr, sampling_occasion, Net_ID
### ARRANGE DETECTOR INPUT FILE FOR 'secr'###

trap<-trap[order(trap$Net_ID, decreasing=F),c(1,21,20,5:19)]
trap[is.na(trap)]<-0


#### LOOK AT DATA ###########
trap
head(ESP)


#######################################################################################
############### Data Preparation for SPACECAP #########################################
#######################################################################################

### IMPORTANT: FOR SPACECAP INPUT OF DETECTORS MUST BE NUMERIC, alphanumeric is not allowed!!!
### TRAP COORDINATES MUST BE IN UTM!!
### DEPLOYMENT DETAILS IN SEPARATE COLUMNS
### HEADERS MUST BE INCLUDED AND EXACTLY LIKE IN HELP FILE!!!!!
### ORDER OF COLUMNS IS IMPORTANT AND MUST BE: detector, animal ID, occasion


#### INPUT FILE 1: Animal Capture Details
#### requires Ring_Nr to be sequential from 1:N, with N=number of individuals 

ESP_spacecap<-ESP[,c(4,2,3)]
names(ESP_spacecap)<-c('LOC_ID','ANIMAL_ID','SO')
birds<-levels(as.factor(ESP_spacecap$ANIMAL_ID))
ESP_spacecap$ANIMAL_ID<-match(ESP_spacecap$ANIMAL_ID,birds) 			## this replaces the original CAT_ID with the position along the vector of unique CAT_IDs
head(ESP_spacecap)


#### INPUT FILE 2: Trap Deployment Details
#### Each net location must be given a unique identification number, ranging from 1: J, where J is the total number of camera trap locations used in the survey
#### Sampling occasions must be numbered from 1:x
trap_SPACECAP<-trap
names(trap_SPACECAP)[1:3]<-c('LOC_ID','X_Coord','Y_Coord')
trap_SPACECAP[,4:dim(trap_SPACECAP)[2]]<-ifelse(trap_SPACECAP[,4:dim(trap_SPACECAP)[2]]>0,1,0)



#### INPUT FILE 3: Potential Home-Range Centers (=Nest Sites)

## create home-range center input files
minX<-min(trap[,2]-150)
minY<-min(trap[,3]-150)
maxX<-max(trap[,2]+150)
maxY<-max(trap[,3]+150)
a<-seq(minX,maxX, 5)
b<-seq(minY,maxY,5)
MH <- expand.grid(x=a, y=b)

## turn the grid into a spatial points object
MHgsps = SpatialPoints(MH)

### TRUNCATE GRID TO REMOVE GRID POINTS THAT ARE OUTSIDE THE TRAPPING AREAS ####

inside<-over(MHgsps, Filflapoly)		### this does the actual overlay and produces a vector of 1 and NA
MH[,3]<-ifelse(inside$Name=="Filfla",1,0)
MH<-data.frame(X_COORD=as.integer(MH[,1]),Y_COORD=as.integer(MH[,2]),HABITAT=as.integer(MH[,3]))
#MH<-subset(MH, HABITAT==1)
MH[is.na(MH)]<-0
head(MH)

#### CHECK WHETHER OVERLAY WAS SUCCESSFUL
plot(Filflapoly)
points(MH[,1:2], pch=3, col=(MH$HABITAT+1))
points(trap[,2:3], pch=1, col="green")

save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Malta\\Analysis\\ESP_abundance_survival\\ESP_SPACECAP_INPUT.RData")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Malta\\Analysis\\ESP_abundance_survival\\ESP_SPACECAP_INPUT.RData")


################################################################################################################
###################### SAVE THE DATA IN THE PROPER FORMAT ######################################################
################################################################################################################

write.table(ESP, "ESP_secr_input.txt", row.names=F, col.names=F, sep=" ")
write.table(trap, "ESP_trap_locations_secr.txt", row.names=F, col.names=F, sep=" ", quote=F)

write.table(ESP_spacecap,"ESP_captures_SPACECAP.csv", sep=',', row.names=F)
write.table(trap_SPACECAP,"traps_SPACECAP.csv", row.names=F, sep=",")
write.table(MH, "Pot_nests.csv", row.names=F, sep=",")





################################################################################################################
###################### SECR ANALYSIS ###########################################################################
################################################################################################################


### IMPORTING DATA (the ones you just created by calling the database) ####
# NOTE: the bird data input file must be in this order: island, ring_NR, sampling_occasion, net
# the two files must be in the directory you specified above

ESP_CH <- read.capthist("ESP_secr_input.txt", "ESP_trap_locations_secr.txt", detector="proximity", verify=T, binary.usage = FALSE)
summary(ESP_CH)




####################### CREATE MASKS AND PLOT THE DATA #####################################################


### CREATE A HABITAT MASK TO REDUCE THE AREA OVER WHICH PROBABILITIES ARE CALCULATED ###
trap<- traps(ESP_CH)   #extracts the trap layout from the capthist object and turns it into a traps object
mask<-make.mask(trap, type='trapbuffer', spacing=10, buffer=80, poly=Filflapoly)

plot(ESP_CH, tracks=T, varycol=FALSE)   # each colour represents an individual, tracks=T specifies that you want lines to connect subsequent sightings of the same individual at different detectors
plot(mask, add=T)
plot(trap, add=T)


### CHECK FOR POPULATION CLOSURE

closure.test(ESP_CH, SB=T)


### CONVERT TO COUNT DETECTOR INPUT IF CLOSURE TEST INDICATES OPEN POPULATION

#ESP.countCH<-reduce(ESP_CH, columns = list(1:20), outputdetector = 'count')
#summary(ESP.countCH)



########################################################################################################################################
############################################## 'secr' analysis       ###################################################################
########################################################################################################################################
### removed hazard rate because too unrealistic estimates

halfnorm.poisson<-secr.fit (ESP_CH, buffer=80, detectfn=0, model=list(D~1, g0~1, sigma~1), details=list(binomN=0), trace=F, mask=mask)  #assumes poisson distribution
exp.poisson<-secr.fit (ESP_CH, buffer=80, detectfn=2, model=list(D~1, g0~1, sigma~1), details=list(binomN=0), trace=F, mask=mask)  #assumes poisson distribution

# AIC table to determine the best of these simple models
model_table1<-AIC(exp.poisson, halfnorm.poisson)
model_table1  ## detectfn and details=list(binomN=...) from top model must then be passed on to next model comparison


### OUTPUT FROM THE TOP MODEL
simple.summary<-derived(halfnorm.poisson) #lists the effective sampling area and the Density estimate


# Model averaged density estimates
model_averaged_estimates<-as.data.frame(model.average(halfnorm.poisson, halfnorm.poisson))   ### provides model-averaged real parameter estimates
model_averaged_estimates<-model_averaged_estimates[,1:4]
model_averaged_estimates
plogis(log(model_averaged_estimates[2,1]))	## reduce.capthist causes default link to switch to 'log' rather than 'logit', hence back-transformed parameter estimates for g0 are >1
plogis(log(model_averaged_estimates[2,3]))
plogis(log(model_averaged_estimates[2,4]))


### ESTIMATE TOTAL NUMBER OF MADEIRAN STORM PETRELS BY MULTIPLYING DENSITY FIGURE WITH AREA
island_area<-3.84			### directly derived from ArcGIS
total_ESP_estimates<-model_averaged_estimates[,1:4]*island_area
for (c in 1:4){
total_ESP_estimates[2,c]<-plogis(log(model_averaged_estimates[2,c]))
}
total_ESP_estimates[3,]<-(qchisq(0.95,2)^0.5)*(model_averaged_estimates[3,])		### equation from Noss et al. 2012 - result gives home range radius
total_ESP_estimates<-rbind(total_ESP_estimates, simple.summary[,1:4])
total_ESP_estimates


write.table(total_ESP_estimates,"ESP_Filfla_abundance_estimate.csv", sep=",", row.names=F)
write.table(simple.summary,"ESP_densities.csv", sep=",", row.names=F)


#### PLOT DETECTION FUNCTION #######
### capture probabilities derived from a normal distribution with sigma, and scaled up to the peak which = g0
pdf("ESP_detection_function.pdf", height=5, width=5)
par(mar=c(5,6,0,0), oma=c(0,0,0,0))
capt_prob<-(dnorm(0:100,mean=0,sd=model_averaged_estimates[3,1]))*(total_ESP_estimates[2,1]/dnorm(0,mean=0,sd=model_averaged_estimates[3,1]))
plot(c(0:100),capt_prob, type='l', ylim=c(0,total_ESP_estimates[2,4]+0.002),axes=F, xlab="Distance from mistnet (m)", ylab="Capture probability", cex.lab=1.4, cex.axis=1.4, bty='n', mgp=c(4,0.5,0))
par(new=T)
capt_prob<-(dnorm(0:100,mean=0,sd=model_averaged_estimates[3,3]))*(total_ESP_estimates[2,1]/dnorm(0,mean=0,sd=model_averaged_estimates[3,3]))
plot(c(0:100),capt_prob, type='l', ylim=c(0,total_ESP_estimates[2,4]+0.002),lty=2, axes=F, xlab="", ylab="", cex.lab=1.4, cex.axis=1.4, bty='n')
par(new=T)
capt_prob<-(dnorm(0:100,mean=0,sd=model_averaged_estimates[3,4]))*(total_ESP_estimates[2,1]/dnorm(0,mean=0,sd=model_averaged_estimates[3,4]))
plot(c(0:100),capt_prob, type='l', ylim=c(0,total_ESP_estimates[2,4]+0.002),lty=2, axes=F, xlab="", ylab="", cex.lab=1.4, cex.axis=1.4, bty='n')
axis(1, at=seq(0,100,10), cex.axis=1.4)
axis(2, at=seq(0,total_ESP_estimates[2,4]+0.002, 0.005), labels=T, cex.axis=1.4, las=1, mgp=c(4,0.5,0))
dev.off()



save.image("C:\\STEFFEN\\RSPB\\Malta\\Analysis\\ESP_abundance_survival\\ESP_SECR.RData")


#####################################################################################################################################################################################
###### PLOT A DENSITY MAP #######################################
#####################################################################################################################################################################################
map<-secr.fit (ESP_CH, buffer=80, detectfn=0, model=list(D~x, g0~1, sigma~1), details=list(binomN=0), trace=F, mask=mask, method='BFGS')  #assumes poisson distribution

MASPmap <- predictDsurface(map, cl.D = TRUE)
plot(MASPmap, plottype = 'shaded', polycol = 'blue', border = 100)
plot(traps(ESP_CH), detpar = list(col = 'black'), add = TRUE)
head(MASPmap)




#######################################################################################
############### Analysis with SPACECAP ################################################
#######################################################################################

SPACECAP()


# area = 0.0025
# niter=80000
# burnin=20000
# thin=3
# N=150



out<-read.table("S:/ConSci/DptShare/SteffenOppel/RSPB/Malta/Analysis/ESP_Abundance_survival/output_132429/param_val_132429.csv", header=T, sep=",")
hist(out$Nsuper, breaks=100, xlim=c(10000,40000), xlab="Number of European Storm Petrels on Filfla", ylab="Probability", main="")


