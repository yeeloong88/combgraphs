# Testing combgraphs function
library(msm)
library(adegenet)

######################################################################
# Getting Dummy Distance matrix
######################################################################
set.seed(24)
nrand=10
gen.dat <- round(c(rtnorm(nrand,mean=0,sd=1,lower=0),rtnorm(nrand,mean=100,sd=1,lower=0),rtnorm(nrand,mean=10,sd=1,lower=0)))
tem.dat <- round(c(rtnorm(nrand,mean=0,sd=2,lower=0),rtnorm(nrand,mean=25,sd=2,lower=0),rtnorm(nrand,mean=50,sd=2,lower=0)))
spa.dat <- round(c(rtnorm(nrand,mean=10,sd=10,lower=0),rtnorm(nrand,mean=200,sd=50,lower=0),rtnorm(nrand,mean=100,sd=30,lower=0)))
dat=data.frame(id=1:length(gen.dat),grp=c(rep("a",nrand),rep("b",nrand),rep("c",nrand)),genetic=gen.dat,temporal=tem.dat,spatial=spa.dat)

gen.dist <- as.matrix(dist(gen.dat))
tem.dist <- as.matrix(dist(tem.dat))
spa.dist <- as.matrix(dist(spa.dat))

#Clustered graphs
gen.g <- gengraph.matrix(gen.dist,ngrp=5)
tem.g <- gengraph.matrix(tem.dist,ngrp=5)  
spa.g <- gengraph.matrix(spa.dist,ngrp=5) 
dat$gen.clust <- gen.g$clust$membership
dat$tem.clust <- tem.g$clust$membership
dat$spa.clust <- spa.g$clust$membership

#Theoretical grouping
par(mfrow=c(1,3))
plot(dat$genetic,pch=20,col=c("blue","red","green")[dat$grp],ylab="Genetic Distance")
plot(dat$temporal,pch=20,col=c("blue","red","green")[dat$grp],ylab="Temporal Distance")
plot(dat$spatial,pch=20,col=c("blue","red","green")[dat$grp],ylab="Spatial Distance")
par(mfrow=c(1,1))


######################################################################
#Testing combgraph function
######################################################################
newlist <- list(gen.dist,tem.dist,spa.dist)
newlist2 <- list(dist(gen.dat),dist(tem.dat),dist(spa.dat))
combgraphs.list(newlist,ngrp=5,plotAll=T)   #Testing if plotAll=T works --> should replicate gen.g, tem.g, spa.g
system.time(combgraphs.list(newlist,ngrp=5,plotAll=F,method="union"))  #How long the function takes to generate union of graphs
system.time(combgraphs.list(newlist,ngrp=5,plotAll=F,method="intersect"))  #How long the function takes to generate intersect of graphs

# Plotting all 5 graphs for comparison
par(mfrow=c(2,3))
combgraphs.list(newlist,ngrp=5,plotAll=T) #Select "1", then "2", then "3", and quit. The "union" graph is also plotted
combgraphs.list(newlist,ngrp=5,plotAll=F,method="intersect")
par(mfrow=c(1,1))

g<-combgraphs.list(newlist,ngrp=5,plotAll=F)$graph
tkplot(g,layout=layout.circle)   #interactive interface to shift vertices around




#######################
## Using ToyOutbreak ##
#######################



######################################################################
######################################################################
# Examples: How to get Distance matrices?
######################################################################
######################################################################
#Genetic distance matrix
dna <- fasta2DNAbin("alignment.fa")   #requires "ape"
D <- dist.dna(dna,model="N")

#Temporal distance matrix
# Temporal data, x, is the number of days since the oldest dated entry in the data
dates <- as.Date(cases$collec.dates)
days <- as.integer(difftime(dates,min(dates),unit="days"))
temporal.x <- dist(dates)

#Spatial distance matrix
filename <- readShapeSpatial("",IDvar="code")  #import shapefile
pts.filename <- fortify(filename,region="code")		#convert shapefile into dataframe
filename.df <- merge(pts.filename, filename, by.x="id",by.y="code")		#merge data

gDistance(filename,byid=TRUE) #Package:rgeos; cartesian distance for geometric objects 	
coord <- coordinate(filename)   #Package: sp
temp2 <- dist1(coord)   #Package: SpatialTools; Euclidean distance matrix of coordinates