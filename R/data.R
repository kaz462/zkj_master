#### read data
setwd("data")

I=83
J=14
unemploy <- matrix(0, nrow = I, ncol = J)
mi.unem90 <- read.csv("laucnty90.csv")
unemploy[, 1] <- as.vector(mi.unem90$Unemployment_rate[mi.unem90$state == 26])
mi.unem91 <- read.csv("laucnty91.csv")
unemploy[, 2] <- as.vector(mi.unem91$Unemployment_rate[mi.unem91$state == 26])
mi.unem92 <- read.csv("laucnty92.csv")
unemploy[, 3] <- as.vector(mi.unem92$Unemployment_rate[mi.unem92$state == 26])
mi.unem93 <- read.csv("laucnty93.csv")
unemploy[, 4] <-  as.vector(mi.unem93$Unemployment_rate[mi.unem93$state == 26])
mi.unem94 <- read.csv("laucnty94.csv")
unemploy[, 5] <-  as.vector(mi.unem94$Unemployment_rate[mi.unem94$state == 26])
mi.unem95 <- read.csv("laucnty95.csv")
unemploy[, 6] <-  as.vector(mi.unem95$Unemployment_rate[mi.unem95$state == 26])
mi.unem96 <- read.csv("laucnty96.csv")
unemploy[, 7] <-  as.vector(mi.unem96$Unemployment_rate[mi.unem96$state == 26])
mi.unem98 <- read.csv("laucnty98.csv")
unemploy[, 8] <-  as.vector(mi.unem98$Unemployment_rate[mi.unem98$state == 26])
mi.unem99 <- read.csv("laucnty99.csv")
unemploy[, 9] <-  as.vector(mi.unem99$Unemployment_rate[mi.unem99$state == 26])
mi.unem00 <- read.csv("laucnty00.csv")
unemploy[, 10] <- as.vector(mi.unem00$Unemployment_rate[mi.unem00$state == 26])
mi.unem01 <- read.csv("laucnty01.csv")
unemploy[, 11] <- as.vector(mi.unem01$Unemployment_rate[mi.unem01$state == 26])
mi.unem02 <- read.csv("laucnty02.csv")
unemploy[, 12] <- as.vector(mi.unem02$Unemployment_rate[mi.unem02$state == 26])
mi.unem03 <- read.csv("laucnty03.csv")
unemploy[, 13] <- as.vector(mi.unem03$Unemployment_rate[mi.unem03$state == 26])
mi.unem04 <- read.csv("laucnty04.csv")
unemploy[, 14] <- as.vector(mi.unem04$Unemployment_rate[mi.unem04$state == 26])

mi.gdl <- read.csv("mi_GDL1617.csv")
mi.gdl.96 <- subset(mi.gdl, YEAR == 1996)
mi.county <- map("county", "michigan", fill = T, plot = F)
mi.IDs <- sapply(strsplit(mi.county$names, ":"), function(x)  x[1])

## Here we transform the dataset in the maps database into a SpatialPolygon object
mi.poly.county <- map2SpatialPolygons(mi.county,IDs=mi.IDs,proj4string=CRS("+proj=longlat + datum=wgs84"))
car.values <- addSpatial(mi.poly.county)
adj <- car.values$adj
num <- car.values$num
sumNumNeigh <- sum(num)

injury <- read.table("injury.txt",header=T)
mi.gdl <- read.table("mi_GDL1617.txt",header=T)
mi.gdl <- as.data.frame(mi.gdl)
y=O<-matrix(0,nrow=I,ncol=J)
x1 <- matrix(0,nrow=I,ncol=J)
x2 <- matrix(0,nrow=I,ncol=J)
x3 <- matrix(0,nrow=I,ncol=J)
T1 <- matrix(0,nrow=I,ncol=J)
T2 <- matrix(0,nrow=I,ncol=J)
Time <- matrix(0,nrow=I,ncol=J)

for(i in 1:I){
  Time[i,] <- c(rep(1,7),rep(0,7))
}

for(i in 1:I){
  y[i,] <- injury$LOGFATAL[which(injury$COUNTY==2*i-1)][c(1:7,9:15)]
  O[i,] <- mi.gdl$count_tn[which(mi.gdl$COUNTY==2*i-1)][c(1:7,9:15)]
}

for(j in 1:J){
  x1[,j] <- as.numeric(unemploy[,j])/100
}

for(i in 1:I){
  x2[i,] <- injury$RURALITY[which(injury$COUNTY==2*i-1)][c(1:7,9:15)]
}

for(i in 1:I){
  T1[i,] <- injury$T1[which(injury$COUNTY==2*i-1)][c(1:7,9:15)]
}
for(i in 1:I){
  T2[i,] <- injury$T2[which(injury$COUNTY==2*i-1)][c(1:7,9:15)]
}

#Measurement error-real data
mi_teen_pop <- read.table("mi_teen_pop.txt",header=T)
mi.teen <- as.data.frame(mi_teen_pop)
popT <- matrix(0,nrow=20,ncol=1)
for(i in 1:20){
  popT[i,1] <- sum(mi.gdl$pop_tn[mi.gdl$YEAR==(1990+i-1)])
}
rate <- matrix(0,nrow=20,ncol=1)
for(i in 1:20){
  rate[i,1] <- mi.teen[i,2]/popT[i,1]
}
r <- log(rate)[c(1:7,9:15)]
rate <- rate[c(1:7,9:15)]

d_rate <- matrix(0,nrow=I,ncol=J)
for(i in 1:I){   
  for(j in 1:J){	
    d_rate[i,j] <- rnorm(1,mean=rate[j],sd=0.005)
  }
}
dd_rate <- as.vector(d_rate)

setwd("..")

## transfer function 
mungeCARdata4stan <- function(adjBUGS,numBUGS) {
  N = length(numBUGS);
  nn = numBUGS;
  N_edges = length(adjBUGS) / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  iAdj = 0;
  iEdge = 0;
  for (i in 1:N) {
    for (j in 1:nn[i]) {
      iAdj = iAdj + 1;
      if (i < adjBUGS[iAdj]) {
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adjBUGS[iAdj];
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}
nbs <- mungeCARdata4stan(adj, num);

## get the car parameters used in stan 
N <- nbs$N;
node1 <- nbs$node1;
node2 <- nbs$node2;
N_edges <- nbs$N_edges;
node_1 <- cbind(node1,node2)
node_2 <- cbind(node2,node1)
W <- matrix(0, I, I)
W[node_1] <- 1
W[node_2] <- 1
D <- diag(num)
W_n <- length(node1)
Eigen <- eigen(solve(D)%*%W) #eigen(solve(D)^(0.5)%*%W%*%solve(D)^(0.5))  0<a<1


yy <- as.vector(O) #original count data
county <-rep(1:83,14)
xx1 <- as.vector(x1)
xx2 <- c(rep(x2[,1],14))
TT1 <- as.vector(-T1)
TT2 <- as.vector(T2)
TTime <- c(rep(1,7*83),rep(0,7*83))
TTime1 <- c(TT1[1:(1162/2)],TT2[(1162/2+1):1162])
x_lambda <- as.matrix(cbind(rep(1,83*14),scale(xx1),scale(xx2)))
x_theta <- as.matrix(cbind(rep(1,83*14),scale(xx1)))
pop_tn <- mi.gdl$pop_tn[mi.gdl$YEAR!=1997 & mi.gdl$YEAR!=2005 & mi.gdl$YEAR!=2006 & mi.gdl$YEAR!=2007 & mi.gdl$YEAR!=2008 & mi.gdl$YEAR!=2009]
