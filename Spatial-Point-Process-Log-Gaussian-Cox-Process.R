library(maptools)
library(sp)
library(spatstat)
library(spatstat.data)
library(dplyr)
library(RandomFields)
library(RandomFieldsUtils)
#**************************input data Indonesia***************************************#
Data<-read.csv("E:\\Magister ITS\\THESIS\\Data Gempa\\Data input R 2009 - 2018 M5.csv",header = T)
gunungapi<-read.csv("E:\\Magister ITS\\THESIS\\Data Gempa\\gunung api Indonesia.csv",header=T)
#deteksi point yang sama#
data2<-distinct(Data,latitude,longitude,.keep_all=TRUE)
data2

#***********************************Data Sulawesi Maluku**********************************#
#Data point gempa#
ii<-with(data2,which(latitude>5.7,arr.ind=TRUE))
datasul1<-data2[-ii,]
jj<-with(datasul1,which(latitude<(-7.5),arr.ind=TRUE))
datasul2<-datasul1[-jj,]
kk<-with(datasul2,which(longitude<118,arr.ind=TRUE))
datasul3<-datasul2[-kk,]
ll<-with(datasul3,which(longitude>136,arr.ind=TRUE))
datasulawesimaluku<-datasul3[-ll,]

#Data Kovariat Sulawesi Maluku#
#==============================
sesarSM<-read.csv("E:\\Magister ITS\\THESIS\\Data Gempa\\sesar maluku sulawesi input r.csv",header=T)
subducSM<-read.csv("E:\\Magister ITS\\THESIS\\Data Gempa\\Subduksi Sulawesi Maluku.csv",header=T)
#pilih gunung sulawesi#
ism<-with(gunungapi,which(Latitude>5.7,arr.ind=TRUE))
gunungsul1<-gunungapi[-ism,]
jsm<-with(gunungsul1,which(Latitude<(-7.5),arr.ind=TRUE))
gunungsul2<-gunungsul1[-jsm,]
ksm<-with(gunungsul2,which(Longitude<118,arr.ind=TRUE))
gunungsul3<-gunungsul2[-ksm,]
lsm<-with(gunungsul3,which(Longitude>136,arr.ind=TRUE))
gunungapiSM<-gunungsul3

#**************************Membuat Window Sulawesi Maluku***************************#
windowgempaSM<-owin(xrange=c(118,136),yrange=c(-7.5,5.7))
area.window<-area.owin(windowgempaSM)

#***********************Membuat Planar Point Pattern gempa bumi********************#
long<-datasulawesimaluku$longitude
lat<-datasulawesimaluku$latitude

####*****Maps World******###
library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(119,136), ylim = c(-7.5,5.7), asp=0)
points(datasulawesimaluku$longitude, datasulawesimaluku$latitude, col = "red", cex = 0.6)

#### PPP Gempa Sulawesi Maluku ####
gempaSM<-ppp(long,lat,window = windowgempaSM) #point pattern gempa sulawesi#
plot(newmap, xlim = c(119,136), ylim = c(-7.5,5.7), asp=0)
plot(gempaSM, add=TRUE,cex=1)

#**********************Membuat PSP Variabel Kovariat Sulawesi************************#
##Gunung Api Aktif Sulawesi Maluku##
latgSM<-gunungapiSM$Latitude
longgSM<-gunungapiSM$Longitude
vulcanoSM<-ppp(longgSM,latgSM,window = windowgempaSM)

##Sesar Aktif Sulawesi Maluku##
x0s<-sesarSM$x0
y0s<-sesarSM$y0
x1s<-sesarSM$x1
y1s<-sesarSM$y1
faultSM<-psp(y0s,x0s,y1s,x1s,window = windowgempaSM)

##Subduksi Sulawesi Maluku##
x0sub<-subducSM$x0
y0sub<-subducSM$y0
x1sub<-subducSM$x1
y1sub<-subducSM$y1
subductionSM<-psp(y0sub,x0sub,y1sub,x1sub,window = windowgempaSM)

#*******************Eksploratory Data Sulawesi Maluku*****************************####
windowgempaSM1<-owin(xrange=c(118*111,136*111),yrange=c(-7.5*111,5.7*111))
area.window<-area.owin(windowgempaSM)
nSM<-gempaSM$n #banyaknya gempa pada window sulawesi Maluku#
intensitySM<-nSM/area.window #intensitas (rho) gempa sulawesi Maluku#

### Stationarity Test Sulawesi Maluku###
gridgempaSM<-quadratcount(gempaSM,8,13)
ujihomogenitas<-quadrat.test(gridgempaSM)
ujihomogenitas

### Uji Distribution Poisson
#Uji Anderson Darling#
library(goftest)
datatest=as.data.frame(ygempaSM)
adtest = ad.test(ygempaSM,"ppois",lambda=mean(ygempaSM), estimated = T)
adtest

##Plot Grid Gempa Sulawesi Maluku##
#plot(gempaSM)
#plot(gridgempaSM, add=T)
#datagridSM=as.data.frame(gridgempaSM)
#ygempaSM=cbind(datagridSM$Freq)

### Independency Test ###
ind.test<-Kinhom(gempaSM,normpower = 2, nlarge = 2000, correction = "border")
independence_test=ind.test
plot(independence_test)
plot(ind.test)

### Deteksi Model Mixture ##
densitygempa<-density(gridgempaSM)
plot(hist(densitygempa, col = "red"))
plot(densitygempa, add=T)
hist(ygempaSM, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "temp",
     ylim = c(0,0.08),
     xlim = c(0,80),
     main = "Density Gempa")
lines(density(ygempaSM), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")

######################################################################################
#                           Function Jarak Kovariat                                  #
######################################################################################
## Fungsi jarak Fault ##
gempa.dfault<-with(gempaSM,distfun(faultSM))
plot(gempa.dfault)

## Fungsi Jarak Subduction ##
gempa.dsub<-with(gempaSM,distfun(subductionSM))
plot(gempa.dsub)

## Fungsi Jarak Volcano ##
gempa.dvol<-with(gempaSM,distfun(vulcanoSM))
plot(gempa.dvol)

##Data Kovariat Jarak##
#Jarak Points ke Faults Terdekat
#cara1
zf<-distfun(faultSM)
dzf<-zf(gempaSM)
#cara2
d.pfSM<-project2segment(gempaSM,faultSM)
x1<-d.pfSM$d
#Jarak Points ke Subduction Terdekat
#cara1
zs<-distfun(subductionSM)
dzs<-zs(gempaSM)
#cara2
d.psSM<-project2segment(gempaSM,subductionSM)
x2<-d.psSM$d
#Jarak Points ke Gunung Api Terdekat
#cara1
zv<-distfun(vulcanoSM)
dzv<-zv(gempaSM)
#cara2
d.pvSM<-nncross(gempaSM,vulcanoSM)
x3<-d.pvSM$dist

#Transformasi Pixel Kovariat
#jarak points ke faults
marksSM<-cbind.data.frame(x1,x2,x3)
gempaSMmarks<-ppp(long,lat,window = windowgempaSM, marks = marksSM)
Z<-Smooth.ppp(gempaSMmarks, sigma=2, at="pixels")

#Plot Contour gempa
pppcontour<-contour(density(gempaSM), nlevels=30, labcex=0.7)
#plot contour fault
x1marks<-ppp(long,lat,window = windowgempaSM, marks = marksSM$x1)
plot.ppp(x1marks, cols="red", border=F, main=NULL)
plot(newmap, xlim = c(118,136), ylim = c(-7.5,5.7), asp=0, add=T)
contour.ssf(x1marks, sigma=0.5)

#plot contour subduction
x2marks<-ppp(long,lat,window = windowgempaSM, marks = marksSM$x2)
plot.ppp(x2marks, cols="red", border=F)
plot(newmap, xlim = c(119,136), ylim = c(-7.5,5.7), asp=0, add=T)
contour.ssf(x2marks, sigma=0.5)

#plot contour volcano
x3marks<-ppp(long,lat,window = windowgempaSM, marks = marksSM$x3)
plot.ppp(x3marks, cols="red", border=F)
plot(newmap, xlim = c(119,136), ylim = c(-7.5,5.7), asp=0, add=T)
contour.ssf(x3marks, sigma=0.5)

#plot smooth fault
plot(Z$x1)
plot(faultSM,add=TRUE)
#plot smooth subduksi
plot(Z$x2)
plot(subductionSM,add=TRUE)
#plot smooth gunung
plot(Z$x3)
plot(vulcanoSM,add=TRUE)
plot.ppp(gempaSMmarks)
smothxmarks_point<-Smooth.ppp(gempaSMmarks, sigma=2, at="points")
mat_imSMx1<-as.matrix.im(Z$x1)
mat_imSMx2<-as.matrix.im(Z$x2)
mat_imSMx3<-as.matrix.im(Z$x3)
write.csv(smothxmarks_point,'E:\\Magister ITS\\THESIS\\Data Gempa\\MatriksMarksatPoints.csv',row.names = FALSE)
write.csv(mat_imSMx1,'E:\\Magister ITS\\THESIS\\Data Gempa\\MatriksMarksX11.csv',row.names = FALSE)
write.csv(mat_imSMx2,'E:\\Magister ITS\\THESIS\\Data Gempa\\MatriksMarksX21.csv',row.names = FALSE)
write.csv(mat_imSMx3,'E:\\Magister ITS\\THESIS\\Data Gempa\\MatriksMarksX31.csv',row.names = FALSE)


##########################################
        #===========================
        #### GEMPA SM LGCP ####
        #===========================
##########################################
library(INLA)
library(Matrix)
library(parallel)

###########################
#BASED ON DISTANCE
##==========================
## Covariate gempa SM ##
#########################
#####################
faultcov1 = Z$x1
subcov2 = Z$x2
volcov3 = Z$x3

# Choosing grid resolution, here using a VERY coarse resolution to speed up calculations
nrow = 8
ncol = 13
n = nrow*ncol

# Vector giving the area for each grid cell
AreaSM = rep(max(windowgempaSM$xrange)*max(windowgempaSM$yrange)/n, n)

# Discretize the observation window
x.gridSM = quadrats(windowgempaSM, ncol, nrow)

# Count the number of points in each grid cell, starting at upper left corner and downwards
ySM = as.vector(quadratcount(gempaSM, tess = x.gridSM))

func.centersSM <- function(windowgempaSM, nrow, ncol)
{
  eps = c((max(windowgempaSM$x)-min(windowgempaSM$x))/ncol/2, (max(windowgempaSM$y)-min(windowgempaSM$y))/nrow/2)
  x.x = matrix(rep(seq(min(windowgempaSM$x)+eps[1], max(windowgempaSM$x)-eps[1],length=ncol),nrow),nrow,ncol,byrow=T)
  x.y = matrix(rep(seq(max(windowgempaSM$y)-eps[2], min(windowgempaSM$y)+eps[2],length=nrow),ncol),nrow,ncol)
  x.ppp = ppp(x.x, x.y, window=windowgempaSM)
  return(x.ppp)
}
centersSM = func.centersSM(windowgempaSM, nrow, ncol)
faultcov1 = lookup.im(faultcov1, centersSM$x, centersSM$y)
subcov2 = lookup.im(subcov2, centersSM$x, centersSM$y)
volcov3 = lookup.im(volcov3, centersSM$x, centersSM$y)
y = ySM
faultcov1 = scale(c(log(faultcov1)))
subcov2 = scale(c(log(subcov2)))
volcov3 = scale(c(log(volcov3)))

# Prepare dataset
data = list(y=ySM, cov1=faultcov1, cov2=subcov2, cov3=volcov3, index=seq(1:n))

# Choose parameters for the PC priors for hyperparameters
#run inla 1 sigma
sig1=1
alpha.sig1=0.01
u.phi1=0.5
alpha.phi1=2/3
form1=y ~ cov1 + cov2 + cov3 +
  f(index, nrow=nrow, ncol=ncol, model="rw2diid", scale.model=T,
    hyper=list(prec=list(prior="pc.prec", param=c(sig1, alpha.sigma)),
               phi =list(prior="pc", param=c(U.phi, alpha.phi))))

result1 = inla(form1, family="poisson", data=data, E=AreaSM, 
               verbose=F, control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE))

sumrres1=summary(result1)
sumrres1

#model without insignificant covariate
######################################
form2=y ~ cov1 + cov3 +
  f(index, nrow=nrow, ncol=ncol, model="rw2diid", scale.model=T,
    hyper=list(prec=list(prior="pc.prec", param=c(sig1, alpha.sigma)),
               phi =list(prior="pc", param=c(U.phi, alpha.phi))))

result2 = inla(form2, family="poisson", data=data, E=AreaSM, 
               verbose=F, control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE))

sumrres2=summary(result2)
sumrres2

#Plot Density
tmpb0 = inla.tmarginal(function(x) x, result1$marginals.fixed[[1]]) 
plot(tmpb0, type = "l", xlab = paste("Fixed effect marginal", i, ":", result1$names.fixed[1]), ylab = "Density")

tmpb1 = inla.tmarginal(function(x) x, result1$marginals.fixed[[2]]) 
plot(tmpb1, type = "l", xlab = paste("Fixed effect marginal", i, ":", result1$names.fixed[2]), ylab = "Density")

tmpb2 = inla.tmarginal(function(x) x, result1$marginals.fixed[[3]]) 
plot(tmpb2, type = "l", xlab = paste("Fixed effect marginal", i, ":", result1$names.fixed[3]), ylab = "Density")

tmpb3 = inla.tmarginal(function(x) x, result1$marginals.fixed[[4]]) 
plot(tmpb3, type = "l", xlab = paste("Fixed effect marginal", i, ":", result1$names.fixed[4]), ylab = "Density")

tmppres = inla.tmarginal(function(x) x, result1$marginals.hyperpar$`Precision for index`) 
plot(tmppres, type = "l", xlab = paste("Fixed effect marginal", i, ":", "precision"), ylab = "Density")

tmpphi = inla.tmarginal(function(x) x, result1$marginals.hyperpar$`Phi for index`) 
plot(tmpphi, type = "l", xlab = paste("Fixed effect marginal", i, ":", "phi"), ylab = "Density")

tmpb0_2 = inla.tmarginal(function(x) x, result2$marginals.fixed[[1]]) 
plot(tmpb0_2, type = "l", xlab = paste("Fixed effect marginal", i, ":", result1$names.fixed[1]), ylab = "Density")

tmpb1_2 = inla.tmarginal(function(x) x, result2$marginals.fixed[[2]]) 
plot(tmpb1_2, type = "l",xlab = paste("Fixed effect marginal", i, ":", result1$names.fixed[2]), ylab = "Density")

mm=result2$marginals.fixed$cov1
plot(inla.smarginal(mm))
tmpb3_2 = inla.smarginal(function(x) exp(x), result2$marginals.fixed[[3]]) 
plot(tmpb3_2, type = "l", xlab = paste("Fixed effect marginal", i, ":", result1$names.fixed[4]), ylab = "Density")

tmppres_2 = inla.tmarginal(function(x) x, result1$marginals.hyperpar$`Precision for index`) 
plot(tmppres_2, type = "l", xlab = paste("Fixed effect marginal", i, ":", "precision"), ylab = "Density")

tmpphi_2 = inla.tmarginal(function(x) x, result1$marginals.hyperpar$`Phi for index`) 
plot(tmpphi_2, type = "l", xlab = paste("Fixed effect marginal", i, ":", "phi"), ylab = "Density")

################################
#Initial Value prior GLM Poisson
################################
library(flexmix)
y = ySM
x1 = faultcov1
x2 = subcov2
x3 = volcov3
ypoiss<-glm(y~1+x1+x2+x3, family = poisson("log"))
summary(ypoiss)

#GLM Mixture Poisson
ypoiss_mix<-flexmix(y~x1+x2+x3,k=2, model=FLXglm(family="poisson"))
summary(ypoiss_mix)
rmodel<-refit(ypoiss_mix)
summary(rmodel)
#################################
