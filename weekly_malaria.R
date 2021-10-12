library(dlnm)
library(imputeTS)
library(AER)
library(ggplot2)
library(RCurl)

##El Bagre
url_path_ElB = "https://raw.githubusercontent.com/juandavidgutier/weekly_malaria_and_environmental_variables/master/ElBagre.csv"
malaria_ElB <- read.csv(url_path_ElB)
malaria_ElB$incidence <- malaria_ElB$Cases*10000/malaria_ElB$total_population 


#imputation for NAs in EVI
statsNA(malaria_ElB$EVI)
malaria_ElB$EVIimputed <- na.kalman(malaria_ElB$EVI, model = "StructTS")
g = ggplot_na_imputations(malaria_ElB$EVI, malaria_ElB$EVIimputed, legend_size = 5, size_points = 2.5, size_imputations = 3.5,
                          title = "El Bagre",  subtitle = "")
g+theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        title=element_text(size=25))

#EVI*10
for (i in 1:length(malaria_ElB$EVIimputed)) {
        malaria_ElB$EVIimputed[i] <- malaria_ElB$EVIimputed[i]/1000
}

#ts
incidence <- malaria_ElB$incidence
ts_incidence <- ts(malaria_ElB$incidence, start = c(2007,1), frequency = 52)
plot(ts_incidence, main="A",  xlab="Time", ylab="Incidence")

temp <- malaria_ElB$Temperature
ts_temp <- ts(malaria_ElB$Temperature, start = c(2007,1), frequency = 52)
rain <- malaria_ElB$Rainfall
ts_rain <- ts(malaria_ElB$Rainfall, start = c(2007,1), frequency = 52)
EVI <- malaria_ElB$EVIimputed
ts_EVI <- ts(malaria_ElB$EVIimputed, start = c(2007,1), frequency = 52)


#crossbasis funtions
cb.temp_ElB <- crossbasis(malaria_ElB$Temperature, lag=14, argvar=list(fun="bs", degree=2), 
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.rain_ElB <- crossbasis(malaria_ElB$Rainfall, lag=14, argvar=list(fun="bs", degree=2),
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.EVI_ElB <- crossbasis(malaria_ElB$EVIimputed, lag=14, argvar=list(fun="bs", degree=2), 
                     arglag=list(fun="ns",intercept=FALSE)) 

model_ElB <- glm(Cases ~ cb.temp_ElB + cb.rain_ElB + cb.EVI_ElB + lag(Cases, k=1) + year + epiweek,
             family=quasipoisson(), malaria_ElB)


#EVI
pred.EVI_ElB <- crosspred(cb.EVI_ElB, model_ElB, by=0.1, bylag=0.5, cen=0, cumul = FALSE, at=0:4)

#surface
plot(pred.EVI_ElB, xlab="EVI index *10", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_EVI_ElB <- crossreduce(cb.EVI_ElB, model_ElB, type="lag", value=1, cen=0, at=0:4)
plot(predlag1_EVI_ElB, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag4_EVI_ElB <- crossreduce(cb.EVI_ElB, model_ElB, type="lag", value=4, cen=0, at=0:4)
plot(predlag4_EVI_ElB, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag12_EVI_ElB <- crossreduce(cb.EVI_ElB, model_ElB, type="lag", value=12, cen=0, at=0:4)
plot(predlag12_EVI_ElB, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)


#####with values of EVI
pred0.5_EVI_ElB <- crossreduce(cb.EVI_ElB, model_ElB, type="var", value=0.5, cen=0, at=0:4)
plot(pred0.5_EVI_ElB, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred2.5_EVI_ElB <- crossreduce(cb.EVI_ElB, model_ElB, type="var", value=2.5, cen=0, at=0:4)
plot(pred2.5_EVI_ElB, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred4_EVI_ElB <- crossreduce(cb.EVI_ElB, model_ElB, type="var", value=4, cen=0, at=0:4)
plot(pred4_EVI_ElB, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#temperature
pred.temp_ElB <- crosspred(cb.temp_ElB, model_ElB, by=0.5, bylag=0.5, cen=25, cumul = FALSE, at=23.55:27.5)

#surface
plot(pred.temp_ElB, xlab="Temperature °C", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_temp_ElB <- crossreduce(cb.temp_ElB, model_ElB, type="lag", value=1, cen=25, at=23.55:27.5)
plot(predlag1_temp_ElB, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag4_temp_ElB <- crossreduce(cb.temp_ElB, model_ElB, type="lag", value=4, cen=25, at=23.55:27.5)
plot(predlag4_temp_ElB, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag12_temp_ElB <- crossreduce(cb.temp_ElB, model_ElB, type="lag", value=12, cen=25, at=23.55:27.5)
plot(predlag12_temp_ElB, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)


######with values de temperature
pred23_temp_ElB <- crossreduce(cb.temp_ElB, model_ElB, type="var", value=23.55, cen=25, at=23.55:27.5)
plot(pred23_temp_ElB, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred25_temp_ElB <- crossreduce(cb.temp_ElB, model_ElB, type="var", value=25.5, cen=25, at=23.55:27.5)
plot(pred25_temp_ElB, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred27_temp_ElB <- crossreduce(cb.temp_ElB, model_ElB, type="var", value=27.5, cen=25, at=23.55:27.5)
plot(pred27_temp_ElB, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#rainfall
#prediction rainfall
pred.rain_ElB <- crosspred(cb.rain_ElB, model_ElB, by=10, bylag=0.5, cumul=TRUE, cen = 5, at=5:180)


#surface
plot(pred.rain_ElB, xlab="Rainfall (mm/week)", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_rain_ElB <- crossreduce(cb.rain_ElB, model_ElB, type="lag", value=1, cen=5, at=5:180)
plot(predlag1_rain_ElB, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag4_rain_ElB <- crossreduce(cb.rain_ElB, model_ElB, type="lag", value=4, cen=5, at=5:180)
plot(predlag4_rain_ElB, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag12_rain_ElB <- crossreduce(cb.rain_ElB, model_ElB, type="lag", value=12, cen=5, at=5:180)
plot(predlag12_rain_ElB, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)


######with values of rainfall
pred20_rain_ElB <- crossreduce(cb.rain_ElB, model_ElB, type="var", value=20, cen=5, at=5:180)
plot(pred20_rain_ElB, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred100_rain_ElB <- crossreduce(cb.rain_ElB, model_ElB, type="var", value=100, cen=5, at=5:180)
plot(pred100_rain_ElB, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred180_rain_ElB <- crossreduce(cb.rain_ElB, model_ElB, type="var", value=180, cen=5, at=5:180)
plot(pred180_rain_ElB, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)



###################################################################################
##Tierraalta
url_path_Tie = "https://raw.githubusercontent.com/juandavidgutier/weekly_malaria_and_environmental_variables/master/Tierralta.csv"
malaria_Tie <- read.csv(url_path_Tie)
malaria_Tie$incidence <- malaria_Tie$Cases*10000/malaria_Tie$total_population

#imputation for NAs in EVI
statsNA(malaria_Tie$EVI)
malaria_Tie$EVIimputed <- na.kalman(malaria_Tie$EVI, model = "StructTS")
g = ggplot_na_imputations(malaria_Tie$EVI, malaria_Tie$EVIimputed, legend_size = 5, size_points = 2.5, size_imputations = 3.5,
                          title = "Tierralta",  subtitle = "")
g+theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        title=element_text(size=25))

#EVI*10
for (i in 1:length(malaria_Tie$EVIimputed)) {
        malaria_Tie$EVIimputed[i] <- malaria_Tie$EVIimputed[i]/1000
}


#ts
incidence <- malaria_Tie$incidence
ts_incidence <- ts(malaria_Tie$incidence, start = c(2007,1), frequency = 52)
plot(ts_incidence, main="B",  xlab="Time", ylab="Incidence")


temp <- malaria_Tie$Temperature
ts_temp <- ts(malaria_Tie$Temperature, start = c(2007,1), frequency = 52)
rain <- malaria_Tie$Rainfall
ts_rain <- ts(malaria_Tie$Rainfall, start = c(2007,1), frequency = 52)
EVI <- malaria_Tie$EVIimputed
ts_EVI <- ts(malaria_Tie$EVIimputed, start = c(2007,1), frequency = 52)



#crossbasis funtions
cb.temp_Tie <- crossbasis(malaria_Tie$Temperature, lag=14, argvar=list(fun="bs", degree=2),  
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.rain_Tie <- crossbasis(malaria_Tie$Rainfall, lag=14, argvar=list(fun="bs", df=2), 
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.EVI_Tie <- crossbasis(malaria_Tie$EVIimputed, lag=14, argvar=list(fun="bs", df=2), 
                     arglag=list(fun="ns",intercept=FALSE)) 


model_Tie <- glm(Cases ~ cb.temp_Tie + cb.rain_Tie + cb.EVI_Tie + lag(Cases, k=1) + year + epiweek,
              family=quasipoisson, malaria_Tie)


#EVI
pred.EVI_Tie <- crosspred(cb.EVI_Tie, model_Tie, by=0.1, bylag=0.5, cen=0, cumul = FALSE, at=0:4.5)

#surface
plot(pred.EVI_Tie, xlab="EVI index *10", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_EVI_Tie <- crossreduce(cb.EVI_Tie, model_Tie, type="lag", value=1, cen=0, at=0:4.5)
plot(predlag1_EVI_Tie, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag4_EVI_Tie <- crossreduce(cb.EVI_Tie, model_Tie, type="lag", value=4, cen=0, at=0:4.5)
plot(predlag4_EVI_Tie, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag12_EVI_Tie <- crossreduce(cb.EVI_Tie, model_Tie, type="lag", value=12, cen=0, at=0:4.5)
plot(predlag12_EVI_Tie, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)


######with values of EVI
pred0.5_EVI_Tie <- crossreduce(cb.EVI_Tie, model_Tie, type="var", value=0.5, cen=0, at=0:4)
plot(pred0.5_EVI_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred2.5_EVI_Tie <- crossreduce(cb.EVI_Tie, model_Tie, type="var", value=2.5, cen=0, at=0:4)
plot(pred2.5_EVI_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred4_EVI_Tie <- crossreduce(cb.EVI_Tie, model_Tie, type="var", value=4, cen=0, at=0:4)
plot(pred4_EVI_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#temperature
pred.temp_Tie <- crosspred(cb.temp_Tie, model_Tie, by=0.5, bylag=0.5, cen=25, cumul = FALSE, at=23.55:27.5)

#surface
plot(pred.temp_Tie, xlab="Temperature °C", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_temp_Tie <- crossreduce(cb.temp_Tie, model_Tie, type="lag", value=1, cen=25, at=23.55:27.5)
plot(predlag1_temp_Tie, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag4_temp_Tie <- crossreduce(cb.temp_Tie, model_Tie, type="lag", value=4, cen=25, at=23.55:27.5)
plot(predlag4_temp_Tie, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag12_temp_Tie <- crossreduce(cb.temp_Tie, model_Tie, type="lag", value=12, cen=25, at=23.55:27.5)
plot(predlag12_temp_Tie, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)


######with values of temperature
pred23_temp_Tie <- crossreduce(cb.temp_Tie, model_Tie, type="var", value=23.55, cen=25, at=23.55:27.5)
plot(pred23_temp_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred25_temp_Tie <- crossreduce(cb.temp_Tie, model_Tie, type="var", value=25.5, cen=25, at=23.55:27.5)
plot(pred25_temp_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred27_temp_Tie <- crossreduce(cb.temp_Tie, model_Tie, type="var", value=27.28, cen=25, at=23.55:27.5)
plot(pred27_temp_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#rainfall
#prediction rainfall
pred.rain_Tie <- crosspred(cb.rain_Tie, model_Tie, by=20, bylag=0.5, cen = 5, cumul = TRUE, at=5:180)

#surface
plot(pred.rain_Tie, xlab="Rainfall (mm/week)", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_rain_Tie <- crossreduce(cb.rain_Tie, model_Tie, type="lag", value=1, cen=5, at=5:180)
plot(predlag1_rain_Tie, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag4_rain_Tie <- crossreduce(cb.rain_Tie, model_Tie, type="lag", value=4, cen=5, at=5:180)
plot(predlag4_rain_Tie, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag12_rain_Tie <- crossreduce(cb.rain_Tie, model_Tie, type="lag", value=12, cen=5, at=5:180)
plot(predlag12_rain_Tie, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)


######with values of rainfall
pred20_rain_Tie <- crossreduce(cb.rain_Tie, model_Tie, type="var", value=20, cen=5, at=5:180)
plot(pred20_rain_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred100_rain_Tie <- crossreduce(cb.rain_Tie, model_Tie, type="var", value=100, cen=5, at=5:180)
plot(pred100_rain_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred180_rain_Tie <- crossreduce(cb.rain_Tie, model_Tie, type="var", value=180, cen=5, at=5:180)
plot(pred180_rain_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)




####################################################################################
##Pto Libertador
url_path_PtoL = "https://raw.githubusercontent.com/juandavidgutier/weekly_malaria_and_environmental_variables/master/PtoLibertador.csv"
malaria_PtoL <- read.csv(url_path_PtoL)
malaria_PtoL$incidence <- malaria_PtoL$Cases*10000/malaria_PtoL$total_population 


#imputation for NAs in EVI
statsNA(malaria_PtoL$EVI)
malaria_PtoL$EVIimputed <- na.kalman(malaria_PtoL$EVI, model = "StructTS")
g = ggplot_na_imputations(malaria_PtoL$EVI, malaria_PtoL$EVIimputed, legend_size = 5, size_points = 2.5, size_imputations = 3.5,
                          title = "Puerto Libertador",  subtitle = "")
g+theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        title=element_text(size=25))

#EVI*10
for (i in 1:length(malaria_PtoL$EVIimputed)) {
        malaria_PtoL$EVIimputed[i] <- malaria_PtoL$EVIimputed[i]/1000
}


#ts
incidence <- malaria_PtoL$incidence
ts_incidence <- ts(malaria_PtoL$incidence, start = c(2007,1), frequency = 52)
plot(ts_incidence, main="C",  xlab="Time", ylab="Incidence")


temp <- malaria_PtoL$Temperature
ts_temp <- ts(malaria_PtoL$Temperature, start = c(2007,1), frequency = 52)
rain <- malaria_PtoL$Rainfall
ts_rain <- ts(malaria_PtoL$Rainfall, start = c(2007,1), frequency = 52)
EVI <- malaria_PtoL$EVIimputed
ts_EVI <- ts(malaria_PtoL$EVIimputed, start = c(2007,1), frequency = 52)


#crossbasis funtions
cb.temp_PtoL <- crossbasis(malaria_PtoL$Temperature, lag=14, argvar=list(fun="bs", degree=2),  
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.rain_PtoL <- crossbasis(malaria_PtoL$Rainfall, lag=14, argvar=list(fun="bs", degree=2), 
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.EVI_PtoL <- crossbasis(malaria_PtoL$EVIimputed, lag=14, argvar=list(fun="bs", degree=2),
                     arglag=list(fun="ns",intercept=FALSE)) 


model_PtoL <- glm(Cases ~ cb.temp_PtoL + cb.rain_PtoL + cb.EVI_PtoL + lag(Cases, k=1) + year + epiweek,
              family=quasipoisson, malaria_PtoL)


#EVI
pred.EVI_PtoL <- crosspred(cb.EVI_PtoL, model_PtoL, by=0.1, bylag=0.5, cen=0, cumul = FALSE, at=0:4.5)

#surface
plot(pred.EVI_PtoL, xlab="EVI index *10", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_EVI_PtoL <- crossreduce(cb.EVI_PtoL, model_PtoL, type="lag", value=1, cen=0, at=0:4.5)
plot(predlag1_EVI_PtoL, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag4_EVI_PtoL <- crossreduce(cb.EVI_PtoL, model_PtoL, type="lag", value=4, cen=0, at=0:4.5)
plot(predlag4_EVI_PtoL, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag12_EVI_PtoL <- crossreduce(cb.EVI_PtoL, model_PtoL, type="lag", value=12, cen=0, at=0:4.5)
plot(predlag12_EVI_PtoL, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)


######with values of EVI
pred0.5_EVI_PtoL <- crossreduce(cb.EVI_PtoL, model_PtoL, type="var", value=0.5, cen=0, at=0:4)
plot(pred0.5_EVI_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred2.5_EVI_PtoL <- crossreduce(cb.EVI_PtoL, model_PtoL, type="var", value=2.5, cen=0, at=0:4)
plot(pred2.5_EVI_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred4_EVI_PtoL <- crossreduce(cb.EVI_PtoL, model_PtoL, type="var", value=4, cen=0, at=0:4)
plot(pred4_EVI_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#temperature
pred.temp_PtoL <- crosspred(cb.temp_PtoL, model_PtoL, by=0.5, bylag=0.5, cen=25, cumul = FALSE, at=23.55:27.5)

#surface
plot(pred.temp_PtoL, xlab="Temperature °C", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_temp_PtoL <- crossreduce(cb.temp_PtoL, model_PtoL, type="lag", value=1, cen=25, at=23.55:27.5)
plot(predlag1_temp_PtoL, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag4_temp_PtoL <- crossreduce(cb.temp_PtoL, model_PtoL, type="lag", value=4, cen=25, at=23.55:27.5)
plot(predlag4_temp_PtoL, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag12_temp_PtoL <- crossreduce(cb.temp_PtoL, model_PtoL, type="lag", value=12, cen=25, at=23.55:27.5)
plot(predlag12_temp_PtoL, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)


######with values of temperature
pred23_temp_PtoL <- crossreduce(cb.temp_PtoL, model_PtoL, type="var", value=23.55, cen=25, at=23.55:27.5)
plot(pred23_temp_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred25_temp_PtoL <- crossreduce(cb.temp_PtoL, model_PtoL, type="var", value=25.5, cen=25, at=23.55:27.5)
plot(pred25_temp_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred27_temp_PtoL <- crossreduce(cb.temp_PtoL, model_PtoL, type="var", value=27.5, cen=25, at=23.55:27.5)
plot(pred27_temp_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#rainfall
#prediction rainfall
pred.rain_PtoL <- crosspred(cb.rain_PtoL, model_PtoL, by=10, bylag=0.5, cen = 5, cumul = TRUE, at=5:180)

#surface
plot(pred.rain_PtoL, xlab="Rainfall (mm/week)", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_rain_PtoL <- crossreduce(cb.rain_PtoL, model_PtoL, type="lag", value=1, cen=5, at=5:180)
plot(predlag1_rain_PtoL, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag4_rain_PtoL <- crossreduce(cb.rain_PtoL, model_PtoL, type="lag", value=4, cen=5, at=5:180)
plot(predlag4_rain_PtoL, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag12_rain_PtoL <- crossreduce(cb.rain_PtoL, model_PtoL, type="lag", value=12, cen=5, at=5:180)
plot(predlag12_rain_PtoL, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)


######with values of rainfall
pred20_rain_PtoL <- crossreduce(cb.rain_PtoL, model_PtoL, type="var", value=20, cen=5, at=5:180)
plot(pred20_rain_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred100_rain_PtoL <- crossreduce(cb.rain_PtoL, model_PtoL, type="var", value=100, cen=5, at=5:180)
plot(pred100_rain_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred180_rain_PtoL <- crossreduce(cb.rain_PtoL, model_PtoL, type="var", value=180, cen=5, at=5:180)
plot(pred180_rain_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)




################################################################################
##Zaragoza
url_path_Zar = "https://raw.githubusercontent.com/juandavidgutier/weekly_malaria_and_environmental_variables/master/Zaragoza.csv"
malaria_Zar <- read.csv(url_path_Zar)
malaria_Zar$incidence <- malaria_Zar$Cases*10000/malaria_Zar$total_population 


#imputation for NAs in EVI
statsNA(malaria_Zar$EVI)
malaria_Zar$EVIimputed <- na.kalman(malaria_Zar$EVI, model = "StructTS")
g = ggplot_na_imputations(malaria_Zar$EVI, malaria_Zar$EVIimputed, legend_size = 5, size_points = 2.5, size_imputations = 3.5,
                          title = "Zaragoza",  subtitle = "")
g+theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        title=element_text(size=25))


#EVI*10
for (i in 1:length(malaria_Zar$EVIimputed)) {
        malaria_Zar$EVIimputed[i] <- malaria_Zar$EVIimputed[i]/1000
}

#ts
incidence <- malaria_Zar$incidence
ts_incidence <- ts(malaria_Zar$incidence, start = c(2007,1), frequency = 52)
plot(ts_incidence, main="D",  xlab="Time", ylab="Incidence")


temp <- malaria_Zar$Temperature
ts_temp <- ts(malaria_Zar$Temperature, start = c(2007,1), frequency = 52)
rain <- malaria_Zar$Rainfall
ts_rain <- ts(malaria_Zar$Rainfall, start = c(2007,1), frequency = 52)
EVI <- malaria_Zar$EVIimputed
ts_EVI <- ts(malaria_Zar$EVIimputed, start = c(2007,1), frequency = 52)


#crossbasis funtions
cb.temp_Zar <- crossbasis(malaria_Zar$Temperature, lag=14, argvar=list(fun="bs", degree=2), 
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.rain_Zar <- crossbasis(malaria_Zar$Rainfall, lag=14, argvar=list(fun="bs", degree=2), 
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.EVI_Zar <- crossbasis(malaria_Zar$EVIimputed, lag=14, argvar=list(fun="bs", degree=2), 
                     arglag=list(fun="ns",intercept=FALSE)) 


model_Zar <- glm(Cases ~ cb.temp_Zar + cb.rain_Zar + cb.EVI_Zar + lag(Cases, k=1) + year + epiweek,
              family=quasipoisson, malaria_Zar)


#EVI
pred.EVI_Zar <- crosspred(cb.EVI_Zar, model_Zar, by=0.1, bylag=0.5, cen=0, cumul = FALSE, at=0:4.5)

#surface
plot(pred.EVI_Zar, xlab="EVI index *10", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_EVI_Zar <- crossreduce(cb.EVI_Zar, model_Zar, type="lag", value=1, cen=0, at=0:4.5)
plot(predlag1_EVI_Zar, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag4_EVI_Zar <- crossreduce(cb.EVI_Zar, model_Zar, type="lag", value=4, cen=0, at=0:4.5)
plot(predlag4_EVI_Zar, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag12_EVI_Zar <- crossreduce(cb.EVI_Zar, model_Zar, type="lag", value=12, cen=0, at=0:4.5)
plot(predlag12_EVI_Zar, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)


######with values of EVI
pred0.5_EVI_Zar <- crossreduce(cb.EVI_Zar, model_Zar, type="var", value=0.5, cen=0, at=0:4)
plot(pred0.5_EVI_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred2.5_EVI_Zar <- crossreduce(cb.EVI_Zar, model_Zar, type="var", value=2.5, cen=0, at=0:4)
plot(pred2.5_EVI_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred4_EVI_Zar <- crossreduce(cb.EVI_Zar, model_Zar, type="var", value=4, cen=0, at=0:4)
plot(pred4_EVI_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#temperature
pred.temp_Zar <- crosspred(cb.temp_Zar, model_Zar, by=0.5, bylag=0.5, cen=25, cumul = FALSE, at=23.55:27.5)

#surface
plot(pred.temp_Zar, xlab="Temperature °C", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_temp_Zar <- crossreduce(cb.temp_Zar, model_Zar, type="lag", value=1, cen=25, at=23.55:27.5)
plot(predlag1_temp_Zar, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag4_temp_Zar <- crossreduce(cb.temp_Zar, model_Zar, type="lag", value=4, cen=25, at=23.55:27.5)
plot(predlag4_temp_Zar, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag12_temp_Zar <- crossreduce(cb.temp_Zar, model_Zar, type="lag", value=12, cen=25, at=23.55:27.5)
plot(predlag12_temp_Zar, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)


######with values of temperature
pred23_temp_Zar <- crossreduce(cb.temp_Zar, model_Zar, type="var", value=23.55, cen=25, at=23.55:27.5)
plot(pred23_temp_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred25_temp_Zar <- crossreduce(cb.temp_Zar, model_Zar, type="var", value=25.5, cen=25, at=23.55:27.5)
plot(pred25_temp_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred27_temp_Zar <- crossreduce(cb.temp_Zar, model_Zar, type="var", value=27.5, cen=25, at=23.55:27.5)
plot(pred27_temp_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#rainfall
#prediction rainfall
pred.rain_Zar <- crosspred(cb.rain_Zar, model_Zar, by=10, bylag=0.5, cen = 5, cumul = TRUE, at=5:180)

#surface
plot(pred.rain_Zar, xlab="Rainfall (mm/week)", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_rain_Zar <- crossreduce(cb.rain_Zar, model_Zar, type="lag", value=1, cen=5, at=5:180)
plot(predlag1_rain_Zar, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag4_rain_Zar <- crossreduce(cb.rain_Zar, model_Zar, type="lag", value=4, cen=5, at=5:180)
plot(predlag4_rain_Zar, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag12_rain_Zar <- crossreduce(cb.rain_Zar, model_Zar, type="lag", value=12, cen=5, at=5:180)
plot(predlag12_rain_Zar, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)


######with values of rainfall
pred20_rain_Zar <- crossreduce(cb.rain_Zar, model_Zar, type="var", value=20, cen=5, at=5:180)
plot(pred20_rain_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred100_rain_Zar <- crossreduce(cb.rain_Zar, model_Zar, type="var", value=100, cen=5, at=5:180)
plot(pred100_rain_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred180_rain_Zar <- crossreduce(cb.rain_Zar, model_Zar, type="var", value=180, cen=5, at=5:180)
plot(pred180_rain_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)





################################################################################
##Tumaco
url_path_Tum = "https://raw.githubusercontent.com/juandavidgutier/weekly_malaria_and_environmental_variables/master/Tumaco.csv"
malaria_Tum <- read.csv(url_path_Tum)
malaria_Tum$incidence <- malaria_Tum$Cases*10000/malaria_Tum$total_population 


#imputation for NAs in EVI
statsNA(malaria_Tum$EVI)
malaria_Tum$EVIimputed <- na.kalman(malaria_Tum$EVI, model = "StructTS")
g = ggplot_na_imputations(malaria_Tum$EVI, malaria_Tum$EVIimputed, legend_size = 5, size_points = 2.5, size_imputations = 3.5,
                          title = "Tumaco",  subtitle = "")
g+theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        title=element_text(size=25))

#EVI*10
for (i in 1:length(malaria_Tum$EVIimputed)) {
        malaria_Tum$EVIimputed[i] <- malaria_Tum$EVIimputed[i]/1000
}

#ts
incidence <- malaria_Tum$incidence
ts_incidence <- ts(malaria_Tum$incidence, start = c(2007,1), frequency = 52)
plot(ts_incidence, main="E",  xlab="Time", ylab="Incidence")


temp <- malaria_Tum$Temperature
ts_temp <- ts(malaria_Tum$Temperature, start = c(2007,1), frequency = 52)
rain <- malaria_Tum$Rainfall
ts_rain <- ts(malaria_Tum$Rainfall, start = c(2007,1), frequency = 52)
EVI <- malaria_Tum$EVIimputed
ts_EVI <- ts(malaria_Tum$EVIimputed, start = c(2007,1), frequency = 52)



#crossbasis funtions
cb.temp_Tum <- crossbasis(malaria_Tum$Temperature, lag=14, argvar=list(fun="bs", degree=2), 
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.rain_Tum <- crossbasis(malaria_Tum$Rainfall, lag=14, argvar=list(fun="bs", degree=2), 
                      arglag=list(fun="ns",intercept=FALSE)) 
cb.EVI_Tum <- crossbasis(malaria_Tum$EVIimputed, lag=14, argvar=list(fun="bs", degree=2), 
                     arglag=list(fun="ns",intercept=FALSE)) 


model_Tum <- glm(Cases ~ cb.temp_Tum + cb.rain_Tum + cb.EVI_Tum + lag(Cases, k=1) + year + epiweek,
              family=quasipoisson, malaria_Tum)


#EVI
pred.EVI_Tum <- crosspred(cb.EVI_Tum, model_Tum, by=0.1, bylag=0.5, cen=0, cumul = FALSE, at=0:4.5)

#surface
plot(pred.EVI_Tum, xlab="EVI index *10", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_EVI_Tum <- crossreduce(cb.EVI_Tum, model_Tum, type="lag", value=1, cen=0, at=0:4.5)
plot(predlag1_EVI_Tum, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag4_EVI_Tum <- crossreduce(cb.EVI_Tum, model_Tum, type="lag", value=4, cen=0, at=0:4.5)
plot(predlag4_EVI_Tum, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)

predlag12_EVI_Tum <- crossreduce(cb.EVI_Tum, model_Tum, type="lag", value=12, cen=0, at=0:4.5)
plot(predlag12_EVI_Tum, ci= "area", xlab="EVI index *10", ylab="RR", main="", col=2, lwd = 3)



######with values of EVI
pred0.5_EVI_Tum <- crossreduce(cb.EVI_Tum, model_Tum, type="var", value=0.5, cen=0, at=0:4)
plot(pred0.5_EVI_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred2.5_EVI_Tum <- crossreduce(cb.EVI_Tum, model_Tum, type="var", value=2.5, cen=0, at=0:4)
plot(pred2.5_EVI_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred4_EVI_Tum <- crossreduce(cb.EVI_Tum, model_Tum, type="var", value=4, cen=0, at=0:4)
plot(pred4_EVI_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#temperature
pred.temp_Tum <- crosspred(cb.temp_Tum, model_Tum, by=0.2, bylag=0.5, cen=25, cumul = FALSE, at=23.55:27.5)

#surface
plot(pred.temp_Tum, xlab="Temperature °C", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_temp_Tum <- crossreduce(cb.temp_Tum, model_Tum, type="lag", value=1, cen=25, at=23.55:27.5)
plot(predlag1_temp_Tum, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag4_temp_Tum <- crossreduce(cb.temp_Tum, model_Tum, type="lag", value=4, cen=25, at=23.55:27.5)
plot(predlag4_temp_Tum, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

predlag12_temp_Tum <- crossreduce(cb.temp_Tum, model_Tum, type="lag", value=12, cen=25, at=23.55:27.5)
plot(predlag12_temp_Tum, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)



######with values of temperature
pred23_temp_Tum <- crossreduce(cb.temp_Tum, model_Tum, type="var", value=23.55, cen=25, at=23.55:27.5)
plot(pred23_temp_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred25_temp_Tum <- crossreduce(cb.temp_Tum, model_Tum, type="var", value=25.5, cen=25, at=23.55:27.5)
plot(pred25_temp_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred27_temp_Tum <- crossreduce(cb.temp_Tum, model_Tum, type="var", value=27.5, cen=25, at=23.55:27.5)
plot(pred27_temp_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


#rainfall
#prediction rainfall
pred.rain_Tum <- crosspred(cb.rain_Tum, model_Tum, by=10, bylag=0.5, cen = 5, cumul = TRUE, at=5:180)

#surface
plot(pred.rain_Tum, xlab="Rainfall (mm/week)", zlab="RR", theta=200, phi=50, lphi=20, col = "gray",
     main="")

predlag1_rain_Tum <- crossreduce(cb.rain_Tum, model_Tum, type="lag", value=1, cen=5, at=5:180)
plot(predlag1_rain_Tum, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag4_rain_Tum <- crossreduce(cb.rain_Tum, model_Tum, type="lag", value=4, cen=5, at=5:180)
plot(predlag4_rain_Tum, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

predlag12_rain_Tum <- crossreduce(cb.rain_Tum, model_Tum, type="lag", value=12, cen=5, at=5:180)
plot(predlag12_rain_Tum, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)


######with values of rainfall
pred20_rain_Tum <- crossreduce(cb.rain_Tum, model_Tum, type="var", value=20, cen=5, at=5:180)
plot(pred20_rain_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred100_rain_Tum <- crossreduce(cb.rain_Tum, model_Tum, type="var", value=100, cen=5, at=5:180)
plot(pred100_rain_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

pred180_rain_Tum <- crossreduce(cb.rain_Tum, model_Tum, type="var", value=180, cen=5, at=5:180)
plot(pred180_rain_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)




####################################################################################################
#Figures 2-7
#The rows correspond to: El Bagre, Tierralta, Puerto Libertador, Zaragoza and Tumaco

Fig2 <- layout(matrix(1:15, ncol=3))
layout.show(Fig2)
par(mar = c(4, 5, 2, 2) + 0.1)

plot(pred0.5_EVI_ElB, ci= "area", xlab="lag", ylab="RR", main="EVI*10=0.5", col=2, lwd = 3)
plot(pred0.5_EVI_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred0.5_EVI_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred0.5_EVI_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred0.5_EVI_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

plot(pred2.5_EVI_ElB, ci= "area", xlab="lag", ylab="RR", main="EVI*10=2.5", col=2, lwd = 3)
plot(pred2.5_EVI_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred2.5_EVI_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred2.5_EVI_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred2.5_EVI_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

plot(pred4_EVI_ElB, ci= "area", xlab="lag", ylab="RR", main="EVI*10=4", col=2, lwd = 3)
plot(pred4_EVI_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred4_EVI_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred4_EVI_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred4_EVI_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


Fig3 <- layout(matrix(1:15, ncol=3))
layout.show(Fig3)
par(mar = c(4, 5, 2, 2) + 0.1)
plot(predlag1_EVI_ElB, ci= "area", xlab="EVI*10", ylab="RR", main="Lag=1", col=2, lwd = 3, lwd = 3)
plot(predlag1_EVI_Tie, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)
plot(predlag1_EVI_PtoL, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)
plot(predlag1_EVI_Zar, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)
plot(predlag1_EVI_Tum, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)

plot(predlag4_EVI_ElB, ci= "area", xlab="EVI*10", ylab="RR", main="Lag=4", col=2, lwd = 3)
plot(predlag4_EVI_Tie, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)
plot(predlag4_EVI_PtoL, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)
plot(predlag4_EVI_Zar, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)
plot(predlag4_EVI_Tum, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)

plot(predlag12_EVI_ElB, ci= "area", xlab="EVI*10", ylab="RR", main="Lag=12", col=2, lwd = 3)
plot(predlag12_EVI_Tie, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)
plot(predlag12_EVI_PtoL, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)
plot(predlag12_EVI_Zar, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)
plot(predlag12_EVI_Tum, ci= "area", xlab="EVI*10", ylab="RR", main="", col=2, lwd = 3)


Fig4 <- layout(matrix(1:15, ncol=3))
layout.show(Fig4)
par(mar = c(4, 5, 2, 2) + 0.1)

plot(pred23_temp_ElB, ci= "area", xlab="lag", ylab="RR", main="Temperature=23.5", col=2, lwd = 3)
plot(pred23_temp_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred23_temp_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred23_temp_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred23_temp_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

plot(pred25_temp_ElB, ci= "area", xlab="lag", ylab="RR", main="Temperature=25.5", col=2, lwd = 3)
plot(pred25_temp_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred25_temp_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred25_temp_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred25_temp_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

plot(pred27_temp_ElB, ci= "area", xlab="lag", ylab="RR", main="Temperature=27.5", col=2, lwd = 3)
plot(pred27_temp_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred27_temp_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred27_temp_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred27_temp_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


Fig5 <- layout(matrix(1:15, ncol=3))
layout.show(Fig5)
par(mar = c(4, 5, 2, 2) + 0.1)
plot(predlag1_temp_ElB, ci= "area", xlab="°C", ylab="RR", main="Lag=1", col=2, lwd = 3)
plot(predlag1_temp_Tie, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)
plot(predlag1_temp_PtoL, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)
plot(predlag1_temp_Zar, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)
plot(predlag1_temp_Tum, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

plot(predlag4_temp_ElB, ci= "area", xlab="°C", ylab="RR", main="Lag=4", col=2, lwd = 3)
plot(predlag4_temp_Tie, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)
plot(predlag4_temp_PtoL, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)
plot(predlag4_temp_Zar, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)
plot(predlag4_temp_Tum, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)

plot(predlag12_temp_ElB, ci= "area", xlab="°C", ylab="RR", main="Lag=12", col=2, lwd = 3)
plot(predlag12_temp_Tie, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)
plot(predlag12_temp_PtoL, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)
plot(predlag12_temp_Zar, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)
plot(predlag12_temp_Tum, ci= "area", xlab="°C", ylab="RR", main="", col=2, lwd = 3)


Fig6 <- layout(matrix(1:15, ncol=3))
layout.show(Fig6)
par(mar = c(4, 5, 2, 2) + 0.1)

plot(pred20_rain_ElB, ci= "area", xlab="lag", ylab="RR", main="Rain=20 mm", col=2, lwd = 3)
plot(pred20_rain_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred20_rain_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred20_rain_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred20_rain_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

plot(pred100_rain_ElB, ci= "area", xlab="lag", ylab="RR", main="Rain=100 mm", col=2, lwd = 3)
plot(pred100_rain_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred100_rain_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred100_rain_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred100_rain_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)

plot(pred180_rain_ElB, ci= "area", xlab="lag", ylab="RR", main="Rain=180 mm", col=2, lwd = 3)
plot(pred180_rain_Tie, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred180_rain_PtoL, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred180_rain_Zar, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)
plot(pred180_rain_Tum, ci= "area", xlab="lag", ylab="RR", main="", col=2, lwd = 3)


Fig7 <- layout(matrix(1:15, ncol=3))
layout.show(Fig7)
par(mar = c(4, 5, 2, 2) + 0.1)

plot(predlag1_rain_ElB, ci= "area", xlab="mm", ylab="RR", main="Lag=1", col=2, lwd = 3)
plot(predlag1_rain_Tie, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)
plot(predlag1_rain_PtoL, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)
plot(predlag1_rain_Zar, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)
plot(predlag1_rain_Tum, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)

plot(predlag4_rain_ElB, ci= "area", xlab="mm", ylab="RR", main="Lag=4", col=2, lwd = 3)
plot(predlag4_rain_Tie, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)
plot(predlag4_rain_PtoL, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)
plot(predlag4_rain_Zar, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)
plot(predlag4_rain_Tum, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)


plot(predlag12_rain_ElB, ci= "area", xlab="mm", ylab="RR", main="Lag=12", col=2, lwd = 3)
plot(predlag12_rain_Tie, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)
plot(predlag12_rain_PtoL, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)
plot(predlag12_rain_Zar, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)
plot(predlag12_rain_Tum, ci= "area", xlab="mm", ylab="RR", main="", col=2, lwd = 3)






