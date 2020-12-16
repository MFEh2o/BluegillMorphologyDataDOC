#Create log-normally distributed x data
x <- rlnorm(50,0,1)
#See, it's log-normally distributed
hist(x)
#Generate y values according to a known regression relationship between y and x
#This regression relationship is log(y) = log(2) + 1.2*log(x) + e
y <- (2*x^1.2)*rlnorm(50,0,0.05)



#Plot y~x on natural scale
plot(y~x)
#Plot with log axes
plot(y~x,log='xy')


#Fit regression model on log-transformed data
lm1 <- lm(log(y)~log(x))
#Summary of fitted model. Note that fitted parameter estimates match the known parameter estimates
summary(lm1)



#Calculate model-predicted values across a range of x values
#First create the range of x values
newdata <- data.frame(x=seq(0.05,10,0.1))
#Now calculate predicted values
hats <- predict.lm(lm1,newdata)
#Exponentiate the predictions to get them back in natural units (i.e. to go from predictions of log(y) to predictions of y)
expHats <- exp(hats)
#Plot the predicted values
points(expHats~newdata$x,type='l')


setwd("~/")
eyeWidthDat<-read.csv("eyewidthsFINAL.csv")

dfeye<-data.frame(eyewidth<-eyeWidthDat$eyeWidth,
                  fishIDeye<-eyeWidthDat$fishID,
                  lakeID<-eyeWidthDat$lakeID,
                  fishLength<-eyeWidthDat$fishStdLength,
                  DOClevel<-eyeWidthDat$DOClevel,
                  DOC<-eyeWidthDat$DOC,
                  eyewidth.ss<-eyeWidthDat$sizeStandardize,
                  lakeSS<-eyeWidthDat$lakeSS,basin<-eyeWidthDat$basin, 
                  area<-eyeWidthDat$area,maxDepth<-eyeWidthDat$max_Depth, fitted<-eyeWidthDat$fitted)
hist(eyewidth.ss)
plot(eyewidth.ss~DOC)
plot(eyewidth.ss~DOC,log='xy')
lm1 <- lmer(log(eyewidth.ss) ~ 1 + log(DOC) + basin + basin:log(DOC) + (1|lakeID),data = dfeye)
summary(lm1)
hats<-fitted(lm1)
exphats<-exp(hats)
plot(eyewidth.ss~DOC)
plot(eyewidth.ss~DOC,log='xy')
plot(hats~DOC, log='xy')
plot(exphats~DOC,log='xy')
plot(exphats~DOC,log='y')
