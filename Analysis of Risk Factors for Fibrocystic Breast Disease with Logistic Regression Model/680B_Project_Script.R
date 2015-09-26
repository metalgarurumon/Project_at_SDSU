##Loading required packages
library(xlsx)
library(survival)

#Loading Dataset

x<-read.xlsx(file="F:/Academics/2015 Spring/STAT680B/Final Project/bbdm13.xls",sheetName="BBDM13")

names(x)
attach(x)

####Descriptive Statistics
#Age at interview
mean(AGMT)
median(AGMT)
sd(AGMT)
#Matched cannot use
#Final diagnosis
table(FNDX)
#Degree
table(DEG)
#No missing
#Regular Medical Check ups
table(CHK)
#No missing values
#Age at first pregnancy
#Has missing values do not use.
#Age at Menarche
mean(AGMN)
median(AGMN)
sd(AGMN)
#Number of stillbirths or carriages
table(NLV)
##Either numerical or categorical
#Has missing
#Number of live births
table(LIV)
#Has missing values
#Either numerical or categorical
#Weight
mean(WT)
median(WT)
sd(WT)
#Age at last menstrual cycle
mean(AGLP)
median(AGLP)
sd(AGLP)
#Marital Status
table(MST)
#No missing
#Highest Grade in School
mean(HIGD)
#Contains 7 good variables
#Just enought
#Use HIGD DEG CHK AGMN WT AGLP MST


#####Real Descriptive Statistics
xC<-subset(x,FNDX==1)
xCont<-subset(x,FNDX==0)


##Cases
#Age
mean(xC$AGMT)
sd(xC$AGMT)
#High School High school age
min(xC$HIGD)
max(xC$HIGD)
mean(xC$HIGD)
sd(xC$HIGD)
#Degree completed
min(xC$DEG)
max(xC$DEG)
mean(xC$DEG)
sd(xC$DEG)
table(xC$DEG)
#Regular Medical Check up
table(xC$CHK)
#Age at Menarche
min(xC$AGMN)
max(xC$AGMN)
mean(xC$AGMN)
sd(xC$AGMN)
#Weight
min(xC$WT)
max(xC$WT)
mean(xC$WT)
sd(xC$WT)
#Age at Last Menstrual Period
min(xC$AGLP)
max(xC$AGLP)   
mean(xC$AGLP)
sd(xC$AGLP)
#Marital Status
table(xC$MST)

##Controls
#Age
mean(xCont$AGMT)
sd(xCont$AGMT)
#High School High school age
min(xCont$HIGD)
max(xCont$HIGD)
mean(xCont$HIGD)
sd(xCont$HIGD)
#Degree completed
min(xCont$DEG)
max(xCont$DEG)
mean(xCont$DEG)
sd(xCont$DEG)
table(xCont$DEG)
#Regular Medical Check up
table(xCont$CHK)
#Age at Menarche
min(xCont$AGMN)
max(xCont$AGMN)
mean(xCont$AGMN)
sd(xCont$AGMN)
#Weight
min(xCont$WT)
max(xCont$WT)
mean(xCont$WT)
sd(xCont$WT)
#Age at Last Menstrual Period
min(xCont$AGLP)
max(xCont$AGLP)
mean(xCont$AGLP)
sd(xCont$AGLP)
#Marital Status
table(xCont$MST)


#check if degree needs to be categorical or continuous
z1<-clogit(FNDX~1+DEG+strata(STR),data=x,method="exact")
z2<-clogit(FNDX~1+factor(DEG)+strata(STR),data=x,method="exact")
LRT=-2*(logLik(z1)-logLik(z2))
1-pchisq(LRT,4)

z3<-clogit(FNDX~1+HIGD+strata(STR),data=x,method="exact")
z4<-clogit(FNDX~1+factor(HIGD)+strata(STR),data=x,method="exact")
LRT=-2*(logLik(z3)-logLik(z4))
1-pchisq(LRT,14)
#Keep continuous

#check if HIGD needs to be categorical or continuous
z3<-clogit(FNDX~1+HIGD+strata(STR),data=x,method="exact")
z4<-clogit(FNDX~1+factor(HIGD)+strata(STR),data=x,method="exact")
LRT=-2*(logLik(z3)-logLik(z4))
1-pchisq(LRT,15)
#Keep continuous

###Model Destroying
#Eliminate at a 0.05 level
y<-clogit(FNDX~HIGD+factor(CHK)+AGMN+WT+AGLP+factor(MST)+DEG+strata(STR),method="exact",data=x)
y
logLik(y)
#Eliminate Divorced
x$sep<-ifelse(x$MST==3,1,0)
x$Wid<-ifelse(x$MST==4,1,0)
x$Nev<-ifelse(x$MST==5,1,0)
y1<-clogit(FNDX~HIGD+factor(CHK)+AGMN+WT+AGLP+factor(sep)+factor(Wid)+factor(Nev)+DEG+strata(STR),method="exact",data=x)
y1
#Eliminate Seperated
y2<-clogit(FNDX~HIGD+factor(CHK)+AGMN+WT+AGLP+factor(Wid)+factor(Nev)+DEG+strata(STR),method="exact",data=x)
y2
#Eliminate Widowed
y3<-clogit(FNDX~HIGD+factor(CHK)+AGMN+WT+AGLP+factor(Nev)+DEG+strata(STR),method="exact",data=x)
y3
#Highest Grade Completed
y4<-clogit(FNDX~factor(CHK)+AGMN+WT+AGLP+factor(Nev)+DEG+strata(STR),method="exact",data=x)
y4
#Eliminate Degree
y5<-clogit(FNDX~factor(CHK)+AGMN+WT+AGLP+factor(Nev)+strata(STR),method="exact",data=x)
y5
#Eliminate Age at last menstrual period
y6<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x)
y6
AICc=function(object){
  n=length(object$y)
  r=length(object$coef)
  AICc=AIC(object)+2*r*(r+1)/(n-r-1)
  list(AIC=AIC(object), AICc=AICc, BIC=BIC(object))
}
AIC(y1,y2,y3,y4,y5,y6)
AICc(y1)
AICc(y2)
AICc(y3)
AICc(y4)
AICc(y5)
AICc(y6)
BIC(y1,y2,y3,y4,y5,y6)
#Final model is y6
summary(y6)
names(y6)
y6$coefficients
exp(cbind(OR = coef(y6), confint(y6)))

#LR Test

LRT=-2*(logLik(y6)-logLik(y1))
1-pchisq(LRT,5)

##Loading excel file from SAS

x_d<-read.csv(file="F:/Academics/2015 Spring/STAT680B/Final Project/diag_680.csv")
names(x_d)

ZC<-subset(x_d,FNDX==1)
ZCont<-subset(x_d,FNDX==0)
#Leverages
plot(x_d$leve~x_d$fitt,xlab="Predicted Values",ylab="Leverages", main="Leverages vs. Predicted Values")
plot(ZC$leve~ZC$fitt,xlab="Predicted Values",ylab="Leverages")
points(ZCont$fitt,ZCont$leve,type="p",pch=4)
legend(x="topright", pch=c(1,4), legend=c("Case", "Control"))
#0.13
#betas
plot(x_d$betas~x_d$fitt,xlab="Predicted Values",ylab="Cook's Distance", main="Cook's Distance vs. Predicted Values")
plot(ZC$betas~ZC$fitt,xlab="Predicted Values",ylab="Leverages")
points(ZCont$fitt,ZCont$betas,type="p",pch=4)
legend(x="topright", pch=c(1,4), legend=c("Case", "Control"))
#0.18
#Pearson Chisquare
plot(x_d$chisqu~x_d$fitt,xlab="Predicted Values",ylab="Change in Pearson x2", main="Change in Pearson x2 vs. Predicted Values")
plot(ZC$chisqu~ZC$fitt,xlab="Predicted Values",ylab="Leverages")
points(ZCont$fitt,ZCont$chisqu,type="p",pch=4)
legend(x="topright", pch=c(1,4), legend=c("Case", "Control"))
#Use 13

#Creating indicator variables
x_d$leverage=ifelse(x_d$leve>0.13,1,0)
sum(x_d$leverage)
#Leverages says 5 observations to worry about
x_d$cbetas=ifelse(x_d$betas>0.18,1,0)
sum(x_d$cbetas)
#betas says to look at 5 observations
x_d$cchi=ifelse(x_d$chisqu>13,1,0)
sum(x_d$cchi)
#3 observations

#Stratas to review 10
#39,31,26,24,18,17,14,12,10,4
#39
x_39<-subset(x,STR!=39)
y_39<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_39)
y_39
y6

#31
x_31<-subset(x,STR!=31)
y_31<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_31)
y_31
y6

#26
x_26<-subset(x,STR!=26)
y_26<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_26)
y_26
y6

#24
x_24<-subset(x,STR!=24)
y_24<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_24)
y_24
y6

#18
x_18<-subset(x,STR!=18)
y_18<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_18)
y_18
y6

#17
x_17<-subset(x,STR!=17)
y_17<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_17)
y_17
y6

#14
x_14<-subset(x,STR!=14)
y_14<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_14)
y_14
y6

#12
x_12<-subset(x,STR!=12)
y_12<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_12)
y_12
y6

#10
x_10<-subset(x,STR!=10)
y_10<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_10)
y_10
y6

#4
x_4<-subset(x,STR!=4)
y_4<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x_4)
y_4
y6


#Not much change in the model, I recommend keeping all the stratum


#testing interactions
y6_int<-clogit(FNDX~factor(CHK)+AGMN+WT+factor(Nev)+strata(STR),method="exact",data=x)
#tested all possible combinations
#No interactions


####################
# Diagnostic plots #
####################
# User written R function: examine.logistic.reg(), which takes a glm object as required argument

one.fourth.root=function(x){
  x^0.25
}

# use source function or copy and paste the function examine.logistic.reg
source(file="http://statistics.unl.edu/faculty/bilder/categorical/Chapter5/Examine.logistic.reg.R")
save1=examine.logistic.reg(y6, identify.points=T, scale.n=one.fourth.root, scale.cookd=sqrt)








