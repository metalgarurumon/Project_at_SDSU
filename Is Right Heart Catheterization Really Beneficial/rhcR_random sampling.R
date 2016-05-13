setwd("F:/Academics/2016 Spring/Survival Analysis/Final Project")
rhc <- read.csv("rhc.csv")
summary(rhc)
attach(rhc)
library(survival)

#################################
# data cleaning
#################################

colSums( is.na(rhc) )
#In dthdte, replace NA with value from lstctdte to find days

rhc$dthdte1<-ifelse(is.na(rhc$dthdte), rhc$lstctdte, rhc$dthdte)
rhc$time<-rhc$dthdte1-rhc$sadmdte

rhc$cat12=ifelse(cat1=="CHF",1,0)
rhc$cat13=ifelse(cat1=="Cirrhosis",1,0)
rhc$cat14=ifelse(cat1=="Colon Cancer",1,0)
rhc$cat15=ifelse(cat1=="Coma",1,0)
rhc$cat16=ifelse(cat1=="COPD",1,0)
rhc$cat17=ifelse(cat1=="Lung Cancer",1,0)
rhc$cat18=ifelse(cat1=="MOSF w/Malignancy",1,0)
rhc$cat19=ifelse(cat1=="MOSF w/Sepsis",1,0)

rhc$caN=ifelse(ca=="No",1,0)
rhc$caY=ifelse(ca=="Yes",1,0)


rhc$ninsclas2=ifelse(ninsclas=="Medicare",1,0)
rhc$ninsclas3=ifelse(ninsclas=="Medicare & Medicaid",1,0)
rhc$ninsclas4=ifelse(ninsclas=="No insurance",1,0)
rhc$ninsclas5=ifelse(ninsclas=="Private",1,0)
rhc$ninsclas6=ifelse(ninsclas=="Private & Medicare",1,0)

rhc$black=ifelse(race=="black",1,0)
rhc$white=ifelse(race=="white",1,0)

rhc$income2=ifelse(income=="$25-$50k",1,0)
rhc$income3=ifelse(income=="> $50k",1,0)
rhc$income4=ifelse(income=="Under $11k",1,0)

rhc$death1<-ifelse(death=="Yes", 1, 0)

#mean(adld3p,na.rm=TRUE)
#rhc$adld3p1=ifelse(is.na(adld3p),1.182071,adld3p)
#mean(urin1,na.rm=TRUE)
#rhc$urin11=ifelse(is.na(urin1),2192.454,urin1)

quantile(rhc$aps1, 0.5)
rhc$aps1_grp <-ifelse(rhc$aps1<54, 1, 2)
table(rhc$aps1_grp)


# Randomly sample 3823 data as training set,
# the remaining 1912 data are used as validation set.

set.seed(1)
train <- sample(1:5735,3823)
rhc.training <- rhc[train,]
dim(rhc.training)
validation <- setdiff(1:5735,train)
rhc.validation <- rhc[validation,]
dim(rhc.validation)


##########################
# Unajusted survival curve
##########################

fitKM <- survfit( Surv(time,death1)~swang1, type='kaplan', 
                conf.type='log-log', data=rhc.training)
summary(fitKM)
plot(fitKM, xlab='Days', ylab='survival probability', col=c(2,4),lty=c(2,4),
     main="Undjusted Survival Estimation")
legend("topright",c("No RHC","RHC"), col=c(2,4),lty=c(2,4))

fit.rhc <- coxph(Surv(time,death1)~swang1, 
              method='breslow', data=rhc)
summary(fit.rhc)


######################
# Model Selection
######################


fit1 <- coxph(Surv(time,death1)~cat12+cat13+cat14+cat15+cat16+cat17+cat18+cat19
              +caN+caY+cardiohx+chfhx+dementhx
              +psychhx+chrpulhx+renalhx+liverhx+gibledhx+malighx
              +immunhx+transhx+amihx+age+sex+edu     
              +das2d3pc+aps1+scoma1  
              +meanbp1+wblc1+hrt1+resp1+temp1+pafi1 
              +alb1+hema1+bili1+crea1+sod1+pot1   
              +paco21+ph1+swang1+wtkilo1+dnr1
              +resp+card+neuro+gastr+renal+meta+hema
              +seps+trauma+ortho+ninsclas2+ninsclas3
              +ninsclas4+ninsclas5+ninsclas6+black+white
              +income2+income3+income4, 
              method='breslow', data=rhc.training)
summary(fit1)


step(fit1)  # Stepwise AIC
step(fit1, k = log(3823))  # Stepwise BIC


# backward BIC
fit1a <- coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                 caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                 swang1 + dnr1 + hema, data = rhc.training, method = "breslow")
summary(fit1a)


#################################
# Transformation
#################################

#age
plot(rhc.training$age, resid(fit1a), xlab='age', ylab='Residual')  
lines(lowess(rhc.training$age, resid(fit1a), f=.5), col="red") 

fit1a.age <- coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                     caN + liverhx + age +I(age^2) + das2d3pc + aps1 + scoma1 + bili1 + 
                     swang1 + dnr1 + hema, data = rhc.training, method = "breslow")
summary(fit1a.age)
# Reduced Model is preferred.

#das2d3pc
plot(rhc.training$das2d3pc, resid(fit1a), xlab='das2d3pc', ylab='Residual')  
lines(lowess(rhc.training$das2d3pc, resid(fit1a), f=.5), col="red") 

fit1a.das2d3pc <- coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                          caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                          swang1 + dnr1 + hema, data = rhc.training, method = "breslow")
summary(fit1a.das2d3pc)
anova(fit1a.das2d3pc, fit1a)
AIC(fit1a.das2d3pc, fit1a)
BIC(fit1a.das2d3pc, fit1a)
# Reduced Model is preferred.

#aps1
plot(rhc.training$aps1, resid(fit1a), xlab='aps1', ylab='Residual')  
lines(lowess(rhc.training$aps1, resid(fit1a), f=.5), col="red") 

fit1a.aps1<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                    caN + liverhx + age + das2d3pc + aps1 +I(aps1^2) + scoma1 + bili1 + 
                    swang1 + dnr1 + hema, data = rhc.training, method = "breslow")
summary(fit1a.aps1)
# Reduced Model is preferred.

#scoma1
plot(rhc.training$scoma1, resid(fit1a), xlab='scoma1', ylab='Residual')  
lines(lowess(rhc.training$scoma1, resid(fit1a), f=.5), col="red") 

fit1a.scoma1<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                      caN + liverhx + age + das2d3pc + aps1 + scoma1 +I(scoma1^2) + bili1 + 
                      swang1 + dnr1 + hema, data = rhc.training, method = "breslow")
summary(fit1a.scoma1)
anova(fit1a.scoma1, fit1a)
AIC(fit1a.scoma1, fit1a)
BIC(fit1a.scoma1, fit1a)
# Reduced Model is preferred.

#bili1
plot(rhc.training$bili1, resid(fit1a), xlab='bili1', ylab='Residual')  
lines(lowess(rhc.training$bili1, resid(fit1a), f=.5), col="red") 

fit1a.bili1<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                     caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 +I(bili1^2)+ 
                     swang1 + dnr1 + hema, data = rhc.training, method = "breslow")
summary(fit1a.bili1)
anova(fit1a.bili1, fit1a)
AIC(fit1a.bili1, fit1a)
BIC(fit1a.bili1, fit1a)
# Reduced Model is preferred.

# Conclusion:  Model with linear forms is preferred.


#############################
# Interaction with RHC
#############################

# Interaction between cat1 and RHC
fit.with.cat1<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                       caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                       swang1 + dnr1 + hema + cat15*swang1 + cat17*swang1 + cat18*swang1, 
                     data = rhc.training, method = "breslow")
summary(fit.with.cat1)
#No interaction

# Interaction between caN and RHC
fit.with.caN<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                      caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                      swang1 + dnr1 + hema + caN*swang1, 
                    data = rhc.training, method = "breslow")
summary(fit.with.caN)
anova(fit.with.caN, fit1a.das2d3pc)
AIC(fit.with.caN, fit1a.das2d3pc)
BIC(fit.with.caN, fit1a.das2d3pc)
#No interaction

# Interaction between age and RHC
fit.with.age<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                      caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                      swang1 + dnr1 + hema + age*swang1, 
                    data = rhc.training, method = "breslow")
summary(fit.with.age)
#No interaction

# Interaction between liverhx and RHC
fit.with.liverhx<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                      caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                      swang1 + dnr1 + hema + liverhx*swang1, 
                    data = rhc.training, method = "breslow")
summary(fit.with.liverhx)
#No interaction

# Interaction between das2d3pc and RHC
fit.with.das2d3pc<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                           caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                           swang1 + dnr1 + hema + das2d3pc*swang1, 
                         data = rhc.training, method = "breslow")
summary(fit.with.das2d3pc)
#No interaction

# Interaction between aps1 and RHC
fit.with.aps1<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                       caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                       swang1 + dnr1 + hema + aps1*swang1, 
                     data = rhc.training, method = "breslow")
summary(fit.with.aps1)
#No interaction

# Interaction between scoma1 and RHC
fit.with.scoma1<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                         caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                         swang1 + dnr1 + hema + scoma1*swang1, 
                       data = rhc.training, method = "breslow")
summary(fit.with.scoma1)
#No interaction

# Interaction between bili1 and RHC
fit.with.bili1<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                        caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                        swang1 + dnr1 + hema + bili1*swang1, 
                      data = rhc.training, method = "breslow")
summary(fit.with.bili1)
#No interaction

# Interaction between hema and RHC
fit.with.hema<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                       caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                       swang1 + dnr1 + hema + hema*swang1, 
                     data = rhc.training, method = "breslow")
summary(fit.with.hema)
#No interaction

#No interaction terms are included in the model.
#Model with linear forms and no interaction terms is the best model after model building.


############################
# Cox PH Assumption
############################

## using schoenfeld residuals

temp<-cox.zph(fit1a)
temp

## log-log survival curves
par(mfrow=c(1,1))

# dnr1
fit_dnr1<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                  caN + liverhx + age + das2d3pc + aps1 + scoma1 + bili1 + 
                  swang1 + strata(dnr1) + hema, data = rhc.training, method = "breslow")
fit_dnr11<-survfit(fit_dnr1)
plot(fit_dnr11,fun='cloglog',lty=1:2,col=1:2)
legend("topleft",c("No DNR", "DNR"),lty=1:2,col=1:2, cex=0.8)

temp<-cox.zph(fit_dnr1)
temp
# do not stratify on dnr1


# aps1
fit_aps1<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                  caN + liverhx + age + das2d3pc + strata(aps1_grp) + scoma1 + bili1 + 
                  swang1 + dnr1 + hema, data = rhc.training, method = "breslow")
fit_aps11<-survfit(fit_aps1)
plot(fit_aps11,fun='cloglog',lty=1:2,col=1:2,
     main="Log-Log survival curve", xlab="Days", ylab="Log-Log survival")
legend("topleft",c("Low APACHE score", "High APACHE score"),lty=1:2,col=1:2)

temp<-cox.zph(fit_aps1)
temp


# use stratified model (aps1 as strata)


############################
# Outliers
############################

fit_training <- coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                        caN + liverhx + age + das2d3pc + strata(aps1_grp) + scoma1 + bili1 + 
                        swang1 + dnr1 + hema, data = rhc.training, method = "breslow")
summary(fit_training)

dresid<-resid(fit_training,'dev')
length(dresid)

par(mfrow=c(1,1))
plot(dresid, ylim=c(-3.5, 3.5), main="Outliers",
     xlab="observations", ylab="Deviance Residuals")

head(sort(dresid))
tail(sort(dresid))


##################################
##Influential Points
##################################

sresid<-resid(fit_training,'score')
dim(sresid) 

apply(sresid, 2, which.max)

data.frame(rhc$cat15, rhc$cat17, rhc$cat18, rhc$caN, rhc$liverhx, 
           rhc$age, rhc$das2d3pc, rhc$scoma1, rhc$bili1, rhc$swang1, 
           rhc$dnr1, rhc$hema)[c(5653, 4624, 1317),]
dresid[1317]

##################################
## Model Validation
##################################


#  Comparing coef. estimates

fit_validation<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                        caN + liverhx + age + das2d3pc + strata(aps1_grp) + scoma1 + bili1 + 
                        swang1 + dnr1 + hema, data = rhc.validation, method = "breslow")
summary(fit_validation)

fit_all<-coxph(Surv(time, death1) ~ cat15 + cat17 + cat18 + 
                        caN + liverhx + age + das2d3pc + strata(aps1_grp) + scoma1 + bili1 + 
                        swang1 + dnr1 + hema, data = rhc, method = "breslow")
summary(fit_all)

tablecomp<-round( cbind(fit_training$coef, confint(fit_training),
                        fit_validation$coef, confint(fit_validation) ), 2)
tablecomp


# Comparing observed vs. expected survival

rs <- cbind(rhc$cat15, rhc$cat17, rhc$cat18, rhc$caN, rhc$liverhx, 
            rhc$age, rhc$das2d3pc, rhc$scoma1, rhc$bili1, rhc$swang1, 
            rhc$dnr1, rhc$hema)[validation,] %*% coef(fit_all)
dim(rs)


## stratum: Low APACHE score 

rs_lowaps <- rs[rhc.validation$aps1_grp==1]

quantile(rs_lowaps, c(.33,.67))               # 33% and 67% quantiles of risk scores

riskgrp_lowaps<-ifelse(rs_lowaps<=1.217248, 1, 2)
riskgrp_lowaps<-ifelse(rs_lowaps>=1.569997, 3, riskgrp_lowaps)
table(riskgrp_lowaps)

rhc.validation.lowaps <- rhc.validation[rhc.validation$aps1_grp==1,]
dim(rhc.validation.lowaps)

# observed

fit_lowaps <- survfit(Surv(time,death1)~riskgrp_lowaps, data=rhc.validation.lowaps)


# expected
zvalues <- data.frame(rhc$cat15, rhc$cat17, rhc$cat18, rhc$caN, rhc$liverhx, 
                      rhc$age, rhc$das2d3pc, rhc$scoma1, rhc$bili1, rhc$swang1, 
                      rhc$dnr1, rhc$hema)[validation,]
dim(zvalues)

zvalues_lowaps <- zvalues[rhc.validation$aps1_grp==1,]
dim(zvalues_lowaps)
names(zvalues_lowaps)

detach()
library(survival)
attach(rhc.validation.lowaps)
expect_lowaps <- survexp(~riskgrp_lowaps+ratetable(cat15=cat15, cat17=cat17, cat18=cat18,
                                            caN=caN, liverhx=liverhx, age=age,
                                            das2d3pc=das2d3pc, aps1_grp=1,
                                            scoma1=scoma1, bili1=bili1, swang1=swang1, 
                                            dnr1=dnr1, hema=hema),
                  data=zvalues_lowaps, ratetable=fit_training)

par(mfrow=c(1,1))
plot(fit_lowaps,xlab="Days", ylab="Survival", main="Model Validation \nStratum: Low APACHE score", 
     col=2:4, lty=1)
lines(expect_lowaps,col=2:4,lty=2)
legend("topright",paste( rep(c("observed", "expected"), 3),
                         rep(c("Low risk", "Medium risk", "High risk"), each=2) ),
       col=rep(2:4, each=2), lty=rep(1:2,3), cex=1 )



## stratum: High APACHE score 

rs_highaps <- rs[rhc.validation$aps1_grp==2]

quantile(rs_highaps, c(.33,.67))               # 33% and 67% quantiles of risk scores

riskgrp_highaps<-ifelse(rs_highaps<=1.359365, 1, 2)
riskgrp_highaps<-ifelse(rs_highaps>=1.812699, 3, riskgrp_highaps)
table(riskgrp_highaps)

rhc.validation.highaps <- rhc.validation[rhc.validation$aps1_grp==2,]
dim(rhc.validation.highaps)

# observed

fit_highaps<-survfit(Surv(time,death1)~riskgrp_highaps, data=rhc.validation.highaps)


# expected

zvalues_highaps <- zvalues[rhc.validation$aps1_grp==2,]
dim(zvalues_highaps)
names(zvalues_highaps)

detach()
library(survival)
attach(rhc.validation.highaps)
expect_highaps <- survexp(~riskgrp_highaps+ratetable(cat15=cat15, cat17=cat17, cat18=cat18,
                                             caN=caN, liverhx=liverhx, age=age,
                                             das2d3pc=das2d3pc, aps1_grp=2,
                                             scoma1=scoma1, bili1=bili1, swang1=swang1, 
                                             dnr1=dnr1, hema=hema),
                  data=zvalues_highaps, ratetable=fit_training)

plot(fit_highaps,xlab="Days", ylab="Survival", main="Model Validation \nStratum: High APACHE score", 
     col=2:4, lty=1)
lines(expect_highaps,col=2:4,lty=2)
legend("topright",paste( rep(c("observed", "expected"), 3),
                         rep(c("Low risk", "Medium risk", "High risk"), each=2) ),
       col=rep(2:4, each=2), lty=rep(1:2,3), cex=1)


detach()

