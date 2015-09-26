
# Procedure 1: generate grouped data set and summary of the data
library(nlme)
raw_NCTD <- read.table("F:/Academics/2014 Fall/Linear Mixed Models/Final Project/NCTD.txt", head=T)
raw_NCTD$drug.f <- factor(raw_NCTD$drug)
NCTD <- groupedData(cell ~ time | drug.f, data=raw_NCTD)
NCTD
summary(NCTD)
plot(NCTD, main='Growth curve of cancer cells')


# Procedure 2: lme model selection

# Top-down Strategy
# Step1: random effects

# We drop the random intercepts because of the experimental design
model.lme.2 <- lme(cell~time+I(time^2)+drug+I(drug^2), random=~time+I(time^2)-1, data=NCTD)

# H0: I(time^2) can be removed from random effects
model.lme.3 <- update(model.lme.2, random=~time-1)
anova(model.lme.2, model.lme.3)
# Conclusion: I(time^2) cannot be removed random effects
# model.lme.2 is the winner

# H0: Heterogeneous residual variances among groups
model.lme.2.hete <- update(model.lme.2, weights=varIdent(form=~1|drug.f) )
anova(model.lme.2, model.lme.2.hete)
# Conclusion: There are heterogeneous residual variances among groups
# model.lme.2.hete is the winner


# Step2: fixed effects

# The growth curve seems like parabola, so we keep I(time^2) in fixed effects

# H0: I(drug^2) can be removed from fixed effects
model.lme.5.hete.ml <- update(model.lme.2.hete.ml, fixed=cell~time+I(time^2)+drug)
anova(model.lme.2.hete.ml, model.lme.5.hete.ml)
# Conclusion: I(drug^2) cannot be removed from fixed effects
# model.lme.2.hete.ml is the winner

# top-down lme winner: model.lme.2.hete


# Step-up Strategy
model.lme.6 <- lme(cell~1, random=~time+I(time^2)-1, data=NCTD)

# Step1: random effects
# H0: I(time^2) can be removed from random effects
model.lme.7 <- update(model.lme.6, random=~time-1)
anova(model.lme.6, model.lme.7)
# Conclusion: I(time^2) cannot be removed random effects
# model.lme.6 is the winner

# H0: Heterogeneous residual variances among groups
model.lme.6.hete <- update(model.lme.6, weights=varIdent(form=~1|drug.f) )
anova(model.lme.6, model.lme.6.hete)
# Conclusion: There are heterogeneous residual variances among groups
# model.lme.6.hete is the winner

# Step2: fixed effects
# H0: I(time^2) dose not need to be added to fixed effects
model.lme.6.hete.ml <- update(model.lme.6.hete, fixed=cell~time, method='ML')
model.lme.7.hete.ml <- update(model.lme.6.hete.ml, fixed=cell~time+I(time^2))
anova(model.lme.6.hete.ml, model.lme.7.hete.ml)
# Conclusion: I(time^2) should be added to fixed effects
# model.lme.7.hete.ml is the winner

# H0: I(drug^2) dose not need to be added to fixed effects
model.lme.8.hete.ml <- update(model.lme.7.hete.ml, fixed=cell~time+I(time^2)+drug)
model.lme.9.hete.ml <- update(model.lme.7.hete.ml, fixed=cell~time+I(time^2)+drug+I(drug^2))
anova(model.lme.8.hete.ml, model.lme.9.hete.ml)
# Conclusion: I(drug^2) should be added to fixed effects
# model.lme.9.hete.ml is the winner

# step-up lme winner: model.lme.9.hete
model.lme.9.hete <- update(model.lme.9.hete.ml, method='REML')

# model.lme.8.hete are the same with model.lme.2.hete
anova(model.lme.2.hete, model.lme.9.hete)

# model.lme.2.hete is the final lme winner
anova(model.lme.2.hete)
summary(model.lme.2.hete)
plot(model.lme.2.hete, resid(., type='p')~fitted(.) | drug.f, abline=0, lty=2)
qqnorm(resid(model.lme.2.hete))
qqline(resid(model.lme.2.hete))
plot(density(resid(model.lme.2.hete)))
plot(augPred(model.lme.2.hete))


# model.lme.5.hete is an alternatve model
model.lme.5.hete <- update(model.lme.5.hete.ml, method='REML')
summary(model.lme.5.hete)
plot(model.lme.5.hete, resid(., type='p')~fitted(.) | drug.f, abline=0, lty=2)
qqnorm(resid(model.lme.5.hete))
qqline(resid(model.lme.5.hete))
hist(resid(model.lme.5.hete))
plot(density(resid(model.lme.5.hete)))
plot(augPred(model.lme.5.hete))

# The diagnostics plot and the fitted value plot are similar for the two models
# we finally choose model.lme.2.hete based on the background 
# that the drug concentration and its cytotoxicity is usually not a simple linear relationship


# Procedure 3: nlme model selection

# model: cell~N*(a+b*drug)^time #
model.nlme.1 <- nlme(cell~N*(a+b*drug)^time, data=NCTD,
                     fixed=N+a+b~1, random=N+a+b~1, start=c(N=200, a=1.05, b=-0.0001) )

# random effects
# H0: N can be removed from random effects
model.nlme.2 <- update(model.nlme.1, random=a+b~1)
anova(model.nlme.1, model.nlme.2)
# Conclusion: N cannot be removed from random effects
# model.nlme.1 is the winner #

# H0: a can be removed from random effects
model.nlme.3 <- update(model.nlme.1, random=N+b~1)
anova(model.nlme.1, model.nlme.3)
# Conclusion: a cannot be removed from random effects
# model.nlme.1 is the winner #

# H0: b can be removed from random effects
model.nlme.4 <- update(model.nlme.1, random=N+a~1)
anova(model.nlme.1, model.nlme.4)
# Conclusion: b should be removed from random effects
# model.nlme.4 is the winner #


# H0: Heterogeneous residual variances among groups
model.nlme.4.hete <- update(model.nlme.4, weights=varIdent(form=~1|drug.f) )
anova(model.nlme.4, model.nlme.4.hete)
# Conclusion: There are heterogeneous residual variances among groups
# model.nlme.4.hete is the final winner #


summary(model.nlme.4.hete)
plot(model.nlme.4.hete, resid(., type='p')~fitted(.) | drug.f, abline=0, lty=2)
qqnorm(resid(model.nlme.4.hete))
qqline(resid(model.nlme.4.hete))
plot(density(resid(model.nlme.4.hete)))
plot(augPred(model.nlme.4.hete))



# Procedure 4: final game between model.lme.2.hete and model.nlme.1.hete2
# compare info. criteria between LMM and NLMM #
AIC(model.lme.2.hete, model.nlme.4.hete)
BIC(model.lme.2.hete, model.nlme.4.hete)

# LMM model.lme.2.hete is the final winner


