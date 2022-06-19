# Lecture 1
# Linear Regression

# Read Data
# kalama = read.table(file = file.choose(), header = T)
# head(kalama, 10)

# Covariance
# cov(kalama$age, kalama$height)

# Correlation
# cor(kalama$age, kalama$height)
# see all of the correlations
# cor(kalama)
# correlation graphs
# plot(kalama)

# Fitting a linear model
# res <- lm(height~age, data=kalama)
# kalama.anova <- anova(res)
# kalama.summary<-summary(res)

# Predict a new observation
# newdata = data.frame(age = 43, anxiety=2.7)
# pred_plim <- predict(linear_model, newdata, interval="predict")
# pred_clim <- predict(linear_model, newdata, interval="confidence")

# Lecture 2: Logistic Regression
# Read/prepare data
# donner <- read.table...
# Remove rows with an NA value
# donner.na<-na.omit(subset(donner,select=c('Age','Outcome','Sex')))
# donner.na$fem = as.numeric(donner.na$Sex=="Female")

# Fit the model
# donner.log<-glm(Outcome ~ Age + fem, data = donner.na, family=binomial(link="logit"))
# summary(donner.log)

# Odds ratios
# exp(donner.log$coefficients)
# Confidence interval
# exp(confint(donner.log))
# exp(cbind(OR = donner.log$coefficients, confint(donner.log)))
# Odds ratio 10 year increase
# exp(donner.log$coefficients*10)
# exp(c(OR = donner.log$coefficients[2]*10,confint(donner.log)[2,]*10))

# Plot logit curve
# logit<-function(x){log((x/1-x))}
# ilogit<- function(x,a,b){exp(a+b*x)/(1+exp(a+b*x))}

# Predicted probs of survival
# newdata2 <- data.frame(fem=1, Age=mean(donner.na$Age))
# newdata2$greP<-predict(donner.log, newdata=newdata2,type="response")

# Interaction model
# m4<-glm(Outcome ~ Age*fem,data=donner.na,family=binomial(link="logit"))

# Comparing models AIC
# donner.list=list()
# donner.list[[1]]=glm(Outcome ~ Age,data=donner.na,family=binomial(link="logit"))
# ... do the same for the other models, (changing the number)
# donner.modnames = c("Age", etc etc.)
# donner.aictab=aictab(cand.set=donner.list,modnames = donner.modnames)

# Lecture 3 Multilevel Models

#Trellis graph
# Who cares!

# lmer stuff
# install.packages("lme4")
# install.packages("arm")
# install.packages("pbkrtest")

# library(lme4)
# library(lattice)
# library(arm)
# library(car)
# library(pbkrtest)

# early.int1$age0<-early.int1$age-1

# Fitting the model
# early.lmer1<-lmer(cog~1+age0*program+(1+age0|id), REML = FALSE, data=early.int1)
# summary(early.lmer1)
# display(early.lmer1)

# Confidence intervals for fixed effects
# Wald, bootstrap, profile likelihood

# confint(early.lmer1,par=5:8,method="Wald",oldNames=FALSE)
# confint(early.lmer1,method="boot",boot.type="perc",oldNames=FALSE,nsim=500)

# Get p values
# confint(early.lmer1, level = 0.95, method="profile", oldNames = FALSE)
# early.lmer1.df.KR <- get_ddf_Lb(early.lmer1, fixed(early.lmer1))
# early.lmer1.coef=coef(summary(early.lmer1))
# early.lmer1.p.KR <- cbind(early.lmer1.coef,2 * (1 - pt(abs(early.lmer1.coef[,3], early.lmer1.df.KR))))
# early.lmer1.p.KR

# likelihood ratio tests between diff models
# early.lmer.noprog(lmer(cog~1+age0+(1 + age|id), REML = FALSE, data=early.int1))
# early.lmer1.intprog<-lmer(cog~1+age0+program+(1 + age0|id), REML = FALSE, data=early.int1)
# anova(early.lmer1.noprog, early.lmer1.intprog,early.lmer1)

# Lecture 4 Multilevel Models Continued

# summary(early.lmer1)
# random effects covariance matrix
# D.early=unclass(VarCorr(early.lmer1))$id
# D.early

# Predict random effects
# early.lmer1.re=anef(early.lmer1)$id
# plot(early.lmer1.re, + main="Random intercept (b0i) versus random slope)

# Creating subject specific intercepts and slopes
# ind.coef=coef(early.lmer1)$id
# head(ind.coef)

# Lecture 4: Cluster data
# ratpup <- read.table(...)
# ratpup$sex1[ratpup$sex == "Female"] <- 1
# ratpup$sex1[ratpup$sex == "Male"] <- 0
# attach(ratpup)

# table
# g <- function(x)c(N=length(x), Mean=mean(x,na.rm=TRUE), + SD=sd(x,na.rm=TRUE), Min=min(x,na.rm=TRUE),Max=max(x,na.rm=TRUE))
# summarize(weight, by=list(treatment,sex),g)

# comparing the distributions of birth weights
# library(lattice)
# library(grid)
# bwplot(weight ~ sex|treatment, data=ratpup, aspect = 2, 
# + ylab="Birth Weights", xlab="SEX",
# + main = "Boxplots of birth weights for levels of treatment by sex")

# Fitting a homoscedastic model
# library(nlme)
# meanfull.hom <- lme(weight~treatment + sex1 + litsize + treatment:sex1,
# random = ~1 | litterid, ratpup, method = "REML")
# lme() treats the lowest level of a factor as the ref category.
# to change use:
# treatment=relevel(treatment,ref="High")
# summary(meanfull.hom)
# anova(meanfull.hom)
# anova is meaningless in REML, comparing models with diff fized effects

# Display random effects
# random.effects(meanfull.hom)

# Fit a heteroscedastic model
# meanfull.het <- lme(weight~treatment + sex1 + litsize + treatment:sex1,
# random = ~1 | litterid, ratpup, method = "REML",
# weights = varIdent(form = ~1 | treatment))
# summary(meanfull.het)

# hetero vs homoscedastic
# anova(meanfull.hom,meanfull.het)

# High=low dose: Equal residual variance
# ratpup$trtgrp[treatment=="Control"] <- 1
# ratpup$trtgrp[treatment == "Low" | treatment == "High"] <- 2

# meanfull.hilo <- lme(weight~treatment + sex1 + litsize + treatment:sex1,
# random = ~1 | litterid, ratpup, method = "REML",
# weights = varIdent(form = ~1 | trtgrp))
# summary(meanfull.hilo)
# anova(meanfull.hilo)

# meanfull.hilo.nolitter <-gls(weight~treatment + sex1 + litsize + treatment:sex1,
# data = ratpup, weights = varIdent(form=~1|trtgrp))

# Lecture 5: Missing Data
# library(mice)
# library(lattice)
# library(VIM)
# library(aod)
# library(BaM)

# Read titanic data
# titanic.missing <- read.table(...)
# head(titanic.missing,10)

# Exploring missingness
# titanic.missing.agg = aggr(titanic.missing,numbers=TRUE,
# + prop=FALSE, ylab=c("Histogram of missing data","Pattern"))
# titanic.missing.agg

# aggr(titanic.missing, combined=TRUE, numbers = TRUE,
# prop = TRUE, cex.numbers=0.87, verheight = FALSE)

# Amount of missingness in age for each survived group
# barMiss(titanic.missing[,c9"survived","age")])

# Amount in each sex
# barMiss(titanic.missing[c,("sex","age")]
# histMiss(titanic.missing))

# Fitting logistic regression for complete cases
# titanic.logistic.omit<-glm(survived~pclass + sex + age, family=binomial,data = titanic,missing)
# summary(titanic.logisitc.omit)

# wald.test(b=coef(titanic.logistic.omit), Sigma=vcov(titanic.logistic.omit).
# + Terms=2:3)

# Odds ratios
# exp(cbind(OR = titanic.logistic.omit$coefficients,
# + confint(titanic.logistic.omit)))

# Multiple imputation
# Study the pattern of missingness
# pattern=md.pattern(titanic.missing)
# pattern

# pairs=md.pairs(titanic.missing)
# pairs

# Imputing missing values
# imp <- mice(titanic.missing, m=100)
# imp

# Checking imputed values
# imp$imp$age[1:10,1:5]

# Can speify the method for mice
# ex: imp <-mice(titanic.missing, meth= c("","","logreg", "pmm"), m=100)

# Analyze imputed data
# fit <-with(data=imp, exp=glm(survived~pclass + sex + age, family=binomial))

# MI.matrix<-matrix(0,100,5)
# for(k in 1:100) MI.matrix[k,]<-coefficients(fit$analyses[[k]])
# MI.results=data.frame(Intercept=MI.matrix[,1],pclass2=MI.matrix[,2].
# MI.results[1:10,] 
#+ pclass3=MI.matrix[,3], sex=MI.matrix[,4], age=MI.matrix[,5])

# combine results using Rubin's rule
# est <- pool(fit)
# summary(est)

# Inverse Probability Weighting (IPW)
# create missing data indicator variable r
# titanic.missing$r<-as.numeric(!is.na(titanic,missing$age))*as.numeric(!is.na(titanic,missing$sex))
# head(titanic.missing, 15)

# Fitting the logistic regression model to calculate the probs of beign complete
# titanic.ipw.glm<-glm(r ~ pclass + survived, data=titanic.missing,family=binomial)
# summary(titanic.ipw.glm)

# Calculating the weights
# titanic.missing$w<-1/fitted(titanic.ipw.glm)

# Final IPW results
# titanic.results.ipw<- glm(survived~pclass + sex + age, data=titanic.missing, weights=titanic.missing$w,
# + family=binomial)
# summary(titanic.results.ipw)


# Part II
# Lecture 1
# I-1 Bootstrap and Cross Validation
# install.packages("ISLR")
# install.packages("boot")
# library("boot")
# library("ISLR")

# set.seed(2)
# dim(Auto)
# train=sample(392,196)
# train
# 
# lm.fit=lm(mpg~horsepower,data=Auto,subset=train)
# 
# test=Auto[-train,]
# predictions=predict(lm.fit,test)
# mean((test$mpg-predictions)^2)
# 
# set.seed(1)
# testMSE=rep(0,6)
# trainMSE=rep(0,6)
# for(i in 1:6){
#   lm.fit=lm(mpg~poly(horsepower,i),data=Auto,subset=train)
#   testMSE[i]=mean((test$mpg-predict(lm.fit,test))^2)
#   trainMSE[i]=mean((Auto$mpg[train]-predict(lm.fit))^2)
# }
# plot(testMSE,type="b",ylim=c(min(testMSE,trainMSE),max(testMSE,trainMSE)))
# lines(trainMSE, type="b",col=2)
# 
# loocv.error=rep(0,6)
# cv10.error=rep(0,6)
# for (i in 1:6){
#   glm.fit=glm(mpg~poly(horsepower,i),data=Auto)
#   loocv.error[i]=cv.glm(Auto,glm.fit)$delta[1]
#   cv10.error[i]=cv.glm(Auto,glm.fit,K=10)$delta[1]
# }
# 
# #Bootstrap
# x <- seq(-4, 4, length=100)
# hx <- dt(x,1.1)
# hx2 <- dnorm(x,sd=1)
# plot(hx2~x,type="l")
# lines(hx~x,col=2)
# set.seed(3)#make sure we get the same results
# x = rnorm(1000)#some values, normal distribution around 0, sd=1
# y = 1 + x + rt(1000,1.1) #simple relation, added noise is broad tailed!
# plot(x,y) #broad tails give outliers
# model = lm(y~x)
# summary(model)#not significant in my case
# 
# bootfun = function(data,b,formula){
#   # b are the bootstrap sample indices
#   d = data[b,]
#   return(lm(d[,1]~d[,2], data = d)$coef[2]) # thats the beta1 coefficient
# }
# 
# result=boot(data=data.frame(y,x), statistic=bootfun, R=1000)
# #R is number of boostrap samples
# result
# plot(result)
# boot.ci(result, index=1, type=c("bca"))#estimated confidence interval for bootstrap
# confint(model)#confidence interval for linear model