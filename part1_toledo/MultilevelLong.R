## Multilevel Models: Longitudinal

# Home

setwd("C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Bioinformatics\\Multilevel-Models\\Data-R-code")



# Packages

install.packages("foreign")
install.packages("UsingR")
install.packages("gplots")
install.packages("xtable")
install.packages("Rcmdr")
install.packages("stats")
install.packages("asbio")
install.packages("mvtnorm")
install.packages("multcomp")
install.packages("lawstat")
install.packages("gmodels")
install.packages("doBy")
install.packages("faraway")



install.packages("corrplot") 

install.packages("lme4") 
install.packages("arm")
install.packages("pbkrtest")
install.packages("LMERConvenienceFunctions")
install.packages("languageR")
install.packages("CorrMixed")
install.packages("lsmeans")

## Practical Regression and Anova using R by Julian J. Faraway

# install.packages("faraway")

library('corrplot')

library(CorrMixed)


library(doBy)
library(foreign)
library(UsingR)
library(gplots)
library(xtable)
library(stats)
library("asbio")
library(graphics)
library(mvtnorm)
library(matrixcalc)
library(multcomp)
library(lawstat)
library(MASS)
library(faraway) # library of the book of linearodels in r. You have this book.
library(Rcmdr) # This one opens a windows

library(lme4)
library(lattice)
library("lsmeans")
library(arm)
library(car)
library(pbkrtest)
library(LMERConvenienceFunctions)
library(coda)
library(languageR)

##############################################################
#           General function to plot error bars              #
##############################################################

errbar=function(x,y,height,width,lty=1,col="black")
{arrows(x,y,x,y+height,angle=90,length=width,lty=lty,
col=col)
arrows(x,y,x,y-height,angle=90,length=width,lty=lty,
col=col)}

#############################################################################
# General to create an Spaghetti Plot with the mean o median function on it #
#############################################################################

## tto is the treatment variable or any other variable defining the  groups
## if bytto is true then the means or medians are calculated for  every group
## colme is a  vector given the colors for every mean/median curve

Spaghetti.Plot.new=
function (Dataset, Outcome, Time, Id, tto, Add.Profiles = TRUE, Add.Mean = TRUE, 
    Add.Median = FALSE, Col = 8, colme, Lwd.Me = 3, xlim, ylim, bytto=FALSE,...) 
{
    Object <- x <- Dataset
    Outcome <- x[, paste(substitute(Outcome))]
    Time <- x[, paste(substitute(Time))]
    Id <- x[, paste(substitute(Id))]
    tto<-x[, paste(substitute(tto))]
    Data <- data.frame(cbind(Outcome, Time, Id,tto))
    ttovalue=unique(tto)
    max_val_y <- max(Data$Outcome, na.rm = TRUE)
    max_val_x <- max(Data$Time, na.rm = TRUE)
    if (missing(xlim) == TRUE) {
        xlim = c(0, max_val_x)
    }
    if (missing(ylim) == TRUE) {
        ylim = c(0, max_val_y)
    }
    plot(y = Outcome, x = Time, type = "n", ylim = ylim, xlim = xlim, 
        ...)
    if (Add.Profiles == TRUE) {
        for (i in 1:length(unique(Id))) {
            lines(y = Data$Outcome[Data$Id == unique(Data$Id)[i]], 
                x = Data$Time[Data$Id == unique(Data$Id)[i]], 
                col = Col)
        }
    }
    if (Add.Mean == TRUE) 
	{
    	 if (bytto==TRUE)
		{
	       #mean.tto=rep(0,length(ttovalue))
    		 for (i in 1:length(ttovalue)) {
		 mean.tto <- tapply(Data$Outcome[Data$tto==ttovalue[i]], INDEX = Data$Time[Data$tto==ttovalue[i]], FUN = mean, 
                               na.rm = TRUE)
            lines(mean.tto, x = unique(Data$Time), lwd = Lwd.Me, col=colme[i])
		}}else
		  {
    		   mean <- tapply(Data$Outcome, INDEX = Data$Time, FUN = mean, 
                       na.rm = TRUE)
              lines(mean, x = unique(Data$Time), lwd = Lwd.Me)
		  }
    }
    if (Add.Median == TRUE) 
	{
    	 if (bytto==TRUE)
		{
	       #median.tto=rep(0,length(ttovalue))
    		 for (i in 1:length(ttovalue)) {
		 median.tto <- tapply(Data$Outcome[Data$tto==ttovalue[i]], INDEX = Data$Time[Data$tto==ttovalue[i]], FUN = median, 
                                 na.rm = TRUE)
            lines(median.tto, x = unique(Data$Time), lwd = Lwd.Me,col=colme[i])
		}}else
		  {
    		   median<- tapply(Data$Outcome, INDEX = Data$Time, FUN = median, 
                       na.rm = TRUE)
              lines(median, x = unique(Data$Time), lwd = Lwd.Me)
		  }
    }
}


########################################################################
######################### Examples #####################################
########################################################################

#############################################################################################
############################## Early dietary intervention ###################################
#############################################################################################


#########################
############################## Generating the data. This create data that leads to problems.
#########################


## Reading in the early.int data

early.sim <- read.table("earlyint_pp.txt", header=T, sep=",")


## Mean:
early.mean.sim=tapply(cog,list(age,program),mean)
early.mean.sim[,2][1]=early.mean.sim[,1][1]

## Reshaping the data into a wide form

early.int2.sim <- reshape(early.sim[,2:5], 
  timevar = "age", idvar = c("id", "program"), direction = "wide")

## Sigma matrix and mu

mu.prog0=early.mean.sim[,1]
mu.prog1=early.mean.sim[,2]

Sigma.prog0<-var(early.int2.sim[early.int2.sim$program==0,3:5])
Sigma.prog1<-var(early.int2.sim[early.int2.sim$program==1,3:5])

Sigma.prog=var(early.int2.sim[,3:5])

set.seed(1)

## This generate data with different covariance sttructure per program. NOT USE

data.prog1=mvrnorm(n = 58, mu.prog1, Sigma.prog1)
data.prog0=mvrnorm(n = 45, mu.prog0, Sigma.prog0)

## This generate data with the same covariance sttructure per program. USE!!!

data.prog1=mvrnorm(n = 58, mu.prog1, Sigma.prog)
data.prog0=mvrnorm(n = 45, mu.prog0, Sigma.prog)

early.int1.sim=as.data.frame (rbind(data.prog1,data.prog0))
early.int1.sim$id=1:103
early.int1.sim$program=rep(c(1,0), c(58,45))

early.long.sim <- reshape(early.int1.sim, varying = c("1", "1.5", "2"), 
           v.names = "cog", timevar = "age", times = c("1", "1.5", "2"), 
           new.row.names = 1:309,direction = "long")

early.long.sim=early.long.sim[order(early.long.sim$id,early.long.sim$age),] 

write.table(early.long.sim, file="earlyint.txt", sep=",")

#########################
############################## Generating the data. Done one time. This is the one USE!!!
#########################

set.seed(1)

## Sigma matrix and mu

age=c(1,1.5,2)

mu.prog0=107.20-20.02*(age-1)
mu.prog1=early.mean.sim[,2]

Sigma.prog0<-var(early.int2.sim[early.int2.sim$program==0,3:5])
Sigma.prog1<-var(early.int2.sim[early.int2.sim$program==1,3:5])

Sigma.prog=var(early.int2.sim[,3:5])

set.seed(1)

## This generate data with different covariance sttructure per program. NOT USE

data.prog1=mvrnorm(n = 58, mu.prog1, Sigma.prog1)
data.prog0=mvrnorm(n = 45, mu.prog0, Sigma.prog0)

## This generate data with the same covariance sttructure per program. USE!!!

data.prog1=mvrnorm(n = 58, mu.prog1, Sigma.prog)
data.prog0=mvrnorm(n = 45, mu.prog0, Sigma.prog)

early.int1.sim=as.data.frame (rbind(data.prog1,data.prog0))
early.int1.sim$id=1:103
early.int1.sim$program=rep(c(1,0), c(58,45))

early.long.sim <- reshape(early.int1.sim, varying = c("1", "1.5", "2"), 
           v.names = "cog", timevar = "age", times = c("1", "1.5", "2"), 
           new.row.names = 1:309,direction = "long")

early.long.sim=early.long.sim[order(early.long.sim$id,early.long.sim$age),] 

write.table(early.long.sim, file="earlyint.txt", sep=",")















## Generating random effects

D=matrix(c(108.84,-30.419,-30.419,9.43),nrow=2,ncol=2,byrow=TRUE)
b=mvrnorm(n = 103, c(0,0), D)

## Generating the error

error=rnorm(103,mean=0,sd=8.61)

## Generating program

program=rep(c(1,0), c(58,45))

## Generating slope and intercept per subject

inti=107.20+b[,1]
slopei=-20.02+7.64*program+b[,2]

## Generating response per age and subject

age1.temp=inti+slopei*1
age1.5.temp=

#########################

#########################
############################## Analysis of the simulated data.
#########################

## Reading the data

early.int1 <- read.table("earlyint.txt", header=T, sep=",")

early.int1.table <- xtable(early.int1[1:24,])
print(early.int1.table)

## Attach data to the search path

attach(early.int1)

## Spaghettiplot

n=length(unique(id))
interaction.plot(age,id,cog, xlab="Agein years", ylab="IQ", legend=F) 

# Plot individual profiles + mean with function

Spaghetti.Plot.new(Dataset=early.int1, Outcome=cog, Time=age, Id=id, tto=program,
			  Add.Profiles = TRUE, Add.Mean = TRUE, Add.Median = FALSE, Col = 8, 
			  colme=c(2,1), Lwd.Me = 3, xlim=c(1,2), ylim=c(40,160), bytto=TRUE) 


## Adding a loess curve per group. There is a problem probably due to small n
with (early.int1, {
lines (loess.smooth (age[program==0], cog[program==0], family = "gaussian"),lty = 3,col=4,lwd = 4)
lines (loess.smooth (age[program==1], cog[program==1], family = "gaussian"),lty = 4,col=5,lwd = 4)
})

# Plot individual profiles + mean  by hand

plot(age, cog, type = "n")

for (i in 1:length(unique(id)))
{ lines(y =early.int1$cog[early.int1$id == unique(early.int1$id)[i]], 
                x =early.int1$age[early.int1$id == unique(early.int1$id)[i]], 
                col = 8)
}

mean.prog0=tapply(cog[program==0], INDEX =age[program==0], FUN = mean, na.rm = TRUE)
lines(mean.prog0, x = unique(age), lwd = 3)

mean.prog1=tapply(cog[program==1], INDEX =age[program==1], FUN = mean, na.rm = TRUE)
lines(mean.prog1, x = unique(age), lwd = 4,col=2)


## Descriptives

## Mean:
early.mean=tapply(cog,list(age,program),mean)

## Standard deviation:
early.sd=tapply(cog,list(age,program),sd)

## Variance:
early.var=tapply(cog,list(age,program),var)

## Frequency:
early.n=table(age,program)

## Boxplots:

boxplot(cog~age,xlab="Age (in years)",ylab="IQ", col = "gray")

## Boxplots per program

par(mfrow=c(2,1))
boxplot(cog[program==0]~age[program==0],main="No intervention",main="No intervention",xlab="Age (in years)",ylab="IQ", col = "gray")
boxplot(cog[program==1]~age[program==1],main="Intervention",main="No intervention",xlab="Age (in years)",ylab="IQ", col = "gray")

## Plotting mean evolutions

plot(age[id==1],early.mean[,1],type="b",xlim=c(1,2),
ylim=c(40,160),xlab="Age (in years)",ylab="IQ",axes=F,
main="Mean evolution (with 1 SE intervals)")
axis(side=1,at=c(1,1.5,2),labels=c(1,1.5,2))
axis(side=2,at=seq(40,160,20))

box()
points(age[id==1],early.mean[,2],type="b",col="red")
errbar(age[id==1]-.005,early.mean[,1],early.sd[,1],.1)
errbar(age[id==1]+.005,early.mean[,2],early.sd[,2],.1,col="red")

## Correlations

## Reshaping the data into a wide form

early.int2 <- reshape(early.int1, 
  timevar = "age", idvar = c("id", "program"), direction = "wide")

## Correlation between the IQ scores  at  different ages

cor(early.int2[,3:5])
cor(early.int2[early.int2$program==1,3:5])
cor(early.int2[early.int2$program==0,3:5])

## Linear regression per person + histograms

## Creating the time variable

early.int1$age0<-early.int1$age-1

## Displaying the linear regression per person:

cf<-sapply(early.int1$id, function(x) 
    coef(lm(cog~age0, data=subset(early.int1, id==x))))

Sx<-reorder(early.int1$id, cf[1,])

xyplot(cog ~ age0|Sx,groups=program,data=early.int1,
 type=c('p','r'),auto.key=T,aspect="xy",
 par.settings=list(axis.text=list(cex=0.6),
 fontsize=list(text=8, points=10))
)

## Linear regression per participant of cog on age


## Coefficients

lin.reg.coef <- by(early.int1, early.int1$id, 
          function(data) coef(lm(cog ~ age0, data=data)))
lin.reg.coef1 <- unlist(lin.reg.coef)
names(lin.reg.coef1) <- NULL 
lin.reg.coef2=matrix(lin.reg.coef1,length(lin.reg.coef1)/2,2,byrow = TRUE)

## R squared

lin.reg.r.squared <- by(early.int1, early.int1$id, 
          function(data) summary(lm(cog ~ age, data=data))$r.squared )
lin.reg.r.squared1<- as.vector(unlist(lin.reg.r.squared))

## Histograms

par(mfrow=c(3,1))
hist(lin.reg.coef2[,1],xlab="Intercept",col="lightblue",main="Histogram of individual intercepts")
hist(lin.reg.coef2[,2],xlab="Slope",col="lightblue",main="Histogram of individual slopes")
hist(lin.reg.r.squared1,xlab="R squared",col="lightblue",main="Histogram of individual R squared")

##   Correlations

int.slope.corr=cor(lin.reg.coef2[,1],lin.reg.coef2[,2])
plot(lin.reg.coef2[,1],lin.reg.coef2[,2],xlab="Intercept", ylab="Slope", main="Intercept versus Slope")

## Plotting individual regression lines per group

reg.coef=cbind(lin.reg.coef2, early.int1[early.int1$age==1,]$program)

mean.int<-tapply(reg.coef[,1],reg.coef[,3],mean)
mean.slope<-tapply(reg.coef[,2],reg.coef[,3],mean)

par(mfrow=c(1,2))
plot(early.int1$age0,early.int1$cog,type="n",xlim=c(0,1),ylim=c(40,160),main="No intervention",xlab="Age (in years)",ylab="IQ",axes=F)
axis(side=1,at=c(0,0.5,1),labels=c(0,0.5,1))
axis(side=2,at=seq(40,160,20))
box()
for (i in 1:103)
{if (reg.coef[i,3]==0) 
{curve(cbind(1,x)%*%reg.coef[i,1:2],add=T,col="gray")}}
curve(cbind(1,x)%*%c(mean.int[1],mean.slope[1]),add=T,lwd=2)

plot(early.int1$age0,early.int1$cog,type="n",xlim=c(0,1),ylim=c(40,160),main="Intervention",xlab="Age (in years)",ylab="IQ",axes=F)
axis(side=1,at=c(0,0.5,1),labels=c(0,0.5,1))
axis(side=2,at=seq(40,160,20))
box()
for (i in 1:103)
{if (reg.coef[i,3]==1) 
{curve(cbind(1,x)%*%reg.coef[i,1:2],add=T,col="gray")}}
curve(cbind(1,x)%*%c(mean.int[2],mean.slope[2]),add=T,lwd=2)


## Fitting the model with ML

## Different random intercept and slope

early.lmer1<-lmer(cog~1+age0*program+(1 + age0|id), REML = FALSE, data=early.int1)
mcp.fnc(early.lmer1)

summary(early.lmer1)
display(early.lmer1)
anova(early.lmer1)

## Predicted random effects

early.lmer1.re=ranef(early.lmer1)[["id"]]
head(early.lmer1.re)
cor(early.lmer1.re[,1],early.lmer1.re[,2])
cor(early.lmer1.re[,1],early.lmer1.re[,2])

plot(early.lmer1.re[,1],lin.reg.coef2[,1])
plot(early.lmer1.re[,2],lin.reg.coef2[,2])

## Paradox for the random effects:
## The random effects from the LMM are positively correlated with the intercept and slope
## obtained from fitting fixed linear regressions to each subject. However, the random effects
## from the LMM are positively correlated (correlation 1) and the estimated intercept and slope
## obtained from the fixed linear regressions per subject are negatively correlated!!!

par(mfrow=c(2,2))
plot(lin.reg.coef2[,1],lin.reg.coef2[,2])
plot(early.lmer1.re[,1],early.lmer1.re[,2])
plot(early.lmer1.re[,1],lin.reg.coef2[,1])
plot(early.lmer1.re[,2],lin.reg.coef2[,2])

int.slope.corr.prg0=cor(reg.coef[,1][reg.coef[,3]==0],reg.coef[,2][reg.coef[,3]==0])
int.slope.corr.prg1=cor(reg.coef[,1][reg.coef[,3]==1],reg.coef[,2][reg.coef[,3]==1])

## Same random intercept and slope

early.lmer1.sameintslope<-lmer(cog~1+age0*program+(-1 + age0|id), REML = FALSE, data=early.int1)

summary(early.lmer1.sameintslope)
display(early.lmer1.sameintslope)
anova(early.lmer1.sameintslope)

## Profiling the  likelihood

pr01.sameintslope<- profile(early.lmer1.sameintslope)
xyplot(pr01.sameintslope, aspect = 1.3)
splom(pr01.sameintslope)

## Predicting the random effects

ranef(early.lmer1.sameintslope)
pred.re.sameintslope=ranef(early.lmer1.sameintslope,condVar = TRUE)
dotplot(pred.re.sameintslope)
qqmath(pred.re.sameintslope)


## Estimating the fixed effects via bootstrap

fixed.boot=bootMer(early.lmer1,  fixef, use.u = TRUE, nsim = 250)
fixed.boot
summary(fixed.boot)

#help(pvalues)

## Calculating confidence intervals for the fixed effects via bootstrap

confint(early.lmer1.sameintslope, level = 0.95,method="profile",oldNames = FALSE)
confint(early.lmer1.sameintslope,method="Wald",oldNames = FALSE)
confint(early.lmer1.sameintslope,method="boot",boot.type ="perc",oldNames = FALSE,nsim=500)
confint(early.lmer1.sameintslope,method="boot",boot.type ="basic",oldNames = FALSE,nsim=500)

## Get the KR-approximated degrees of freedom

early.lmer1.df.KR <- get_ddf_Lb(early.lmer1.sameintslope, fixef(early.lmer1))

## Get p-values from the t-distribution using the t-values and approximated
## degrees of freedom

early.lmer1.coef=coef(summary(early.lmer1.sameintslope))
early.lmer1.p.KR <- cbind(early.lmer1.coef,2 * (1 - pt(abs(early.lmer1.coef[,3]), early.lmer1.df.KR)))
early.lmer1.p.KR

## Likelihood ratio tests

early.lmer1.noprog<-lmer(cog~1+age0+(1 + age0|id), REML = FALSE, data=early.int1)
early.lmer1.int<-lmer(cog~1+age0+program+(1 + age0|id), REML = FALSE, data=early.int1)
anova(early.lmer1.noprog,early.lmer1.int,early.lmer1)

early.lmer1.noint<-lmer(cog~1+age0+age0:program+(1 + age0|id), REML = FALSE, data=early.int1)
anova(early.lmer1.noint,early.lmer1)
anova(early.lmer1.noint,early.lmer1.noprog)


confint(early.lmer1, level = 0.99)
confint(early.lmer1,method="Wald")
confint(early.lmer1,method="boot",boot.type ="perc")


warnings()


RLRsim(early.lmer1)



###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################






str(alcohol)
dim(alcohol)

alcohol.table <- xtable(alcohol[1:10,])
print(alcohol.table)

## Exploring the alcohol data


# Boxplots and mean-plot per alcohol group

boxplot(alcohol$Attractiveness~alcohol$Alcohol, horizontal=F, col = "gray")
plotmeans(Attractiveness~Alcohol,xlab="Alcohol", ylab="Attractiveness", main="Mean Plot\nwith 95% CI", data = alcohol) 


# Density plot per group

simple.densityplot(Attractiveness~Alcohol, data = alcohol)

# Means per factors

plot.design(Attractiveness ~ Alcohol*Gender, data = alcohol, xlab = "")
coplot(Attractiveness~Alcohol | Gender, data=alcohol, panel=panel.smooth, xlab="Alcohol")

# Plotting the data per alcohol group

stripchart(alcohol$Attractiveness~alcohol$Alcohol,vertical=F, xlab='Attractiveness')
meanvec = tapply(alcohol$Attractiveness,alcohol$Alcohol,mean)
abline(h=1:3,lty=2,col="black")
points(meanvec,1:3 , pch=17, col="red", cex=2)

####################
#################### Analyzing the alcohol data: One-Way ANOVA
####################

attract.fit<-aov(Attractiveness~Alcohol,data = alcohol)
summary(attract.fit)

## Creating tables in Latex for the anova output

attract.fit.table <- xtable(attract.fit)
print(attract.fit.table)

## Testing normality

hist(rstandard(attract.fit))
attract.fit.shapiro<-shapiro.test(alcohol$Attractiveness)
ks.test(alcohol$Attractiveness, "plnorm")

## Creating tables in Latex for the shapiro-test

attract.fit.shapiro[1]
attract.fit.shapiro.table<-xtable(as.data.frame(W=attract.fit.shapiro$W[1] ,p.value=attract.fit.shapiro$p.value))


## Testing homocedasticity

# library(lawstat)
levene.test(alcohol$Attractiveness,alcohol$Alcohol,location="median")


## Ploting residuals

plot(attract.fit) # escape to get each graph separated

## Studentized residuals

plot(fitted.values(attract.fit),rstandard(attract.fit), xlab="Fitted values",
     ylab="Studentized residuals",
     main="Residuals vs fitted values plot")
     abline(h=0,lty="dashed")


plot(alcohol$Alcohol,rstandard(attract.fit), xlab="Observed values",
     ylab="Studentized residuals",
     main="Residuals vs factor levels plot")

abline(h=0,lty="dashed")

## Time trends

plot(rstandard(attract.fit)[-c(1)],rstandard(attract.fit)[-c(dim(alcohol)[1])],
     xlab="Studentized residuals at 1 lag",
     ylab="Studentized residuals",
     main="Sequence plot")
abline(a=0,b=1,lty="dashed")

## QQ-plot

qqnorm(residuals(attract.fit))
qqline(residuals(attract.fit))

## Studentized residuals stratified by Gender

alcohol$rs<-rstandard(attract.fit)
alcohol$fit<-fitted.values(attract.fit)

# Create new column filled with default colour
alcohol$Colour="black"

# Set new column values to appropriate colours

alcohol$Colour[alcohol$Gender=="Female"]="red"
alcohol$Colour[alcohol$Gender=="Male"]="blue"

# Plot all points at once, using newly generated colours

plot(alcohol$fit,alcohol$rs, col=alcohol$Colour, pch=16, xlab="Fitted values",
     ylab="Studentized residuals",
     main="Residuals vs fitted values plot per gender")
abline(h=0,lty="dashed")

plot(jitter(alcohol$fit),jitter(alcohol$rs), col=alcohol$Colour, pch=16, xlab="Fitted values",
     ylab="Studentized residuals",
     main="Residuals vs fitted values plot per gender")
abline(h=0,lty="dashed")

legend(53,2.2,pch=c(19,19),lty=c(1,1),col=c("red","blue"),legend=c("Female","Male"))

## Outliers

pvalue_outliers = NULL
r.al=3
nT.al<-dim(alcohol)[1]
for(i in 1:nT.al){
pvalue_outliers[i]=1-pt(abs(rstudent(attract.fit)[i]),nT.al-r.al-1)
}
pvalue_outliers[pvalue_outliers>(0.05/(2*nT.al))]=1

out.alcohol<-data.frame(Stud.Deleted.Res=rstudent(attract.fit),Outlier.p.value=pvalue_outliers)

## Normality

shapiro.test(residuals(attract.fit))

ks.test(residuals(attract.fit),"pnorm",alternative="two.sided")

## Box-Cox

attract.boxcox<-boxcox(attract.fit,lambda=seq(0,4,by=.1))
lambda.attract.max<-max(attract.boxcox$y)
attract.boxcox$x[attract.boxcox$y==lambda.attract.max]

 ## Kruskal-Wallis test

kruskal.test(Attractiveness~Alcohol,data = alcohol, na.action=na.omit)

####################
#################### Analyzing the alcohol data: Two-Way ANOVA
####################

# Exploring the data

op0 = par()    # Get current graphical parameters
op1 = op0$mar  # Get current margins in lines
op1[1] = 8     # Modify bottom margins to 20 (lines)
op1 
par(mar = op1)

boxplot(Attractiveness~Alcohol*Gender,data = alcohol, main="Boxplots of Attractiveness Data", col = "gray", las=3)


par(mar = op0)
with(alcohol, interaction.plot(x.factor=Alcohol, trace.factor=Gender, response=Attractiveness, fun=mean, type="b", legend=T,
                 ylab="Attractiveness Means", main="Interaction Plot", pch=c(1,19)))

# Fitting the model

attract.fit.two<-aov(Attractiveness~Alcohol*Gender,data = alcohol)
attract.fit.two.summary<-summary(attract.fit.two)

# Fitting the model with other codes for the factors

alcohol$Al<-'0'
alcohol$Al[alcohol$Alcohol=='2 Pints']<-'2'
alcohol$Al[alcohol$Alcohol=='4 Pints']<-'4'

alcohol$Ge<-'M'
alcohol$Ge[alcohol$Gender=='Female']<-'F'

alcohol$Al<-as.factor(alcohol$Al)
alcohol$Ge<-as.factor(alcohol$Ge)

# Creating tables in Latex for the anova output. Sometimes if xtable cannot read a format it is useful to use the unclass statement like this example shows.

table.mean.alcohol.temp<-unclass(model.tables(attract.fit.two, type="means", se=T)$tables$`Alcohol:Gender`)
table.mean.alcohol<-xtable(table.mean.alcohol.temp , include.rownames=FALSE, include.colnames=TRUE)
print(table.mean.alcohol)


# Multiple comparison


attract.fit.two.tukey<-TukeyHSD(attract.fit.two,which=c("Alcohol:Gender"), conf.level=.95)

op0 = par()    # Get current graphical parameters
op1 = op0$mar  # Get current margins in lines
op1[2] = 13     # Modify bottom margins to 20 (lines)
op1 
par(mar = op1)
plot(attract.fit.two.tukey, las=1)


# Creating tables in Latex for the multiple comparison

attract.fit.two.tukey.temp<-xtable(attract.fit.two.tukey$`Alcohol:Gender`, include.rownames=FALSE, include.colnames=TRUE)
print(attract.fit.two.tukey.temp)


plot(TukeyHSD(aov(Attractiveness~Al*Ge,data = alcohol),which=c("Al:Ge"), conf.level=.95))

summary(aov(Attractiveness~Al*Ge,data = alcohol))

## New analysis

## Exploring the data

attach(alcohol)

## Interaction plot

interaction.plot(Alcohol,Gender,Attractiveness, xlab="Treatment")

## Other plots

plot.design(Attractiveness~Alcohol+Gender, data = alcohol , fun = mean)

## Estimating the models and performing the F tests:

## define the contrasts

contrasts(alcohol$Alcohol) = contr.sum
contrasts(alcohol$Gender) = contr.sum

## Fit the initial model

alcohol.fit2= lm(Attractiveness~Alcohol*Gender,data=alcohol)
summary(alcohol.fit2)

## Performing the anova tests

anova(alcohol.fit2)  # interaction effect is significant

## Performing multiple comparisons

## Tukey all treatments

alcohol$treatments =alcohol$Alcohol:alcohol$Gender
alcohol.fit2.treatments= lm(Attractiveness~treatments,data=alcohol)
tukey_treatments = glht(alcohol.fit2.treatments, linfct = mcp(treatments = "Tukey"))

tukey_ci_treatments = confint(tukey_treatments)
tukey_ci_treatments

op0 = par()    # Get current graphical parameters
op1 = op0$mar  # Get current margins in lines
op1[2] = 13     # Modify bottom margins to 20 (lines)
op1 
par(mar = op1)
plot(tukey_ci_treatments)

## Estimable functions

## Analysis under different restrictions

library(car)
alcohol.estimable.sum=lm(Attractiveness~Alcohol*Gender,contrasts=list(Alcohol='contr.sum',Gender='contr.sum'),data = alcohol)
summary(alcohol.estimable.sum)
Anova(alcohol.estimable.sum)

alcohol.estimable.treatment=lm(Attractiveness~Alcohol*Gender,contrasts=list(Alcohol='contr.treatment',Gender='contr.treatment'),data = alcohol)
summary(alcohol.estimable.treatment)
Anova(alcohol.estimable.treatment)

detach(alcohol)

#####################################################################################
############################ Productivity Improvement ###############################
#####################################################################################

data = read.table("CH16PR07.txt",col.names=c("y","group","j_index")) 
str(data)

## Plotting the observations

attach(data)
group_means = c(mean(data$y[data$group==1]),mean(data$y[data$group==2]),mean(data$y[data$group==3]))


plot(data$group,data$y,xlab="Expenditures for research and development",
             ylab="Productivity improvement",axes=F)
points(c(1,2,3),group_means,pch=19,cex=1.2)
axis(side=1, at = c(1,2,3), labels = c("Low","Moderate","High"))
axis(side=2)

## Estimating $\mu_i$, $\sigma^2$, calculating SSE,SSTR and SSTO

# means

mu_i = group_means

# SSE,SSTR and SSTO

sum_i1_sstr = length(y[group==1])*((mu_i[1]-mean(y))^2)
sum_i2_sstr = length(y[group==2])*((mu_i[2]-mean(y))^2)
sum_i3_sstr = length(y[group==3])*((mu_i[3]-mean(y))^2)
sum_i_sstr = c(sum_i1_sstr,sum_i2_sstr,sum_i3_sstr)
sum_i_sstr 
[1] 10.3827160  0.3952263  9.3472428
sstr = sum(sum_i_sstr )

sum_i1_sse = sum((y[group==1]-mu_i[1])^2)
sum_i2_sse = sum((y[group==2]-mu_i[2])^2)
sum_i3_sse = sum((y[group==3]-mu_i[3])^2)
sum_i_sse = c(sum_i1_sse,sum_i2_sse,sum_i3_sse)

sse = sum(sum_i_sse)

ssto = sum((y-mean(y))^2)
ssto 
[1] 35.48741

# $\sigma^2$

sigma2 = sse/(length(y)-3)   # we'll see it later
sigma2

# Fitting the one-way ANOVA model

data.fit<-aov(y~as.factor(group),data = data)
summary(data.fit)

# line plot
 
plot(mu_i,c(1,1,1), axes=F,pch=19,cex=3, ylim=c(0,2.100),xlim=c(6,10),
     xlab="Productivity improvement",  ylab=" ",col="grey")
axis(1)
abline(h=1)
points(mu_i,c(1,1,1),pch=19,cex=3,col="grey")
text(6.85,1.3,"Low")
text(8.1,1.3,"Moderate")
text(9.18,1.3,"High")

# Bar  plot
barplot(height=mu_i,width=1,space=1.2,  ylim=c(0,10), names.arg=c("Low","Moderate","High"),
    xlab="Expenditures for research and development",
    ylab="Productivity improvement",col="grey")

# main effects plot

plot(c(1,2,3),mu_i, axes=F, ylim=c(0,10), xlab="Expenditures for research and development",
    ylab="Productivity improvement",type="l")
abline(h=mean(mu_i),lty="dashed")
points(c(1,2,3),mu_i,pch=19,cex=3,col="grey")
axis(1,at=c(1,2,3),labels=c("Low","Moderate","High"))
axis(2)

# Factor Level Means: t-test

## Single factor level mean

alpha = 0.05
n_T = length(data$y)
r = length(mu_i)
n_1 = sum(data$group==1)
c = 8
 
# confidence interval

ci_mu1_single =  c(mu_i[1]+qt(alpha/2,n_T-r)*sqrt(sigma2/n_1),mu_i[1]+qt(1-alpha/2,n_T-r)*sqrt(sigma2/n_1))
ci_mu1_single


# two-sided test (H_1: mu_1 \neq 8)
t.test(y[as.factor(group)=="1"],mu=8,conf.level = 0.95)

# one-sided test (H_1: mu_1 < 8)
t.test(y[as.factor(group)=="1"],mu=8,conf.level = 0.95,alternative="less")

# two-sided test (H_1: mu_2 - mu_1 \neq 0)
t.test(y[as.factor(group)=="2"],y[as.factor(group)=="1"],conf.level = 0.95,var.equal=T)

#  # one-sided test (H_1: mu_2 - mu_1 > 0)
t.test(y[as.factor(group)=="2"],y[as.factor(group)=="1"],conf.level = 0.95,var.equal=T,alternative="greater")

# Tukey Method I

tukey.t.data<-TukeyHSD(data.fit)
tukey.t.data 
plot(tukey.t.data)

# Tukey Method II

data$group1 = factor(data$group,labels=c("Low","Moderate","High"))
data.fit1<-aov(y~group1,data = data)
data.fit1.summary<-summary(data.fit1)
str(summary(data.fit1))

tukey = glht(data.fit1, linfct = mcp(group1= "Tukey"))
tukey_ci = confint(tukey)

## Sheffé

all_pairwise = rbind("2 - 1" = c(-1, 1, 0),
                     "3 - 1" = c(-1, 0, 1), 
                     "3 - 2" = c(0, -1, 1))
n_i=c(length(y[group==1]),length(y[group==2]),length(y[group==3]))
y = data$y
tested_contrasts = t(all_pairwise)
n_T=sum(n_i)
r=data.fit1.summary[[1]]$Df[1]+1
sse=data.fit1.summary[[1]]$'Sum Sq'[2] 
sheffe = sheffe_fun(y,mu_i,tested_contrasts,sse,n_T,r,n_i,alpha=0.05)

# Bonferoni Method

bonferroni = glht(data.fit1, linfct = mcp(group1 = all_pairwise), test = adjusted("bonferroni"))
summary(bonferroni, test = adjusted("bonferroni"))
bonferroni_ci = confint(bonferroni)

pairwise.t.test(y, as.factor(group), p.adj = "bonf")

## Graphical comparison between the three methods

plot(tukey_ci$confint[,1]~c(0.9,1.9,2.9),xlab="Expenditures for research and development",
	ylab="Productivity improvement",axes=F,
	ylim=c(min(tukey_ci$confint,sheffe$ci_sheffe,bonferroni_ci$confint),
	max(tukey_ci$confint,sheffe$ci_sheffe,bonferroni_ci$confint)),
	xlim=c(0.5,3.5))

group_means = c(mean(y[group==1]),mean(y[group==2]),mean(y[group==3]))
points(tukey_ci$confint[,1]~c(0.9,1.9,2.9),pch=19,cex=1.2)
points(sheffe$ci_sheffe[,2]~c(1,2,3),pch=19,cex=1.2,col=2)
points(bonferroni_ci$confint[,1]~c(1.1,2.1,3.1),pch=19,cex=1.2,col=3)

axis(side=1, at = c(1,2,3), labels = c("Low","Moderate","High"))
axis(side=2)
abline(h=0,lty="dashed")

segments(0.9,tukey_ci$confint[1,2],0.9,tukey_ci$confint[1,3])
segments(1,sheffe$ci_sheffe[1,1],1,sheffe$ci_sheffe[1,3],col=2)
segments(1.1,bonferroni_ci$confint[1,2],1.1,bonferroni_ci$confint[1,3],col=3)
segments(1.9,tukey_ci$confint[2,2],1.9,tukey_ci$confint[2,3])
segments(2,sheffe$ci_sheffe[2,1],2,sheffe$ci_sheffe[2,3],col=2)
segments(2.1,bonferroni_ci$confint[2,2],2.1,bonferroni_ci$confint[2,3],col=3)
segments(2.9,tukey_ci$confint[3,2],2.9,tukey_ci$confint[3,3])
segments(3,sheffe$ci_sheffe[3,1],3,sheffe$ci_sheffe[3,3],col=2)
segments(3.1,bonferroni_ci$confint[3,2],3.1,bonferroni_ci$confint[3,3],col=3)

legend(0.5,3.5,pch=c(19,19,19),lty=c(1,1,1),col=c(1,2,3),legend=c("Tukey","Sheffe","Bonferroni"))

## Residuals and remedial measures

## Analysis of residuals
##  Residuals
 
n_groups = c(rep(n_i[1],sum(data$group==1)), rep(n_i[2],sum(data$group==2)), rep(n_i[3],sum(data$group==3)))

e_star = residuals(data.fit1)/sqrt(sse/(n_T-r))
r_stud = residuals(data.fit1)/sqrt((sse/(n_T-r))*(1-1/n_groups))
t_deleted = residuals(data.fit1)*sqrt((n_T-r-1)/(sse*(1-1/n_groups)-residuals(data.fit1)^2))

res.all<- data.frame("e"=data$y-fitted.values(data.fit1),"e with R"=residuals(data.fit1),"e star"=e_star,
                     "r stud"=r_stud,"r stud with R"=rstandard(data.fit1),
                     "t deleted"=t_deleted,"t deleted with R"=rstudent(data.fit1))

## Studentized residuals

plot(fitted.values(data.fit1),rstandard(data.fit1), xlab="Fitted values",
     ylab="Studentized residuals",
     main="Residuals vs fitted values plot")
     abline(h=0,lty="dashed")


plot(data$group,rstandard(data.fit1), xlab="Observed values",
     ylab="Studentized residuals",
     main="Residuals vs factor levels plot",axes=F)

axis(side=1, at = c(1,2,3), labels = c("Low","Moderate","High"))
axis(side=2)
abline(h=0,lty="dashed")

plot(data$group1,rstandard(data.fit1), xlab="Observed values",
     ylab="Studentized residuals",
     main="Residuals vs factor levels plot")
abline(h=0,lty="dashed")

## Time trends

plot(rstandard(data.fit1)[-c(1)],rstandard(data.fit1)[-c(n_T)],
     xlab="Studentized residuals at 1 lag",
     ylab="Studentized residuals",
     main="Sequence plot")
abline(a=0,b=1,lty="dashed")

## QQ-plots

qqnorm(residuals(data.fit1))
qqline(residuals(data.fit1))

## Testing homocedasticity

# library(lawstat)
levene.test(data$y,data$group1,location="median")

## Outliers

pvalue_outliers = NULL
#aaa = NULL
for(i in 1:n_T){
pvalue_outliers[i]=1-pt(abs(rstudent(data.fit1)[i]),n_T-r-1)
#pvalue_outliers[i] = ((1-pt(abs(rstudent(data.fit1)[i]),n_T-r-1))*2)*n_T
}
pvalue_outliers[pvalue_outliers>(0.05/(2*n_T))]=1
#pvalue_outliers[pvalue_outliers>1] = 1

out.data<-data.frame(Stud.Deleted.Res=rstudent(data.fit1),Outlier.p.value=pvalue_outliers)

## Normality

shapiro.test(residuals(data.fit1))

ks.test(residuals(data.fit1),"pnorm",alternative="two.sided")

## Box-Cox transformation

data.boxcox<-boxcox(data.fit1,lambda=seq(0,4,by=.1))
lambda.max<-max(data.boxcox$y)
data.boxcox$x[data.boxcox$y==lambda.max]

## Kruskal-Wallis test

kruskal.test(y~group1, data = data, na.action=na.omit)

## Regression approach

## How to define the contrasts in R. Model matrixes per type of contrast

data$group1 = factor(data$group,labels=c("Low","Moderate","High"))

model.matrix(~ group1, data);

model.matrix(~ group1, data, contrasts = list(group1="contr.treatment"))

model.matrix(~ -1+group1, data)

model.matrix(~ group1, data, contrasts = list(group1="contr.sum"))

## Contrast matrixes per type of contrast

## First group reference

contrasts(data$group1) = contr.treatment
contrasts(data$group1)

## Unweighted means

contrasts(data$group1) = contr.sum
contrasts(data$group1)

## Weighted means

n_i=c(length(y[group==1]),length(y[group==2]),length(y[group==3]))
C = matrix(c(1,0,-(n_i[1]/n_i[3]),0,1,-(n_i[2]/n_i[3])),byrow=F,nrow=3)
contrasts(data$group1) = C
model.matrix(~ group1, data)

## Factor effects model
## Unweighted mean
##Define the contrasts

contrasts(data$group1) = contr.sum
mod_fem.um = lm(y~group1,data=data)
summary(mod_fem.um)

## Regression approach Unweighted mean

mu_i = c(mean(data$y[data$group==1]),mean(data$y[data$group==2]),mean(data$y[data$group==3]))

data.frame("lm function"=c(coef(mod_fem.um)[1],coef(mod_fem.um)[1]+coef(mod_fem.um)[2],
         coef(mod_fem.um)[1]+coef(mod_fem.um)[3],
         coef(mod_fem.um)[1]-coef(mod_fem.um)[2]-coef(mod_fem.um)[3]),
         "by hand"=c(mean(mu_i),mu_i[1],mu_i[2],mu_i[3]),row.names=c("mu","mu_1","mu_2","mu_3"))

## Perform the anova test

## Model under H_0

mod_H0 = lm(y~1,data=data)   

anova(mod_fem.um,mod_H0)


## Weighted mean
## Define the contrasts

C = matrix(c(1,0,-(n_i[1]/n_i[3]),0,1,-(n_i[2]/n_i[3])),byrow=F,nrow=3)
contrasts(data$group1) = C

## Fit the model

mod_fem.wm = lm(y~group1,data=data)
summary(mod_fem.wm)

## Regression approach Weighted mean

data.frame("lm function"=c(coef(mod_fem.wm)[1],coef(mod_fem.wm)[1]+coef(mod_fem.wm)[2],
         coef(mod_fem.wm)[1]+coef(mod_fem.wm)[3],
         coef(mod_fem.wm)[1]-coef(mod_fem.wm)[2]-coef(mod_fem.wm)[3]),
         "by hand"=c(mean(data$y),mu_i[1],mu_i[2],mu_i[3]),row.names=c("mu","mu_1","mu_2","mu_3"))

## Perform the anova test

anova(mod_fem.wm,mod_H0)

##  Cell means model
##  Define the contrasts

contrasts(data$group1) = contr.treatment

## Fit the model

mod_cmm = lm(y~-1+group1,data=data)
summary(mod_cmm)

## Regression approach cell mean

data.frame("lm function"=c(coef(mod_cmm)[1],coef(mod_cmm)[2],coef(mod_cmm)[3]),
         "by hand"=c(mu_i[1],mu_i[2],mu_i[3]),row.names=c("mu_1","mu_2","mu_3"))
     lm.function  by.hand

## Perform the anova test

anova(mod_cmm,mod_H0)

## Estimable functions

## Means per group

tapply(y, group1, mean)

## Mean of the means

sum(tapply(y, group1, mean))/3

## Analysis under different restrictions

library(car)
prod.estimable.sum=lm(y~group1,contrasts=list(group1='contr.sum'),data = data)
summary(prod.estimable.sum)
Anova(prod.estimable.sum)

prod.estimable.treatment=lm(y~group1,contrasts=list(group1='contr.treatment'),data = data)
summary(prod.estimable.treatment)
Anova(prod.estimable.treatment)



detach(data)

#############################################################################################
##################################  Alertness level  ########################################
#############################################################################################

## Reading and saving the data

datafilename="http://personality-project.org/r/datasets/R.appendix2.data"
data.ex2=read.table(datafilename,header=T)    #read the data into a table
data.ex2                                      #show the data
summary(data.ex2)

# FPath.alertness  <- "C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\ANOVA-Ariel\\Data-Codes\\alertness.txt" 
# write.table(data.ex2, file ="alertness.txt", sep = "\t")

data.ex2 = read.table("alertness.txt") 

## Exploring the data

boxplot(Alertness~Gender*Dosage,data=data.ex2, main="Boxplots of Alertness Data", col = "gray")
with(data.ex2, interaction.plot(x.factor=Dosage, trace.factor=Gender, response=Alertness,
 fun=mean, type="b", legend=T, ylab="Alertness Means",
 main="Interaction Plot", pch=c(1,19)))

## Fitting the model

aov.ex2 = aov(Alertness~Gender*Dosage,data=data.ex2) #do the analysis of variance
summary(aov.ex2)                                     #show the summary table
print(model.tables(aov.ex2,"means"),digits=3)        #report the means and the number of subjects/cell



#####################################################################################
################################# Toxicity Rats #####################################
#####################################################################################

data.ex2 = read.table("rats.txt", header=TRUE) 
str(data.ex2)
attach(data.ex2)

summary(data.ex2)

## Boxplot

boxplot(time~treat,data=data.ex2, main="Survival versus treatment", col = "gray")
boxplot(time~poison,data=data.ex2, main="Survival versus poison", col = "gray")

## Interaction plot

interaction.plot(treat,poison,time, xlab="Treatment")

## Other plots

plot.design(time~poison+treat, data = data.ex2, fun = mean)


plot(time~treat + poison, data=data.ex2)	

## Estimating the models and performing the F tests:

# Define the contrasts

contrasts(data.ex2$treat) = contr.sum
contrasts(data.ex2$poison) = contr.sum

## Fit the initial model

aov.lm.ex2= lm(time~treat*poison,data=data.ex2)
summary(aov.lm.ex2)

## Perform the anova test

anova(aov.lm.ex2)  # interaction effect is not significant

## Performing multiple comparisons

library(multcomp)

## Tukey for factor treat excluding factor poison

toxicity_ex2_treat= lm(time~treat,data=data.ex2) # ignoring factor A
tukey_treat= glht(toxicity_ex2_treat, linfct = mcp(treat= "Tukey"))
tukey_ci_treat= confint(tukey_treat)
tukey_ci_treat

## Tukey for factor treat including factor poison+treat+poison:treat

tukey_treat_all= glht(aov.lm.ex2, linfct = mcp(treat= "Tukey"))
tukey_ci_treat_all= confint(tukey_treat_all)
tukey_ci_treat_all
plot(tukey_ci_treat_all)

## Tukey for factor poison including factor poison+treat+poison:treat

tukey_poison_all= glht(aov.lm.ex2, linfct = mcp(poison= "Tukey"))
tukey_ci_poison_all= confint(tukey_poison_all)
tukey_ci_poison_all
plot(tukey_ci_poison_all)

## Dropping the interaction

aov.lm.ni.ex2= lm(time~treat+poison,data=data.ex2)
anova(aov.lm.ni.ex2)

detach(rats)

#####################################################################################
######################### Synthetic growth hormone ##################################
#####################################################################################

####################
#################### Synthetic growth hormone: Unbalanced
####################

growth.unbalanced<-read.table("CH23TA01.txt", col.names=c("Growth", "SexO", "BDO", "ID"))

growth.unbalanced$Sex<-'Male'
growth.unbalanced$Sex[growth.unbalanced$SexO==2]<-'Female'
growth.unbalanced$Sex<-as.factor(growth.unbalanced$Sex)

growth.unbalanced$BD<-'Mildly'
growth.unbalanced$BD[growth.unbalanced$BDO==2]<-'Moderately'
growth.unbalanced$BD[growth.unbalanced$BDO==3]<-'Severely'
growth.unbalanced$BD<-as.factor(growth.unbalanced$BD)

attach(growth.unbalanced)


# Summary tables

freq.unbalanced<-table(growth.unbalanced$Sex,growth.unbalanced$BD)
mean.freq.unbalanced<-tapply(growth.unbalanced$Growth, list(growth.unbalanced$Sex,growth.unbalanced$BD), mean)

freq.table.unbalanced<-xtable(freq.unbalanced)
print(freq.table.unbalanced)

mean.freq.table.unbalanced<-xtable(mean.freq.unbalanced)
print(mean.freq.table.unbalanced)

# Interaction plot

interaction.plot(BD,Sex,Growth, xlab="Treatment")

## Other plots

plot.design(Growth~Sex+BD, data = growth.unbalanced, fun = mean)


#### Regression approach

## Model matrixes per type of contrast

model.matrix(~ -1+Sex:BD, growth.unbalanced)

## Regression model

growth.lm.unbalanced=lm(Growth~-1+Sex:BD, data = growth.unbalanced) 
summary(growth.lm.unbalanced)

## Contrasts for main effect of Sex (Factor A)

cmsex.unbalanced <-  matrix(c(1,-1,1,-1,1,-1), 1)
Sex.const.unbalanced <- glht(growth.lm.unbalanced, linfct = cmsex.unbalanced)
summary(Sex.const.unbalanced,test =Ftest())

## Contrasts for main effect of Bone Density (Factor B)

cmbd.unbalanced <- matrix(c(1,1,-1,-1,0,0,
                 1,1,0,0,-1,-1),  nrow = 2, ncol = 6, byrow =T)
BD.const.unbalanced <- glht(growth.lm.unbalanced, linfct = cmbd.unbalanced)
summary(BD.const.unbalanced,test =Ftest())

## Contrasts for effect of Interaction

cmint.unbalanced <- matrix(c(-1,1,1,-1,0,0,
                 0,0,-1,1,1,-1),  nrow = 2, ncol = 6, byrow =T)
Int.const.unbalanced <- glht(growth.lm.unbalanced, linfct = cmint.unbalanced)
summary(Int.const.unbalanced,test =Ftest())

## Analysis Type I SS

growth.fit.SBD.unbalanced<-aov(Growth~Sex*BD,data = growth.unbalanced)
growth.fit.SBD.summary.unbalanced<-summary(growth.fit.SBD.unbalanced)

growth.fit.SBD.summary.table.unbalanced<-xtable(growth.fit.SBD.summary.unbalanced)
print(growth.fit.SBD.summary.table.unbalanced)

## Type I, II and III

library(car)
growth.fit.type2.unbalanced<-Anova(lm(Growth~Sex*BD,contrasts=list(Sex='contr.sum', BD='contr.sum'),data = growth.unbalanced),type='II')
growth.fit.type2.table.unbalanced<-xtable(growth.fit.type2.unbalanced)
print(growth.fit.type2.table.unbalanced)

growth.fit.type3.unbalanced<-Anova(lm(Growth~Sex*BD,contrasts=list(Sex='contr.sum', BD='contr.sum'),data = growth.unbalanced),type='III')
growth.fit.type3.table.unbalanced<-xtable(growth.fit.type3.unbalanced)
print(growth.fit.type3.table.unbalanced)

## Another way of getting Type III SS

growth.fit3.unbalanced<-aov(Growth~Sex*BD,data = growth.unbalanced,contrasts=list(Sex='contr.sum', BD='contr.sum'))
drop1(growth.fit3.unbalanced,~.,test="F")

## Factor level model: Regression approach

## Creating the indicator variables

growth.regre.data.unbalanced=as.data.frame(model.matrix(~ Sex*BD, growth.unbalanced, contrasts = list(Sex="contr.sum",BD="contr.sum"))[,2:6])
growth.regre.data.unbalanced$Growth=growth.unbalanced$Growth

## Models

growth.full.unbalanced=lm(Growth~Sex1+BD1+BD2+Sex1:BD1+Sex1:BD2, data = growth.regre.data.unbalanced)
growth.NoInt.unbalanced=lm(Growth~Sex1+BD1+BD2, data = growth.regre.data.unbalanced)
growth.NoA.unbalanced=lm(Growth~BD1+BD2+Sex1:BD1+Sex1:BD2, data = growth.regre.data.unbalanced)
growth.NoB.unbalanced=lm(Growth~Sex1+Sex1:BD1+Sex1:BD2, data = growth.regre.data.unbalanced)

anova(growth.full.unbalanced)
anova(growth.NoInt.unbalanced)

## Test for the interaction

anova(growth.NoInt.unbalanced,growth.full.unbalanced)

## Test for factor A

anova(growth.NoA.unbalanced,growth.full.unbalanced)

## Test for factor B

anova(growth.NoB.unbalanced,growth.full.unbalanced)


## Tukey Method

tukey.growth.unbalanced<-TukeyHSD(growth.fit.SBD.unbalanced) 

op0 = par()    # Get current graphical parameters
op1 = op0$mar  # Get current margins in lines
op1[2] = 15     # Modify bottom margins to 20 (lines)
op1 
par(mar = op1)
plot(tukey.growth.unbalanced,las=1)

## Factor level model: Regression approach

## Correlation matrix  of the predictors

library('corrplot')
dummy.cor.unbalanced=cor(growth.regre.data.unbalanced[,1:5])
corrplot(dummy.cor.unbalanced, method = "circle")

## Type I by hand

## Models

## Empty model

growth.unbalanced.empty=lm(Growth~1, data = growth.regre.data.unbalanced)
summary(growth.unbalanced.empty)

## Model for  factor A(Sex)

growth.unbalanced.A=lm(Growth~Sex1, data = growth.regre.data.unbalanced)
summary(growth.unbalanced.A)

## Model for factor A(Sex) and B(BD)

growth.unbalanced.AB=lm(Growth~Sex1+BD1+BD2, data = growth.regre.data.unbalanced)
summary(growth.unbalanced.AB)
anova(growth.unbalanced.AB)

## Full model

growth.unbalanced.full=lm(Growth~Sex1+BD1+BD2+Sex1:BD1+Sex1:BD2, data = growth.regre.data.unbalanced)
summary(growth.unbalanced.full)
anova(growth.unbalanced.full)

## Computing Type I SS

anova(growth.unbalanced.empty,growth.unbalanced.A,growth.unbalanced.AB,growth.unbalanced.full)

## Computing Type I SS for Sex as a contrats proportional to the sample size of the cells

cmsex.typeI.unbalanced <-  matrix(c(1/7,-3/7,3/7,-2/7,3/7,-2/7), 1)
Sex.const.typeI.unbalanced <- glht(growth.lm.unbalanced, linfct = cmsex.typeI.unbalanced)
summary(Sex.const.typeI.unbalanced,test =Ftest())

## Type II by hand

## Model for  factor B(BD)

growth.unbalanced.B=lm(Growth~BD1+BD2, data = growth.regre.data.unbalanced)
summary(growth.unbalanced.B)
anova(growth.unbalanced.B)

## Type II for factor A(Sex)

anova(growth.unbalanced.B,growth.unbalanced.AB,test="F")

## Type II for factor B(BD)

anova(growth.unbalanced.A,growth.unbalanced.AB,test="F")

## Type II for interaction

anova(growth.unbalanced.AB,growth.unbalanced.full,test="F")

## Type III by hand

## Type III for factor A(Sex)

anova(growth.NoA.unbalanced,growth.unbalanced.full,test="F")

## Type III for factor B(BD)

anova(growth.NoB.unbalanced,growth.unbalanced.full,test="F")

## Type III for interaction)

anova(growth.NoInt.unbalanced,growth.unbalanced.full,test="F")

detach(growth.unbalanced)

####################
#################### Synthetic growth hormone: Balanced
####################

## data

growth.balanced<-read.table("CH23TA01-balanced.txt", col.names=c("Sex","BD","Growth"))

attach(growth.balanced)

# Summary tables

freq.balanced<-table(Sex,BD)
mean.freq.balanced<-tapply(Growth, list(Sex,BD), mean)

freq.table.balanced<-xtable(freq.balanced)
print(freq.table.balanced)

mean.freq.table.balanced<-xtable(mean.freq.balanced)
print(mean.freq.table.balanced)

# Interaction plot

interaction.plot(BD,Sex,Growth, xlab="Treatment")

## Factor level model: Regression approach

## Creating the indicator variables

growth.regre.data.balanced=as.data.frame(model.matrix(~ Sex*BD, growth.balanced, contrasts = list(Sex="contr.sum",BD="contr.sum"))[,2:6])
growth.regre.data.balanced$Growth=growth.balanced$Growth

## Correlation matrix  of the predictors

library('corrplot')
dummy.cor.balanced=cor(growth.regre.data.balanced[,1:5])
corrplot(dummy.cor.balanced, method = "circle")

## Models

## Model for  factor A(Sex)

growth.balanced.A=lm(Growth~Sex1, data = growth.regre.data.balanced)
summary(growth.balanced.A)

## Model for  factor A(Sex) and B(BD)

growth.balanced.AB=lm(Growth~Sex1+BD1+BD2, data = growth.regre.data.balanced)
summary(growth.balanced.AB)

## Full model

growth.balanced.full=lm(Growth~Sex1+BD1+BD2+Sex1:BD1+Sex1:BD2, data = growth.regre.data.balanced)
summary(growth.balanced.full)

detach(growth.balanced)


growth.balanced.full=lm(Growth~Sex1+BD1+BD2+Sex1:BD1+Sex1:BD2, data = growth.regre.data.balanced)
growth.balanced.NoInt=lm(Growth~Sex1+BD1+BD2, data = growth.regre.data.balanced)
growth.balanced.NoA=lm(Growth~BD1+BD2+Sex1:BD1+Sex1:BD2, data = growth.regre.data.balanced)
growth.balanced.NoB=lm(Growth~Sex1+Sex1:BD1+Sex1:BD2, data = growth.regre.data.balanced)

anova(growth.balanced.full)
anova(growth.balanced.NoInt)

## Test for the interaction

anova(growth.balanced.NoInt,growth.balanced.full)

## Test for factor A

anova(growth.balanced.NoA,growth.balanced.full)

## Test for factor B

anova(growth.balanced.NoB,growth.balanced.full)


#####################################################################################
#####################################################################################

################# Explaining variability in the ANOVA
x <- seq(0, 8, 0.1)

# Plot an empty chart with tight axis boundaries, and axis lines on bottom and left

plot(x, type="n", xaxs="i", yaxs="i", xlim=c(0, 8), ylim=c(0, 6),
     bty="l", xaxt="n", yaxt="n", xlab="", ylab="",)
axis(2, at=c(3,6), labels=c(expression(bar(Y)[..]),"Y"), las = 2,cex=10.5)
abline(h=3, col="red", lwd=2)
points(7, 5, col = "black", pch=16, cex=1.5)
points(1, 1, col = "black", pch=16, cex=1.5)
#loc<-locator(1)
text(6.977431,5.304798, expression(Y[45]),col="black", cex=1.5)
text(1.002315,0.7, expression(Y[11]),col="black", cex=1.5)

arrows(7, 5,7, 3,lwd=3, col="blue")
arrows(1, 1,1, 3,lwd=3, col="blue")

axis(2, at=4, labels=c(expression(bar(Y)[4.])), las = 2,cex=10.5)
abline(h=4, col="black", lwd=2,  lty=2)
axis(2, at=1.2, labels=c(expression(bar(Y)[1.])), las = 2,cex=10.5)
abline(h=1.4, col="black", lwd=2,  lty=2)

segments(7, 4,7, 3,lwd=3, col="blue")
segments(1, 1.4,1, 3,lwd=3, col="blue")

segments(7, 5,7, 4,lwd=3, col="blue")
segments(1, 1,1, 1.4,lwd=3, col="blue")

par(mfrow = c(1, 3)) 

################# Illustrating the Multiple Comparisons Problem

B <- 2500 # number of simulation replications
N<-500 # number of subjects per group
p<-6 # number of groups
mup<-c(10, rep(0,p-1)) # means of the groups
Group<-rep(1:p,N)[order(rep(1:p,N))]
anova.pvalue <- matrix(0, B, 1)
con.unadj <- 0
con.adj <- 0
sigmaerror<-1
mult.unadj<-matrix(0, p-1, p-1)
mult.adj<-matrix(0, p-1, p-1)


for (i in 1:B) 
{

	set.seed(i)
      res.temp1<-rmvnorm(N, mean =mup, sigma =sigmaerror*diag(p))
	data.simu<-data.frame(resp=vec(res.temp1),group=as.factor(Group))
	anova.fit<-aov(resp~group,data = data.simu)
	anova.fit.summary<-summary(anova.fit)
	anova.pvalue[i]<-anova.fit.summary[[1]]$'Pr(>F)'[1]
	
	if (anova.pvalue[i]<0.05) 
	{
	pairwise.nocorr<-pairwise.t.test(data.simu$resp, data.simu$group, "none")
	mult.unadj.temp<-1*(pairwise.nocorr$p.value[,]<0.05)
	mult.unadj<-mult.unadj+mult.unadj.temp

	con.unadj.temp<-sum(mult.unadj.temp[,2:(p-1)][lower.tri(mult.unadj.temp[,2:(p-1)], diag=F)])
	#print(con.unadj.temp)
	con.unadj<-con.unadj+(con.unadj.temp!=0)*1

	mult.adj.temp <- matrix(NA,p-1,p-1)
    	mult.adj.temp[lower.tri(mult.adj.temp, diag=TRUE)] <- TukeyHSD(anova.fit)$group[,4]
	mult.adj.temp<-1*(mult.adj.temp[,]<0.05)
	mult.adj<-mult.adj+mult.adj.temp

	con.adj.temp<-sum(mult.adj.temp[,2:(p-1)][lower.tri(mult.adj.temp[,2:(p-1)], diag=F)])
	con.adj<-con.adj+(con.adj.temp!=0)*1
	}
} 

mult.unadj/B
mult.adj/B

con.unadj/B
con.adj/B	

################# Multiple Comparisons

# A t-test for each comparison separated, no correction

CN2<-t.test(alcohol$Attractiveness[alcohol$Alcohol=="None"],alcohol$Attractiveness[alcohol$Alcohol=="2 Pints"])
CN4<-t.test(alcohol$Attractiveness[alcohol$Alcohol=="None"],alcohol$Attractiveness[alcohol$Alcohol=="4 Pints"])
V24<-t.test(alcohol$Attractiveness[alcohol$Alcohol=="2 Pints"],alcohol$Attractiveness[alcohol$Alcohol=="4 Pints"])

CN2$estimate[1]-CN2$estimate[2]
CN2$conf.int

CN4$estimate[1]-CN4$estimate[2]
CN4$conf.int

V24$estimate[1]-V24$estimate[2]
V24$conf.int

# Tukey Method

tukey.t<-TukeyHSD(attract.fit) 
plot(tukey.t)

tukey.t.table <- xtable(tukey.t$Alcohol)
print(tukey.t.table)

help(xtable)
methods(xtable)

############### Making the interaction plots of the slides

# First interaction plot for the lecture, no interaction

lev.names <- list(c('M', 'F'), c('4-Pints', '2-Pints', 'None'))
plot.2by2(11,13,18,7,9,14, lev.names)


# Second interaction plot for the lecture, interaction

lev.names <- list(c('M', 'F'), c('4-Pints', '2-Pints', 'None'))
plot.2by2(9,12,18,9,10,14, lev.names)

# Third interaction plot for the lecture, interaction: 2x2 design

lev.names <- list(c('M', 'F'), c('None', 'Alcohol', NA))
plot.2by2(3,7,NA,7,3,NA, lev.names, leg.loc=c(1,5.3))

# Unimportant interaction plot for the lecture, interaction: 2x2 design

lev.names <- list(c('M', 'F'), c('None', 'Alcohol', NA))
plot.2by2(3.5,3,NA,7,7,NA, lev.names, leg.loc=c(1,5.3))

# Interaction plot for the secon seminar

par(mfrow = c(2, 2)) 

lev.names <- list(c('B1', 'B2'), c('A1', 'A2', 'A3'))
plot.2by2(3,3,3,3,3,3, lev.names, leg.loc=c(1.8,2.95),factor.labels=c('Factor B', 'Factor A'), main="Graph A")

lev.names <- list(c('B1', 'B2'), c('A1', 'A2', 'A3'))
plot.2by2(1,3,2,3,5,4, lev.names, leg.loc=c(0.8,5),factor.labels=c('Factor B', 'Factor A'), main="Graph B")

lev.names <- list(c('B1', 'B2'), c('A1', 'A2', 'A3'))
plot.2by2(5,4.6,2.8,1,3,2.8, lev.names, leg.loc=c(1,4),factor.labels=c('Factor B', 'Factor A'), main="Graph C")

lev.names <- list(c('B1', 'B2'), c('A1', 'A2', 'A3'))
plot.2by2(3,2,4,3,2,4, lev.names,factor.labels=c('Factor B', 'Factor A'), leg.loc=c(1,4), main="Graph D")











