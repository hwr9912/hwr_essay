# 动脉瘤LASSO模型代码
## 1 数据处理 ----
rm(list = ls())
library(readxl)
windowsFonts(A=windowsFont("Arial"))
sourcedata <- read_excel("training.xlsx")
str(sourcedata)
sourcedata <- na.omit(sourcedata)

sourcedata$sex <- factor(sourcedata$sex, 
                         levels = c(0,1), labels = c("female", "male"))
sourcedata$drinking <- factor(sourcedata$drinking,
                              levels = c(0,1), labels = c("good", "poor"))
sourcedata$smoking <- factor(sourcedata$smoking,
                             levels = c(0,1), labels = c("good", "poor"))
# sourcedata$hypertension <- factor(sourcedata$hypertension,
#                                   levels = c(0,1), labels = c("good", "poor"))
sourcedata$diabetes <- factor(sourcedata$diabetes, 
                              levels = c(0,1), labels = c("good", "poor"))
sourcedata$hyperlipidemia <- factor(sourcedata$hyperlipidemia, 
                                    levels = c(0,1), labels = c("good", "poor"))
sourcedata$CHD <- factor(sourcedata$CHD, 
                         levels = c(0,1), labels = c("good", "poor"))
# sourcedata$location <- factor(sourcedata$location, 
#                               levels = c(0,1,2), labels = c("left", "right", "bilateral"))
sourcedata$irrugular <- factor(sourcedata$irrugular, 
                               levels = c(0,1), labels = c("regular", "irregular"))
# sourcedata$rupture <- factor(sourcedata$rupture,
#                              levels = c(0,1), labels = c("unruptured", "ruptured"))

cluster0 <- na.omit(sourcedata)
cluster0$OSI <- cluster0$OSI*10
#数据纳入
x <- model.matrix(rupture~., cluster0[1:23])[,-1]
#x1 <- model.matrix(rupture~., cluster0)
y <- cluster0$rupture

## 2 训练模型 ----
library(glmnet)
library(ggplot2)
#k <- 10
grid<-10^seq(1,-4,length = 400)
#设置对数梯度以便挑选λ
set.seed(1)
cv.fit<-glmnet(x, y, family = "binomial",
               alpha = 1, lambda = grid)
#以上的lambda即LASSO中的λ，数值越大惩罚力度越强。
#LASSO回归的α=1，=2时进化为岭回归

## cross-validation
#set.seed(k)
cv.out_class<-cv.glmnet(x, y, type.measure = "class",
                        family='binomial',
                        lambda = grid, nfolds=10)
cv.out_auc<-cv.glmnet(x, y, type.measure = "class",
                      family='binomial',
                      lambda = grid, nfolds=10)

##组成广义线性模型
l.coef1<-coef(cv.out_class$glmnet.fit,s=cv.out_class$lambda.min,exact = F)
l.coef2<-coef(cv.out_class$glmnet.fit,s=cv.out_class$lambda.1se,exact = F)
l.coef1
l.coef2
#根据挑选的变量建立模型
mod<-glm(rupture~OSI+AR+hypertension+WSS,
         family="binomial",data = cluster0)
summary(mod)
#求入选变量的OR和95%CI
exp(confint(mod))
exp(coef(mod))

## 作图 AR\hypertension\OSI\WSS
plot(cv.out_class, xlab = "Log(λ)", ylab = "Misclassification Error")

plot(cv.out_auc, xlab = "Log(λ)", ylab = "Area Under ROC")

plot(cv.fit, xvar = "lambda", label = TRUE)
abline(v = log10(cv.out_class$lambda.1se), lty = 2)
abline(v = log10(cv.out_class$lambda.min), lty = 2)

## 3 作出nomogram ----
library(foreign)
library(survival)
library(rms)

cluster0$OSI <- cluster0$OSI/10
ia_sort <- datadist(cluster0)
options(datadist="ia_sort")

formula1<-as.formula(rupture ~ OSI + AR + hypertension + WSS)
f<-lrm(formula1,data=sourcedata,x=T,y=T)

##构建列线图
#nom1 <- nomogram(f, fun = function(x) surv(0, x), lp=F,
#                 funlabel = "rupture probability")
#plot(nom1, xfrac = 0.15)
nom1<-nomogram(f,
               fun=function(x)1/(1+exp(-x)),
               lp=F,
               fun.at = c(0.1,0.3,0.5,0.7,0.9),
               funlabel = "Risk")
plot(nom1)

##内部验证
validate(f, method="boot", B=1000, dxy=T)
#rcorrcens可计算C-index
rcorrcens(rupture ~ predict(f), data = sourcedata)

##建立校准曲线并绘图
cal1<-calibrate(f,method = "boot",B=1000)
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Nomogram Predicted Rupture", ylab = "Actual Rupture")
## 4 ROC曲线作图 ----
library(pROC)
library(ggplot2)

names(cluster0)

# OSI+AR+LSA+hypertension+WSS
roc_OSI <- roc(cluster0$rupture, cluster0$OSI, ci=T, auc = T)
roc_AR <- roc(cluster0$rupture, cluster0$AR)
roc_hyp <- roc(cluster0$rupture, cluster0$hypertension)
roc_WSS <- roc(cluster0$rupture, cluster0$WSS)

set.seed(1)
train<-sample(1:nrow(x),nrow(x)/2)
pr <- predict(f, data = cluster0[train,])
roc_training<- roc(cluster0$rupture~pr)

g2 <- ggroc(list(OSI=roc_OSI, AR=roc_AR, hypertension=roc_hyp,
                  WSS=roc_WSS, Nomogram=roc_training),
            size = 1)
g2+
  annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1, linetype = 1, size = 1)+
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.1),
        legend.justification = c(0.9,0.1),
        legend.key = element_blank(),
        axis.title=element_text(size=12, family = "A"),
        axis.text=element_text(size=12, family = "A"))+
  guides(color=guide_legend(title=NULL))+
  scale_colour_discrete(breaks = c("OSI", "AR", "hypertension", "WSS", "Nomogram"),
                        labels = c("OSI                  AUC=0.71", 
                                   "AR                   AUC=0.77", 
                                   "hypertension AUC=0.58",
                                   "WSS                AUC=0.77", 
                                   "Nomogram    AUC=0.87"))

## 5 作出DCA曲线与CIC曲线 ----
library(rmda)
dca1<- decision_curve(rupture ~ AR,data = sourcedata, 
                      family = binomial(link ='logit'),#模型类型，这里是二分类
                      thresholds= seq(0,1, by = 0.01),
                      confidence.intervals =0.95,#95可信区间
                      study.design = 'cohort')#研究类型，这里是队列研究
dca4<- decision_curve(rupture ~ OSI,data = sourcedata, 
                      family = binomial(link ='logit'),
                      thresholds= seq(0,1, by = 0.01),
                      confidence.intervals =0.95,
                      study.design = 'cohort')
All <- decision_curve(rupture ~ OSI + AR + hypertension + WSS,data = sourcedata, 
                      family = binomial(link ='logit'),
                      thresholds= seq(0,1, by = 0.01),
                      confidence.intervals =0.95,
                      study.design = 'cohort')

List <- list(dca1,dca4,All)
plot_decision_curve(List,curve.names= c("AR","OSI","Nomogram"),
                    cost.benefit.axis =FALSE,col = c('red','blue','green'),
                    confidence.intervals =FALSE,standardize = FALSE)
plot_clinical_impact(All,curve.names= c("Nomogram"))
## 6 Hosmer-Lemeshow拟合优度检验及其他数据 ----
# install.packages("ResourceSelection")
library(ResourceSelection)
hoslem.test(mod$y, fitted(mod), g=10)

## 7 外部验证 ----
# 读取数据集
certificate <- read_excel("certification.xlsx")

### 7.1 校准曲线作图 ----
fc<-lrm(formula1,data=certificate,x=T,y=T)
# validate(f, method="boot", B=1000, dxy=T)
#rcorrcens可计算C-index
rcorrcens(rupture ~ predict(fc), data = certificate)
cal1<-calibrate(fc,method = "boot",B=1000)
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Nomogram Predicted Rupture", ylab = "Actual Rupture")

### 7.2 ROC曲线作图 ----
library(pROC)
library(ggplot2)

set.seed(1)
pr <- predict(fc, data = certificate)
roc_certi<- roc(certificate$rupture~pr)
g3 <- ggroc(list(Nomogram=roc_certi),
            size = 1)
# ggroc(rocobj, alpha = 0.5, colour = "red", linetype = 2, size = 2)
g3+
  geom_segment(x = -1, y = 0, xend = 0, yend = 1, linetype = 1, colour = "black",size = 1)+
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.9,0.1),
        legend.justification = c(0.85,0.25),
        legend.key = element_blank(),
        axis.title=element_text(size=12, family = "A"),
        axis.text=element_text(size=12, family = "A"))+
  guides(color=guide_legend(title=NULL))+
  scale_colour_discrete(breaks = c("Nomogram"), labels = c("Nomogram  AUC=0.87"))
  
### 7.3 Hosmer-Lemeshow拟合优度检验及其他数据 ----
library(ResourceSelection)
hoslem.test(fc$y, fitted(fc), g=10)
