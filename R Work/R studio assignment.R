##John Salber

#This is importing my data from python work
link='https://github.com/jsalber/542/raw/main/First%20Assignment/Final%20Data/FinalData.csv'
myFile=url(link)
fromPy=read.csv(file = myFile)
row.names(fromPy)=NULL
str(fromPy)

# I only used the per capita data
selection=c("Country","PerCapitaEmissions", "PerCapitaGDP.Thousands.","KilowattsPerCapita")


dataToCluster=fromPy[,selection]
row.names(dataToCluster)=dataToCluster$Country
dataToCluster$Country=NULL
boxplot(dataToCluster,horizontal = T, las=2,cex.axis=2)
as.data.frame(scale(dataToCluster))
log(dataToCluster)
library(factoextra)

#finding out how many clusters to use. I ended up using 4 because it was in the middle
#between 9 and 2
set.seed(999)
library(cluster)
dataToCluster_DM=daisy(x=dataToCluster, metric = "gower")
fviz_nbclust(dataToCluster, 
             pam,
             diss=dataToCluster_DM,
             method = "gap_stat",
             k.max = 10,verbose = F)
fviz_nbclust(dataToCluster, 
             hcut,
             diss=dataToCluster_DM,
             method = "gap_stat",
             k.max = 10,
             verbose = F,
             hc_func = "agnes")
fviz_nbclust(dataToCluster, 
             hcut,
             diss=dataToCluster_DM,
             method = "gap_stat",
             k.max = 10,
             verbose = F,
             hc_func = "diana")


NumberOfClusterDesired=4

# Partitioning technique
res.pam = pam(x=dataToCluster_DM,
              k = NumberOfClusterDesired,
              cluster.only = F)

# Hierarchical technique- agglomerative approach

#library(factoextra)
res.agnes= hcut(dataToCluster_DM, 
                k = NumberOfClusterDesired,
                isdiss=TRUE,
                hc_func='agnes',
                hc_method = "ward.D2")

# Hierarchical technique- divisive approach
res.diana= hcut(dataToCluster_DM, 
                k = NumberOfClusterDesired,
                isdiss=TRUE,
                hc_func='diana',
                hc_method = "ward.D2")


fromPy$pam=as.factor(res.pam$clustering)
fromPy$agn=as.factor(res.agnes$cluster)
fromPy$dia=as.factor(res.diana$cluster)

aggregate(data=fromPy,
         HDI~pam,
          FUN=mean)

aggregate(data=fromPy,
          HDI~agn,
          FUN=mean)

aggregate(data=fromPy,
          HDI~dia,
          FUN=mean)

library(dplyr)
##Picking the right cluster technique

fromPy$pam=dplyr::recode_factor(fromPy$pam, 
                                `1` = '4',`2`='3',`3`='2',`4`='1')
fromPy$agn=dplyr::recode_factor(fromPy$agn, 
                                `1` = '4',`2`='3',`3`='2',`4`='1')
fromPy$dia=dplyr::recode_factor(fromPy$dia, 
                                `1` = '4',`2`='3',`3`='2',`4`='1')
fviz_silhouette(res.pam)

fviz_silhouette(res.agnes)

library(factoextra)
fviz_silhouette(res.diana)

head(data.frame(res.pam$silinfo$widths),10)

pamEval=data.frame(res.pam$silinfo$widths)
agnEval=data.frame(res.agnes$silinfo$widths)
diaEval=data.frame(res.diana$silinfo$widths)

pamPoor=rownames(pamEval[pamEval$sil_width<0,])
agnPoor=rownames(agnEval[agnEval$sil_width<0,])
diaPoor=rownames(diaEval[diaEval$sil_width<0,])

library("qpcR") 


bap_Clus=as.data.frame(qpcR:::cbind.na(sort(pamPoor), sort(agnPoor),sort(diaPoor)))
names(bap_Clus)=c("pam","agn","dia")
bap_Clus

projectedData = cmdscale(dataToCluster_DM, k=2)
#Going with Agnes because its has the least negatives and the highest silhouette width


fromPy$dim1 = projectedData[,1]
fromPy$dim2 = projectedData[,2]



fromPy[,c('dim1','dim2')][1:10,]

base= ggplot(data=fromPy,
             aes(x=dim1, y=dim2,
                 label=Country)) 
base + geom_text(size=2)

pamPlot=base + labs(title = "PAM") + geom_point(size=2,
                                                aes(color=pam),
                                                show.legend = T) 

diaPlot=base + labs(title = "DIANA") + geom_point(size=2,
                                                  aes(color=dia),
                                                  show.legend = T) 
                                                  
agnPlot=base + labs(title = "AGNES") + geom_point(size=2,
                                                  aes(color=agn),
                                                  show.legend = T)
library(ggpubr)

ggarrange(pamPlot, agnPlot, diaPlot,ncol = 3,common.legend = T)

fviz_dend(res.agnes,k=NumberOfClusterDesired, cex = 0.45, horiz = T,main = "AGNES approach")

selection=c("Country","PerCapitaEmissions", "PerCapitaGDP.Thousands.","KilowattsPerCapita")

dataForFA=fromPy[,selection]

names(dataForFA)

library(lavaan)

#this is my model data. It is all in per capita so that the numbers work together
model='energymoney=~PerCapitaEmissions + PerCapitaGDP.Thousands.+KilowattsPerCapita'

fit<-cfa(model, data = dataForFA,std.lv=TRUE)
indexCFA=lavPredict(fit)

indexCFA[1:10]


library(scales)
indexCFANorm=rescale(as.vector(indexCFA), 
                     to = c(0, 10))
indexCFANorm[1:10]

fromPy$demo_FA=indexCFANorm

base=ggplot(data=fromPy,
            aes(x=demo_FA,y=HDI))
base+geom_point()

evalCFA1=parameterEstimates(fit, standardized =TRUE)

evalCFA1[evalCFA1$op=="=~",c('rhs','std.all','pvalue')]

evalCFA2=as.list(fitMeasures(fit))

evalCFA2[c("chisq", "df", "pvalue")] 

evalCFA2$tli

evalCFA2[c( 'rmsea.ci.lower','rmsea','rmsea.ci.upper')] 

library(semPlot)

semPaths(fit, what='std', nCharNodes=0, sizeMan=12,
         edge.label.cex=1.5, fade=T,residuals = F)

#picking the right model for my data
hypo1=formula(HDI~ KilowattsPerCapita)


hypo2=formula(HDI~ PerCapitaEmissions+KilowattsPerCapita+PerCapitaGDP.Thousands. )

gauss1=glm(hypo1,
           data = fromPy,
           family = 'gaussian')

gauss2=glm(hypo2,
           data = fromPy,
           family = 'gaussian')

summary(gauss1)

summary(gauss2)

#this shows that the second model is better
anova(gauss1,gauss2,test="Chisq")



library(rsq)
rsq(gauss2,adj=T)


plot(gauss2,1)

plot(gauss2,2)

#p-value is much below 0.05 so data is not normal
shapiro.test(gauss2$residuals)


plot(gauss2, 3)

library(lmtest)

#P-value is above 0.05 so we can assume homosedasticity
bptest(gauss2) 

library(car)

vif(gauss2)

plot(gauss2,5)

gaussInf=as.data.frame(influence.measures(gauss2)$is.inf)
gaussInf[gaussInf$cook.d,]

library(sjPlot)
plot_models(gauss2,vline.color = "grey")

library(caret)

set.seed(123)

selection = createDataPartition(fromPy$HDI,
                                p = 0.75,
                                list = FALSE)

trainGauss = fromPy[ selection, ]
#
testGauss  = fromPy[-selection, ]


ctrl = trainControl(method = 'cv',number = 5)


gauss2CV = train(hypo2,
                 data = trainGauss, 
                 method = 'glm',
                 trControl = ctrl)
##there is a moderate amount of correlation
gauss2CV


# this tells us that its an ok predictor since r^2 is above 0.5
#the data isnt a perfect perdictor but its better than a random guess
predictedVal<-predict(gauss2CV,testGauss)

postResample(obs = testGauss$HDI,
             pred=predictedVal)


