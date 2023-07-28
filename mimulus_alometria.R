####Associated to Gangas, Reyes and Carvallo Manuscript##
####sent to Plant Species Biology - 2023

R.Version()

#dataset
alometria<-read.table(file="Mimulus-alometria.csv", header=T, sep=",", check.names=T, dec=".") #, row.names=1)
alometria$species_abrev<-as.factor(alometria$species_abrev)
alometria$location<-as.factor(alometria$location)

# Descriptive statistics
library(psych)
names(alometria)
describeBy(alometria, c("species_abrev","location"), skew=F) # descriptive statistical by species

#######UNIVARIATE ANALYSIS########
####Differences among taxa########
#m: models, change traits
m0<-glm(pedicele~species_abrev,family=gaussian(),data=alometria)
summary(m0)
dropterm(m0,test="F") 
library(multcomp)
summary(glht(m0, mcp(species_abrev="Tukey")))

######MULTIVARIATE ANALYSYS########
####Difference among traits#######
#Multi-dimensional scaling (MDS)
library(vegan)
#1. First, determine the subset of variables used in analyses
names(alometria)
rasgos<-alometria[,c(6:11,13)]
names(rasgos)
#2. distance among taxa
d_Bray<-vegdist(rasgos, method="euclidean",na.rm=T)
#3. MDS, using two dimenssions
mds<-cmdscale(d_Bray, k=2)
#4. Adding the factor (taxa names),1 agregamos la info del factor a cada score (especie/localidad)
taxa<-factor(alometria[,3]) # c(3:20,22:27)
#5. Assessing the statistical differences among taxa 
an<-anosim(d_Bray,taxa, permutations=999,distance="euclidean")
summary(an)
adonis(d_Bray~taxa,method="euclidean")
#6. Estimating mean distance among taxa
meandist(d_Bray, taxa)

#######
#Principal components analysis (PCA)
######
#####incorporate variables according the scenarios described in the text
names(alometria)
datos.pca<-as.data.frame(alometria[,c(1,3,6,7,8,13)]) 
#####incorporate variables names according last funtion
colnames(datos.pca) <- c("ID","Taxa", "Corolla throat width", 
                         "Corolla tube lenght","Corolla throat height",
                         "Corolla area")#"Corolla tube lenght", #"Corolla tube length",
####preparing data to PCA
row.names(datos.pca)<-paste(datos.pca$ID)
datos.pca$Taxa<-NULL
datos.pca$ID<-NULL
head(datos.pca)
datos.pca
######pca analysis
pca.alometria<-prcomp(na.omit(datos.pca), center = TRUE,scale. = TRUE)
summary(pca.alometria) #PCA Result
pca.alometria

#######
#contribution of each variable in the PCA (%)
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}

loadings <- pca.alometria$rotation
sdev <- pca.alometria$sdev
var.coord <- t(apply(loadings, MARGIN = 1, var_coord_func, sdev)) 
var.cos2 <- var.coord^2
comp.cos2 <- apply(var.cos2, MARGIN = 2, FUN = sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2, MARGIN = 1, contrib, comp.cos2))
(var.contrib)

#PCA plot
df_out <- as.data.frame(pca.alometria$x)
df_out
datos.pca<-datos.pca[complete.cases(datos.pca),]
datos.pca
df_out$group<-sapply(strsplit(as.character(row.names(datos.pca)), "_"), "[[", 1)
df_out$observation <- 1:nrow(df_out) 
head(df_out)
df_out[c(1:27),6]<-"E. depressa" # last value (5) depict the column with the grouping factor (taxa)
df_out[c(28:42),6]<-"E. glabrata"
df_out[c(43:60),6]<-"Hybrids"
df_out[c(61:110),6]<-"E. lutea"
names(df_out)
df_out

library(grid)
library(gridExtra)
library(factoextra)
g<-fviz_pca_biplot(pca.alometria, geom=c("point"),pointsize=5,
                   labelsize=10,
                   col.ind = df_out$observation, geom.ind="point",
                   addEllipses = TRUE, label = "var",
                   fill.ind = df_out$observation,
                   col.var = "black", repel = TRUE,
                   legend.title = "Species",alpha.var=1,
                   title="")   +
  scale_color_manual(values=c("blueviolet","red", "darkgreen","coral4"))+
  scale_fill_manual(values=c("blueviolet","red", "darkgreen","coral4"))+
  scale_shape_manual(values=c(24,22,21,23))+
  theme(text = element_text(size = 40),
        axis.title = element_text(size = 40),
        axis.text = element_text(size = 40))

g$layers[[2]]$aes_params$size <- 4 #número de variables usadas(?)
g$layers[[2]]$geom_params$arrow$length <- unit(35, units = "points")
g

###################################
####To regress PCA with ovule number
df_out #show PCA scores, previously obtained
alometria<-alometria[order(alometria$ID), ]
rm.alometria<-alometria[c(-26,-31,-36,-96,-97,-99),] #deleting gla012, gla017,  lut070,lut073
head(rm.alometria)
head(df_out)
names(rm.alometria)
df_out<-cbind(df_out, ID2=rm.alometria$ID, species=rm.alometria$species,
              species_abrev=rm.alometria$species_abrev,ovules=rm.alometria$ovule_number,
              ovules_per_ovaryvolume=rm.alometria$ovules_per_ovaryvolume)
names(df_out)
head(df_out)

#linear_regression
library(lme4)
fit<-lm(ovules ~ PC1, data=df_out)
fits<-lmList(ovules ~ PC1 | species_abrev, data=df_out)
summary(fit)
summary(fits)

######plot regression by group 
#(not included in the manuscript)
library(ggplot2)
ggplot(df_out, aes(x = PC1, y = ovules, color = species_abrev)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "PC1 scores",
       y = "Ovule number") +
  theme_minimal()

#####################
#Allometric models###
#####################
# pedicel vs tube length
library(smatr)
names(alometria)
m1.1<-(sma(pedicele ~ corolla_tube_length+species_abrev,
           data=alometria,log="xy",slope.test = 0))
summary(m1.1)
plot(m1.1)
plot(m1.1, which= "residual" )
plot(m1.1, which= "qq" )

# corolla_width vs corolla_tube_length
names(alometria)
m1.2<-sma(corolla_width ~ corolla_tube_length+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m1.2)
plot(m1.2)
plot(m1.2, which= "residual" )
plot(m1.2, which= "qq" )

# corola_throat_height vs corolla_tube_length
names(alometria)
m1.3<-sma(corolla_height_lateral ~ corolla_tube_length+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m1.3)
plot(m1.3)
plot(m1.3, which= "residual" )
plot(m1.3, which= "qq" )

# corolla area vs corolla_tube_length
names(alometria)
m1.4<-sma(corolla_area ~ corolla_tube_length+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m1.4)
plot(m1.4)
plot(m1.4, which= "residual" )
plot(m1.4, which= "qq" )

# style length vs corolla_tube_length
names(alometria)
m1.5<-sma(style_lenght ~ corolla_tube_length+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m1.5)
plot(m1.5)
plot(m1.5, which= "residual" )
plot(m1.5, which= "qq" )

# ovary volume vs corolla_tube_length
names(alometria)
m1.6<-sma(ovary_volume ~ corolla_tube_length+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m1.6)
plot(m1.6)
plot(m1.6, which= "residual" )
plot(m1.6, which= "qq")

# pedicel vs ovule number
names(alometria)
m2.1<-sma(ovule_number ~ pedicele+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m2.1)
plot(m2.1)
plot(m2.1, which= "residual" )
plot(m2.1, which= "qq")

# corolla throat width vs ovule number
names(alometria)
m2.2<-sma(ovule_number ~ corolla_width+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m2.2)
plot(m2.2)
plot(m2.2, which= "residual" )
plot(m2.2, which= "qq")

# corolla throat heigth vs ovule number
names(alometria)
m2.3<-sma(ovule_number ~ corolla_height_lateral+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m2.3)
plot(m2.3)
plot(m2.3, which= "residual" )
plot(m2.3, which= "qq")

# corolla length vs ovule number
names(alometria)
m2.4<-sma(ovule_number ~ corolla_tube_length+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m2.4)
plot(m2.4)
plot(m2.4, which= "residual" )
plot(m2.4, which= "qq")

# corolla area vs ovule number
names(alometria)
m2.5<-sma(ovule_number ~ corolla_area+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m2.5)
plot(m2.5)
plot(m2.5, which= "residual" )
plot(m2.5, which= "qq")

# style lenght vs ovule number
names(alometria)
m2.6<-sma(ovule_number ~ style_lenght+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m2.6)
plot(m2.6)
plot(m2.6, which= "residual" )
plot(m2.6, which= "qq")

# ovary_volume vs ovule number
names(alometria)
m2.7<-sma(ovule_number ~ ovary_volume+species_abrev,
          data=alometria,log="xy",slope.test = 1)
summary(m2.7)
plot(m2.7)
plot(m2.7, which= "residual" )
plot(m2.7, which= "qq")

# Formatted plots for the manuscript
m2.7
x_expression <- expression("Ovary volume (log"[10]*")")
y_expression <- expression("Ovule number (log"[10]*")")
library("scales")
par(mgp=c(0.25,2.5,0),mar = c(12.5,14,1,1)+0.1)
plot(m2.7, xlab="", ylab="",
     col=alpha(c("blueviolet","coral4","red", "darkgreen"),.7),cex=5,
     cex.axis=5., cex.lab=5.,pch=c(17,18,15,16),las=1,lwd=5) # hybrids,lutea, glabrata, depressa
abline(a=0.8,b=1,col="gray",lwd=3,lty=1)
mtext(text=y_expression,side=2,line=9,cex=5)
mtext(text = x_expression,
      side = 1,
      line = 9,cex=5)
