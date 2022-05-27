R.Version()
rm(list=ls()) # Borra todos los objetos
setwd("/Users/user/Downloads/mimulus_Francisca") #ESCOGE EL DIRECTORIO
dir()

# datos obtenidos desde el archivo Drive. Guardar como .xls y posterior transformación a 
# archivo .csv. Eliminar primera fila en excel
mimulus<-read.table(file="Morfologia_Mimulus_2021.csv", header=T, sep=",", check.names=T, dec=".") #, row.names=1)
mimulus
alometria<-read.table(file="Mimulus-alometria.csv", header=T, sep=",", check.names=T, dec=".") #, row.names=1)
alometria

names(alometria)
alometria<-alometria[c(1:117),c(1:20)]
#datos[fila,columnas] #formato de lectura de los datos
summary(as.factor(alometria$species_abrev))
summary(as.factor(alometria$species))
alometria$species_abrev

#incluir sólo luteus
#mimulus_lut<-as.data.frame(mimulus[mimulus$especie %in% c('LUT', 
                                                    #"VAR"), ])
#summary(as.factor(mimulus_lut$especie))
#summary(as.factor(mimulus_lut$poblacion))

library(psych)
library(DT)
?describeBy
describeBy(alometria, "species_abrev", skew=F) # para describir las poblaciones por el factor "poblacion"
# de la función anterior se obtienen los tamaños muestreales (N) para cada población

#mimulus_lut$latitud
#mimulus_lut$latitud<-(mimulus_lut$latitud)*-1

# Graficos de frecuencia de los datos
library(ggplot2)
library(ggridges)
#install.packages("stringi",type="mac.binary")
#library(stringi)
library(tidyverse)

# Grafico de distribución de variables. Modificar el nombre de la variable
  #en x=_____ para ir cambiándola. El nombre la variable se cambia en 
  # xlab=______
summary(as.factor(alometria$species_abrev))
summary(as.factor(alometria$species))
names(alometria)
alometria$species_abrev

mimulus_meaned <- alometria %>%
  group_by(species_abrev) %>% 
  summarise_all(funs(mean), na.rm=T)
mimulus_meaned

names(alometria)
alometria$species_abrev

alometria %>%
  mutate(species = fct_rev(fct_relevel(species, 
                            "M. depressus", "Hybrid", "M. luteus",
                            "M. glabratus"))) %>%
ggplot(aes(x = pedicele, y = species)) +
  theme(legend.position = "None", 
        axis.title.x=element_text(hjust = 0))+
  geom_density_ridges(aes(fill = species),jittered_points = TRUE, scale = .95,
                      rel_min_height = .01,
                      point_shape = "|", point_size = 2, size = 0.5,
                      position = position_points_jitter(height = 0),
                      quantile_lines=TRUE,
                      quantile_fun=mean, show_guides=F)+ #para sacar legenda show_guide=F
  xlab("Pedicele length [mm]")+
  ylab("")+
   #se debe usar para aquellas variables que no pueden tomar valores negativos
  scale_fill_viridis_c(option = "turbo",direction=1,trans="reverse") + #colores para variables discretas
  coord_cartesian(clip = "off")+
  theme_ridges(font_size = 13, grid = TRUE)+
  scale_fill_manual(values = c("darkblue","yellow","darkgreen","red"))
  #geom_text(x=1.1, y=1.5, label="N = 26")+ #rio tepu
  #xlim(-0.4,1.2)

?scale_fill_viridis_c

# Modelos
  #Los modelos son de la forma y = x donde y es la variable a estudiar y x es
  #la latitud. Por lo tanto, se está evaluando la variación de cada rasgo en
  #el gradiente latitudinal.
  #Las familias (family) debe modificarse de acuerdo a la estructura de los
  #datos que estamos estudiando: gaussian() para datos continuos (la mayoría)
  #y poisson() para datos discretos (ejemplo, el número de manchas)

names(mimulus_lut)
mimulus_lut<-cbind(mimulus_lut,lat_sd=((mimulus_lut$lat)/41.2329))
mimulus_lut$lat_sd  
mimulus_lut<-cbind(mimulus_lut,alt_sd=((mimulus_lut$alt)/3411))
mimulus_lut<-cbind(mimulus_lut,alt_lat_sd=(mimulus_lut$alt_sd*mimulus_lut$lat_sd))

#diferencias entre grupos (TABLA 1)
names(alometria)
alometria$species_abrev <- factor(alometria$species_abrev)
m0<-glm(ovule_number~species_abrev,family=poisson(),data=alometria) # formula, familia y datos
summary(m0)
library(MASS)
dropterm(m0,test="F") #el menor valor de AIC es el mejor modelo
library(multcomp)
summary(glht(m0, mcp(species_abrev="Tukey")))
?glht()

names(alometria)
my_data <- alometria[, c(3,8:16)]
my_data

library(dplyr)
library(tidyr)
library(purrr)

#correlación
library("PerformanceAnalytics")
names(my_data)
my_data$species_abrev
(chart.Correlation(my_data[c(1:28),-1], histogram=TRUE, pch=19))
colnames(my_data)[1]<-"taxa"
colnames(my_data)[2]<-"PL"
colnames(my_data)[3]<-"CW"
colnames(my_data)[4]<-"CTL"
colnames(my_data)[5]<-"CTH"
colnames(my_data)[6]<-"OL"
colnames(my_data)[7]<-"OW"
colnames(my_data)[8]<-"SL"
colnames(my_data)[9]<-"ON"
colnames(my_data)[10]<-"CA"
library(corrplot)
cor.mtest <- function(mat, ...)
{
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) 
  {
    for (j in (i + 1):n)
    {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
names(my_data)
my_data
p.mat <- cor.mtest(my_data[c(65:117),-1])
M<-cor(my_data[c(65:117),-1],method = "pearson", use = "pairwise.complete.obs")
corrplot(M,type="upper",order="alphabet",method = "circle",
         p.mat=p.mat, sig.level=0.05, insig="blank",
         diag=F, addCoef.col ='black', number.cex = 0.8)
?corrplot

#Relaciones alométricas
# estos son utilizados para construir los gráficos
library(smatr)
names(alometria)
model<-(sma(ovule_number ~ ovary_volume+species_abrev,
         data=alometria,log="xy",slope.test = 1))
x_expression <- expression("Ovary volume (log"[10]*")")
y_expression <- expression("Ovule number (log"[10]*")")
library("scales")
par(mgp=c(0,2.5,0),mar = c(13,14,5,5)+0.1)
plot(model, xlab="", ylab="",
     col=alpha(c("blueviolet","coral4","red", "darkgreen"),.7),cex=5,
     cex.axis=4, cex.lab=4,pch=c(17,18,15,16),las=1,lwd=5) # hybrids,lutea, glabrata, depressa
abline(a=0.7,b=1,col="gray",lwd=3,lty=1)
mtext(text = x_expression,
      side = 1,
      line = 9,cex=4)
mtext(text=y_expression,side=2,line=9,cex=4)


# pedicelo vs tube length
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


#Escalamiento multidimensional (MDS)
#1. permite hacer comparaci?n entre distintos grupos. Escogeremos los datos a analizar
names(rasgos)
#2. se debe estimar una medida de distancia entre los grupos, en este caso utilizo "Bray-Curtis"
library(vegan)
rasgos
names(alometria)
rasgos<-alometria[,c(9:11,16)]
summary(rasgos)
d_Bray<-vegdist(rasgos, method="euclidean",na.rm=T)
#3. se realiza el MDS, se define a priori el n?mero de dimensiones a la cual se reducen nuestras variables (en nuestro caso, cada especie es una variable)
mds<-cmdscale(d_Bray, k=2)
#3,1 agregamos la info del factor a cada score (especie/localidad)
mds1<-as.data.frame(mds)
mds1
#names(alometria)
#mds1<-cbind(alometria[,3],mds1)
#mds1
#4. determinamos un factor de agrupamiento, en nuestro ejemplo, el tipo de sitio
names(alometria)
taxa<-factor(alometria[,3]) # c(3:20,22:27)
taxa

#5. Finalmente evaluamos estad?sticamente si las distancias observadas entre grupos son significativas 
an<-anosim(d_Bray,taxa, permutations=999,distance="euclidean")
summary(an)
#?adonis
adonis(d_Bray~taxa,method="euclidean")

# para estimar las distancias medias
meandist(d_Bray, taxa)

# principal components
names(alometria)
head(alometria)
datos.pca<-as.data.frame(alometria[,c(3,9:11,16)]) #no olvidar agregar variable que separa grupos, aunque en la linea 249 se suprime
head(datos.pca)
names(datos.pca)
colnames(datos.pca) <- c("Taxa", "Corolla throat width", 
                         "Corolla tube lenght","Corolla throat height",
                                                    "Corolla area")#"Corolla tube lenght", #"Corolla tube length",
row.names(datos.pca)<-paste(datos.pca$Taxa,row.names(datos.pca),sep="_")
datos.pca$Taxa<-NULL
head(datos.pca)
datos.pca

pca.alometria<-prcomp(na.omit(datos.pca), center = TRUE,scale. = TRUE)
summary(pca.alometria) #RESULTADO
pca.alometria

#contribución de cada variable a PCA en % (porcentaje)
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
var.contrib


#PCA plot
df_out <- as.data.frame(pca.alometria$x)
df_out
datos.pca<-datos.pca[complete.cases(datos.pca),]
datos.pca
df_out$group<-sapply(strsplit(as.character(row.names(datos.pca)), "_"), "[[", 1)
df_out$observation <- 1:nrow(df_out) 
df_out
df_out[c(1:27),6]<-"E. depressa" # el último valor corresponde a la columna donde se agrega el factor de agrupamiento
df_out[c(28:45),6]<-"Hybrids"
df_out[c(46:61),6]<-"E. glabrata"
df_out[c(62:111),6]<-"E. lutea"
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

#relacionando el nº ovulos con el PCA
names(df_out)
names(alometria)
alometria
df_out$Ovules<-alometria[c(1:116),15]
df_out$OvulesperPCA1<-abs(df_out$Ovules/(df_out$PC1)) #se ha puesto en positivo
df_out

#plot
names(df_out)
summary(df_out)
df_out$group<-as.factor(df_out$group)
library(dplyr)

df_out
names(alometria)
df.summary<-alometria%>%
  group_by(species_abrev)%>%
  summarise(sd = sd(ovules_ovaryvolume, na.rm = T),
    mean = mean(ovules_ovaryvolume,na.rm=T))
class(df.summary)

ggplot(df.summary, aes(species_abrev, mean)) +
  geom_jitter(position = position_jitter(1), color = "darkgray") + 
  geom_pointrange(aes(ymin = mean-sd, ymax = mean+sd),data = df.summary)

#modelos
names(df_out)
m1.1<-glm(Ovules~PC1*group,family=gaussian(),data=df_out) # formula, familia y datos
summary(m1.1)
dropterm(m1.1,test="F")

summary(as.factor(df_out$group))
df.glm.dep<-subset(df_out,group=="dep") #depresssa
m1.2<-glm(Ovules~PC1,family=gaussian(),data=df.glm.dep) # formula, familia y datos
summary(m1.2)
dropterm(m1.2,test="F")

summary(as.factor(df_out$group))
df.glm.gla<-subset(df_out,group=="gla") #glabrata
m1.2<-glm(Ovules~PC1,family=gaussian(),data=df.glm.gla) # formula, familia y datos
summary(m1.2)
dropterm(m1.2,test="F")

df.glm.lut<-subset(df_out,group=="gla") #lut
m1.3<-glm(Ovules~PC1,family=gaussian(),data=df.glm.lut) # formula, familia y datos
summary(m1.3)
dropterm(m1.3,test="F")

df.glm.hyb<-subset(df_out,group=="hyb") #hybridos
m1.4<-glm(Ovules~PC1,family=gaussian(),data=df.glm.hyb) # formula, familia y datos
summary(m1.4)
dropterm(m1.4,test="F")

#grafico entre variables
names(df_out)

library(ggplot2)
ggplot(df_out, aes(x=(PC1), y=(Ovules), color=group)) +
  geom_point(aes(fill=group, shape=group)) + 
  geom_smooth(method=lm, aes(color=group,fill=group),
              fullrange=T, alpha=0.2)+
  scale_color_manual(values=c('red','forestgreen', 'purple','goldenrod1'))+
  scale_fill_manual(values=c('red','forestgreen', 'purple','goldenrod1'))+
  scale_y_continuous(name="Ovule number", expand=c(0,0))+
  scale_x_continuous(name="PC1",expand=c(0,0))+
  scale_shape_manual(values=c(21, 22, 23,24))+
  theme(axis.text.x=element_text(hjust=0.5,size=18,colour="black"), 
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=18,colour="black"),
        panel.border=element_rect(colour = "black", 
                                  fill=NA, size=1.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())#,
        #legend.position = "none")#,
#coord_cartesian(xlim=c(0,0))#, ylim=c(0,2000))

#filtrado de datos
mimulus = filter(mimulus, Species != "var_x_nai" & Species != "var") # esta función filtra los datos

mimulus %>%
  ggplot(aes(x=log10(Largo.Ovario..mm.), 
             y=log10(Ancho.ovario..mm.), 
             color=Species))+
  geom_point()+xlab("Largo ovario (mm) [log10]") + ylab("Ancho ovario [log10]")+
  geom_smooth(method="lm")

