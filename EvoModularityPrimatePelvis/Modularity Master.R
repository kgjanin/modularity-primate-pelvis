## code: Katrien G. Janin 
## Evolutionary Modularity of the Primate Pelvis


## update R version/ roll back if necessary
#if(!require(installr)) {
#install.packages("installr"); require(installr)}
#library(installr)
#updateR()

################# R VERSION 3.5.2 ###################
## clear R workspace
#rm(list=ls())
## check working directory
getwd()
setwd("/Users/trientje/Dropbox/R!/EvoModularityPrimatePelvis") #change to your work location

################# Evolutionary Modularity Primate pelvic Girdle   ###################

## set working directories
input.loc <- "/Users/trientje/Dropbox/R!/EvoModularityPrimatePelvis/input"  #change to your work location

## load libraries
#library(RRPP)
#library(rgl)
library(geomorph) ## VERSION 3.2.1  
#library(phytools) ## VERSION 0.7-20
#library(ape) ## VERSION ape_5.3 
#library(maps) ## VERSION maps_3.3.0 
#library(png) ## VERSION 0.1-7
#library(Morpho) ## VERSION Morpho_2.8
#library("ggplot2") ## VERSION gplot2_3.1.0
#library(geiger) ## VERSION geiger_2.0.6.1
library("viridis") 

### check versions
sessionInfo()


##############################################################
########### import initial data - primates no human ##########
##############################################################
setwd(input.loc)
classifiers.NoHom <- read.csv("classifiers.csv", header = T)
Csize_NoHom <- read.delim("Csize.txt", row.names=1,header = T)

data_NoHom <- read.delim("primate_data.txt", row.names=1, header = T) #shape data
sym.data.NoHom <- arrayspecs(data_NoHom, 153, 3)
dim(sym.data.NoHom)

phy.NoHom <- read.nexus("phylo_gen_con.nex")  #10K website genetic consenus tree
phy.NoHom$tip.label <- gsub('_', ' ', phy.NoHom$tip.label)

partition <- read_excel("partition.xlsx") #landmark region partitions


###########################################
#### Organising data in R environment #####
###########################################
viridis(5)
### (> viridis(5) [1] (1)- HOM "#440154FF" (2)- LEM "#3B528BFF" (3)- LOR "#21908CFF"  (4)- NWM "#5DC863FF" - OWM(5)"#FDE725FF") 

## by group
group.nh.cols <- rep("", length(classifiers.NoHom$Group))
names(group.nh.cols) <- classifiers.NoHom$Group

LEM.nh.grp <- as.numeric(which(classifiers.NoHom$Group == "LEM"))
LOR.nh.grp <- as.numeric(which(classifiers.NoHom$Group == "LOR"))
NWM.nh.grp <- as.numeric(which(classifiers.NoHom$Group == "NWM"))
OWM.nh.grp <- as.numeric(which(classifiers.NoHom$Group == "OWM"))
HOM.nh.grp <- as.numeric(which(classifiers.NoHom$Group == "HOM"))

group.nh.cols[LEM.nh.grp] <- "#3B528BFF" 
group.nh.cols[LOR.nh.grp] <- "#21908CFF" 
group.nh.cols[NWM.nh.grp] <- "#5DC863FF"
group.nh.cols[OWM.nh.grp] <- "#FDE725FF" 
group.nh.cols[HOM.nh.grp] <- "#440154FF" 


## for phylomorphospaces group # gives warnings, safe to disregard
grp.nh.phy <- rep("black", length(phy.NoHom$tip.label) + phy.NoHom$Nnode) 
control=list(col.node=grp.nh.phy)
names(grp.nh.phy) <- 1:length(grp.nh.phy)
grp.nh.phy[getDescendants(phy.NoHom,LEM.nh.grp)]<-"#3B528BFF" 
grp.nh.phy[getDescendants(phy.NoHom,LOR.nh.grp)]<-"#21908CFF" 
grp.nh.phy[getDescendants(phy.NoHom,NWM.nh.grp)]<-"#5DC863FF" 
grp.nh.phy[getDescendants(phy.NoHom,OWM.nh.grp)]<-"#FDE725FF" 
grp.nh.phy[getDescendants(phy.NoHom,HOM.nh.grp)]<-"#440154FF" 


primates.pca <- plotTangentSpace(sym.data.NoHom, axis1 = 1, axis2 = 2, warpgrids=T,verbose=T, groups = group.nh.cols)


lm <-1:62 #landmarks
lm.L <-1:33 #left landmarks + midline
slm.L <-c(1:33,64:108) #left landmarks + midline + semilandmarks

###################################
#### Phy signal test #####
##################################

phy.sig.symdata.nh <- physignal(sym.data.NoHom, phy.NoHom, iter = 999, print.progress = T)
phy.sig.symdata.nh #Observed Phylogenetic Signal (K):  0.4002 - P-value: 0.001 effect size 9.4763
plot(phy.sig.symdata.nh)


LCS <-Csize_NoHom[2]
phy.sig.LCS <-physignal(LCS, phy.NoHom, iter = 999, print.progress = T)
phy.sig.LCS#Observed Phylogenetic Signal (K):  1.3935 - P-value: 0.001 effect size 8.2511
plot(phy.sig.LCS)

###################################
#### Creating  data sets   #####
##################################


#PIC data
gen.tree <- read.nexus("phylo_gen_con.nex")
coords.2d <- two.d.array(sym.data.NoHom)
PIC.data <- apply(coords.2d, 2, FUN=pic, phy=multi2di(gen.tree)) %>% arrayspecs(., p=dim(sym.data.NoHom)[1], k = 3)



# allo data - evolutionary allometry
gdf <- geomorph.data.frame(coords = sym.data.NoHom, Csize = Csize_NoHom$Log.Centroid.Size) # sym data and logCsize
allo <- procD.pgls(coords ~ Csize, phy.NoHom, iter = 999, data = gdf, verbose = T)
summary(allo) # Rsq 0.111   - 0.001 **

allo.primate.plot <- plot(allo, type = "regression", reg.type = "RegScore", predictor = Csize_NoHom$Log.Centroid.Size, 
                      col = group.nh.cols, pch= 16, cex = 1.5) 

allo.resid.array <- arrayspecs(allo$pgls.residuals, p=153, k=3)# allometry-adjusted residuals
allo.shape.resid <- allo.resid.array + array(mshape(sym.data.NoHom), dim(allo.resid.array))# allometry-free shapes
allo.data <- gpagen(allo.shape.resid, Proj = TRUE)
allo.data <-allo.data$coords


####################################################################################
#################        MODULARITY TEST  / primates order       ###################
####################################################################################

####################### shape data  ###################### 


sym.nh.h1 <-  modularity.test(sym.data.NoHom, partition.gp = partition$Dev.Mod5 , iter = 999, print.progress = T, CI = F)
summary(sym.nh.h1) 

sym.nh.h2<-  modularity.test(sym.data.NoHom, partition.gp = partition$Dev.Mod4, iter = 999, print.progress = T, CI = F)
summary(sym.nh.h2)

sym.nh.h3 <-  modularity.test(sym.data.NoHom, partition.gp = partition$H3, iter = 999, print.progress = T, CI = F)
summary(sym.nh.h3)

sym.nh.h4 <-  modularity.test(sym.data.NoHom, partition.gp = partition$Dev.Mod3, iter = 999, print.progress = T, CI = F)
summary(sym.nh.h4)

sym.nh.h5<-  modularity.test(sym.data.NoHom, partition.gp = partition$Fun.Mod2, iter = 999, print.progress = T, CI = F)
summary(sym.nh.h5) 

shape.mod <- compare.CR(sym.nh.h1,
                        sym.nh.h2,
                        sym.nh.h3, 
                        sym.nh.h4, 
                        sym.nh.h5,
                        CR.null=TRUE)
shape.mod
###################### PIC data ############################

PIC.gen.h1 <-  modularity.test(PIC.data, partition.gp = partition$Dev.Mod5 , iter = 999, print.progress = T, CI = F)
summary(PIC.gen.h1) 

PIC.gen.h2 <-  modularity.test(PIC.data, partition.gp = partition$Dev.Mod4 , iter = 999, print.progress = T, CI = F)
summary(PIC.gen.h2) 

PIC.gen.h3 <-  modularity.test(PIC.data, partition.gp = partition$H3 , iter = 999, print.progress = T, CI = F)
summary(PIC.gen.h3) 

PIC.gen.h4 <-  modularity.test(PIC.data, partition.gp = partition$Dev.Mod3 , iter = 999, print.progress = T, CI = F)
summary(PIC.gen.h4) 

PIC.gen.h5 <-  modularity.test(PIC.data, partition.gp = partition$Fun.Mod2 , iter = 999, print.progress = T, CI = F)
summary(PIC.gen.h5) 


PIC.mod <- compare.CR(PIC.gen.h1,
                      PIC.gen.h2, 
                      PIC.gen.h3, 
                      PIC.gen.h4, 
                      PIC.gen.h5, 
                      CR.null=TRUE)
PIC.mod

####################### allo data  ###################### 


allo.nh.h1 <-  modularity.test(allo.data, partition.gp = partition$Dev.Mod5 , iter = 999, print.progress = T, CI = F)
summary(allo.nh.h1) 

allo.nh.h2<-  modularity.test(allo.data, partition.gp = partition$Dev.Mod4, iter = 999, print.progress = T, CI = F)
summary(allo.nh.h2)

allo.nh.h3 <-  modularity.test(allo.data, partition.gp = partition$H3, iter = 999, print.progress = T, CI = F)
summary(allo.nh.h3)

allo.nh.h4 <-  modularity.test(allo.data, partition.gp = partition$Dev.Mod3, iter = 999, print.progress = T, CI = F)
summary(allo.nh.h4)

allo.nh.h5<-  modularity.test(allo.data, partition.gp = partition$Fun.Mod2, iter = 999, print.progress = T, CI = F)
summary(allo.nh.h5) 


allo.mod <- compare.CR(allo.nh.h1,
                       allo.nh.h2,
                       allo.nh.h3,
                       allo.nh.h4,
                       allo.nh.h5,
                      CR.null=TRUE)
allo.mod



##########################################################################################################################
#################    5 BONES MODULARITY TEST / primate order / no new gp to preserve relative position of bones    ###################
##########################################################################################################################


#LM's associated per bone
L.Ilium <- as.numeric(which(partition$Dev.Mod5_bones == "L_il")) #LM's associated with L.Ilium
L.Ischium <-as.numeric(which(partition$Dev.Mod5_bones == "L_is")) #LM's associated with L.Ilischium
L.Pubis <- as.numeric(which(partition$Dev.Mod5_bones == "L_pu")) #LM's associated with L.Pubis
L.Acetabulum <-as.numeric(which(partition$Dev.Mod5_bones == "L_ac")) #LM's associated with L.Ac
Sacrum <-as.numeric(which(partition$Dev.Mod5_bones == "sa")) #LM's associated with sacrum 

# setting partitions
il.is.part <-c(rep(1, times=29), rep(2, times=20))
il.pu.part <-c(rep(1, times=29), rep(2, times=10))
il.ac.part <-c(rep(1, times=29), rep(2, times=11))
il.sa.part <-c(rep(1, times=29), rep(2, times=13))
is.pu.part <-c(rep(1, times=20), rep(2, times=10))
is.ac.part <-c(rep(1, times=20), rep(2, times=11))
is.sa.part <-c(rep(1, times=20), rep(2, times=13))
pu.ac.part <-c(rep(1, times=10), rep(2, times=11))
pu.sa.part <-c(rep(1, times=10), rep(2, times=13))
ac.sa.part <-c(rep(1, times=11), rep(2, times=13))


##### SYM BONE DATA

# no new gpa to preserve relative positions of the bones within the overal primate pelvic girdle
# SYM.nh : creating bone pair data sets

Sym.Ilium.nh <- sym.data.NoHom[L.Ilium,,] #29
Sym.Ischium.nh <-sym.data.NoHom[L.Ischium,,] #20
Sym.Pubis.nh <- sym.data.NoHom[L.Pubis,,] #10
Sym.Acetabulum.nh <-sym.data.NoHom[L.Acetabulum,,] #11
Sym.Sacrum.nh <-sym.data.NoHom[Sacrum,,] #13

a.nh<-two.d.array(Sym.Ilium.nh)
b.nh<-two.d.array(Sym.Ischium.nh)
c.nh<-two.d.array(Sym.Pubis.nh)
d.nh<-two.d.array(Sym.Acetabulum.nh)
e.nh<-two.d.array(Sym.Sacrum.nh)

il.is.sym.nh <-cbind(a.nh,b.nh)
il.pu.sym.nh <-cbind(a.nh,c.nh)
il.ac.sym.nh <-cbind(a.nh,d.nh)
il.sa.sym.nh <-cbind(a.nh,e.nh)
is.pu.sym.nh <-cbind(b.nh,c.nh)
is.ac.sym.nh <-cbind(b.nh,d.nh)
is.sa.sym.nh <-cbind(b.nh,e.nh)
pu.ac.sym.nh <-cbind(c.nh,d.nh)
pu.sa.sym.nh <-cbind(c.nh,e.nh)
ac.sa.sym.nh <-cbind(d.nh,e.nh)

il.is.sym.nh <-arrayspecs(il.is.sym.nh,49,3)
il.pu.sym.nh <-arrayspecs(il.pu.sym.nh,39,3)
il.ac.sym.nh <-arrayspecs(il.ac.sym.nh,40,3)
il.sa.sym.nh <-arrayspecs(il.sa.sym.nh,42,3)
is.pu.sym.nh <-arrayspecs(is.pu.sym.nh,30,3)
is.ac.sym.nh <-arrayspecs(is.ac.sym.nh,31,3)
is.sa.sym.nh <-arrayspecs(is.sa.sym.nh,33,3)
pu.ac.sym.nh <-arrayspecs(pu.ac.sym.nh,21,3)
pu.sa.sym.nh <-arrayspecs(pu.sa.sym.nh,23,3)
ac.sa.sym.nh <-arrayspecs(ac.sa.sym.nh,24,3)


# PIC nh 5 bones
Phy.Ilium.nh <- PIC.data[L.Ilium,,]
Phy.Pubis.nh <- PIC.data[L.Pubis,,]
Phy.Ischium.nh <-PIC.data[L.Ischium,,]
Phy.Acetabulum.nh <-PIC.data[L.Acetabulum,,]
Phy.Sacrum.nh <-PIC.data[Sacrum,,]

p.a.nh<-two.d.array(Phy.Ilium.nh)
p.b.nh<-two.d.array(Phy.Ischium.nh)
p.c.nh<-two.d.array(Phy.Pubis.nh)
p.d.nh<-two.d.array(Phy.Acetabulum.nh)
p.e.nh<-two.d.array(Phy.Sacrum.nh)

il.is.phy.nh <-cbind(p.a.nh,p.b.nh)
il.pu.phy.nh <-cbind(p.a.nh,p.c.nh)
il.ac.phy.nh <-cbind(p.a.nh,p.d.nh)
il.sa.phy.nh <-cbind(p.a.nh,p.e.nh)
is.pu.phy.nh <-cbind(p.b.nh,p.c.nh)
is.ac.phy.nh <-cbind(p.b.nh,p.d.nh)
is.sa.phy.nh <-cbind(p.b.nh,p.e.nh)
pu.ac.phy.nh <-cbind(p.c.nh,p.d.nh)
pu.sa.phy.nh <-cbind(p.c.nh,p.e.nh)
ac.sa.phy.nh <-cbind(p.d.nh,p.e.nh)

il.is.phy.nh <-arrayspecs(il.is.phy.nh,49,3)
il.pu.phy.nh <-arrayspecs(il.pu.phy.nh,39,3)
il.ac.phy.nh <-arrayspecs(il.ac.phy.nh,40,3)
il.sa.phy.nh <-arrayspecs(il.sa.phy.nh,42,3)
is.pu.phy.nh <-arrayspecs(is.pu.phy.nh,30,3)
is.ac.phy.nh <-arrayspecs(is.ac.phy.nh,31,3)
is.sa.phy.nh <-arrayspecs(is.sa.phy.nh,33,3)
pu.ac.phy.nh <-arrayspecs(pu.ac.phy.nh,21,3)
pu.sa.phy.nh <-arrayspecs(pu.sa.phy.nh,23,3)
ac.sa.phy.nh <-arrayspecs(ac.sa.phy.nh,24,3)

# evo ALLO nh 5 bones
Allo.Ilium.nh <- allo.data[L.Ilium,,]
Allo.Pubis.nh <- allo.data[L.Pubis,,]
Allo.Ischium.nh <-allo.data[L.Ischium,,]
Allo.Acetabulum.nh <-allo.data[L.Acetabulum,,]
Allo.Sacrum.nh <-allo.data[Sacrum,,]

a.a.nh<-two.d.array(Allo.Ilium.nh)
a.b.nh<-two.d.array(Allo.Ischium.nh)
a.c.nh<-two.d.array(Allo.Pubis.nh)
a.d.nh<-two.d.array(Allo.Acetabulum.nh)
a.e.nh<-two.d.array(Allo.Sacrum.nh)

il.is.allo.nh <-cbind(a.a.nh,a.b.nh)
il.pu.allo.nh <-cbind(a.a.nh,a.c.nh)
il.ac.allo.nh <-cbind(a.a.nh,a.d.nh)
il.sa.allo.nh <-cbind(a.a.nh,a.e.nh)
is.pu.allo.nh <-cbind(a.b.nh,a.c.nh)
is.ac.allo.nh <-cbind(a.b.nh,a.d.nh)
is.sa.allo.nh <-cbind(a.b.nh,a.e.nh)
pu.ac.allo.nh <-cbind(a.c.nh,a.d.nh)
pu.sa.allo.nh <-cbind(a.c.nh,a.e.nh)
ac.sa.allo.nh <-cbind(a.d.nh,a.e.nh)

il.is.allo.nh <-arrayspecs(il.is.allo.nh,49,3)
il.pu.allo.nh <-arrayspecs(il.pu.allo.nh,39,3)
il.ac.allo.nh <-arrayspecs(il.ac.allo.nh,40,3)
il.sa.allo.nh <-arrayspecs(il.sa.allo.nh,42,3)
is.pu.allo.nh <-arrayspecs(is.pu.allo.nh,30,3)
is.ac.allo.nh <-arrayspecs(is.ac.allo.nh,31,3)
is.sa.allo.nh <-arrayspecs(is.sa.allo.nh,33,3)
pu.ac.allo.nh <-arrayspecs(pu.ac.allo.nh,21,3)
pu.sa.allo.nh <-arrayspecs(pu.sa.allo.nh,23,3)
ac.sa.allo.nh <-arrayspecs(ac.sa.allo.nh,24,3)

#### shape nh per bone pairwise modularity test

# il.is
sym.il.is.mod5.nh <- modularity.test(il.is.sym.nh, il.is.part,  iter = 999, print.progress = T, CI=F)
summary(sym.il.is.mod5.nh)
plot(sym.il.is.mod5.nh)

# il.pu 
sym.il.pu.mod5.nh <- modularity.test(il.pu.sym.nh, il.pu.part, iter = 999, print.progress = T, CI=F)
summary(sym.il.pu.mod5.nh)
plot(sym.il.pu.mod5.nh)

# il.ac 
sym.il.ac.mod5.nh <- modularity.test(il.ac.sym.nh, il.ac.part, iter = 999, print.progress = T, CI=F)
summary(sym.il.ac.mod5.nh)
plot(sym.il.ac.mod5.nh)

# il.sa 
sym.il.sa.mod5.nh <- modularity.test(il.sa.sym.nh, il.sa.part, iter = 999, print.progress = T, CI = F)
summary(sym.il.sa.mod5.nh)
plot(sym.il.sa.mod5.nh)

# is.pu
sym.is.pu.mod5.nh <- modularity.test(is.pu.sym.nh, is.pu.part, iter = 999, print.progress = T, CI = F)
summary(sym.is.pu.mod5.nh)
plot(sym.is.pu.mod5.nh)

# is.ac 
sym.is.ac.mod5.nh <- modularity.test(is.ac.sym.nh, is.ac.part, iter = 999, print.progress = T, CI = F)
summary(sym.is.ac.mod5.nh)
plot(sym.is.ac.mod5.nh)

# is.sa 
sym.is.sa.mod5.nh <- modularity.test(is.sa.sym.nh, is.sa.part, iter = 999, print.progress = T, CI = F)
summary(sym.is.sa.mod5.nh)
plot(sym.is.sa.mod5.nh)

# pu.ac 
sym.pu.ac.mod5.nh <- modularity.test(pu.ac.sym.nh, pu.ac.part, iter = 999, print.progress = T, CI = F)
summary(sym.pu.ac.mod5.nh)
plot(sym.pu.ac.mod5.nh)

# pu.sa 
sym.pu.sa.mod5.nh <- modularity.test(pu.sa.sym.nh, pu.sa.part, iter = 999, print.progress = T, CI= F)
summary(sym.pu.sa.mod5.nh)
plot(sym.pu.sa.mod5.nh)

# ac.sa 
sym.ac.sa.mod5.nh <- modularity.test(ac.sa.sym.nh, ac.sa.part, iter = 999, print.progress = T, CI =F)
summary(sym.ac.sa.mod5.nh)
plot(sym.ac.sa.mod5.nh)


# CR
CR.5mod.shape <- c(sym.il.is.mod5.nh$CR, sym.il.pu.mod5.nh$CR, sym.il.ac.mod5.nh$CR, sym.il.sa.mod5.nh$CR, sym.is.pu.mod5.nh$CR, 
                   sym.is.ac.mod5.nh$CR, sym.is.sa.mod5.nh$CR, sym.pu.ac.mod5.nh$CR, sym.pu.sa.mod5.nh$CR, sym.ac.sa.mod5.nh$CR)
round(CR.5mod.shape,digits = 3)

#compare
sym.bones.mod5.z.nh <-compare.CR(sym.il.is.mod5.nh, sym.il.pu.mod5.nh, sym.il.ac.mod5.nh, sym.il.sa.mod5.nh, sym.is.pu.mod5.nh, 
                                 sym.is.ac.mod5.nh, sym.is.sa.mod5.nh, sym.pu.ac.mod5.nh, sym.pu.sa.mod5.nh, sym.ac.sa.mod5.nh)
round(sym.bones.mod5.z.nh$sample.z, digits = 3) 



####  PIC nh per bone 

# il.is
phy.il.is.mod5.nh <- modularity.test(il.is.phy.nh, il.is.part,  iter = 999, print.progress = T, CI=F)
summary(phy.il.is.mod5.nh)
plot(phy.il.is.mod5.nh)

# il.pu 
phy.il.pu.mod5.nh <- modularity.test(il.pu.phy.nh, il.pu.part, iter = 999, print.progress = T, CI=F)
summary(phy.il.pu.mod5.nh)
plot(phy.il.pu.mod5.nh)

# il.ac 
phy.il.ac.mod5.nh <- modularity.test(il.ac.phy.nh, il.ac.part, iter = 999, print.progress = T, CI=F)
summary(phy.il.ac.mod5.nh)
plot(phy.il.ac.mod5.nh)

# il.sa 
phy.il.sa.mod5.nh <- modularity.test(il.sa.phy.nh, il.sa.part, iter = 999, print.progress = T, CI = F)
summary(phy.il.sa.mod5.nh)
plot(phy.il.sa.mod5.nh)

# is.pu
phy.is.pu.mod5.nh <- modularity.test(is.pu.phy.nh, is.pu.part, iter = 999, print.progress = T, CI = F)
summary(phy.is.pu.mod5.nh)
plot(phy.is.pu.mod5.nh)

# is.ac 
phy.is.ac.mod5.nh <- modularity.test(is.ac.phy.nh, is.ac.part, iter = 999, print.progress = T, CI = F)
summary(phy.is.ac.mod5.nh)
plot(phy.is.ac.mod5.nh)

# is.sa 
phy.is.sa.mod5.nh <- modularity.test(is.sa.phy.nh, is.sa.part, iter = 999, print.progress = T, CI = F)
summary(phy.is.sa.mod5.nh)
plot(phy.is.sa.mod5.nh)

# pu.ac 
phy.pu.ac.mod5.nh <- modularity.test(pu.ac.phy.nh, pu.ac.part, iter = 999, print.progress = T, CI = F)
summary(phy.pu.ac.mod5.nh)
plot(phy.pu.ac.mod5.nh)

# pu.sa 
phy.pu.sa.mod5.nh <- modularity.test(pu.sa.phy.nh, pu.sa.part, iter = 999, print.progress = T, CI= F)
summary(phy.pu.sa.mod5.nh)
plot(phy.pu.sa.mod5.nh)

# ac.sa 
phy.ac.sa.mod5.nh <- modularity.test(ac.sa.phy.nh, ac.sa.part, iter = 999, print.progress = T, CI =F)
summary(phy.ac.sa.mod5.nh)
plot(phy.ac.sa.mod5.nh)


# CR
CR.5mod.pic <- c(phy.il.is.mod5.nh$CR, phy.il.pu.mod5.nh$CR, phy.il.ac.mod5.nh$CR, phy.il.sa.mod5.nh$CR, phy.is.pu.mod5.nh$CR, 
                phy.is.ac.mod5.nh$CR, phy.is.sa.mod5.nh$CR, phy.pu.ac.mod5.nh$CR, phy.pu.sa.mod5.nh$CR, phy.ac.sa.mod5.nh$CR)
round(CR.5mod.pic,digits = 3)


#compare
pic.bones.mod5.z.nh <-compare.CR(phy.il.is.mod5.nh, phy.il.pu.mod5.nh, phy.il.ac.mod5.nh, phy.il.sa.mod5.nh, phy.is.pu.mod5.nh, 
                                 phy.is.ac.mod5.nh, phy.is.sa.mod5.nh, phy.pu.ac.mod5.nh, phy.pu.sa.mod5.nh, phy.ac.sa.mod5.nh)
round(pic.bones.mod5.z.nh$sample.z, digits = 3) 


####  ALLO nh per bone 

# il.is
allo.il.is.mod5.nh <- modularity.test(il.is.allo.nh, il.is.part,  iter = 999, print.progress = T, CI=F)
summary(allo.il.is.mod5.nh)
plot(allo.il.is.mod5.nh)

# il.pu 
allo.il.pu.mod5.nh <- modularity.test(il.pu.allo.nh, il.pu.part, iter = 999, print.progress = T, CI=F)
summary(allo.il.pu.mod5.nh)
plot(allo.il.pu.mod5.nh)

# il.ac 
allo.il.ac.mod5.nh <- modularity.test(il.ac.allo.nh, il.ac.part, iter = 999, print.progress = T, CI=F)
summary(allo.il.ac.mod5.nh)
plot(allo.il.ac.mod5.nh)

# il.sa 
allo.il.sa.mod5.nh <- modularity.test(il.sa.allo.nh, il.sa.part, iter = 999, print.progress = T, CI = F)
summary(allo.il.sa.mod5.nh)
plot(allo.il.sa.mod5.nh)

# is.pu
allo.is.pu.mod5.nh <- modularity.test(is.pu.allo.nh, is.pu.part, iter = 999, print.progress = T, CI = F)
summary(allo.is.pu.mod5.nh)
plot(allo.is.pu.mod5.nh)

# is.ac 
allo.is.ac.mod5.nh <- modularity.test(is.ac.allo.nh, is.ac.part, iter = 999, print.progress = T, CI = F)
summary(allo.is.ac.mod5.nh)
plot(allo.is.ac.mod5.nh)

# is.sa 
allo.is.sa.mod5.nh <- modularity.test(is.sa.allo.nh, is.sa.part, iter = 999, print.progress = T, CI = F)
summary(allo.is.sa.mod5.nh)
plot(allo.is.sa.mod5.nh)

# pu.ac 
allo.pu.ac.mod5.nh <- modularity.test(pu.ac.allo.nh, pu.ac.part, iter = 999, print.progress = T, CI = F)
summary(allo.pu.ac.mod5.nh)
plot(allo.pu.ac.mod5.nh)

# pu.sa 
allo.pu.sa.mod5.nh <- modularity.test(pu.sa.allo.nh, pu.sa.part, iter = 999, print.progress = T, CI= F)
summary(allo.pu.sa.mod5.nh)
plot(allo.pu.sa.mod5.nh)

# ac.sa 
allo.ac.sa.mod5.nh <- modularity.test(ac.sa.allo.nh, ac.sa.part, iter = 999, print.progress = T, CI =F)
summary(allo.ac.sa.mod5.nh)
plot(allo.ac.sa.mod5.nh)


# CR
CR.5mod.allo <- c(allo.il.is.mod5.nh$CR, allo.il.pu.mod5.nh$CR, allo.il.ac.mod5.nh$CR, allo.il.sa.mod5.nh$CR, allo.is.pu.mod5.nh$CR, 
                  allo.is.ac.mod5.nh$CR, allo.is.sa.mod5.nh$CR, allo.pu.ac.mod5.nh$CR, allo.pu.sa.mod5.nh$CR, allo.ac.sa.mod5.nh$CR)
round(CR.5mod.allo,digits = 3)


#compare
allo.bones.mod5.z.nh <-compare.CR(allo.il.is.mod5.nh, allo.il.pu.mod5.nh, allo.il.ac.mod5.nh, allo.il.sa.mod5.nh, allo.is.pu.mod5.nh, 
                                  allo.is.ac.mod5.nh, allo.is.sa.mod5.nh, allo.pu.ac.mod5.nh, allo.pu.sa.mod5.nh, allo.ac.sa.mod5.nh)
round(allo.bones.mod5.z.nh$sample.z, digits = 3) 



####################################################################################
#################      MODULARITY TEST  / primates group    ###################
####################################################################################

### import GPA data
setwd(input.loc)
# HOM
data_NoHom.HOM <- read.delim("sym.HOM.txt", row.names=1, header = TRUE)
classifiers.HOM <- read.csv("classifiers.HOM.csv", header = T)
Csize.HOM <- read.delim("Csize.HOM.txt", row.names=1, header = T)
sym.HOM <- arrayspecs(data_NoHom.HOM , 153, 3)


## by HOM.genus (5)
viridis(5)
HOM.genus.cols <- rep("", length(classifiers.HOM$Genus))
names(HOM.genus.cols) <- classifiers.HOM$Genus

Gor.grp <- as.numeric(which(classifiers.HOM$Genus == "Gor"))
Hyl.grp <- as.numeric(which(classifiers.HOM$Genus == "Hyl"))
Pan.grp <- as.numeric(which(classifiers.HOM$Genus == "Pan"))
Pon.grp <- as.numeric(which(classifiers.HOM$Genus == "Pon"))
Sym.grp <- as.numeric(which(classifiers.HOM$Genus == "Sym"))

HOM.genus.cols[Gor.grp] <- "#440154FF"
HOM.genus.cols[Hyl.grp] <- "#3B528BFF" 
HOM.genus.cols[Pan.grp] <- "#21908CFF" 
HOM.genus.cols[Pon.grp] <- "#5DC863FF"
HOM.genus.cols[Sym.grp] <- "#FDE725FF"
#View(HOM.genus.cols)

sym.HOM.pca <-plotTangentSpace(sym.HOM, warpgrids=T, groups = HOM.genus.cols)

allo.HOM.gdf <- geomorph.data.frame(coords = sym.HOM, Csize = Csize.HOM$Centroid.Size )
allo.HOM.procD <-procD.lm(coords ~ log(Csize), iter = 999, RRPP = T, data =  allo.HOM.gdf, verbose=T, print.progress = T)
summary(allo.HOM.procD) 

allo.hom.plot <- plot(allo.HOM.procD, type = "regression", reg.type = "RegScore", predictor = allo.HOM.gdf$Csize, 
                      col = HOM.genus.cols, pch= 16, cex = 1.5) 

allo_HOM.resid <- arrayspecs(allo.HOM.procD$residuals,p=153, k=3) # allometry-adjusted residuals
allo.hom.free.shapes.resid <-allo_HOM.resid + array(mshape(sym.HOM), dim(allo_HOM.resid)) # allometry-free data
allo.hom.free.data <- gpagen(allo.hom.free.shapes.resid)
allo.HOM <- allo.hom.free.data$coords


# LEM
setwd(input.loc)

data_LEM <- read.delim("sym.LEM.txt", row.names=1, header = TRUE)
classifiers.LEM <- read.csv("classifiers.LEM.csv", header = T)
Csize.LEM <- read.delim("Csize.LEM.txt", row.names=1, header = T)
sym.LEM <- arrayspecs(data_LEM , 153, 3)

LEM.genus.cols <- rep("", length(classifiers.LEM$Genus))
names(LEM.genus.cols) <- classifiers.LEM$Genus

Dau.grp <- as.numeric(which(classifiers.LEM$Genus == "Dau"))
Eul.grp <- as.numeric(which(classifiers.LEM$Genus == "Eul"))
Ind.grp <- as.numeric(which(classifiers.LEM$Genus == "Ind"))
Lem.grp <- as.numeric(which(classifiers.LEM$Genus == "Lem"))
Var.grp <- as.numeric(which(classifiers.LEM$Genus == "Var"))

LEM.genus.cols[Dau.grp] <- "#440154FF"
LEM.genus.cols[Eul.grp] <- "#3B528BFF" 
LEM.genus.cols[Ind.grp] <- "#21908CFF" 
LEM.genus.cols[Lem.grp] <- "#5DC863FF"
LEM.genus.cols[Var.grp] <- "#FDE725FF"


sym.LEM.pca <- plotTangentSpace(sym.LEM, axis1 = 1, axis2 = 2, warpgrids=T,
                                verbose=T, groups = LEM.genus.cols)

allo.LEM.gdf <- geomorph.data.frame(coords = sym.LEM, Csize = Csize.LEM$Centroid.Size )
allo.LEM.procD <-procD.lm(coords ~ log(Csize), iter = 999, RRPP = T, data =  allo.LEM.gdf, verbose=T, print.progress = T)
summary(allo.LEM.procD) #stat sign only at 0.05 ! 

allo.lem.plot <- plot(allo.LEM.procD, type = "regression", reg.type = "RegScore", predictor = allo.LEM.gdf$Csize, 
                      col = LEM.genus.cols, pch= 16, cex = 1.5) 

allo_LEM.resid <- arrayspecs(allo.LEM.procD$residuals,p=153, k=3) # allometry-adjusted residuals
allo.lem.free.shapes.resid <-allo_LEM.resid + array(mshape(sym.LEM), dim(allo_LEM.resid)) # allometry-free data
allo.lem.free.data <- gpagen(allo.lem.free.shapes.resid)
allo.LEM <- allo.lem.free.data$coords


# NWM
setwd(input.loc)

data_NWM <- read.delim("sym.NWM.txt", row.names=1, header = TRUE)
classifiers.NWM <- read.csv("classifiers.NWM.csv", header = T)
Csize.NWM <- read.delim("Csize.NWM.txt", row.names=1, header = T)
sym.NWM<- arrayspecs(data_NWM , 153, 3)

NWM.genus.cols <- rep("", length(classifiers.NWM$Genus))
names(NWM.genus.cols) <- classifiers.NWM$Genus

Alo.grp <- as.numeric(which(classifiers.NWM$Genus == "Alo"))
Ate.grp <- as.numeric(which(classifiers.NWM$Genus == "Ate"))
Ceb.grp <- as.numeric(which(classifiers.NWM$Genus == "Ceb"))
Lag.grp <- as.numeric(which(classifiers.NWM$Genus == "Lag"))


NWM.genus.cols[Alo.grp] <- "#440154FF"
NWM.genus.cols[Ate.grp] <- "#31688EFF"
NWM.genus.cols[Ceb.grp] <- "#35B779FF" 
NWM.genus.cols[Lag.grp] <- "#FDE725FF"
#View(NWM.genus.cols)


sym.NWM.pca <- plotTangentSpace(sym.NWM, axis1 = 1, axis2 = 2, warpgrids=T,verbose=T, groups = NWM.genus.cols)

allo.NWM.gdf <- geomorph.data.frame(coords = sym.NWM, Csize = Csize.NWM$Centroid.Size)
allo.NWM.procD <-procD.lm(coords ~ log(Csize), iter = 999, RRPP = T, data =  allo.NWM.gdf, verbose=T, print.progress = T)
summary(allo.NWM.procD) #R2 at 0.019 (0.01*)

allo.nwm.plot <- plot(allo.NWM.procD, type = "regression", reg.type = "RegScore", predictor = allo.NWM.gdf$Csize, 
                      col = NWM.genus.cols, pch= 16, cex = 1.5) 

allo_NWM.resid <- arrayspecs(allo.NWM.procD$residuals,p=153, k=3) # allometry-adjusted residuals
allo.nwm.free.shapes.resid <-allo_NWM.resid + array(mshape(sym.NWM), dim(allo_NWM.resid)) # allometry-free data
allo.NWM.free.data <- gpagen(allo.nwm.free.shapes.resid)

allo.NWM <-allo.NWM.free.data$coords

# OWM
setwd(input.loc)

data_OWM <- read.delim("sym.OWM.txt", row.names=1, header = TRUE)
classifiers.OWM <- read.csv("classifiers.OWM.csv", header = T)
Csize.OWM <- read.delim("Csize.OWM.txt", row.names=1, header = T)
sym.OWM <- arrayspecs(data_OWM , 153, 3)

OWM.genus.cols <- rep("", length(classifiers.OWM$Genus))
names(OWM.genus.cols) <- classifiers.OWM$Genus

Cer.grp <- as.numeric(which(classifiers.OWM$Genus == "Cer"))
Chl.grp <- as.numeric(which(classifiers.OWM$Genus == "Chl"))
Col.grp <- as.numeric(which(classifiers.OWM$Genus == "Col"))
Ery.grp <- as.numeric(which(classifiers.OWM$Genus == "Ery"))
Lop.grp <- as.numeric(which(classifiers.OWM$Genus == "Lop"))
Mac.grp <- as.numeric(which(classifiers.OWM$Genus == "Mac"))
Man.grp <- as.numeric(which(classifiers.OWM$Genus == "Man"))
Pap.grp <- as.numeric(which(classifiers.OWM$Genus == "Pap"))
Pil.grp <- as.numeric(which(classifiers.OWM$Genus == "Pil"))
Pre.grp <- as.numeric(which(classifiers.OWM$Genus == "Pre"))
Pyg.grp <- as.numeric(which(classifiers.OWM$Genus == "Pyg"))
Sem.grp <- as.numeric(which(classifiers.OWM$Genus == "Sem"))
The.grp <- as.numeric(which(classifiers.OWM$Genus == "The"))

OWM.genus.cols[Cer.grp] <- "#440154FF"
OWM.genus.cols[Chl.grp] <- "#481F70FF"
OWM.genus.cols[Col.grp] <- "#443A83FF"
OWM.genus.cols[Ery.grp] <- "#3B528BFF" 
OWM.genus.cols[Lop.grp] <- "#31688EFF" 
OWM.genus.cols[Mac.grp] <- "#287C8EFF"
OWM.genus.cols[Man.grp] <-  "#21908CFF"
OWM.genus.cols[Pap.grp] <-  "#20A486FF"
OWM.genus.cols[Pil.grp] <-  "#35B779FF"
OWM.genus.cols[Pre.grp] <-  "#5DC863FF"
OWM.genus.cols[Pyg.grp] <-  "#8FD744FF"
OWM.genus.cols[Sem.grp] <- "#C7E020FF" 
OWM.genus.cols[The.grp] <- "#FDE725FF"
#View(OWM.genus.cols)


sym.OWM.pca <- plotTangentSpace(sym.OWM, axis1 = 1, axis2 = 2, warpgrids=T, verbose=T, groups = OWM.genus.cols)

allo.OWM.gdf <- geomorph.data.frame(coords = sym.OWM, Csize = Csize.OWM$Centroid.Size)
allo.OWM.procD <-procD.lm(coords ~ log(Csize), iter = 999, RRPP = T, data =  allo.OWM.gdf, verbose=T, print.progress = T)
summary(allo.OWM.procD)

allo.owm.plot <- plot(allo.OWM.procD, type = "regression", reg.type = "RegScore", predictor = allo.OWM.gdf$Csize, 
                      col= OWM.genus.cols,  pch= 16, cex = 1.5) 

allo_OWM.resid <- arrayspecs(allo.OWM.procD$residuals,p=153, k=3) # allometry-adjusted residuals
allo.owm.free.shapes.resid <-allo_OWM.resid + array(mshape(sym.OWM), dim(allo_OWM.resid)) # allometry-free data
allo.OWM.free.data <- gpagen(allo.owm.free.shapes.resid)

allo.OWM <- allo.OWM.free.data$coords

##################################################################################################################################
###### HOM #####

##shape data
sym.hom.data.H1 <-modularity.test(sym.HOM, partition.gp = partition$Dev.Mod5, iter = 999, print.progress = T, CI = F)
sym.hom.data.H1
sym.hom.data.H2 <-modularity.test(sym.HOM , partition.gp = partition$Dev.Mod4_bones, iter = 999, print.progress = T, CI = F)
sym.hom.data.H2
sym.hom.data.H3 <-modularity.test(sym.HOM, partition.gp = partition$H3_bones, iter = 999, print.progress = T, CI = F)
sym.hom.data.H3
sym.hom.data.H4 <-modularity.test(sym.HOM, partition.gp = partition$Dev.Mod3_bones, iter = 999, print.progress = T, CI = F)
sym.hom.data.H4
sym.hom.data.H5 <-modularity.test(sym.HOM, partition.gp = partition$Fun.Mod2, iter = 999, print.progress = T, CI = F)
sym.hom.data.H5

sym.hom.data.CR <- c(sym.hom.data.H1$CR, sym.hom.data.H2$CR, sym.hom.data.H3$CR, sym.hom.data.H4$CR, sym.hom.data.H5$CR)
round(sym.hom.data.CR, digits = 4)

sym.hom.data.mod.z <-compare.CR(sym.hom.data.H1, sym.hom.data.H2, sym.hom.data.H3, sym.hom.data.H4, sym.hom.data.H5, CR.null = TRUE)
sym.hom.data.mod.z

##allo data
allo.hom.data.H1 <-modularity.test(allo.HOM, partition.gp = partition$Dev.Mod5_bones, iter = 999, print.progress = T, CI = F)
allo.hom.data.H1
allo.hom.data.H2 <-modularity.test(allo.HOM, partition.gp = partition$Dev.Mod4_bones, iter = 999, print.progress = T, CI = F)
allo.hom.data.H2
allo.hom.data.H3 <-modularity.test(allo.HOM, partition.gp = partition$H3_bones, iter = 999, print.progress = T, CI = F)
allo.hom.data.H3
allo.hom.data.H4 <-modularity.test(allo.HOM, partition.gp = partition$Dev.Mod3_bones, iter = 999, print.progress = T, CI = F)
allo.hom.data.H4
allo.hom.data.H5 <-modularity.test(allo.HOM, partition.gp = partition$Fun.Mod2_bones, iter = 999, print.progress = T, CI = F)
allo.hom.data.H5

allo.hom.data.CR <- c(allo.hom.data.H1$CR, allo.hom.data.H2$CR, allo.hom.data.H3$CR, allo.hom.data.H4$CR, allo.hom.data.H5$CR)
round(allo.hom.data.CR, digits = 4)

allo.hom.mod.z <-compare.CR(allo.hom.data.H1, allo.hom.data.H2, allo.hom.data.H3, allo.hom.data.H4, allo.hom.data.H5, CR.null = TRUE)
allo.hom.mod.z


###### LEM #####

##shape data
sym.LEM.data.H1 <-modularity.test(sym.LEM, partition.gp = partition$Dev.Mod5, iter = 999, print.progress = T, CI = F)
sym.LEM.data.H1
sym.LEM.data.H2 <-modularity.test(sym.LEM , partition.gp = partition$Dev.Mod4_bones, iter = 999, print.progress = T, CI = F)
sym.LEM.data.H2
sym.LEM.data.H3 <-modularity.test(sym.LEM, partition.gp = partition$H3_bones, iter = 999, print.progress = T, CI = F)
sym.LEM.data.H3
sym.LEM.data.H4 <-modularity.test(sym.LEM, partition.gp = partition$Dev.Mod3_bones, iter = 999, print.progress = T, CI = F)
sym.LEM.data.H4
sym.LEM.data.H5 <-modularity.test(sym.LEM, partition.gp = partition$Fun.Mod2, iter = 999, print.progress = T, CI = F)
sym.LEM.data.H5

sym.LEM.data.CR <- c(sym.LEM.data.H1$CR, sym.LEM.data.H2$CR, sym.LEM.data.H3$CR, sym.LEM.data.H4$CR, sym.LEM.data.H5$CR)
round(sym.LEM.data.CR, digits = 4)

sym.LEM.data.mod.z <-compare.CR(sym.LEM.data.H1, sym.LEM.data.H2, sym.LEM.data.H3, sym.LEM.data.H4, sym.LEM.data.H5, CR.null = TRUE)
sym.LEM.data.mod.z

##allo data
allo.LEM.data.H1 <-modularity.test(allo.LEM, partition.gp = partition$Dev.Mod5_bones, iter = 999, print.progress = T, CI = F)
allo.LEM.data.H1
allo.LEM.data.H2 <-modularity.test(allo.LEM, partition.gp = partition$Dev.Mod4_bones, iter = 999, print.progress = T, CI = F)
allo.LEM.data.H2
allo.LEM.data.H3 <-modularity.test(allo.LEM, partition.gp = partition$H3_bones, iter = 999, print.progress = T, CI = F)
allo.LEM.data.H3
allo.LEM.data.H4 <-modularity.test(allo.LEM, partition.gp = partition$Dev.Mod3_bones, iter = 999, print.progress = T, CI = F)
allo.LEM.data.H4
allo.LEM.data.H5 <-modularity.test(allo.LEM, partition.gp = partition$Fun.Mod2_bones, iter = 999, print.progress = T, CI = F)
allo.LEM.data.H5

allo.LEM.data.CR <- c(allo.LEM.data.H1$CR, allo.LEM.data.H2$CR, allo.LEM.data.H3$CR, allo.LEM.data.H4$CR, allo.LEM.data.H5$CR)
round(allo.LEM.data.CR, digits = 4)

allo.LEM.mod.z <-compare.CR(allo.LEM.data.H1, allo.LEM.data.H2, allo.LEM.data.H3, allo.LEM.data.H4, allo.LEM.data.H5, CR.null = TRUE)
allo.LEM.mod.z


###### NWM #####

##shape data
sym.NWM.data.H1 <-modularity.test(sym.NWM, partition.gp = partition$Dev.Mod5, iter = 999, print.progress = T, CI = F)
sym.NWM.data.H1
sym.NWM.data.H2 <-modularity.test(sym.NWM , partition.gp = partition$Dev.Mod4_bones, iter = 999, print.progress = T, CI = F)
sym.NWM.data.H2
sym.NWM.data.H3 <-modularity.test(sym.NWM, partition.gp = partition$H3_bones, iter = 999, print.progress = T, CI = F)
sym.NWM.data.H3
sym.NWM.data.H4 <-modularity.test(sym.NWM, partition.gp = partition$Dev.Mod3_bones, iter = 999, print.progress = T, CI = F)
sym.NWM.data.H4
sym.NWM.data.H5 <-modularity.test(sym.NWM, partition.gp = partition$Fun.Mod2, iter = 999, print.progress = T, CI = F)
sym.NWM.data.H5

sym.NWM.data.CR <- c(sym.NWM.data.H1$CR, sym.NWM.data.H2$CR, sym.NWM.data.H3$CR, sym.NWM.data.H4$CR, sym.NWM.data.H5$CR)
round(sym.NWM.data.CR, digits = 4)

sym.NWM.data.mod.z <-compare.CR(sym.NWM.data.H1, sym.NWM.data.H2, sym.NWM.data.H3, sym.NWM.data.H4, sym.NWM.data.H5, CR.null = TRUE )
sym.NWM.data.mod.z

##allo data
allo.NWM.data.H1 <-modularity.test(allo.NWM, partition.gp = partition$Dev.Mod5_bones, iter = 999, print.progress = T, CI = F)
allo.NWM.data.H1
allo.NWM.data.H2 <-modularity.test(allo.NWM, partition.gp = partition$Dev.Mod4_bones, iter = 999, print.progress = T, CI = F)
allo.NWM.data.H2
allo.NWM.data.H3 <-modularity.test(allo.NWM, partition.gp = partition$H3_bones, iter = 999, print.progress = T, CI = F)
allo.NWM.data.H3
allo.NWM.data.H4 <-modularity.test(allo.NWM, partition.gp = partition$Dev.Mod3_bones, iter = 999, print.progress = T, CI = F)
allo.NWM.data.H4
allo.NWM.data.H5 <-modularity.test(allo.NWM, partition.gp = partition$Fun.Mod2_bones, iter = 999, print.progress = T, CI = F)
allo.NWM.data.H5

allo.NWM.data.CR <- c(allo.NWM.data.H1$CR, allo.NWM.data.H2$CR, allo.NWM.data.H3$CR, allo.NWM.data.H4$CR, allo.NWM.data.H5$CR)
round(allo.NWM.data.CR, digits = 4)

allo.NWM.mod.z <-compare.CR(allo.NWM.data.H1, allo.NWM.data.H2, allo.NWM.data.H3, allo.NWM.data.H4, allo.NWM.data.H5, CR.null = TRUE)
allo.NWM.mod.z





###### OWM #####

##shape data
sym.OWM.data.H1 <-modularity.test(sym.OWM, partition.gp = partition$Dev.Mod5, iter = 999, print.progress = T, CI = F)
sym.OWM.data.H1
sym.OWM.data.H2 <-modularity.test(sym.OWM , partition.gp = partition$Dev.Mod4_bones, iter = 999, print.progress = T, CI = F)
sym.OWM.data.H2
sym.OWM.data.H3 <-modularity.test(sym.OWM, partition.gp = partition$H3_bones, iter = 999, print.progress = T, CI = F)
sym.OWM.data.H3
sym.OWM.data.H4 <-modularity.test(sym.OWM, partition.gp = partition$Dev.Mod3_bones, iter = 999, print.progress = T, CI = F)
sym.OWM.data.H4
sym.OWM.data.H5 <-modularity.test(sym.OWM, partition.gp = partition$Fun.Mod2, iter = 999, print.progress = T, CI = F)
sym.OWM.data.H5

sym.OWM.data.CR <- c(sym.OWM.data.H1$CR, sym.OWM.data.H2$CR, sym.OWM.data.H3$CR, sym.OWM.data.H4$CR, sym.OWM.data.H5$CR)
round(sym.OWM.data.CR, digits = 4)

sym.OWM.data.mod.z <-compare.CR(sym.OWM.data.H1, sym.OWM.data.H2, sym.OWM.data.H3, sym.OWM.data.H4, sym.OWM.data.H5, CR.null = TRUE)
sym.OWM.data.mod.z

##allo data
allo.OWM.data.H1 <-modularity.test(allo.OWM, partition.gp = partition$Dev.Mod5_bones, iter = 999, print.progress = T, CI = F)
allo.OWM.data.H1
allo.OWM.data.H2 <-modularity.test(allo.OWM, partition.gp = partition$Dev.Mod4_bones, iter = 999, print.progress = T, CI = F)
allo.OWM.data.H2
allo.OWM.data.H3 <-modularity.test(allo.OWM, partition.gp = partition$H3_bones, iter = 999, print.progress = T, CI = F)
allo.OWM.data.H3
allo.OWM.data.H4 <-modularity.test(allo.OWM, partition.gp = partition$Dev.Mod3_bones, iter = 999, print.progress = T, CI = F)
allo.OWM.data.H4
allo.OWM.data.H5 <-modularity.test(allo.OWM, partition.gp = partition$Fun.Mod2_bones, iter = 999, print.progress = T, CI = F)
allo.OWM.data.H5

allo.OWM.data.CR <- c(allo.OWM.data.H1$CR, allo.OWM.data.H2$CR, allo.OWM.data.H3$CR, allo.OWM.data.H4$CR, allo.OWM.data.H5$CR)
round(allo.OWM.data.CR, digits = 4)

allo.OWM.mod.z <-compare.CR(allo.OWM.data.H1, allo.OWM.data.H2, allo.OWM.data.H3, allo.OWM.data.H4, allo.OWM.data.H5, CR.null = TRUE)
allo.OWM.mod.z
