
##################################################################
##################################################################
##################################################################

library(phyloseq)
library(metagMisc) #devtools::install_github("vmikk/metagMisc")
library(ade4)
library(ggd)
library(MicEco)# install_github("Russel88/MicEco")

##################################################################
##################################################################
##################################################################

setwd("/Users/francoisrigal/Documents/THESE_TOURBIERES/BERNADOUZE/R")

load("bernadouze.Rdata")

### Indices de diversité phylogenetique ###

MPD <- metagMisc::phyloseq_phylo_div(bernadouze, measures = c("MPD"))

PD <- metagMisc::phyloseq_phylo_div(bernadouze, measures = c("PD"))

### standardized effect size du modele nul avec pool complet (180 échantillons) et 
### technique de la réassignation aléatoire des noms de tips. Méthode à choisir par défaut

sesMPD <- metagMisc::phyloseq_phylo_div(bernadouze, measures = c("SES.MPD"))

sesPD <- metagMisc::phyloseq_phylo_div(bernadouze, measures = c("SES.PD"))


# avec picante.....

sesMPD.pic <- picante::ses.mpd(t(bernadouze@otu_table), cophenetic(bernadouze@phy_tree), 
                               abundance.weighted=T, null.model = "taxa.labels", runs=9)

sesPD.pic <- picante::ses.pd(t(bernadouze@otu_table), bernadouze@phy_tree, 
                             abundance.weighted=T, null.model = "taxa.labels", runs=9)


### Utilisation de la fonction pour randomizeMatrix pour redistribuer aléatoirement la présence et 
### l’abondance des espèces au sein de la matrice de communauté en imposant des contraintes 


MPD.swap <- list()

for (i in 1:100)
{
  rd.mat <- randomizeMatrix(t(bernadouze@otu_table), null.model = "independentswap")
  
  MPD.swap[[i]] <- PhyloMeasures::mpd.query(bernadouze@phy_tree, rd.mat)
  
  print(i)
  
}


tab.MPD.swap <- do.call(cbind, MPD.swap)

ses.MPD.swap <- (MPD$MPD - apply(tab.MPD.swap, 1, mean))/apply(tab.MPD.swap, 1, sd)

df.plot.div.phylo <- data.frame(MPD, PD, sesMPD, sesPD, ses.MPD.swap, bernadouze@sam_data)



# plot par zone, par habitat, par profondeur, avec PH

library(RColorBrewer)

df.plot.div.phylo$PROF <- factor(df.plot.div.phylo$PROF, levels = c("Bas", "Milieu", "Haut"))

ggplot(df.plot.div.phylo, aes(x=ZONE, y=SES.MPD, color=ZONE)) + geom_boxplot() + 
  theme_bw() + scale_color_brewer(palette="Dark2", name="Zones") 

ggplot(df.plot.div.phylo, aes(x=HAB, y=SES.MPD, color=HAB)) + geom_boxplot()+ 
  theme_bw() + scale_color_manual(values = c("green4", "green", "goldenrod4"), name="Habitats") 

ggplot(df.plot.div.phylo, aes(x=PROF, y=SES.MPD, color=PROF)) + geom_boxplot()+ 
  theme_bw() + scale_color_manual(values = c("grey80", "grey60", "black"), name="Depths") 

ggplot(df.plot.div.phylo, aes(x=PH, y=SES.MPD)) + geom_point() + geom_smooth(method="lm") + 
  theme_bw()

ggplot(df.plot.div.phylo, aes(x=PH, y=SES.MPD)) + geom_point() + geom_smooth(method="loess") + 
  theme_bw()

summary(lm(SES.MPD ~ PH+I(PH^2), data=df.plot.div.phylo))



# on va selectionner un groupe de sites pour l'exemple #

# les échantillons du bas de la butte dans la zone 3 #

but_haut_3 <- filter_taxa(subset_samples(bernadouze, HAB=="Butte" & PROF=="Haut" & ZONE==3), function(x) sum(x) > 0, TRUE) # uniquement les échantillons des buttes

sesMPD[rownames(but_haut_3@sam_data),]

mean(sesMPD[rownames(but_haut_3@sam_data),])


# definition des pools d'espèces #

buttes <- filter_taxa(subset_samples(bernadouze, HAB=="Butte" ), function(x) sum(x) > 0, TRUE) # uniquement les échantillons des buttes

prof_haut <- filter_taxa(subset_samples(bernadouze, PROF=="Haut"), function(x) sum(x) > 0, TRUE)# uniquement les échantillons de surfaces

buttes_haut <- filter_taxa(subset_samples(bernadouze, HAB=="Butte" & PROF=="Haut"), function(x) sum(x) > 0, TRUE) # uniquement les échantillons des buttes

zone_3 <- filter_taxa(subset_samples(bernadouze, ZONE==3), function(x) sum(x) > 0, TRUE) # uniquement les échantillons de la zone ombrotrophe 
 
zone_3_Haut <- filter_taxa(subset_samples(bernadouze, ZONE==3 & PROF=="Haut"), function(x) sum(x) > 0, TRUE) # uniquement les échantillons de la zone ombrotrophe 



# selectionner par affinités envir.

colnames( bernadouze@sam_data)

environment <- bernadouze@sam_data[,c("PH", "FREE_PHENOL","CN_ratio", "NUT", "DOC", "TEMP", "WFPS")]

km_env <- kmeans(scale(environment), centers = 3)

pca_env <- dudi.pca(as.matrix(environment), scannf=F, nf=2)

s.corcircle(pca_env$co)

plot(pca_env$li[,1:2]);ordispider(pca_env$li[,1:2], groups = km_env$cluster, label = T)

points(pca_env$li[,1:2][rownames(but_bas_3@sam_data),], col="red") # ou sont mes 4 points ?

bernadouze@sam_data$cluster_env <- km_env$cluster # ajouter le facteur "cluster" au data

cluster <- filter_taxa(subset_samples(bernadouze, cluster_env==2), function(x) sum(x) > 0, TRUE) # uniquement les échantillons du cluster 1


# regardons les SES par type de pool #



ses.buttes <- metagMisc::phyloseq_phylo_div(buttes, measures = c("SES.MPD"))

ses.prof_haut <- metagMisc::phyloseq_phylo_div(prof_haut, measures = c("SES.MPD"))

ses.zone_3 <- metagMisc::phyloseq_phylo_div(zone_3, measures = c("SES.MPD"))

ses.cluster <- metagMisc::phyloseq_phylo_div(cluster, measures = c("SES.MPD"))

ses.buttes_haut <- metagMisc::phyloseq_phylo_div(buttes_haut, measures = c("SES.MPD"))

ses.zone_3_Haut <- metagMisc::phyloseq_phylo_div(zone_3_Haut, measures = c("SES.MPD"))

ses.but_haut_3 <- metagMisc::phyloseq_phylo_div(but_haut_3, measures = c("SES.MPD"))


sesMPD[rownames(but_haut_3@sam_data),]%>%mean

ses.buttes[rownames(but_haut_3@sam_data),]%>%mean

ses.prof_haut[rownames(but_haut_3@sam_data),]%>%mean

ses.zone_3[rownames(but_haut_3@sam_data),]%>%mean

ses.cluster[rownames(but_haut_3@sam_data),]%>%mean

ses.buttes_haut[rownames(but_haut_3@sam_data),]%>%mean

ses.zone_3_Haut[rownames(but_haut_3@sam_data),]%>%mean




