#source("http://bioconductor.org/biocLite.R")
library(vegan)
library(ggplot2)
library(grid)
library(ggrepel)
library(fossil)
library(reshape2)
library(metagenomeSeq)
library(coda.base)
library(zCompositions)
library(ALDEx2)
library(Tjazi)
library(iNEXT)
library(patchwork)
library(lme4)
library(rstatix)
library("ggsignif")
library("ggforce")
library("magick")
library(ggpubr)
library(rsvg)
library(grImport2)
library(png)
library(grid)
library(cowplot)
library(magick)
library(pdftools)

#library(devtools)
#devtools::install_github(repo = "thomazbastiaanssen/Tjazi")

options(stringsAsFactors = F)
getwd()
setwd("/home/thomaz/Documents/PhD/FMT aging/")

counts    <- read.delim("genus_table_from_dada2_Aging_FMT.csv", sep=",", row.names = 1) #mind the encoding UTF-7 vs UTF-8 issues
counts    = counts[-209,]

GBM       <- read.delim("piphillin/20200219160155__Thomaz_keggPiphillin/GBMs_labeled_FMT_AGING_2020.csv", sep = ",", row.names = 1)

metadata  <- read.delim("metadata_Pool3_included.csv", sep = ",")
#metadata  <- metadata[metadata$Treatment != "Young.CTR" & metadata$Treatment!= "Aged.CTR",]
metadata  <- metadata[metadata$Experiment != "GF",]
metadata = metadata[metadata$Mouse_ID %in% intersect(metadata[metadata$Timepoint == "Pre",]$Mouse_ID,
                                                     metadata[metadata$Timepoint == "Post",]$Mouse_ID),  ]

metadata$KW_legend <- metadata$Legend
metadata$KW_legend[metadata$KW_legend != "Aged oFMT"] <- "All yFMT"


#We should have 75 samples remaining per timepoint
sum(metadata$Timepoint == "Pre")
sum(metadata$Timepoint == "Post")

metadata$order_for_volatility  <- paste(metadata$Timepoint, metadata$Mouse_ID)
metadata  <- metadata[order(metadata$order_for_volatility),]
counts    <- counts  [,metadata$Sample_ID]
GBM       <- floor(GBM[,metadata$Sample_ID])

alphadiv <- get_asymptotic_alpha(species = counts)
alphadiv$Legend <- factor(paste(metadata$Timepoint, metadata$Legend), levels = unique(paste(metadata$Timepoint, metadata$Legend))[c(6, 3, 5, 2, 4, 1)])

species   <- counts



species   <- apply(species,c(1,2),function(x) as.numeric(as.character(x)))
species   <- species[apply(species == 0, 1, sum) <= (ncol(species) * 0.90), ]   #remove rows with 2 or fewer hits
#rspecies  <- apply(species, 2, function(i) i/sum(i))                        #normalize to relative abundance
#species   <- species[rowMax(rspecies) > 0.01,]


#species  <- t(species)

####
####Calculate Aitchison distance and prepare for PCA
####
species.clr <- aldex.clr(species, mc.samples = 1000, denom="all", verbose=TRUE, useMC=TRUE) 

species.exp <- Tjazi::get_aldex_exp(species.clr)

data.a.pca  <- prcomp(t(species.exp))




pc1 <- round(data.a.pca$sdev[1]^2/sum(data.a.pca$sdev^2),4) *100
pc2 <- round(data.a.pca$sdev[2]^2/sum(data.a.pca$sdev^2),4) *100
pc3 <- round(data.a.pca$sdev[3]^2/sum(data.a.pca$sdev^2),4) *100
pc4 <- round(data.a.pca$sdev[4]^2/sum(data.a.pca$sdev^2),4) *100

pca  = data.frame(PC1 = data.a.pca$x[,1], 
                  PC2 = data.a.pca$x[,2], 
                  PC3 = data.a.pca$x[,3], 
                  PC4 = data.a.pca$x[,4])
#metadata          <- metadata[metadata$Mapping_file == "Adolescence_002",]
pca$ID                   = metadata$Sample_ID
pca$Longname             = metadata$Mouse_ID
pca$Treatment            = metadata$Treatment
pca$Legend               = metadata$Legend
pca$Timepoint            = metadata$Timepoint
pca$group                = paste(metadata$Timepoint, metadata$Legend)


metadata

ta.pca <- (t(data.a.pca$x ))* ((data.a.pca$sdev^2)/sum(data.a.pca$sdev^2))
a.pca  <- t(ta.pca)
loadings <- as.data.frame(data.a.pca$rotation) 
loadings.sub = loadings[which(abs(loadings$PC1)>0.2 | abs(loadings$PC2)>0.2),]


a = ggplot(pca, aes(x=PC1, y=PC2, fill = Legend, label = ID, shape = Timepoint, group = paste(Timepoint, Legend ))) + 
#  geom_line(size = 1.5, alpha = 0.1) +
  stat_ellipse(geom = "polygon", alpha = 1/4)+
  #stat_ellipse()+
  
  geom_point(size=4, stroke = 1) +
  #geom_label_repel()+
  xlab(paste("PC1: ", pc1,  "%", sep="")) + 
  ylab(paste("PC2: ", pc2,  "%", sep="")) + 
  theme_bw() + 
  scale_shape_manual( values = c("Pre"  = 24, 
                                 "Post" = 21, 
                                 "OMT"  = 23, 
                                 "YMT"  = 23)) +
  scale_fill_manual( values = c("Aged oFMT"  = "#b2182b",
                                "Aged yFMT"  = "#ef8a62",
                                "Young yFMT" = "#2166ac",
                                "oFMT"       = "#b2182b", 
                                "yFMT"       = "#386cb0",
                                "Young CTR"  = "#b3b3b3", 
                                "Aged CTR"   = "#666666"
  ))+
  
  scale_color_manual(values = c( "Pre"  = "#666666", 
                                 "Post" = "#000000", 
                                 "OMT"  = "#fdae61", 
                                 "YMT"  = "#fdae61"
  )
  )+ 
  theme(text = element_text(size = 20) , legend.position = "none") + #aspect.ratio = 1, 
  guides(fill = guide_legend(override.aes = list(shape= 24))) 
# scale_shape_manual(values = c("1" = 21, 
geom_segment(data = loadings.sub, aes(x=0,y=0,xend=PC1*50,yend=PC2*50),arrow=arrow(length=unit(0.1,"cm")), color = "black")+
  geom_text(data = loadings.sub, aes(x=PC1*50, y=PC2*50, label=row.names(loadings.sub)),color="black") 


a

means_per_group = function(PCA, groups, relevant.dims = ncol(PCA)){
  means_df = data.frame(1:relevant.dims)
  row.names(means_df) = colnames(PCA)[1:relevant.dims]
  for(group in 1:length(unique(groups))){
    for(n in 1:relevant.dims){
      means_df[n,group] = mean(PCA[which(groups == unique(groups)[group]),n])
    }
    colnames(means_df)[group] = as.character(unique(groups)[group])
  }
  return(means_df)
}

means_pg <- means_per_group(PCA = data.a.pca$x, groups = paste(metadata$Timepoint, metadata$Legend))


pcams  = data.frame(t(means_pg[1:4,]))
#View(pcams)
#metadata          <- metadata[metadata$Mapping_file == "Adolescence_002",]

pcams $Legend               = row.names(pcams)


b = a +  
  geom_curve(aes(x = pcams$PC1[6],     y = pcams$PC2[6], 
                 xend = pcams$PC1[3], yend = pcams$PC2[3]), 
             colour = "black", size = 2.5, arrow = arrow()) +
  geom_curve(aes(x = pcams$PC1[5],     y = pcams$PC2[5], 
                 xend = pcams$PC1[2], yend = pcams$PC2[2]), 
             colour = "black", size = 2.5, arrow = arrow()) +
  geom_curve(aes(x = pcams$PC1[4],     y = pcams$PC2[4], 
                 xend = pcams$PC1[1], yend = pcams$PC2[1]), 
             colour = "black", size = 2.5, arrow = arrow(), curvature = -0.5) +
  geom_curve(aes(x = pcams$PC1[6],     y = pcams$PC2[6],       #Young yFMT
                 xend = pcams$PC1[3], yend = pcams$PC2[3]), 
             colour = "#2166ac", size = 1, arrow = arrow()) +
  geom_curve(aes(x = pcams$PC1[5],     y = pcams$PC2[5], 
                 xend = pcams$PC1[2], yend = pcams$PC2[2]), 
             colour = "#ef8a62", size = 1, arrow = arrow()) +
  geom_curve(aes(x = pcams$PC1[4],     y = pcams$PC2[4], 
                 xend = pcams$PC1[1], yend = pcams$PC2[1]), 
             colour = "#b2182b", size = 1, arrow = arrow(), curvature = -0.5) 

b 

pca_leg <- image_read_svg('for_pca_legend.svg', width = 194, height = 150)

ggpca_leg = ggplot() + 
  draw_image(pca_leg) + theme_nothing()  + 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x = NULL, y = NULL) 
ggpca_leg
?draw_image
res_df_mw <- pairwise_DA_wrapper(reads = species, 
                                 groups = paste(metadata$Timepoint, metadata$KW_legend), 
                                 comparisons = data.frame(
                                   a = c("Post Aged oFMT"),
                                   b =  c("Post All yFMT")), 
                                 parametric = F, 
                                 ignore.posthoc = F)

View(res_df_mw[res_df_mw$`Post Aged oFMT vs Post All yFMT BH.adjusted.p.value` < 0.2,])

gg_df <- data.frame(microbe = unlist(species.exp["Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus",]), 
                    Legend = factor(paste(metadata$Timepoint, metadata$Legend), levels = unique(paste(metadata$Timepoint, metadata$Legend))[c(6, 3, 5, 2, 4, 1)]))

ent <- ggplot(gg_df, aes(x = Legend, y = microbe,  fill = Legend, shape = Legend)) +
  geom_violin(alpha = 1/2, trim = F) + 
  #geom_boxplot(coef = 1000, width = 0.3, alpha = 1/2)+
  geom_dotplot(binaxis  = 'y',  
               stackdir = "center",  
               dotsize  = 2.2) +
  
  scale_fill_manual( values = c("Pre Aged oFMT"  = "#b2182b",
                                "Pre Aged yFMT"  = "#ef8a62",
                                "Pre Young yFMT" = "#2166ac",
                                "Post Aged oFMT"  = "#b2182b",
                                "Post Aged yFMT"  = "#ef8a62",
                                "Post Young yFMT" = "#2166ac"
  ))  + theme_bw()+ guides(fill = "none", shape = "none") + 
  ylab("Enterococcus abundance (clr)") + 
  xlab("") + 
  ggtitle("Enterococcus colonizes in Aged microbiome after yFMT") +
  theme(text = element_text(size = 12)) + geom_hline(yintercept = -1, col = "red", linetype = "dashed") +
  annotate(geom = "text", x = Inf, y = -Inf, 
                          label = paste("q = ", 
                                        round(res_df_mw[res_df_mw$microbe == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus",]$`Post Aged oFMT vs Post All yFMT BH.adjusted.p.value`, 
                                              digits = 3) ), hjust = 1, vjust = -0.25) +
  annotate(geom = "text", x = 5.5, y = Inf, label = "***", vjust = 2)

ent
species["Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus",]  
species.exp["Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus",]  

ent  


dim(a.pca)
#ndim = 1:10
disait = sqrt(rowSums( (a.pca[1:75,1:10] - a.pca[76:150,1:10])^2) ) #save for future comparison with BC

gg_ait <- data.frame(disait = disait, 
                     Legend = factor(metadata$Legend[1:75], levels = unique(metadata$Legend[1:75])[c(3, 2, 1)]))

vola <- ggplot(data = gg_ait, aes(x = Legend, y = disait, fill = Legend)) +   
  geom_violin(alpha = 1/2, trim = F, draw_quantiles = 0.5) + 
  #geom_boxplot( width = 0.3, alpha = 1/2, coef = 10000)+
  geom_dotplot(binaxis  = 'y',  
               stackdir = "center",  
               dotsize  = 2.5) +
  
  scale_fill_manual( values = c("Aged oFMT"   = "#b2182b",
                                "Aged yFMT"   = "#ef8a62",
                                "Young yFMT"  = "#2166ac"
  ))  + theme_bw()+ guides(fill = "none") + xlab("") + ylab("Aitchison distance") + 
  theme(text = element_text(size = 12)) + ggtitle("Microbiome Change During Treatment") +
  annotate(geom = "text", x = 1.5, y = Inf, label = "*", vjust = 2) +
  annotate(geom = "text", x = 2.5, y = Inf, label = "#", vjust = 2) 


vola
summary(aaa)
aaa <- aov(data = gg_ait, disait~Legend)
TukeyHSD(aaa)
a



chao <- ggplot(alphadiv, aes(x = Legend, y = chao1,  fill = Legend)) +
  geom_violin(alpha = 1/2, trim = F, draw_quantiles = 0.5) + 
 # geom_boxplot(coef = 1000, width = 0.3, alpha = 1/2)+
  geom_dotplot(binaxis  = 'y',  
               stackdir = "center",  
               dotsize  = 2) +
  
  scale_fill_manual( values = c("Pre Aged oFMT"  = "#b2182b",
                                "Pre Aged yFMT"  = "#ef8a62",
                                "Pre Young yFMT" = "#2166ac",
                                "Post Aged oFMT"  = "#b2182b",
                                "Post Aged yFMT"  = "#ef8a62",
                                "Post Young yFMT" = "#2166ac"
  ))  + theme_bw()+ guides(fill = "none") + ylab("") + xlab("Young yFMT     Aged oFMT     Aged yFMT") + ggtitle("Chao1") + 
  scale_x_discrete(labels=c("Pre Aged oFMT"   = "Pre",
                            "Pre Aged yFMT"   = "Pre",
                            "Pre Young yFMT"  = "Pre",
                            "Post Aged oFMT"  = "Post",
                            "Post Aged yFMT"  = "Post",
                            "Post Young yFMT" = "Post")) +
  annotate(geom = "text", x = 1.5, y = Inf, label = "***", vjust = 2) +
  annotate(geom = "text", x = 3.5, y = Inf, label = "***", vjust = 2) +
  annotate(geom = "text", x = 5.5, y = Inf, label = "***", vjust = 2) 


chao


simps <- ggplot(alphadiv, aes(x = Legend, y = asymptotic_simps,  fill = Legend)) +
  geom_violin(alpha = 1/2, trim = F, draw_quantiles = 0.5) + 
  #geom_boxplot(coef = 1000, width = 0.3, alpha = 1/2)+
  geom_dotplot(binaxis  = 'y',  
               stackdir = "center",  
               dotsize  = 2) +
  
  scale_fill_manual( values = c("Pre Aged oFMT"  = "#b2182b",
                                "Pre Aged yFMT"  = "#ef8a62",
                                "Pre Young yFMT" = "#2166ac",
                                "Post Aged oFMT"  = "#b2182b",
                                "Post Aged yFMT"  = "#ef8a62",
                                "Post Young yFMT" = "#2166ac"
  ))  + theme_bw()+ guides(fill = "none") + ylab("")  + xlab("Young yFMT     Aged oFMT     Aged yFMT")  + ggtitle("Simpson Index") + 
  scale_x_discrete(labels=c("Pre Aged oFMT"   = "Pre",
                            "Pre Aged yFMT"   = "Pre",
                            "Pre Young yFMT"  = "Pre",
                            "Post Aged oFMT"  = "Post",
                            "Post Aged yFMT"  = "Post",
                            "Post Young yFMT" = "Post"))

simps




shan <- ggplot(alphadiv, aes(x = Legend, y = asymptotic_shannon,  fill = Legend)) +
  geom_violin(alpha = 1/2, trim = F, draw_quantiles = 0.5) + 
#  geom_boxplot(coef = 1000, width = 0.3, alpha = 1/2)+
  geom_dotplot(binaxis  = 'y',  
               stackdir = "center",  
               dotsize  = 2) +
  
  scale_fill_manual( values = c("Pre Aged oFMT"  = "#b2182b",
                                "Pre Aged yFMT"  = "#ef8a62",
                                "Pre Young yFMT" = "#2166ac",
                                "Post Aged oFMT"  = "#b2182b",
                                "Post Aged yFMT"  = "#ef8a62",
                                "Post Young yFMT" = "#2166ac"
  ))  + theme_bw()+ guides(fill = "none")  + ylab("") + xlab("Young yFMT     Aged oFMT     Aged yFMT")  + ggtitle("Shannon Index") + 
  scale_x_discrete(labels=c("Pre Aged oFMT"   = "Pre",
                            "Pre Aged yFMT"   = "Pre",
                            "Pre Young yFMT"  = "Pre",
                            "Post Aged oFMT"  = "Post",
                            "Post Aged yFMT"  = "Post",
                            "Post Young yFMT" = "Post"))
shan

comp <- read.delim("rel_comp_FMT_effect.csv", header = F, sep = ",")
comp

species_df <- pairwise_DA_wrapper(reads = species, 
                              groups = paste(metadata$Timepoint, metadata$Legend, sep = " "), 
                              comparisons = comp,
                              paired.test = T)

res_df <- species_df
row.names(res_df) <- res_df$microbe

cols_containing_pvals <- 1:((ncol(res_df))/2)*2
cols_containing_evals <- (1:((ncol(res_df))/2)*2) +1
change_NA_to_1      <- function(x) { replace(x, is.na(x), 1) }
change_NA_to_0      <- function(x) { replace(x, is.na(x), 0) }

res_df[,cols_containing_pvals] <- change_NA_to_1(res_df[,cols_containing_pvals])
res_df[,cols_containing_evals] <- change_NA_to_0(res_df[,cols_containing_evals])

#############
#############Set up filter to pick displayed microbes based on p-values
#############
View(head(res_df))

signig_filter = rep(F, nrow(res_df))
limit_filter = 0.1
limit_e      = 0.5
for(number in (1:nrow(res_df))){
  for(pval in (cols_containing_pvals)){
    if(res_df[number, pval] <= limit_filter & abs(res_df[number, pval + 1]) >= limit_e){
      signig_filter[number] <- T
    }
  }
}

#############
#############Shorten names of microbes while keeping taxonomic order
#############

res_df$microbe <- sub(".*ales_", "", res_df$microbe )

#############
#############Split effect size and p-value into two different data.frames and define when to add a "star" for significance
#############
plot.m <- melt(res_df[signig_filter,c(1, cols_containing_evals)], id.vars = "microbe")

pval<- melt(res_df[signig_filter,c(1, cols_containing_pvals)], id.vars = "microbe")


stars <- pval$value
for(number in 1:nrow(pval)){
  if(pval$value[number] <= 0.1){
    stars[number] <- "*"
  }
  else{
    stars[number] <- ""
  }
  if(pval$value[number] <= 0.01){
    stars[number] <- "**"
  }
  if(pval$value[number] <= 0.001){
    stars[number] <- "***"
  }
}

plot.m$value <- as.numeric(plot.m$value)

plot.m$stars <- stars


plot.m$microbe <- sub("*; D_5__", " ", plot.m$microbe)
plot.m$variable = factor(plot.m$variable, levels = c("Pre Young yFMT vs Post Young yFMT", 
                                                     "Pre Aged oFMT vs Post Aged oFMT", 
                                                     "Pre Aged yFMT vs Post Aged yFMT"))


p <- ggplot(plot.m, aes(variable, microbe), height=20, width=20) + geom_tile(aes(fill = value), colour = "white") + 
  scale_fill_gradientn(colours = c("darkblue"   , "darkblue"  ,
                                   "blue"       , "blue"      , "blue"   ,
                                   "white"      , "white"     , "white"  ,
                                   "orange"     , "orange"    , "orange" ,
                                   "orange"     , "orange"    , 
                                   "orange"     ,  "orange"   ,    
                                   "darkorange" , "darkorange", 
                                   "red"        , "red"       ), limits = c(-1.5, 3)) +
  geom_text(aes(label=stars)) + xlab("") + ylab("") + labs(fill = "Effect Size") +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 330, hjust = 0, vjust = 0.5,  size = 10),
        axis.text.y = element_text(size = 10)) +  
  scale_y_discrete(position = "right") + 
  scale_x_discrete(labels=c("Delta Young yFMT", "Delta Aged oFMT", "Delta Aged yFMT"))

p



GBM_df <- pairwise_DA_wrapper(reads = GBM, 
                              groups = paste(metadata$Timepoint, metadata$Legend, sep = " "), 
                              comparisons = comp,
                              paired.test = T)

res_df <- GBM_df
row.names(res_df) <- res_df$microbe

cols_containing_pvals <- 1:((ncol(res_df))/2)*2
cols_containing_evals <- (1:((ncol(res_df))/2)*2) +1
change_NA_to_1      <- function(x) { replace(x, is.na(x), 1) }
change_NA_to_0      <- function(x) { replace(x, is.na(x), 0) }

res_df[,cols_containing_pvals] <- change_NA_to_1(res_df[,cols_containing_pvals])
res_df[,cols_containing_evals] <- change_NA_to_0(res_df[,cols_containing_evals])

#############
#############Set up filter to pick displayed microbes based on p-values
#############
View(head(res_df))

signig_filter = rep(F, nrow(res_df))
limit_filter = 0.1
limit_e      = 0.5
for(number in (1:nrow(res_df))){
  for(pval in (cols_containing_pvals)){
    if(res_df[number, pval] <= limit_filter & abs(res_df[number, pval + 1]) >= limit_e){
      signig_filter[number] <- T
    }
  }
}

#############
#############Shorten names of microbes while keeping taxonomic order
#############

res_df$microbe <- sub(".*ales_", "", res_df$microbe )

#############
#############Split effect size and p-value into two different data.frames and define when to add a "star" for significance
#############
plot.m <- melt(res_df[signig_filter,c(1, cols_containing_evals)], id.vars = "microbe")

pval<- melt(res_df[signig_filter,c(1, cols_containing_pvals)], id.vars = "microbe")


stars <- pval$value
for(number in 1:nrow(pval)){
  if(pval$value[number] <= 0.1){
    stars[number] <- "*"
  }
  else{
    stars[number] <- ""
  }
  if(pval$value[number] <= 0.01){
    stars[number] <- "**"
  }
  if(pval$value[number] <= 0.001){
    stars[number] <- "***"
  }
}

plot.m$value <- as.numeric(plot.m$value)

plot.m$stars <- stars


plot.m$microbe <- sub("*; D_5__", " ", plot.m$microbe)
plot.m$variable = factor(plot.m$variable, levels = c("Pre Young yFMT vs Post Young yFMT", 
                                                     "Pre Aged oFMT vs Post Aged oFMT", 
                                                     "Pre Aged yFMT vs Post Aged yFMT"))


pGBM <- ggplot(plot.m, aes(variable, microbe), height=20, width=20) + geom_tile(aes(fill = value), colour = "white") + 
  scale_fill_gradientn(colours = c("darkblue"   , "darkblue"  ,
                                   "blue"       , "blue"      , "blue"   ,
                                   "white"      , "white"     , "white"  ,
                                   "orange"     , "orange"    , "orange" ,
                                   "orange"     , "orange"    , 
                                   "orange"     ,  "orange"   ,    
                                   "darkorange" , "darkorange", 
                                   "red"        , "red"       ), limits = c(-1.5, 3)) +
  geom_text(aes(label=stars)) + xlab("") + ylab("") + labs(fill = "Effect Size") +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 330, hjust = 0, vjust = 0.5,  size = 10),
        axis.text.y = element_text(size = 10)) +  
  scale_y_discrete(position = "right") + 
  scale_x_discrete(labels=c("Delta Young yFMT", "Delta Aged oFMT", "Delta Aged yFMT"))

pGBM










anova(lm(alphadiv$chao1 ~ metadata$Timepoint + metadata$Legend + metadata$Timepoint*metadata$Legend))
anova(lm(alphadiv$asymptotic_simps ~ metadata$Timepoint + metadata$Legend + metadata$Timepoint*metadata$Legend))
anova(lm(alphadiv$asymptotic_shannon ~ metadata$Timepoint + metadata$Legend + metadata$Timepoint*metadata$Legend))



fig1 <- (b |(chao /shan/simps)) / 
  (vola | ent)/
  (p|pGBM) & theme(text = element_text(size = 12))

fig1 + 
  plot_annotation(
    tag_levels = "A")  

alphadiv




(b | (chao / shan / simps)) /
  (vola | ent)/
  (p + guides(fill = "none") | (guide_area()/  pGBM )+ plot_layout(guides = "collect") ) 





(b | (chao / shan / simps)) /
  (vola | ent)/
  (p + guides(fill = "none") | ( (plot_spacer() /pGBM  + guides(fill = "none")))) 



layout <- "
AAAABB##
AAAABB##
AAAACCDD
AAAACCDD
"

b + chao  + simps + shan +#+ vola + ent + p + guides(fill = "none") + pGBM + guides(fill = "none") +
  plot_layout(design = layout)   + 
  plot_annotation(
    tag_levels = "A")  & theme(text = element_text(size = 12)) 





####This is the one for now
timeline <- image_read_svg('Fig_1_timeline.svg', width = 1460)

ggtimeline = ggplot() + 
  draw_image(timeline) + theme_nothing()  + 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x = NULL, y = NULL) 

top = ((((b + draw_image(magick::image_read_pdf('for_pca_legend.pdf', density = 600), x = 14, y = 13,  width = 13, height = 13) + labs(tag = "b", title = "Principal Component Analysis")) |(((chao + labs(tag = "c")) /simps/shan))) + plot_layout(widths = c(3, 1)) & theme(text = element_text(size = 12)))) #+ plot_annotation(tag_levels = "A")

  
  
  bottom = p  +  theme(legend.position = c(4, 0.2), legend.background = element_rect(colour = "black")) + labs(tag = "d", title = "Differential Abundance of Microbial Genera")  | ((ent + labs(tag = "e", title = "Enterococcus Colonizes Aged Microbiome after FMT"))
                                           /
                                             (pGBM + guides(fill = 'none') + labs(tag = "f", title = "Differential Abundance of GBMs")  + vola + labs(tag = "g", title = "Microbiome Change During Treatment")) )   & theme(text = element_text(size = 12))   
  
bottom
fig1 = (ggtimeline + labs(tag = "a", title = "Timeline of Studies")) / top / bottom  & theme(plot.title = element_text(face = "bold"), 
                                       plot.tag = element_text(face = "bold"))

fig1
