#source("http://bioconductor.org/biocLite.R")
#options(scipen=1) #To get rid of scientific notation in ggplot2
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
library(dplyr)
library(tidyverse)
library(qvalue)
library(omixerRpm)
library(goeveg)
library(scales)
library("gridExtra")
library(cowplot)
library("genefilter")
library(patchwork)
library(stringr) 
library("DescTools")
#library(devtools)
#devtools::install_github(repo = "thomazbastiaanssen/Tjazi")

options(stringsAsFactors = F)
getwd()
setwd("/home/thomaz/Documents/PhD/FMT aging/")

metabol <- read.delim("Metabolomics/FMT_Aging _Metabolomics_All_t.csv", sep = ",", row.names = 1)
metabol = t(metabol) #legacy reasons
metadata  <- read.delim("metadata_Pool3_included.csv", sep = ",")
intersect(metadata$Metabolon_ID, rownames(metabol))

metadata <- metadata[metadata$Timepoint == "Post"  & metadata$Legend != "Aged CTR" & metadata$Legend != "Young CTR",]
metadata <- metadata[metadata$Metabolon_ID %in% row.names(metabol) ,]
metadata$KW_legend <- metadata$Legend
metadata$KW_legend[metadata$KW_legend != "Aged oFMT"] <- "All yFMT"
metabol <- metabol[metadata$Metabolon_ID,]



#cvs <- apply(X = metabol, MARGIN = 2, FUN = goeveg::cv, na.rm = T)

#metabol = metabol[,cvs < 0.25]
dim(metabol)
colnames(metabol)
NAs <- apply(metabol, 2, function(x) sum(is.na(x)))

#metabol = metabol[,NAs < nrow(metabol)/10]
metabol[is.na(metabol)] <- floor(min(metabol, na.rm = T)*0.95)
dim(metabol)
#metabol <- t(varFilter(t(metabol), var.func=IQR, var.cutoff=0.150, filterByQuantile=T))



metabol    <- metabol[metadata$Metabolon_ID,]

#get_asymptotic_alpha(counts = species, metadata = paste(metadata$Timepoint, metadata$Legend, sep = " "))
species   <- data.frame(t(metabol))


species   <- apply(species,c(1,2),function(x) as.numeric(as.character(x)))

dim(species)
#species  <- t(species)
####
####Calculate Aitchison distance and prepare for PCA
####
conds       <- c(rep("A", ncol(species)-10 ), rep("B", 10)) #If you have less than 12 animals, adjust!
species.clr <- aldex.clr(species, conds, mc.samples = 1000, denom="all", verbose=TRUE, useMC=TRUE) 
species.eff <- aldex.effect(species.clr, verbose = TRUE, include.sample.summary = TRUE)

species.exp <- (species.eff[,c(4:(ncol(species.eff)-4))]) #remove the useless t-test-like results

data.a.pca  <- prcomp(t(species.exp))


loadings <- as.data.frame(data.a.pca$rotation) 
loadings.sub = loadings[which(abs(loadings$PC1)>0.2 | abs(loadings$PC2)>0.2),]

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
pca$Legend               = factor(metadata$Legend, levels = c("Young yFMT", "Aged oFMT", "Aged yFMT"))
pca$Timepoint            = metadata$Timepoint
pca$group                = metadata$Legend



a = ggplot() + 
  #geom_line(size = 1.5) +
  stat_ellipse(data = pca, aes(group = Legend, x=PC1, y=PC3, fill = Legend),alpha = 1/4,  geom = "polygon")+
  geom_point(size=4, stroke = 1, shape = 21, data = pca, aes(x=PC1, y=PC3, fill = Legend, group = group)) +
  geom_label_repel()+
  xlab(paste("PC1: ", pc1,  "%", sep="")) + 
  ylab(paste("PC3: ", pc3,  "%", sep="")) + 
  theme_bw() + 
  scale_fill_manual( values = c("Young yFMT" = "#2166ac",
                                "Aged oFMT"  = "#b2182b",
                                "Aged yFMT"  = "#ef8a62",
                                "oFMT"       = "#b2182b", 
                                "yFMT"       = "#386cb0",
                                "Young CTR"  = "#b3b3b3", 
                                "Aged CTR"   = "#666666"
  )) + ggtitle("Principal Component Analysis")  + 
  theme(legend.position = c(0.85, 0.1), legend.background = element_rect(colour = "black"))
#geom_segment(data = loadings.sub, aes(x=0,y=0,xend=PC1*50,yend=PC2*50),arrow=arrow(length=unit(0.1,"cm")), color = "black")+
#geom_text(data = loadings.sub, aes(x=PC1*50, y=PC2*50, label=row.names(loadings.sub)),color="black") 

# scale_shape_manual(values = c("1" = 21, 


a



out_df <- pairwise_DA_wrapper(reads          = species, 
                              groups         = metadata$Legend, 
                              comparisons    = data.frame(a = c("Aged oFMT",  "Young yFMT", "Young yFMT"),
                                                          b = c("Aged yFMT", "Aged yFMT"  , "Aged oFMT")), 
                              parametric     = F, 
                              ignore.posthoc = T)#, "Aged yFMT",  "Aged oFMT")), parametric = T)

res_df_mw <- pairwise_DA_wrapper(reads = species, 
                                 groups = metadata$KW_legend, comparisons = data.frame(a = c("Aged oFMT"),
                                                                                       b =  c("All yFMT")), 
                                 parametric = F, 
                                 ignore.posthoc = T)


res_df = cbind(out_df, res_df_mw[,-1])


#View(res_df)
res_df$`Aged oFMT vs Aged yFMT q.value`  <- qvalue(p = res_df$`Aged oFMT vs Aged yFMT p.value` )$qvalues
res_df$`Young yFMT vs Aged yFMT q.value` <- qvalue(p = res_df$`Young yFMT vs Aged yFMT p.value`)$qvalues
res_df$`Young yFMT vs Aged oFMT q.value` <- qvalue(p = res_df$`Young yFMT vs Aged oFMT p.value`)$qvalues
res_df$`Aged oFMT vs All yFMT q.value` <- qvalue(p = res_df$`Aged oFMT vs All yFMT p.value`)$qvalues
#res_df$`Aged CTR vs Young CTR q.value`   <- qvalue(p = res_df$`Aged CTR vs Young CTR p.value`)$qvalues


#up_in_young <- res_df[((res_df$`Aged CTR vs Young CTR`) >= 0.95),]$microbe 
#up_in_aged  <- res_df[((res_df$`Aged CTR vs Young CTR`) <= -0.95),]$microbe 
#dif_in_ctr  <- res_df[abs(res_df$`Aged CTR vs Young CTR`) >= 0.95,]$microbe


res_df$`Aged oFMT vs Aged yFMTpq`        <- rep("ns", nrow(res_df))
res_df$`Aged oFMT vs Aged yFMTpq`[res_df$`Aged oFMT vs Aged yFMT p.value` < 0.05] <- "p"
res_df$`Aged oFMT vs Aged yFMTpq`[res_df$`Aged oFMT vs Aged yFMT p.value` < 0.05 &
                                    res_df$`Aged oFMT vs Aged yFMT q.value`< 0.2 ] <- "pq"

res_df$`Young yFMT vs Aged yFMTpq`        <- rep("ns", nrow(res_df))
res_df$`Young yFMT vs Aged yFMTpq`[res_df$`Young yFMT vs Aged yFMT p.value` < 0.05] <- "p"
res_df$`Young yFMT vs Aged yFMTpq`[res_df$`Young yFMT vs Aged yFMT p.value` < 0.05 &
                                     res_df$`Young yFMT vs Aged yFMT q.value`< 0.2 ] <- "pq"



res_df$`Young yFMT vs Aged oFMTpq`        <- rep("ns", nrow(res_df))
res_df$`Young yFMT vs Aged oFMTpq`[res_df$`Young yFMT vs Aged oFMT p.value` < 0.05] <- "p"
res_df$`Young yFMT vs Aged oFMTpq`[res_df$`Young yFMT vs Aged oFMT p.value` < 0.05 &
                                     res_df$`Young yFMT vs Aged oFMT q.value`< 0.2 ] <- "pq"


res_df$`Aged oFMT vs All yFMTpq`        <- rep("ns", nrow(res_df))
res_df$`Aged oFMT vs All yFMTpq`[res_df$`Aged oFMT vs All yFMT p.value` < 0.05] <- "p"
res_df$`Aged oFMT vs All yFMTpq`[res_df$`Aged oFMT vs All yFMT p.value` < 0.05 &
                                   res_df$`Aged oFMT vs All yFMT q.value`< 0.2 ] <- "pq"



#View(res_df)


metab_names <- read.delim("metadata_metabolomics.csv", sep = ",")
res_df$microbe



res_df$`Metabolite Type` <- metab_names$SUPER_PATHWAY
res_df$label <- metab_names$Full_abbrev


reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}






metab_agedo_vs_agedy <- 
  ggplot()+
  scale_y_continuous(trans=reverselog_trans(10)) +
  scale_x_continuous(limits = c(-1.25, 1.25))+
  scale_alpha_manual(values = c("ns" = 0.5,
                                "p"  = 0.75, 
                                "pq" = 1), )+
  scale_size_manual(values =  c("ns" = 3,
                                "p"  = 3, 
                                "pq" = 3.5))+
  scale_discrete_manual(    aesthetics = "stroke",
                            values =  c("ns" = 0.5,
                                        "p"  = 0.5, 
                                        "pq" = 1))+
  scale_fill_brewer(palette = "Accent")+
  geom_point(data = res_df, shape = 21,  
             aes(
               x       = `Aged oFMT vs Aged yFMT`, 
               y       = `Aged oFMT vs Aged yFMT p.value`, 
               fill    = `Metabolite Type`, 
               alpha   = `Aged oFMT vs Aged yFMTpq`, 
               size    = `Aged oFMT vs Aged yFMTpq`,
               stroke  = `Aged oFMT vs Aged yFMTpq`)) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", alpha = 0.75)+
  theme_bw() + theme(axis.text.y = element_text(size = 5)) +
  ylab("") + xlab("Effect Size") + 
  guides(alpha = FALSE, size = FALSE, stroke = FALSE, fill = guide_legend(override.aes = list(size= 4)))

metab_agedo_vs_agedy

metab_agedo_vs_agedy = metab_agedo_vs_agedy +  
  geom_text_repel(data = subset(res_df, `Aged oFMT vs Aged yFMT p.value` < 0.05 & 
                                  (`Aged oFMT vs Aged yFMT` >  quantile((res_df$`Aged oFMT vs Aged yFMT`), 95/100) |
                                     `Aged oFMT vs Aged yFMT` <  quantile((res_df$`Aged oFMT vs Aged yFMT`), 5/100) )), 
                  aes(  x       = `Aged oFMT vs Aged yFMT`, 
                        y       = `Aged oFMT vs Aged yFMT p.value`, 
                        label   = label), 
                  size    = 3, box.padding = 0.5) + ggtitle("Aged oFMT vs Aged yFMT")

metab_agedo_vs_agedy






metab_youngy_vs_agedy <- 
  ggplot()+
  scale_y_continuous(trans=reverselog_trans(10)) +
  scale_x_continuous(limits = c(-2.5, 5.5))+
  scale_alpha_manual(values = c("ns" = 0.5,
                                "p"  = 0.75, 
                                "pq" = 1), )+
  scale_size_manual(values =  c("ns" = 3,
                                "p"  = 3, 
                                "pq" = 3.5))+
  scale_discrete_manual(    aesthetics = "stroke",
                            values =  c("ns" = 0.5,
                                        "p"  = 0.5, 
                                        "pq" = 1))+
  scale_fill_brewer(palette = "Accent")+
  geom_point(data = res_df, shape = 21,  
             aes(
               x       = `Young yFMT vs Aged yFMT`, 
               y       = `Young yFMT vs Aged yFMT p.value`, 
               fill    = `Metabolite Type`, 
               alpha   = `Young yFMT vs Aged yFMTpq`, 
               size    = `Young yFMT vs Aged yFMTpq`,
               stroke  = `Young yFMT vs Aged yFMTpq`)) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", alpha = 0.75)+
  theme_bw() + theme() +
  ylab("p-value") + xlab("Effect Size") + 
  guides(alpha = FALSE, size = FALSE, stroke = FALSE, fill = guide_legend(override.aes = list(size= 4)))

metab_youngy_vs_agedy

metab_youngy_vs_agedy = metab_youngy_vs_agedy +  
  geom_text_repel(data = subset(res_df, `Young yFMT vs Aged yFMT p.value` < 0.05 &   
                                  (`Young yFMT vs Aged yFMT` >  quantile((res_df$`Young yFMT vs Aged yFMT`), 95/100) |
                                     `Young yFMT vs Aged yFMT` <  quantile((res_df$`Young yFMT vs Aged yFMT`), 5/100)  )),
                  
                  
                  
                  aes(  x       = `Young yFMT vs Aged yFMT`, 
                        y       = `Young yFMT vs Aged yFMT p.value`, 
                        label   = label), 
                  size    = 3, box.padding = 0.5)  + ggtitle("Young yFMT vs Aged yFMT")

metab_youngy_vs_agedy








metab_youngy_vs_agedo <- 
  ggplot()+
  scale_y_continuous(trans=reverselog_trans(10)) +
  scale_x_continuous(limits = c(-2, 5))+
  scale_alpha_manual(values = c("ns" = 0.5,
                                "p"  = 0.75, 
                                "pq" = 1), )+
  scale_size_manual(values =  c("ns" = 3,
                                "p"  = 3, 
                                "pq" = 3.5))+
  scale_discrete_manual(    aesthetics = "stroke",
                            values =  c("ns" = 0.5,
                                        "p"  = 0.5, 
                                        "pq" = 1))+
  scale_fill_brewer(palette = "Accent")+
  geom_point(data = res_df, shape = 21,  
             aes(
               x       = `Young yFMT vs Aged oFMT`, 
               y       = `Young yFMT vs Aged oFMT p.value`, 
               fill    = `Metabolite Type`, 
               alpha   = `Young yFMT vs Aged oFMTpq`, 
               size    = `Young yFMT vs Aged oFMTpq`,
               stroke  = `Young yFMT vs Aged oFMTpq`)) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", alpha = 0.75)+
  theme_bw() + theme(axis.text.y = element_text(size = 5)) +
  ylab("p-value") + xlab("Effect Size") + 
  guides(alpha = FALSE, size = FALSE, stroke = FALSE, fill = guide_legend(override.aes = list(size= 4)))

metab_youngy_vs_agedo

metab_youngy_vs_agedo = metab_youngy_vs_agedo +  
  geom_text_repel(data = subset(res_df, `Young yFMT vs Aged oFMT p.value` < 0.05 & 
                                  (`Young yFMT vs Aged oFMT` >  quantile((res_df$`Young yFMT vs Aged oFMT`), 95/100) |
                                     `Young yFMT vs Aged oFMT` <  quantile((res_df$`Young yFMT vs Aged oFMT`), 5/100)  )), 
                  aes(  x       = `Young yFMT vs Aged oFMT`, 
                        y       = `Young yFMT vs Aged oFMT p.value`, 
                        label   = label), 
                  size    = 3, box.padding = 0.5)  + ggtitle("Young yFMT vs Aged oFMT")

metab_youngy_vs_agedo




metab_reverse_aging <- 
  ggplot()+
  scale_y_continuous(trans=reverselog_trans(10)) +
  scale_x_continuous(limits = c(-1.50, 1.25))+
  scale_alpha_manual(values = c("ns" = 0.5,
                                "p"  = 0.75, 
                                "pq" = 1), )+
  scale_size_manual(values =  c("ns" = 3,
                                "p"  = 3, 
                                "pq" = 3.5))+
  scale_discrete_manual(    aesthetics = "stroke",
                            values =  c("ns" = 0.5,
                                        "p"  = 0.5, 
                                        "pq" = 1))+
  scale_fill_brewer(palette = "Accent")+
  geom_point(data = res_df, shape = 21,  
             aes(
               x       = `Aged oFMT vs All yFMT`, 
               y       = `Aged oFMT vs All yFMT p.value`, 
               fill    = `Metabolite Type`, 
               alpha   = `Aged oFMT vs All yFMTpq`, 
               size    = `Aged oFMT vs All yFMTpq`,
               stroke  = `Aged oFMT vs All yFMTpq`)) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", alpha = 0.75)+
  theme_bw() + theme(aspect.ratio = 1) +
  ylab("p-value") + xlab("Effect Size") + 
  guides(alpha = FALSE, size = FALSE, stroke = FALSE, fill = guide_legend(override.aes = list(size= 4)))

metab_reverse_aging

metab_reverse_aging = metab_reverse_aging +  
  geom_text_repel(data = subset(res_df, `Aged oFMT vs All yFMT p.value` < 0.05 & 
                                  (`Aged oFMT vs All yFMT` >  quantile((res_df$`Aged oFMT vs All yFMT`), 95/100) |
                                     `Aged oFMT vs All yFMT` <  quantile((res_df$`Aged oFMT vs All yFMT`), 5/100) )), 
                  aes(  x       = `Aged oFMT vs All yFMT`, 
                        y       = `Aged oFMT vs All yFMT p.value`, 
                        label   = label), 
                  size    = 3) + ggtitle("Aged oFMT vs All yFMT")

metab_reverse_aging





#species.exp[inter,]

#View(metab_names)

labs_res_df <- subset(res_df, (`Aged oFMT vs Aged yFMT p.value` < 0.05 &   
                        (`Aged oFMT vs Aged yFMT` >  quantile((res_df$`Aged oFMT vs Aged yFMT`), 95/100) |
                           `Aged oFMT vs Aged yFMT` <  quantile((res_df$`Aged oFMT vs Aged yFMT`), 5/100)  ))
                     #   |(`Aged oFMT vs All yFMT p.value` < 0.05 & 
                    #      (`Aged oFMT vs All yFMT` >  quantile((res_df$`Aged oFMT vs All yFMT`), 95/100) |
                    #         `Aged oFMT vs All yFMT` <  quantile((res_df$`Aged oFMT vs All yFMT`), 5/100) ))
                        |(`Young yFMT vs Aged oFMT p.value` < 0.05 & 
                        (`Young yFMT vs Aged oFMT` >  quantile((res_df$`Young yFMT vs Aged oFMT`), 95/100) |
                           `Young yFMT vs Aged oFMT` <  quantile((res_df$`Young yFMT vs Aged oFMT`), 5/100)  ))
                          )


dim(labs_res_df)

ggplot(data = metab_names[metab_names$Full_abbrev %in%labs_res_df$label, ], aes(x = SUPER_PATHWAY, y =Super_number, fill = SUPER_PATHWAY )) + geom_point(shape = 21)


#View( metab_names[metab_names$Full_abbrev %in%labs_res_df$label, ])


gg_leg = metab_names[metab_names$Full_abbrev %in% labs_res_df$label, ]
gg_leg = gg_leg[order(paste(gg_leg$SUPER_PATHWAY, gg_leg$Full_abbrev)),]
dim(gg_leg)

gg_leg$y = rep(nrow(gg_leg):1, 1)
gg_leg$x = rep(1, each = nrow(gg_leg))

gg_leg[gg_leg$Metabolite == "fructose 1,6-diphosphate/glucose 1,6-diphosphate/myo-inositol diphosphates",]$Metabolite = "fructose 1,6-/glucose 1,6-/myo-inositol - diphosphate"

gg_leg$test_labels = paste(gg_leg$Full_abbrev, gg_leg$Metabolite, sep = " ")


legplot = ggplot(data = gg_leg, aes(x = x, y = y, fill = SUPER_PATHWAY)) + 
  geom_point(shape= 22, size = 3) + 
  geom_text(aes(label = test_labels),nudge_x = 0.003, hjust =0, vjust = 0.5, size = 2.5)+  theme_void() + guides(fill = "none") + scale_x_continuous(limits = c(1, 1.15)) +  
  scale_fill_manual(values = c("Amino Acid"   = "#7fc97f", 
                                                                                                              "Carbohydrate" = "#beaed4", 
                                                                                                              "Cofactors and Vitamins" = "#fdc086", 
                                                                                                              "Energy" = "#ffff99", 
                                                                                                              "Lipid"  = "#386cb0", 
                                                                                                              "Nucleotide" = "#f0027f", 
                                                                                                              "Peptide"  = "#bf5b17", 
                                                                                                              "Xenobiotics" = "#666666")) #+ 


legplot


#(metab_youngy_vs_agedo  | metab_agedo_vs_agedy & guides(fill = "none")) | legplot   #+ plot_layout(guides = 'collect')



(metab_youngy_vs_agedo  +  metab_agedo_vs_agedy   & guides(fill = "none")) + legplot  + theme(panel.spacing.y =  unit(-5, "pt")) + plot_layout(widths = c(1, 1, 0.4)) + labs(caption = paste0("Produced by Thomaz ", Sys.time()))

inter <- ((res_df[res_df$`Aged oFMT vs All yFMT q.value`   < 0.2 & 
                  res_df$`Aged oFMT vs Aged yFMT p.value`  < 0.1 & 
                  res_df$`Young yFMT vs Aged oFMT p.value` < 0.1,]))$microbe




layout ="
######CC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
AAABBBCC
######CC"


top = (metab_youngy_vs_agedo + theme(legend.position = c(0.8275, 0.178), legend.background = element_rect(colour = "black"))+ labs(tag = "a") +  
         guides(fill = guide_legend(override.aes = list(stroke= 1, size = 3.5)))+ 
         metab_agedo_vs_agedy + guides(fill = 'none') + labs(tag = "b") + 
         theme(
           axis.ticks.y = element_blank(),
           axis.title.y = element_blank()) ) + 
  legplot + 
  plot_layout(design = layout)


top


firstup <- function(x) {
  substr(x, which((unlist(strsplit(c(x), ""))) %in% c(letters, LETTERS))[1],
         which((unlist(strsplit(c(x), ""))) %in% c(letters, LETTERS))[1]) <- toupper(substr(x, 
                                    which((unlist(strsplit(c(x), ""))) %in% c(letters, LETTERS))[1], 
                                    which((unlist(strsplit(c(x), ""))) %in% c(letters, LETTERS))[1]
))
  x
} ###Remember to use sapply on this one

gg_boxplots = data.frame(t(species.exp[inter,]))
colnames(gg_boxplots)  = sapply(X = annot_ao_vs_ay$microbe, FUN = firstup )
gg_boxplots$Legend = metadata$Legend

gg_box = gg_boxplots %>%
  pivot_longer(!Legend)
gg_box$Legend = factor(gg_box$Legend, levels = c("Young yFMT", "Aged oFMT", "Aged yFMT"))
gg_box

break_fct <-  function(x){
  round(seq(min(x), max(x), length = 5))
}



annot_ao_vs_ay <- res_df[res_df$microbe %in% inter,]
annot_ao_vs_ay$name = sapply(X = annot_ao_vs_ay$microbe, FUN = firstup )
annot_ao_vs_ay$stars = "ERROR"
annot_ao_vs_ay$x = 2.5
annot_ao_vs_ay$y = Inf
annot_ao_vs_ay$stars[annot_ao_vs_ay$`Aged oFMT vs Aged yFMT p.value` < 0.1] <- "#"
annot_ao_vs_ay$stars[annot_ao_vs_ay$`Aged oFMT vs Aged yFMT p.value` < 0.05] <- "*"
annot_ao_vs_ay$stars[annot_ao_vs_ay$`Aged oFMT vs Aged yFMT p.value` < 0.01] <- "**"
annot_ao_vs_ay$stars[annot_ao_vs_ay$`Aged oFMT vs Aged yFMT p.value` < 0.001] <- "***"
annot_ao_vs_ay$stars



annot_ao_vs_y <- res_df[res_df$microbe %in% inter,]
annot_ao_vs_y$name = sapply(X = annot_ao_vs_ay$microbe, FUN = firstup )
annot_ao_vs_y$stars = "ERROR"
annot_ao_vs_y$x = 1.5
annot_ao_vs_y$y = Inf
annot_ao_vs_y$stars[annot_ao_vs_y$`Young yFMT vs Aged oFMT p.value` < 0.1] <- "#"
annot_ao_vs_y$stars[annot_ao_vs_y$`Young yFMT vs Aged oFMT p.value`  < 0.05] <- "*"
annot_ao_vs_y$stars[annot_ao_vs_y$`Young yFMT vs Aged oFMT p.value`  < 0.01] <- "**"
annot_ao_vs_y$stars[annot_ao_vs_y$`Young yFMT vs Aged oFMT p.value`  < 0.001] <- "***"
annot_ao_vs_y$stars



bottom = ggplot(data = gg_box) + 
  geom_violin(aes(group =Legend, x = Legend, y = value, fill = Legend ), alpha = 1/2,  col = "black", draw_quantiles = 0.5)+ 
  #stat_boxplot(coef = 100000,  aes(fill = Legend, alpha = 1/4))+ 
  geom_dotplot( aes(x = Legend, y = value, fill = Legend), 
               binaxis  = 'y',  
               stackdir = "center",  
               dotsize  = 2.5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = expand_scale(mult = c(0.15, .4))) +
  theme_bw() +   
  scale_fill_manual( values = c("Aged oFMT"  = "#b2182b",
                                "Aged yFMT"  = "#ef8a62",
                                "Young yFMT" = "#2166ac")) + 
  facet_wrap(.~name, scales = "free_y", nrow = 7,) + guides(fill = "none", 
                                                            alpha = "none") +
  ylab("Metabolite concentration CLR(pg/mL)") + xlab("") + 
  ggtitle("Metabolites Altered by Ageing and Reversed by yFMT") +
  theme(axis.text.x  = element_text(size = 6),
                              strip.text.x = element_text(size = 6.8))+
  geom_text(data = annot_ao_vs_ay, 
            mapping = aes(x = x, y = y, label = stars), vjust = 2) +
  geom_text(data = annot_ao_vs_y, 
            mapping = aes(x = x, y = y, label = stars), vjust = 2) + 
  geom_text(data = annot_ao_vs_ay, 
            mapping = aes(x = Inf, y = -Inf, 
                          label = paste("q = ", round(annot_ao_vs_ay$`Aged oFMT vs All yFMT q.value`, digits = 3) )), hjust = 1, vjust = -0.25, size = 2) 



bottom  




metaban = read.delim("pathway_results.csv", sep = ",")
metaban$X[3] =   "Alanine, aspartate and\nglutamate metabolism"

metaboanalyst = ggplot(data = metaban, aes(x    = Impact, 
                           y    = X.log.p., 
                           size = Impact, 
                           alpha = FDR < 0.2 ,
                           label = X)) +
  geom_point(shape = 21, fill = "#0c2c84", stroke = 1.1) + theme_bw() + 
  geom_text_repel(data = subset(metaban,FDR < 0.2 ), size = 3, box.padding =1, direction = "y", seed = 1)+
  scale_alpha_manual(values = c("TRUE"  = 1, 
                                "FALSE" = 0.5)) + guides(alpha = "none", 
                                                         size = "none") +
  ylab("-log10(p-value)") + 
  xlab("Pathway Impact Score") +
  ggtitle("Metabolomic Pathway Analysis")






#top / 
#  ((((a + guides(fill = "none"))/metaboanalyst) | bottom) + plot_layout(widths = c(2, 5))) + plot_layout(heights = c(2, 3)) 


fig3 = top / 
  ((((a + labs(tag = "c"))/metaboanalyst + labs(tag = "e")) | bottom + labs(tag = "d")) + 
     plot_layout(widths = c(2, 5))) + 
  plot_layout(heights = c(2, 3)) & theme(plot.title = element_text(face = "bold"), 
                                         plot.tag = element_text(face = "bold"))

fig3
###width  1460
###height 2000
write.csv(res_df, file = "stats_metabolomics_and_info_for_figures.csv", quote = T)
write.csv(gg_box, file = "data_used_for_boxplots.csv", quote = T)

