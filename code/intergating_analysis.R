# Intergating analysis for longevity
# On DEC 14 2023
# By Jun
#====================================
library(data.table)
# Basic information about samples of each dataset
# merge sample_id.csv and run_id.csv
s_id <- fread('sample_id.csv')
r_id <- fread('run_id.csv')

sra_sampleid <- merge(s_id,r_id, by = "sample_id")

sra_sampleid$age <- as.numeric(sra_sampleid$age)

# GROUP SAMPLE BASED ON AGE 
sra_sampleid <- sra_sampleid %>% mutate(category = case_when(age < 62  ~ "W",
                                                             age < 89.5 & age > 62 ~ "R",
                                                             age >89.5  ~ "C",
                                                             TRUE ~"Not listed"))

write.csv(sra_sampleid, file = "sra_sampleid.csv")


#-------------------------------------------------------------------------------
# plot 1.0sample age_gender
sra_sampleid <- fread("sra_sampleid.csv")


# plot using ggplot
g_age_gender <- sra_sampleid[,c("gender","country","category")]


df2 <- g_age_gender %>% group_by(gender,country,category) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

df2_m<- df2[which(df2$gender %in% "M"),]
df2_f <- df2[which(df2$gender %in% "F"),]
df2_f$total_count <- df2_f$total_count *-1


# Male population data graph
pdf(paste("F1.age_gender", "pdf", sep="."),width=12, height=9)
ggplot()+
  geom_bar(data =df2_m,aes(category,total_count,fill=country),stat = "identity",position="stack",color="black",width = 0.7,size=0.25)+
  
  # Famale population data graph
  
  
  geom_bar(data =df2_f,aes(category,y= total_count *(-1),fill=country),stat = "identity",position="stack",color="black",width = 0.7,size=0.25)+
  coord_flip()+
  scale_y_continuous(limits = c(-275,125),
                     breaks = seq(-275,125,25),
                     labels = abs(seq(-275,125,25)))+
  labs(x = "", y = "")+
  scale_fill_manual(values=c("#C2C2C2","#d9544d","#a6cee3","#b2df8a"),label = c("China","Italy_1","Italy_2","Japan"))+
  # display themes
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 24, face = "bold"),  # fonts/sizes
    axis.title = element_text(size = 24, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5), 
    plot.caption = element_text(hjust=0, size=11, face = "italic")
  )
dev.off()
#====================================
# prepare taxonomy information for grouping by bacteria and fungal virus 
prepareDatabase('accessionTaxa.sql') # only need to run once 


# obtain microbes profile of each sample

setwd("D:/OneDrive/dataanalysis/longevity/data/kraken")

kraken_files <- list.files(full.names = T, recursive = FALSE)

file_names = ""
for (name in kraken_files) {
  file_names <- c(file_names, unlist(strsplit(name, split='/', fixed=TRUE))[2])}
file_names <- file_names[-1]

# READ IN FILES
# loading the kraken table
y <- 0
for (x in file_names) {
  y <- y + 1
  if (y == 1) {
    kraken_table <- read.table(file = x, header = T, quote = "", sep = "\t")
    kraken_table <- kraken_table[,c(2,6)]
    n <- str_remove_all(x, "_d_k2_brackentaxa.tsv")
    colnames(kraken_table)[2] <- n
  }
  if (y > 1) {
    temp_table <- read.table(file = x, header = T, quote = "", sep = "\t")
    temp_table <- temp_table[,c(2,6)]
    n <- str_remove_all(x, "_d_k2_brackentaxa.tsv")
    colnames(temp_table)[2] <- n
    kraken_table <- merge(kraken_table, temp_table, by = "taxonomy_id", all = T)  }
}
kraken_table[is.na(kraken_table)] <- 0

# obtain taxid in kraken results
taxid <- data.frame()
for (x in file_names) {
  kraken_taxid <- read.table(file = x, header = T, quote = "", sep = "\t")
  kraken_taxid <- kraken_taxid[,c(1,2)]
  taxid <- rbind(taxid,kraken_taxid)
}
taxid <- unique(taxid)
setwd("D:/OneDrive/dataanalysis/longevity/data")
taxa<-getTaxonomy(taxid$taxonomy_id,'accessionTaxa.sql')
taxa <- as.data.frame(taxa)
taxa$taxonomy_id <- as.integer(rownames(taxa))
kraken_taxa <- merge(taxa , kraken_table, by = "taxonomy_id")
save(kraken_taxa, file = "kraken_taxa.rda")

#===============================================================================
# 1. Reduce batch effect for each taxa read counts among batchs
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)
rownames(sra_sampleid) <- sra_sampleid$run_id
# At this step, always do factor() for discrete variables, and do droplevels() to drop unused factor levels.

sra_sampleid$bioproject_id <- as.factor(sra_sampleid$bioproject_id)
sra_sampleid$location <- as.factor(sra_sampleid$country)
sra_sampleid$gender <- as.factor(sra_sampleid$gender)
sra_sampleid$category <- as.factor(sra_sampleid$category)   
sra_sampleid$reads <- as.factor(sra_sampleid$reads)


batchid = sra_sampleid[, 'location']
summary(batchid)
#-------------------------------
# 1.1 Assigned counts

load("kraken_taxa.rda")

# Convert reads to TPM
rownames(kraken_taxa) <- kraken_taxa$taxonomy_id
taxa <- kraken_taxa[,c(9:554)]

pdf(paste("F1.1_rawcounts", "pdf", sep="."))
Plot_PCoA(TAX=t(taxa), factor=batchid, main="Before Correction, Bray-Curtis")
dev.off()

# calcualte permanova R2

PERMANOVA_R2(t(taxa),batchid,covar,1) # category

# $tab_count
# standard sqrt.dist=T      add=T
# batch 0.1318531  0.07335235 0.06395214
# key   0.0407508  0.02435416 0.02082305
# 
# $tab_rel
# standard sqrt.dist=T      add=T
# batch 0.11371626  0.06327618 0.11371626
# key   0.02019498  0.01220531 0.02019498

PERMANOVA_R2(t(taxa),batchid,covar,2) # gender
# $tab_count
# standard sqrt.dist=T       add=T
# batch 0.131853115 0.073352346 0.063952140
# key   0.009293103 0.006054958 0.005284969
# 
# $tab_rel
# standard sqrt.dist=T       add=T
# batch 0.113716261 0.063276177 0.113716261
# key   0.005557061 0.003756617 0.005557061


#-----------------------------
# 1.2 Normalization (for each sample ) # is not tpm IS RPM
taxa_tpm <- sweep(taxa,2,colSums(taxa)/(10^6),`/`)
taxa_tpm <-apply(taxa_tpm,1:2,as.integer)
taxa_tpm <- as.data.frame(taxa_tpm)

pdf(paste("F1.1_norm", "pdf", sep="."))
Plot_PCoA(TAX=t(taxa_tpm), factor=batchid, main="Before Correction, Bray-Curtis")
dev.off()

# calcualte permanova R2
# tpm counts

PERMANOVA_R2(t(taxa_tpm),batchid,covar,1)

# $tab_count
# standard sqrt.dist=T      add=T
# batch 0.07407425  0.04253402 0.02848703
# key   0.05018312  0.02947762 0.01925959
# 
# $tab_rel
# standard sqrt.dist=T      add=T
# batch 0.07789612  0.04254993 0.07789612
# key   0.02464458  0.01466518 0.02464458
PERMANOVA_R2(t(taxa_tpm),batchid,covar,2)
# $tab_count
# standard sqrt.dist=T       add=T
# batch 0.07407425 0.042534020 0.028487033
# key   0.01177986 0.007446648 0.005168121
# 
# $tab_rel
# standard sqrt.dist=T       add=T
# batch 0.077896116 0.042549928 0.077896116
# key   0.007144107 0.004467701 0.007144107


#-------------------------------------------
# 2 Correction batch effects
covar = sra_sampleid[, c('category','gender')]
#-------------------------------
# 2.1 after corrected assigned counts
options(warn=-1)
taxa_corrected_japan = ConQuR(tax_tab=t(taxa), batchid=batchid, covariates=covar, batch_ref="Japan")
#taxa_corrected1 = ConQuR(tax_tab=t(taxa_tpm), batchid=batchid, covariates=covar, batch_ref="PRJNA675598")
save(taxa_corrected_japan,file = "taxa_corrected_japan.rda")
load("taxa_corrected_japan.rda")
pdf(paste("F1.1_rawcount_correctedjapan", "pdf", sep="."))
Plot_PCoA(TAX=taxa_corrected_japan, factor=batchid, main="After Correction, Bray-Curtis")
dev.off()


# calcualte permanova R2
# after corrected assigned reads 
PERMANOVA_R2(taxa_corrected_japan,batchid,covar,1)

# $tab_count
# standard sqrt.dist=T       add=T
# batch 0.006736281 0.006993911 0.005965544
# key   0.060248361 0.034455614 0.024844099
# 
# $tab_rel
# standard sqrt.dist=T       add=T
# batch 0.005601828 0.006150905 0.005601828
# key   0.025865081 0.014926701 0.025865081

PERMANOVA_R2(taxa_corrected_japan,batchid,covar,2)
# $tab_count
# standard sqrt.dist=T       add=T
# batch 0.006736281 0.006993911 0.005965544
# key   0.011614509 0.007370217 0.005494864
# 
# $tab_rel
# standard sqrt.dist=T       add=T
# batch 0.005601828 0.006150905 0.005601828
# key   0.006625387 0.004249841 0.006625387
#----------------
## 2.2 Normalization after correction taxa 
TAXA <- taxa_corrected_japan
taxa_tpm <- sweep(TAXA,2,colSums(TAXA)/(10^6),`/`)
taxa_tpm <- sweep(TAXA,2,colSums(taxa)/(10^6),`/`)
taxa_tpm <-apply(taxa_tpm,1:2,as.integer)
taxa_tpm <- as.data.frame(taxa_tpm)

pdf(paste("F1.1_tpmOFcorrectedreads", "pdf", sep="."))
Plot_PCoA(TAX=t(taxa_tpm), factor=batchid, main="Tpm of corrected, Bray-Curtis")
dev.off()

PERMANOVA_R2(taxa_tpm,batchid,covar,1)
# $tab_count
# standard sqrt.dist=T       add=T
# batch 0.005788867  0.00590518 0.005649412
# key   0.037591141  0.02162411 0.020950787
# 
# $tab_rel
# standard sqrt.dist=T       add=T
# batch 0.003987434 0.005278255 0.003987434
# key   0.026490666 0.015267442 0.026490666

PERMANOVA_R2(taxa_tpm,batchid,covar,2)
# $tab_count
# standard sqrt.dist=T       add=T
# batch 0.005788867  0.00590518 0.005649412
# key   0.007498410  0.00490223 0.004720123
# 
# $tab_rel
# standard sqrt.dist=T       add=T
# batch 0.003987434 0.005278255 0.003987434
# key   0.006038507 0.003960754 0.006038507

#----------------------------------
# 2.3 correction of RMP taxa reads 

options(warn=-1)
taxa_corrected_japan_tpm = ConQuR(tax_tab=t(taxa_tpm), batchid=batchid, covariates=covar, batch_ref="Japan")
#taxa_corrected1 = ConQuR(tax_tab=t(taxa_tpm), batchid=batchid, covariates=covar, batch_ref="PRJNA675598")
save(taxa_corrected_japan_tpm,file = "taxa_corrected_japan_tmp.rda")
load("taxa_corrected_japan_tmp.rda")
pdf(paste("F1.1japan_tmp", "pdf", sep="."))
Plot_PCoA(TAX=taxa_corrected_japan_tpm, factor=batchid, main="After Correction, Bray-Curtis")
dev.off()

PERMANOVA_R2(taxa_corrected_japan_tpm,batchid,covar,1)

# $tab_count
# standard sqrt.dist=T       add=T
# batch 0.002527346 0.004929286 0.004556971
# key   0.046609480 0.026944931 0.017336871
# 
# $tab_rel
# standard sqrt.dist=T       add=T
# batch 0.004336678 0.005210447 0.004336678
# key   0.023131264 0.013903585 0.023131264
PERMANOVA_R2(taxa_corrected_japan_tpm,batchid,covar,2)

# $tab_count
# standard sqrt.dist=T       add=T
# batch 0.002527346 0.004929286 0.004556971
# key   0.011416083 0.007166692 0.004884436
# 
# $tab_rel
# standard sqrt.dist=T       add=T
# batch 0.004336678 0.005210447 0.004336678
# key   0.006594819 0.004178288 0.006594819



###=========
# 4.figure reflect the assiged reads related to raw counts
#

# table(sra_sampleid$bioproject_id)
# 4.1 boxplot for distribution of reads before and after assignment
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)



load("kraken_taxa.rda")
# normalized reads based on kingdom. 
# Convert reads to TPM
rownames(kraken_taxa) <- kraken_taxa$taxonomy_id
taxa <- kraken_taxa[,c(9:554)]


sample_rawcounts <- as.data.frame(apply(taxa, 2, sum))
sample_rawcounts$run_id <- rownames(sample_rawcounts)
colnames(sample_rawcounts) <- c("raw_counts","run_id")

lb_rc <- merge(sample_rawcounts,sra_sampleid, by ="run_id")
rownames(lb_rc) <- lb_rc$run_id
lb_rc <- lb_rc[,c("raw_counts","reads","country")]
colnames(lb_rc) <- c("Assigned reads","Raw counts","Country")
lb_rc <- melt(lb_rc)
lb_rc$value <- log2(lb_rc$value+1)
save(lb_rc, file = "lb_rc.rda")
##plot figure
load("lb_rc.rda")
pdf(paste("F1.1readnumbers", "pdf", sep="."),width=12, height=9)
ggplot(lb_rc,aes(variable, value, fill = factor(Country)))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(size=24,color = "black"),axis.text.y=element_text(size=24,color = "black"),axis.title =element_text(size=24,face="bold"))+
  labs(x = "", y = expression("Log"["2"]("reads+1")))+
  scale_fill_manual(values=c("#C2C2C2","#d9544d","#a6cee3","#b2df8a"),label = c("China","Italy_1","Italy_2","Japan"))
dev.off()

#-----------------------------
# 4.2 Relationship read count before and after assignment
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)

load("kraken_taxa.rda")
# normalized reads based on kingdom. 
# Convert reads to TPM
rownames(kraken_taxa) <- kraken_taxa$taxonomy_id
taxa <- kraken_taxa[,c(9:554)]


sample_rawcounts <- as.data.frame(apply(taxa, 2, sum))
sample_rawcounts$run_id <- rownames(sample_rawcounts)
colnames(sample_rawcounts) <- c("raw_counts","run_id")

lb_rc <- merge(sample_rawcounts,sra_sampleid, by ="run_id")
rownames(lb_rc) <- lb_rc$run_id
lb_rc <- lb_rc[,c("run_id","raw_counts","reads","country")]
colnames(lb_rc) <- c("run_id","Assigned reads","Raw counts","Country")


j <- lb_rc[which(lb_rc$Country %in% "Japan"),]
j$`Assigned reads` <-log2(j$`Assigned reads`+1)
j$`Raw counts` <-log2(j$`Raw counts`+1)

wilcox.test(j$`Assigned reads`,j$`Raw counts`)
cor(j$`Assigned reads`,j$`Raw counts`, method = 'pearson')
# 0.9614

i1 <-  lb_rc[which(lb_rc$Country %in% "Italy_1"),]
i1$`Assigned reads` <-log2(i1$`Assigned reads`+1)
i1$`Raw counts` <-log2(i1$`Raw counts`+1)

wilcox.test(i1$`Assigned reads`,i1$`Raw counts`)
cor(i1$`Assigned reads`,i1$`Raw counts`, method = 'pearson')
# 0.69

i2 <-  lb_rc[which(lb_rc$Country %in% "Italy_2"),]
i2$`Assigned reads` <-log2(i2$`Assigned reads`+1)
i2$`Raw counts` <-log2(i2$`Raw counts`+1)

wilcox.test(i2$`Assigned reads`,i2$`Raw counts`)
cor(i2$`Assigned reads`,i2$`Raw counts`, method = 'pearson')
# 0.77
c <-  lb_rc[which(lb_rc$Country %in% "China"),]
c$`Assigned reads` <-log2(c$`Assigned reads`+1)
c$`Raw counts` <-log2(c$`Raw counts`+1)

wilcox.test(i2$`Assigned reads`,i2$`Raw counts`)
cor(c$`Assigned reads`,c$`Raw counts`, method = 'pearson')
# 0.34

# make plots
pdf(paste("F1.4.2readrelation", "pdf", sep="."),width=12, height=9)
par(mfrow=c(2,2))
plot(j[,2:3],main= unique(j$Country))
plot(i1[,2:3],main= unique(i1$Country))
plot(i2[,2:3],main= unique(i2$Country))
plot(c[,2:3],main= unique(c$Country))
dev.off()

###=========
# remove the participants  with illness 


###=========
# fraction of microbio after corrected 
load("kraken_taxa.rda")
load("taxa_corrected_japan_tmp.rda")
taxa_corrected <- t(taxa_corrected_japan_tpm)
taxa_corrected <- as.data.frame(taxa_corrected)
taxa_corrected$taxonomy_id <- rownames(taxa_corrected)
taxa_species <- kraken_taxa[,c("taxonomy_id","superkingdom","phylum","genus","species")]
taxa_corrected <- merge(taxa_species, taxa_corrected, by = "taxonomy_id")


save(taxa_corrected, file = "kraken_taxa_corrected.rda")


# rownames(kraken_taxa) <- kraken_taxa$taxonomy_id

sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)
unique(sra_sampleid$category)
# for group C
c_runid <- sra_sampleid[which(sra_sampleid$category %in% "C"),"run_id"]
c_taxa <- taxa_corrected[,c("superkingdom",c_runid)]
c_taxa$sumreads <- rowSums(c_taxa[,-1])
c_taxa <- c_taxa[c_taxa$sumreads >0,]
table(c_taxa$superkingdom)

# Archaea  Bacteria Eukaryota   Viruses 
# 265      7126        88       639 

# for group W
w_runid <- sra_sampleid[which(sra_sampleid$category %in% "W"),"run_id"]
w_taxa <- taxa_corrected[,c("superkingdom",w_runid)]
w_taxa$sumreads <- rowSums(w_taxa[,-1])
w_taxa <- w_taxa[w_taxa$sumreads >0,]
table(w_taxa$superkingdom)
#
# Archaea  Bacteria Eukaryota   Viruses 
# 201      6790        87       291 



# for group R
r_runid <- sra_sampleid[which(sra_sampleid$category %in% "R"),"run_id"]
r_taxa <- taxa_corrected[,c("superkingdom",r_runid)]
r_taxa$sumreads <- rowSums(r_taxa[,-1])
r_taxa <- r_taxa[r_taxa$sumreads >0,]
table(r_taxa$superkingdom)


#
# Archaea  Bacteria Eukaryota   Viruses 
# 240      7085        91       501


# create a dataframe about number of m
group <- c("C","C","C","C","R","R","R","R","W","W","W","W")
microbiome <- c("Archaea","Bacteria","Eukaryota","Virues","Archaea","Bacteria","Eukaryota","Virues","Archaea","Bacteria","Eukaryota","Virues")
value <- c(265,7126,88,639,240,7085,91,501,201,6790,87,291)
m_number <- data.frame(group,microbiome,value)

# total microbe in each category
# C:8118; R:7917;W:7369
# THE LABLE OF EACH circle is W,R,C from outer circle to inner circle


##plot figure 
pdf(paste("F1.1microbes_corrected", "pdf", sep="."),width=12, height=9)
ggplot(m_number,aes(group,value,fill = microbiome))+
  geom_bar(stat = 'identity',position = 'fill',width = 1,color = "black")+
  coord_polar(theta = "y")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.text = element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=c("#ffffb3","#bebada","#fb8072","#80b1d3"))+
  labs(x = "", y = "case")
dev.off()

###=========
# fraction of microbe before corrected 

load("kraken_taxa.rda")

sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)
unique(sra_sampleid$category)
# for group C
c_runid <- sra_sampleid[which(sra_sampleid$category %in% "C"),"run_id"]
c_taxa <- kraken_taxa[,c("superkingdom",c_runid)]
c_taxa$sumreads <- rowSums(c_taxa[,-1])
c_taxa <- c_taxa[c_taxa$sumreads >0,]
table(c_taxa$superkingdom)

# Archaea  Bacteria Eukaryota   Viruses 
# 315      7370       112       714 

# for group W
w_runid <- sra_sampleid[which(sra_sampleid$category %in% "W"),"run_id"]
w_taxa <- kraken_taxa[,c("superkingdom",w_runid)]
w_taxa$sumreads <- rowSums(w_taxa[,-1])
w_taxa <- w_taxa[w_taxa$sumreads >0,]
table(w_taxa$superkingdom)
#
# Archaea  Bacteria Eukaryota   Viruses 
# 279      7142       112       364 



# for group R
r_runid <- sra_sampleid[which(sra_sampleid$category %in% "R"),"run_id"]
r_taxa <- kraken_taxa[,c("superkingdom",r_runid)]
r_taxa$sumreads <- rowSums(r_taxa[,-1])
r_taxa <- r_taxa[r_taxa$sumreads >0,]
table(r_taxa$superkingdom)

# Archaea  Bacteria Eukaryota   Viruses 
# 300      7300       112       567




# create a dataframe about number of m
group <- c("C","C","C","C","R","R","R","R","W","W","W","W")
microbiome <- c("Archaea","Bacteria","Eukaryota","Virues","Archaea","Bacteria","Eukaryota","Virues","Archaea","Bacteria","Eukaryota","Virues")
value <- c(315,7370,112,714,300,7300,112,567,279,7142,112,364)
m_number <- data.frame(group,microbiome,value)

# total microbe in each category
# C:8511; R:8274;W:7897
# THE LABLE OF EACH circle is W,R,C from outer circle to inner circle


##plot figure 
pdf(paste("F1.1microbes", "pdf", sep="."),width=12, height=9)
ggplot(m_number,aes(group,value,fill = microbiome))+
  geom_bar(stat = 'identity',position = 'fill',width = 1,color = "black")+
  coord_polar(theta = "y")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.text = element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=c("#ffffb3","#bebada","#fb8072","#80b1d3"))+
  labs(x = "", y = "case")
dev.off()


### # total microbes before and after corrected batch effects 
group <- c("C","C","R","R","W","W")
value <- c(8511,8118,8274,7917,7897,7369)
varants <- c("Non-corrected","Corrected","Non-corrected","Corrected","Non-corrected","Corrected")
non_corrected_mircobes <- data.frame(group,value,varants)

pdf(paste("F1.1totalmicrobes", "pdf", sep="."),width=12, height=9)
ggplot(non_corrected_mircobes,aes(varants,value,fill = group))+
  geom_bar(stat = "identity",width = 0.7)+labs(x= "", y="")+ coord_flip()+theme_classic()+
  scale_fill_manual(values=c("#6a3d9a","#fdb462","#33a02c"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+guides(fill="none")
dev.off()

#===============================================================================
#
#2. Microbiome communities  diversity analysis_FOR LOCATIONS 
#
#===============================================================================

########################
#
# alpha diversity and richness analysis 
#
########################
load("kraken_taxa_corrected.rda")
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)

# ----------------------------------For bacteria -------------------------------
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]

#--------------------------
# alpha diversity based on shannon index

diversity_shannon_b <- diversity(t(b_taxa_a), index = "shannon", MARGIN = 1)
# diversity_shannon_b <- diversity(t(b_taxa_a), index = "simpson", MARGIN = 1)

diversity_shannon_b <- as.data.frame(diversity_shannon_b)
diversity_shannon_b$run_id <- rownames(diversity_shannon_b) 
diversity_shannon_b <- merge(diversity_shannon_b,sra_sampleid, by="run_id")
bb <- diversity_shannon_b
kruskal.test(diversity_shannon_b~country, data = diversity_shannon_b) # alpha diversity based on location 
# 0.07314

compaired <- list(c("China","Italy_1"),c("China","Italy_2"),c("China","Japan"),c("Italy_1","Italy_2"),c("Italy_1","Japan"),c("Italy_2","Japan"))
pdf(paste("F2.2bacteria", "pdf", sep="."),width=12, height=9)
ggboxplot(diversity_shannon_b,x="country",y="diversity_shannon_b",fill = "country",add="jitter",size=0.5)+
  #stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()


# shannon diversity index density analysis 
pdf(paste("F3.1bacteria", "pdf", sep="."),width=12, height=9)
ggdensity(diversity_shannon_b,x='diversity_shannon_b',color ="country",fill = "country", add="mean", rug = TRUE)+
  
  #  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  #  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  geom_density()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  # scale_fill_manual(label = c("China","Italy_1","Italy_2","Japan"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")#+guides(fill="none")+
dev.off()




# ----------------------------------For Archaea -------------------------------
load("kraken_taxa_corrected.rda")
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)

b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Archaea"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]
 diversity_shannon_a <- diversity(t(b_taxa_a), index = "shannon", MARGIN = 1)
#diversity_shannon_a <- diversity(t(b_taxa_a), index = "simpson", MARGIN = 1)

diversity_shannon_a <- as.data.frame(diversity_shannon_a)
diversity_shannon_a$run_id <- rownames(diversity_shannon_a) 
diversity_shannon_a <- merge(diversity_shannon_a,sra_sampleid, by="run_id")

# 
# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_a by country


kruskal.test(diversity_shannon_a~country, data = diversity_shannon_a)

# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_a by country
# Kruskal-Wallis chi-squared = 8.134, df = 3, p-value = 0.04332

compaired <- list(c("China","Italy_1"),c("China","Italy_2"),c("China","Japan"),c("Italy_1","Italy_2"),c("Italy_1","Japan"),c("Italy_2","Japan"))

pdf(paste("F2.2Archaea", "pdf", sep="."),width=12, height=9)
ggboxplot(diversity_shannon_a,x="country",y="diversity_shannon_a",fill = "country",add="jitter",size=0.5)+
#  stat_compare_means(comparisons = compaired, method = "wilcox.test")+ 
  stat_compare_means( method = "kruskal.test")+
  theme_light()+
  scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()

#-----------------------------
# shannon diversity index density analysis 
pdf(paste("F3.1archaea", "pdf", sep="."),width=12, height=9)
ggdensity(diversity_shannon_a,x='diversity_shannon_a',color ="country",fill = "country", add="mean", rug = TRUE)+
  
  #  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  #  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  geom_density()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  # scale_fill_manual(label = c("China","Italy_1","Italy_2","Japan"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")#+guides(fill="none")+
dev.off()

#------------------------------
# find a cutoff for remove sample that leading difference of shannon diversity among cohorts
summary(diversity_shannon_a$diversity_shannon_a)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.2952  1.2558  1.2958  2.1391  3.6843

for (x in sort(diversity_shannon_a$diversity_shannon_a)) {
  #x <- x+0.000001
  aa <- diversity_shannon_a[which(diversity_shannon_a$diversity_shannon_a > x),]
  if ( kruskal.test(diversity_shannon_a~country, data = aa)$p.value >0.05){
    print(x)
    break
  }
}
kruskal.test(diversity_shannon_a~country, data = aa)
# p-value = 0.05049
# 0.6892943

# ----------------------------------For Eukaryota -------------------------------
load("kraken_taxa_corrected.rda")
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)

b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Eukaryota"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]
b_taxa_a <- b_taxa_a[-which(rownames(b_taxa_a) %in% "9606"),]
diversity_shannon_e <- diversity(t(b_taxa_a), index = "shannon", MARGIN = 1)
# diversity_shannon_e <- diversity(t(b_taxa_a), index = "simpson", MARGIN = 1)

diversity_shannon_e <- as.data.frame(diversity_shannon_e)
diversity_shannon_e$run_id <- rownames(diversity_shannon_e) 
diversity_shannon_e <- merge(diversity_shannon_e,sra_sampleid, by="run_id")

# 
# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_e by country


kruskal.test(diversity_shannon_e~country, data = diversity_shannon_e)

# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_e by country
# Kruskal-Wallis chi-squared = 11.372, df = 3, p-value = 0.009877

compaired <- list(c("China","Italy_1"),c("China","Italy_2"),c("China","Japan"),c("Italy_1","Italy_2"),c("Italy_1","Japan"),c("Italy_2","Japan"))

pdf(paste("F2.2Eukaryota", "pdf", sep="."),width=12, height=9)
ggboxplot(diversity_shannon_e,x="country",y="diversity_shannon_e",fill = "country",add="jitter",size=0.5)+
#  stat_compare_means(comparisons = compaired, method = "wilcox.test")+ 
  stat_compare_means( method = "kruskal.test")+
  theme_light()+
  scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()

#-----------------------------
# shannon diversity index density analysis 
pdf(paste("F3.1fungi", "pdf", sep="."),width=12, height=9)
ggdensity(diversity_shannon_e,x='diversity_shannon_e',color ="country",fill = "country", add="mean", rug = TRUE)+
  
  #  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  #  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  geom_density()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  # scale_fill_manual(label = c("China","Italy_1","Italy_2","Japan"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")#+guides(fill="none")+
dev.off()

# 
# aa <- top_n(diversity_shannon_e,-10, diversity_shannon_e$diversity_shannon_e)
summary(diversity_shannon_e$diversity_shannon_e)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1685  0.8404  0.9585  1.5624  3.2999 
for (x in sort(diversity_shannon_e$diversity_shannon_e)) {
  #x <- x+0.000001
  ff <- diversity_shannon_e[which(diversity_shannon_e$diversity_shannon_e > x),]
  if ( kruskal.test(diversity_shannon_e~country, data = ff)$p.value >0.05){
    print(x)
    break
  }
}
kruskal.test(diversity_shannon_e~country, data = ff)
# p-value = 0.8941
# 0

# ----------------------------------For Viruses -------------------------------
load("kraken_taxa_corrected.rda")
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)

b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Viruses"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]
 diversity_shannon_v <- diversity(t(b_taxa_a), index = "shannon", MARGIN = 1)
# diversity_shannon_v <- diversity(t(b_taxa_a), index = "simpson", MARGIN = 1)

diversity_shannon_v <- as.data.frame(diversity_shannon_v)
diversity_shannon_v$run_id <- rownames(diversity_shannon_v) 
diversity_shannon_v <- merge(diversity_shannon_v,sra_sampleid, by="run_id")

# 
# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_v by country


kruskal.test(diversity_shannon_v~country, data = diversity_shannon_v)

# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_v by country
# Kruskal-Wallis chi-squared = 53.511, df = 3, p-value = 1.427e-11


compaired <- list(c("China","Italy_1"),c("China","Italy_2"),c("China","Japan"),c("Italy_1","Italy_2"),c("Italy_1","Japan"),c("Italy_2","Japan"))

pdf(paste("F2.2viruses", "pdf", sep="."),width=12, height=9)
ggboxplot(diversity_shannon_v,x="country",y="diversity_shannon_v",fill = "country",add="jitter",size=0.5)+
#  stat_compare_means(comparisons = compaired, method = "wilcox.test")+ 
  stat_compare_means( method = "kruskal.test")+
  theme_light()+
  scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()


# shannon diversity index density analysis 
pdf(paste("F3.1virus", "pdf", sep="."),width=12, height=9)
ggdensity(diversity_shannon_v,x='diversity_shannon_v',color ="country",fill = "country", add="mean", rug = TRUE)+
  
  #  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  #  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  geom_density()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  # scale_fill_manual(label = c("China","Italy_1","Italy_2","Japan"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")#+guides(fill="none")+
dev.off()
#-------------------------------
summary(diversity_shannon_v$diversity_shannon_v)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.7181  1.0652  1.1508  1.5270  3.2464

 
for (x in sort(diversity_shannon_v$diversity_shannon_v)) {
  #x <- x+0.000001
  vv <- diversity_shannon_v[which(diversity_shannon_v$diversity_shannon_v > x),]
  if ( kruskal.test(diversity_shannon_v~country, data = vv)$p.value >0.05){
    print(x)
    break
  }
}
kruskal.test(diversity_shannon_v~country, data = vv)
# p-value = 0.05121
# 1.112332


### combine aa ff vv and diversity shannon v to find common samples 
# 1.1 save runid of microbe groups based on shannon diversity
b_samples <- sra_sampleid[which(sra_sampleid$run_id %in% bb$run_id ),]
f_samples <- sra_sampleid[which(sra_sampleid$run_id %in% ff$run_id ),]
a_samples <- sra_sampleid[which(sra_sampleid$run_id %in% aa$run_id ),]
v_samples <- sra_sampleid[which(sra_sampleid$run_id %in% vv$run_id ),]
save(b_samples,file = "b_samples.rda") # 546
save(f_samples,file = "f_samples.rda") # 424
save(a_samples,file = "a_samples.rda") # 356
save(v_samples,file = "v_samples.rda") # 261
load("b_samples.rda")
load("f_samples.rda")
load("a_samples.rda")
load("v_samples.rda")
# use to make venn figure

X <- list(Bacteria =b_samples$run_id,
          Fungi =f_samples$run_id,
          Archaea = a_samples$run_id,
          Viruses = v_samples$run_id)


pdf(paste("F4.1", "pdf", sep="."),width=12, height=9)
ggVennDiagram(X,label_alpha = 0)+ggplot2::scale_fill_gradient(low="lightblue",high = 'yellow')
dev.off()


a <- as.data.frame(aa$run_id)
colnames(a) <- "afvb"
f <- as.data.frame(ff$run_id)
colnames(f) <- "afvb"
v <- as.data.frame(vv$run_id)
colnames(v) <- "afvb"
b <- as.data.frame(diversity_shannon_b$run_id)
colnames(b) <- "afvb"
afvb <- rbind(a,f,v,b)
afvb <- as.data.frame(table(afvb))
afvb <- afvb[which(afvb$Freq >3),]

comsamples <- sra_sampleid[which(sra_sampleid$run_id %in% afvb$afvb ),] # 161
save(comsamples, file = "comsamples.rda")
load('comsamples.rda')
comsamples <- comsamples %>% mutate(category = case_when(age < 60  ~ "Y",
                                                             age < 90 & age > 65 ~ "E",
                                                             age > 95  ~ "C",
                                                             TRUE ~"Not listed"))
comsample_id <- comsamples[-which(comsamples$category %in% "Not listed"),]
table(comsample_id$category)
save(comsample_id,file = 'comsample_id.rda')


# plot using ggplot
g_age_gender <- comsample_id[,c("gender","country","category")]


# df2 <- g_age_gender %>% group_by(gender,country,category) %>% 
#   summarise(total_count=n(),.groups = 'drop') %>%
#   as.data.frame()

df2 <- g_age_gender %>% group_by(gender,category) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()


df2_m<- df2[which(df2$gender %in% "M"),]
df2_f <- df2[which(df2$gender %in% "F"),]
df2_f$total_count <- df2_f$total_count *-1




pdf(paste("4.2comage_gender", "pdf", sep="."),width=12, height=9)
ggplot()+
  geom_bar(data =df2_m,aes(category,total_count,fill=gender),stat = "identity",position="stack",color="black",width = 0.4,size=0.25)+
  
  # Famale population data graph
  geom_bar(data =df2_f,aes(category,y= total_count,fill=gender),stat = "identity",position="stack",color="black",width = 0.4,size=0.25)+
  
  #geom_bar(data =df2_f,aes(category,y= total_count *(-1),fill=country),stat = "identity",position="stack",color="black",width = 0.7,size=0.25)+
  #coord_flip()+
  scale_y_continuous(limits = c(-50,30),
                     breaks = seq(-50,30,25),
                     labels = abs(seq(-50,30,25)))+
  labs(x = "", y = "")+
  scale_fill_manual(values=c("#FFA07A","#87CEEB"),label = c("M","F"))+
  # display themes
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 24, face = "bold"),  # fonts/sizes
    axis.title = element_text(size = 24, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5), 
    plot.caption = element_text(hjust=0, size=11, face = "italic")
  )

dev.off()







####################################################################################
#
# PCOA analysis for population and age group
#
###################################################################################


#------------------------------------------------#
#                                                #
#           Beta for bacteria                    # 
#                                                #
#------------------------------------------------#


set.seed(2023)

load("kraken_taxa_corrected.rda")
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
count_data_zeroReplaced <- b_taxa[,-1:-5]
# Define pseudocount
pseudocount <- 0.5

# Impute zero
count_data_zeroReplaced[count_data_zeroReplaced == 0] <- pseudocount

# CLR transformation
count_data_CLR <- t(apply(count_data_zeroReplaced, 1, function(x) log(x) - mean(log(x))))

# Calculate Euclidean distances
Euclidean_CLR <- vegdist(t(count_data_CLR), method = "euclidean")


# 
# # dist_b <- vegdist(t(b_taxa_a),method = "bray")
# dist_b <- dist(t(b_taxa_a), method = "euclidean")
# # pcoa_b <-wcmdscale(bc_b) 

# print(pcoa_b)

pcoa <- cmdscale (Euclidean_CLR, eig = TRUE,add = TRUE) # add make eig convert to positive value 
positions <- pcoa$points
colnames(positions) <- c("pcoa1","pcoa2")
percent_explained <- 100*pcoa$eig / sum(pcoa$eig)
prett_pe <- format(round(percent_explained[1:2],digits=1),nsmall = 1,trim=TRUE)
labs <- c(glue("PCo 1 ({prett_pe[1]}%)"),
          glue("PCo 2 ({prett_pe[2]}%)"))
positions_samples <- positions %>% as_tibble(rownames="run_id")
metadata_gray <- inner_join(positions_samples,sra_sampleid)

# variance is country

centroid_country <- metadata_gray %>% group_by(country) %>%
  summarize(pcoa1 =mean(pcoa1),pcoa2 =mean(pcoa2))

## p_PERMANOVA
value_adnois <-adonis(Euclidean_CLR~metadata_gray$country,permutations = 999,method = "euclidean")
value_adnois$aov.tab  
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# metadata_gray$country   3      7259  2419.8 0.37844 0.00209      1
# Residuals             542   3465633  6394.2         0.99791       
# Total                 545   3472892                 1.00000

# BETADISPER
anova(betadisper(Euclidean_CLR,metadata_gray$country)) # or use permutest(betadisper(Euclidean_CLR,metadata_gray$country))
# Analysis of Variance Table
# 
# Response: Distances
# Df Sum Sq Mean Sq F value Pr(>F)
# Groups      3    805  268.32  0.8102 0.4885


pdf(paste("F2.4_bacteria_country", "pdf", sep="."),width=14, height=12)
 ggplot(metadata_gray,aes(pcoa1,pcoa2,color = country))+
  geom_point()+
  scale_color_manual( values = c("#ffcb4dff",  "#282554ff", "#de2524ff",
                                 "#6b9064ff"))+

  stat_ellipse(show.legend = F)+
  geom_point(data = centroid_country, size=5,shape=21,color="black",aes(fill=country))+
  scale_fill_manual( values = c("#ffcb4dff",  "#282554ff", "#de2524ff",
                                "#6b9064ff"))+
  theme_light()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = labs[1], y = labs[2])+guides(fill="none")
dev.off()



#------------------------------------------------#
#                                                #
#           Beta for Archaea                     # 
#                                                #
#------------------------------------------------#


load("kraken_taxa_corrected.rda")
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Archaea"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
count_data_zeroReplaced <- b_taxa[,-1:-5]
count_data_zeroReplaced <- count_data_zeroReplaced[,which(colSums(count_data_zeroReplaced)>0)]
# Define pseudocount
pseudocount <- 0.5

# Impute zero
count_data_zeroReplaced[count_data_zeroReplaced == 0] <- pseudocount

# CLR transformation
count_data_CLR <- t(apply(count_data_zeroReplaced, 1, function(x) log(x) - mean(log(x))))

# Calculate Euclidean distances
Euclidean_CLR <- vegdist(t(count_data_CLR), method = "euclidean")

# plot pcoa
pcoa <- cmdscale (Euclidean_CLR, eig = TRUE,add = TRUE) # add make eig convert to positive value 
positions <- pcoa$points
colnames(positions) <- c("pcoa1","pcoa2")
percent_explained <- 100*pcoa$eig / sum(pcoa$eig)
prett_pe <- format(round(percent_explained[1:2],digits=1),nsmall = 1,trim=TRUE)
labs <- c(glue("PCo 1 ({prett_pe[1]}%)"),
          glue("PCo 2 ({prett_pe[2]}%)"))
positions_samples <- positions %>% as_tibble(rownames="run_id")
metadata_gray <- inner_join(positions_samples,sra_sampleid)

# variance is country

centroid_country <- metadata_gray %>% group_by(country) %>%
  summarize(pcoa1 =mean(pcoa1),pcoa2 =mean(pcoa2))

## p_PERMANOVA
set.seed(2023)
value_adnois <-adonis(Euclidean_CLR~metadata_gray$country,permutations = 999,method = "euclidean")
value_adnois$aov.tab  
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# metadata_gray$country   3     170.7  56.905  1.1585 0.00672  0.215
# Residuals             514   25247.9  49.120         0.99328       
# Total                 517   25418.6                 1.00000 

# BETADISPER
anova(betadisper(Euclidean_CLR,metadata_gray$country)) # or use permutest(betadisper(Euclidean_CLR,metadata_gray$country))
# Analysis of Variance Table
# 
# Response: Distances
# Df Sum Sq Mean Sq F value Pr(>F)
# Groups      3   22.1  7.3554  0.8109 0.4882
# Residuals 514 4662.5  9.0710  

pdf(paste("F2.4_Archaea_country", "pdf", sep="."),width=14, height=12)
ggplot(metadata_gray,aes(pcoa1,pcoa2,color = country))+
  geom_point()+
  scale_color_manual( values = c("#ffcb4dff",  "#282554ff", "#de2524ff",
                                 "#6b9064ff"))+
  
  stat_ellipse(show.legend = F)+
  geom_point(data = centroid_country, size=5,shape=21,color="black",aes(fill=country))+
  scale_fill_manual( values = c("#ffcb4dff",  "#282554ff", "#de2524ff",
                                "#6b9064ff"))+
  theme_light()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = labs[1], y = labs[2])+guides(fill="none")
dev.off()


#------------------------------------------------#
#                                                #
#           Beta for Eukaryota                   # 
#                                                #
#------------------------------------------------#


load("kraken_taxa_corrected.rda")
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Eukaryota"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
count_data_zeroReplaced <- b_taxa[,-1:-5]
count_data_zeroReplaced <- count_data_zeroReplaced[,which(colSums(count_data_zeroReplaced)>0)]
# Define pseudocount
pseudocount <- 0.5

# Impute zero
count_data_zeroReplaced[count_data_zeroReplaced == 0] <- pseudocount

# CLR transformation
count_data_CLR <- t(apply(count_data_zeroReplaced, 1, function(x) log(x) - mean(log(x))))

# Calculate Euclidean distances
Euclidean_CLR <- vegdist(t(count_data_CLR), method = "euclidean")

# plot pcoa
pcoa <- cmdscale (Euclidean_CLR, eig = TRUE,add = TRUE) # add make eig convert to positive value 
positions <- pcoa$points
colnames(positions) <- c("pcoa1","pcoa2")
percent_explained <- 100*pcoa$eig / sum(pcoa$eig)
prett_pe <- format(round(percent_explained[1:2],digits=1),nsmall = 1,trim=TRUE)
labs <- c(glue("PCo 1 ({prett_pe[1]}%)"),
          glue("PCo 2 ({prett_pe[2]}%)"))
positions_samples <- positions %>% as_tibble(rownames="run_id")
metadata_gray <- inner_join(positions_samples,sra_sampleid)

# variance is country

centroid_country <- metadata_gray %>% group_by(country) %>%
  summarize(pcoa1 =mean(pcoa1),pcoa2 =mean(pcoa2))

## p_PERMANOVA
set.seed(2023)
value_adnois <-adonis(Euclidean_CLR~metadata_gray$country,permutations = 999,method = "euclidean")
value_adnois$aov.tab  
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# metadata_gray$country   3     114.9  38.295  2.2373 0.01223  0.001 ***
#   Residuals             542    9277.3  17.117         0.98777           
# Total                 545    9392.2                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 

# BETADISPER
set.seed(2023)
anova(betadisper(Euclidean_CLR,metadata_gray$country)) # or use permutest(betadisper(Euclidean_CLR,metadata_gray$country))

# Analysis of Variance Table
# 
# Response: Distances
# Df  Sum Sq Mean Sq F value Pr(>F)
# Groups      3   11.68  3.8935   0.894  0.444
# Residuals 542 2360.45  4.3551  


pdf(paste("F2.4_Eukaryota_country", "pdf", sep="."),width=14, height=12)
ggplot(metadata_gray,aes(pcoa1,pcoa2,color = country))+
  geom_point()+
  scale_color_manual( values = c("#ffcb4dff",  "#282554ff", "#de2524ff",
                                 "#6b9064ff"))+
  
  stat_ellipse(show.legend = F)+
  geom_point(data = centroid_country, size=5,shape=21,color="black",aes(fill=country))+
  scale_fill_manual( values = c("#ffcb4dff",  "#282554ff", "#de2524ff",
                                "#6b9064ff"))+
  theme_light()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = labs[1], y = labs[2])+guides(fill="none")
dev.off()

#------------------------------------------------#
#                                                #
#           Beta for Viruses                     # 
#                                                #
#------------------------------------------------#


load("kraken_taxa_corrected.rda")
sra_sampleid <- fread("sra_sampleid.csv")
sra_sampleid <- as.data.frame(sra_sampleid)
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Viruses"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
count_data_zeroReplaced <- b_taxa[,-1:-5]
count_data_zeroReplaced <- count_data_zeroReplaced[,which(colSums(count_data_zeroReplaced)>0)]
# Define pseudocount
pseudocount <- 0.5

# Impute zero
count_data_zeroReplaced[count_data_zeroReplaced == 0] <- pseudocount

# CLR transformation
count_data_CLR <- t(apply(count_data_zeroReplaced, 1, function(x) log(x) - mean(log(x))))

# Calculate Euclidean distances
Euclidean_CLR <- vegdist(t(count_data_CLR), method = "euclidean")

# plot pcoa
pcoa <- cmdscale (Euclidean_CLR, eig = TRUE,add = TRUE) # add make eig convert to positive value 
positions <- pcoa$points
colnames(positions) <- c("pcoa1","pcoa2")
percent_explained <- 100*pcoa$eig / sum(pcoa$eig)
prett_pe <- format(round(percent_explained[1:2],digits=1),nsmall = 1,trim=TRUE)
labs <- c(glue("PCo 1 ({prett_pe[1]}%)"),
          glue("PCo 2 ({prett_pe[2]}%)"))
positions_samples <- positions %>% as_tibble(rownames="run_id")
metadata_gray <- inner_join(positions_samples,sra_sampleid)

# variance is country

centroid_country <- metadata_gray %>% group_by(country) %>%
  summarize(pcoa1 =mean(pcoa1),pcoa2 =mean(pcoa2))

## p_PERMANOVA
set.seed(2023)
value_adnois <-adonis(Euclidean_CLR~metadata_gray$country,permutations = 999,method = "euclidean")
value_adnois$aov.tab  
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# metadata_gray$country   3      5579 1859.80  10.187 0.05394  0.001 ***
#   Residuals             536     97854  182.56         0.94606           
# Total                 539    103434                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#----
# BETADISPER
set.seed(2023)
anova(betadisper(Euclidean_CLR,metadata_gray$country)) # or use permutest(betadisper(Euclidean_CLR,metadata_gray$country))
# Analysis of Variance Table
# 
# Response: Distances
# Df  Sum Sq Mean Sq F value Pr(>F)
# Groups      3   187.2  62.417  2.0874 0.1009
# Residuals 536 16026.9  29.901  
pdf(paste("F2.4_Viruses_country", "pdf", sep="."),width=14, height=12)
ggplot(metadata_gray,aes(pcoa1,pcoa2,color = country))+
  geom_point()+
  scale_color_manual( values = c("#ffcb4dff",  "#282554ff", "#de2524ff",
                                 "#6b9064ff"))+
  
  stat_ellipse(show.legend = F)+
  geom_point(data = centroid_country, size=5,shape=21,color="black",aes(fill=country))+
  scale_fill_manual( values = c("#ffcb4dff",  "#282554ff", "#de2524ff",
                                "#6b9064ff"))+
  theme_light()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = labs[1], y = labs[2])+guides(fill="none")
dev.off()




# ------------------------------------------------------------
# Microbiome community alpha diversity analysis in three age groups
# ------------------------------------------------------------

########################
#
# alpha diversity and richness analysis 
#
###########################

load("kraken_taxa_corrected.rda")
load('comsample_id.rda')

# ----------------------------------For bacteria -------------------------------
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]

#--------------------------
# alpha diversity 

diversity_shannon_b <- diversity(t(b_taxa_a), index = "shannon", MARGIN = 1)
# diversity_shannon_b <- diversity(t(b_taxa_a), index = "simpson", MARGIN = 1)

diversity_shannon_b <- as.data.frame(diversity_shannon_b)
diversity_shannon_b$run_id <- rownames(diversity_shannon_b) 
diversity_shannon_b <- merge(diversity_shannon_b,comsample_id, by="run_id")

kruskal.test(diversity_shannon_b~category, data = diversity_shannon_b) # alpha diversity based on age category

# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_b by category
# Kruskal-Wallis chi-squared = 1.6871, df = 2, p-value = 0.4302

compaired <- list(c("Y","E"),c("Y","C"),c("E","C"))
pdf(paste("F4.3bacteria_agegroup", "pdf", sep="."),width=12, height=9)
ggboxplot(diversity_shannon_b,x="category",y="diversity_shannon_b",fill = "category",add="jitter",size=0.5)+
  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  stat_compare_means(method = "kruskal.test")+
  theme_light()+
#   scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  scale_fill_manual(values=c("#FFFFFF","#FFFFFF","#FFFFFF"),label = c("Y","W","C"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()


#-----------------------------
# Richness analysis 

richness_b <- estimateR(t(b_taxa_a))
richness_b <- as.data.frame(t(richness_b))
richness_b$run_id <- rownames(richness_b)
richness_b <- merge(richness_b,comsample_id, by="run_id")
kruskal.test(S.chao1~category, data = richness_b)


compaired <- list(c("Y","E"),c("Y","C"),c("E","C"))
pdf(paste("F4.3bacteria_agegroup_chao1", "pdf", sep="."),width=12, height=9)
ggboxplot(richness_b,x="category",y="S.chao1",fill = "category",add="jitter",size=0.5)+
  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  scale_fill_manual(values=c("#F0F8FF","#F0F8FF","#F0F8FF"),label = c("Y","W","C"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()


# ----------------------------------For Archaea -------------------------------

b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Archaea"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]


#--------------------------
# alpha diversity 

diversity_shannon_b <- diversity(t(b_taxa_a), index = "shannon", MARGIN = 1)
# diversity_shannon_b <- diversity(t(b_taxa_a), index = "simpson", MARGIN = 1)

diversity_shannon_b <- as.data.frame(diversity_shannon_b)
diversity_shannon_b$run_id <- rownames(diversity_shannon_b) 
diversity_shannon_b <- merge(diversity_shannon_b,comsample_id, by="run_id")

kruskal.test(diversity_shannon_b~category, data = diversity_shannon_b) # alpha diversity based on age category

# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_b by category
# Kruskal-Wallis chi-squared = 11.191, df = 2, p-value = 0.003715

compaired <- list(c("Y","E"),c("Y","C"),c("E","C"))
pdf(paste("F4.3archaea_agegroup", "pdf", sep="."),width=12, height=9)
ggboxplot(diversity_shannon_b,x="category",y="diversity_shannon_b",fill = "category",add="jitter",size=0.5)+
  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  #   scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  scale_fill_manual(values=c("#FFFFFF","#FFFFFF","#FFFFFF"),label = c("Y","W","C"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()


#-----------------------------
# Richness analysis 

richness_b <- estimateR(t(b_taxa_a))
richness_b <- as.data.frame(t(richness_b))
richness_b$run_id <- rownames(richness_b)
richness_b <- merge(richness_b,comsample_id, by="run_id")
kruskal.test(S.chao1~category, data = richness_b)


compaired <- list(c("Y","E"),c("Y","C"),c("E","C"))
pdf(paste("F4.3archaea_agegroup_chao1", "pdf", sep="."),width=12, height=9)
ggboxplot(richness_b,x="category",y="S.chao1",fill = "category",add="jitter",size=0.5)+
  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  scale_fill_manual(values=c("#F0F8FF","#F0F8FF","#F0F8FF"),label = c("Y","W","C"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()


# ----------------------------------For Eukaryota -------------------------------
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Eukaryota"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]

#--------------------------
# alpha diversity 

diversity_shannon_b <- diversity(t(b_taxa_a), index = "shannon", MARGIN = 1)
# diversity_shannon_b <- diversity(t(b_taxa_a), index = "simpson", MARGIN = 1)

diversity_shannon_b <- as.data.frame(diversity_shannon_b)
diversity_shannon_b$run_id <- rownames(diversity_shannon_b) 
diversity_shannon_b <- merge(diversity_shannon_b,comsample_id, by="run_id")

kruskal.test(diversity_shannon_b~category, data = diversity_shannon_b) # alpha diversity based on age category

# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_b by category
# Kruskal-Wallis chi-squared = 3.9209, df = 2, p-value = 0.1408

#compaired <- list(c("Y","W"),c("Y","C"),c("W","C"))
pdf(paste("F4.3fungi_agegroup", "pdf", sep="."),width=12, height=9)
ggboxplot(diversity_shannon_b,x="category",y="diversity_shannon_b",fill = "category",add="jitter",size=0.5)+
  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  #   scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  scale_fill_manual(values=c("#FFFFFF","#FFFFFF","#FFFFFF"),label = c("Y","W","C"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()


#-----------------------------
# Richness analysis 

richness_b <- estimateR(t(b_taxa_a))
richness_b <- as.data.frame(t(richness_b))
richness_b$run_id <- rownames(richness_b)
richness_b <- merge(richness_b,comsample_id, by="run_id")
kruskal.test(S.chao1~category, data = richness_b)


compaired <- list(c("Y","E"),c("Y","C"),c("E","C"))
pdf(paste("F4.3bfungi_agegroup_chao1", "pdf", sep="."),width=12, height=9)
ggboxplot(richness_b,x="category",y="S.chao1",fill = "category",add="jitter",size=0.5)+
  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  scale_fill_manual(values=c("#F0F8FF","#F0F8FF","#F0F8FF"),label = c("Y","W","C"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()



# ----------------------------------For Viruses -------------------------------
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Viruses"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]

#--------------------------
# alpha diversity 

diversity_shannon_b <- diversity(t(b_taxa_a), index = "shannon", MARGIN = 1)
# diversity_shannon_b <- diversity(t(b_taxa_a), index = "simpson", MARGIN = 1)

diversity_shannon_b <- as.data.frame(diversity_shannon_b)
diversity_shannon_b$run_id <- rownames(diversity_shannon_b) 
diversity_shannon_b <- merge(diversity_shannon_b,comsample_id, by="run_id")

kruskal.test(diversity_shannon_b~category, data = diversity_shannon_b) # alpha diversity based on age category

# Kruskal-Wallis rank sum test
# 
# data:  diversity_shannon_b by category
# Kruskal-Wallis chi-squared = 4.7685, df = 2, p-value = 0.09216

compaired <- list(c("Y","E"),c("Y","C"),c("E","C"))
pdf(paste("F4.3viruses_agegroup", "pdf", sep="."),width=12, height=9)
ggboxplot(diversity_shannon_b,x="category",y="diversity_shannon_b",fill = "category",add="jitter",size=0.5)+
  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  #   scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  scale_fill_manual(values=c("#FFFFFF","#FFFFFF","#FFFFFF"),label = c("Y","W","C"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()


#-----------------------------
# Richness analysis 

richness_b <- estimateR(t(b_taxa_a))
richness_b <- as.data.frame(t(richness_b))
richness_b$run_id <- rownames(richness_b)
richness_b <- merge(richness_b,comsample_id, by="run_id")
kruskal.test(S.chao1~category, data = richness_b)


compaired <- list(c("Y","E"),c("Y","C"),c("E","C"))
pdf(paste("F4.3viruses_agegroup_chao1", "pdf", sep="."),width=12, height=9)
ggboxplot(richness_b,x="category",y="S.chao1",fill = "category",add="jitter",size=0.5)+
  stat_compare_means(comparisons = compaired, method = "wilcox.test")+
  stat_compare_means(method = "kruskal.test")+
  theme_light()+
  scale_fill_manual(values=c("#F0F8FF","#F0F8FF","#F0F8FF"),label = c("Y","W","C"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")+guides(fill="none")
dev.off()



########################
#
# beta diveristy for age groups
#
########################

#------------------------------------------------#
#                                                #
#           Beta for bacteria                    # 
#                                                #
#------------------------------------------------#


set.seed(2023)

load("kraken_taxa_corrected.rda")
load('comsample_id.rda')

# ----------------------------------For bacteria -------------------------------
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
count_data_zeroReplaced <- b_taxa[,-1:-5]
# Define pseudocount
pseudocount <- 0.5

# Impute zero
count_data_zeroReplaced[count_data_zeroReplaced == 0] <- pseudocount

# CLR transformation
count_data_CLR <- t(apply(count_data_zeroReplaced, 1, function(x) log(x) - mean(log(x))))
count_data_CLR <- count_data_CLR[,comsample_id$run_id]
# Calculate Euclidean distances
Euclidean_CLR <- vegdist(t(count_data_CLR), method = "euclidean")


# 
# # dist_b <- vegdist(t(b_taxa_a),method = "bray")
# dist_b <- dist(t(b_taxa_a), method = "euclidean")
# # pcoa_b <-wcmdscale(bc_b) 

# print(pcoa_b)

pcoa <- cmdscale (Euclidean_CLR, eig = TRUE,add = TRUE) # add make eig convert to positive value 
positions <- pcoa$points
colnames(positions) <- c("pcoa1","pcoa2")
percent_explained <- 100*pcoa$eig / sum(pcoa$eig)
prett_pe <- format(round(percent_explained[1:2],digits=1),nsmall = 1,trim=TRUE)
labs <- c(glue("PCo 1 ({prett_pe[1]}%)"),
          glue("PCo 2 ({prett_pe[2]}%)"))
positions_samples <- positions %>% as_tibble(rownames="run_id")
metadata_gray <- inner_join(positions_samples,comsample_id)

# variance is country

centroid_country <- metadata_gray %>% group_by(category) %>%
  summarize(pcoa1 =mean(pcoa1),pcoa2 =mean(pcoa2))

## p_PERMANOVA
value_adnois <-adonis(Euclidean_CLR~metadata_gray$category,permutations = 999,method = "euclidean")
value_adnois$aov.tab  
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# metadata_gray$category   2     27125 13562.6  2.2792 0.02949  0.004 **
#   Residuals              150    892579  5950.5         0.97051          
# Total                  152    919704                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# BETADISPER
anova(betadisper(Euclidean_CLR,metadata_gray$category)) # or use permutest(betadisper(Euclidean_CLR,metadata_gray$country))
# Analysis of Variance Table
# 
# Response: Distances
# Df  Sum Sq Mean Sq F value  Pr(>F)  
# Groups      2  1026.7  513.34  2.5591 0.08075 .
# Residuals 150 30089.1  200.59                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf(paste("F5.1_bacteria_age", "pdf", sep="."),width=14, height=12)
ggplot(metadata_gray,aes(pcoa1,pcoa2,color = category))+
  geom_point()+
  scale_color_manual( values = c("#ADFF2F",  "#F08080", "#778899"
                                 ))+
  
  stat_ellipse(show.legend = F)+
  geom_point(data = centroid_country, size=5,shape=21,color="black",aes(fill=category))+
  scale_fill_manual( values = c("#ADFF2F",  "#F08080", "#778899"))+
  theme_light()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = labs[1], y = labs[2])+guides(fill="none")
dev.off()

#==============================================
# Density of beta diversity analysis for population 

Euclidean_CLR <- as.matrix(Euclidean_CLR)
# Convert to three columns 
Euclidean_CLR <- as.data.frame(as.table(Euclidean_CLR))

# add category column to dist_b matrix
run_id_Y <- comsample_id[which(comsample_id$category %in% "Y"),"run_id"]
run_id_E <- comsample_id[which(comsample_id$category %in% "E"),"run_id"]
run_id_C <- comsample_id[which(comsample_id$category %in% "C"),"run_id"]


dist_b_country <- Euclidean_CLR %>% mutate(category = case_when(Var1 %in% run_id_Y & Var2 %in% run_id_Y  ~ "Y",
                                                                Var1 %in% run_id_E & Var2 %in% run_id_E ~ "E",
                                                                Var1 %in% run_id_C & Var2 %in% run_id_C  ~ "C",
                                                                TRUE ~"Not listed"))

dist_b_country <- dist_b_country[-which(dist_b_country$category %in% "Not listed"),] 

kruskal.test(Freq~category, data = dist_b_country)
# Kruskal-Wallis rank sum test
# 
# data:  Freq by category
# Kruskal-Wallis chi-squared = 533.4, df = 2, p-value < 2.2e-16

pdf(paste("F5.2bacteria", "pdf", sep="."),width=12, height=9)
ggdensity(dist_b_country,x='Freq',color ="category",fill = "category", add="mean", rug = FALSE,  palette = c("#ADFF2F",  "#F08080", "#778899"))+
  theme_light()+
  geom_density()+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")
dev.off()



#------------------------------------------------#
#                                                #
#           Beta for Archaea                     # 
#                                                #
#------------------------------------------------#


b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Archaea"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
count_data_zeroReplaced <- b_taxa[,-1:-5]
# Define pseudocount
pseudocount <- 0.5

# Impute zero
count_data_zeroReplaced[count_data_zeroReplaced == 0] <- pseudocount

# CLR transformation
count_data_CLR <- t(apply(count_data_zeroReplaced, 1, function(x) log(x) - mean(log(x))))
count_data_CLR <- count_data_CLR[,comsample_id$run_id]
# Calculate Euclidean distances
Euclidean_CLR <- vegdist(t(count_data_CLR), method = "euclidean")

pcoa <- cmdscale (Euclidean_CLR, eig = TRUE,add = TRUE) # add make eig convert to positive value 
positions <- pcoa$points
colnames(positions) <- c("pcoa1","pcoa2")
percent_explained <- 100*pcoa$eig / sum(pcoa$eig)
prett_pe <- format(round(percent_explained[1:2],digits=1),nsmall = 1,trim=TRUE)
labs <- c(glue("PCo 1 ({prett_pe[1]}%)"),
          glue("PCo 2 ({prett_pe[2]}%)"))
positions_samples <- positions %>% as_tibble(rownames="run_id")
metadata_gray <- inner_join(positions_samples,comsample_id)

# variance is country

centroid_country <- metadata_gray %>% group_by(category) %>%
  summarize(pcoa1 =mean(pcoa1),pcoa2 =mean(pcoa2))

## p_PERMANOVA
value_adnois <-adonis(Euclidean_CLR~metadata_gray$category,permutations = 999,method = "euclidean")
value_adnois$aov.tab  
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# metadata_gray$category   2     127.8  63.893  1.5093 0.01973  0.015 *
#   Residuals              150    6349.8  42.332         0.98027         
# Total                  152    6477.6                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# BETADISPER
anova(betadisper(Euclidean_CLR,metadata_gray$category)) # or use permutest(betadisper(Euclidean_CLR,metadata_gray$country))
# Analysis of Variance Table
# 
# Response: Distances
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Groups      2  29.75 14.8741  2.3789 0.09614 .
# Residuals 150 937.89  6.2526                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf(paste("F5.1_archaea_age", "pdf", sep="."),width=14, height=12)
ggplot(metadata_gray,aes(pcoa1,pcoa2,color = category))+
  geom_point()+
  scale_color_manual( values = c("#ADFF2F",  "#F08080", "#778899"
  ))+
  
  stat_ellipse(show.legend = F)+
  geom_point(data = centroid_country, size=5,shape=21,color="black",aes(fill=category))+
  scale_fill_manual( values = c("#ADFF2F",  "#F08080", "#778899"))+
  theme_light()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = labs[1], y = labs[2])+guides(fill="none")
dev.off()

#==============================================
# Density of beta diversity analysis for population 

Euclidean_CLR <- as.matrix(Euclidean_CLR)
# Convert to three columns 
Euclidean_CLR <- as.data.frame(as.table(Euclidean_CLR))

# add category column to dist_b matrix
run_id_Y <- comsample_id[which(comsample_id$category %in% "Y"),"run_id"]
run_id_E <- comsample_id[which(comsample_id$category %in% "E"),"run_id"]
run_id_C <- comsample_id[which(comsample_id$category %in% "C"),"run_id"]


dist_b_country <- Euclidean_CLR %>% mutate(category = case_when(Var1 %in% run_id_Y & Var2 %in% run_id_Y  ~ "Y",
                                                                Var1 %in% run_id_E & Var2 %in% run_id_E ~ "E",
                                                                Var1 %in% run_id_C & Var2 %in% run_id_C  ~ "C",
                                                                TRUE ~"Not listed"))

dist_b_country <- dist_b_country[-which(dist_b_country$category %in% "Not listed"),] 

kruskal.test(Freq~category, data = dist_b_country)
# Kruskal-Wallis rank sum test
# 
# data:  Freq by category
# Kruskal-Wallis chi-squared = 533.4, df = 2, p-value < 2.2e-16

pdf(paste("F5.2archaea", "pdf", sep="."),width=12, height=9)
ggdensity(dist_b_country,x='Freq',color ="category",fill = "category", add="mean", rug = FALSE,  palette = c("#ADFF2F",  "#F08080", "#778899"))+
  theme_light()+
  geom_density()+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")
dev.off()



#------------------------------------------------#
#                                                #
#           Beta for Eukaryota                     # 
#                                                #
#------------------------------------------------#


b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Eukaryota"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
count_data_zeroReplaced <- b_taxa[,-1:-5]
# Define pseudocount
pseudocount <- 0.5

# Impute zero
count_data_zeroReplaced[count_data_zeroReplaced == 0] <- pseudocount

# CLR transformation
count_data_CLR <- t(apply(count_data_zeroReplaced, 1, function(x) log(x) - mean(log(x))))
count_data_CLR <- count_data_CLR[,comsample_id$run_id]
# Calculate Euclidean distances
Euclidean_CLR <- vegdist(t(count_data_CLR), method = "euclidean")

pcoa <- cmdscale (Euclidean_CLR, eig = TRUE,add = TRUE) # add make eig convert to positive value 
positions <- pcoa$points
colnames(positions) <- c("pcoa1","pcoa2")
percent_explained <- 100*pcoa$eig / sum(pcoa$eig)
prett_pe <- format(round(percent_explained[1:2],digits=1),nsmall = 1,trim=TRUE)
labs <- c(glue("PCo 1 ({prett_pe[1]}%)"),
          glue("PCo 2 ({prett_pe[2]}%)"))
positions_samples <- positions %>% as_tibble(rownames="run_id")
metadata_gray <- inner_join(positions_samples,comsample_id)

# variance is country

centroid_country <- metadata_gray %>% group_by(category) %>%
  summarize(pcoa1 =mean(pcoa1),pcoa2 =mean(pcoa2))

## p_PERMANOVA
value_adnois <-adonis(Euclidean_CLR~metadata_gray$category,permutations = 999,method = "euclidean")
value_adnois$aov.tab  
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# metadata_gray$category   2     88.31  44.157    2.27 0.02938  0.001 ***
#   Residuals              150   2917.93  19.453         0.97062           
# Total                  152   3006.25                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# BETADISPER
anova(betadisper(Euclidean_CLR,metadata_gray$category)) # or use permutest(betadisper(Euclidean_CLR,metadata_gray$country))
# Analysis of Variance Table
# 
# Response: Distances
# Df Sum Sq Mean Sq F value Pr(>F)
# Groups      2   2.32  1.1597  0.4559 0.6347
# Residuals 150 381.56  2.5437   


pdf(paste("F5.1_Eukaryota_age", "pdf", sep="."),width=14, height=12)
ggplot(metadata_gray,aes(pcoa1,pcoa2,color = category))+
  geom_point()+
  scale_color_manual( values = c("#ADFF2F",  "#F08080", "#778899"
  ))+
  
  stat_ellipse(show.legend = F)+
  geom_point(data = centroid_country, size=5,shape=21,color="black",aes(fill=category))+
  scale_fill_manual( values = c("#ADFF2F",  "#F08080", "#778899"))+
  theme_light()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = labs[1], y = labs[2])+guides(fill="none")
dev.off()

#==============================================
# Density of beta diversity analysis for population 

Euclidean_CLR <- as.matrix(Euclidean_CLR)
# Convert to three columns 
Euclidean_CLR <- as.data.frame(as.table(Euclidean_CLR))

# add category column to dist_b matrix
run_id_Y <- comsample_id[which(comsample_id$category %in% "Y"),"run_id"]
run_id_E <- comsample_id[which(comsample_id$category %in% "E"),"run_id"]
run_id_C <- comsample_id[which(comsample_id$category %in% "C"),"run_id"]


dist_b_country <- Euclidean_CLR %>% mutate(category = case_when(Var1 %in% run_id_Y & Var2 %in% run_id_Y  ~ "Y",
                                                                Var1 %in% run_id_E & Var2 %in% run_id_E ~ "E",
                                                                Var1 %in% run_id_C & Var2 %in% run_id_C  ~ "C",
                                                                TRUE ~"Not listed"))

dist_b_country <- dist_b_country[-which(dist_b_country$category %in% "Not listed"),] 

kruskal.test(Freq~category, data = dist_b_country)
# Kruskal-Wallis rank sum test
# 
# data:  Freq by category
# Kruskal-Wallis chi-squared = 64.358, df = 2, p-value = 1.059e-14

pdf(paste("F5.2fungi", "pdf", sep="."),width=12, height=9)
ggdensity(dist_b_country,x='Freq',color ="category",fill = "category", add="mean", rug = FALSE,  palette = c("#ADFF2F",  "#F08080", "#778899"))+
  theme_light()+
  geom_density()+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")
dev.off()

#------------------------------------------------#
#                                                #
#           Beta for Viruses                     # 
#                                                #
#------------------------------------------------#


b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Viruses"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
count_data_zeroReplaced <- b_taxa[,-1:-5]
# Define pseudocount
pseudocount <- 0.5

# Impute zero
count_data_zeroReplaced[count_data_zeroReplaced == 0] <- pseudocount

# CLR transformation
count_data_CLR <- t(apply(count_data_zeroReplaced, 1, function(x) log(x) - mean(log(x))))
count_data_CLR <- count_data_CLR[,comsample_id$run_id]
# Calculate Euclidean distances
Euclidean_CLR <- vegdist(t(count_data_CLR), method = "euclidean")

pcoa <- cmdscale (Euclidean_CLR, eig = TRUE,add = TRUE) # add make eig convert to positive value 
positions <- pcoa$points
colnames(positions) <- c("pcoa1","pcoa2")
percent_explained <- 100*pcoa$eig / sum(pcoa$eig)
prett_pe <- format(round(percent_explained[1:2],digits=1),nsmall = 1,trim=TRUE)
labs <- c(glue("PCo 1 ({prett_pe[1]}%)"),
          glue("PCo 2 ({prett_pe[2]}%)"))
positions_samples <- positions %>% as_tibble(rownames="run_id")
metadata_gray <- inner_join(positions_samples,comsample_id)

# variance is country

centroid_country <- metadata_gray %>% group_by(category) %>%
  summarize(pcoa1 =mean(pcoa1),pcoa2 =mean(pcoa2))

## p_PERMANOVA
value_adnois <-adonis(Euclidean_CLR~metadata_gray$category,permutations = 999,method = "euclidean")
value_adnois$aov.tab  
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# metadata_gray$category   2      1069  534.53  2.2901 0.02963  0.001 ***
#   Residuals              150     35012  233.41         0.97037           
# Total                  152     36081                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# BETADISPER
anova(betadisper(Euclidean_CLR,metadata_gray$category)) # or use permutest(betadisper(Euclidean_CLR,metadata_gray$country))
# Analysis of Variance Table
# 
# Response: Distances
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Groups      2  222.2 111.102  3.1461 0.04587 *
#   Residuals 150 5297.2  35.315                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf(paste("F5.1_Viruses_age", "pdf", sep="."),width=14, height=12)
ggplot(metadata_gray,aes(pcoa1,pcoa2,color = category))+
  geom_point()+
  scale_color_manual( values = c("#ADFF2F",  "#F08080", "#778899"
  ))+
  
  stat_ellipse(show.legend = F)+
  geom_point(data = centroid_country, size=5,shape=21,color="black",aes(fill=category))+
  scale_fill_manual( values = c("#ADFF2F",  "#F08080", "#778899"))+
  theme_light()+
  #scale_fill_manual(values=c("#a6cee3","#C2C2C2","#b2df8a","#d9544d"),label = c("Italy_2","China","Japan","Italy_1"))+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = labs[1], y = labs[2])+guides(fill="none")
dev.off()



#==============================================
# Density of beta diversity analysis for population 

Euclidean_CLR <- as.matrix(Euclidean_CLR)
# Convert to three columns 
Euclidean_CLR <- as.data.frame(as.table(Euclidean_CLR))

# add category column to dist_b matrix
run_id_Y <- comsample_id[which(comsample_id$category %in% "Y"),"run_id"]
run_id_E <- comsample_id[which(comsample_id$category %in% "E"),"run_id"]
run_id_C <- comsample_id[which(comsample_id$category %in% "C"),"run_id"]


dist_b_country <- Euclidean_CLR %>% mutate(category = case_when(Var1 %in% run_id_Y & Var2 %in% run_id_Y  ~ "Y",
                                                               Var1 %in% run_id_E & Var2 %in% run_id_E ~ "E",
                                                               Var1 %in% run_id_C & Var2 %in% run_id_C  ~ "C",
                                                               TRUE ~"Not listed"))

dist_b_country <- dist_b_country[-which(dist_b_country$category %in% "Not listed"),] 

kruskal.test(Freq~category, data = dist_b_country)
# Kruskal-Wallis rank sum test
# 
# data:  Freq by category
# Kruskal-Wallis chi-squared = 520.87, df = 2, p-value < 2.2e-16

pdf(paste("F5.2viruses", "pdf", sep="."),width=12, height=9)
ggdensity(dist_b_country,x='Freq',color ="category",fill = "category", add="mean", rug = FALSE,  palette = c("#ADFF2F",  "#F08080", "#778899"))+
  theme_light()+
  geom_density()+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=24,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"))+
  labs(x = "", y = "")
dev.off()


#========================================================END

