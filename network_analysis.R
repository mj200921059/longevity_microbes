# Microbes composition and network building

Cb64k <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
           "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
           "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
           "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
           "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
           "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
           "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
           "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
           "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
           "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
           "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
           "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
           "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")






#================================================
#$
#
# microbes composition 
#
#================================================
load("kraken_taxa_corrected.rda")
load('comsample_id.rda')
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Archaea"),]
#b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
sid <- comsample_id[,c("run_id","category")]
colnames(sid) <- c("variable","category")
b_taxa <- b_taxa[,-1:-2]
b_taxa_m <- melt(b_taxa)
b_taxa_m_s <- merge(b_taxa_m,sid, by="variable")

# # mean of taxon_rel_abund 
# taxon_rel_abund <- b_taxa_m_s %>% 
#   group_by(variable) %>%
#   mutate(rel_abund = value / sum(value)) %>%
#   ungroup() %>%
#   group_by(category,variable, phylum) %>%
#   summarize(rel_abund = sum(rel_abund),.groups = "drop") %>%
#   group_by(category,phylum) %>%
#   summarize(mean_rel_abund = 100*mean(rel_abund),.groups = "drop")

# medain of taxon_rel_abund 
taxon_rel_abund <- b_taxa_m_s %>% 
  group_by(variable) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  group_by(category,variable, phylum) %>%
  summarize(rel_abund = sum(rel_abund),.groups = "drop") %>% # sum every phylum abundance in samples of each group.
  group_by(category,phylum) %>%
  summarize(median_rel_abund = 100*median(rel_abund),.groups = "drop") # median phylum abundance among samples in each group.


# 
taxon_rel_abund %>%
  group_by(phylum) %>%
  summarize(max=max(median_rel_abund)) %>%
  arrange(desc(max))


taxon_pool <- taxon_rel_abund %>%
  group_by(phylum) %>%
  #summarize(max=max(mean_rel_abund)) %>%
  #arrange(desc(max))
  summarize(pool = max(median_rel_abund) <0.5,
            mean = median(median_rel_abund),
            .groups = "drop" )



inner_join(taxon_rel_abund,taxon_pool,by="phylum") %>%
  mutate(taxon= if_else(pool, "Other", phylum)) %>%
  group_by(category, taxon) %>%
  summarize(median_rel_abund =sum(median_rel_abund),
            mean= min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon), # REORDER LEGEND
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n=1)) %>%  #pull(taxon) %>% levels() # TO CHECK THE TAXON ORDER 
  
  ggplot(aes(x=category,y=median_rel_abund,fill=taxon))+
  geom_col(position="fill")+
  #scale_fill_discrete(name=NULL)+
  scale_fill_manual(name=NULL,
                    breaks=c('Euryarchaeota','Candidatus Thermoplasmatota','Thermoproteota',
                             'Other'),
                    values= c(brewer.pal(3,'Dark2'),"gray"))+
  scale_y_continuous(expand = c(0,0),labels = percent_format())+
  labs(x=NULL,
       y="Median Relative Abundance")+
  theme_classic()+
  theme(legend.text = element_text(face = 'italic'),
        legend.key.size = unit(10,"pt"))

ggsave("Fig6.1archaea_abundance_phlyum.pdf", width = 5,height = 4)



# ---------------------------Bacteria
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]


#b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
sid <- comsample_id[,c("run_id","category")]
colnames(sid) <- c("variable","category")
b_taxa <- b_taxa[,-1:-2]
b_taxa_m <- melt(b_taxa)
b_taxa_m_s <- merge(b_taxa_m,sid, by="variable")

# # mean of taxon_rel_abund 
# taxon_rel_abund <- b_taxa_m_s %>% 
#   group_by(variable) %>%
#   mutate(rel_abund = value / sum(value)) %>%
#   ungroup() %>%
#   group_by(category,variable, phylum) %>%
#   summarize(rel_abund = sum(rel_abund),.groups = "drop") %>%
#   group_by(category,phylum) %>%
#   summarize(mean_rel_abund = 100*mean(rel_abund),.groups = "drop")

# medain of taxon_rel_abund 
taxon_rel_abund <- b_taxa_m_s %>% 
  group_by(variable) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  group_by(category,variable, phylum) %>%
  summarize(rel_abund = sum(rel_abund),.groups = "drop") %>% # sum every phylum abundance in samples of each group.
  group_by(category,phylum) %>%
  summarize(median_rel_abund = 100*median(rel_abund),.groups = "drop") # median phylum abundance among samples in each group.


# 
taxon_rel_abund %>%
  group_by(phylum) %>%
  summarize(max=max(median_rel_abund)) %>%
  arrange(desc(max))


taxon_pool <- taxon_rel_abund %>%
  group_by(phylum) %>%
  #summarize(max=max(mean_rel_abund)) %>%
  #arrange(desc(max))
  summarize(pool = max(median_rel_abund) <0.5,
            mean = median(median_rel_abund),
            .groups = "drop" )


#set.seed(2025)
inner_join(taxon_rel_abund,taxon_pool,by="phylum") %>%
  mutate(taxon= if_else(pool, "Other", phylum)) %>%
  group_by(category, taxon) %>%
  summarize(median_rel_abund =sum(median_rel_abund),
            mean= min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon), # REORDER LEGEND
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n=1)) %>%  #pull(taxon) %>% levels() # TO CHECK THE TAXON ORDER 
  
  ggplot(aes(x=category,y=median_rel_abund,fill=taxon))+
  geom_col(position="fill")+
  #scale_fill_discrete(name=NULL)+
  scale_fill_manual(name=NULL,
                    breaks=c('Bacillota','Bacteroidota','Actinomycetota',
                             'Pseudomonadota',
                             'Other'),
                    values= c(sample(Cb64k,4),"gray"))+
  scale_y_continuous(expand = c(0,0),labels = percent_format())+
  labs(x=NULL,
       y="Median Relative Abundance")+
  theme_classic()+
  theme(legend.text = element_text(face = 'italic'),
        legend.key.size = unit(10,"pt"))

ggsave("Fig6.1bacteria_abundance_phlyum.pdf", width = 5,height = 4)







# ---------------------------Eukaryota
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Eukaryota"),]


#b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
sid <- comsample_id[,c("run_id","category")]
colnames(sid) <- c("variable","category")
b_taxa <- b_taxa[,-1:-2]
b_taxa_m <- melt(b_taxa)
b_taxa_m_s <- merge(b_taxa_m,sid, by="variable")

# # mean of taxon_rel_abund 
# taxon_rel_abund <- b_taxa_m_s %>% 
#   group_by(variable) %>%
#   mutate(rel_abund = value / sum(value)) %>%
#   ungroup() %>%
#   group_by(category,variable, phylum) %>%
#   summarize(rel_abund = sum(rel_abund),.groups = "drop") %>%
#   group_by(category,phylum) %>%
#   summarize(mean_rel_abund = 100*mean(rel_abund),.groups = "drop")

# medain of taxon_rel_abund 
taxon_rel_abund <- b_taxa_m_s %>% 
  group_by(variable) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  group_by(category,variable, phylum) %>%
  summarize(rel_abund = sum(rel_abund),.groups = "drop") %>% # sum every phylum abundance in samples of each group.
  group_by(category,phylum) %>%
  summarize(median_rel_abund = 100*median(rel_abund),.groups = "drop") # median phylum abundance among samples in each group.


# 
taxon_rel_abund %>%
  group_by(phylum) %>%
  summarize(max=max(median_rel_abund)) %>%
  arrange(desc(max))


taxon_pool <- taxon_rel_abund %>%
  group_by(phylum) %>%
  #summarize(max=max(mean_rel_abund)) %>%
  #arrange(desc(max))
  summarize(pool = max(median_rel_abund) <0.5,
            mean = median(median_rel_abund),
            .groups = "drop" )


set.seed(203)
inner_join(taxon_rel_abund,taxon_pool,by="phylum") %>%
  mutate(taxon= if_else(pool, "Other", phylum)) %>%
  group_by(category, taxon) %>%
  summarize(median_rel_abund =sum(median_rel_abund),
            mean= min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon), # REORDER LEGEND
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n=1)) %>%  #pull(taxon) %>% levels() # TO CHECK THE TAXON ORDER 
  
  ggplot(aes(x=category,y=median_rel_abund,fill=taxon))+
  geom_col(position="fill")+
  #scale_fill_discrete(name=NULL)+
  scale_fill_manual(name=NULL,
                    breaks=c('Chordata','Ascomycota',
                             'Other'),
                    values= c(sample(Cb64k,2),"gray"))+
  scale_y_continuous(expand = c(0,0),labels = percent_format())+
  labs(x=NULL,
       y="Median Relative Abundance")+
  theme_classic()+
  theme(legend.text = element_text(face = 'italic'),
        legend.key.size = unit(10,"pt"))

ggsave("Fig6.1fungi_abundance_phlyum.pdf", width = 5,height = 4)




# ---------------------------Viruses
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Viruses"),]


#b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
sid <- comsample_id[,c("run_id","category")]
colnames(sid) <- c("variable","category")
b_taxa <- b_taxa[,-1:-2]
b_taxa_m <- melt(b_taxa)
b_taxa_m_s <- merge(b_taxa_m,sid, by="variable")

# # mean of taxon_rel_abund 
# taxon_rel_abund <- b_taxa_m_s %>% 
#   group_by(variable) %>%
#   mutate(rel_abund = value / sum(value)) %>%
#   ungroup() %>%
#   group_by(category,variable, phylum) %>%
#   summarize(rel_abund = sum(rel_abund),.groups = "drop") %>%
#   group_by(category,phylum) %>%
#   summarize(mean_rel_abund = 100*mean(rel_abund),.groups = "drop")

# medain of taxon_rel_abund 
taxon_rel_abund <- b_taxa_m_s %>% 
  group_by(variable) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  group_by(category,variable, phylum) %>%
  summarize(rel_abund = sum(rel_abund),.groups = "drop") %>% # sum every phylum abundance in samples of each group.
  group_by(category,phylum) %>%
  summarize(median_rel_abund = 100*median(rel_abund),.groups = "drop") # median phylum abundance among samples in each group.


# 
taxon_rel_abund %>%
  group_by(phylum) %>%
  summarize(max=max(median_rel_abund)) %>%
  arrange(desc(max))


taxon_pool <- taxon_rel_abund %>%
  group_by(phylum) %>%
  #summarize(max=max(mean_rel_abund)) %>%
  #arrange(desc(max))
  summarize(pool = max(median_rel_abund) <0.5,
            mean = median(median_rel_abund),
            .groups = "drop" )


set.seed(2034)
inner_join(taxon_rel_abund,taxon_pool,by="phylum") %>%
  mutate(taxon= if_else(pool, "Other", phylum)) %>%
  group_by(category, taxon) %>%
  summarize(median_rel_abund =sum(median_rel_abund),
            mean= min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon), # REORDER LEGEND
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n=1)) %>%  #pull(taxon) %>% levels() # TO CHECK THE TAXON ORDER 
  
  ggplot(aes(x=category,y=median_rel_abund,fill=taxon))+
  geom_col(position="fill")+
  #scale_fill_discrete(name=NULL)+
  scale_fill_manual(name=NULL,
                    breaks=c('Uroviricota',
                             'Other'),
                    values= c("#7A87A1","gray"))+
  scale_y_continuous(expand = c(0,0),labels = percent_format())+
  labs(x=NULL,
       y="Median Relative Abundance")+
  theme_classic()+
  theme(legend.text = element_text(face = 'italic'),
        legend.key.size = unit(10,"pt"))

ggsave("Fig6.1virus_abundance_phlyum.pdf", width = 5,height = 4)


#===============================================================
#
# Figure for species
#
#===============================================================


#================================================
load("kraken_taxa_corrected.rda")
load('comsample_id.rda')
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Archaea"),]
#b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
sid <- comsample_id[,c("run_id","category")]
colnames(sid) <- c("variable","category")
b_taxa <- b_taxa[,-1:-2]
b_taxa_m <- melt(b_taxa)
b_taxa_m_s <- merge(b_taxa_m,sid, by="variable")

# taxon_rel_abund 
taxon_rel_abund <- b_taxa_m_s %>% 
  group_by(variable) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  group_by(category,variable, species) %>%
  summarize(rel_abund = 100*sum(rel_abund),.groups = "drop") 

# aa <- taxon_rel_abund[which(taxon_rel_abund$rel_abund >0),]
# 
# 
# # 
# taxon_rel_abund %>%
#   group_by(species) %>%
#   summarize(max=max(rel_abund)) %>%
#   arrange(desc(max))

aa <- taxon_rel_abund %>%
  group_by(category,species) %>%
  summarize(mean=mean(rel_abund),.groups = 'drop') %>%
  group_by(species) %>%
  summarize(pool = max(mean)) %>%
  arrange(desc(pool))



taxon_pool <- taxon_rel_abund %>%
  group_by(category,species) %>%
  summarize(mean=mean(rel_abund),.groups = 'drop') %>%
  group_by(species) %>%
  summarize(pool = max(mean) <1.7,
            mean = mean(mean),
            .groups = "drop" )

sample_order <- taxon_rel_abund %>%
  filter(species == "Methanobrevibacter smithii") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(variable,order)

g_data <- inner_join(taxon_rel_abund,taxon_pool,by="species") %>%
  mutate(taxon= if_else(pool, "Other", species)) %>%
  group_by(variable,category, taxon) %>%
  summarize(rel_abund =sum(rel_abund),
            mean= min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon), # REORDER LEGEND
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n=1))  %>%  #pull(taxon) %>% levels() # TO CHECK THE TAXON ORDER 
  inner_join(.,sample_order,by="variable") %>%
  mutate(variable = factor(variable),
         variable =fct_reorder(variable, order))
g_data <- g_data[-which(g_data$taxon %in% "Other"),]
g_data %>% 
  ggplot(aes(x=variable,y=rel_abund,fill=taxon))+
  geom_col(width = 1)+
  #scale_fill_discrete(name=NULL)+
  # scale_fill_manual(name=NULL,
  #                   breaks=c('Euryarchaeota','Candidatus Thermoplasmatota','Thermoproteota',
  #                            'Other'),
  #                   values= c(brewer.pal(3,'Dark2'),"gray"))+
  scale_fill_manual(name=NULL, values = Cb64k[32:52])+
  scale_y_continuous(expand = c(0,0))+
  # facet_wrap(~ category, nrow=1, scale='free_x')+
  facet_grid(~ category, scale='free_x',space = "free",switch = "x")+
  labs(x=NULL,
       y="Relative Abundance(%)")+
  theme_classic()+
  theme(axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = 'italic'),
        legend.key.size = unit(10,"pt"),
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave("Fig6.2archaea_abundance_speices.pdf", width = 5,height = 4)





#================================================Bacteria
load("kraken_taxa_corrected.rda")
load('comsample_id.rda')
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
sid <- comsample_id[,c("run_id","category")]
colnames(sid) <- c("variable","category")
b_taxa <- b_taxa[,-1:-2]
b_taxa_m <- melt(b_taxa)
b_taxa_m_s <- merge(b_taxa_m,sid, by="variable")

# taxon_rel_abund 
taxon_rel_abund <- b_taxa_m_s %>% 
  group_by(variable) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  group_by(category,variable, species) %>%
  summarize(rel_abund = 100*sum(rel_abund),.groups = "drop") 

# aa <- taxon_rel_abund[which(taxon_rel_abund$rel_abund >0),]
# 
# 
# # 
# taxon_rel_abund %>%
#   group_by(species) %>%
#   summarize(max=max(rel_abund)) %>%
#   arrange(desc(max))

aa <- taxon_rel_abund %>%
  group_by(category,species) %>%
  summarize(mean=mean(rel_abund),.groups = 'drop') %>%
  group_by(species) %>%
  summarize(pool = max(mean)) %>%
  arrange(desc(pool))



taxon_pool <- taxon_rel_abund %>%
  group_by(category,species) %>%
  summarize(mean=mean(rel_abund),.groups = 'drop') %>%
  group_by(species) %>%
  summarize(pool = max(mean) <1.68,
            mean = mean(mean),
            .groups = "drop" )

sample_order <- taxon_rel_abund %>%
  filter(species == "Faecalibacterium prausnitzii") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(variable,order)

g_data <- inner_join(taxon_rel_abund,taxon_pool,by="species") %>%
  mutate(taxon= if_else(pool, "Other", species)) %>%
  group_by(variable,category, taxon) %>%
  summarize(rel_abund =sum(rel_abund),
            mean= min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon), # REORDER LEGEND
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n=1))  %>%  #pull(taxon) %>% levels() # TO CHECK THE TAXON ORDER 
  inner_join(.,sample_order,by="variable") %>%
  mutate(variable = factor(variable),
         variable =fct_reorder(variable, order))
g_data <- g_data[-which(g_data$taxon %in% "Other"),]
g_data %>% 
  ggplot(aes(x=variable,y=rel_abund,fill=taxon))+
  geom_col(width = 1)+
  #scale_fill_discrete(name=NULL)+
  # scale_fill_manual(name=NULL,
  #                   breaks=c('Euryarchaeota','Candidatus Thermoplasmatota','Thermoproteota',
  #                            'Other'),
  #                   values= c(brewer.pal(3,'Dark2'),"gray"))+
  scale_fill_manual(name=NULL, values = Cb64k)+
  scale_y_continuous(expand = c(0,0))+
  # facet_wrap(~ category, nrow=1, scale='free_x')+
  facet_grid(~ category, scale='free_x',space = "free",switch = "x")+
  labs(x=NULL,
       y="Relative Abundance(%)")+
  theme_classic()+
  theme(axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = 'italic'),
        legend.key.size = unit(10,"pt"),
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave("Fig6.2bacteria_abundance_speices.pdf", width = 5,height = 4)


#================================================Eukaryota

b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Eukaryota"),]

sid <- comsample_id[,c("run_id","category")]
colnames(sid) <- c("variable","category")
b_taxa <- b_taxa[,-1:-2]
b_taxa_m <- melt(b_taxa)
b_taxa_m_s <- merge(b_taxa_m,sid, by="variable")
b_taxa_m_s <- b_taxa_m_s[-which(b_taxa_m_s$species %in% "Homo sapiens"),]
# taxon_rel_abund 
taxon_rel_abund <- b_taxa_m_s %>% 
  group_by(variable) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  group_by(category,variable, species) %>%
  summarize(rel_abund = 100*sum(rel_abund),.groups = "drop") 

# aa <- taxon_rel_abund[which(taxon_rel_abund$rel_abund >0),]
# 
# 
# # 
# taxon_rel_abund %>%
#   group_by(species) %>%
#   summarize(max=max(rel_abund)) %>%
#   arrange(desc(max))

aa <- taxon_rel_abund %>%
  group_by(category,species) %>%
  summarize(mean=mean(rel_abund),.groups = 'drop') %>%
  group_by(species) %>%
  summarize(pool = max(mean)) %>%
  arrange(desc(pool))



taxon_pool <- taxon_rel_abund %>%
  group_by(category,species) %>%
  summarize(mean=mean(rel_abund),.groups = 'drop') %>%
  group_by(species) %>%
  summarize(pool = max(mean) <2,
            mean = mean(mean),
            .groups = "drop" )

sample_order <- taxon_rel_abund %>%
  filter(species == "Nakaseomyces glabratus") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(variable,order)

g_data <- inner_join(taxon_rel_abund,taxon_pool,by="species") %>%
  mutate(taxon= if_else(pool, "Other", species)) %>%
  group_by(variable,category, taxon) %>%
  summarize(rel_abund =sum(rel_abund),
            mean= min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon), # REORDER LEGEND
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n=1))  %>%  #pull(taxon) %>% levels() # TO CHECK THE TAXON ORDER 
  inner_join(.,sample_order,by="variable") %>%
  mutate(variable = factor(variable),
         variable =fct_reorder(variable, order))
g_data <- g_data[-which(g_data$taxon %in% "Other"),]
g_data %>% 
  ggplot(aes(x=variable,y=rel_abund,fill=taxon))+
  geom_col(width = 1)+
  #scale_fill_discrete(name=NULL)+
  # scale_fill_manual(name=NULL,
  #                   breaks=c('Euryarchaeota','Candidatus Thermoplasmatota','Thermoproteota',
  #                            'Other'),
  #                   values= c(brewer.pal(3,'Dark2'),"gray"))+
  scale_fill_manual(name=NULL, values = Cb64k[20:40])+
  scale_y_continuous(expand = c(0,0))+
  # facet_wrap(~ category, nrow=1, scale='free_x')+
  facet_grid(~ category, scale='free_x',space = "free",switch = "x")+
  labs(x=NULL,
       y="Relative Abundance(%)")+
  theme_classic()+
  theme(axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = 'italic'),
        legend.key.size = unit(10,"pt"),
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave("Fig6.2fungi_abundance_speices.pdf", width = 5,height = 4)



#================================================Viruses

b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Viruses"),]

sid <- comsample_id[,c("run_id","category")]
colnames(sid) <- c("variable","category")
b_taxa <- b_taxa[,-1:-2]
b_taxa_m <- melt(b_taxa)
b_taxa_m_s <- merge(b_taxa_m,sid, by="variable")

# taxon_rel_abund 
taxon_rel_abund <- b_taxa_m_s %>% 
  group_by(variable) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  group_by(category,variable, species) %>%
  summarize(rel_abund = 100*sum(rel_abund),.groups = "drop") 

# aa <- taxon_rel_abund[which(taxon_rel_abund$rel_abund >0),]
# 
# 
# # 
# taxon_rel_abund %>%
#   group_by(species) %>%
#   summarize(max=max(rel_abund)) %>%
#   arrange(desc(max))

aa <- taxon_rel_abund %>%
  group_by(category,species) %>%
  summarize(mean=mean(rel_abund),.groups = 'drop') %>%
  group_by(species) %>%
  summarize(pool = max(mean)) %>%
  arrange(desc(pool))



taxon_pool <- taxon_rel_abund %>%
  group_by(category,species) %>%
  summarize(mean=mean(rel_abund),.groups = 'drop') %>%
  group_by(species) %>%
  summarize(pool = max(mean) <1.7,
            mean = mean(mean),
            .groups = "drop" )

sample_order <- taxon_rel_abund %>%
  filter(species == "CrAss-like virus sp.") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(variable,order)

g_data <- inner_join(taxon_rel_abund,taxon_pool,by="species") %>%
  mutate(taxon= if_else(pool, "Other", species)) %>%
  group_by(variable,category, taxon) %>%
  summarize(rel_abund =sum(rel_abund),
            mean= min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon), # REORDER LEGEND
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n=1))  %>%  #pull(taxon) %>% levels() # TO CHECK THE TAXON ORDER 
  inner_join(.,sample_order,by="variable") %>%
  mutate(variable = factor(variable),
         variable =fct_reorder(variable, order))
g_data <- g_data[-which(g_data$taxon %in% "Other"),]
g_data %>% 
  ggplot(aes(x=variable,y=rel_abund,fill=taxon))+
  geom_col(width = 1)+
  #scale_fill_discrete(name=NULL)+
  # scale_fill_manual(name=NULL,
  #                   breaks=c('Euryarchaeota','Candidatus Thermoplasmatota','Thermoproteota',
  #                            'Other'),
  #                   values= c(brewer.pal(3,'Dark2'),"gray"))+
  scale_fill_manual(name=NULL, values = Cb64k[50:70])+
  scale_y_continuous(expand = c(0,0))+
  # facet_wrap(~ category, nrow=1, scale='free_x')+
  facet_grid(~ category, scale='free_x',space = "free",switch = "x")+
  labs(x=NULL,
       y="Relative Abundance(%)")+
  theme_classic()+
  theme(axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = 'italic'),
        legend.key.size = unit(10,"pt"),
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave("Fig6.2viruses_abundance_speices.pdf", width = 5,height = 4)


#===============================================================
#
# Figure for difference abundance analysis
#
#===============================================================



#----------------------------------------------------------------------------
#
# Maaslin2 for archaea 
#----------------------------------------------------------------------------

set.seed(2023)
load("kraken_taxa_corrected.rda")
load('comsample_id.rda')
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Archaea"),]

rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]

# ===========================Y-E
YE_id <- comsample_id[which(comsample_id$category %in% c('Y','E')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_archaea_ye <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_archaea_ye",
  transform = "LOG",
  fixed_effects = "category",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,Y",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)

# ===========================Y-C
YE_id <- comsample_id[which(comsample_id$category %in% c('Y','C')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_archaea_yc <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_archaea_yc",
  transform = "LOG",
  fixed_effects = "category",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,Y",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)

# ===========================C-E
YE_id <- comsample_id[which(comsample_id$category %in% c('E','C')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_archaea_ce <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_archaea_ce",
  transform = "LOG",
  fixed_effects = c("category"),
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,E",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)





#----------------------------------------------------------------------------
#
# Maaslin2 FOR fungi
#----------------------------------------------------------------------------
set.seed(2023)
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Eukaryota"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]

# ===========================Y-E
YE_id <- comsample_id[which(comsample_id$category %in% c('Y','E')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_fungi_ye <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_fungi_ye",
  transform = "LOG",
  fixed_effects = "category",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,Y",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)






# ===========================Y-C
YE_id <- comsample_id[which(comsample_id$category %in% c('Y','C')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_fungi_yc <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_fungi_yc",
  transform = "LOG",
  fixed_effects = "category",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,Y",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)



# ===========================E-C
YE_id <- comsample_id[which(comsample_id$category %in% c('E','C')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_fungi_ce <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_fungi_ce",
  transform = "LOG",
  fixed_effects = c("category"),
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,E",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)


#----------------------------------------------------------------------------
#
# Maaslin2 for viruses
#----------------------------------------------------------------------------
set.seed(2023)
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Viruses"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]

# ===========================Y-E
YE_id <- comsample_id[which(comsample_id$category %in% c('Y','E')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_viruses_ye <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_viruses_ye",
  transform = "LOG",
  fixed_effects = "category",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,Y",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)



# ===========================Y-C
YE_id <- comsample_id[which(comsample_id$category %in% c('Y','C')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_viruses_yc <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_viruses_yc",
  transform = "LOG",
  fixed_effects = "category",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,Y",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)



# ===========================E-C
YE_id <- comsample_id[which(comsample_id$category %in% c('E','C')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_viruses_ce <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_viruses_ce",
  transform = "LOG",
  fixed_effects = c("category"),
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,E",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)



#----------------------------------------------------------------------------
#
# Maaslin2 for bacteria
#----------------------------------------------------------------------------

set.seed(2023)
b_taxa <- taxa_corrected[which(taxa_corrected$superkingdom %in% "Bacteria"),]
rownames(b_taxa) <- b_taxa$taxonomy_id
b_taxa_a <- b_taxa[,-1:-5]

# ===========================Y-E
YE_id <- comsample_id[which(comsample_id$category %in% c('Y','E')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_bacteria_ye <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_bacteria_ye",
  transform = "LOG",
  fixed_effects = "category",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,Y",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)



# ===========================Y-C
YE_id <- comsample_id[which(comsample_id$category %in% c('Y','C')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_bacteria_yc <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_bacteria_yc",
  transform = "LOG",
  fixed_effects = "category",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,Y",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)



# ===========================E-C
YE_id <- comsample_id[which(comsample_id$category %in% c('E','C')),]
taxa_ye <- b_taxa_a[,YE_id$run_id]
rownames(YE_id) <- YE_id$run_id
masslin2_bacteria_ce <- Maaslin2(
  as.data.frame(t(taxa_ye)),
  #taxa_ye,
  YE_id,
  # A folder will be created that is called like the below specified output.
  # It contains also figures to visualize the difference between genera 
  # for the significant ones.
  output = "masslin2_bacteria_ce",
  transform = "LOG",
  fixed_effects = c("category"),
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "category,E",  
  normalization = "TSS",# total-sum scaling
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)



#==============================================================================
#
# Merge all significant difference microbes from masslin2 folders into a excel file
# 
# HEATMAP
#==============================================================================





set.seed(2023)
load("kraken_taxa_corrected.rda")
load('comsample_id.rda')

# Check value
# masslin2_bacteria_ce$results %>%
#   filter(qval < 0.25) %>%
#   knitr::kable()




bacteria_ce <- masslin2_bacteria_ce$results %>%
  filter(qval < 0.25) 
bacteria_ce$micro <- "bacteria"
bacteria_ce$group <- "E_C"
bacteria_ye <- masslin2_bacteria_ye$results %>%
  filter(qval < 0.25) 
bacteria_ye$micro <- "bacteria"
bacteria_ye$group <- "Y_E"
bacteria_yc <- masslin2_bacteria_yc$results %>%
  filter(qval < 0.25) 
bacteria_yc$micro <- "bacteria"
bacteria_yc$group <- "Y_C"


archaea_ce <- masslin2_archaea_ce$results %>%
  filter(qval < 0.25) 
archaea_ce$micro <- "archaea"
archaea_ce$group <- "E_C"
archaea_ye <- masslin2_archaea_ye$results %>%
  filter(qval < 0.25) 
archaea_ye$micro <- "archaea"
archaea_ye$group <- "Y_E"
archaea_yc <- masslin2_archaea_yc$results %>%
  filter(qval < 0.25) 
archaea_yc$micro <- "archaea"
archaea_yc$group <- "Y_C"





fungi_ce <- masslin2_fungi_ce$results %>%
  filter(qval < 0.25) 
fungi_ce$micro <- "fungi"
fungi_ce$group <- "E_C"
fungi_ye <- masslin2_fungi_ye$results %>%
  filter(qval < 0.25) 
fungi_ye$micro <- "fungi"
fungi_ye$group <- "Y_E"
fungi_yc <- masslin2_fungi_yc$results %>%
  filter(qval < 0.25) 
fungi_yc$micro <- "fungi"
fungi_yc$group <- "Y_C"


viruses_ce <- masslin2_viruses_ce$results %>%
  filter(qval < 0.25) 
viruses_ce$micro <- "viruses"
viruses_ce$group <- "E_C"
viruses_ye <- masslin2_viruses_ye$results %>%
  filter(qval < 0.25) 
viruses_ye$micro <- "viruses"
viruses_ye$group <- "Y_E"
viruses_yc <- masslin2_viruses_yc$results %>%
  filter(qval < 0.25) 
viruses_yc$micro <- "viruses"
viruses_yc$group <- "Y_C"

# 
# masslin2_viruses_ye$results %>%
#   filter(qval < 0.25) %>%
#   knitr::kable()


hmap <- rbind(bacteria_ce,bacteria_ye,bacteria_yc,
              archaea_ye,archaea_yc,
              fungi_ce,fungi_ye,fungi_yc,
              viruses_ce,viruses_ye,viruses_yc)
save(hmap, file = "hmap.rda")





# number of DA plot using ggplot
# GROUP SAMPLE BASED ON AGE 
hmap_da <- hmap %>% mutate(category = case_when(coef < 0  ~ "Negative",
                                                             coef >0  ~ "Positive",
                                                             TRUE ~"Not listed"))


hmap_da_num <- hmap_da[,c("group","micro","category")]


df2 <- hmap_da_num %>% group_by(group,micro,category) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

df2_m<- df2[which(df2$category %in% "Positive"),]
df2_f <- df2[which(df2$category %in% "Negative"),]
# df2_f$total_count <- df2_f$total_count *-1


# Male population data graph
pdf(paste("F7.1num_da", "pdf", sep="."),width=12, height=9)
ggplot()+
  geom_bar(data =df2_m,aes(group,total_count,fill=micro),stat = "identity",position="stack",color="black",width = 0.7,size=0.25)+
  
  # Famale population data graph
  
  
  geom_bar(data =df2_f,aes(group,y= total_count *(-1),fill=micro),stat = "identity",position="stack",color="black",width = 0.7,size=0.25)+
  coord_flip()+
  scale_y_continuous(limits = c(-600,500),
                     breaks = seq(-600,500,100),
                     labels = abs(seq(-600,500,100)))+
  labs(x = "", y = "")+
  # scale_fill_manual(values=Cb64k[45:49],label = c("Archaea","Bacteria","Fungi","Viruses"))+
  scale_fill_manual(values = c("#4a2377","#8cc5e3","#f55f74","#0d7d87"),label = c("Archaea","Bacteria","Fungi","Viruses"))+
  #scale_fill_manual(values = c("#3594cc","#8cc5e3","#ea801c","#f0b077"),label = c("Archaea","Bacteria","Fungi","Viruses"))+
  # display themes
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 24, face = "bold",colour = 'black'),  # fonts/sizes
    axis.title = element_text(size = 24, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5), 
    plot.caption = element_text(hjust=0, size=11, face = "italic")
  )
dev.off()

# use to make venn figure for da of bacteria 

X <- list(Y_E =bacteria_ye$feature,
          Y_C =bacteria_yc$feature,
          E_C = bacteria_ce$feature)


pdf(paste("F7.2", "pdf", sep="."),width=12, height=9)
ggVennDiagram(X,label_alpha = 0)+ggplot2::scale_fill_gradient(low="lightblue",high = 'yellow')
dev.off()

aa <- hmap[which(hmap$micro %in% "bacteria"),]
a <- as.data.frame(table(aa$feature))
# five common bacterium in YE,YC,CE difference abundance analysis
# 1623
# # Ligilactobacillus ruminis
# 1712675
# Turicibacter sp. H121
# 1737424
# #Blautia massiliensis
# 1902
# Streptomyces coelicolor
# 649756
# Anaerostipes hadrus


#--------------------------------------------------------------
# heatmap for selected DA microbes using pheatmap


load("hmap.rda")
hmap$feature <- str_remove(hmap$feature,"X")

rownames(taxa_corrected) <- taxa_corrected$taxonomy_id

bb <- hmap[-which(hmap$micro %in% "bacteria"),]

damicro <- c("1623","1712675","1737424","1902","649756",unique(bb$feature)) # 5 common DA Bacteria+ all other DA microbes

#hmap_abundance <- taxa_corrected[which(taxa_corrected$taxonomy_id %in% damicro),comsample_id$run_id]
cc <- comsample_id[order(comsample_id$category, decreasing = T),] # group sample id with age group for heatmap group
hmap_abundance <- taxa_corrected[damicro,c("species",cc$run_id)]
#hmap_abundance <- taxa_corrected[damicro,cc$run_id]
rownames(hmap_abundance) <- hmap_abundance$species
hmap_abundance <- hmap_abundance[,-1] 
hmap_abund_log2 <- log2(hmap_abundance+1)

# annoation col and row
sample_col <- as.data.frame(comsample_id[,c("category")])
row.names(sample_col) <- comsample_id$run_id
colnames(sample_col) <- "Age group"
#sample_col <- sample_col[,2]

# annoation row


dd <- taxa_corrected[which(taxa_corrected$species %in% rownames(hmap_abund_log2)),c("species","superkingdom")]
rownames(dd) <- dd$species
dd <- dd[rownames(hmap_abund_log2),]

micro_col <- as.data.frame(dd[,"superkingdom"])
row.names(micro_col) <- dd$species
colnames(micro_col) <- "Microbiome"

f_bar_hmap <- hmap[which(hmap$feature %in% damicro), c('feature','coef','group')]
hmap_list <- list(hmap_abund_log2,micro_col,sample_col,f_bar_hmap) # used to make bar plot 

save(hmap_list, file = "hmap_list.rda")
load("hmap_list.rda")
hmap_abund_log2 <-hmap_list[[1]]
micro_col <- hmap_list[[2]]
sample_col <- hmap_list[[3]]
pdf(paste("F7.3heatmap", "pdf", sep="."),width=12, height=14)
pheatmap(hmap_abund_log2, cluster_rows=FALSE, cluster_cols=F,annotation_col = sample_col,annotation_row = micro_col,
         gaps_col = c(41, 105), gaps_row = c(5,7,11),
         #col = brewer.pal(10, 'RdYlGn'),
         fontsize =18,
         show_colnames = F, show_rownames = T,
         annotation_names_row = F, 
         annotation_names_col = F)
dev.off()




f_bar_hmap <- hmap_list[[4]]
f_bar_hmap$group <- factor(f_bar_hmap$group, levels = c('Y_E','Y_C','E_C')) 
f_bar_hmap$feature <- factor(f_bar_hmap$feature, levels = rev(damicro))
pdf(paste("F7.4hmapgroup", "pdf", sep="."),width=12, height=14)
ggbarplot(f_bar_hmap,x="feature",y="coef",orientation="horiz",facet.by = 'group')+
  #theme_light()+
  theme(axis.text.y =element_text(size=24,face = "bold",color = "black"),axis.text.x=element_text(size=20,face="bold",color = "black"),axis.title =element_text(size=26,face="bold",color = "black"),
        strip.text = element_text(
          size = 15, color = "black"))+
  labs(x = "", y = expression(Log["2"](Fold-Change)))
dev.off()




#===============================================================
#
#   Build Bacterial co-concurrence networks to calculate degree of each node # SpiecEasi took too long time for build network
#   according to Difficulty in inferring microbial community structure based on co-occurrence network approaches
#   trandition method pearson exhit simmlar performance as spieceasi algorithm 
#===============================================================
load("kraken_taxa_corrected.rda")
load('comsample_id.rda')
rownames(taxa_corrected) <- taxa_corrected$taxonomy_id
id_sample <- comsample_id[which(comsample_id$category %in% "Y"),"run_id"]
taxa_abundance <- taxa_corrected[,id_sample]


# Number of expressed bacteria among samples based on CV<=2 and log2(counts+1)>1

# set CV threshold is 2
# calculate coefficient of variations
CV <- function(x){
  (sd(x)/mean(x))
}


sp_abundance<- log2(taxa_abundance+1)
dim(sp_abundance) # 194
# select bacteria based on CV
cc<-apply(sp_abundance,1,CV)
cc<-as.data.frame(cc)
cc$specis<-rownames(cc)
cc<-cc[which(cc$cc <= 2),]
spabund_cv<-sp_abundance[which(rownames(sp_abundance)%in% cc$specis),]



# 170
y_Net<- cor(t(spabund_cv), method=c("pearson"), use = "pairwise.complete.obs")
diag(y_Net) <- NA
y_Net3col<-as.data.frame(as.table(y_Net))
# summary(abs(y_Net3col$Freq))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   0.094   0.194   0.217   0.316   0.999    3519 


# remove weakly connections 
ii <- abs(y_Net) < 0.194 #  based on summary analysis for y_Net3col
y_Net[ii] <- NA

y_Net3col<-as.data.frame(as.table(y_Net))
length(unique(y_Net3col$Var1)) #3519
y_Net3col<-y_Net3col[-which(is.na(y_Net3col$Freq)),]
length(unique(y_Net3col$Var1)) #6988

save(y_Net3col,file = "y_Net3col.RData")

y_Netdegree <- apply(y_Net, 2, function(x) {sum(!is.na(x) & x!=0)})
y_Netdegree<- as.data.frame(y_Netdegree)
y_Netdegree$y_Netdegree<-as.numeric(y_Netdegree$y_Netdegree)
save(y_Netdegree, file = "y_Netdegree.RData")





library(igraph)
#  Nor network structure
set.seed(2023)
load("y_Net3col.RData")

g2 <- graph.data.frame(y_Net3col, directed=FALSE)

cw_hall <- cluster_walktrap(g2)
h_cw<- as.data.frame(sizes(cw_hall))  # 8 big clusters 
#save(cw_hall,file = "cw_hall.RData")
# ADD COLOR BASED ON THE CLUSTERS
components<- cw_hall$membership
colours = sample ( rainbow ( max ( components )
                             + 1) )
V (g2)$ color = colours [ components +1]

# Change Layout Structure in IGraph Plot based on Community
weights_H <- ifelse(crossing(cw_hall, g2), 1, 50)
coords<- layout_with_fr(g2, weights=weights_H)
png ('Fig8.1Ynetwork.png ', width = 1000 , height
     =900)
plot (g2 , layout =coords , vertex.label =NA ,edge.color="white",
      vertex.size =3)

dev.off()














