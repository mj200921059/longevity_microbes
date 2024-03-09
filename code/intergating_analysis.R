# Intergating analysis for longevity
# On DEC 14 2023
# By Jun
#====================================
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
