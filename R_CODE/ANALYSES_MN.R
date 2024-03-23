##############################################################################################################
#  Analyses for "Benefits of fencing farm dams to exclude livestock on seasonal total methane emissions and water quality" #
##############################################################################################################


rm(list = ls(all = TRUE))

#General libarary
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)
library(readxl)
library(writexl)
library(lubridate)
library(sf)
library(car)
library(GGally)
library(psych)
library(Hmisc)
library(PerformanceAnalytics)
library(corrplot)
library(lme4)
library(nlme)
library(MuMIn)
library(car)
library(effects)
library(tidyr)
library(lsmeans)
library(cowplot)
library(visreg)
library(ggsignif)
library(scales)
library(ggmap)
library(mgcv)
library(remef)
library(patchwork)


##############
# CHANGE REQUIRED
# Set where the R code is located
dir.script = "C:/Users/odebs/Desktop/GHG DATA AND ANALYSIS/BCL_Paper1_Analysis_data_code"
setwd(dir.script)
##############

#ANALYSIS TO CONDUCT
#(A) Organize the data
#(B) Fencing effects on greenhouse gas emissions across seasons in farm dams 
#(C) Fencing effects on nutrients and water quality
#(D) Effects of environmental factors on greenhouse gas emissions
#(E) Plot all results and data



#(A) Organise the data
#read the data---the csv fie
Data.imp = read.csv("GHG_Data.csv")

str(Data.imp)

if(FALSE){
  
  Data = Data.imp %>%
    mutate(SiteNameTreat = as.factor(SiteNameTreat),
           SiteNameRep = as.factor(SiteNameRep),
           SiteName = as.factor(SiteName), 
           Owner = as.factor(gsub('[[:digit:]]+', '', SiteName)), # Remove numbers (e.g., WebbWare1 becomes WebbWare)
           Site = as.factor(Site),
           Gas = as.factor(Gas),
           FileName = as.factor(FileName),
           Date_R = as.Date(Date_R, format = "%d/%m/%Y"),
           SlopeID = as.factor(SlopeID),
           Treatment = as.factor(Treatment),
           SplitID = as.factor(SplitID),
           Season = as.factor(Season),
           ppm_range = ppm_max - ppm_min
    )
  
  Data$Nitrogen[Data$Nitrogen == "<2.5"] = "2.5"
  Data$Phosphorus[Data$Phosphorus == "<0.10"] = "0.1"
  
  Data$Nitrogen = as.numeric(as.character(Data$Nitrogen))
  Data$Phosphorus = as.numeric(as.character(Data$Phosphorus))
  
  dim(Data)

### Separate Data averages or mean (Data_sum) for CH4 and CO2 
  Data_sum_CH4 = Data %>%
    filter(Gas == "CH4" & mgm2day_mean > 0) %>%
    group_by(SiteNameTreat, SiteName, Owner, Treatment, Season) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE) %>% 	
    rename(mgm2day_mean_CH4 = mgm2day_mean, mgm2day_sd_CH4 = mgm2day_sd, ppm_min_CH4 = ppm_min, ppm_max_CH4 = ppm_max, ppm_range_CH4 = ppm_range)
  
  
  Data_sum_CO2 = Data %>%
    filter(Gas == "CO2") %>%
    group_by(SiteNameTreat, SiteName, Owner, Treatment, Season) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE) %>% 	
    ungroup() %>%
    select(SiteNameTreat, mgm2day_mean, mgm2day_sd, ppm_min, ppm_max,ppm_range) %>%
    group_by(SiteNameTreat) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
    rename(mgm2day_mean_CO2 = mgm2day_mean, mgm2day_sd_CO2 = mgm2day_sd, ppm_min_CO2 = ppm_min, ppm_max_CO2 = ppm_max, ppm_range_CO2 = ppm_range)
  
  
  Data_sum = Data_sum_CH4 %>%
    merge(Data_sum_CO2, by = "SiteNameTreat", all = T) %>%
    mutate(mgm2day_mean_CO2eq = mgm2day_mean_CO2 + mgm2day_mean_CH4*28)
  
  dim(Data_sum)
  
  save(Data,Data_sum, file = "Data.RData")
}


#(B) Fencing effects on greenhouse gas emissions across seasons in farm dams 
##Begin the analysis for the gases and also add CO2eq
Data.CH4 = subset(Data, Gas == "CH4" & mgm2day_mean > 0)	
Data.CO2 = subset(Data, Gas == "CO2")

#EFFECTS OF TREATMENT AND SEASON ON GASES
#(1) Fixed-effect models for CH4
Data.CH4.WinterSummer = Data.CH4 %>% 
  filter(Owner == "WebbWare" | Season == "Summer")  %>% 
  filter(PondiRep > 10)

if(FALSE){
  
  Data.CH4.WinterSummer.lme1 = lme(data = Data.CH4.WinterSummer, 
                                   log10(mgm2day_mean) ~ Treatment*Season,
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varIdent(form = ~1|Treatment)
  )
  
  Data.CH4.WinterSummer.lme2= lme(data = Data.CH4.WinterSummer, 
                                  log10(mgm2day_mean) ~ Treatment*Season,
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  Data.CH4.WinterSummer.lme3= lme(data = Data.CH4.WinterSummer, 
                                  log10(mgm2day_mean) ~ Treatment*Season +
                                  scale(Surface_Ar) +scale(Temp_mean) +
                                  scale(Humidity_mean) + scale(log10(Rainfall_mean+1)),# +
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CH4.WinterSummer.lme4= lme(data = Data.CH4.WinterSummer, 
                                  log10(mgm2day_mean) ~ Treatment*Season +
                                    scale(log10(Rainfall_mean+1)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CH4.WinterSummer.lme5= lme(data = Data.CH4.WinterSummer, 
                                  log10(mgm2day_mean) ~ Treatment*Season +
                                    scale(log10(Rainfall_mean+1)) + 
                                    scale(Temp_mean),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CH4.WinterSummer.lme6= lme(data = Data.CH4.WinterSummer, 
                                  log10(mgm2day_mean) ~ Treatment*Season +
                                    scale(log10(Rainfall_mean+1))*scale(log10(Surface_Ar)) + 
                                    scale(Temp_mean)*scale(log10(Surface_Ar)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CH4.WinterSummer.lme7= lme(data = Data.CH4.WinterSummer, 
                                  log10(mgm2day_mean) ~ Treatment*Season +
                                    #scale(Temp_mean)*scale(log10(Surface_Ar)) +
                                    scale(log10(Rainfall_mean+1))*scale(log10(Surface_Ar)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CH4.WinterSummer.lme8= lme(data = Data.CH4.WinterSummer, 
                                  log10(mgm2day_mean) ~ Treatment*Season + scale(log10(Surface_Ar)) +
                                    #scale(Temp_mean)*scale(log10(Surface_Ar)) +
                                    scale(log10(Rainfall_mean+1)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CH4.WinterSummer.lme9= lme(data = Data.CH4.WinterSummer, 
                                  log10(mgm2day_mean) ~ Treatment*scale(Temp_mean) +
                                    scale(log10(Rainfall_mean+1)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CH4.WinterSummer.lme10= lme(data = Data.CH4.WinterSummer, 
                                   log10(mgm2day_mean) ~ Treatment + 
                                     scale(log10(Temp_mean)) +
                                     scale(log10(Rainfall_mean+1)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CH4.WinterSummer.lme11= lme(data = Data.CH4.WinterSummer, 
                                   log10(mgm2day_mean) ~ Treatment * 
                                     scale(log10(Temp_mean)) *
                                     scale(log10(Rainfall_mean+1)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CH4.WinterSummer.lme12= lme(data = Data.CH4.WinterSummer, 
                                   log10(mgm2day_mean) ~ Treatment*Season +
                                     scale(log10(Rainfall_mean+1)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Owner), varIdent(form = ~1|Season))
  )
  
  
  Data.CH4.WinterSummer.lme13 = lme(data = Data.CH4.WinterSummer, 
                                    log10(mgm2day_mean) ~ Treatment*Season +
                                      scale(log10(Temp_mean))*Treatment +
                                      scale(log10(Rainfall_mean+1))*Treatment,
                                    random = ~ 1|Owner/SiteNameTreat,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CH4.WinterSummer.lme14 = lme(data = Data.CH4.WinterSummer, 
                                    log10(mgm2day_mean) ~ Treatment*Season +
                                      scale(log10(Temp_mean))*Treatment +
                                      scale(log10(Rainfall_mean+1)),#*Treatment,
                                    random = ~ 1|Owner/SiteNameTreat,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CH4.WinterSummer.lme15 = lme(data = Data.CH4.WinterSummer, 
                                    log10(mgm2day_mean) ~ Treatment+Season +
                                      scale(log10(Temp_mean))*Treatment +
                                      scale(log10(Rainfall_mean+1)),#*Treatment,
                                    #		weight = varIdent(form = ~1|Treatment),
                                    #		weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season)),
                                    random = ~ 1+scale(log10(Temp_mean))|Owner/SiteNameTreat
  )
  
  
  Data.CH4.WinterSummer.lme16= lme(data = Data.CH4.WinterSummer, 
                                   log10(mgm2day_mean) ~ Treatment*Season +
                                     scale(log10(Rainfall_mean+1)) + 
                                     scale(log10(Surface_Ar)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  Data.CH4.WinterSummer.lme17= lme(data = Data.CH4.WinterSummer, 
                                   log10(mgm2day_mean) ~ Treatment*Season +
                                     (log10(Temp_mean))+
                                     (log10(Rainfall_mean+1)),
                                   # (log10(Surface_Ar)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  Data.CH4.WinterSummer.lme18= lme(data = Data.CH4.WinterSummer, 
                                   log10(mgm2day_mean) ~ Treatment*Season +
                                     scale(log10(Temp_mean)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  AICc(Data.CH4.WinterSummer.lme1,Data.CH4.WinterSummer.lme2,Data.CH4.WinterSummer.lme3,
       Data.CH4.WinterSummer.lme4,Data.CH4.WinterSummer.lme5,Data.CH4.WinterSummer.lme6,Data.CH4.WinterSummer.lme7,
       Data.CH4.WinterSummer.lme8,Data.CH4.WinterSummer.lme9,Data.CH4.WinterSummer.lme10,Data.CH4.WinterSummer.lme11,
       Data.CH4.WinterSummer.lme12,Data.CH4.WinterSummer.lme13,Data.CH4.WinterSummer.lme14,Data.CH4.WinterSummer.lme15,
       Data.CH4.WinterSummer.lme16,Data.CH4.WinterSummer.lme17, Data.CH4.WinterSummer.lme18)
  Anova(Data.CH4.WinterSummer.lme4)
  plot(allEffects(Data.CH4.WinterSummer.lme4))
  Data.CH4.WinterSummer.Best = Data.CH4.WinterSummer.lme4
  #lsmeans(object=Data.CH4.WinterSummer.Best, pairwise ~ Treatment, adjust= "tukey")
  save(Data.CH4.WinterSummer.Best, file = "RESULTS/Data.CH4.WinterSummer.Best.RData")
  

} 

# Check assumptions
if(FALSE){
  
  plot(allEffects(Data.CH4.WinterSummer.Best))
  
  # Check assumptions
  plot(Data.CH4.WinterSummer.Best)
  scatterplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ scale(log10(Data.CH4.WinterSummer$Temp_mean)))
  scatterplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ scale(Data.CH4.WinterSummer$Humidity_mean))
  scatterplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ scale(log10(Data.CH4.WinterSummer$Rainfall_mean+1)))
  boxplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ Data.CH4.WinterSummer$Treatment)
  boxplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ Data.CH4.WinterSummer$SiteNameTreat)
  boxplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ Data.CH4.WinterSummer$SiteName)
  boxplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ Data.CH4.WinterSummer$Owner)
  boxplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ Data.CH4.WinterSummer$Season)
  
  scatterplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ scale(Data.CH4.WinterSummer$Surface_Ar))
  scatterplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ scale(Data.CH4.WinterSummer$Vegetation.20m.Area.m2.))
  scatterplot(resid(Data.CH4.WinterSummer.Best, type = "pearson") ~ scale(Data.CH4.WinterSummer$Depth))
}

# Plots
if(FALSE){
  
  Data.CH4.WinterSummer.Best.df = as.data.frame(allEffects(Data.CH4.WinterSummer.Best))$`Treatment:Season`
  
  # Convert from mg/m2/day to kg.ha.year
  Data.CH4.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")] = 365.25*10000*(10^Data.CH4.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")])/1e+6
  
  Data.CH4.WinterSummer.Best.df$Treatment = factor(Data.CH4.WinterSummer.Best.df$Treatment, levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))
  
  ggplot(data = Data.CH4.WinterSummer.Best.df, aes(x = Season, y = fit)) +
    geom_bar(aes(fill = Treatment), position=position_dodge(1), stat="identity") +
    geom_linerange(aes(col = Treatment ,ymin = lower, ymax = upper), position=position_dodge(1)) +
    theme_bw() +
    labs(y = "Methane emissions from farm dams (kg.ha.year)") +
    annotate(geom = "text", x = c(1,2), y = c(500, 500), label = c("-72%", "-92%"), size = 12)
  
  ggsave(width = 11.666666, height = 5.194444, filename = "RESULTS/Data.CH4.WinterSummer.Best.plot.pdf", device = "pdf")
  ggsave(width = 11.666666, height = 5.194444, filename = "RESULTS/Data.CH4.WinterSummer.Best.plot.jpeg", device = "jpeg")
  
}


#(2) Fixed-effect models for CO2
Data.CO2.WinterSummer = Data.CO2 %>% 
  filter(Owner == "WebbWare" | Season == "Summer")  %>% 
  filter(PondiRep > 10) %>%
  dplyr::mutate(log.mgm2day_mean = log10(mgm2day_mean + 1170))#1170 to remove negatives


if(FALSE){
  
  Data.CO2.WinterSummer.lme1 = lme(data = Data.CO2.WinterSummer, 
                                   log.mgm2day_mean ~ Treatment*Season,
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varIdent(form = ~1|Treatment)
  )
  
  Data.CO2.WinterSummer.lme2= lme(data = Data.CO2.WinterSummer, 
                                  log.mgm2day_mean ~ Treatment*Season,
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  Data.CO2.WinterSummer.lme3= lme(data = Data.CO2.WinterSummer, 
                                  log.mgm2day_mean ~ Treatment*Season +
                                  scale(Surface_Ar) + scale(Temp_mean) +
                                  scale(Humidity_mean) + scale(log10(Rainfall_mean+1)),# +
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2.WinterSummer.lme4= lme(data = Data.CO2.WinterSummer, 
                                  log.mgm2day_mean ~ Treatment*Season +
                                    scale(log10(Rainfall_mean+1)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2.WinterSummer.lme5= lme(data = Data.CO2.WinterSummer, 
                                  log.mgm2day_mean ~ Treatment*Season +
                                    scale(log10(Rainfall_mean+1)) + 
                                    scale(Temp_mean),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2.WinterSummer.lme6= lme(data = Data.CO2.WinterSummer, 
                                  log.mgm2day_mean ~ Treatment*Season +
                                    scale(log10(Rainfall_mean+1))*scale(log10(Surface_Ar)) + 
                                    scale(Temp_mean)*scale(log10(Surface_Ar)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2.WinterSummer.lme7= lme(data = Data.CO2.WinterSummer, 
                                  log.mgm2day_mean ~ Treatment*Season +
                                    #scale(Temp_mean)*scale(log10(Surface_Ar)) +
                                    scale(log10(Rainfall_mean+1))*scale(log10(Surface_Ar)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2.WinterSummer.lme8= lme(data = Data.CO2.WinterSummer, 
                                  log.mgm2day_mean ~ Treatment*Season + scale(log10(Surface_Ar)) +
                                    #scale(Temp_mean)*scale(log10(Surface_Ar)) +
                                    scale(log10(Rainfall_mean+1)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2.WinterSummer.lme9= lme(data = Data.CO2.WinterSummer, 
                                  log.mgm2day_mean ~ Treatment*scale(Temp_mean) +
                                    scale(log10(Rainfall_mean+1)),
                                  random = ~ 1|Owner/SiteNameTreat,
                                  weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CO2.WinterSummer.lme10= lme(data = Data.CO2.WinterSummer, 
                                   log.mgm2day_mean ~ Treatment + 
                                     scale(log10(Temp_mean)) +
                                     scale(log10(Rainfall_mean+1)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CO2.WinterSummer.lme11= lme(data = Data.CO2.WinterSummer, 
                                   log.mgm2day_mean ~ Treatment * 
                                     scale(log10(Temp_mean)) *
                                     scale(log10(Rainfall_mean+1)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CO2.WinterSummer.lme12= lme(data = Data.CO2.WinterSummer, 
                                   log.mgm2day_mean ~ Treatment*Season +
                                     scale(log10(Rainfall_mean+1)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Owner), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2.WinterSummer.lme13 = lme(data = Data.CO2.WinterSummer, 
                                    log.mgm2day_mean ~ Treatment*Season +
                                      scale(log10(Temp_mean))*Treatment +
                                      scale(log10(Rainfall_mean+1))*Treatment,
                                    random = ~ 1|Owner/SiteNameTreat,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2.WinterSummer.lme14 = lme(data = Data.CO2.WinterSummer, 
                                    log.mgm2day_mean ~ Treatment*Season +
                                      scale(log10(Temp_mean))*Treatment +
                                      scale(log10(Rainfall_mean+1)),#*Treatment,
                                    random = ~ 1|Owner/SiteNameTreat,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2.WinterSummer.lme15= lme(data = Data.CO2.WinterSummer, 
                                   log.mgm2day_mean ~ Treatment*Season +
                                     scale(log10(Rainfall_mean+1)) + 
                                     scale(log10(Surface_Ar)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  Data.CO2.WinterSummer.lme16= lme(data = Data.CO2.WinterSummer, 
                                   log.mgm2day_mean ~ Treatment*Season +
                                     (log10(Temp_mean))+
                                     (log10(Rainfall_mean+1)),
                                   # (log10(Surface_Ar)),
                                   random = ~ 1|Owner/SiteNameTreat,
                                   weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  AICc(Data.CO2.WinterSummer.lme1,Data.CO2.WinterSummer.lme2,Data.CO2.WinterSummer.lme3,
       Data.CO2.WinterSummer.lme4,Data.CO2.WinterSummer.lme5,Data.CO2.WinterSummer.lme6,
       Data.CO2.WinterSummer.lme7,Data.CO2.WinterSummer.lme8,Data.CO2.WinterSummer.lme9,
       Data.CO2.WinterSummer.lme10,Data.CO2.WinterSummer.lme11,Data.CO2.WinterSummer.lme12,
       Data.CO2.WinterSummer.lme13,Data.CO2.WinterSummer.lme14,Data.CO2.WinterSummer.lme15,Data.CO2.WinterSummer.lme16)
  
  Data.CO2.WinterSummer.Best = Data.CO2.WinterSummer.lme12
  Anova(Data.CO2.WinterSummer.Best)
  plot(allEffects(Data.CO2.WinterSummer.Best))
  save(Data.CO2.WinterSummer.Best, file = "RESULTS/Data.CO2.WinterSummer.Best.RData")
  
}

# Check assumptions
if(FALSE){
  
  plot(allEffects(Data.CO2.WinterSummer.Best))
  
  # Check assumptions
  plot(Data.CO2.WinterSummer.Best)
  scatterplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ scale(log10(Data.CO2.WinterSummer$Temp_mean)))
  scatterplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ scale(Data.CO2.WinterSummer$Humidity_mean))
  scatterplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ scale(log10(Data.CO2.WinterSummer$Rainfall_mean+1)))
  boxplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ Data.CO2.WinterSummer$Treatment)
  boxplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ Data.CO2.WinterSummer$SiteNameTreat)
  boxplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ Data.CO2.WinterSummer$SiteName)
  boxplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ Data.CO2.WinterSummer$Owner)
  boxplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ Data.CO2.WinterSummer$Season)
  
  scatterplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ scale(Data.CO2.WinterSummer$Surface_Ar))
  scatterplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ scale(Data.CO2.WinterSummer$Vegetation.20m.Area.m2.))
  scatterplot(resid(Data.CO2.WinterSummer.Best, type = "pearson") ~ scale(Data.CO2.WinterSummer$Depth))
}


# Plots
if(FALSE){
  
  Data.CO2.WinterSummer.Best.df = as.data.frame(allEffects(Data.CO2.WinterSummer.Best))$`Treatment:Season`
  
  # Convert from mg/m2/day to kg.ha.year
  Data.CO2.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")] = 365.25*10000*(10^Data.CO2.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")])/1e+6
  
  Data.CO2.WinterSummer.Best.df$Treatment = factor(Data.CO2.WinterSummer.Best.df$Treatment, levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))
  
  ggplot(data = Data.CO2.WinterSummer.Best.df, aes(x = Season, y = fit)) +
    geom_bar(aes(fill = Treatment), position=position_dodge(1), stat="identity") +
    geom_linerange(aes(col = Treatment ,ymin = lower, ymax = upper), position=position_dodge(1)) +
    theme_bw() +
    labs(y = "Carbon dioxide emissions from farm dams (kg.ha.year)") +
    annotate(geom = "text", x = c(1,2), y = c(5500, 5500), label = c("4.32%", "1.55%"), size = 12)
  
  ggsave(width = 11.666666, height = 5.194444, filename = "RESULTS/Data.CO2.WinterSummer.Best.plot.pdf", device = "pdf")
  ggsave(width = 11.666666, height = 5.194444, filename = "RESULTS/Data.CO2.WinterSummer.Best.plot.jpeg", device = "jpeg")
  
}

#(3) Fixed-effect models for CO2eq
Data.CO2eq.WinterSummer = Data_sum %>% 
  filter(Owner == "WebbWare" | Season == "Summer")  %>% 
  mutate(log.mgm2day_mean_CO2eq = log10(mgm2day_mean_CO2eq))

hist(Data.CO2eq.WinterSummer$log.mgm2day_mean_CO2eq)

ggplot(data = Data.CO2eq.WinterSummer, aes(x = log.mgm2day_mean_CO2eq)) +
  geom_histogram() + facet_grid(Treatment~.)


if(FALSE){
  
  Data.CO2eq.WinterSummer.lme1 = lme(data = Data.CO2eq.WinterSummer, 
                                     log.mgm2day_mean_CO2eq ~ Treatment*Season,
                                     random = ~ 1|Owner,
                                     weight = varIdent(form = ~1|Treatment)
  )
  
  Data.CO2eq.WinterSummer.lme2= lme(data = Data.CO2eq.WinterSummer, 
                                    log.mgm2day_mean_CO2eq ~ Treatment*Season,
                                    random = ~ 1|Owner,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  Data.CO2eq.WinterSummer.lme3= lme(data = Data.CO2eq.WinterSummer, 
                                    log.mgm2day_mean_CO2eq ~ Treatment*Season +
                                    scale(Surface_Ar) + scale(Temp_mean) +
                                    scale(Humidity_mean) + scale(log10(Rainfall_mean+1)),# +
                                    random = ~ 1|Owner,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme4= lme(data = Data.CO2eq.WinterSummer, 
                                    log.mgm2day_mean_CO2eq ~ Treatment*Season +
                                      scale(log10(Rainfall_mean+1)),
                                    random = ~ 1|Owner,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme5= lme(data = Data.CO2eq.WinterSummer, 
                                    log.mgm2day_mean_CO2eq ~ Treatment*Season +
                                      scale(log10(Rainfall_mean+1)) + 
                                      scale(Temp_mean),
                                    random = ~ 1|Owner,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme6= lme(data = Data.CO2eq.WinterSummer, 
                                    log.mgm2day_mean_CO2eq ~ Treatment*Season +
                                      scale(log10(Rainfall_mean+1))*scale(log10(Surface_Ar)) + 
                                      scale(Temp_mean)*scale(log10(Surface_Ar)),
                                    random = ~ 1|Owner,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme7= lme(data = Data.CO2eq.WinterSummer, 
                                    log.mgm2day_mean_CO2eq ~ Treatment*Season +
                                      #scale(Temp_mean)*scale(log10(Surface_Ar)) +
                                      scale(log10(Rainfall_mean+1))*scale(log10(Surface_Ar)),
                                    random = ~ 1|Owner,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme8= lme(data = Data.CO2eq.WinterSummer, 
                                    log.mgm2day_mean_CO2eq ~ Treatment*Season + scale(log10(Surface_Ar)) +
                                      #scale(Temp_mean)*scale(log10(Surface_Ar)) +
                                      scale(log10(Rainfall_mean+1)),
                                    random = ~ 1|Owner,
                                    weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme9= lme(data = Data.CO2eq.WinterSummer, 
                                    log.mgm2day_mean_CO2eq ~ Treatment*scale(Temp_mean) +
                                      scale(log10(Rainfall_mean+1)),
                                    random = ~ 1|Owner,
                                    weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CO2eq.WinterSummer.lme10= lme(data = Data.CO2eq.WinterSummer, 
                                     log.mgm2day_mean_CO2eq ~ Treatment*Season + 
                                       scale(log10(Temp_mean)) +
                                       scale(log10(Rainfall_mean+1)),
                                     random = ~ 1|Owner,
                                     weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CO2eq.WinterSummer.lme11= lme(data = Data.CO2eq.WinterSummer, 
                                     log.mgm2day_mean_CO2eq ~ Treatment * 
                                       scale(log10(Temp_mean)) *
                                       scale(log10(Rainfall_mean+1)),
                                     random = ~ 1|Owner,
                                     weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CO2eq.WinterSummer.lme12= lme(data = Data.CO2eq.WinterSummer, 
                                     log.mgm2day_mean_CO2eq ~ Treatment*Season +
                                       scale(log10(Rainfall_mean+1)),
                                     random = ~ 1|Owner
  )
  
  
  Data.CO2eq.WinterSummer.lme13 = lme(data = Data.CO2eq.WinterSummer, 
                                      log.mgm2day_mean_CO2eq ~ Treatment*Season +
                                        scale(log10(Temp_mean))*Treatment +
                                        scale(log10(Rainfall_mean+1))*Treatment,
                                      random = ~ 1|Owner,
                                      weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme14 = lme(data = Data.CO2eq.WinterSummer, 
                                      log.mgm2day_mean_CO2eq ~ Treatment*Season +
                                        scale(log10(Temp_mean))*Treatment +
                                        scale(log10(Rainfall_mean+1)),#*Treatment,
                                      random = ~ 1|Owner,
                                      weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme15= lme(data = Data.CO2eq.WinterSummer, 
                                     log.mgm2day_mean_CO2eq ~ Treatment*Season +
                                       scale(log10(Rainfall_mean+1)) + 
                                       scale(log10(Surface_Ar)),
                                     random = ~ 1|Owner,
                                     weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme16= lme(data = Data.CO2eq.WinterSummer, 
                                     log.mgm2day_mean_CO2eq ~ Treatment + 
                                       scale(log10(Temp_mean)) *
                                       scale(log10(Rainfall_mean+1)),
                                     random = ~ 1|Owner,
                                     weight = varIdent(form = ~1|Treatment)
  )
  
  
  
  Data.CO2eq.WinterSummer.lme17= lme(data = Data.CO2eq.WinterSummer, 
                                     log.mgm2day_mean_CO2eq ~ Treatment*Season + 
                                       scale(log10(Temp_mean)) +
                                       scale(log10(Rainfall_mean+1)) +
                                       scale(log10(Surface_Ar)),
                                     random = ~ 1|Owner,
                                     weight = varIdent(form = ~1|Treatment)
  )
  
  
  Data.CO2eq.WinterSummer.lme18= lme(data = Data.CO2eq.WinterSummer, 
                                     log.mgm2day_mean_CO2eq ~ Treatment*Season + 
                                       scale(log10(Temp_mean)) +
                                       scale(log10(Rainfall_mean+1)) +
                                       scale(log10(Surface_Ar)),
                                     random = ~ 1|Owner,
                                     weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  Data.CO2eq.WinterSummer.lme19= lme(data = Data.CO2eq.WinterSummer, 
                                     log.mgm2day_mean_CO2eq ~ Treatment + 
                                       scale(log10(Temp_mean)) +
                                       scale(log10(Rainfall_mean+1)) +
                                       scale(log10(Surface_Ar)),
                                     random = ~ 1|Owner,
                                     weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  
  AICc(Data.CO2eq.WinterSummer.lme1,Data.CO2eq.WinterSummer.lme2,Data.CO2eq.WinterSummer.lme3,
       Data.CO2eq.WinterSummer.lme4,Data.CO2eq.WinterSummer.lme5,Data.CO2eq.WinterSummer.lme6,
       Data.CO2eq.WinterSummer.lme7,Data.CO2eq.WinterSummer.lme8,Data.CO2eq.WinterSummer.lme9,
       Data.CO2eq.WinterSummer.lme10,Data.CO2eq.WinterSummer.lme11,Data.CO2eq.WinterSummer.lme12,
       Data.CO2eq.WinterSummer.lme13,Data.CO2eq.WinterSummer.lme14,Data.CO2eq.WinterSummer.lme15,
       Data.CO2eq.WinterSummer.lme16,Data.CO2eq.WinterSummer.lme17,Data.CO2eq.WinterSummer.lme18,
       Data.CO2eq.WinterSummer.lme19)
  
  Data.CO2eq.WinterSummer.Best = Data.CO2eq.WinterSummer.lme1 
  Anova(Data.CO2eq.WinterSummer.lme1)
  plot(allEffects(Data.CO2eq.WinterSummer.Best))
  save(Data.CO2eq.WinterSummer.Best, file = "RESULTS/Data.CO2eq.WinterSummer.Best.RData")
  
  
}


# Check assumptions
if(FALSE){
  
  plot(allEffects(Data.CO2eq.WinterSummer.Best))
  
  # Check assumptions
  plot(Data.CO2eq.WinterSummer.Best)
  scatterplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ scale(log10(Data.CO2eq.WinterSummer$Temp_mean)))
  scatterplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ scale(Data.CO2eq.WinterSummer$Humidity_mean))
  scatterplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ scale(log10(Data.CO2eq.WinterSummer$Rainfall_mean+1)))
  boxplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Treatment)
  boxplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$SiteNameTreat)
  boxplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$SiteName)
  boxplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Owner)
  boxplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Season)
  
  scatterplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ scale(Data.CO2eq.WinterSummer$Surface_Ar))
  scatterplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ scale(Data.CO2eq.WinterSummer$Vegetation.20m.Area.m2.))
  scatterplot(resid(Data.CO2eq.WinterSummer.Best, type = "pearson") ~ scale(Data.CO2eq.WinterSummer$Depth))
}

# Plots
if(FALSE){
  
  Data.CO2eq.WinterSummer.Best.df = as.data.frame(allEffects(Data.CO2eq.WinterSummer.Best))$`Treatment:Season`
  
  # Convert from mg/m2/day to kg.ha.year
  Data.CO2eq.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")] = 365.25*10000*(10^Data.CO2eq.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")])/1e+6
  
  Data.CO2eq.WinterSummer.Best.df$Treatment = factor(Data.CO2eq.WinterSummer.Best.df$Treatment, levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))
  
  ggplot(data = Data.CO2eq.WinterSummer.Best.df, aes(x = Season, y = fit)) +
    geom_bar(aes(fill = Treatment), position=position_dodge(1), stat="identity") +
    geom_linerange(aes(col = Treatment ,ymin = lower, ymax = upper), position=position_dodge(1)) +
    theme_bw() +
    labs(y = "CO2eq emissions from farm dams (kg.ha.year)") +
    annotate(geom = "text", x = c(1,2), y = c(14000, 14000), label = c("-59.22%", "-73.18%"), size = 12)
  
  ggsave(width = 11.666666, height = 5.194444, filename = "RESULTS/Data.CO2eq.WinterSummer.Best.plot.pdf", device = "pdf")
  ggsave(width = 11.666666, height = 5.194444, filename = "RESULTS/Data.CO2eq.WinterSummer.Best.plot.jpeg", device = "jpeg")
  
}


#(C) Fencing effects on nutrients and water quality
##We used total GHG data(CO2-eq)
## Effects of treatment and season on Dissolved Oxygen##

if(FALSE){
  
  Data.DO.WinterSummer.lme1 = lme(data = Data.CO2eq.WinterSummer, 
                                  DO ~ Treatment*Season,
                                  random = ~ 1|Owner,
                                  weight = varIdent(form = ~1|Treatment)
                                  
  )
  
  Data.DO.WinterSummer.lme2 = lme(data = Data.CO2eq.WinterSummer, 
                                  DO ~ Treatment*Season,
                                  random = ~ 1|Owner,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season))
  )
  
  AICc(Data.DO.WinterSummer.lme1,Data.DO.WinterSummer.lme2)
  Data.DO.WinterSummer.Best = Data.DO.WinterSummer.lme2
  Anova(Data.DO.WinterSummer.Best)
  plot(allEffects(Data.DO.WinterSummer.Best))
  save(Data.DO.WinterSummer.Best, file = "RESULTS/Data.DO.WinterSummer.Best.RData")
  
}

# Check assumptions
if(FALSE){
  scatterplot(residuals(Data.DO.WinterSummer.Best, type = "pearson")~predict(Data.DO.WinterSummer.Best))
  scatterplot(residuals(Data.DO.WinterSummer.Best, type = "pearson")~scale(log10(Data.CO2eq.WinterSummer$DO)))
  boxplot(resid(Data.DO.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Treatment)
  boxplot(resid(Data.DO.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Season)
}


##Effect of treatment and season on pH##
if(FALSE){
  
  Data.pH.WinterSummer.lme1 = lme(data = Data.CO2eq.WinterSummer, 
                                  pH ~ Treatment*Season,
                                  random = ~ 1|Owner,
                                  weight = varIdent(form = ~1|Treatment),
                                  na.action = na.exclude
  )
  
  Data.pH.WinterSummer.lme2 = lme(data = Data.CO2eq.WinterSummer, 
                                  pH ~ Treatment*Season,
                                  random = ~ 1|Owner,
                                  weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season)),
                                  na.action = na.exclude
  )
  
  AICc(Data.pH.WinterSummer.lme1,Data.pH.WinterSummer.lme2)
  Data.pH.WinterSummer.Best = Data.pH.WinterSummer.lme2
  Anova(Data.pH.WinterSummer.Best)
  plot(allEffects(Data.pH.WinterSummer.Best))
  save(Data.pH.WinterSummer.Best, file = "RESULTS/Data.pH.WinterSummer.Best.RData")
  
}

# Check assumptions
if(FALSE){
  scatterplot(residuals(Data.pH.WinterSummer.Best, type = "pearson")~predict(Data.pH.WinterSummer.Best))
  scatterplot(residuals(Data.pH.WinterSummer.Best, type = "pearson")~scale(log10(Data.CO2eq.WinterSummer$pH)))
  boxplot(resid(Data.pH.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Treatment)
  boxplot(resid(Data.pH.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Season)
}


##Effect of treatment on Water Temperature##
if(FALSE){
  Data.WaterT.WinterSummer.lme1 = lme(data = Data.CO2eq.WinterSummer, 
                                      Water_Temp ~ Treatment*Season,
                                      random = ~ 1|Owner,
                                      weight = varIdent(form = ~1|Treatment),
                                      na.action = na.exclude
  )
  
  Data.WaterT.WinterSummer.lme2 = lme(data = Data.CO2eq.WinterSummer, 
                                      Water_Temp ~ Treatment*Season,
                                      random = ~ 1|Owner,
                                      weight = varComb(varIdent(form = ~1|Treatment), varIdent(form = ~1|Season)),
                                      na.action = na.exclude
  )
  
  AICc(Data.WaterT.WinterSummer.lme1,Data.WaterT.WinterSummer.lme2)
  Data.WaterT.WinterSummer.Best = Data.WaterT.WinterSummer.lme2
  Anova(Data.WaterT.WinterSummer.Best)
  plot(allEffects(Data.WaterT.WinterSummer.Best))
  save(Data.WaterT.WinterSummer.Best, file = "RESULTS/Data.WaterT.WinterSummer.Best.RData")
  
}

# Check assumptions
if(FALSE){
  scatterplot(residuals(Data.WaterT.WinterSummer.Best, type = "pearson")~predict(Data.WaterT.WinterSummer.Best))
  scatterplot(residuals(Data.WaterT.WinterSummer.Best, type = "pearson")~scale(log10(Data.CO2eq.WinterSummer$Water_Temp)))
  boxplot(resid(Data.WaterT.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Treatment)
  boxplot(resid(Data.WaterT.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Season)
}


##Effects of treatment on Nitrogen/Note we only analysed for winter##
if(FALSE){
  
  Data.N.WinterSummer.lme1 = lme(data = Data.CO2eq.WinterSummer[!is.na(Data.CO2eq.WinterSummer$Nitrogen),], 
                                 log10(Nitrogen) ~ Treatment,
                                 random = ~ 1|Owner,
                                 weight = varIdent(form = ~1|Treatment)
  )
  
  Data.N.WinterSummer.lme2 = lme(data = Data.CO2eq.WinterSummer[!is.na(Data.CO2eq.WinterSummer$Nitrogen),], 
                                 log10(Nitrogen) ~ Treatment,
                                 random = ~ 1|Owner
  )
  
  AICc(Data.N.WinterSummer.lme1,Data.N.WinterSummer.lme2)
  Data.N.WinterSummer.Best = Data.N.WinterSummer.lme2
  Anova(Data.N.WinterSummer.Best)
  plot(allEffects(Data.N.WinterSummer.Best))
  save(Data.N.WinterSummer.Best, file = "RESULTS/Data.N.WinterSummer.Best.RData")
  
}

# Check assumptions
if(FALSE){
  scatterplot(residuals(Data.N.WinterSummer.Best, type = "pearson")~predict(Data.N.WinterSummer.Best))
  #scatterplot(residuals(Data.N.WinterSummer.Best, type = "pearson")~Data.CO2eq.WinterSummer$Nitrogen) #remove NAs before plot
  #boxplot(resid(Data.N.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Treatment) #remove NAs before plot
}


##Effects of treatment on Phosphorus/Note we only analysed for winter##
if(FALSE){
  
  Data.P.WinterSummer.lme1 = lme(data = Data.CO2eq.WinterSummer[!is.na(Data.CO2eq.WinterSummer$Phosphorus),],
                                 scale(log10(Phosphorus)) ~ Treatment,
                                 random = ~ 1|Owner,
                                 #weight = varIdent(form = ~1|Treatment)
  )
  
  Data.P.WinterSummer.lme2 = lme(data = Data.CO2eq.WinterSummer[!is.na(Data.CO2eq.WinterSummer$Phosphorus),], 
                                 log10(Phosphorus) ~ Treatment,
                                 random = ~ 1|Owner
  )
  
  AICc(Data.P.WinterSummer.lme1,Data.P.WinterSummer.lme2)
  Data.P.WinterSummer.Best = Data.P.WinterSummer.lme2
  Anova(Data.P.WinterSummer.Best)
  plot(allEffects(Data.P.WinterSummer.Best))
  save(Data.P.WinterSummer.Best, file = "RESULTS/Data.P.WinterSummer.Best.RData")
  
}

# Check assumptions
if(FALSE){
  scatterplot(residuals(Data.P.WinterSummer.Best, type = "pearson")~predict(Data.P.WinterSummer.Best))
  #scatterplot(residuals(Data.P.WinterSummer.Best, type = "pearson")~Data.CO2eq.WinterSummer$Phosphorus #remove NAs before plot
  #boxplot(resid(Data.P.WinterSummer.Best, type = "pearson") ~ Data.CO2eq.WinterSummer$Treatment) #remove NAs before plot
}



#(D) Effects of environmental factors on greenhouse gas emissions
#We used the data averages (Data_sum)

# Check on the effects of factors on methane
Data_sum.CH4 = Data_sum %>% 
  filter(Owner == "WebbWare" | Season == "Summer") 
  
if(FALSE){
  
  Data.CH4.Temp.1 = lme(data = Data_sum.CH4, 
                        log10(mgm2day_mean_CH4) ~ scale(log10(Temp_mean)) + 
                          scale(log10(Rainfall_mean+1)) +
                          scale(log10(Humidity_mean)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
  )
  
  Data.CH4.Temp.2 = lme(data = Data_sum.CH4, 
                        log10(mgm2day_mean_CH4) ~ scale(log10(Temp_mean)) + 
                          scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
  )
  
  Data.CH4.Temp.3 = lme(data = Data_sum.CH4, 
                        log10(mgm2day_mean_CH4) ~ scale(log10(Temp_mean)),# + 
                        #scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
                        na.action = na.exclude,
  )
  
  Data.CH4.Temp.4 = lme(data = Data_sum.CH4, 
                        log10(mgm2day_mean_CH4) ~ scale(log10(Temp_mean)),# +  
                        #scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
                        #method = "ML",
                        na.action = na.exclude,
  )
  
  Data.CH4.Temp.5 = lme(data = Data_sum.CH4, 
                        log10(mgm2day_mean_CH4) ~ scale(Temp_mean),# + 
                        #scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        #weight = varExp(form = ~Temp_mean),
                        random = ~ 1|Owner,
  )
  
  Data.CH4.Temp.6 = lme(data = Data_sum.CH4, 
                        log10(mgm2day_mean_CH4) ~ scale(Temp_mean),# + 
                        #scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        weight = varPower(form = ~Temp_mean),
                        random = ~ 1|Owner,
  )
  
  
  Data.CH4.Temp.7 = lme(data = Data_sum.CH4, 
                        log10(mgm2day_mean_CH4) ~ scale(Temp_mean),# + 
                        #scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        weight = varConstProp(form = ~Temp_mean),
                        random = ~ 1|Owner,
  )
  
  
  Data.CH4.Temp.8 = lme(data = Data_sum.CH4, 
                        log10(mgm2day_mean_CH4) ~ scale(log10(Water_Temp))+  
                          scale(log10(Rainfall_mean+1))* (Treatment),# +
                        #scale(log10(Humidity_mean)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
                        #method = "ML",
                        na.action = na.exclude,
  )
  
  
  Data.CH4.Temp.9 = lme(data = Data_sum.CH4, 
                        log10(mgm2day_mean_CH4) ~ scale(log10(Water_Temp)) +
                          scale(log10(Temp_mean)),
                        #scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
                        #method = "ML",
                        na.action = na.exclude,
  )
  
  
  Data.CH4.Temp.10 = lme(data = Data_sum.CH4, 
                         log10(mgm2day_mean_CH4) ~ scale(log10(Water_Temp))+
                           scale(log10(Temp_mean))+
                           scale(log10(Rainfall_mean+1))+# +
                           scale(log10(Humidity_mean))+
                           scale(log10(Depth))+
                           scale(log10(Surface_Ar)),
                         #scale(log10(Vegetation.20m.Area.m2.)),
                         #weight = varIdent(form = ~1|Treatment),
                         random = ~ 1|Owner,
                         #method = "ML",
                         na.action = na.exclude,
                         
  )
  
  AICc(Data.CH4.Temp.1,Data.CH4.Temp.2,Data.CH4.Temp.3,Data.CH4.Temp.4,Data.CH4.Temp.5,
       Data.CH4.Temp.6,Data.CH4.Temp.7, Data.CH4.Temp.8, Data.CH4.Temp.9, Data.CH4.Temp.10)
  Data.CH4.Temp.Best = Data.CH4.Temp.9
  Anova(Data.CH4.Temp.Best)
  plot(allEffects(Data.CH4.Temp.Best))
  #vif(Data.CH4.Temp.Best)
  save(Data.CH4.Temp.Best, file = "RESULTS/Data.CH4.Temp.Best.RData")
  
}

# MODEL VALIDATION
if(FALSE){
  Data_sum.CH4$Resid = residuals(Data.CH4.Temp.Best, type = "pearson")
  Data_sum.CH4$Fitted = predict(Data.CH4.Temp.Best)
  
  scatterplot(data = Data_sum.CH4, Resid~Fitted)
  
  # Verify the final model
  scatterplot(data = Data_sum.CH4, Resid~Temp_mean)
  boxplot(data = Data_sum.CH4, Resid~Owner)
  
}

##FIT THE BEST MODEL FOR THE OVERALL PLOT LATER
# First, remove NA in Stock from the dataset because makes issues below
Data_sum.CH4_noNA = Data_sum.CH4[!is.na(Data_sum.CH4$Water_Temp),]


# Fit an equivalent model with all raw values
Data.CH4.Temp.Best.plot = lme(data = Data_sum.CH4_noNA, 
                              mgm2day_mean_CH4 ~ Water_Temp +
                                Temp_mean,
                              random = ~ 1|Owner,
                              
)

Anova(Data.CH4.Temp.Best.plot)


# 1) Plot Temp_mean EFFECT CH4
Data.CH4.Temp.Best.res = Effect("Temp_mean", partial.residuals= T, Data.CH4.Temp.Best.plot)
Temp.plotCH4 = plot(Data.CH4.Temp.Best.res, smooth.residuals=F, colors = c('black', 'black'),
                    xlab="Air Temperature (°C)", ylab="CH4 Emissions (kg ha^-1 year^-1)",
                    main='')
print(Temp.plotCH4)
dev.copy(pdf, "RESULTS/Temp.plotCH4.pdf", width = 6.24, height = 5.06)
dev.off()
dev.copy(jpeg, "RESULTS/Temp.plotCH4.jpeg", width = 6.24 * 300, height = 5.06 * 300, units = "px", res = 300)
dev.off()


# 2) Plot Water_temp EFFECT CH4
Data.CH4.Temp.Best.res = Effect("Water_Temp", partial.residuals= T, Data.CH4.Temp.Best.plot)
Water_Temp.plotCH4 = plot(Data.CH4.Temp.Best.res, smooth.residuals=F, colors = c('black', 'black'),
                          xlab="Water Temperature (°C)", ylab="CH4 Emissions (kg ha^-1 year^-1)",
                          main='')
print(Water_Temp.plotCH4)
dev.copy(pdf, "RESULTS/Water_Temp.plotCH4.pdf", width = 6.24, height = 5.06)
dev.off()
dev.copy(jpeg, "RESULTS/Water_Temp.plotCH4.jpeg", width = 6.24 * 300, height = 5.06 * 300, units = "px", res = 300)
dev.off()


# Check on the effects of factors on Carbon dioxide
Data_sum.CO2 = Data_sum %>% 
  filter(Owner == "WebbWare" | Season == "Summer")


if(FALSE){
  
  Data.CO2.Temp.1 = lme(data = Data_sum.CO2, 
                        log10(mgm2day_mean_CO2) ~ scale(log10(Temp_mean)) + 
                          scale(log10(Rainfall_mean+1)) +
                          scale(log10(Humidity_mean)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
  )
  
  Data.CO2.Temp.2 = lme(data = Data_sum.CO2, 
                        log10(mgm2day_mean_CO2) ~ scale(log10(Temp_mean)) + 
                          scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
  )
  
  Data.CO2.Temp.3 = lme(data = Data_sum.CO2, 
                        log10(mgm2day_mean_CO2) ~ scale(log10(Temp_mean)) + 
                          scale(log10(Rainfall_mean+1)) +
                          scale(log10(Humidity_mean)),
                        weight = varIdent(form = ~Temp_mean),
                        random = ~ 1|Owner,
  )
  
  Data.CO2.Temp.4 = lme(data = Data_sum.CO2, 
                        log10(mgm2day_mean_CO2) ~ scale(Temp_mean)+
                          scale(Water_Temp)+ 
                          scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        weight = varExp(form = ~Temp_mean),
                        random = ~ 1|Owner,
                        na.action = na.exclude,
  )
  
  Data.CO2.Temp.5 = lme(data = Data_sum.CO2, 
                        log10(mgm2day_mean_CO2) ~ scale(log10(Temp_mean))+ 
                          scale(log10(Water_Temp))+
                          scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        weight = varPower(form = ~Temp_mean),
                        random = ~ 1|Owner,
                        na.action = na.exclude,
  )
  
  
  Data.CO2.Temp.6 = lme(data = Data_sum.CO2, 
                        log10(mgm2day_mean_CO2) ~ scale(Temp_mean),# + 
                        #scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        weight = varConstProp(form = ~Temp_mean),
                        random = ~ 1|Owner,
  )
  
  
  Data.CO2.Temp.7 = lme(data = Data_sum.CO2, 
                        log10(mgm2day_mean_CO2) ~ scale(log10(Water_Temp))+  
                          scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
                        #method = "ML",
                        na.action = na.exclude,
  )
  
  
  Data.CO2.Temp.8 = lme(data = Data_sum.CO2, 
                        log10(mgm2day_mean_CO2) ~ scale(log10(Water_Temp))+
                          scale(log10(Temp_mean)),
                        #scale(log10(Rainfall_mean+1)),# +
                        #scale(log10(Humidity_mean)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
                        #method = "ML",
                        na.action = na.exclude,
  )
  
  
  Data.CO2.Temp.9 = lme(data = Data_sum.CO2, 
                        log10(mgm2day_mean_CO2) ~ scale(log10(Water_Temp))+
                          scale(log10(Temp_mean))+
                          scale(log10(Rainfall_mean+1))+# +
                          scale(log10(Humidity_mean))+
                          scale(log10(Depth))+
                          scale(log10(Surface_Ar)),
                        #scale(log10(Vegetation.20m.Area.m2.)),
                        #weight = varIdent(form = ~1|Treatment),
                        random = ~ 1|Owner,
                        #method = "ML",
                        na.action = na.exclude,
                        
  )
  
  AICc(Data.CO2.Temp.1,Data.CO2.Temp.2,Data.CO2.Temp.3,Data.CO2.Temp.4,Data.CO2.Temp.5,Data.CO2.Temp.6,
       Data.CO2.Temp.7,Data.CO2.Temp.8, Data.CO2.Temp.9)
  Data.CO2.Temp.Best = Data.CO2.Temp.8 
  Anova(Data.CO2.Temp.Best)
  plot(Data.CO2.Temp.8)
  save(Data.CO2.Temp.Best, file = "RESULTS/Data.CO2.Temp.Best.RData")
  
}


# MODEL VALIDATION
if(FALSE){
  Data_sum.CO2$Resid = residuals(Data.CO2.Temp.Best, type = "pearson")
  Data_sum.CO2$Fitted = predict(Data.CO2.Temp.Best)
  scatterplot(data = Data_sum.CO2, Resid~Fitted)
  
  # Verify the final model
  scatterplot(data = Data_sum.CO2, Resid~Temp_mean)
  boxplot(data = Data_sum.CO2, Resid~Owner)
  
}

##FIT THE BEST MODEL FOR THE OVERALL PLOT LATER
# First, remove NA in Stock from the dataset because makes issues below
Data_sum.CO2_noNA = Data_sum.CO2[!is.na(Data_sum.CO2$Water_Temp),]


# Fit an equivalent model with all raw values
Data.CO2.Temp.Best.plot = lme(data = Data_sum.CO2_noNA, 
                              mgm2day_mean_CO2 ~ Water_Temp +
                                Temp_mean,
                              random = ~ 1|Owner,
                              
)

Anova(Data.CO2.Temp.Best.plot)


# 1) Plot Temp_mean EFFECT CO2
Data.CO2.Temp.Best.res = Effect("Temp_mean", partial.residuals= T, Data.CO2.Temp.Best.plot)
Temp.plotCO2 = plot(Data.CO2.Temp.Best.res, smooth.residuals=F, colors = c('black', 'black'),
                    xlab="Air Temperature (°C)", ylab="CO2 Emissions (kg ha^-1 year^-1)",
                    main='')
print(Temp.plotCO2)
dev.copy(pdf, "RESULTS/Temp.plotCO2.pdf", width = 6.24, height = 5.06)
dev.off()
dev.copy(jpeg, "RESULTS/Temp.plotCO2.jpeg", width = 6.24 * 300, height = 5.06 * 300, units = "px", res = 300)
dev.off()


# 2) Plot water_temp EFFECT CO2
Data.CO2.Temp.Best.res = Effect("Water_Temp", partial.residuals= T, Data.CO2.Temp.Best.plot)
Water_Temp.plotCO2 = plot(Data.CO2.Temp.Best.res, smooth.residuals=F, colors = c('black', 'black'),
                          xlab="Water Temperature (°C)", ylab="CO2 Emissions (kg ha^-1 year^-1)",
                          main='')
print(Water_Temp.plotCO2)
dev.copy(pdf, "RESULTS/Water_Temp.plotCO2.pdf", width = 6.24, height = 5.06)
dev.off()
dev.copy(jpeg, "RESULTS/Water_Temp.plotCO2.jpeg", width = 6.24 * 300, height = 5.06 * 300, units = "px", res = 300)
dev.off()


# Check on the effects of factors on CO2eq
if(FALSE){
  
  Data.CO2eq.Temp.1 = lme(data = Data.CO2eq.WinterSummer, 
                          log.mgm2day_mean_CO2eq ~ scale(log10(Temp_mean)) + 
                            scale(log10(Rainfall_mean+1)) +
                            scale(log10(Humidity_mean)),
                          #weight = varIdent(form = ~1|Treatment),
                          random = ~ 1|Owner,
  )
  
  Data.CO2eq.Temp.2 = lme(data = Data.CO2eq.WinterSummer, 
                          log.mgm2day_mean_CO2eq ~ scale(log10(Temp_mean)) + 
                            scale(log10(Rainfall_mean+1)),# +
                          #scale(log10(Humidity_mean)),
                          #weight = varIdent(form = ~1|Treatment),
                          random = ~ 1|Owner,
  )
  
  Data.CO2eq.Temp.3 = lme(data = Data.CO2eq.WinterSummer, 
                          log.mgm2day_mean_CO2eq ~ scale(log10(Temp_mean)),# + 
                          #scale(log10(Rainfall_mean+1)),# +
                          #scale(log10(Humidity_mean)),
                          #weight = varIdent(form = ~1|Treatment),
                          random = ~ 1|Owner,
  )
  
  Data.CO2eq.Temp.4 = lme(data = Data.CO2eq.WinterSummer, 
                          log.mgm2day_mean_CO2eq ~ scale(log10(Temp_mean))+ 
                            scale(log10(Rainfall_mean+1))+
                            scale(log10(Humidity_mean))+
                            scale(log10(Water_Temp)),
                          #scale(log10(DO))+
                          #scale(log10(pH)),
                          weight = varExp(form = ~Temp_mean),
                          random = ~ 1|Owner,
                          na.action = na.exclude,
                          
  )
  
  Data.CO2eq.Temp.5 = lme(data = Data.CO2eq.WinterSummer, 
                          log.mgm2day_mean_CO2eq ~ scale(Temp_mean),# + 
                          #scale(log10(Rainfall_mean+1)),# +
                          #scale(log10(Humidity_mean)),
                          weight = varPower(form = ~Temp_mean),
                          random = ~ 1|Owner,
                          
  )
  
  
  Data.CO2eq.Temp.6 = lme(data = Data.CO2eq.WinterSummer, 
                          log.mgm2day_mean_CO2eq ~ scale(Temp_mean),# + 
                          #scale(log10(Rainfall_mean+1)),# +
                          #scale(log10(Humidity_mean)),
                          weight = varConstProp(form = ~Temp_mean),
                          random = ~ 1|Owner,
                          
  )
  
  
  Data.CO2eq.Temp.7 = lme(data = Data.CO2eq.WinterSummer, 
                          log.mgm2day_mean_CO2eq ~ scale(log10(Temp_mean))+#*(Treatment) + 
                            scale(log10(Rainfall_mean+1)) +
                            scale(log10(Humidity_mean)),
                          weight = varPower(form = ~Temp_mean),
                          random = ~ 1|Owner,
                          
  )
  
  
  Data.CO2eq.Temp.8 = lme(data = Data.CO2eq.WinterSummer, 
                          log.mgm2day_mean_CO2eq ~ scale(log10(Water_Temp)),# +  
                          #scale(log10(Rainfall_mean+1)),# +
                          #scale(log10(Humidity_mean)),
                          #weight = varIdent(form = ~1|Treatment),
                          random = ~ 1|Owner,
                          #method = "ML",
                          na.action = na.exclude,
  )
  
  
  Data.CO2eq.Temp.9 = lme(data = Data.CO2eq.WinterSummer, 
                          log.mgm2day_mean_CO2eq ~ scale(log10(Water_Temp))+#* (Treatment) +
                            scale(log10(Temp_mean)),
                          #scale(log10(Rainfall_mean+1)),# +
                          #scale(log10(Humidity_mean)),
                          #weight = varIdent(form = ~1|Treatment),
                          random = ~ 1|Owner,
                          #method = "ML",
                          na.action = na.exclude,
  )
  
  Data.CO2eq.Temp.10 = lme(data = Data.CO2eq.WinterSummer, 
                           log.mgm2day_mean_CO2eq ~ scale(log10(Water_Temp))+
                             scale(log10(Temp_mean))+
                             scale(log10(Rainfall_mean+1))+# +
                             scale(log10(Humidity_mean))+
                             scale(log10(Depth))+
                             scale(log10(Surface_Ar)),
                           #scale(log10(Vegetation.20m.Area.m2.)),
                           #weight = varIdent(form = ~1|Treatment),
                           random = ~ 1|Owner,
                           #method = "ML",
                           na.action = na.exclude,
                           
  )
  
  AICc(Data.CO2eq.Temp.1,Data.CO2eq.Temp.2,Data.CO2eq.Temp.3,Data.CO2eq.Temp.4,Data.CO2eq.Temp.5,Data.CO2eq.Temp.6,
       Data.CO2eq.Temp.7,Data.CO2eq.Temp.8,Data.CO2eq.Temp.9,Data.CO2eq.Temp.10)
  Data.CO2eq.Temp.Best = Data.CO2eq.Temp.9
  Anova(Data.CO2eq.Temp.Best)
  plot(Data.CO2eq.Temp.Best)
  save(Data.CO2eq.Temp.Best, file = "RESULTS/Data.CO2eq.Temp.Best.RData")
  
}

# MODEL VALIDATION
if(FALSE){
  Data.CO2eq.WinterSummer$Resid = residuals(Data.CO2eq.Temp.Best, type = "pearson")
  Data.CO2eq.WinterSummer$Fitted = predict(Data.CO2eq.Temp.Best)
  
  scatterplot(data = Data.CO2eq.WinterSummer, Resid~Fitted)
  
  # Verify the final model
  scatterplot(data = Data.CO2eq.WinterSummer, Resid~Temp_mean)
  boxplot(data = Data.CO2eq.WinterSummer, Resid~Owner)
  
}

##FIT THE BEST MODEL FOR THE OVERALL PLOT LATER
# First, remove NA in Stock from the dataset because makes issues below
Data_sum.CO2eq_noNA = Data.CO2eq.WinterSummer[!is.na(Data.CO2eq.WinterSummer$Water_Temp),]

# Fit an equivalent model with all raw values
Data.CO2eq.Temp.Best.plot = lme(data = Data_sum.CO2eq_noNA, 
                                mgm2day_mean_CO2eq ~ Water_Temp +
                                  Temp_mean,
                                random = ~ 1|Owner,
                                
)

Anova(Data.CO2eq.Temp.Best.plot)


# 1) Plot Temp_mean EFFECT CO2eq
Data.CO2eq.Temp.Best.res = Effect("Temp_mean", partial.residuals= T, Data.CO2eq.Temp.Best.plot)
Temp.plotCO2eq = plot(Data.CO2eq.Temp.Best.res, smooth.residuals=F, colors = c('black', 'black'),
                      xlab="Air Temperature (°C)", ylab="CO2eq Emissions (kg ha^-1 year^-1)",
                      main='')
print(Temp.plotCO2eq)
dev.copy(pdf, "RESULTS/Temp.plotCO2eq.pdf", width = 6.24, height = 5.06)
dev.off()
dev.copy(jpeg, "RESULTS/Temp.plotCO2eq.jpeg", width = 6.24 * 300, height = 5.06 * 300, units = "px", res = 300)
dev.off()


# 2) Plot water_temp EFFECT CO2eq
Data.CO2eq.Temp.Best.res = Effect("Water_Temp", partial.residuals= T, Data.CO2eq.Temp.Best.plot)
Water_Temp.plotCO2eq = plot(Data.CO2eq.Temp.Best.res, smooth.residuals=F, colors = c('black', 'black'),
                            xlab="Water Temperature (°C)", ylab="CO2eq Emissions (kg ha^-1 year^-1)",
                            main='')
print(Water_Temp.plotCO2eq)
dev.copy(pdf, "RESULTS/Water_Temp.plotCO2eq.pdf", width = 6.24, height = 5.06)
dev.off()
dev.copy(jpeg, "RESULTS/Water_Temp.plotCO2eq.jpeg", width = 6.24 * 300, height = 5.06 * 300, units = "px", res = 300)
dev.off()

##GENERAL COMMENT ON SECTION D
#We noticed that temperature(both air and water temperature) were the best-- 
#combination of environmental factors to explain seasonal GHG emissions


#(E) Plot all results and data
#############
# All plots #
#############

# # Figure 1 used in the manuscript
load(file = "RESULTS/Data.CH4.WinterSummer.Best.RData")  # Treat/season on CH4 (Data.CH4.WinterSummer.Best)
load(file = "RESULTS/Data.CO2.WinterSummer.Best.RData")   # Treat/season  on CO2 (Data.CO2.WinterSummer.Best)
load(file = "RESULTS/Data.CO2eq.WinterSummer.Best.RData") # Treat/season  on CO2eq (Data.CO2eq.WinterSummer.Best)


# Figure 2 used in the manuscript
load(file = "RESULTS/Data.N.WinterSummer.Best.RData")  # Treat on N (Data.N.WinterSummer.Best)
load(file = "RESULTS/Data.P.WinterSummer.Best.RData")  # Treat on P (Data.P.WinterSummer.Best)
load(file = "Results/Data.DO.WinterSummer.Best.RData") # Treat/season on DO (Data.DO.WinterSummer.Best)
load(file = "Results/Data.pH.WinterSummer.Best.RData")# Treat/season on pH (Data.pH.WinterSummer.Best)

# Figure 3 used in the manuscript
load(file = "RESULTS/Data.CH4.Temp.Best.RData")   # Temp on CH4 (Data.CH4.Temp.Best)
load(file = "RESULTS/Data.CO2.Temp.Best.RData")   # Temp on CO2 (Data.CO2.Temp.Best)
load(file = "RESULTS/Data.CO2eq.Temp.Best.RData") # Temp on CO2eq (Data.CO2eq.Temp.Best)


#Create plot functions for easy plotting
if(FALSE){
  
  #FIGURE 1
  # Adjust the dodge width if necessary to align with your bar width
  dodge_width <- 0.9
  
  PlotTemplate <- list(
    geom_bar(aes(fill = Treatment), position=position_dodge(dodge_width), stat="identity"),
    geom_linerange(aes(ymin = lower, ymax = upper, group = Treatment),
                   position=position_dodge(dodge_width), color="black", size=0.5),
    theme_bw(),
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
    )
  )
  
  PercCalc = function(model.df) {
    #model.df = Data.CH4.WinterSummer.Best.df
    model.df %>% group_by(Season) %>%
      summarise(Perc = 100*(fit[Treatment == "Fenced"]-fit[Treatment == "Unfenced"])/fit[Treatment == "Unfenced"]) %>%
      select(Perc) %>% unlist(c()) %>% round(1) %>% paste0(.,"%")
  }


# Model: Data.CH4.WinterSummer.Best
Data.CH4.WinterSummer.Best.df = as.data.frame(allEffects(Data.CH4.WinterSummer.Best))$`Treatment:Season`
# Convert from mg/m2/day to kg.ha.year
Data.CH4.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")] = 365.25*10000*(10^Data.CH4.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")])/1e+6
Data.CH4.WinterSummer.Best.df$Treatment = factor(Data.CH4.WinterSummer.Best.df$Treatment, levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))

####Plot
Plot.Data.CH4.WinterSummer.Best.df = 
  ggplot(data = Data.CH4.WinterSummer.Best.df, aes(x = Season, y = fit, fill = Treatment)) +
  PlotTemplate +
  labs(y = expression(paste("Methane emissions (kg CH"[4]," ha"^-1, " year"^-1, ")"))) +
  annotate(geom = "text", x = c(1,2), y = c(500, 50), label = PercCalc(Data.CH4.WinterSummer.Best.df), size = 6) +
  scale_y_log10() +
  theme(legend.position = c(0.80, 0.95), # Place the legend inside the plot
        legend.background = element_blank(), # Remove background of legend
        legend.key = element_blank(), # Remove key background
        legend.text = element_text(size = 12), # Reduce size of legend text
        legend.key.size = unit(0.5, 'cm'),  # Reduce size of legend keys
        #legend.title = element_text(size = 9)) #Reduce the legend tittle
        legend.title = element_blank(), # Remove the legend title
        axis.title.x = element_blank(), # Remove x-axis title
        axis.text.x = element_blank())  # Remove x-axis text

print(Plot.Data.CH4.WinterSummer.Best.df)

# Model: Data.CO2.WinterSummer.Best
Data.CO2.WinterSummer.Best.df = as.data.frame(allEffects(Data.CO2.WinterSummer.Best))$`Treatment:Season`
# Convert from mg/m2/day to kg.ha.year
Data.CO2.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")] = 365.25*10000*(10^Data.CO2.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")])/1e+6
Data.CO2.WinterSummer.Best.df$Treatment = factor(Data.CO2.WinterSummer.Best.df$Treatment, levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))

Plot.Data.CO2.WinterSummer.Best.df = 
  ggplot(data = Data.CO2.WinterSummer.Best.df, aes(x = Season, y = fit)) +
  PlotTemplate +
  labs(y = expression(paste("Carbon dioxide emissions (kg CO"[2]," ha"^-1, " year"^-1, ")"))) +
  scale_y_log10() +
  annotate(geom = "text", x = c(1.5), y = c(7000), label = c("no effect"), size = 6)+
  theme(legend.position = "none",#c(0.92, 0.97), # Place the legend inside the plot
        legend.background = element_blank(), # Remove background of legend
        legend.key = element_blank(), # Remove key background
        legend.text = element_text(size = 10), # Reduce size of legend text
        legend.key.size = unit(0.4, 'cm'),  # Reduce size of legend keys
        #legend.title = element_text(size = 9)) #Reduce the legend tittle
        legend.title = element_blank(), # Remove the legend title
        axis.title.x = element_blank(), # Remove x-axis title
        axis.text.x = element_blank())  # Remove x-axis text

print(Plot.Data.CO2.WinterSummer.Best.df)

# Model: Data.CO2eq.WinterSummer.Best
Data.CO2eq.WinterSummer.Best.df = as.data.frame(allEffects(Data.CO2eq.WinterSummer.Best))$`Treatment:Season`
# Convert from mg/m2/day to kg.ha.year
Data.CO2eq.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")] = 365.25*10000*(10^Data.CO2eq.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")])/1e+6
Data.CO2eq.WinterSummer.Best.df$Treatment = factor(Data.CO2eq.WinterSummer.Best.df$Treatment, 
                                                   levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))

Plot.Data.CO2eq.WinterSummer.Best.df = 
  ggplot(data = Data.CO2eq.WinterSummer.Best.df, aes(x = Season, y = fit)) +
  PlotTemplate +
  labs(y = expression(paste("Carbon dioxide-eq. emissions (kg CO"[2],"-eq ha"^-1, " year"^-1, ")"))) +
  annotate(geom = "text", x = c(1,2), y = c(17500, 3000), label = PercCalc(Data.CO2eq.WinterSummer.Best.df), size = 6) +
  scale_y_log10()+#(limits = c(0.1, max(Data.CO2eq.WinterSummer.Best.df$upper)))+
  #scale_y_continuous(limits = c(100, max(Data.CO2eq.WinterSummer.Best.df$upper)))+
  theme(legend.position = "none",#c(0.90, 0.97), # Place the legend inside the plot
        legend.background = element_blank(), # Remove background of legend
        legend.key = element_blank(), # Remove key background
        legend.text = element_text(size = 10), # Reduce size of legend text
        legend.key.size = unit(0.4, 'cm'),  # Reduce size of legend keys
        #legend.title = element_text(size = 9)) #Reduce the legend tittle
        legend.title = element_blank()) # Remove the legend title

print(Plot.Data.CO2eq.WinterSummer.Best.df)

##JOIN FIGURES INTO ONE DATAFRAME FOR MANUSCRIPT USING PATCHWORK OR COWPLOT
# Create with three plots side by side
Fig1_combined_plot = (Plot.Data.CH4.WinterSummer.Best.df + 
                      Plot.Data.CO2.WinterSummer.Best.df + 
                      Plot.Data.CO2eq.WinterSummer.Best.df) +
  plot_layout(ncol = 3, widths = c(1, 1, 1)) + 
  plot_annotation(tag_levels = c('A', 'B', 'C'),tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag.position = "topleft", plot.tag = element_text(face = "bold"))

# Print the combined plot
print(Fig1_combined_plot)

# When you print or save the plot, you can specify the dimensions to control the overall size
ggsave("RESULTS/Fig1_combined_plot.png", Fig1_combined_plot, width = 15, height = 5) # Adjust 'width' and 'height' as needed
ggsave("RESULTS/Fig1_combined_plot.pdf", Fig1_combined_plot, width = 15, height = 5)



#FIGURE 2
# Model: Data.DO.WinterSummer.Best
Data.DO.WinterSummer.Best.df = as.data.frame(allEffects(Data.DO.WinterSummer.Best))$`Treatment:Season`
Data.DO.WinterSummer.Best.df$Treatment = factor(Data.DO.WinterSummer.Best.df$Treatment,
                                                levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))

Plot.Data.DO.WinterSummer.Best.df = 
  ggplot(data = Data.DO.WinterSummer.Best.df, aes(x = Season, y = fit)) +
  PlotTemplate +
  labs(y = expression(paste("Dissolved Oxygen (mg L"^-1, ")"))) +
  #scale_y_log10() +
  annotate(geom = "text", x = c(1,2), y = c(10,13), label = PercCalc(Data.DO.WinterSummer.Best.df), size = 6)+
  theme(legend.position = "none",#c(0.92, 0.97), # Place the legend inside the plot
        legend.background = element_blank(), # Remove background of legend
        legend.key = element_blank(), # Remove key background
        legend.text = element_text(size = 10), # Reduce size of legend text
        legend.key.size = unit(0.4, 'cm'),  # Reduce size of legend keys
        #legend.title = element_text(size = 9)) #Reduce the legend tittle
        legend.title = element_blank(), # Remove the legend title
        axis.title.x = element_text(size = 14),  # Increase font size of x axis label
        axis.title.y = element_text(size = 14),  # Increase font size of y axis label
        axis.text.x = element_text(size = 11),   # X axis tick label font size
        axis.text.y = element_text(size = 11))   # Y axis tick label font size

print(Plot.Data.DO.WinterSummer.Best.df)


# Model: Data.pH.WinterSummer.Best
Data.pH.WinterSummer.Best.df = as.data.frame(allEffects(Data.pH.WinterSummer.Best))$`Treatment:Season`
Data.pH.WinterSummer.Best.df$Treatment = factor(Data.pH.WinterSummer.Best.df$Treatment,
                                                levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))
Plot.Data.pH.WinterSummer.Best.df = 
  ggplot(data = Data.pH.WinterSummer.Best.df, aes(x = Season, y = fit)) +
  PlotTemplate +
  labs(y = expression(paste("pH"))) +
  #scale_y_log10() +
  annotate(geom = "text", x = c(1,2), y = c(7.5,7.6), label = PercCalc(Data.pH.WinterSummer.Best.df), size = 6)+
  theme(legend.position = "none",#c(0.92, 0.97), # Place the legend inside the plot
        legend.background = element_blank(), # Remove background of legend
        legend.key = element_blank(), # Remove key background
        legend.text = element_text(size = 10), # Reduce size of legend text
        legend.key.size = unit(0.4, 'cm'),  # Reduce size of legend keys
        #legend.title = element_text(size = 9)) #Reduce the legend tittle
        legend.title = element_blank(), # Remove the legend title
        axis.title.x = element_text(size = 14),  # Increase font size of x axis label
        axis.title.y = element_text(size = 14),  # Increase font size of y axis label
        axis.text.x = element_text(size = 11),   # X axis tick label font size
        axis.text.y = element_text(size = 11))   # Y axis tick label font size

print(Plot.Data.pH.WinterSummer.Best.df)

#plotting nutrient (i.e.,Nitrogen and Phosphorus)
PercCalc2 = function(model.df) {
  #model.df = Data.N.WinterSummer.Best.df
  model.df %>% #group_by(Season) %>%
    summarise(Perc = 100*(fit[Treatment == "Fenced"]-fit[Treatment == "Unfenced"])/fit[Treatment == "Unfenced"]) %>%
    select(Perc) %>% unlist(c()) %>% round(1) %>% paste0(.,"%")
}

# Model: Data.N.WinterSummer.Best
Data.N.WinterSummer.Best.df = as.data.frame(allEffects(Data.N.WinterSummer.Best))$`Treatment`
# Convert from log10 to arithmentic
Data.N.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")] = 10^Data.N.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")]
Data.N.WinterSummer.Best.df$Treatment = factor(Data.N.WinterSummer.Best.df$Treatment, levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))

Plot.Data.N.WinterSummer.Best.df =
  ggplot(data = Data.N.WinterSummer.Best.df, aes(x = Treatment, y = fit, fill = Treatment)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_linerange(aes(ymin = lower, ymax = upper), position=position_dodge(0.85)) +
  theme_bw() +
  labs(y = expression(paste("Total nitrogen (mg L"^-1, ")"))) +
  #scale_y_log10() +
  annotate(geom = "text", x = c(1.5), y = c(15), label = PercCalc2(Data.N.WinterSummer.Best.df), size = 6) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),  # Added to remove grid lines
        legend.position = c(0.09, 0.95), # Place the legend inside the plot
        legend.background = element_blank(), # Remove background of legend
        legend.key = element_blank(), # Remove key background
        legend.text = element_text(size = 12), # Reduce size of legend text
        legend.key.size = unit(0.5, 'cm'),  # Reduce size of legend keys
        #legend.title = element_text(size = 9)) #Reduce the legend tittle
        legend.title = element_blank(), # Remove the legend title
        axis.title.x = element_text(size = 14),  # Increase font size of x axis label
        axis.title.y = element_text(size = 14),  # Increase font size of y axis label
        axis.text.x = element_text(size = 11),   # X axis tick label font size
        axis.text.y = element_text(size = 11))  # Y axis tick label font size

print(Plot.Data.N.WinterSummer.Best.df)

# Model: Data.P.WinterSummer.Best
Data.P.WinterSummer.Best.df = as.data.frame(allEffects(Data.P.WinterSummer.Best))$`Treatment`
# Convert from log10 to arithmentic
Data.P.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")] = 10^Data.P.WinterSummer.Best.df[,c("fit", "se", "lower", "upper")]
Data.P.WinterSummer.Best.df$Treatment = factor(Data.P.WinterSummer.Best.df$Treatment, levels = c("UD", "FD"), labels = c("Unfenced", "Fenced"))

Plot.Data.P.WinterSummer.Best.df =
  ggplot(data = Data.P.WinterSummer.Best.df, aes(x = Treatment, y = fit, fill = Treatment)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_linerange(aes(ymin = lower, ymax = upper), position=position_dodge(0.85)) +
  theme_bw() +
  labs(y = expression(paste("Total phosphorous (mg L"^-1, ")"))) +
  #scale_y_log10() +
  annotate(geom = "text", x = c(1.5), y = c(0.4), label = PercCalc2(Data.P.WinterSummer.Best.df), size = 6) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),  # Added to remove grid lines
        legend.position = "none",#c(0.09, 0.95), # Place the legend inside the plot
        legend.background = element_blank(), # Remove background of legend
        legend.key = element_blank(), # Remove key background
        legend.text = element_text(size = 12), # Reduce size of legend text
        legend.key.size = unit(0.5, 'cm'),  # Reduce size of legend keys
        #legend.title = element_text(size = 9)) #Reduce the legend tittle
        legend.title = element_blank(), # Remove the legend title
        axis.title.x = element_text(size = 14),  # Increase font size of x axis label
        axis.title.y = element_text(size = 14),  # Increase font size of y axis label
        axis.text.x = element_text(size = 11),   # X axis tick label font size
        axis.text.y = element_text(size = 11))  # Y axis tick label font size

print(Plot.Data.P.WinterSummer.Best.df)

##JOIN FIGURES INTO ONE DATAFRAME FOR MANUSCRIPT USING PATCHWORK OR COWPLOT
# Create with two plots each on two rows side by side
# Set a common theme to reduce plot margins for all plots
common_theme = theme(plot.margin = margin(5, 5, 5, 5))

Fig2_combined_plot = (Plot.Data.N.WinterSummer.Best.df + common_theme) +
  (Plot.Data.P.WinterSummer.Best.df + common_theme) +
  (Plot.Data.DO.WinterSummer.Best.df + common_theme) +
  (Plot.Data.pH.WinterSummer.Best.df + common_theme) +
  #plot_layout(ncol = 2, nrow = 2) & # Use & to apply globally
  plot_layout(ncol = 2, nrow = 2, widths = c(1, 1), heights = c(1, 1))&
  plot_annotation(tag_levels = 'A',tag_prefix = "(", tag_suffix = ")")&
  theme(plot.tag = element_text(face = "bold"))

print(Fig2_combined_plot)

# When you print or save the plot, you can specify the dimensions to control the overall size
ggsave("RESULTS/Fig2_combined_plot.png", Fig2_combined_plot, width = 13.5, height = 9) # Adjust 'width' and 'height' as needed
ggsave("RESULTS/Fig2_combined_plot.pdf", Fig2_combined_plot, width = 13.5, height = 9)



#FIGURE 3 
Anova(Data.CH4.Temp.Best.plot)
Anova(Data.CO2.Temp.Best.plot)
Anova(Data.CO2eq.Temp.Best.plot)

# Create a function to get the dataset to plot the effects in arithmetic scale
# -------------------
PlotArithmetic = function(mod){
  # mod = CH4_D.env.best.plot
  
  Terms = names(fixef(mod))[-1]  # Make sure this line is correctly obtaining the terms.
  print(Terms)
  
  
  PredictedValues.all = c()
  Residuals.all = c()
  
  for(t in 1:length(Terms))	{
    #t = 1
    
    # Associate a term
    term = Terms[t]
    
    # Get the effects
    options("na.action")
    Effects = Effect(term, partial.residuals=T, mod, na.action = na.exclude)
    
    # Get the predicted values
    PredictedValues = as.data.frame(Effects)
    names(PredictedValues)[1] = "x"
    
    # Add the name of the character
    PredictedValues$Val = as.factor(term)
    
    # Create a new dataset with the partial residuals
    Residuals = data.frame("y" = c(predict(mod, level = 0) + resid(Effects)))
    Residuals$y_raw = Effects$data[,Effects$response]
    Residuals$x = as.numeric(unlist(c(mod$data[,term])))
    Residuals$Val = as.factor(term)
    
    # Residuals from this post: https://stackoverflow.com/questions/43950459/use-ggplot-to-plot-partial-effects-obtained-with-effects-library
    closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)
    x.fit <- unlist(Effects$x.all)
    x <- data.frame(lower = Effects$lower, upper = Effects$upper, fit = Effects$fit, term = Effects$x[,term])
    xy <- data.frame(x = x.fit, y = x$fit[closest(x.fit, x$term)] + Effects$residuals)
    Residuals$y_NewRes = xy$y
    
    # Merge all values together	
    PredictedValues.all = rbind(PredictedValues.all, PredictedValues)
    Residuals.all = rbind(Residuals.all, Residuals)
    
  }
  
  return(list(Predicted = PredictedValues.all, Residuals = Residuals.all))
  
}

# Run the function
PlotData_CH4 = PlotArithmetic(Data.CH4.Temp.Best.plot)
PlotData_CO2 = PlotArithmetic(Data.CO2.Temp.Best.plot)
PlotData_CO2eq = PlotArithmetic(Data.CO2eq.Temp.Best.plot)

# Add the info on which model
PlotData_CH4$Predicted$Model = "CH4"
PlotData_CH4$Residuals$Model = "CH4"

PlotData_CO2$Predicted$Model = "CO2"
PlotData_CO2$Residuals$Model = "CO2"

PlotData_CO2eq$Predicted$Model = "CO2eq"
PlotData_CO2eq$Residuals$Model = "CO2eq"

# Combine together
AllPlotData_Residuals = rbind(PlotData_CH4$Residuals,PlotData_CO2$Residuals,PlotData_CO2eq$Residuals)
AllPlotData_Predicted = rbind(PlotData_CH4$Predicted,PlotData_CO2$Predicted,PlotData_CO2eq$Predicted)

# Save p values
Pvalues_CH4 = as.data.frame(Anova(Data.CH4.Temp.Best.plot))
Pvalues_CH4$levels = rownames(Pvalues_CH4)

Pvalues_CO2 = as.data.frame(Anova(Data.CO2.Temp.Best.plot))
Pvalues_CO2$levels = rownames(Pvalues_CO2)

Pvalues_CO2eq = as.data.frame(Anova(Data.CO2eq.Temp.Best.plot))
Pvalues_CO2eq$levels = rownames(Pvalues_CO2eq)

# Add info on models
Pvalues_CH4$Model = "CH4"
Pvalues_CO2$Model = "CO2"
Pvalues_CO2eq$Model = "CO2eq"

# Combine
Pvalues = rbind(Pvalues_CH4,Pvalues_CO2,Pvalues_CO2eq)
rownames(Pvalues) = NULL

names(Pvalues)[3] = "p"
Pvalues$p_round = round(Pvalues$p,3)
Pvalues$p_round2 = ifelse(Pvalues$p_round < 0.01, "p < 0.01", paste0("p = ", Pvalues$p_round))

# Finally we identify and remove the non-sign. fitted lines
AllPlotData_Predicted$ModelVal = interaction(AllPlotData_Predicted$Model, AllPlotData_Predicted$Val)
ModelVal_NonSign = c("CH4.Water_Temp", "CO2.Water_Temp", "CO2.Temp_mean", "CO2eq.Water_Temp")
AllPlotData_Predicted$fit_2 = ifelse(AllPlotData_Predicted$ModelVal %in% ModelVal_NonSign, NA, AllPlotData_Predicted$fit)
AllPlotData_Predicted$lower_2 = ifelse(AllPlotData_Predicted$ModelVal %in% ModelVal_NonSign, NA, AllPlotData_Predicted$lower)
AllPlotData_Predicted$upper_2 = ifelse(AllPlotData_Predicted$ModelVal %in% ModelVal_NonSign, NA, AllPlotData_Predicted$upper)

# ----------------
# Labelling - Log scale
# ----------------

# Create the parsed lables in the right order
ParseLabeled_log = c(expression("Air Temperature ("*~degree*C*")"),
                     expression("Water Temperature ("*~degree*C*")")
                     
)

ParseModel_log = c(expression(paste("CH"[4], " flux (Kg ha"^"-1", " y"^-1, ")")),
                   expression(paste("CO"[2], " flux (Kg ha"^"-1", " y"^-1, ")")),
                   expression(paste("CO"[2],"-eq flux (Kg ha"^"-1", " y"^-1, ")"))
)


# Label the residuals
AllPlotData_Residuals$Val_parsed_log = factor(AllPlotData_Residuals$Val, 
                                              levels = c("Temp_mean","Water_Temp"),
                                              labels = ParseLabeled_log
)

# Label the predicted values
AllPlotData_Predicted$Val_parsed_log = factor(AllPlotData_Predicted$Val, 
                                              levels = c("Temp_mean","Water_Temp"),
                                              labels = ParseLabeled_log
)

# Label the models
AllPlotData_Predicted$Model_parsed_log = factor(AllPlotData_Predicted$Model, 
                                                levels = c("CH4","CO2","CO2eq"),
                                                labels = ParseModel_log
)

AllPlotData_Residuals$Model_parsed_log = factor(AllPlotData_Residuals$Model, 
                                                levels = c("CH4","CO2","CO2eq"),
                                                labels = ParseModel_log
)

# Label the p-values
Pvalues$Val_parsed_log = factor(Pvalues$levels, 
                                levels = c("Temp_mean","Water_Temp"),
                                labels = ParseLabeled_log
)

Pvalues$Model_parsed_log = factor(Pvalues$Model, 
                                  levels = c("CH4","CO2","CO2eq"),
                                  labels = ParseModel_log
)

# Turn from mg to kg/ha
conversion_factor <- 365 / 100

AllEffectPlots_g_log10 <- ggplot() +
  geom_point(data = AllPlotData_Residuals, aes(x = x, y = (y_NewRes * conversion_factor))) +
  geom_line(data = AllPlotData_Predicted, aes(x = x, y = (fit_2 * conversion_factor))) +
  geom_ribbon(data = AllPlotData_Predicted, aes(x = x, ymin = (lower_2 * conversion_factor), ymax = (upper_2 * conversion_factor)), alpha = 0.2) +
  facet_grid(Model_parsed_log ~ Val_parsed_log, scales = "free", labeller = labeller(Model_parsed_log = label_parsed, Val_parsed_log = label_parsed), switch = "y") +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14)) +
  geom_text(data = Pvalues, aes(x = -Inf, y = -Inf, label = p_round2), parse = FALSE, inherit.aes = FALSE, hjust = -0.25, vjust = -1, size = 6, color = "dark grey") +
  scale_y_continuous(position = "right",labels = scales::label_number())  # Assuming your data does not need a log scale anymore
print(AllEffectPlots_g_log10)

ggsave(AllEffectPlots_g_log10, file = "RESULTS/AllEffectPlots_g_log10_3mods.pdf", width = 10, height  = 8)
ggsave(AllEffectPlots_g_log10, file = "RESULTS/AllEffectPlots_g_log10_3mods.png", width = 10, height  = 8)


}

#End of analyses
