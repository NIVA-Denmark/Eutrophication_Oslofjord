library(data.table)
library(extrafont)
library(extrafontdb)
library(hms)
library(lubridate)
library(patchwork)
library(readxl)
library(sf)
library(terra)
library(tidyverse)

output_folder <- "C:/your_folder/output/"
input_folder <- "C:/your_folder/input/"
data <- data_with_thresholds


# ---------- vfk information  ---------------
shp <- sf::st_read(paste0(input_folder, "gis/oslofjord_waterbodies.shp")) 
shp$geometry <- NULL

pts_coords <- shp %>%
  select(Vannforeko, A_group, Vannfore_1, x_m, y_m)%>%
  mutate(vfk = paste0(Vannfore_1, " [", Vannforeko, "]")) %>%
  arrange(A_group, desc(y_m), desc(x_m))

list_vfk <- pts_coords$vfk

rm(shp)

xlimits <- c(9.245,11.535)
ylimits <- c(58.781,59.964)

IndicatorLevels<- c("Softbottom fauna", "MSMDI","summer klorofyll a",
                    "po4-p, sommer", "po4-p, vinter", "no3+no2-n, sommer", "no3+no2-n, vinter", "totn, sommer", "totn, vinter", "totp, sommer","totp, vinter", "nh4-n, sommer", "nh4-n, vinter","siktdyp, sommer", "oksygen","oksygenmetning")

StatusLevels <- c("High","Good","Mod","Poor","Bad")
color_palette <- c("High"="#2C96F6", 
                   "Good"="#8FFF07",
                   "Mod"="#FFFF00" ,
                   "Poor"="#FAB70E",
                   "Bad"="#FF2121")

areas<-data.frame(A_group= c(1:6),
                  Basin_name=c("1 Inner fjord",
                               "2 Middle fjord",
                               "3 Outer fjord",
                               "4 Drammensfjord",
                               "5 West sidefjords",
                               "6 East sidefjords"))

area_f<-unique(areas$Basin_name)



source(paste0(input_folder, "../EQR_functions.R"))

# Heatmap by basin ----
  df <- data %>%
  rename("value"="Avg")
  
  
  df<-df%>%
    CalcEQR()
  # get quality element names
  df <- df %>%
    GetQEs()
  
  
  df<-df%>%
    mutate(HEATCat= ifelse(Parameter=="summer klorofyll a","Pri",HEATCat ),
           HEATCat= ifelse(Parameter=="Softbottom fauna","Sec",HEATCat ))%>%
    mutate(QE= ifelse(Parameter=="summer klorofyll a","Phytoplankton",QE),
           QE= ifelse(Parameter=="Softbottom fauna","Bundfauna",QE))%>%
    mutate(QEtype= ifelse(Parameter=="summer klorofyll a","Bio",QEtype ),
           QEtype= ifelse(Parameter=="Softbottom fauna","Bio",QEtype ))%>%
    left_join(pts_coords,by=c("Vannforeko","Vannfore_1"))
  
  
  
  df<-df%>%
    filter(!is.na(A_group))
  
  # calculate HEAT eutrophication ratio - using the G/M boundary as threshold
  df <- df %>%
    CalcER()
  df1<-df%>%
    mutate(ER=EQR)
  
  #HEAT Heatmap
  # do aggregation and overall classification using HEAT
  dfHEAT <- aggregateHEAT(df1, var_cat="HEATCat",group_vars=c("A_group","Year","Parameter"))
  
  dfHEAT <- dfHEAT %>% 
    left_join(areas,by="A_group")%>%
    mutate(Basin=factor(Basin_name,levels=area_f))
  
  years<-min(dfHEAT$Year,na.rm=T):max(dfHEAT$Year,na.rm=T)
  
  dfHEAT <- tibble(Year=years) %>% 
    merge(tibble(Basin=area_f),all=T) %>%
    left_join(dfHEAT, by=c("Year","Basin")) 
  
  dfHEAT$Year <- factor(dfHEAT$Year, levels=years)
  dfHEAT$Class <- factor(dfHEAT$HEATClass, levels = c("Bad","Poor", "Mod","Good","High"))
  
  
  dfHEAT<-dfHEAT%>%
    mutate(Criteria= case_when(
      Worst =="Nut" ~ "Nutrient levels",
      Worst =="Pri" ~ "Direct effects",
      Worst =="Sec" ~ "Indirect effects",
      TRUE ~ NA
    ))%>%mutate(Year=as.character(Year),Year=as.numeric(Year))
  
  dfHEAT1 <- tibble(Year=years) %>% 
    merge(tibble(Basin=area_f),all=T) %>%
    merge(tibble(Parameter=IndicatorLevels),all=T)%>%
    left_join(dfHEAT, by=c("Year","Basin","Parameter")) %>%
    mutate(Parameter=factor(Parameter,levels=rev(IndicatorLevels)))%>%
    group_by(Basin)%>%
    fill(A_group,.direction = "updown")%>%ungroup()
  
  dfHEAT1<-dfHEAT1%>%filter(Year>=1960)
  for(i in c(1:6)){
    dfHEATx<-dfHEAT1%>%filter(A_group==i)
    p<-ggplot(dfHEATx, aes(Year, Parameter, fill=Class))+ 
      geom_tile(color="grey", show.legend=F) +
      labs(x="",y="",title=unique(dfHEATx$Basin))+
      scale_fill_manual(values=color_palette,drop=F,breaks=StatusLevels,na.value = "white") +
      scale_x_continuous(expand=c(0,0),
                         breaks=c(1960,1970,1980,1990,2000,2010,2020),
                         labels= c("1960", "1970", "1980", "1990", "2000", "2010", "2020"))+
      scale_y_discrete(
        labels=c(
          "summer klorofyll a"  = "Chlorophyll, summer",
          "siktdyp, sommer" = "Secchi depth",
          "Softbottom fauna" = "Softbottom fauna",
          "MSMDI" = "MSMDI Hardbottom fauna",
          "oksygen" = "Bottom oxygen",
          "oksygenmetning" = "Oxygen saturation",
          "po4-p, sommer" = "Phosphate, summer",
          "po4-p, vinter" = "Phosphate, winter",
          "no3+no2-n, sommer" = "Sum NOx, summer",
          "no3+no2-n, vinter" = "Sum NOx, winter",
          "nh4-n, sommer"= "Ammonium, summer",
          "nh4-n, vinter"= "Ammonium, winter",
          "totn, sommer" = "Total N, summer",
          "totn, vinter" = "Total N, winter",
          "totp, sommer" = "Total P, summer",
          "totp, vinter" = "Total P, winter"))+
      theme_bw() +
      theme(
        text=element_text(family="Calibri"),
        axis.text.x=element_text(size=18,colour="#000000",hjust=0,vjust=0.5),
        axis.text.y=element_text(size=18,vjust = 0.2,colour="#000000"),
        axis.ticks.y=element_line(linewidth=0.4),
        axis.ticks.x=element_line(linewidth=0.4),
        title = element_text(size=19),
        panel.grid.major.y = element_blank(),
        legend.text = element_text(size=24),
        legend.title = element_text(size=24),
        legend.position = "bottom",
        plot.title = element_text(hjust=0.5),
        plot.margin = margin(r=20,t=10)
      ) +  
      geom_hline(yintercept = c(14.5,10.5),linetype="dashed")
    p
    ggsave(p, file=paste0(output_folder,"Heatmap basin ",i,".png"), height=15, width=40, units="cm")
  }
  
  
# Heatmap by vannforekomst ----
  
  df <- data %>%
    rename("value"="Avg")
  
  
  df<-df%>%
    CalcEQR()
  # get quality element names
  df <- df %>%
    GetQEs()
  
  
  df<-df%>%
    mutate(HEATCat= ifelse(Parameter=="summer klorofyll a","Pri",HEATCat ),
           HEATCat= ifelse(Parameter=="Softbottom fauna","Sec",HEATCat ))%>%
    mutate(QE= ifelse(Parameter=="summer klorofyll a","Phytoplankton",QE),
           QE= ifelse(Parameter=="Softbottom fauna","Bundfauna",QE))%>%
    mutate(QEtype= ifelse(Parameter=="summer klorofyll a","Bio",QEtype ),
           QEtype= ifelse(Parameter=="Softbottom fauna","Bio",QEtype ))%>%
    left_join(pts_coords,by=c("Vannforeko","Vannfore_1"))
  
  
  df<-df%>%
    filter(!is.na(A_group))
  
  # calculate HEAT eutrophication ratio - using the G/M boundary as threshold
  df <- df %>%
    CalcER()
  df1<-df%>%
    mutate(ER=EQR)
  
  #HEAT Heatmap
  # do aggregation and overall classification using HEAT
  dfHEAT <- aggregateHEAT(df1, var_cat="HEATCat",group_vars=c("Vannforeko","Year","Parameter"))
  
  dfHEAT <- dfHEAT %>%
    left_join(pts_coords,by="Vannforeko")%>%
    mutate(Basin=factor(vfk,levels=list_vfk))
  
  years<-min(dfHEAT$Year,na.rm=T):max(dfHEAT$Year,na.rm=T)
  
  dfHEAT <- tibble(Year=years) %>%
    merge(tibble(Basin=list_vfk),all=T) %>%
    left_join(dfHEAT, by=c("Year","Basin"))
  
  dfHEAT$Year <- factor(dfHEAT$Year, levels=years)
  dfHEAT$Class <- factor(dfHEAT$HEATClass, levels = c("Bad","Poor", "Mod","Good","High"))
  
  
  dfHEAT<-dfHEAT%>%
    mutate(Criteria= case_when(
      Worst =="Nut" ~ "Nutrient levels",
      Worst =="Pri" ~ "Direct effects",
      Worst =="Sec" ~ "Indirect effects",
      TRUE ~ NA
    ))%>%mutate(Year=as.character(Year),Year=as.numeric(Year))
  
  dfHEAT1 <- tibble(Year=years) %>%
    merge(tibble(Basin=list_vfk),all=T) %>%
    merge(tibble(Parameter=IndicatorLevels),all=T)%>%
    left_join(dfHEAT, by=c("Year","Basin","Parameter")) %>%
    mutate(Parameter=factor(Parameter,levels=rev(IndicatorLevels)))
  
  add_basin<-shp%>%select(A_group,Vannforeko,Vannfore_1)
  add_basin$geometry<-NULL
  add_basin<-add_basin%>%
    merge(tibble(Basin=list_vfk),all=T)%>%
    rowwise()%>%
    mutate(keep= ifelse(grepl(Vannforeko,Basin),0,1))%>%
    ungroup()%>%
    filter(keep==0)%>%
    select(-keep,-Vannforeko,-Vannfore_1)
  
  dfHEAT1<-dfHEAT1%>%
    select(-A_group)%>%
    left_join(add_basin, by="Basin")
  
  dfHEAT1<-dfHEAT1%>%filter(Year>=1960)
  
  for(vannforeko in list_vfk){
    dfHEATx<-dfHEAT1%>%filter(Basin==vannforeko)
    bas<-unique(dfHEATx$A_group)
    
    p<-ggplot(dfHEATx, aes(Year, Parameter, fill=Class))+
      geom_tile(color="grey", show.legend=F) +
      labs(x="",y="",title=unique(dfHEATx$Basin))+
      scale_fill_manual(values=color_palette,drop=F,breaks=StatusLevels,na.value = "white") +
      scale_x_continuous(expand=c(0,0),
                         breaks=c(1960,1970,1980,1990,2000,2010,2020),
                         labels= c("1960", "1970", "1980", "1990", "2000", "2010", "2020"))+
      scale_y_discrete(
        labels=c(
          "Softbottom fauna" = "Softbottom fauna",
          "MSMDI" = "MSMDI Hardbottom fauna",
          "summer klorofyll a"  = "Chlorophyll, summer",
          "po4-p, sommer" = "Phosphate, summer",
          "po4-p, vinter" = "Phosphate, winter",
          "no3+no2-n, sommer" = "Sum NOx, summer",
          "no3+no2-n, vinter" = "Sum NOx, winter",
          "totn, sommer" = "Total N, summer",
          "totn, vinter" = "Total N, winter",
          "totp, sommer" = "Total P, summer",
          "totp, vinter" = "Total P, winter",
          "nh4-n, sommer"= "Ammonium, summer",
          "nh4-n, vinter"= "Ammonium, winter",
          "siktdyp, sommer" = "Secchi depth",
          "oksygen" = "Bottom oxygen",
          "oksygenmetning" = "Oxygen saturation"))+
      theme_bw() +
      theme(
        text=element_text(family="Calibri"),
        axis.text.x=element_text(size=18,colour="#000000",hjust=0,vjust=0.5),
        axis.text.y=element_text(size=18,vjust = 0.2,colour="#000000"),
        axis.ticks.y=element_line(linewidth=0.4),
        axis.ticks.x=element_line(linewidth=0.4),
        title = element_text(size=19),
        panel.grid.major.y = element_blank(),
        legend.text = element_text(size=24),
        legend.title = element_text(size=24),
        legend.position = "bottom",
        plot.title = element_text(hjust=0.5)
      ) +
      geom_hline(yintercept = c(14.5,10.5),linetype="dashed")
    p
    ggsave(p, file=paste0(output_folder,"Vannforekomst/Heatmap basin ", bas," ",vannforeko,".png"), height=15, width=40, units="cm")

  }
  

# Criteria Integrated by basin -----
  df <- data %>%
    rename("value"="Avg")
  
  df<-df%>%
    CalcEQR()
  # get quality element names
  df <- df %>%
    GetQEs()
  df$Parameter<-tolower(df$Parameter)
  
  df<-df%>%
    mutate(HEATCat= ifelse(Parameter=="summer klorofyll a","Pri",HEATCat ),
           HEATCat= ifelse(Parameter=="softbottom fauna","Sec",HEATCat ))%>%
    mutate(QE= ifelse(Parameter=="summer klorofyll a","Phytoplankton",QE),
           QE= ifelse(Parameter=="softbottom fauna","Bundfauna",QE))%>%
    mutate(QEtype= ifelse(Parameter=="summer klorofyll a","Bio",QEtype ),
           QEtype= ifelse(Parameter=="softbottom fauna","Bio",QEtype ))%>%
    left_join(pts_coords,by=c("Vannforeko","Vannfore_1"))
  
  
  df<-df%>%
    filter(!is.na(A_group))
  
  # calculate HEAT eutrophication ratio - using the G/M boundary as threshold
  df <- df %>%
    CalcER()
  df1<-df%>%
    mutate(ER=EQR)
  #integrated by basin + confidence interval
  
  
  random_indicator_values <- function(avg, se, logavg, logse, param="", nsim=1000){
    
    list_log <- c("po4-p, sommer","no3-n, sommer","totn, sommer","totp, sommer",
                  "no3-n, vinter","po4-p, vinter","totn, vinter","totp, vinter",
                  "nh4-n, sommer","nh4-n, vinter","summer klorofyll a")

    
    se <- ifelse(is.na(se), 0, se) 
    
    se <- ifelse(param %in% list_log, logse , se)
    avg <- ifelse(param %in% list_log, logavg, avg)
    
    vals <- rnorm(n=nsim, mean=avg, sd=se)
    param <- rep(param, nsim)
    vals <- ifelse(param %in% list_log, exp(vals)-0.001, vals)
    
    return(vals)
  }

  #we need to keep the "true" averaged values so we can compare them to the CV
  df_mean <- df1 %>%
    select(A_group, Type, Year,QE, HEATCat, Parameter, Indicator, mean=value, ER, Class)%>%filter(!is.na(A_group))
  
  #Using the function "random_indicator_values" we apply a Monte Carlo simulation. This simulation takes the actual value and its SE
  #and calculates 1000 random values that fit with the ones provided. 
  df_sim <- df1 %>%
    filter(!is.na(A_group))%>%
    rowwise() %>%
    mutate(val_sim=list(random_indicator_values(value, SE, logAvg, logSE, Parameter))) %>%
    ungroup() %>%
    select(A_group, Type, Year,QE, HEATCat, Parameter, Indicator, MatchValue, Match, Ref, HG, GM, MP, PB, Worst, Resp, val_sim)
  
  #this simulation gives you 1000 simulated values in a list, we need to unpack it. Therefore, for x value, now we have 1000 values.
  #the simulation is made by row, in another words, for each type of waterbody, for each year, and each parameter, we will get 1000 simulations.
  
  df_sim <- tidyr::unnest(df_sim, cols=val_sim) %>%
    rename(value=val_sim) %>%
    group_by(A_group, Type, Year, QE,HEATCat, Parameter, Indicator, MatchValue, Match, Ref, HG, GM, MP, PB, Worst, Resp) %>%
    mutate(simID=row_number()) %>%
    ungroup()
  
  #now, from the simulated values, we calculate a simulated ER value.
  df_sim <- df_sim %>%
    CalcEQR(getClass=F)%>%
    mutate(ER=EQR)
  
  
  
  df_sim$Class <- factor(df_sim$Class, levels = c("Bad","Poor", "Mod","Good","High"))
  
  #We use the aggregation function, aggregating by year and simID
  df_sim<-df_sim%>%
    mutate(Criteria= case_when(
      HEATCat =="Nut" ~ "Nutrients",
      HEATCat =="Pri" ~ "Direct effects",
      HEATCat =="Sec" ~ "Indirect effects",
      TRUE ~ NA
    ))
  
  df_mean<-df_mean%>%
    mutate(Criteria= case_when(
      HEATCat =="Nut" ~ "Nutrients",
      HEATCat =="Pri" ~ "Direct effects",
      HEATCat =="Sec" ~ "Indirect effects",
      TRUE ~ NA
    ))
  
  dfHEATsim <- aggregateHEAT(df_sim, var_cat="HEATCat",group_vars=c("Year","simID","A_group","Criteria"))
  dfHEATmean <- aggregateHEAT(df_mean, var_cat="HEATCat",group_vars=c("Year","A_group","Criteria"))
  
  rm(df_sim,df_mean)
  gc()
  dfHEATmean<-dfHEATmean%>%
    left_join(areas,by="A_group")%>%mutate(Basin_name=factor(Basin_name,levels=area_f))
  
  #to avoid joining years too far away from each other we make a new df with all the years, to have the gaps.
  dfYears <- data.frame(YearAvg=(min(dfHEATsim$Year):max(dfHEATsim$Year)))
  dfCriteria<-data.frame(Criteria= c("Direct effects","Indirect effects","Nutrients"))
  
  dfCombined<-expand.grid(YearAvg=dfYears$YearAvg,Criteria=dfCriteria$Criteria)
  
  dfHEATsim_avg<- dfCombined %>%
    merge(dfHEATsim, all=T)
  
  #now we have all combinations possible, i.e. for 1933 we have a row for years from 1933 to 2023.
  #Of course, we do not need all this combinations,we are only interested in the year, and its previous 
  # 5 years (year + 4 previous years) so we can calculate the standard error for that year.
  dfHEATsim_avgx<- dfHEATsim_avg %>%
    filter(Year<=YearAvg, Year>=(YearAvg-4)) 
  
  dfHEATsim_avgx<-dfHEATsim_avgx%>%
    group_by(YearAvg,A_group,Criteria) %>%
    summarise(ERmean=mean(ER), SD=sd(ER), SE=sd(ER)/sqrt(n()), .groups="drop")
  
  dfHEATsim_avgx<- split(dfHEATsim_avgx,dfHEATsim_avgx$A_group)
  
  dfHEATsim_avgx <- dfHEATsim_avgx %>%
    map(~right_join(.x, dfCombined, by = c("YearAvg","Criteria")) %>%
          mutate(A_group = unique(.x$A_group[!is.na(.x$A_group)]))
    )%>%
    bind_rows()%>%
    rename(Year=YearAvg)
  
  #to plot the bands for the confidence interval, we need a minimum and a max
  dfHEATsim_avgx <- dfHEATsim_avgx %>%
    mutate(ERmax=ERmean+SD, ERmin=ERmean-SD)
  dfHEATsim_avgx <- dfHEATsim_avgx%>%
    left_join(areas,by="A_group")%>%mutate(Basin_name=factor(Basin_name,levels=area_f))
  
  dfHEATsim_avgx <- dfHEATsim_avgx%>%filter(Year>=1960)
  
  #we do the same process for the true values, so we can plot the simulated one and the original ones.
  dfHEATmean5yr <- dfCombined %>%
    merge(dfHEATmean, all=T)
  
  
  dfHEATmean5yr<- dfHEATmean5yr %>%
    filter(Year<=YearAvg, Year>=(YearAvg-4)) %>%
    group_by(YearAvg,A_group,Criteria) %>%
    summarise(ERmean=mean(ER), SD=sd(ER), SE=sd(ER)/sqrt(n()), n=n(), .groups="drop")  
  
  dfHEATmean5yr<- split(dfHEATmean5yr,dfHEATmean5yr$A_group)
  
  dfHEATmean5yr <- dfHEATmean5yr %>%
    map(~right_join(.x, dfCombined, by = c("YearAvg","Criteria")) %>%
          mutate(A_group = unique(.x$A_group[!is.na(.x$A_group)]))
    )%>%
    bind_rows()%>%
    rename(Year=YearAvg)
  
  dfHEATmean5yr <- dfHEATmean5yr%>%
    left_join(areas,by="A_group")%>%mutate(Basin_name=factor(Basin_name,levels=area_f))
  
  dfHEATmean5yr <- dfHEATmean5yr%>%filter(Year>=1960)
  #we plot both, the simulation and the true values
  dfHEATmean<-dfHEATmean%>%filter(Year>=1960)
  for(i in c(1:6)){
    bas<-dfHEATsim_avgx%>%
      filter(A_group==i)%>%
      distinct(as.character(Basin_name))
    bas<-bas[[1]]
    p<-ggplot(subset(dfHEATsim_avgx,A_group==i)) + 
      # geom_rect(aes(ymin = 0.8, ymax = Inf, xmin = -Inf, xmax = Inf), fill="#2C96F6",color=NA, alpha = 0.2) +
      # geom_rect(aes(ymin = 0.6, ymax = 0.8, xmin = -Inf, xmax = Inf), fill = "#8FFF07",color=NA, alpha = 0.2) +
      # geom_rect(aes(ymin = 0.4, ymax = 0.6, xmin = -Inf, xmax = Inf), fill = "#FFFF00",color=NA, alpha = 0.2) +
      # geom_rect(aes(ymin = 0.2, ymax = 0.4, xmin = -Inf, xmax = Inf), fill = "#FAB70E",color=NA, alpha = 0.2) +
      # geom_rect(aes(ymin = -Inf, ymax = 0.2, xmin = -Inf, xmax = Inf), fill = "#FF2121",color=NA, alpha = 0.2) +
      labs(title =paste0(bas),y="EQR",x="Year")+
      geom_ribbon(aes(x=Year, ymin=ERmin, ymax=ERmax), alpha=0.2) +
      geom_line(data=subset(dfHEATmean5yr,A_group==i), aes(x=Year, y=ERmean), colour="black") +
      geom_point(data=subset(dfHEATmean,A_group==i), aes(x=Year, y=ER),shape=1)+
      scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
      scale_x_continuous(limits = c(1955,2025),
                         labels= c("1960","1980","2000","2020"),
                         breaks = c(1960,1980,2000,2020))+
      geom_hline(yintercept = c(0.2,0.4,0.6,0.8),linetype="dashed")+
      facet_wrap(~Criteria,ncol=3)+
      theme_minimal()+
      theme(
        plot.background = element_rect(fill="white",colour = "darkgrey",linewidth=1.5),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey",fill=NA),
        text = element_text(size=20,family = "Calibri"),
        plot.title= element_text(size=35,hjust=0.5),
        strip.text = element_text(size=30),
        axis.title = element_text(size=30),
        axis.text = element_text(size=25),
        axis.ticks.x = element_line()) +
      coord_cartesian(ylim=c(0,1),expand = T)
    
    
    ggsave(p, file=paste0(output_folder,"Criteria basin ",i,".png"), height=15, width=40, units="cm")
  }
