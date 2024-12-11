library(tidyverse)
library(readxl)
library(lubridate)
library(hms)
library(terra)
library(sf)
library(patchwork)
library(data.table)
library(extrafont)
library(extrafontdb)
  
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
  
# Fig 1 - Catchment area ----
  
  #This has been done in ArcGis

# Fig 3 - Integrated assessment by basins ----

  
  # calculate EQR values and status class
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
    select(A_group, Type, Year,QE, HEATCat, Parameter, Indicator, mean=value, ER, Class)%>%
    filter(!is.na(A_group))
  
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
    CalcEQR()
  df_sim <- df_sim %>%
    mutate(ER=EQR)
  
  
  
  df_sim$Class <- factor(df_sim$Class, levels = c("Bad","Poor", "Mod","Good","High"))
  
  
  #We use the aggregation function, aggregating by year and simID
  
  dfHEATsim <- aggregateHEAT(df_sim, var_cat="HEATCat",group_vars=c("Year","simID","A_group"))
  dfHEATmean <- aggregateHEAT(df_mean, var_cat="HEATCat",group_vars=c("Year","A_group"))
  dfHEATmean<-dfHEATmean%>%
    left_join(areas,by="A_group")%>%mutate(Basin_name=factor(Basin_name,levels=area_f))
  
  #to avoid joining years too far away from each other we make a new df with all the years, to have the gaps.
  dfYears <- data.frame(YearAvg=(min(dfHEATsim$Year):max(dfHEATsim$Year)))
  
  dfHEATsim_avg<- dfYears %>%
    merge(dfHEATsim, all=T)
  
  #now we have all combinations possible, i.e. for 1933 we have a row for years from 1933 to 2023.
  #Of course, we do not need all this combinations,we are only interested in the year, and its previous 
  # 5 years (year + 4 previous years) so we can calculate the standard error for that year.
  dfHEATsim_avgx<- dfHEATsim_avg %>%
    filter(Year<=YearAvg, Year>=(YearAvg-4)) %>%
    group_by(YearAvg,A_group) %>%
    summarise(ERmean=mean(ER), SD=sd(ER), SE=sd(ER)/sqrt(n()), .groups="drop")
  
  dfHEATsim_avgx<- split(dfHEATsim_avgx,dfHEATsim_avgx$A_group)
  
  dfHEATsim_avgx <- dfHEATsim_avgx %>%
    map(~right_join(.x, dfYears, by = "YearAvg") %>%
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
  dfHEATmean5yr <- dfYears %>%
    merge(dfHEATmean, all=T)
  
  dfHEATmean5yr<- dfHEATmean5yr %>%
    filter(Year<=YearAvg, Year>=(YearAvg-4)) %>%
    group_by(YearAvg,A_group) %>%
    summarise(ERmean=mean(ER), SD=sd(ER), SE=sd(ER)/sqrt(n()), n=n(), .groups="drop")  
  
  dfHEATmean5yr<- split(dfHEATmean5yr,dfHEATmean5yr$A_group)
  
  dfHEATmean5yr <- dfHEATmean5yr %>%
    map(~right_join(.x, dfYears, by = "YearAvg") %>%
          mutate(A_group = unique(.x$A_group[!is.na(.x$A_group)]))
    )%>%
    bind_rows()%>%
    rename(Year=YearAvg)
  dfHEATmean5yr <- dfHEATmean5yr%>%
    left_join(areas,by="A_group")%>%mutate(Basin_name=factor(Basin_name,levels=area_f))
  
  dfHEATmean5yr <- dfHEATmean5yr%>%filter(Year>=1960)
  
  
  #arrange it 
  
  dfHEATmean5yrx<- dfHEATmean5yr%>%filter(A_group==1)
  dfHEATsim_avgy<- dfHEATsim_avgx%>%filter(A_group==1)
  dfHEATmeanx<-dfHEATmean%>%filter(A_group==1)
  bas<-unique(as.character(dfHEATmeanx$Basin_name))
  #we plot both, the simulation and the true values
  p1<-ggplot(dfHEATsim_avgy) + 
    geom_rect(aes(ymin = 0.8, ymax = Inf, xmin = -Inf, xmax = Inf), fill="#2C96F6",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.6, ymax = 0.8, xmin = -Inf, xmax = Inf), fill = "#8FFF07",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.4, ymax = 0.6, xmin = -Inf, xmax = Inf), fill = "#FFFF00",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.2, ymax = 0.4, xmin = -Inf, xmax = Inf), fill = "#FAB70E",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = -Inf, ymax = 0.2, xmin = -Inf, xmax = Inf), fill = "#FF2121",color=NA, alpha = 0.2) +
    labs(title =paste0(bas),y="EQR",x="Year")+
    geom_ribbon(aes(x=Year, ymin=ERmin, ymax=ERmax), alpha=0.2) +
    geom_line(data=dfHEATmean5yrx, aes(x=Year, y=ERmean), colour="black") +
    geom_point(data=dfHEATmeanx, aes(x=Year, y=ER),shape=1)+
    scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
    scale_x_continuous(limits = c(1960,2025),
                       labels= c("1960","1980","2000","2020"),
                       breaks = c(1960,1980,2000,2020))+
    theme_minimal()+
    theme(
      plot.background = element_rect(fill="white",colour = "white"),
      panel.grid = element_blank(),
      panel.spacing = unit(1.5,"lines"),
      panel.border = element_rect(colour = "grey",fill=NA),
      text = element_text(size=20,family = "Calibri"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_text(size=24),
      axis.ticks.x = element_line(),
      strip.text = element_text(size=20),
      strip.background = element_rect(colour = NA,linewidth = 6))+
    coord_cartesian(ylim=c(0,1),expand = T)
  
  dfHEATmean5yrx<- dfHEATmean5yr%>%filter(A_group==2)
  dfHEATsim_avgy<- dfHEATsim_avgx%>%filter(A_group==2)
  dfHEATmeanx<-dfHEATmean%>%filter(A_group==2)
  bas<-unique(as.character(dfHEATmeanx$Basin_name))
  #we plot both, the simulation and the true values
  p2<-ggplot(dfHEATsim_avgy) + 
    geom_rect(aes(ymin = 0.8, ymax = Inf, xmin = -Inf, xmax = Inf), fill="#2C96F6",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.6, ymax = 0.8, xmin = -Inf, xmax = Inf), fill = "#8FFF07",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.4, ymax = 0.6, xmin = -Inf, xmax = Inf), fill = "#FFFF00",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.2, ymax = 0.4, xmin = -Inf, xmax = Inf), fill = "#FAB70E",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = -Inf, ymax = 0.2, xmin = -Inf, xmax = Inf), fill = "#FF2121",color=NA, alpha = 0.2) +
    labs(title =paste0(bas),y="EQR",x="Year")+
    geom_ribbon(aes(x=Year, ymin=ERmin, ymax=ERmax), alpha=0.2) +
    geom_line(data=dfHEATmean5yrx, aes(x=Year, y=ERmean), colour="black") +
    geom_point(data=dfHEATmeanx, aes(x=Year, y=ER),shape=1)+
    scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
    scale_x_continuous(limits = c(1960,2025),
                       labels= c("1960","1980","2000","2020"),
                       breaks = c(1960,1980,2000,2020))+
    theme_minimal()+
    theme(
      plot.background = element_rect(fill="white",colour = "white"),
      panel.grid = element_blank(),
      panel.spacing = unit(1.5,"lines"),
      panel.border = element_rect(colour = "grey",fill=NA),
      text = element_text(size=20,family = "Calibri"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_text(size=24),
      axis.ticks.x = element_line(),
      strip.text = element_text(size=20),
      strip.background = element_rect(colour = NA,linewidth = 6))+
    coord_cartesian(ylim=c(0,1),expand = T)
  
  dfHEATmean5yrx<- dfHEATmean5yr%>%filter(A_group==3)
  dfHEATsim_avgy<- dfHEATsim_avgx%>%filter(A_group==3)
  dfHEATmeanx<-dfHEATmean%>%filter(A_group==3)
  bas<-unique(as.character(dfHEATmeanx$Basin_name))
  #we plot both, the simulation and the true values
  p3<-ggplot(dfHEATsim_avgy) + 
    geom_rect(aes(ymin = 0.8, ymax = Inf, xmin = -Inf, xmax = Inf), fill="#2C96F6",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.6, ymax = 0.8, xmin = -Inf, xmax = Inf), fill = "#8FFF07",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.4, ymax = 0.6, xmin = -Inf, xmax = Inf), fill = "#FFFF00",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.2, ymax = 0.4, xmin = -Inf, xmax = Inf), fill = "#FAB70E",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = -Inf, ymax = 0.2, xmin = -Inf, xmax = Inf), fill = "#FF2121",color=NA, alpha = 0.2) +
    labs(title =paste0(bas),y="EQR",x="Year")+
    geom_ribbon(aes(x=Year, ymin=ERmin, ymax=ERmax), alpha=0.2) +
    geom_line(data=dfHEATmean5yrx, aes(x=Year, y=ERmean), colour="black") +
    geom_point(data=dfHEATmeanx, aes(x=Year, y=ER),shape=1)+
    scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
    scale_x_continuous(limits = c(1960,2025),
                       labels= c("1960","1980","2000","2020"),
                       breaks = c(1960,1980,2000,2020))+
    theme_minimal()+
    theme(
      plot.background = element_rect(fill="white",colour = "white"),
      panel.grid = element_blank(),
      panel.spacing = unit(1.5,"lines"),
      panel.border = element_rect(colour = "grey",fill=NA),
      text = element_text(size=20,family = "Calibri"),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size=24,margin = margin(t=10,b=5,unit="mm")),
      axis.title.y = element_blank(),
      axis.text = element_text(size=24),
      axis.ticks.x = element_line(),
      strip.text = element_text(size=20),
      strip.background = element_rect(colour = NA,linewidth = 6))+
    coord_cartesian(ylim=c(0,1),expand = T)
  
  dfHEATmean5yrx<- dfHEATmean5yr%>%filter(A_group==4)
  dfHEATsim_avgy<- dfHEATsim_avgx%>%filter(A_group==4)
  dfHEATmeanx<-dfHEATmean%>%filter(A_group==4)
  bas<-unique(as.character(dfHEATmeanx$Basin_name))
  #we plot both, the simulation and the true values
  p4<-ggplot(dfHEATsim_avgy) + 
    geom_rect(aes(ymin = 0.8, ymax = Inf, xmin = -Inf, xmax = Inf), fill="#2C96F6",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.6, ymax = 0.8, xmin = -Inf, xmax = Inf), fill = "#8FFF07",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.4, ymax = 0.6, xmin = -Inf, xmax = Inf), fill = "#FFFF00",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.2, ymax = 0.4, xmin = -Inf, xmax = Inf), fill = "#FAB70E",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = -Inf, ymax = 0.2, xmin = -Inf, xmax = Inf), fill = "#FF2121",color=NA, alpha = 0.2) +
    labs(title =paste0(bas),y="EQR",x="Year")+
    geom_ribbon(aes(x=Year, ymin=ERmin, ymax=ERmax), alpha=0.2) +
    geom_line(data=dfHEATmean5yrx, aes(x=Year, y=ERmean), colour="black") +
    geom_point(data=dfHEATmeanx, aes(x=Year, y=ER),shape=1)+
    scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
    scale_x_continuous(limits = c(1960,2025),
                       labels= c("1960","1980","2000","2020"),
                       breaks = c(1960,1980,2000,2020))+
    theme_minimal()+
    theme(
      plot.background = element_rect(fill="white",colour = "white"),
      panel.grid = element_blank(),
      panel.spacing = unit(1.5,"lines"),
      panel.border = element_rect(colour = "grey",fill=NA),
      text = element_text(size=20,family = "Calibri"),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text(size=24,margin = margin(r=10,l=5,unit="mm")),
      axis.title.x = element_blank(),
      axis.text = element_text(size=24),
      axis.ticks.x = element_line(),
      strip.text = element_text(size=20),
      strip.background = element_rect(colour = NA,linewidth = 6))+
    coord_cartesian(ylim=c(0,1),expand = T)
  
  dfHEATmean5yrx<- dfHEATmean5yr%>%filter(A_group==5)
  dfHEATsim_avgy<- dfHEATsim_avgx%>%filter(A_group==5)
  dfHEATmeanx<-dfHEATmean%>%filter(A_group==5)
  bas<-unique(as.character(dfHEATmeanx$Basin_name))
  #we plot both, the simulation and the true values
  p5<-ggplot(dfHEATsim_avgy) + 
    geom_rect(aes(ymin = 0.8, ymax = Inf, xmin = -Inf, xmax = Inf), fill="#2C96F6",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.6, ymax = 0.8, xmin = -Inf, xmax = Inf), fill = "#8FFF07",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.4, ymax = 0.6, xmin = -Inf, xmax = Inf), fill = "#FFFF00",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.2, ymax = 0.4, xmin = -Inf, xmax = Inf), fill = "#FAB70E",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = -Inf, ymax = 0.2, xmin = -Inf, xmax = Inf), fill = "#FF2121",color=NA, alpha = 0.2) +
    labs(title =paste0(bas),y="EQR",x="Year")+
    geom_ribbon(aes(x=Year, ymin=ERmin, ymax=ERmax), alpha=0.2) +
    geom_line(data=dfHEATmean5yrx, aes(x=Year, y=ERmean), colour="black") +
    geom_point(data=dfHEATmeanx, aes(x=Year, y=ER),shape=1)+
    scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
    scale_x_continuous(limits = c(1960,2025),
                       labels= c("1960","1980","2000","2020"),
                       breaks = c(1960,1980,2000,2020))+
    theme_minimal()+
    theme(
      plot.background = element_rect(fill="white",colour = "white"),
      panel.grid = element_blank(),
      panel.spacing = unit(1.5,"lines"),
      panel.border = element_rect(colour = "grey",fill=NA),
      text = element_text(size=20,family = "Calibri"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_text(size=24),
      axis.ticks.x = element_line(),
      strip.text = element_text(size=20),
      strip.background = element_rect(colour = NA,linewidth = 6))+
    coord_cartesian(ylim=c(0,1),expand = T)
  
  dfHEATmean5yrx<- dfHEATmean5yr%>%filter(A_group==6)
  dfHEATsim_avgy<- dfHEATsim_avgx%>%filter(A_group==6)
  dfHEATmeanx<-dfHEATmean%>%filter(A_group==6)
  bas<-unique(as.character(dfHEATmeanx$Basin_name))
  #we plot both, the simulation and the true values
  p6<-ggplot(dfHEATsim_avgy) + 
    geom_rect(aes(ymin = 0.8, ymax = Inf, xmin = -Inf, xmax = Inf), fill="#2C96F6",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.6, ymax = 0.8, xmin = -Inf, xmax = Inf), fill = "#8FFF07",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.4, ymax = 0.6, xmin = -Inf, xmax = Inf), fill = "#FFFF00",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.2, ymax = 0.4, xmin = -Inf, xmax = Inf), fill = "#FAB70E",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = -Inf, ymax = 0.2, xmin = -Inf, xmax = Inf), fill = "#FF2121",color=NA, alpha = 0.2) +
    labs(title =paste0(bas),y="EQR",x="Year")+
    geom_ribbon(aes(x=Year, ymin=ERmin, ymax=ERmax), alpha=0.2) +
    geom_line(data=dfHEATmean5yrx, aes(x=Year, y=ERmean), colour="black") +
    geom_point(data=dfHEATmeanx, aes(x=Year, y=ER),shape=1)+
    scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
    scale_x_continuous(limits = c(1960,2025),
                       labels= c("1960","1980","2000","2020"),
                       breaks = c(1960,1980,2000,2020))+
    theme_minimal()+
    theme(
      plot.background = element_rect(fill="white",colour = "white"),
      panel.grid = element_blank(),
      panel.spacing = unit(1.5,"lines"),
      panel.border = element_rect(colour = "grey",fill=NA),
      text = element_text(size=20,family = "Calibri"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_text(size=24),
      axis.ticks.x = element_line(),
      strip.text = element_text(size=20),
      strip.background = element_rect(colour = NA,linewidth = 6))+
    coord_cartesian(ylim=c(0,1),expand = T)
  
  layout <- "
  #1#
  42#
  536
  "
  p<-p1+p2+p3+p4+p5+p6+plot_layout(design = layout,guides = "collect")+plot_annotation(theme = theme(
    axis.title.x = element_text(margin = margin(t = 1), vjust = 0),  # Add margin to x-axis title
    axis.title.y = element_text(margin = margin(r = 1)),
    panel.border = element_rect(colour = "grey",fill=NA),
    plot.background = element_rect(fill="white",colour = "darkgrey",linewidth=1.5)# Add margin to y-axis title
  )
  )
  
  ggsave(p, file=paste0(output_folder,"integrated_HEAT_basins.png"), height=30, width=40, units="cm")
  
  p<-ggplot(dfHEATsim_avgx) + 
    geom_rect(aes(ymin = 0.8, ymax = Inf,xmin = -Inf, xmax = Inf), fill="#2C96F6",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.6, ymax = 0.8, xmin = -Inf, xmax = Inf), fill = "#8FFF07",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.4, ymax = 0.6, xmin = -Inf, xmax = Inf), fill = "#FFFF00",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.2, ymax = 0.4, xmin = -Inf, xmax = Inf), fill = "#FAB70E",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = -Inf, ymax = 0.2, xmin = -Inf, xmax = Inf), fill = "#FF2121",color=NA, alpha = 0.2) +
    labs(y="EQR",x="Year")+
    geom_ribbon(aes(x=Year, ymin=ERmin, ymax=ERmax), alpha=0.2) +
    geom_line(data=dfHEATmean5yr, aes(x=Year, y=ERmean), colour="black") +
    geom_point(data=dfHEATmean, aes(x=Year, y=ER),shape=1)+
    facet_wrap(~Basin_name,nrow=2,scales="fixed")+
    scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
    scale_x_continuous(limits = c(1960,2025),
                       labels= c("1960","1980","2000","2020"),
                       breaks = c(1960,1980,2000,2020))+
    theme_minimal()+
    theme(
      plot.background = element_rect(fill="white",colour = "darkgrey",linewidth=1.5),
      panel.grid = element_blank(),
      panel.spacing = unit(1.5,"lines"),
      plot.margin = margin(r=20,l=10,t=10,b=10),
      panel.border = element_rect(colour = "grey",fill=NA),
      text = element_text(size=20,family = "Calibri"),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 8), vjust = 1),  # Add margin to x-axis title
      axis.title.y = element_text(margin = margin(r = -8)),
      axis.text = element_text(size=24),
      axis.ticks.x = element_line(),
      strip.text = element_text(size=20),
      strip.background = element_rect(colour = NA,linewidth = 6))
  ggsave(p, file=paste0(output_folder,"integrated_HEAT_basins_x.png"), height=25, width=40, units="cm")

  
# Fig 4 - Heatmap by basins ----
 
  # calculate EQR values and status class
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
  df <- df %>%
    CalcER()
  df1<-df%>%
    mutate(ER=EQR)
  
  #HEAT Heatmap
  # do aggregation and overall classification using HEAT
  dfHEAT <- aggregateHEAT(df1, var_cat="HEATCat",group_vars=c("A_group","Year"))
  
  dfHEAT <- dfHEAT %>% 
    left_join(areas,by="A_group")%>%
    mutate(Basin=factor(Basin_name,levels=area_f))%>%
    filter(Year>=1960)
  
  years<-1960:max(dfHEAT$Year,na.rm=T)
  
  dfHEATx <- tibble(Year=years) %>% 
    merge(tibble(Basin=area_f),all=T) %>%
    left_join(dfHEAT, by=c("Year","Basin"))
  
  dfHEATx$Year <- factor(dfHEATx$Year, levels=years)
  dfHEATx$Class <- factor(dfHEATx$HEATClass, levels = c("Bad","Poor", "Mod","Good","High"))
  dfHEATx$Basin <- factor(dfHEATx$Basin, levels= (area_f))
  
  
  fig1<-ggplot(dfHEATx, aes(Basin, Year, fill=Class)) +
    geom_tile(color="grey",show.legend = F) +
    labs(x="",y="",title="")+
    scale_fill_manual(values=color_palette,drop=F,breaks=StatusLevels,na.value = "white") +
    scale_x_discrete(position="top")+
    scale_y_discrete(expand=c(0,0),
                     breaks=c("2020","2010","2000","1990","1980","1970","1960"))+
    theme_bw() +
    theme(
      text= element_text(family = "Calibri"),
      axis.text.x=element_text(size=12,colour="#000000",angle=90,hjust=0,vjust=0.5),
      axis.text.y=element_text(size=12,vjust = 0.2,colour="#000000"),
      axis.ticks.y=element_line(linewidth=0.4),
      axis.ticks.x=element_blank(),
      title = element_text(size=12),
      panel.grid=element_blank(),
      legend.text = element_text(size=11),
      plot.title = element_text(hjust=0.5)
    ) +  
    geom_hline(yintercept = c("2020","2010","2000","1990","1980","1970"),color="darkgrey",linetype="dashed")
  
  fig1
  ggsave(fig1, file=paste0(output_folder,"heatmap HEAT basins.png"), height=22, width=7, units="cm")
  
  
  
# Fig 5 - Integrated assessment ----
  df <- data %>%
    filter(!(Parameter=="90_PCT_klfa"))%>%rename("value"="Avg") #we will be using the summer mean
  
  
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
  
  df <- df %>%
    CalcER()
  df1<-df%>%
    mutate(ER=EQR)
  
  #HEAT Heatmap
  
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
    select(A_group, Type, Year, HEATCat, Parameter, Indicator, mean=value, EQR,ER, Class)%>%
    filter(!is.na(A_group))
  
  #Using the function "random_indicator_values" we apply a Monte Carlo simulation. This simulation takes the actual value and its SE
  #and calculates 1000 random values that fit with the ones provided. 
  df_sim <- df1 %>%
    rowwise() %>%
    mutate(val_sim=list(random_indicator_values(value, SE, logAvg, logSE, Parameter))) %>%
    ungroup() %>%
    select(A_group, Type, Year, HEATCat, Parameter, Indicator, MatchValue, Match, Ref, HG, GM, MP, PB, Worst, Resp, val_sim)
  
  #this simulation gives you 1000 simualted values in a list, we need to unpack it. Therefore, for x value, now we have 1000 values.
  #the simulation is made by row, in another words, for each type of waterbody, for each year, and each parameter, we will get 1000 simulations.
  
  df_sim <- tidyr::unnest(df_sim, cols=val_sim) %>%
    rename(value=val_sim) %>%
    group_by(A_group, Type, Year, HEATCat, Parameter, Indicator, MatchValue, Match, Ref, HG, GM, MP, PB, Worst, Resp) %>%
    mutate(simID=row_number()) %>%
    ungroup()
  
  #now, from the simulated values, we calculate a simulated ER value.
  df_sim <- df_sim %>%
    CalcEQR()%>%
    mutate(ER=EQR)
  
  
  df_sim$Class <- factor(df_sim$Class, levels = c("Bad","Poor","Mod","Good","High"))
  
  
  #We use the aggregation function, aggregating by year and simID
  
  dfHEATsim <- aggregateHEAT(df_sim, var_cat="HEATCat",group_vars=c("Year","simID"))
  dfHEATmean <- aggregateHEAT(df_mean, var_cat="HEATCat",group_vars=c("Year"))
  
  
  #to avoid joining years too far away from each other we make a new df with all the years, to have the gaps.
  dfYears <- data.frame(YearAvg=(min(dfHEATsim$Year):max(dfHEATsim$Year)))
  
  dfHEATsim_avg<- dfYears %>%
    merge(dfHEATsim, all=T)
  
  #now we have all combinations possible, i.e. for 1933 we have a row for years from 1933 to 2023.
  #Of course, we do not need all this combinations,we are only intereseted in the year, and its previous 
  # 5 years (year + 4 previous years) so we can calculate the standard error for that year.
  dfHEATsim_avgx<- dfHEATsim_avg %>%
    filter(Year<=YearAvg, Year>=(YearAvg-4)) %>%
    group_by(YearAvg) %>%
    summarise(ERmean=mean(ER), SD=sd(ER), SE=sd(ER)/sqrt(n()), .groups="drop")
  
  dfHEATsim_avgx <- dfYears %>%
    left_join(dfHEATsim_avgx, by="YearAvg") %>%
    rename(Year=YearAvg)
  
  #to plot the bands for the confidence interval, we need a minimum and a max
  dfHEATsim_avgx <- dfHEATsim_avgx %>%
    mutate(ERmax=ERmean+SD, ERmin=ERmean-SD)
  
  dfHEATsim_avgx<-dfHEATsim_avgx%>%filter(Year>=1960)
  #we do the same process for the true values, so we can plot the simulated one and the original ones.
  dfHEATmean5yr <- dfYears %>%
    merge(dfHEATmean, all=T)
  
  dfHEATmean5yr<- dfHEATmean5yr %>%
    filter(Year<=YearAvg, Year>=(YearAvg-4)) %>%
    group_by(YearAvg) %>%
    summarise(ERmean=mean(ER), SD=sd(ER), SE=sd(ER)/sqrt(n()), n=n(), .groups="drop")  
  
  
  dfHEATmean5yr <- dfYears %>%
    left_join(dfHEATmean5yr, by="YearAvg") %>%
    rename(Year=YearAvg)
  dfHEATmean5yr<-dfHEATmean5yr%>%filter(Year>=1960)
  #we plot both, the simulation and the true values
  p<-ggplot(dfHEATsim_avgx) + 
    geom_rect(aes(ymin = 0.8, ymax = Inf, xmin = -Inf, xmax = Inf), fill="#2C96F6",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.6, ymax = 0.8, xmin = -Inf, xmax = Inf), fill = "#8FFF07",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.4, ymax = 0.6, xmin = -Inf, xmax = Inf), fill = "#FFFF00",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = 0.2, ymax = 0.4, xmin = -Inf, xmax = Inf), fill = "#FAB70E",color=NA, alpha = 0.2) +
    geom_rect(aes(ymin = -Inf, ymax = 0.2, xmin = -Inf, xmax = Inf), fill = "#FF2121",color=NA, alpha = 0.2) +
    labs(title =paste0("HEAT Integrated Assessment"),y="EQR",x="Year")+
    geom_ribbon(aes(x=Year, ymin=ERmin, ymax=ERmax), alpha=0.2) +
    geom_line(aes(x=Year, y=ERmean), linewidth=1) +
    geom_point(data=dfHEATmean, aes(x=Year, y=ER),shape=1)+
    scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
    scale_x_continuous(limits = c(1960,2024),
                       labels= c("1960","1970","1980","1990","2000","2010","2020"),
                       breaks = c(1960,1970,1980,1990,2000,2010,2020))+
    theme_minimal()+
    theme(plot.background = element_rect(fill="white",colour = "darkgrey",linewidth=1.5),
          panel.border = element_rect(colour = "grey",fill=NA),
          panel.grid = element_blank(),
          plot.margin = margin(r=20,l=10,t=10,b=10),
          text = element_text(size=20,family = "Calibri"),
          axis.title = element_text(size=20),
          axis.text = element_text(size=20),
          axis.ticks.x = element_line()) +
    coord_cartesian(ylim=c(0,1),expand = T)
  p
  ggsave(p, file=paste0(output_folder,"integrated HEAT.png"), height=15, width=25, units="cm")
  
# Fig 6 - Distance to target ----
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
    filter(!is.na(WB_group))
  
  # calculate HEAT eutrophication ratio - using the G/M boundary as threshold
  df <- df %>%
    CalcER()
  df1<-df%>%
    mutate(ER=EQR)
  
  
  df_mean <- df1 %>%
    select(A_group, Type, Year,QE, HEATCat, Parameter, Indicator, mean=value, ER, Class)%>%filter(!is.na(A_group))%>%
    mutate(Criteria= case_when(
      HEATCat =="Nut" ~ "Nutrients",
      HEATCat =="Pri" ~ "Direct effects",
      HEATCat =="Sec" ~ "Indirect effects",
      TRUE ~ NA
    ))
  
  df_mean<-df_mean%>%filter(Year>=1960)
  
  dfHEATmean <- aggregateHEAT(df_mean, var_cat="HEATCat",group_vars=c("Criteria","A_group"))
  
  
  dist_HEAT<-dfHEATmean%>%
    mutate(Distance= (ER-0.6)*-1)%>%
    left_join(areas,by="A_group")
  
  p<-ggplot(data=dist_HEAT, aes(x = Basin_name, y = Distance, fill = Criteria)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.4), width = 0.4,color="black") +
    scale_y_continuous(limits = c(-0.25, 0.5)) +
    scale_fill_manual(values = c("Nutrients"="white", "Direct effects"="grey70", "Indirect effects"="grey30")) +
    labs(x = "",
         y = "Distance to target") +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 0.5),
          axis.ticks.x = element_line(colour="black"),
          text= element_text(size=14,family="Calibri"),
          axis.text =element_text(size=14), 
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = c(0.8,0.8),
          legend.key.width = unit(22,"pt"),
          legend.key.height = unit(12,"pt"),
          legend.text = element_text(size=15),
          legend.box.background = element_rect(fill="white",color="black"),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white",colour = "darkgrey",linewidth=1.5),
          panel.border = element_rect(colour = "grey",fill=NA),
          plot.margin = margin(r=20,l=10,t=10,b=10),)
  p
  
  ggsave(p, file=paste0(output_folder,"distance target basins.png"), height=17, width=27, units="cm")
  