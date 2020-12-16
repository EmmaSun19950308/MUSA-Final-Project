setwd("~/Box Sync/GitHub/Final Project")


# 1. Load library
library(tidyverse)
library(sf)
library(QuantPsyc)
library(RSocrata)
library(viridis)
library(caret)
library(spatstat)
library(spdep)
library(FNN)
library(grid)
library(gridExtra)
library(knitr)
library(kableExtra)
library(tidycensus)
library(nnet)
library(reshape2)
library(class)
library(MASS)
library(dplyr)
library(tidyr)
library(lubridate)
library(plotly)
library(ggmap)
library(data.table)
library(gganimate)


# 2. Identify functions
mapTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 15,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    strip.text.x = element_text(size = 14))
}

plotTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(color = "#0277a2", size=15, face="bold"),
    plot.subtitle = element_text(face="italic"),
    plot.caption = element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_line("#E5E5E5", size = 0.1),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    strip.background = element_rect(fill = "#E5E5E5", color = "white"),
    strip.text = element_text(size=12),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(colour = "black", face = "italic"),
    legend.text = element_text(colour = "black", face = "italic"),
    strip.text.x = element_text(size = 14)
  )
}

#3. Load Quantile break functions

qBr <- function(df, variable, rnd) {
  if (missing(rnd)) {
    as.character(quantile(round(df[[variable]],0),
                          c(.01,.2,.4,.6,.8), na.rm=T))
  } else if (rnd == FALSE | rnd == F) {
    as.character(formatC(quantile(df[[variable]]), digits = 3),
                 c(.01,.2,.4,.6,.8), na.rm=T)
  }
}

q5 <- function(variable) {as.factor(ntile(variable, 5))}

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

# for calculating average nearest neighbor distance.

nn_function <- function(measureFrom,measureTo,k) {
  measureFrom_Matrix <- as.matrix(measureFrom)
  measureTo_Matrix <- as.matrix(measureTo)
  nn <-   
    get.knnx(measureTo, measureFrom, k)$nn.dist
  output <-
    as.data.frame(nn) %>%
    rownames_to_column(var = "thisPoint") %>%
    gather(points, point_distance, V1:ncol(.)) %>%
    arrange(as.numeric(thisPoint)) %>%
    group_by(thisPoint) %>%
    dplyr::summarize(pointDistance = mean(point_distance)) %>%
    arrange(as.numeric(thisPoint)) %>% 
    dplyr::select(-thisPoint) %>%
    pull() # pull() is similar to $. It's mostly useful because it looks a little nicer in pipes, it also works with remote data frames, and it can optionally name the output.
  
  return(output)  
}


palette7 <- c('#4c281a','#8f4405',"#ff6d69","#0be7fb","#3182bd","#08519c",'#f9c37a')
palette5 <- c("#eff3ff","#bdd7e7","#6baed6","#3182bd","#08519c")
palette4 <- c("#D2FBD4","#FFA500","#005582",'#B22222')
palette3 <- c("#41c447","#005582", '#B22222')
palette2 <- c("#eadd46","#485123")



# DATA

## 1. train data
train = read.csv("trains_train.csv") 

## 2. station data
stations0 = read.csv("stations.csv")

stations.sf <- 
  stations0 %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant") %>%
  st_transform(31370)
belgium_basemap <- get_stamenmap(bbox = c(left = 2.739258, bottom = 49.563228,
                                          right = 6.551513, top = 51.695260), zoom = 11)

### split URI for the station codes
stations = stations.sf %>%
  separate(URI, c("res", "station_name"), "NMBS/")

### drop unneccesary columns
drops <- c("res","alternative.fr","alternative.nl","alternative.de", "alternative.en")
stations <- stations[ , !(names(stations) %in% drops)]

### Check station distribution
ggmap(belgium_basemap)+
  geom_point(aes(x = longitude, y = latitude), data = stations0, size = 1, color ='darkred') +
  labs(title="Train Station Distribution in Belgium\n",
       caption = 'Figure 1',
       x="Longtitude", 
       y="Latitude")+
  plotTheme()


## 3. Weather data
W701 <- read.csv("weather_data_july_1.csv")
W702 <- read.csv("weather_data_july_2.csv")
W801 <- read.csv("weather_data_aug_1.csv")
W802 <- read.csv("weather_data_aug_2.csv")
W901 <- read.csv("weather_data_sep_1.csv")
W902 <- read.csv("weather_data_sep_2.csv")
W101 <- read.csv("weather_data_oct_1.csv")
W102 <- read.csv("weather_data_oct_2.csv")

weather <- rbind(W701,W702,W801,W802,W901,W902,W101,W102)

weather <- weather %>%
  mutate(interval60 = ymd_hms(date_time))


grid.arrange(
  ggplot(trn_j_w, aes(interval60,humidity)) + geom_line() +
    labs(title="Percipitation", x="Hour", y="Perecipitation") + 
    plotTheme(),
  ggplot(trn_j_w, aes(interval60,windspeed)) + geom_line() +
    labs(title="Wind Speed", x="Hour", y="Wind Speed") + 
    plotTheme(),
  ggplot(trn_j_w, aes(interval60,temperature)) + geom_line() +
    labs(title="Temperature", x="Hour", y="temperature") + plotTheme(),
  top="Weather Data - Belgium - May, 2018")




## 4. Demographic data: population density
### population
### https://ec.europa.eu/eurostat/cache/metadata/en/demo_r_gind3_esms.htm
### topo data
### https://gisco-services.ec.europa.eu/distribution/v1/nuts-2016.html
bel <- st_read("NUTS_RG_60M_2016_4326_LEVL_3.shp") 
popdens <- read.csv("demo_r_d3dens_1_Data.csv")
belpopden <- bel %>%
  left_join(dplyr::select(popdens, c(GEO, Value)), by = c("FID" = "GEO")) %>%
  st_transform(31370)



## 5. Facility data
facility <- read_csv("https://raw.githubusercontent.com/iRail/stations/master/facilities.csv")

facility[is.na(facility)] <- 0
facility <- 
  facility %>% 
  mutate(facility_sum = select(.,ticket_vending_machine:audio_induction_loop) %>% 
           rowSums())

facility <- 
  facility %>% 
  select(name,facility_sum)

stations <- merge(stations,facility,by="name",all=T) 

## 6. Railway data
### https://mapcruzin.com/free-belgium-arcgis-maps-shapefiles.htm
rail <- st_read("belgium-railways-shape/railways.shp") %>%
  st_transform(4326)

dat <- read.table(textConnection("start.date start.time end.date end.time
2012-07-13   15:01:32 2012-07-13 15:02:42
2012-07-05   18:26:31 2012-07-05 18:27:19 
2012-07-14   20:23:21 2012-07-14 20:24:11"), header=TRUE)
starttime = hm(facility$sales_open_monday)
endtime = hm(facility$sales_close_monday)
interval = difftime(endtime,starttime,units = "hours")

facility$sales_close_monday <- 
ifelse(facility$sales_close_monday %in% c(""," ","NA"), 
       "00:00", facility$sales_close_monday)
facility$sales_open_monday <- 
  ifelse(facility$sales_open_monday %in% c(""," ","NA"), 
         "00:00", facility$sales_open_monday)

# Final Dataset

trn_j <-
  read.csv("trains_train.csv") %>%
  mutate(from = as.character(from), to = as.character(to)) %>%
  left_join(dplyr::select(stations, station_name), by = c("from" = "station_name")) %>%
  st_sf() %>% 
  mutate(from.X = st_coordinates(.)[,1], from.Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  left_join(dplyr::select(stations, station_name), by = c("to" = "station_name")) %>%
  st_sf() %>% 
  mutate(to.X = st_coordinates(.)[,1], to.Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  mutate(distance = sqrt((from.X - to.X)^2 + (from.Y - to.Y)^2)) %>%
  arrange(-distance)


trn_j <- trn_j %>%
  mutate(datetime <- paste(date,time),
         datetime = parse_date_time(datetime, '%y-%m-%d %I:%M:%S %p'),
         interval60 = floor_date(ymd_hms(datetime), unit = "hour"),
         interval15 = floor_date(ymd_hms(datetime), unit = "15 mins"),
         week = week(interval60),
         dotw = wday(interval60,label = TRUE, abbr = TRUE))
trn_j$`datetime <- paste(date, time)` <- NULL

trn_j$veh <- gsub("[^a-zA-Z]", "", trn_j$vehicle) 
trn_j$veh[trn_j$veh==""] <- "undefined"
trn_j$veh[trn_j$veh=="ic"] <- "ICE"


trn_j_f <- merge(trn_j, stations,  by.x = "from", by.y = "station_name") %>%
  st_as_sf()

## add weather data 
trn_j_w <- merge(trn_j_f,weather,by.x=c("name","interval60"),by.y=c("station_name","interval60"))




# III. EDA 

## EDA on some features


grid.arrange(
  ggplot(trn_j_w, aes(interval60,humidity)) + geom_line() +
    labs(title="Percipitation", x="Hour", y="Perecipitation") + 
    plotTheme(),
  ggplot(trn_j_w, aes(interval60,windspeed)) + geom_line() +
    labs(title="Wind Speed", x="Hour", y="Wind Speed") + 
    plotTheme(),
  ggplot(trn_j_w, aes(interval60,temperature)) + geom_line() +
    labs(title="Temperature", x="Hour", y="temperature") + plotTheme(),
  top="Weather Data - Belgium - May, 2018")

## Temporal Trend

### count by week
ggplot(trn_j_w %>% 
         mutate(hour = hour(datetime)) %>% 
         filter(occupancy == "high" | occupancy == "medium" | occupancy == "low"), 
       aes(hour, color = dotw))+
  geom_freqpoly(binwidth = 1)+
  facet_wrap(~occupancy, ncol=3) +
  scale_colour_manual(values = palette7) +
  labs(title="Train Shifts in Belgium, by day of the week, 2016\n",
       x="Hour", 
       y="Trip Counts") +
  plotTheme()

### count by day
ggplot(trn_j_w %>% 
         mutate(hour = hour(datetime)))+
  geom_freqpoly(aes(hour, color = occupancy), binwidth = 1)+
  scale_colour_manual(values = palette3) +
  labs(title="Number of Train Shifts in Belgium, by occupancy, 2016\n",
       x="Hour", 
       y="Trip Counts")+
  plotTheme()

### count by weekday vs weekend
ggplot(trn_j_w %>% mutate(hour = hour(datetime),
                             weekend = ifelse(dotw %in% c("7", "6"), "Weekend", "Weekday"))) +
  geom_freqpoly(aes(hour, color = weekend), binwidth = 1)+
  scale_color_manual(values = palette2,
                     name = 'Period') +
  labs(title="Bike share trips in DC - weekend vs weekday, Week 35 - Week 39, 2020\n",
       caption = 'Figure ',
       x="Hour", 
       y="Trip Counts")+
  plotTheme()

### count by weather: visibility
ggplot(trn_j_w %>% 
         mutate(hour = hour(datetime)) %>% 
         filter(occupancy == "high" | occupancy == "medium" | occupancy == "low"), 
       aes(visibility, color = dotw))+
  geom_freqpoly(binwidth = 1)+
  facet_wrap(~occupancy,ncol=3) +
  scale_colour_manual(values = palette3) +
  labs(title="Train Shifts in Belgium, by day of the week, 2016\n",
       caption = 'Figure ',
       x="Hour", 
       y="Trip Counts") +
  plotTheme()

### count by weather: temperature
ggplot(trn_j_w %>% 
         mutate(hour = hour(datetime)) %>% 
         filter(occupancy == "high" | occupancy == "medium" | occupancy == "low"), 
       aes(temperature, color = dotw))+
  geom_freqpoly(binwidth = 1)+
  facet_wrap(~occupancy,ncol=3) +
  scale_colour_manual(values = palette7) +
  labs(title="Train Shifts in Belgium, by day of the week, 2016\n",
       caption = "Figure ",
       x="Hour", 
       y="Trip Counts") +
  plotTheme()


# Spatial Trend
sort(summary(trn_j$vehicle,maxsum=1500),decreasing = TRUE)
freq_v = c("IC1518","IC429","IC1515","IC407","P7305","IC1807","IC3631","1828","8015","IC1831","IC539","P7444","S83978", "IC1507","IC716","L557","S23665","IC3432","IC4317","S53586")

## occupancy by lines: origin station
p1 <- ggplot() +
  geom_sf(data = rail, aes(color = "grey")) + 
  geom_sf(data=subset(trn_j_f, trn_j_f$vehicle %in% freq_v),aes(colour = vehicle),size=2,show.legend = "point")+ 
  labs(title="Top 20 Frequent Lines Plotted to the Departure\n",
       caption = "Figure")+
  mapTheme()+
  plotTheme() +
  theme(legend.position = "bottom") + 
  transition_manual(factor(vehicle, levels = freq_v), cumulative = TRUE) 
 
gganimate::animate(p1,  duration=10,renderer = gifski_renderer())

## occupancy by lines: destination
trn_j_t <- merge(trn_j, stations,  by.x = "to", by.y = "station_name") %>%
  st_as_sf()

p2 <- ggplot()+
  geom_sf(data=rail,aes(color = "grey"))+
  geom_sf(data=subset(trn_j_t, trn_j_t$vehicle %in% freq_v),aes(colour = vehicle),size=2,show.legend = "point")+ 
  labs(title="Top 20 Frequent Lines Plotted to Destination\n",
       caption = 'Figure')+
  mapTheme()+
  plotTheme() +
  theme(legend.position = "bottom") + 
  transition_manual(factor(vehicle, levels = freq_v), cumulative = TRUE)

gganimate::animate(p2,  duration=10,renderer = gifski_renderer())

## occupancy by week days: departure

stations.df <- as.data.frame(stations.sf)
trn_j <- trn_j %>% 
  mutate(occ_num=as.numeric(case_when(
    trn_j$occupancy=="high" ~ 3,
    trn_j$occupancy=="medium" ~ 2,
    trn_j$occupancy=="low" ~ 1,
  )))
trn_j_from <- setDT(trn_j)
trn_j_from <- trn_j_from[,list(from_occ=mean(occ_num)),by=list(from,dotw)] %>% 
  as.data.frame() 
colnames(trn_j_from)[colnames(trn_j_from)=="dotw"]<- "from_dotw"
stations_from <-merge(stations,trn_j_from, by.x="station_name", by.y="from") %>%
  mutate(from_dotw_text=case_when(
    from_dotw==1~"Mon",
    from_dotw==2~"Tue",
    from_dotw==3~"Wed",
    from_dotw==4~"Thu",
    from_dotw==5~"Fri",
    from_dotw==6~"Sat",
    from_dotw==7~"Sun"
  ))

ggplot() +
  geom_sf(data = rail, color = "grey") +  
  geom_sf(data = st_sf(stations_from[stations_from$country.code=='be',]), 
          aes(color = from_occ) ,size = 0.5)+
  labs(title="Occupancy Rate of Departure Station by Week Days\n",
       caption = 'Figure')+
  scale_colour_viridis(direction = -1,
                       discrete = FALSE, option = "plasma",
                       name = 'Occupancy')+
  facet_wrap(~from_dotw_text,nrow =2) +
  mapTheme()+
  plotTheme() 


## Occupancy by week days: Destination

trn_j_to <- setDT(trn_j)
trn_j_to <- trn_j_to[,list(to_occ=mean(occ_num)),by=list(to,dotw)] %>% 
  as.data.frame() 
colnames(trn_j_from)[colnames(trn_j_from)=="dotw"]<- "to_dotw"

stations_to <-merge(stations,trn_j_to,by.x="station_name",by.y="to")

stations_to <-merge(stations,trn_j_to, by.x="station_name", by.y="to") %>%
  mutate(to_dotw_text=case_when(
    dotw==1~"Mon",
    dotw==2~"Tue",
    dotw==3~"Wed",
    dotw==4~"Thu",
    dotw==5~"Fri",
    dotw==6~"Sat",
    dotw==7~"Sun"
  ))


ggplot() +
  geom_sf(data = rail, color = "grey") +  
  geom_sf(data = st_sf(stations_to[stations_to$country.code=='be',]), 
          aes(color = to_occ) ,size = 0.5)+
  labs(title="Occupancy Rate of Destination Station by Week Days\n",
       caption = 'Figure')+ 
  scale_colour_viridis(direction = -1,
                       discrete = FALSE, option = "plasma",
                       name = 'Occupancy')+
  
  facet_wrap(~to_dotw_text,nrow =2) +
  mapTheme()+
  plotTheme() 


## Feature Engineering 

### time lag 
trn_j <- 
  trn_j %>% 
  arrange(from, interval60) %>% 
  mutate(lagHour = dplyr::lag(occ_num,1),
         lag2Hours = dplyr::lag(occ_num,2),
         lag3Hours = dplyr::lag(occ_num,3),
         lag4Hours = dplyr::lag(occ_num,4),
         lag12Hours = dplyr::lag(occ_num,12))

as.data.frame(trn_j) %>%
  group_by(interval60) %>% 
  summarise_at(c(vars(starts_with("lag")), "occ_num"), mean, na.rm = TRUE) %>%
  gather(Variable, Value, -interval60, -occ_num) %>%
  mutate(Variable = factor(Variable, levels=c("lagHour","lag2Hours","lag3Hours","lag4Hours",
                                              "lag12Hours")))%>%
  group_by(Variable) %>%  
  summarize(correlation = round(cor(Value, occ_num),2))




