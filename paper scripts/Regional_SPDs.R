#load packages
library(rcarbon)
library(tidyverse)
library(ggpubr)
library(sf)


#brute force stops scientific numbers in favour of full digits.
options(scipen=999)

#load the c14 dates
mydates = read.csv("Master Paper/csv/arabia_dates_20.09.24.csv")


##Add the subregion filtering here!
#read the study area polygon
subregions= read_sf("Master Paper/shp/Arabia_subregions.gpkg")

#subset the levant region
southwest = subregions[1,1]
southeast = subregions[2,1]
bahrain = subregions[3,1]
northwest = subregions[4,1]


####Turn the sites into a simple feature with geometry.####
# provide the EPSG code using sf::st_crs()
crs_wgs84 = st_crs(4326) # WGS84 has EPSG code 4326

#provide the coordinates
sites_sf = st_as_sf(mydates, coords = c("Longitude", "Latitude"), crs = crs_wgs84)
st_crs(southwest)= crs_wgs84
st_crs(southeast)= crs_wgs84
st_crs(bahrain)= crs_wgs84
st_crs(northwest)= crs_wgs84


# provide the EPSG code using st_crs() - https://epsg.io/32637
UTM_Zone_37N = st_crs(32637) # IGRS_UTM_Zone_38N EPSG: 3891

# transform sites from crs 4326 (WGS84) to 32637(WGS84/UTM zone 37N) so that it is more relevant to the study area.
sites_sf = st_transform(sites_sf, UTM_Zone_37N)
southwest = st_transform(southwest, UTM_Zone_37N)
southeast = st_transform(southeast, UTM_Zone_37N)
bahrain = st_transform(bahrain, UTM_Zone_37N)
northwest = st_transform(northwest, UTM_Zone_37N)

####turn them into polygons####
sw.polygon = st_geometry(obj = southwest)
se.polygon = st_geometry(obj = southeast)
b.polygon = st_geometry(obj = bahrain)
nw.polygon = st_geometry(obj = northwest)

#subset the sites within the relevant area
sites_sf$sw.intersect = lengths(st_intersects(sites_sf$geometry, sw.polygon, sparse = TRUE, prepared = FALSE))
sites_sf$se.intersect = lengths(st_intersects(sites_sf$geometry, se.polygon, sparse = TRUE, prepared = FALSE))
sites_sf$b.intersect = lengths(st_intersects(sites_sf$geometry, b.polygon, sparse = TRUE, prepared = FALSE))
sites_sf$nw.intersect = lengths(st_intersects(sites_sf$geometry, nw.polygon, sparse = TRUE, prepared = FALSE))

#subset just the sites in southwest
sw.dates = sites_sf %>% 
  group_by(sw.intersect) %>%
  filter(sw.intersect >= 1)
#subset just the sites in southeast
se.dates = sites_sf %>% 
  group_by(se.intersect) %>%
  filter(se.intersect >= 1)
#subset just the sites in bahrain
b.dates = sites_sf %>% 
  group_by(b.intersect) %>%
  filter(b.intersect >= 1)
#subset just the sites in northwest
nw.dates = sites_sf %>% 
  group_by(nw.intersect) %>%
  filter(nw.intersect >= 1)


#count the number of dates
n.sw.dates = sw.dates %>% 
  filter(CRA<=5200&CRA>=3200)
n.se.dates = se.dates %>% 
  filter(CRA<=5200&CRA>=3200)
n.b.dates = b.dates %>% 
  filter(CRA<=5200&CRA>=3200)
n.nw.dates = nw.dates %>% 
  filter(CRA<=5200&CRA>=3200)

######SPDs######
###bin prep
#define bin clustering for rcarbon
bins.n = 50
#smoothing of SPDs
runm <- 50

#regional bins
#define the bin object for use below
sw.bins = binPrep(sites=sw.dates$SiteName, ages=sw.dates$CRA, h = bins.n)
se.bins = binPrep(sites=se.dates$SiteName, ages=se.dates$CRA, h = bins.n)
b.bins = binPrep(sites=b.dates$SiteName, ages=b.dates$CRA, h = bins.n)
nw.bins = binPrep(sites=nw.dates$SiteName, ages=nw.dates$CRA, h = bins.n)

#calibrate the southwest dates
sw.caldates=calibrate(x=sw.dates$CRA, errors=sw.dates$Error, calCurves='intcal20', normalised=TRUE, spdnormalised = TRUE, type = 'full')
#run the southwest spd
sw.spd = spd(sw.caldates,bins = sw.bins, runm=runm, timeRange=c(5200,3200), spdnormalised = TRUE)
#extract the values into a tidyverse dataframe and add the BC dates
extracted_sw.spd = as.data.frame(sw.spd$grid) %>% 
  mutate(time = (calBP*-1)+1950)


#calibrate the southeast dates
se.caldates=calibrate(x=se.dates$CRA, errors=se.dates$Error, calCurves='intcal20', normalised=TRUE, spdnormalised = TRUE, type = 'full')
#run the southeast spd
se.spd = spd(se.caldates,bins = se.bins, runm=runm, timeRange=c(5200,3200), spdnormalised = TRUE)
#extract the values into a tidyverse dataframe and add the BC dates
extracted_se.spd = as.data.frame(se.spd$grid) %>% 
  mutate(time = (calBP*-1)+1950)

#calibrate the bahrain dates
b.caldates=calibrate(x=b.dates$CRA, errors=b.dates$Error, calCurves='intcal20', normalised=TRUE, spdnormalised = TRUE, type = 'full')
#run the bahrain spd
b.spd = spd(b.caldates,bins = b.bins, runm=runm, timeRange=c(5200,3200), spdnormalised = TRUE)
#extract the values into a tidyverse dataframe and add the BC dates
extracted_b.spd = as.data.frame(b.spd$grid) %>% 
  mutate(time = (calBP*-1)+1950)

#calibrate the northwest dates
nw.caldates=calibrate(x=nw.dates$CRA, errors=nw.dates$Error, calCurves='intcal20', normalised=TRUE, spdnormalised = TRUE, type = 'full')
#run the northwest spd
nw.spd = spd(nw.caldates,bins = nw.bins, runm=runm, timeRange=c(5200,3200), spdnormalised = TRUE)
#extract the values into a tidyverse dataframe and add the BC dates
extracted_nw.spd = as.data.frame(nw.spd$grid) %>% 
  mutate(time = (calBP*-1)+1950)


joined_regionals = extracted_sw.spd %>% 
  left_join(extracted_se.spd, join_by(calBP == calBP, time == time))%>% 
  left_join(extracted_b.spd, join_by(calBP == calBP, time == time))%>% 
  left_join(extracted_nw.spd, join_by(calBP == calBP, time == time)) %>% 
  mutate(Southwest = PrDens.x,
         PrDens.x = NULL,
         Southeast = PrDens.y,
         PrDens.y = NULL,
         Bahrain = PrDens.x.x,
         PrDens.x.x = NULL,
         Northwest = PrDens.y.y,
         PrDens.y.y = NULL)
  

#####Permutation
#merge the regional dates for permutation:
merged_regional = bind_rows(sw.dates, se.dates, b.dates, nw.dates) %>% 
  pivot_longer(cols = c(sw.intersect,se.intersect,b.intersect,nw.intersect), names_to = "region") %>% 
  filter(value>0)

#prep the permutation
cal.all.dates = calibrate(merged_regional$CRA,merged_regional$Error,normalised=TRUE, spdnormalised = TRUE)
bins.all.dates = binPrep(ages = merged_regional$CRA,sites = merged_regional$SiteName,h=50)
#run it
perm.all.dates=permTest(x=cal.all.dates,marks=merged_regional$region,timeRange=c(5200,3200),backsight = 50,bins=bins.all.dates,nsim=10000,runm=50, spdnormalised = TRUE)
round(perm.all.dates$pValueList,4) #extract p-values
#look at the summary of the permutations to see significant deviations
summary(perm.all.dates)

par(mfrow=c(4,1))
plot(perm.all.dates,focalm = 1,main="Southern Levant")
plot(perm.all.dates,focalm = 2,main="Northern Levant/Upper Mesopotamia")
plot(perm.all.dates,focalm = 3,main="South-Central Anatolia")
plot(perm.all.dates,focalm = 4,main="South-Central Anatolia")

#extract for tidyverse compatability
extracted_perm = as.data.frame(perm.all.dates$envelope) %>% 
  mutate(calBP = 5200:3200) %>% 
  left_join(joined_regionals)

####a bit of wrangling to make sure the ribbon behaves properly for red and blue colouration
extracted_perm.plot = extracted_perm %>% 
  mutate(sw.above = ifelse((Southwest>=sw.intersect.2), Southwest, NA_real_), sw.below = ifelse((Southwest<=sw.intersect.1), Southwest, NA_real_)) %>% 
  mutate(se.above = ifelse((Southeast>=se.intersect.2), Southeast, NA_real_), se.below = ifelse((Southeast<=se.intersect.1), Southeast, NA_real_)) %>% 
  mutate(b.above = ifelse((Bahrain>=b.intersect.2), Bahrain, NA_real_), b.below = ifelse((Bahrain<=b.intersect.1), Bahrain, NA_real_)) %>% 
  mutate(nw.above = ifelse((Northwest>=nw.intersect.2), Northwest, NA_real_), nw.below = ifelse((Northwest<=nw.intersect.1), Northwest, NA_real_))


##plot all of the figures now
Figure_xa = ggplot(extracted_perm.plot, aes(x = calBP)) +
  geom_ribbon(aes(ymin=sw.intersect.1, ymax=sw.intersect.2, fill = "Envelope"), alpha = 0.5)+
  geom_ribbon(aes(ymin=sw.intersect.2, ymax=sw.above, fill = ">95%"), alpha = 0.5)+
  geom_ribbon(aes(ymin=sw.intersect.1, ymax=sw.below, fill = "<5%"), alpha = 0.5)+
  geom_line(aes(y = Southwest, colour = "Southwest SPD"), size=1, show.legend = FALSE)+
  scale_fill_manual(values=c("Envelope"="grey75", ">95%" = "blue", "<5%" = "red"))+
  scale_colour_manual(values=c("Southwest SPD"="#F4C431"))+
  scale_y_continuous(limits = c(0,0.0015), breaks=seq(0,0.0015, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  #plot the number of dates
  annotate("text", x=5200,y=0.0008, label = paste0("n = ",nrow(n.sw.dates)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  #plot the p-value
  annotate("text", x=5200,y=0.0007, label = paste0("global p-value = ",round((perm.all.dates$pValueList[1]), 5)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  #Southwestern
  annotate("text", x=5175,y=0.001475, label = "Southwest", color = "navyblue", 
           angle = 0, hjust = 0.0, size = 4, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, y = 0.001415, xend = 5200, yend = 0.001415, color = "navyblue", linetype = "dotted",linewidth = .5)+
  annotate("segment", x = 3200, y = 0.001415, xend = 5150, yend = 0.001415, color = "navyblue", linetype = "solid",linewidth = .5)+
  #Neolithic end
  annotate("text", x=5025,y=0.00136, label = "Neolithic", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4850, xend = 4850, y = 0.001415, yend = 0.001345,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Bronze age end
  annotate("text", x=4025,y=0.00136, label = "Bronze Age", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3200, xend = 3200, y = 0.001415, yend = 0.001345,  color = "navyblue",linetype="solid",linewidth = .5) +
  ##
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.25, 'cm'),legend.position = c(.11, .75),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.4,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

Figure_xb = ggplot(extracted_perm.plot, aes(x = calBP)) +
  geom_ribbon(aes(ymin=se.intersect.1, ymax=se.intersect.2, fill = "Envelope"), alpha = 0.5)+
  geom_ribbon(aes(ymin=se.intersect.2, ymax=se.above, fill = ">95%"), alpha = 0.5)+
  geom_ribbon(aes(ymin=se.intersect.1, ymax=se.below, fill = "<5%"), alpha = 0.5)+
  geom_line(aes(y = Southeast, colour = "Southeast SPD"), size=1, show.legend = FALSE)+
  scale_fill_manual(values=c("Envelope"="grey75", ">95%" = "blue", "<5%" = "red"))+
  scale_colour_manual(values=c("Southeast SPD"="#F4C431"))+
  scale_y_continuous(limits = c(0,0.0014), breaks=seq(0,0.0014, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  #plot the number of dates
  annotate("text", x=5200,y=0.00075, label = paste0("n = ",nrow(n.se.dates)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  #plot the p-value
  annotate("text", x=5200,y=0.00065, label = paste0("global p-value = ",round((perm.all.dates$pValueList[2]), 5)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  #Southeastern
  annotate("text", x=5175,y=0.00135, label = "Southeast", color = "navyblue", 
           angle = 0, hjust = 0.0, size = 4, fontface = "bold", family = "sans")+
  annotate("segment", x = 3200, y = 0.00129, xend = 5150, yend = 0.00129, color = "navyblue", linetype = "solid",linewidth = .5)+
  annotate("segment", x = 5150, y = 0.00129, xend = 5200, yend = 0.00129, color = "navyblue", linetype = "dotted",linewidth = .5)+
  #Hafit start
  annotate("text", x=4900,y=0.00124, label = "Hafit", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, xend = 5150, y = 0.00129, yend = 0.00123,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Hafit end
  annotate("segment", x = 4650, xend = 4650, y = 0.00129, yend = 0.00123, color = "navyblue",linetype="solid",linewidth = .5) +
  #Umm an-Nar start
  annotate("text", x=4300,y=0.00124, label = "Umm an-Nar", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.00129, yend = 0.00123,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Umm an-Nar end
  annotate("segment", x = 3950, xend = 3950, y = 0.00129, yend = 0.00123,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Wadi Suq start
  annotate("text", x=3775,y=0.00124, label = "Wadi Suq", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.00129, yend = 0.00123,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Wadi Suq end
  annotate("segment", x = 3600, xend = 3600, y = 0.00129, yend = 0.00123,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Late Bronze start
  annotate("text", x=3400,y=0.00124, label = "Late Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3600, xend = 3600, y = 0.00129, yend = 0.00123,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Late Bronze start
  annotate("segment", x = 3200, xend = 3200, y = 0.00129, yend = 0.00123,  color = "navyblue",linetype="solid",linewidth = .5) +
  ##
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.25, 'cm'),legend.position = c(.11, .75),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.4,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

Figure_xc = ggplot(extracted_perm.plot, aes(x = calBP)) +
  geom_line(aes(y = Bahrain, colour = "Bahrain SPD"), size=1, show.legend = FALSE)+
  scale_colour_manual(values=c("Bahrain SPD"="#F4C431"))+
  scale_y_continuous(limits = c(0,0.0024), breaks=seq(0,0.0024, .0006))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  #plot the number of dates
  annotate("text", x=5200,y=0.001, label = paste0("n = ",nrow(n.b.dates)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  #plot the p-value
  annotate("text", x=5200,y=0.000875, label = paste0("global p-value = ",round((perm.all.dates$pValueList[3]), 5)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  #Bahrain
  annotate("text", x=5175,y=0.0023, label = "Bahrain", color = "navyblue", 
           angle = 0, hjust = 0.0, size = 4, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, y = 0.0022, xend = 5200, yend = 0.0022, color = "navyblue", linetype = "dotted",linewidth = .5)+
  annotate("segment", x = 3250, y = 0.0022, xend = 5150, yend = 0.0022, color = "navyblue", linetype = "solid",linewidth = .5)+
  annotate("segment", x = 3250, y = 0.0022, xend = 3200, yend = 0.0022, color = "navyblue", linetype = "dotted",linewidth = .5)+
  #Early bronze start
  annotate("text", x=4700,y=0.002115, label = "Late Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4200, xend = 4200, y = 0.0022, yend = 0.0021,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Early Dilmun start
  annotate("text", x=3875,y=0.002115, label = "Late Dilmun", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3550, xend = 3550, y = 0.0022, yend = 0.0021,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Middle Dilmun start
  annotate("text", x=3375,y=0.002115, label = "Middle Dilmun", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.25, 'cm'),legend.position = c(.11, .7),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.4,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

Figure_xd = ggplot(extracted_perm.plot, aes(x = calBP)) +
  geom_ribbon(aes(ymin=nw.intersect.1, ymax=nw.intersect.2, fill = "Envelope"), alpha = 0.5)+
  #geom_ribbon(aes(ymin=nw.intersect.2, ymax=nw.above, fill = ">95%"), alpha = 0.5)+
  #geom_ribbon(aes(ymin=nw.intersect.1, ymax=nw.below, fill = "<5%"), alpha = 0.5)+
  geom_line(aes(y = Northwest, colour = "Northwest SPD"), size=1, show.legend = FALSE)+
  scale_fill_manual(values=c("Envelope"="grey75", ">95%" = "blue", "<5%" = "red"))+
  scale_colour_manual(values=c("Northwest SPD"="#F4C431"))+
  scale_y_continuous(limits = c(0,0.0018), breaks=seq(0,0.0018, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  #plot the number of dates
  annotate("text", x=5200,y=0.00075, label = paste0("n = ",nrow(n.nw.dates)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  #plot the p-value
  annotate("text", x=5200,y=0.00065, label = paste0("global p-value = ",round((perm.all.dates$pValueList[4]), 5)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  #Northwest
  annotate("text", x=5175,y=0.001675, label = "Northwest", color = "navyblue", 
           angle = 0, hjust = 0.0, size = 4, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, y = 0.0016, xend = 5200, yend = 0.0016, color = "navyblue", linetype = "dotted",linewidth = .5)+
  annotate("segment", x = 3200, y = 0.0016, xend = 5150, yend = 0.0016, color = "navyblue", linetype = "solid",linewidth = .5)+
  #Late EB I end
  annotate("text", x=5025,y=0.001535, label = "Late EB I", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4850, xend = 4850, y = 0.0016, yend = 0.00155,  color = "navyblue",linetype="solid",linewidth = .5) +
  #EB II end
  annotate("text", x=4750,y=0.001535, label = "EB II", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.0016, yend = 0.00155,  color = "navyblue",linetype="solid",linewidth = .5) +
  #EB III end
  annotate("text", x=4450,y=0.001535, label = "EB III", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4250, xend = 4250, y = 0.0016, yend = 0.00155,  color = "navyblue",linetype="solid",linewidth = .5) +
  #EB IV end
  annotate("text", x=4100,y=0.001535, label = "EB IV", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3950, xend = 3950, y = 0.0016, yend = 0.00155,  color = "navyblue",linetype="solid",linewidth = .5) +
  #MBA end
  annotate("text", x=3750,y=0.001535, label = "Middle Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3550, xend = 3550, y = 0.0016, yend = 0.00155,  color = "navyblue",linetype="solid",linewidth = .5) +
  #LBA end
  annotate("text", x=3375,y=0.001535, label = "Late Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3200, xend = 3200, y = 0.0016, yend = 0.00155,  color = "navyblue",linetype="solid",linewidth = .5) +
  ##
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.25, 'cm'),legend.position = c(.11, .6),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.4,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

####connect figure 3a & 3b together into one.
Figure_X = ggarrange(Figure_xa,Figure_xb,Figure_xc,Figure_xd,
                     labels = c("A", "B", "C", "D"),
                     ncol = 1, nrow = 4)

Figure_X

#save figure 3  to file.
ggsave("Figure_X.pdf",  plot = last_plot(), device = "pdf",
       path = "Master Paper/paper pdf/", width = 3508, height = 3500,
       units = 'px', dpi = 300)




ggplot(joined_regionals, aes(x = calBP)) +
  geom_line(aes(y = Southwest, colour = "Southwest"), size=1)+
  geom_line(aes(y = Southeast, colour = "Southeast"), size=1)+
  geom_line(aes(y = Bahrain, colour = "Bahrain"), size=1)+
  geom_line(aes(y = Northwest, colour = "Northwest"), size=1)+
  # scale_colour_manual(values=c("Southwest"="#50786B","Normalised SPD"="#F4C431"))+
  #this is to show the pearsons correlation between normalised and unnormalised SPDs
  #scale_y_continuous(limits = c(0,0.002), breaks=seq(0,0.002, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.5, 'cm'),legend.position = c(.1415, .925),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))


