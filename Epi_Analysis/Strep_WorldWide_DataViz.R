
library(tidyverse)
library(readxl)
library(epitools)

wd = "C:/Users/dac23/Downloads"
wd = "/home/ron/Downloads"
setwd(wd)

# 12F 
# Strep <- read_excel("12F_international_incidence_v2_DC_01.04.2024.xlsx", col_names = T)

# 1
Strep <- read_excel("12F_international_incidence_v2_DC_01.04.2024.xlsx", sheet = "1_incidence_Clean", col_names = T)

# Vaccine
Vaccine <- read_excel("12F_international_incidence_v2_DC_01.04.2024.xlsx", sheet = "Vaccine", col_names = T)

# view(Strep)
glimpse(Strep)
unique(sort(Strep$Region))
unique(sort(Strep$Area))
unique(sort(Strep$Period))
unique(sort(Strep$Demographic))
unique(sort(Strep$Reference))

# Incidence components:
unique(sort(Strep$Count))
unique(sort(Strep$Total))

# Vaccines:
glimpse(Vaccine)
unique(sort(Vaccine$Vaccine))


# Things to consider before incidence analysis: ################################
# 1. Have to check why some data in <Total> have zero value.
# 2. <Region> & <Area> are fine
# 3. <Period> & <Demographic> is quite messy, need to cross-check the reference.

# Some arguments and grumblings ################################################
# Idk. I think it's impossible to group the data considering the various Perd & Demo value?
ZeroTot <- Strep %>% 
  filter(Total == 0) %>% 
  # view() %>% 
  glimpse()
# That's it. The source is located in Hong Kong (2017 & 2019).
# I checked the Ref again and the data's fine & make sense.

IntervalPerd <- Strep %>% 
  filter(grepl("-", Period)) %>% 
  view() %>% 
  glimpse()
# The data that have a messy interval Perd also have messy Age structures:
unique(sort(IntervalPerd$Region))
# [1] "Belgium" "Germany" "Israel"  "Japan"   "Kenya"   "Morocco" "USA"

# Q: Is the Perd-filtered messy data specific to these Region? (manually filter the data)
# If yes, the analysis should be done in country-specific approach.
# Belgium = y (1 Ref)
# Germany = n (3 Refs)
# Israel  = y
# Japan   = n (2 Refs)
# Kenya   = y           --> huge time range (1999-2010 AND 2012-2016)
# Morocco = y           --> in range of 4 OR 5 years (inconsistent)
# USA     = y           --> freaking weird data with age-grouping = "All" -_-)"

unique(sort(IntervalPerd$Demographic))
unique(IntervalPerd$Reference) # ugh

# Decision:
# 1. DataViz as simple as (Period, Incidence, col = ALL_Ages, Children, Adult, Elderly) col according to WHO
# 2. Data is separated based on visualisation:
#   2.1. Worldwide
#   2.2. Area
#   2.3. Region --> Region-based DataViz can be vary based on the original data (Germany & Japan are special cases)
# 3. Keypoints for DataViz components:
#   3.1. Period with intervals (e.g. 2000-2004) will use midpoints instead (e.g. 2002)
#   3.2. Incidence components (Count, Total) will be summed based on the group
#   3.3. Special to Hong Kong: Total with value == 0 will be ommitteed (for CI calculations)
#   3.3. Confidence Intervals then calculated for each dataframe
#   3.4. Vaccine data as geom_vline OR abline(v = PeriodValue)

# The said functions are FUN:
FunYearMid <- function(eyy, delimit = "-") {
  range <- strsplit(eyy, delimit)[[1]]
  
  if (length(range ==2)) {
    numb <- as.numeric(range)
    midp <- mean(numb, na.rm = T)
    
  } else {
    midp <- as.numeric(eyy)
  }
  
  return(midp)
}

Strep <- Strep %>% 
  filter(!is.na(Count),
         !is.na(Total)) %>% # Coz data is considered uncleaned (temporary)
  mutate(Count = as.integer(Count),
         Total = as.integer(Total)) %>% 
  mutate(Demographic2 = case_when(
    Demographic == "<5"  ~ "Toddler",
    Demographic %in% c(">5","5-14","5-17","5-19","5-64",
                       "Children","15-29","15-44","15-59","15-64","<16","<18") ~ "Children",
    
    Demographic %in% c("Adults","≥16","≥18","18-49","19-49","30-49",
                       "45-64","50-59","50-64") ~ "Adults",
    
    Demographic %in% c("≥60","≥65") ~ "Elderly",
    Demographic == "All" ~ "All",
    TRUE ~ "other_value"
  ))


# 2. DataViz ###################################################################
# 2.1. WorldWide ###############################################################

# No grouping by ages
Strep_1ww_ALLages <- Strep %>% 
  # filter(Total != 0) %>% # Only required when we group the data specifically
  mutate(New_Period = sapply(Period, FUN = FunYearMid)) %>% 
  group_by(New_Period) %>% 
  summarise(sum_Count = sum(Count),
            sum_Total = sum(Total)) %>%
  ungroup() %>% 
  mutate(Conf_Int = binom.exact(sum_Count, sum_Total)) %>% 
  # view() %>% 
  glimpse()

max_up <- max(Strep_1ww_ALLages$Conf_Int$upper)+.01
plot(Strep_1ww_ALLages$New_Period, Strep_1ww_ALLages$Conf_Int$proportion,
     ylim = c(0, max_up), cex = Strep_1ww_ALLages$sum_Total/5000,
     main = "The Incidence of Serotype 1 from Publictly-Available Data Worldwide")
segments(Strep_1ww_ALLages$New_Period, Strep_1ww_ALLages$Conf_Int$lower,
         Strep_1ww_ALLages$New_Period, Strep_1ww_ALLages$Conf_Int$upper, col = "black")


# Grouping by ages
Strep_1ww_GRages <- Strep %>% 
  filter(Total != 0) %>% # Only required when we group the data specifically
  mutate(New_Period = sapply(Period, FUN = FunYearMid),
         New_Demographic = "TOBECONTINUED") %>% 
  group_by(New_Period, Demographic2) %>% # Instead of Demographic
  summarise(sum_Count = sum(Count),
            sum_Total = sum(Total)) %>%
  ungroup() %>% 
  mutate(Conf_Int = binom.exact(sum_Count, sum_Total)) %>% 
  # view() %>% 
  glimpse()

# Define the desired colors for each demographic value
col_map1 <- c("<5" = "indianred4", # < 5
              
              ">5" = "indianred3", # more or less 5-65
              "5-14" = "indianred3",
              "5-17" = "indianred3",
              "5-19" = "indianred3",
              "5-64" = "indianred3",
              
              "Children" = "indianred1",
              "15-29" = "indianred1",
              "15-44" = "indianred1",
              "15-59" = "indianred1",
              "15-64" = "indianred1",
              "<16" = "indianred1",
              "<18" = "indianred1",
              
              "Adults" = "seagreen1",
              "≥16" = "seagreen1",
              "≥18" = "seagreen1",
              "18-49" = "seagreen1",
              "19-49" = "seagreen1",
              "30-49" = "seagreen1",
              
              "45-64" = "seagreen3",
              "50-59" = "seagreen3",
              "50-64" = "seagreen3",
              
              "≥60" = "purple3", # > 60
              "≥65" = "purple3",
              "All" = "black")

col_map2 <- c("Toddler" = "indianred4", # < 5
              "Children" = "indianred3", # more or less 5-65
              "Adults" = "seagreen3",
              "Elderly" = "purple3", # > 60
              "All" = "black")

# Create a vector of colors based on the demographic values
col <- col_map2[Strep_1ww_GRages$Demographic2] # Instead of Demographic

max_up <- max(Strep_1ww_GRages$Conf_Int$upper)+.01
plot(Strep_1ww_GRages$New_Period, Strep_1ww_GRages$Conf_Int$proportion,
     pch = 19, col = col,
     ylim = c(0, max_up), cex = 1.5, # cex = Strep_1ww_GRages$sum_Total/800, 
     main = "The Incidence of Serotype 1 from Publictly-Available Data Worldwide")
segments(Strep_1ww_GRages$New_Period, Strep_1ww_GRages$Conf_Int$lower,
         Strep_1ww_GRages$New_Period, Strep_1ww_GRages$Conf_Int$upper, col = col,)

legend("topleft", legend = c("Toddler","Children","Adults","Elderly","All"),
       cex = 1, pch = 19, bty = "n", bg = "transparent",
       col = col_map2)



# 2.2. Facet-wrap by <Area> ####################################################

# No grouping by ages
Strep_2Area_ALLages <- Strep %>% 
  # filter(Total != 0) %>% # Only required when we group the data specifically
  mutate(New_Period = sapply(Period, FUN = FunYearMid)) %>% 
  group_by(Area, New_Period) %>% 
  summarise(sum_Count = sum(Count),
            sum_Total = sum(Total)) %>%
  ungroup() %>% 
  mutate(Conf_Int = binom.exact(sum_Count, sum_Total)) %>% 
  # view() %>% 
  glimpse()

ggplot(Strep_2Area_ALLages, aes(x = New_Period, y = Conf_Int$proportion)) +
  geom_point() +
  geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
                width = .1, color = "black") +
  ggtitle("The Incidence of Serotype 1 Grouped by Area") +
  facet_wrap(~ Area)


# Grouping by ages
Strep_2Area_GRages <- Strep %>% 
  filter(Total != 0) %>% # Only required when we group the data specifically
  mutate(New_Period = sapply(Period, FUN = FunYearMid)) %>% 
  group_by(Area, New_Period, Demographic2) %>% 
  summarise(sum_Count = sum(Count),
            sum_Total = sum(Total)) %>%
  ungroup() %>% 
  mutate(Conf_Int = binom.exact(sum_Count, sum_Total)) %>% 
  # view() %>% 
  glimpse()

ggplot(Strep_2Area_GRages, aes(x = New_Period, y = Conf_Int$proportion,
                               color = factor(Demographic2,
                                                      levels = c("Toddler",
                                                                 "Children",
                                                                 "Adults",
                                                                 "Elderly",
                                                                 "All")))) +
  geom_point() +
  geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
                width = .1) +
  scale_color_manual(values = col_map2, name = "Demographic") +
  ggtitle("The Incidence of Serotype 1 Grouped by Area") +
  facet_wrap(~ Area, scales = "free_y")



# 2.3. Facet-wrap by <Region> ##################################################

# No grouping by ages
Strep_3Region_ALLages <- Strep %>% 
  # filter(Total != 0) %>% # Only required when we group the data specifically
  mutate(New_Period = sapply(Period, FUN = FunYearMid)) %>% 
  group_by(Area, New_Period, Region) %>% 
  summarise(sum_Count = sum(Count),
            sum_Total = sum(Total)) %>%
  ungroup() %>% 
  mutate(Conf_Int = binom.exact(sum_Count, sum_Total)) %>% 
  # view() %>% 
  glimpse()

unique(sort(Vaccine$Vaccine))
Vaccine$Vaccine <- factor(Vaccine$Vaccine,
                          levels = c("PCV7", "PCV10", "PCV13", "PCV10 & PCV13")) # coz of that weird automatic levels

ggplot(Strep_3Region_ALLages, aes(x = New_Period, y = Conf_Int$proportion)) +
  geom_point() +
  geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
                width = .1, color = "black") +
  geom_vline(data = Vaccine, aes(xintercept = Period, colour = Vaccine), linetype = "dashed") +
  scale_color_manual(values = c("PCV7" = "gray80",
                                "PCV10" = "gray60",
                                "PCV13" = "gray20",
                                "PCV10 & PCV13" = "gray10")) +
  ggtitle("The Incidence of Serotype 1 Specific to Regions") +
  facet_wrap(~ Region)


# Grouping by ages
Strep_3Region_GRages <- Strep %>% 
  filter(Total != 0) %>% # Only required when we group the data specifically
  mutate(New_Period = sapply(Period, FUN = FunYearMid)) %>% 
  group_by(Area, New_Period, Region, Demographic2) %>% 
  summarise(sum_Count = sum(Count),
            sum_Total = sum(Total)) %>%
  ungroup() %>% 
  mutate(Conf_Int = binom.exact(sum_Count, sum_Total)) %>% 
  # view() %>% 
  glimpse()

Vaccine$Vaccine <- factor(Vaccine$Vaccine,
                          levels = c("PCV7", "PCV10", "PCV13", "PCV10 & PCV13")) # coz of that weird automatic levels

ggplot(Strep_3Region_GRages, aes(x = New_Period, y = Conf_Int$proportion,
                               color = factor(Demographic2,
                                              levels = c("Toddler",
                                                         "Children",
                                                         "Adults",
                                                         "Elderly",
                                                         "All")))) +
  geom_point() +
  geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
                width = .1) +
  scale_color_manual(values = c(col_map2,
                                "PCV7" = "gray80",
                                "PCV10" = "gray60",
                                "PCV13" = "gray20",
                                "PCV10 & PCV13" = "gray10"),
                     name = "Demographic") +
  geom_vline(data = Vaccine, aes(xintercept = Period,
                                 colour = Vaccine),
             linetype = "dashed") +
  ggtitle("The Incidence of Serotype 1 Specific to Regions") +
  facet_wrap(~ Region, scales = "free_y")
