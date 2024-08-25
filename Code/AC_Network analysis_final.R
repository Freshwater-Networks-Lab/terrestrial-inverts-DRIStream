#### Code for constructing terrestrial food webs along a drought gradient####

#by Alexandra Cooper(adapted from code created by Fredric Windsor and Tom Aspin)

#################################################################

####CONTENTS####

#1. Set up workspace
#2. Create food webs across all 21 channels
#3. Data manipulation
#4. Community analysis
#5. Abundance indices
#6. Calculate sampling completeness of the networks
#7. Calculate topological metrics
#8. Statistical analysis of abundance data
#9. Statistical analysis of food web topology metrics
#10. Species Richness (S) 
#11. Number of Links (L) 
#12. Edge density (LS)
#13. Connectance (C) 
#14. Fraction basal nodes (B)
#15. Fraction of intermediate nodes (I) 
#16. Fraction of top nodes (t)
#17. Food chain length (FCL)
#18. Trophic generality (InDeg)
#19. Plot figures: trophic networks 

#############################################################

####1.Set up workspace ####

# Clear the environment
rm(list=ls())

# Set working directory
setwd

# Load required packages 
library(reshape2); library(ggplot2);library(gridExtra); library(grid); 
library(gtable); library(dplyr); library(cheddar)

# For further information on using cheddar:
vignette('CheddarQuickstart') 
vignette('Community')        
vignette('PlotsAndStats')
vignette('Collections')
vignette('ImportExport')

# Read in the WebBuilder function
source("Code/Functions/WebBuilder.R")

# Load the registry data (trophic links from literature search)
registry <- read.csv('Data/Trophic links database_terrestrial inverts_AC_final.csv')

###################################################################

####2. Create food webs across all 21 channels####

#Channel 1
# Read in the csv 
ch1 <- read.csv("Data/ch1.csv")

# Remove taxa with non positive abundance (NAs)
ch1_clean <- subset(ch1, N >= 0)

# Create a community
ch1_comm <- Community(nodes = ch1_clean, 
                      trophic.links = NULL, 
                      properties = list(title = "Channel 1",
                                        DI = "0.10666",
                                        N.units = "m^-2"))

# Extract the node properties
ch1_nps <- NPS(ch1_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch1_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch1_clean$Consumer.Method

# Store the node data for Channel 1
ch1_nodes <- cbind(ch1_nps[,c(1,3:8)], # Extract the necessary columns
                   minimum.res.method, 
                   minimum.con.method)

# Store the link data for Channel 1 using the registry
ch1_links <- WebBuilder(ch1_nodes, registry,
                        method = c("exact", "genus", 
                                   "subfamily", "family",
                                   "order", "class"))

# Create a web object for Channel 1
ch1_web <- Community(properties = CPS(ch1_comm), 
                     nodes = NPS(ch1_comm),
                     trophic.links = ch1_links)

# Remove objects that are not necessary to keep
rm(ch1_nps, minimum.con.method, minimum.res.method)

#Channel 2 
# Read in the csv 
ch2 <- read.csv("Data/ch2.csv")

# Remove taxa with non positive abundance (NAs)
ch2_clean <- subset(ch2, N >= 0)

# Create a community
ch2_comm <- Community(nodes = ch2_clean, 
                      trophic.links = NULL, 
                      properties = list(title = "Channel 2",
                                        DI = "0.22935",
                                        N.units = "m^-2"))

# Extract the node properties
ch2_nps <- NPS(ch2_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch2_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch2_clean$Consumer.Method

# Store the node data for Channel 2
ch2_nodes <- cbind(ch2_nps[,c(1,3:8)], 
                   minimum.res.method,
                   minimum.con.method)

# Store the link data for Channel 2
ch2_links <- WebBuilder(ch2_nodes, registry, 
                        method = c("exact", "genus", 
                                   "subfamily", "family",
                                   "order", "class"))

# Create a web object for Channel 2
ch2_web <- Community(properties = CPS(ch2_comm), 
                     nodes = NPS(ch2_comm),
                     trophic.links = ch2_links)

# Remove objects that are not necessary to keep
rm(ch2_nps, minimum.con.method, minimum.res.method)

#Channel 3
# Read in the csv 
ch3 <- read.csv("Data/ch3.csv")

# Remove taxa with non positive abundance (NAs)
ch3_clean <- subset(ch3, N >= 0)

# Create a community
ch3_comm <- Community(nodes = ch3_clean, 
                      trophic.links = NULL, 
                      properties = list(title = "Channel 3",
                                         DI = "0.44182", 
                                         N.units = "m^-2"))

# Extract the node properties
ch3_nps <- NPS(ch3_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch3_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch3_clean$Consumer.Method

# Store the node data for Channel 3
ch3_nodes <- cbind(ch3_nps[,c(1,3:8)], minimum.res.method, minimum.con.method)

# Store the link data for Channel 3
ch3_links <- WebBuilder(ch3_nodes, registry, 
                        method = c("exact", "genus", 
                                   "subfamily", "family",
                                   "order", "class"))

# Create a web object for Channel 3
ch3_web <- Community(properties = CPS(ch3_comm),
                     nodes = NPS(ch3_comm),
                     trophic.links = ch3_links)

# Remove objects that are not necessary to keep
rm(ch3_nps, minimum.con.method, minimum.res.method)

#Channel 4
# Read in the csv 
ch4 <- read.csv("Data/ch4.csv")

# Remove taxa with non positive abundance (NAs)
ch4_clean <- subset(ch4, N >= 0)

# Create a community
ch4_comm <- Community(nodes = ch4_clean, trophic.links = NULL, 
                      properties = list(title = "Channel 4", 
                                        DI = "0.37060", 
                                        N.units = "m^-2"))

# Extract the node properties
ch4_nps <- NPS(ch4_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch4_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch4_clean$Consumer.Method

# Store the node data for Channel 4
ch4_nodes <- cbind(ch4_nps[,c(1,3:8)], minimum.res.method, minimum.con.method)

# Store the link data for Channel 4
ch4_links <- WebBuilder(ch4_nodes, registry,
                        method = c("exact", "genus",
                                   "subfamily", "family",
                                   "order", "class"))

# Create a web object for Channel 4
ch4_web <- Community(properties = CPS(ch4_comm),
                     nodes = NPS(ch4_comm),
                     trophic.links = ch4_links)

# Remove objects that are not necessary to keep
rm(ch4_nps, minimum.con.method, minimum.res.method)

#Channel 5
# Read in the csv 
ch5 <- read.csv("Data/ch5.csv")

# Remove taxa with non positive abundance (NAs)
ch5_clean <- subset(ch5, N >= 0)

# Create a community
ch5_comm <- Community(nodes = ch5_clean, trophic.links = NULL, 
                      properties = list(title = "Channel 5", 
                                        DI = "0.62503", 
                                        N.units = "m^-2"))

# Extract the node properties
ch5_nps <- NPS(ch5_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch5_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch5_clean$Consumer.Method

# Store the node data for Channel 5
ch5_nodes <- cbind(ch5_nps[,c(1,3:8)], minimum.res.method, minimum.con.method)

# Store the link data for Channel 5
ch5_links <- WebBuilder(ch5_nodes, registry,
                        method = c("exact", "genus",
                                   "subfamily", "family",
                                   "order", "class"))

# Create a web object for Channel 5
ch5_web <- Community(properties = CPS(ch5_comm), 
                     nodes = NPS(ch5_comm),
                     trophic.links = ch5_links)

# Remove objects that are not necessary to keep
rm(ch5_nps, minimum.con.method, minimum.res.method)

#Channel 6
# Read in the csv 
ch6 <- read.csv("Data/ch6.csv")

# Remove taxa with non positive abundance (NAs)
ch6_clean <- subset(ch6, N >= 0)

# Create a community
ch6_comm <- Community(nodes = ch6_clean, trophic.links = NULL, 
                      properties = list(title = "Channel 6",
                                        DI= "0.00000",
                                        N.units = "m^-2"))

# Extract the node properties
ch6_nps <- NPS(ch6_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch6_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch6_clean$Consumer.Method

# Store the node data for Channel 6
ch6_nodes <- cbind(ch6_nps[,c(1,3:8)],
                   minimum.res.method,
                   minimum.con.method)

# Store the link data for Channel 6
ch6_links <- WebBuilder(ch6_nodes, registry,
                        method = c("exact", "genus", 
                                   "subfamily", "family",
                                   "order", "class"))

# Create a web object for Channel 6
ch6_web <- Community(properties = CPS(ch6_comm),
                     nodes = NPS(ch6_comm),
                     trophic.links = ch6_links)

# Remove objects that are not necessary to keep
rm(ch6_nps, minimum.con.method, minimum.res.method)

#Channel 7
# Read in the csv 
ch7 <- read.csv("Data/ch7.csv")

#Remove taxa with non positive abundance (NAs)
ch7_clean <- subset(ch7, N >= 0)

# Create a community
ch7_comm <- Community(nodes = ch7_clean, trophic.links = NULL, 
                      properties = list(title = "Channel 7",
                                        DI = "0.62073",
                                        N.units = "m^-2"))

# Extract the node properties
ch7_nps <- NPS(ch7_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch7_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch7_clean$Consumer.Method

# Store the node data for Channel 7
ch7_nodes <- cbind(ch7_nps[,c(1,3:8)],
                   minimum.res.method,
                   minimum.con.method)

# Store the link data for Channel 7
ch7_links <- WebBuilder(ch7_nodes, registry, 
                        method = c("exact", "genus",
                                   "subfamily", "family",
                                   "order", "class"))

# Create a web object for Channel 7
ch7_web <- Community(properties = CPS(ch7_comm),
                     nodes = NPS(ch7_comm),
                     trophic.links = ch7_links)

# Remove objects that are not necessary to keep
rm(ch7_nps, minimum.con.method, minimum.res.method)

#Channel 8
# Read in the csv 
ch8 <- read.csv("Data/ch8.csv")

# Remove taxa with non positive abundance (NAs)
ch8_clean <- subset(ch8, N >= 0)

# Create a community
ch8_comm <- Community(nodes = ch8_clean, trophic.links = NULL, 
                      properties = list(title = "Channel 8",
                                        DI = "0.85876", 
                                        N.units = "m^-2"))

# Extract the node properties
ch8_nps <- NPS(ch8_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch8_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch8_clean$Consumer.Method

# Store the node data for Channel 8
ch8_nodes <- cbind(ch8_nps[,c(1,3:8)],
                   minimum.res.method, 
                   minimum.con.method)

# Store the link data for Channel 8
ch8_links <- WebBuilder(ch8_nodes, registry, 
                        method = c("exact", "genus",
                                   "subfamily", "family",
                                    "order", "class"))

# Create a web object for Channel 8
ch8_web <- Community(properties = CPS(ch8_comm), 
                     nodes = NPS(ch8_comm),
                     trophic.links = ch8_links)

# Remove objects that are not necessary to keep
rm(ch8_nps, minimum.con.method, minimum.res.method)

#Channel 9
# Read in the csv 
ch9 <- read.csv("Data/ch9.csv")

# Remove taxa with non positive abundance (NAs)
ch9_clean <- subset(ch9, N >= 0)

# Create a community
ch9_comm <- Community(nodes = ch9_clean, 
                      trophic.links = NULL, 
                      properties = list(title = "Channel 9",
                                        DI = "0.26769", 
                                        N.units = "m^-2"))

# Extract the node properties
ch9_nps <- NPS(ch9_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch9_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch9_clean$Consumer.Method

# Store the node data for Channel 9
ch9_nodes <- cbind(ch9_nps[,c(1,3:8)], 
                   minimum.res.method, 
                   minimum.con.method)

# Store the link data for Channel 9
ch9_links <- WebBuilder(ch9_nodes, registry, 
             method = c("exact", "genus", 
                        "subfamily", "family",
                        "order", "class"))

# Create a web object for Channel 9
ch9_web <- Community(properties = CPS(ch9_comm), 
                     nodes = NPS(ch9_comm),
                     trophic.links = ch9_links)

# Remove objects that are not necessary to keep
rm(ch9_nps, minimum.con.method, minimum.res.method)

#Channel 10
# Read in the csv 
ch10 <- read.csv("Data/ch10.csv")

# Remove taxa with non positive abundance (NAs)
ch10_clean <- subset(ch10, N >= 0)

# Create a community
ch10_comm <- Community(nodes = ch10_clean, trophic.links = NULL, 
                       properties = list(title = "Channel 10",
                                         DI = "0.00443", 
                                         N.units = "m^-2"))

# Extract the node properties
ch10_nps <- NPS(ch10_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch10_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch10_clean$Consumer.Method

# Store the node data for Channel 10
ch10_nodes <- cbind(ch10_nps[,c(1,3:8)],
                    minimum.res.method,
                    minimum.con.method)

# Store the link data for Channel 10
ch10_links <- WebBuilder(ch10_nodes,registry, 
                         method=c("exact","genus", "subfamily", 
                                  "family","order","class"))

# Create a web object for Channel 10
ch10_web <- Community(properties = CPS(ch10_comm), 
                      nodes = NPS(ch10_comm),
                      trophic.links = ch10_links)

# Remove objects that are not necessary to keep
rm(ch10_nps, minimum.con.method, minimum.res.method)

#Channel 11
# Read in the csv 
ch11 <- read.csv("Data/ch11.csv")

# Remove taxa with non positive abundance (NAs)
ch11_clean <- subset(ch11, N >= 0)

# Create a community
ch11_comm <- Community(nodes = ch11_clean, trophic.links = NULL, 
                       properties = list(title = "Channel 11", 
                                         DI = "0.11765",
                                         N.units = "m^-2"))

# Extract the node properties
ch11_nps <- NPS(ch11_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch11_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch11_clean$Consumer.Method

# Store the node data for Channel 11
ch11_nodes <- cbind(ch11_nps[,c(1,3:8)],
                    minimum.res.method,
                    minimum.con.method)

# Store the link data for Channel 11
ch11_links <- WebBuilder(ch11_nodes, registry, 
                         method = c("exact", "genus", 
                                    "subfamily", "family",
                                    "order", "class"))

# Create a web object for Channel 11
ch11_web <- Community(properties = CPS(ch11_comm), 
                      nodes = NPS(ch11_comm),
                      trophic.links = ch11_links)

# Remove objects that are not necessary to keep
rm(ch11_nps, minimum.con.method, minimum.res.method)

#Channel 12
# Read in the csv 
ch12 <- read.csv("Data/ch12.csv")

# Remove taxa with non positive abundance (NAs)
ch12_clean <- subset(ch12, N >= 0)

# Create a community
ch12_comm <- Community(nodes = ch12_clean, trophic.links = NULL, 
                       properties = list(title = "Channel 12",
                                         DI = "0.03578", 
                                         N.units = "m^-2"))

# Extract the node properties
ch12_nps <- NPS(ch12_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch12_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch12_clean$Consumer.Method


# Store the node data for Channel 12
ch12_nodes <- cbind(ch12_nps[,c(1,3:8)],
                    minimum.res.method, 
                    minimum.con.method)

# Store the link data for Channel 12
ch12_links <- WebBuilder(ch12_nodes, registry, 
                         method = c("exact", "genus",
                                    "subfamily", "family",
                                    "order","class"))

# Create a web object for Channel 12
ch12_web <- Community(properties = CPS(ch12_comm),
                      nodes = NPS(ch12_comm),
                      trophic.links = ch12_links)

# Remove objects that are not necessary to keep
rm(ch12_nps, minimum.con.method, minimum.res.method)

#Channel 13
# Read in the csv 
ch13 <- read.csv("Data/ch13.csv")

# Remove taxa with non positive abundance (NAs)
ch13_clean <- subset(ch13, N >= 0)

# Create a community
ch13_comm <- Community(nodes = ch13_clean, trophic.links = NULL, 
                       properties = list(title = "Channel 13",
                                         DI = "0.36487",
                                         N.units = "m^-2"))

# Extract the node properties
ch13_nps <- NPS(ch13_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch13_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch13_clean$Consumer.Method

# Store the node data for Channel 13
ch13_nodes <- cbind(ch13_nps[,c(1,3:8)],
                    minimum.res.method,
                    minimum.con.method)

# Store the link data for Channel 13
ch13_links <- WebBuilder(ch13_nodes, registry,
                         method = c("exact", "genus",
                                    "subfamily", "family",
                                    "order", "class"))

# Create a web object for Channel 13
ch13_web <- Community(properties = CPS(ch13_comm), 
                      nodes = NPS(ch13_comm),
                      trophic.links = ch13_links)

# Remove objects that are not necessary to keep
rm(ch13_nps, minimum.con.method, minimum.res.method)

#Channel 14
# Read in the csv 
ch14 <- read.csv("Data/ch14.csv")

# Remove taxa with non positive abundance (NAs)
ch14_clean <- subset(ch14, N >= 0)

# Create a community
ch14_comm <- Community(nodes = ch14_clean, trophic.links = NULL, 
                       properties = list(title = "Channel 14",
                                         DI = "0.00442", 
                                         N.units = "m^-2"))

# Extract the node properties
ch14_nps <- NPS(ch14_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch14_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch14_clean$Consumer.Method

# Store the node data for Channel 14
ch14_nodes <- cbind(ch14_nps[,c(1,3:8)], 
                    minimum.res.method,
                    minimum.con.method)

# Store the link data for Channel 14
ch14_links <- WebBuilder(ch14_nodes,  registry,
                         method = c("exact", "genus", 
                                    "subfamily", "family",
                                    "order", "class"))

# Create a web object for Channel 14
ch14_web <- Community(properties = CPS(ch14_comm), 
                      nodes = NPS(ch14_comm),
                      trophic.links = ch14_links)

# Remove objects that are not necessary to keep
rm(ch14_nps, minimum.con.method, minimum.res.method)

#Channel 15
# Read in the csv 
ch15 <- read.csv("Data/ch15.csv")

# Remove taxa with non positive abundance (NAs)
ch15_clean <- subset(ch15, N >= 0)

# Create a community
ch15_comm <- Community(nodes = ch15_clean, trophic.links = NULL, 
                       properties = list(title = "Channel 15",
                                         DI = "0.99076", 
                                         N.units = "m^-2"))

# Extract the node properties
ch15_nps <- NPS(ch15_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch15_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch15_clean$Consumer.Method


# Store the node data for Channel 15
ch15_nodes <- cbind(ch15_nps[,c(1,3:8)],
                    minimum.res.method, 
                    minimum.con.method)

# Store the link data for Channel 15
ch15_links <- WebBuilder(ch15_nodes, registry, 
                         method = c("exact", "genus",
                                    "subfamily", "family",
                                    "order", "class"))

# Create a web object for Channel 15
ch15_web <- Community(properties = CPS(ch15_comm),
                      nodes = NPS(ch15_comm),
                      trophic.links = ch15_links)

# Remove objects that are not necessary to keep
rm(ch15_nps, minimum.con.method, minimum.res.method)

#Channel 16
# Read in the csv 
ch16 <- read.csv("Data/ch16.csv")

# Remove taxa with non positive abundance (NAs)
ch16_clean <- subset(ch16, N >= 0)

# Create a community
ch16_comm <- Community(nodes = ch16_clean, 
                       trophic.links = NULL, 
                       properties = list(title = "Channel 16",
                                         DI = "0.55430",
                                         N.units = "m^-2"))

# Extract the node properties
ch16_nps <- NPS(ch16_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch16_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch16_clean$Consumer.Method


# Store the node data for Channel 16
ch16_nodes <- cbind(ch16_nps[,c(1,3:8)], 
                    minimum.res.method,
                    minimum.con.method)

# Store the link data for Channel 16
ch16_links <- WebBuilder(ch16_nodes, registry, 
                         method = c("exact", "genus",
                                    "subfamily", "family",
                                    "order", "class"))

# Create a web object for Channel 16
ch16_web <- Community(properties = CPS(ch16_comm),
                      nodes = NPS(ch16_comm),
                      trophic.links = ch16_links)

# Remove objects that are not necessary to keep
rm(ch16_nps, minimum.con.method, minimum.res.method)

#Channel 17
# Read in the csv 
ch17 <- read.csv("Data/ch17.csv")

# Remove taxa with non positive abundance (NAs)
ch17_clean <- subset(ch17, N >= 0)

# Create a community
ch17_comm <- Community(nodes = ch17_clean, 
                       trophic.links = NULL, 
                       properties = list(title = "Channel 17",
                                         DI = "0.28358",
                                         N.units = "m^-2"))

# Extract the node properties
ch17_nps <- NPS(ch17_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch17_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch17_clean$Consumer.Method

# Store the node data for Channel 17
ch17_nodes <- cbind(ch17_nps[,c(1,3:8)], 
                    minimum.res.method, 
                    minimum.con.method)

# Store the link data for Channel 17
ch17_links <- WebBuilder(ch17_nodes, registry,
                         method=c("exact", "genus",
                                  "subfamily", "family",
                                  "order", "class"))

# Create a web object for Channel 17
ch17_web <- Community(properties = CPS(ch17_comm), 
                      nodes = NPS(ch17_comm),
                      trophic.links = ch17_links)

# Remove objects that are not necessary to keep
rm(ch17_nps, minimum.con.method, minimum.res.method)

#Channel 18
# Read in the csv 
ch18 <- read.csv("Data/ch18.csv")

# Remove taxa with non positive abundance (NAs)
ch18_clean <- subset(ch18, N >= 0)

# Create a community
ch18_comm <- Community(nodes = ch18_clean, 
                       trophic.links = NULL, 
                       properties = list(title = "Channel 18",
                                         DI = "0.40667",
                                         N.units = "m^-2"))

# Extract the node properties
ch18_nps <- NPS(ch18_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch18_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch18_clean$Consumer.Method


# Store the node data for Channel 18
ch18_nodes <- cbind(ch18_nps[,c(1,3:8)],
                    minimum.res.method,
                    minimum.con.method)

# Store the link data for Channel 18
ch18_links <- WebBuilder(ch18_nodes, registry,
                         method = c("exact","genus",
                                    "subfamily", "family",
                                    "order", "class"))

# Create a web object for Channel 18
ch18_web <- Community(properties = CPS(ch18_comm), 
                      nodes = NPS(ch18_comm),
                      trophic.links = ch18_links)

# Remove objects that are not necessary to keep
rm(ch18_nps, minimum.con.method, minimum.res.method)

#Channel 19
# Read in the csv 
ch19 <- read.csv("Data/ch19.csv")

# Remove taxa with non positive abundance (NAs)
ch19_clean <- subset(ch19, N >= 0)

# Create a community
ch19_comm <- Community(nodes = ch19_clean, 
                       trophic.links = NULL, 
                       properties = list(title = "Channel 19",
                                         DI = "0.38453",
                                         N.units = "m^-2"))

# Extract the node properties
ch19_nps <- NPS(ch19_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch19_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch19_clean$Consumer.Method

# Store the node data for Channel 19
ch19_nodes <- cbind(ch19_nps[,c(1,3:8)],
                    minimum.res.method,
                    minimum.con.method)

# Store the link data for Channel 19
ch19_links <- WebBuilder(ch19_nodes, registry,
                         method = c("exact", "genus",
                                    "subfamily", "family",
                                    "order", "class"))

# Create a web object for Channel 19
ch19_web <- Community(properties = CPS(ch19_comm),
                      nodes = NPS(ch19_comm),
                      trophic.links = ch19_links)

# Remove objects that are not necessary to keep
rm(ch19_nps, minimum.con.method, minimum.res.method)

#Channel 20
# Read in the csv 
ch20 <- read.csv("Data/ch20.csv")

# Remove taxa with non-positive abundance (NAs)
ch20_clean <- subset(ch20, N >= 0)

# Create a community
ch20_comm <- Community(nodes = ch20_clean, 
                       trophic.links = NULL, 
                       properties = list(title = "Channel 20",
                                         DI = "0.52676",
                                         N.units = "m^-2"))

# Extract the node properties
ch20_nps <- NPS(ch20_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch20_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch20_clean$Consumer.Method

# Store the node data for Channel 20
ch20_nodes <- cbind(ch20_nps[,c(1,3:8)],
                    minimum.res.method,
                    minimum.con.method)

# Store the link data for Channel 20
ch20_links <- WebBuilder(ch20_nodes, registry,
                         method = c("exact", "genus",
                                    "subfamily", "family",
                                    "order", "class"))

# Create a web object for Channel 20
ch20_web <- Community(properties = CPS(ch20_comm),
                      nodes = NPS(ch20_comm),
                      trophic.links = ch20_links)

# Remove objects that are not necessary to keep
rm(ch20_nps, minimum.con.method, minimum.res.method)

#Channel 21
# Read in the csv 
ch21 <- read.csv("Data/ch21.csv")

# Remove taxa with non-positive abundance (NAs)
ch21_clean <- subset(ch21, N >= 0)

# Create a community
ch21_comm <- Community(nodes = ch21_clean, trophic.links = NULL, 
                       properties = list(title = "Channel 21",
                                         DI = "1", 
                                         N.units = "m^-2"))

# Extract the node properties
ch21_nps <- NPS(ch21_comm)

# List the minimum resolution of the resources 
minimum.res.method <- ch21_clean$Resource.Method

# List the minimum resolution of the consumers
minimum.con.method <- ch21_clean$Consumer.Method

# Store the node data for Channel 21
ch21_nodes <- cbind(ch21_nps[,c(1,3:8)], 
                    minimum.res.method, 
                    minimum.con.method)

# Store the link data for Channel 21
ch21_links <- WebBuilder(ch21_nodes, registry, 
                         method = c("exact", "genus",
                                    "subfamily", "family",
                                    "order","class"))

# Create a web object for Channel 21
ch21_web <- Community(properties = CPS(ch21_comm), 
                      nodes = NPS(ch21_comm),
                      trophic.links = ch21_links)

# Remove objects that are not necessary 
rm(ch21_nps, minimum.con.method, minimum.res.method)

# Remove the final functions that are not needed for the RData file
rm(WebBuilder, registry)

# Save inferred webs as a RData file 
save.image(file = "Data/inferred.foodwebs.RData")

####################################################################

####3. Data manipulation####

# Load packages required:

library(cheddar); library(ggplot2);library(gridExtra); library(grid); 
library(gtable); library(dplyr);library(DescTools); library(purrr); 
library(segmented); library(magrittr); library(ggh4x); library(tidyr);
library(tidyverse); library(vegan);library(igraph); library(bipartite)


# Load in inferred foodwebs into environment (if needed) 
load(file = "Data/inferred.foodwebs.RData")
 
# Re-configure the data into different formats (i.e. adjacency  and incidence
# matrices, and edgelists) to facilitate further analysis.

# Create edgelists (i.e. a graphical representation of the lists of edges/links)
# from the trophic link data

edgelists <- list(ch6_links, ch14_links, ch10_links, 
                  ch12_links, ch1_links, ch11_links, 
                  ch2_links, ch9_links, ch17_links,
                  ch13_links, ch4_links, ch19_links, 
                  ch18_links, ch3_links, ch20_links, 
                  ch16_links, ch7_links, ch5_links,
                  ch8_links, ch15_links, ch21_links)

# Rename edgelists with individual channel number
names(edgelists) <- c("ch6", "ch14", "ch10", 
                      "ch12", "ch1", "ch11", 
                      "ch2",  "ch9", "ch17", 
                      "ch13", "ch4", "ch19", 
                      "ch18", "ch3", "ch20",
                      "ch16", "ch7", "ch5", 
                      "ch8", "ch15", "ch21")

# Create a list of nodes (species sampled from each channel)
nodes <- list(ch6, ch14, ch10, ch12, ch1, ch11, ch2, 
              ch9, ch17, ch13, ch4, ch19, ch18, ch3, 
              ch20, ch16, ch7, ch5, ch8, ch15, ch21)

# Name nodes by channel
names(nodes) <- c("ch6", "ch14", "ch10", "ch12", "ch1", "ch11", "ch2",
                  "ch9", "ch17", "ch13", "ch4", "ch19","ch18", "ch3", 
                  "ch20", "ch16", "ch7", "ch5", "ch8", "ch15", "ch21")

# Simplify the edgelists(to make it easier to convert them into a network later)

edgelists_simple <- lapply(edgelists, 
                    function(x){dplyr::select(x, resource, consumer)})

# Rename the simplified edgelists
names(edgelists_simple) <- c("ch6", "ch14", "ch10", 
                             "ch12", "ch1", "ch11", 
                             "ch2", "ch9", "ch17", 
                             "ch13", "ch4", "ch19", 
                             "ch18", "ch3", "ch20", 
                             "ch16", "ch7", "ch5", 
                             "ch8","ch15", "ch21")

# Create a list of network objects 
networks <- lapply(edgelists_simple, graph_from_data_frame)

# Rename the network objects by channel
names(networks) <- c("ch6", "ch14", "ch10", 
                     "ch12", "ch1", "ch11",
                     "ch2",  "ch9", "ch17",
                     "ch13", "ch4", "ch19", 
                     "ch18", "ch3", "ch20", 
                     "ch16", "ch7", "ch5", 
                     "ch8", "ch15", "ch21")

# Create adjacency matrices 

adj_matrices <- lapply(networks, 
                       function(x){as_adjacency_matrix(x, sparse = F)})

# Rename adjacency matrices by channel
names(adj_matrices) <- c("ch6", "ch14", "ch10",
                         "ch12", "ch1", "ch11",
                         "ch2", "ch9", "ch17",
                         "ch13", "ch4", "ch19",
                         "ch18", "ch3", "ch20", 
                         "ch16", "ch7", "ch5",
                         "ch8", "ch15", "ch21")

# Create incidence matrices (graphical representation of vertices and edges)
matrices <- lapply(adj_matrices, empty)

# Rename adjacency matrices
names(matrices) <- c("ch6", "ch14", "ch10", 
                     "ch12", "ch1", "ch11", 
                     "ch2",   "ch9", "ch17",
                     "ch13", "ch4", "ch19", 
                     "ch18", "ch3", "ch20",
                     "ch16", "ch7", "ch5", 
                     "ch8", "ch15", "ch21")

# Create a list of webs
webs <- list(ch6_web, ch14_web, ch10_web,
             ch12_web, ch1_web, ch11_web,
             ch2_web, ch9_web, ch17_web, 
             ch13_web, ch4_web, ch19_web, 
             ch18_web, ch3_web, ch20_web, 
             ch16_web, ch7_web, ch5_web, 
             ch8_web, ch15_web, ch21_web)

# Rename webs
names(webs) <- c("ch6", "ch14", "ch10", "ch12", "ch1", "ch11", "ch2", 
                 "ch9", "ch17", "ch13", "ch4", "ch19","ch18", "ch3", 
                 "ch20", "ch16", "ch7", "ch5", "ch8", "ch15", "ch21")

#Create a merged dataframe of links and nodes

# Sort the node data into long format
nodes_long <- bind_rows(nodes, .id = "Channel")
nodes_dframe <- dplyr::select(nodes_long, Channel, node, N)

# Sort the link data into long format
links_long <- bind_rows(edgelists, .id = "Channel")
links_dframe <- dplyr::select(links_long, 
                              Channel,
                              resource, 
                              consumer,
                              res.cat,con.cat)

# Bind the two datasets together 
links1 <- left_join(links_dframe, nodes_dframe, 
                    by = c("Channel", "resource" = "node"))

# Rename the columns to indicate that N relates to resource
colnames(links1)[6] <- c("N.res")

# Bind the two datasets together 
tritrophic_links <- left_join(links1, nodes_dframe, 
                              by = c("Channel", "consumer" = "node"))

# Rename the columns to indicate that N relates to consumer
colnames(tritrophic_links)[7] <- c("N.con")

# Order the channels by DI (drought index) and append
tritrophic_links$Channel <- factor(tritrophic_links$Channel, 
                                   levels = c("ch6", "ch14","ch10", 
                                              "ch12","ch1", "ch11",
                                              "ch2", "ch9", "ch17",
                                              "ch13", "ch4","ch19", 
                                              "ch18", "ch3", "ch20", 
                                              "ch16", "ch7", "ch5", 
                                              "ch8",  "ch15", "ch21"))

# Add trophic categories to the node_dframe 
nodes_dframe$cat <- nodes_long$cat
nodes_dframe$Channel <- factor(nodes_dframe$Channel, 
                               levels = c("ch6", "ch14", "ch10",
                                          "ch12", "ch1", "ch11",
                                          "ch2", "ch9", "ch17",
                                          "ch13", "ch4", "ch19",
                                          "ch18", "ch3", "ch20",
                                          "ch16", "ch7", "ch5",
                                          "ch8", "ch15", "ch21"))

# Remove objects and functions that are not  needed 
rm(edgelists_simple, links_long, nodes_long, links_dframe, links1)

#Save as RData file
save.image(file = "Data/Data.manipulation.RData")

#####################################################################

#Descriptive statistics

#remove NAs from node data
nodes_dframe2 <-  na.omit(nodes_dframe)

sum(nodes_dframe2$N) 
var((nodes_dframe2$N))
table(nodes_dframe2$cat) # N of each trophic category

write.csv(nodes_dframe2,'nodes.dframe2.csv')

####4. Community analysis####

# Load the previous file if needed
load(file = "Data/Data.manipulation.RData")

# Create a community dataframe for each channel in the list
comm_list <- lapply(nodes,
                    function(x){pivot_wider(dplyr::select(x, node, N),
                                                   names_from = node, 
                                                   values_from = N)})

# Bind them all together
comm_dframe <- bind_rows(comm_list, .id = "Channel")

# Convert NAs to 0s 
comm_dframe[is.na(comm_dframe)] <- 0

# DI coded as numeric and added to dataframe
comm_dframe$DI <- as.numeric(lapply(webs,
                                    function(x){x$properties$DI}))


# Create the broad categories (flowing, fragmented, dry)
comm_dframe$category <- c(rep("CON", 7), rep("FRAG", 11), 
                          rep("DRY", 3)) 

# Save as RData file.
save.image(file = "Data/Community-composition.RData")

###################################################################
####5. Abundance indices####

# Load the previous file if needed
load(file = "Data/Community-composition.RData")

# Load required packages into R as needed:
install.packages ("vegan")
install.packages("ggplot2")
install.packages ("ggrepel")
library(dplyr)
library(vegan)
library(lattice)
library(vegan)
library(ggplot2)
library(ggrepel)

# Create two separate data frames to add 3 biodiversity metrics, one for plants
# and one for inverts (due to different units of abundance)

# Copy columns containing site and species data (i.e. channel, DI and category) 
#for inverts and plants:

invert.div <- comm_dframe[c(1, 33:249, 258, 259)] 
plant.div <- comm_dframe[c(1:32, 258, 259)] 


#  Total abundance of invert is sum of rows containing invert data: 

invert.div$abund <- rowSums(invert.div[, 2:218])

# and for plants (columns 2-32)
plant.div$abund <-rowSums(plant.div[, 2:32])
 
#Calculate species richness for plants and inverts

invert.div$rich <- specnumber(invert.div[, 2:218])

plant.div$rich <- specnumber(plant.div[, 2:32])

#Calculate evenness (Simpson Diversity; aka Inverse Simpson Index) for plants 
# and inverts separately:

invert.div$even <- diversity(invert.div[, 2:218], index = "invsimpson") / 
  specnumber(invert.div[, 2:218])

plant.div$even<- diversity(plant.div[, 2:32], index = "invsimpson") / 
  specnumber(plant.div[, 2:32])

# Write and save as csv
write.csv(invert.div,'invert.div.csv')
write.csv(plant.div,'plant.div.csv')

####6. Calculate sampling completeness for the networks #### 

# Read in the source code for sampling completeness algorithms, after MacGregor 
# et al. 2017:
source("Code/Functions/sampling_completeness.R")

# Convert the data into a suitable format
matrices_comp <- lapply(matrices, 
                        function(x){cbind.data.frame("Species" = rownames(x),
                                                     x)})

# Create one dataframe that includes all 21 networks
matrices_dframe <- bind_rows(matrices_comp)

# Sum the values for the taxa to calculate frequency of detection

matrices_sum <- matrices_dframe %>% 
  group_by(Species) %>%
  summarise(across(everything(), sum))

# Remove NAs from the database
matrices_sum[is.na(matrices_sum)] <- 0

# Calculate sampling completeness for the consumers and resources separately
completeness_con <- SCw1(matrices_sum[,-1], estimator = "Chao1")
completeness_res <- SCw1(t(matrices_sum[,-1]), estimator = "Chao1")

# Find the median completeness, to give overall sampling completeness across 
# resources and consumers i.e. proportion of potential interactions realised.

(completeness_con[2] + completeness_res[2])/2 

####################################################################

####7. Calculate the topological metrics ####

library (cheddar)

# Remove isolated nodes from web objects

ch1_web.no.isolated <- RemoveIsolatedNodes(ch1_web)
ch2_web.no.isolated <- RemoveIsolatedNodes(ch2_web)
ch3_web.no.isolated <- RemoveIsolatedNodes(ch3_web)
ch4_web.no.isolated <- RemoveIsolatedNodes(ch4_web)
ch5_web.no.isolated <- RemoveIsolatedNodes(ch5_web)
ch6_web.no.isolated <- RemoveIsolatedNodes(ch6_web)
ch7_web.no.isolated <- RemoveIsolatedNodes(ch7_web)
ch8_web.no.isolated <- RemoveIsolatedNodes(ch8_web)
ch9_web.no.isolated <- RemoveIsolatedNodes(ch9_web)
ch10_web.no.isolated <- RemoveIsolatedNodes(ch10_web)
ch11_web.no.isolated <- RemoveIsolatedNodes(ch11_web)
ch12_web.no.isolated <- RemoveIsolatedNodes(ch12_web)
ch13_web.no.isolated <- RemoveIsolatedNodes(ch13_web)
ch14_web.no.isolated <- RemoveIsolatedNodes(ch14_web)
ch15_web.no.isolated <- RemoveIsolatedNodes(ch15_web)
ch16_web.no.isolated <- RemoveIsolatedNodes(ch16_web)
ch17_web.no.isolated <- RemoveIsolatedNodes(ch17_web)
ch18_web.no.isolated <- RemoveIsolatedNodes(ch18_web)
ch19_web.no.isolated <- RemoveIsolatedNodes(ch19_web)
ch20_web.no.isolated <- RemoveIsolatedNodes(ch20_web)
ch21_web.no.isolated <- RemoveIsolatedNodes(ch21_web)


# Remove cannabilistic links from web objects and add title to improve plotting
ch1_web.clean <- RemoveCannibalisticLinks(ch1_web.no.isolated, "Channel 1")
ch2_web.clean <- RemoveCannibalisticLinks(ch2_web.no.isolated, "Channel 2")
ch3_web.clean <- RemoveCannibalisticLinks(ch3_web.no.isolated, "Channel 3")
ch4_web.clean <- RemoveCannibalisticLinks(ch4_web.no.isolated, "Channel 4")
ch5_web.clean <- RemoveCannibalisticLinks(ch5_web.no.isolated, "Channel 5")
ch6_web.clean <- RemoveCannibalisticLinks(ch6_web.no.isolated, "Channel 6")
ch7_web.clean <- RemoveCannibalisticLinks(ch7_web.no.isolated, "Channel 7")
ch8_web.clean <- RemoveCannibalisticLinks(ch8_web.no.isolated, "Channel 8")
ch9_web.clean <- RemoveCannibalisticLinks(ch9_web.no.isolated, "Channel 9")
ch10_web.clean <- RemoveCannibalisticLinks(ch10_web.no.isolated, "Channel 10")
ch11_web.clean <- RemoveCannibalisticLinks(ch11_web.no.isolated, "Channel 11")
ch12_web.clean <- RemoveCannibalisticLinks(ch12_web.no.isolated, "Channel 12")
ch13_web.clean <- RemoveCannibalisticLinks(ch13_web.no.isolated, "Channel 13")
ch14_web.clean <- RemoveCannibalisticLinks(ch14_web.no.isolated, "Channel 14")
ch15_web.clean <- RemoveCannibalisticLinks(ch15_web.no.isolated, "Channel 15")
ch16_web.clean <- RemoveCannibalisticLinks(ch16_web.no.isolated, "Channel 16")
ch17_web.clean <- RemoveCannibalisticLinks(ch17_web.no.isolated, "Channel 17")
ch18_web.clean <- RemoveCannibalisticLinks(ch18_web.no.isolated, "Channel 18")
ch19_web.clean <- RemoveCannibalisticLinks(ch19_web.no.isolated, "Channel 19")
ch20_web.clean <- RemoveCannibalisticLinks(ch20_web.no.isolated, "Channel 20")
ch21_web.clean <- RemoveCannibalisticLinks(ch21_web.no.isolated, "Channel 21")


# Ensure the networks are ordered from low to high drought severity (DI) first.

collection <- CommunityCollection(list(ch6_web.clean, ch14_web.clean, 
                                       ch10_web.clean,ch12_web.clean,
                                       ch1_web.clean, ch11_web.clean,
                                       ch2_web.clean, ch9_web.clean, 
                                       ch17_web.clean,ch13_web.clean, 
                                       ch4_web.clean, ch19_web.clean,
                                       ch18_web.clean, ch3_web.clean,
                                       ch20_web.clean,ch16_web.clean,
                                       ch7_web.clean, ch5_web.clean,
                                       ch8_web.clean, ch15_web.clean,
                                       ch21_web.clean))

# Run topological metrics, and collate into one data frame

topology <- CollectionCPS(collection,
                          c('DI',
                            S = 'NumberOfNodes',
                            L = 'NumberOfTrophicLinks',
                            LS = 'LinkageDensity',
                            C = 'DirectedConnectance',
                            B = 'FractionBasalNodes',
                            I = 'FractionIntermediateNodes',
                            t = 'FractionTopLevelNodes'))
                          
                        
# Extract the channel number (for later)
topology$Channel <- substr(rownames(topology), 9, 10)

# Add in chain averaged trophic level metric (FCL):
topology$FCL <- unlist(lapply(collection, 
                              function(x){max(ChainAveragedTrophicLevel(x),
                                              na.rm = T)}))

# Add in InDegree metric (InDeg):
topology$InDeg <- unlist(lapply(collection, 
                              function(x){max(InDegree(x), na.rm = T)}))

## Remove all objects that are not required and save as RData file.
rm(SCw1, SCw2, species.interaction.completeness, 
   completeness_con,completeness_res, 
   matrices_sum, matrices_dframe, matrices_comp)

save.image(file = "Data/sampling completeness.RData")

#For later analysis of InDeg by trophic category
write.csv(topology,'InDeg.topology.csv')

##################################################################
####8. Statistical analysis of abundance data####

# Note only final models are shown here, All associated link functions were 
#trialled for GAMS/ GLMS as appropriate. The final model was selected on basis 
#that it gave the lowest AIC value. 

library(mgcv) # required for GAMs

# Load plant.csv and invert.csv
plant.div <- read.csv('plant.div.csv')
invert.div <-read.csv('invert.div.csv')

# Does total abundance vary with DI?
hist(invert.div$abund) 
hist(plant.div$abund) 

plot(invert.div$abund ~ invert.div$DI)
plot(plant.div$abund ~ plant.div$DI) 

# Negative binomial (to correct overdispersion of poisson models)

library(MASS) # required for negative binomials

model1<- glm.nb(abund~DI, link = "identity", data=invert.div)
AIC(model1)

summary(model1)
theta <- model1$deviance / model1$df.residual
theta # check overdispersion is around 1.
summary.lm(model1)
plot(model1) # check model assumptions

pseudo.r2<- (model1$null.deviance - model1$deviance) / model1$null.deviance
pseudo.r2 


# Plot predictions of final model (model1)
plot(abund ~ DI, data=invert.div,
     xlab= "Drought index (DI)",
     ylab="Invertebrate abundance")

pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(model1, newdata=pdat, na.rm=T, type="response", se.fit=TRUE)
pred

predframe <-data.frame(pdat, preds=pred$fit, se=pred$se.fit)
lines (predframe$preds ~ predframe$DI, col="red", lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5) # Adds lower SE line


# GAMs with total plant abundance

gam1<-gam(abund ~ s(DI, fx= FALSE, bs="tp", k= -1), 
          family=Gamma (link= "log"),
          data=plant.div) 

AIC(gam1) 
plot(gam1) 
summary.gam(gam1)
gam.check(gam1)

# Check model assumptions

devresid <- resid(gam1, type="deviance")
plot(devresid~gam1$fitted.values)
plot(devresid~plant.div$abund)

# Prediction plotting
plot(abund ~ DI, data=plant.div, 
     xlab= "Drought index (DI)",
     ylab = "Plant abundance" )

pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(gam1, newdata=pdat,
              na.rm=T, 
              type="response",
              se.fit=TRUE)
pred

predframe <-data.frame(pdat, 
                       preds=pred$fit,
                       se=pred$se.fit)

lines (predframe$preds ~ predframe$DI,
       col="red", 
       lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5)

# Simpson index/ Evenness (plants)

library(betareg) 

hist(plant.div$even) #  beta distribution between 0 and 1.
plot(plant.div$even~plant.div$DI) 

# Tranform evenness variable first as it contained 0 and 1s 
plant.div$even.adjusted<-plant.div$even-0.0000001

model2 <- betareg(even.adjusted ~ as.numeric(DI),
                   link = "log", 
                   link.phi="identity",
                   data =plant.div)
summary (model2)
AIC(model2) 

#Evenness inverts

hist(invert.div$even) 
plot(invert.div$even~invert.div$DI)

model3 <- betareg(even ~ as.numeric(DI),
                   link = "cauchit", 
                   link.phi="identity",
                   data = invert.div)
summary (model3)
AIC(model3)
summary(model3)# no significant differences

median(invert.div$even) 
var(invert.div$even)

# Does species richness of plants/ inverts vary with DI?

hist(invert.div$rich) 
hist(plant.div$rich)

plot(invert.div$rich~invert.div$DI) 
plot(plant.div$rich~plant.div$DI)  

#Quasipoisson to correct overdispersion

model4<- gam(rich~s(DI, fx= FALSE, bs="tp", k= -1), 
             family=quasipoisson (link="identity"), 
             data=invert.div)
             
theta <- model4$deviance / model4$df.residual
theta 

summary.gam(model4)
gam.check(model4)

pseudo.r2<- (model4$null.deviance - model4$deviance) / model4$null.deviance
pseudo.r2

devresid <- resid(model4, type="deviance")
plot(devresid~model4$fitted.values)
plot(devresid~invert.div$rich)

# Prediction plotting
plot(rich ~ DI, data=invert.div, 
     xlab= "Drought index (DI)",
     ylab = "Invertebrate richness" )

pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(model4, newdata=pdat,
              na.rm=T, 
              type="response",
              se.fit=TRUE)
pred

predframe <-data.frame(pdat, 
                       preds=pred$fit,
                       se=pred$se.fit)

lines (predframe$preds ~ predframe$DI,
       col="red", 
       lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5)

#GLM  withGamma (plant richness vs DI)

# Models with Gamma error family.
model5<-glm(rich ~ as.numeric(DI),
            family = Gamma (link = "identity"), 
            data = plant.div) 
AIC(model5)

summary(model5)
summary.lm(model5)
plot(model5)

# Prediction plotting
plot(rich ~ DI, data=plant.div, 
     xlab= "Drought index (DI)",
     ylab = "Plant richness" )

pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(model5, newdata=pdat,
              na.rm=T, 
              type="response",
              se.fit=TRUE)
pred

predframe <-data.frame(pdat, 
                       preds=pred$fit,
                       se=pred$se.fit)

lines (predframe$preds ~ predframe$DI,
       col="red", 
       lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5)


####9. Statistical analysis of food web topology metrics ####

load(file = "Data/sampling completeness.RData")

# Explore the effect of DI on:

#S = 'NumberOfNodes', 
#L = 'NumberOfTrophicLinks', 
#LS = 'LinkageDensity',
#C = 'DirectedConnectance',
#B = 'FractionBasalNodes',
#I = 'FractionIntermediateNodes',
#t = 'FractionTopLevelNodes', 
#FCL= 'Food chain length'
#InDeg= 'Trophic Generality',


#Ensure R reads DI as numeric
topology$DI<-as.numeric (topology$DI)

#Visualise the distributions of the variables
hist(topology$DI)
hist(topology$S) 
hist(topology$L)
hist(topology$C) 
hist(topology$LS) 
hist(topology$B) 
hist(topology$I) 
hist(topology$t)  
hist(topology$FCL) 
hist(topology$InDeg)

###################################################################

####10. Network size####
# Does network size (S) vary across foodwebs, as a function of DI?

plot(topology$DI~topology$S) 

# Negative binomial to correct overdispersion

# Load packages if needed
install.packages("MASS")
library(MASS)

model6<- glm.nb(S ~ as.numeric(DI), link = "identity", data = topology)

summary(model6, cor=FALSE)
summary.lm(model6)
theta <- model6$deviance / model6$df.residual
theta # check overdispersion

#Check model assumptions
devresid <- resid(model6, type = "deviance")
plot(devresid ~ model6$fitted.values) 
plot(devresid ~ as.numeric(topology$DI))

plot(model6)

# R squared values are not appropriate for binomial models such as poisson, 
#quasi poisson, and neg binomial, so need to calculate pseudo r- squared.

pseudo.r2 <- (model6$null.deviance - model6$deviance) / model6$null.deviance
pseudo.r2 

# Plot predictions of  model 6
plot(S ~ DI, data=topology, xlab= "Drought index (DI)", ylab="Network size (S)")
pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(model6, newdata=pdat, na.rm=T, type="response", se.fit=TRUE)
pred

predframe <-data.frame(pdat, preds=pred$fit, se=pred$se.fit)
lines (predframe$preds ~ predframe$DI, col="red", lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5) # Adds lower SE line

#####################################################################

####11. Number of links####

#Does the Number of Links (L) vary significantly across foodwebs, 
#as a function of DI?

# Visualise raw data
plot(topology$DI~topology$L) 

# Generalised  additive model (GAM), negative binomial, error family

gam7<-gam(L ~ s(DI, fx= FALSE, bs="tp", k= -1), 
          family=negbin(theta=100>1, link= "log"), 
          data=topology) 

plot(gam7) 
summary(gam7)   
AIC(gam7)

theta <- gam7$deviance / gam7$df.residual
theta # no longer overdispersed

devresid <- resid(gam7, type = "deviance")
plot(devresid ~ gam7$fitted.values)

#Prediction plotting

plot(L ~ DI, data=topology, 
     xlab= "Drought index (DI)", 
     ylab = "Number of links (L)" )

pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(gam7, newdata=pdat, na.rm=T, type="response", se.fit=TRUE)
pred

predframe <-data.frame(pdat, preds=pred$fit, se=pred$se.fit)
lines (predframe$preds ~ predframe$DI, col="red", lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5) # Adds lower SE line

###################################################################

####12. Edge density (LS)####

#Does Edge density (LS) vary significantly across foodwebs as a function of DI?
hist(topology$LS)
plot(topology$DI ~ topology$LS)


#Gamma GAM 
gam8<-gam(LS ~ s(DI, fx= FALSE, bs="tp", k= -1), 
          family=Gamma (link= "log"),
          data=topology) 
AIC(gam8)
plot(gam8) 
summary(gam8)   

# Prediction plotting
plot(LS ~ DI, data=topology, 
     xlab= "Drought index (DI)",
     ylab = "Edge density (LS)" )

pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(gam8, newdata=pdat,
              na.rm=T, 
              type="response",
              se.fit=TRUE)
pred

predframe <-data.frame(pdat, 
                       preds=pred$fit,
                       se=pred$se.fit)

lines (predframe$preds ~ predframe$DI,
       col="red", 
       lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5) # Adds lower SE line

###################################################################

####13. Connectance (C) ####

#Does Connectance (C) vary significantly across foodwebs as a function of DI?
hist(topology$C) # skewed to left data. Trial inverse gaussian and gamma.
plot(topology$C ~topology$DI) 

# GAM with inverse gaussian error family, and inverse link

gam9<-gam(C ~ s(DI, fx= FALSE, bs="tp", k= -1), 
          family=inverse.gaussian (link= "inverse"),
          data=topology) 
AIC(gam9) #-93.72052
plot(gam9)

# Final model
summary.gam(gam9)
AIC(gam9)
gam.check(gam9)

# Prediction plotting

plot(C ~ DI, data=topology, xlab= "Drought intensity (DI)", ylab = "Connectance (C)" )

pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(gam9, newdata=pdat, na.rm=T, type="response", se.fit=TRUE)
pred

predframe <-data.frame(pdat, preds=pred$fit, se=pred$se.fit)
lines (predframe$preds ~ predframe$DI, col="red", lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5) # Adds lower SE line

#####################################################################

####14. Fraction of basal nodes (B)####

#Does the fraction of basal nodes (B)  vary as a function of DI across channels?

# Beta error family given proportions bounded between 0 and 1.
library (betareg) # Load betareg package into library

# GLM with beta error family

model10 <- betareg(IAdjusted ~ as.numeric(DI),
                   link = "probit", 
                   link.phi="identity",
                   data = topology)
summary (model10) 

AIC(model10) 

#Check model assumptions
sresid<- resid (model10, type="pearson")
hist(sresid, freq=F)
lines(density (sresid, adjust=1))

qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)

plot(sresid~model10$fitted.values,
     pch=20,
     cex=2,
     cex.lab=1.5)

 # to get significance 
library(lmtest)

summary(model10)
lrtest(model10)

#####################################################################

####16.  Fraction of top nodes (t)####

#Does the fraction of top nodes (t) vary as a function of DI across channels?
hist(topology$t)
plot(topology$t~topology$DI) 

#GLM with beta error family

model11 <- betareg(t~ as.numeric(DI),
                   link = "log", 
                   link.phi="sqrt",
                   data = topology)
summary (model11)
AIC(model11) #  -46.69559

#Check model assumptions 
sresid<- resid (model11, type="pearson")
hist(sresid, freq=F)
lines(density (sresid, adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)

plot(sresid~model11$fitted.values,
     pch=20,
     cex=2,
     cex.lab=1.5)

summary (model11) # sig
lrtest(model11)

#####################################################################

####17. Food chain length (FCL)####

#Does Food chain length (FCL) vary as a function of DI across channels?
hist(topology$FCL) # assymetrical distribution, zero bounded continuous data suggests Gamma or inverse gaussian error families.

plot(topology$FCL~topology$DI)

# Gamma GAM 
gam12<-gam(FCL+0.0001 ~ s(DI,fx= FALSE, bs="tp",k= -1), 
          family=Gamma (link= "log"),
          data=topology) 
AIC(gam12)

# Final model validation

AIC(gam12)
plot(gam12)
summary(gam12)
summary.gam(gam12)
gam.check(gam12)

# Predication plotting
plot(FCL ~ DI, xlab = "Drought index (DI)", 
     ylab = "Food chain length (FCL)",
     data=topology)

pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(gam12, newdata=pdat, 
              na.rm=T, type="response",
              se.fit=TRUE)
pred

predframe <-data.frame(pdat, 
                       preds=pred$fit,
                       se=pred$se.fit)

lines (predframe$preds ~ predframe$DI, 
       col="red", 
       lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5) # Adds lower SE line

#####################################################################
####18. Trophic generality####

#Does the Trophic generality vary as a function of DI across channels?

hist(topology$InDeg)# in-degree represents generality, that is the number of
#preyconsumed by each predator.


plot(topology$InDeg~topology$DI) 

# GAM with Gamma error family

gam13<-gam(InDeg~s(as.numeric(DI), bs="tp", fx= FALSE, k= -1), 
           family= Gamma (link="identity"), 
           data= topology)
AIC(gam13) 
plot(gam13)
summary.gam(gam13)
gam.check (gam13)

# Predication plotting
plot(InDeg ~ DI, xlab = "Drought index (DI)", ylab = "InDegree", data=topology)

pdat<-expand.grid(DI=seq(0,1,0.1))
pred<-predict(gam13, 
              newdata=pdat,
              na.rm=T, 
              type="response",
              se.fit=TRUE)
pred

predframe <-data.frame(pdat, 
                       preds=pred$fit,
                       se=pred$se.fit)

lines (predframe$preds ~ predframe$DI, 
       col="red", 
       lwd =2)

predframe$upperse <- (predframe$preds + predframe$se)
lines (predframe$upperse ~ predframe$DI,
       lty = 2, 
       col = "blue", 
       lwd = 1.5) # Adds upper SE line

predframe$lowerse <- (predframe$preds - predframe$se)
lines (predframe$lowerse ~ predframe$DI,
       lty = 2, 
       col = "blue",
       lwd = 1.5) # Adds lower SE line





###################################

####19. Plot figures: Trophic networks####

# Plot circular food webs for all channels,  to give an overview of community
# complexity. Specify that the first node is always plotted at the bottom 
#(i.e. at the  6'o clock position), with the remaining nodes plotted 
# counter clockwise, to aid cross referencing between data sets.

library(cheddar) # load cheddar into memory as needed.


# Set own colours for trophic categories using the colour.spec function. 
#This specification will be applied across all 21 webs to aid comparability
#between channels:

colours() # gives an overview of colours available in R

colour.spec<-c('producer'='yellowgreen',
               'predator'='tomato2',
               'detritivore'='tan4', 
               'herbivore'='springgreen4',
               'omnivore'='magenta4',
               'parasitoids'= 'grey', 
               'fungivore'= 'orange',
               'hyperparasite'= 'steelblue',
               'parasite'= 'black')

#Plot circular webs for all channels, to give an overview of complexity. 

par(cex.axis=0.6, mar=c(1,1,1,1)) # Changes  margins of the plot (bottom, 
#left, top, right) to adjust position of legend as needed. 


#Channel 1
#Plot first point at '6 o'clock', and rest of nodes counter clockwise.
PlotCircularWeb(ch1_web.clean,
              origin.degrees=180, 
              clockwise=FALSE, 
              colour.by='cat', 
              colour.spec = colour.spec, 
              pch=19,
              highlight.nodes=NULL)

legend("topright", 
       legend=names(colour.spec), 
       cex= 0.6, 
       title = "Trophic category", 
       pch=19, 
       col=colour.spec,
       bty="n") # remove legend box

#Channel 2 
PlotCircularWeb(ch2_web.clean, 
                origin.degrees=180, 
                clockwise=FALSE, 
                colour.by="cat",
                colour.spec = colour.spec, 
                pch=19)

legend("topright", 
       legend=names(colour.spec), 
       cex=0.39, 
       title = "Trophic category", 
       pch=19,
       col=colour.spec, 
       bty="n") 

#Channel 3
PlotCircularWeb(ch3_web.clean, 
                origin.degrees=180,
                clockwise=FALSE,
                colour.by = "cat",
                colour.spec= colour.spec, 
                pch=19)

legend("topright", 
       legend=names(colour.spec),
       cex=0.39, 
       title = "Trophic category",
       pch=19, 
       col=colour.spec, 
       bty="n") 

#Channel 4
PlotCircularWeb(ch4_web.clean, 
                origin.degrees=180,
                clockwise=FALSE,
                colour.by = "cat",
                colour.spec=colour.spec,
                pch=19)

legend("topright",
       legend=names(colour.spec),
       cex=0.39, 
       title = "Trophic category",
       pch=19, 
       col=colour.spec,
       bty="n") 

#Channel 5
PlotCircularWeb(ch5_web.clean, 
                origin.degrees=180,
                clockwise=FALSE, 
                colour.by = "cat", 
                colour.spec=colour.spec, 
                pch=19)

legend("topright", 
       legend=names(colour.spec), 
       cex=0.39, 
       title = "Trophic category", 
       pch=19, 
       col=colour.spec,
       bty="n") 

#Channel 6
PlotCircularWeb(ch6_web.clean,
                origin.degrees=180,
                clockwise=FALSE,
                colour.by = "cat", 
                colour.spec=colour.spec,
                pch=19)

legend("topright", 
       legend=names(colour.spec),
       cex=0.39,
       title = "Trophic category", 
       pch=19, 
       col=colour.spec, 
       bty="n") 

#Channel 7
PlotCircularWeb(ch7_web.clean,
                origin.degrees=180,
                clockwise=FALSE, 
                colour.by = "cat",
                colour.spec=colour.spec,
                pch=19)

legend("topright",
       legend=names(colour.spec),
       cex=0.39,
       title = "Trophic category", 
       pch=19, 
       col=colour.spec,
       bty="n") 

#Channel 8
PlotCircularWeb(ch8_web.clean,
                origin.degrees=180,
                clockwise=FALSE,
                colour.by = "cat",
                colour.spec=colour.spec,
                pch=19)

legend("topright",
       legend=names(colour.spec), 
       cex=0.39, 
       title = "Trophic category", 
       pch=19, 
       col=colour.spec,
       bty="n") 

#Channel 9
PlotCircularWeb(ch9_web.clean,
                origin.degrees=180,
                clockwise=FALSE, 
                colour.by = "cat", 
                colour.spec=colour.spec, 
                pch=19)

legend("topright",
       legend=names(colour.spec), 
       cex=0.39,
       title = "Trophic category",
       pch=19, 
       col=colour.spec, 
       bty="n") 

#Channel 10
PlotCircularWeb(ch10_web.clean,
                origin.degrees=180, 
                clockwise=FALSE, 
                colour.by = "cat", 
                colour.spec=colour.spec,
                pch=19)

legend("topright",
       legend=names(colour.spec),
       cex=0.39, 
       title = "Trophic category", 
       pch=19, 
       col=colour.spec, 
       bty="n")

#Channel 11
PlotCircularWeb(ch11_web.clean,
                origin.degrees=180,
                clockwise=FALSE, 
                colour.by = "cat",
                colour.spec=colour.spec,
                pch=19)

legend("topright", 
       legend=names(colour.spec), 
       cex=0.39,
       title = "Trophic category", 
       pch=19,
       col=colour.spec, 
       bty="n") 

#Channel 12
PlotCircularWeb(ch12_web.clean,
                origin.degrees=180,
                clockwise=FALSE,
                colour.by = "cat",
                colour.spec=colour.spec,
                pch=19)

legend("topright", 
       legend=names(colour.spec), 
       cex=0.39, 
       title = "Trophic category", 
       pch=19,
       col=colour.spec,
       bty="n") 

#Channel 13
PlotCircularWeb(ch13_web.clean,
                origin.degrees=180,
                clockwise=FALSE, 
                colour.by = "cat",
                colour.spec=colour.spec, 
                pch=19)

legend("topright", 
       legend=names(colour.spec),
       cex=0.39, 
       title = "Trophic category",
       pch=19,
       col=colour.spec, 
       bty="n") 

#Channel 14
PlotCircularWeb(ch14_web.clean, 
                origin.degrees=180,
                clockwise=FALSE, 
                colour.by = "cat",
                colour.spec=colour.spec,
                pch=19)

legend("topright", 
       legend=names(colour.spec),
       cex=0.39, 
       title = "Trophic category",
       pch=19, 
       col=colour.spec,
       bty="n") 

#Channel 15
PlotCircularWeb(ch15_web.clean, 
                origin.degrees=180, 
                clockwise=FALSE,
                colour.by = "cat",
                colour.spec=colour.spec, 
                pch=19)

legend("topright",
       legend=names(colour.spec), 
       cex=0.39, 
       title = "Trophic category", 
       pch=19,
       col=colour.spec, 
       bty="n") 

#Channel 16
PlotCircularWeb(ch16_web.clean,
                origin.degrees=180, 
                clockwise=FALSE, 
                colour.by = "cat",
                colour.spec=colour.spec, 
                pch=19)

legend("topright", 
       legend=names(colour.spec), 
       cex=0.39, 
       title = "Trophic category", 
       pch=19, 
       col=colour.spec, 
       bty="n") 

#Channel 17
PlotCircularWeb(ch17_web.clean, 
                origin.degrees=180,
                clockwise=FALSE, 
                colour.by = "cat",
                colour.spec=colour.spec,
                pch=19)

legend("topright", 
       legend=names(colour.spec),
       cex=0.39,
       title = "Trophic category", 
       pch=19, 
       col=colour.spec, 
       bty="n") 

#Channel 18
PlotCircularWeb(ch18_web.clean,
                origin.degrees=180, 
                clockwise=FALSE,
                colour.by = "cat",
                colour.spec=colour.spec,
                pch=19)

legend("topright", 
       legend=names(colour.spec), 
       cex=0.39, 
       title = "Trophic category",
       pch=19, 
       col=colour.spec, 
       bty="n") 

#Channel 19
PlotCircularWeb(ch19_web.clean, 
                origin.degrees=180, 
                clockwise=FALSE,
                colour.by = "cat", 
                colour.spec=colour.spec, 
                pch=19)

legend("topright",
       legend=names(colour.spec), 
       cex=0.39, 
       title = "Trophic category",
       pch=19,
       col=colour.spec,
       bty="n") 

#Channel 20
PlotCircularWeb(ch20_web.clean, 
                origin.degrees=180, 
                clockwise=FALSE, 
                colour.by = "cat",
                colour.spec=colour.spec, 
                pch=19)

legend("topright", 
       legend=names(colour.spec),
       cex=0.39,
       title = "Trophic category",
       pch=19,
       col=colour.spec,
       bty="n") 

#Channel 21
PlotCircularWeb(ch21_web.clean, 
                origin.degrees=180, 
                clockwise=FALSE, 
                colour.by = "cat",
                colour.spec=colour.spec,
                pch=19)

legend("topright", 
       legend=names(colour.spec), 
       cex=0.39, 
       title = "Trophic category", 
       pch=19, 
       col=colour.spec,
       bty="n") 

dev.off() #restore margins


#########################END####################################################
citation("cheddar")
citation("igraph")
citation("rshape")
citation("dplyr")
citation("gtable")
citation("ggplot2")
