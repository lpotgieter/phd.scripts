library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(nortest)
library(Hmisc)
library(readODS)
setwd("your_directory")
master = read_ods("you master file with all info in different columns")
#Depth per site
ggplot(depth, aes(x=POS, y= MEAN_DEPTH)) +
geom_line() +
facet_wrap(~CHROM)
#Depth per individual
ggplot(master, aes(x=reorder(Isolate.Name, Row.Number), y =Mean.Depth, fill=Host)) +
geom_bar(stat = "identity") +
facet_wrap(~Country, scales = "free_x") +
theme(legend.position = "right",axis.text.x=element_text(angle=90,hjust=1)) +
scale_color_manual(labels = c("Sea beet", "Sugar beet", "Table beet"),
values = c("skyblue2", "darkolivegreen3", "grey"), aesthetics = "fill") +
labs(title="Mean Depth per Individual", x = "Isolate", y = "Mean Depth", color = "Host")
#plot number of snps per individual
ggplot(master, aes(x=reorder(Isolate.Name, Row.Number), y = Variant, fill = Host)) +
geom_bar(stat = "identity") +
facet_wrap(~Country, scales = "free_x") +
theme(legend.position = "right",axis.text.x=element_text(angle=90,hjust=1)) +
scale_color_manual(labels = c("Sea beet", "Sugar beet", "Table beet"),
values = c("skyblue2", "darkolivegreen3", "grey"), aesthetics = "fill") +
labs(title="Number of Variant Positions Per Individual", x = "Isolate", y = "Variant Positions", color
= "Host")

#plot number of missingness per individual
ggplot(master, aes(x=reorder(Isolate.Name, Row.Number), y = Percent.Missing, fill = Host)) +
geom_bar(stat = "identity") +
facet_wrap(~Country, scales = "free_x") +
theme(legend.position = "right",axis.text.x=element_text(angle=90,hjust=1)) +
scale_color_manual(labels = c("Sea beet", "Sugar beet", "Table beet"),
values = c("skyblue2", "darkolivegreen3", "grey"), aesthetics = "fill") +
labs(title="Fraction of Missing Sites Per Individual", x = "Isolate", y = "Fraction of Missing Sites",
color = "Host")
#Prepare environment for pi plot
uk = read.table("uk.recode.vcf.windowed.pi", header = T)
nd = read.table("nd.recode.vcf.windowed.pi", header = T)
ny = read.table("ny.recode.vcf.windowed.pi", header = T)
croatia = read.table("croatia.recode.vcf.windowed.pi", header = T)
italy = read.table("italy.recode.vcf.windowed.pi", header = T)
uk = mutate(uk, Location = "UK")
uk = mutate(uk, Host = "Sea beet")
nd = mutate(nd, Location = "ND")
nd = mutate(nd, Host = "Sugar beet")
ny = mutate(ny, Location = "NY")
ny = mutate(ny, Host = "Table beet")
croatia = mutate(croatia, Location = "Croatia")
croatia = mutate(croatia, Host = "Sea beet")
italy = mutate(italy, Location = "Italy")
italy = mutate(italy, Host = "Sugar beet")
combined = c('uk', 'nd', 'ny', 'croatia', "italy")
x.list <- lapply(combined, get)
combined.pi = do.call(rbind, x.list)
#Prepare environment for Tajima’s D plot
uk = read.table("uk.recode.vcf.Tajima.D", header = T)
nd = read.table("nd.recode.vcf.Tajima.D", header = T)
ny = read.table("ny.recode.vcf.Tajima.D", header = T)
croatia = read.table("croatia.recode.vcf.Tajima.D", header = T)
italy = read.table("italy.recode.vcf.Tajima.D", header = T)
uk = mutate(uk, Location = "UK")
uk = mutate(uk, Host = "Sea beet")
nd = mutate(nd, Location = "ND")
nd = mutate(nd, Host = "Sugar beet")
ny = mutate(ny, Location = "NY")
ny = mutate(ny, Host = "Table beet")
croatia = mutate(croatia, Location = "Croatia")
croatia = mutate(croatia, Host = "Sea beet")
italy = mutate(italy, Location = "Italy")
italy = mutate(italy, Host = "Sugar beet")
combined = c('uk', 'nd', 'ny', 'croatia', "italy")

x.list <- lapply(combined, get)
combined.tajd = do.call(rbind, x.list)
level_order <- c('UK', 'Croatia', 'Italy', "ND", "NY")
#Plot Tajima’s D and pi
ggplot(combined.pi, aes(x= factor(Location, level = level_order), y = PI, fill = Host)) +
geom_violin() +
ggtitle("Pi Average in 5 kb Windows") +
stat_summary(fun=median, geom="point", size=2, color="black")
stat_summary(fun.data="mean_sdl", mult=1,
geom="pointrange", width=0.2 )
ggplot(combined.tajd, aes(x= factor(Location, level = level_order), y = TajimaD, fill = Host)) +
geom_violin() +
ggtitle("Tajima's D Average in 5 kb Windows") +
xlab("Location") + ylab("Tajima's D") +
stat_summary(fun=median, geom="point", size=2, color="black")
stat_summary(fun.data="mean_sdl", mult=1,
geom="pointrange", width=0.2 )
