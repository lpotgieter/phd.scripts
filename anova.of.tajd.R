library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(nortest)
library(Hmisc)
setwd("your_directory")
#Prepare environment for pi
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
#Prepare environment for Tajima’s D
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
group_by(combined.tajd, Location) %>%
summarise(
count = n(),
mean = mean(TajimaD, na.rm = TRUE),
sd = sd(TajimaD, na.rm = TRUE)
)
#HOST
#check outliers
temp = combined.tajd %>%
group_by(Host) %>%
identify_outliers(TajimaD)
#build linear model
model <- lm(TajimaD ~ Host, data = combined.tajd)
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
ad.test(combined.tajd$TajimaD)

ggqqplot(combined.tajd, "TajimaD", facet.by = "Host")
plot(model, 1)
res.aov <- combined.tajd %>% anova_test(TajimaD ~ Host)
res.aov
pwc <- combined.tajd %>% tukey_hsd(TajimaD ~ Host)
pwc
pwc <- pwc %>% add_xy_position(x = "Host")
ggboxplot(combined.tajd, x = "Host", y = "TajimaD") +
stat_pvalue_manual(pwc, hide.ns = TRUE) +
labs(
subtitle = get_test_label(res.aov, detailed = TRUE),
caption = get_pwc_label(pwc)
)
#LOCATION
#check outliers
temp = combined.tajd %>%
group_by(Location) %>%
identify_outliers(TajimaD)
#build linear model
model <- lm(TajimaD ~ Location, data = combined.tajd)
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
ad.test(combined.tajd$TajimaD)
ggqqplot(combined.tajd, "TajimaD", facet.by = "Location")
plot(model, 1)
res.aov <- combined.tajd %>% anova_test(TajimaD ~ Location)
res.aov
pwc <- combined.tajd %>% tukey_hsd(TajimaD ~ Location)
pwc
pwc <- pwc %>% add_xy_position(x = "Location")
ggboxplot(combined.tajd, x = "Location", y = "TajimaD") +
stat_pvalue_manual(pwc, hide.ns = TRUE) +
labs(
subtitle = get_test_label(res.aov, detailed = TRUE),
caption = get_pwc_label(pwc)
)

res.aov2 <- aov(TajimaD ~ Location + Host, data = combined.tajd)
summary(res.aov2)
