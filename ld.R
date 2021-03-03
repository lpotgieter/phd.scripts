#from https://www.biostars.org/p/300381/
library(dplyr)
library(stringr)
library(ggplot2)
dfr <-read.delim("/media/lizel/82AEF34EAEF3396D/chapter2/2_vcf/34_rehh/chr/
uk.recode.vcf.chr2.recode.vcf.biallelic.nomissing.recode.vcf.75k.recode.vcf.summary",sep="",hea
der=F,check.names=F,stringsAsFactors=F)
colnames(dfr) <- c("dist","rsq")
dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=10000))
dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-
e+.]+")),
end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
mid=start+((end-start)/2))
ggplot(dfr1, aes(x=start,y=mean))+
geom_line()
