library(ggplot2)
library(tidyverse)
library(seqinr)
dir_name = "your_directory"
file_names=c("your_polished_assembly")
polishing_algo="arrow"
suffix_cov="-alignment_summary.gff"
#Note: in the new format, the gff gives several information about the coverage (based on the
python code found at https://github.com/PacificBiosciences/pbreports/blob/master/pbreports/
report/summarize_coverage/summarize_coverage.py):
#cov=(min_cov, median_cov, max_cov)))
#cov2=(mean_cov, sd_cov)))
#Because for the old version of SMRT assembly, we had only the average, I will use the average
with the new version too.
suffix_fasta=".fasta"
for ( i in file_names) {
#File names
fasta_input_name = paste(dir_name, i, suffix_fasta, sep = "")
fasta_output_name = paste(dir_name, i, ".filtered_on_cov", suffix_fasta, sep = "")
cov_file_name = paste(dir_name, i, suffix_cov, sep = "")
print(fasta_input_name)
#Get the average coverage data from the gff file
names_T = c("Unitig", "Source", "Feature", "Start",
"End", "Score", "Strand", "Phase", "Attributes")
T = read_tsv(cov_file_name, col_names = names_T, comment = "#")
T_sep = separate(data = T, col = Attributes, into = c("Cov", "Cov2", "Gaps"), sep=";")
T_sep$Cov2 = str_remove(T_sep$Cov2, "cov2=")
T = separate(data = T_sep, col = Cov2, into = c("Mean_coverage", "Sd_cov"), sep=",", convert =
TRUE)
#Extracting the length of the unitig and the average coverage
long = aggregate(T$End, by = list(T$Unitig), FUN = max)
mean = aggregate(T$Mean_coverage, by = list(T$Unitig), FUN = mean)
E = merge(long, mean, by= "Group.1")
names(E) = c("Unitig", "Length", "Mean")
med_value = mean(rep(E$Mean, times=E$Length))
med_max = med_value * 1.5
med_min = med_value / 1.5
kept = E[E$Mean > med_min & E$Mean < med_max,]
E$Pass = ifelse(E$Unitig %in% kept$Unitig, "pass", "fail")
T$Pass = ifelse(T$Unitig %in% kept$Unitig, "pass", "fail")
print(dim(kept))
print(dim(E))

len_kept=sum(kept$Length)
len_all = sum(E$Length)
prop_kept = round(len_kept*100/len_all, 2)
#Now, that we have the thresholds and results, we will print the results in graphs and tables
update_geom_defaults("point", list(colour = NULL))
value_plot = ggplot(data = T, aes(x=Unitig, y = Mean_coverage, col = Pass, fill = Pass))
value_plot + geom_boxplot(outlier.alpha = 0.5, outlier.size = 0.6) +
geom_hline(aes(yintercept=med_min)) +
geom_hline(aes(yintercept=med_max)) +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust= 1)) +
labs(title = paste("Coverage filter on PacBio assembly for", i),
subtitle = paste("The total assembly goes from", len_all,
"to", len_kept, "bp (", prop_kept," kept) and from",
nrow(E), "contigs to", nrow(kept),"."))
ggsave(paste(dir_name, i, "_filter.png", sep=""), width = 11, height = 8)
write.table(E, file = paste(dir_name, i, "_filter.tab", sep=""),
quote = FALSE, sep = "\t", row.names = FALSE)
write.table(kept, file = paste(dir_name, i, "_pass_filter.tab", sep=""),
quote = FALSE, sep = "\t", row.names = FALSE)
#Let's filter the data, if the fasta files are there
if (file.exists(fasta_input_name)) {
fastafile = read.fasta(file = fasta_input_name, seqtype = "DNA",
forceDNAtolower = FALSE, as.string = TRUE, set.attributes = FALSE)
f<-fastafile[names(fastafile) %in% paste(kept$Unitig, "|", polishing_algo, sep = "")]
write.fasta(f, names(f), file.out=fasta_output_name)
}else{
print(paste("WARNING. I can't find the fasta file ", fasta_input_name,
" so I will ignore the filtereing part", sep =""))
}
}
