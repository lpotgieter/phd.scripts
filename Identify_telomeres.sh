# <<>><<>><<>><<>>
# | Needed inputs |
# <<>><<>><<>><<>>
bowtie_dir=”bowtie_directory”
work_dir=”your_directory”
polished_fasta="assembly_Czeina_default.fasta"
sample_name="Czeina"

telom_fn=${polished_fasta%.fasta}_telomeric_repeats
filtered_fn=${polished_fasta%.fasta}.filtered_on_cov.fasta
mt_ref_seq="Pseudocercospora_mori_mitochondrion.fasta" #From NCBI accession MG543071.1
# -----------------------------
# Detecting telomeric repeats
# -----------------------------
${bowtie_dir}bowtie-build ${work_dir}${polished_fasta} ${work_dir}${sample_name}
${bowtie_dir}bowtie \
${work_dir}${sample_name} \
-c CCCTAA \
--all -v 0 \
--threads 4 \
| sort -k 3 -nk 4 \
> ${work_dir}${telom_fn}.txt
#Conversion from bowtie output format to bed format
awk 'BEGIN {FS= "\t"; OFS="\t"} {print $3, $4, $4+length($5), $1, 111, $2}' ${work_dir}$
{telom_fn}.txt \
> ${work_dir}${telom_fn}.bed
#Merge repeats closer than the length of 1 (in case one repeat is mutated or has a sequencing
error, likely at the contig end since the coverage depth drops)
bedtools merge -i ${work_dir}${telom_fn}.bed -d 7 \
> ${work_dir}${telom_fn}_merged.bed
#Keep only blocks of more than 10 repeats (so longer than 60)
awk 'BEGIN {FS= "\t"; OFS="\t"} {diff=($3 - $2); if (diff > 60) print $1,$2,$3,diff} ' \
${work_dir}${telom_fn}_merged.bed \
> ${work_dir}${telom_fn}_merged_long.bed
