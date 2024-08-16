##################################################
##################################################
##################################################
#
# 170329
# Allele-specific ChIP analysis based on RNA-seq pipeline:
# '20170220_RNAseq_reprocessAllLibs_ASEstreamlined_workflow.sh'
#
# ##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="<projDirName>"
WORKDIR="<workDirName>"
REFDIR="<refdataDirName>"

##################################################
# Working folder
mkdir -p ${PROJDIR}/"${WORKDIR}"
cd  ${PROJDIR}/"${WORKDIR}"


##################################################
# Variables

CHROMS=( $( seq 1 19 ) X )
# CHROMS=( X )

cellType=patski

# libIDs generated below:
INPUTS=$( cat ${PROJDIR}/"${WORKDIR}"/libIDs ) 

# DATADIR also generated below:
DATADIR=$( cat ${PROJDIR}/"${WORKDIR}"/DATADIR )

SNPDIR=${REFDIR}/mm10pseudoSpretus






###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
# Parse SNP file for chr SNPs & generate BED files
# *** THIS WAS ALREADY DONE FOR RNA-seq ***
# *** So just link those files          ***

ls -ltrh ${RNASEQDATADIR}/chrSNPbeds
ln -sf ${RNASEQDATADIR}/chrSNPbeds ${THEDATADIR}

########
# Check:

ls -ltrh "${DATADIR}"/chrSNPbeds

for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

VALIDSNPCHRFILE="${SNPDIR}/spretus.${cellType4SNPs}.chr${CURRCHR}.snps.vcf.gz"
zcat "${VALIDSNPCHRFILE}" | wc -l 

OUTDIR="${DATADIR}/chrSNPbeds/chr${CURRCHR}snpBEDs"
zcat "${OUTDIR}"/spretus.${cellType4SNPs}.snps.chr${CURRCHR}.bed.gz | wc -l
zcat "${OUTDIR}"/spretus.${cellType4SNPs}.snps.chr${CURRCHR}.bed.200bp.gz | wc -l 

done









###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
# Must have already called peaks on black (i.e. unsegregated) PE data
# See 'ChIPseq_dataQCandMapping_workflow.sh': Peak calls on black (i.e. unsegregated) PE data
#
# Now partion the peaks into chromosome-specific files:

OUTPUTDIR="peaksPE"
ls -ltrh "${DATADIR}"/"${OUTPUTDIR}"
ls -ltrh "${DATADIR}"/"${OUTPUTDIR}"/*

for INPUT in $( echo ${INPUTS[@]} | tr ' ' '\n' | grep -v input ); do
echo "${INPUT}" 

for CURRCHR in ${CHROMS[@]}; do 
# CURRCHR=19
printf "\t${CURRCHR}\n"

zcat "${DATADIR}"/"${OUTPUTDIR}"/"${INPUT}".black6/"${INPUT}"_peaks.bed.gz | cut -f 1-3 | grep "^chr${CURRCHR}[[:blank:]]" | \
gzip -f > "${DATADIR}"/"${OUTPUTDIR}"/"${INPUT}".black6/"${INPUT}"_peaks.chr${CURRCHR}.bed.gz 

done
done


ls -ltrh "${DATADIR}"/"${OUTPUTDIR}"/*













#################################################
#################################################
#################################################
# Step 0: First must deduplicate sorted, filtered mapped reads.
# Not sufficient to just link the data

# *** Use sorted, filtered PE mapped reads! *** 
# See 'ChIPseq_dataQCandMapping_workflow.sh':
# Key command:
# samtools view -h -F 256 -F 4 -f 2 -q 10 "${INPUTFILE}" | samtools sort -m 64G -n -T "${SAMPLE}"."${ASSEMBLY}".tmp - > "${OUTPUTFILE}"

find "${DATADIR}"/mappedReadsPE -name "L1.filtered.bam"

mkdir "${DATADIR}"/sortedBAMs

#################################################
# Picard way slow... try samtools rmdup
# And sort & index files while at it.

mkdir -p "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/rmdupMappedReads
cd "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/rmdupMappedReads

for INPUT in ${INPUTS[@]}; do
echo "${INPUT}" 

jobFile="${INPUT}.rmdupMappedReads.job"
cat << EOF > "${jobFile}"
#!/bin/bash -x
#\$ -l h_rt=7:59:59
#\$ -l h_vmem=8G
#\$ -cwd
#\$ -j y

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
module load samtools/0.1.19

hostname

printf "\n\nstart: %s\n\n" "\$( date )"

ln -sf "${DATADIR}"/mappedReadsPE/"${INPUT}".black6/L1.filtered.bam "${DATADIR}"/sortedBAMs/"${INPUT}".filtered.nsorted.bam

samtools sort -m 20G "${DATADIR}"/sortedBAMs/"${INPUT}".filtered.nsorted.bam "${DATADIR}"/sortedBAMs/"${INPUT}".filtered.sorted 

samtools rmdup "${DATADIR}"/sortedBAMs/"${INPUT}".filtered.sorted.bam "${DATADIR}"/sortedBAMs/"${INPUT}".filtered.sorted.rmdup.bam

samtools index "${DATADIR}"/sortedBAMs/"${INPUT}".filtered.sorted.rmdup.bam

printf "\n\nfinish: %s\n" "\$( date )"

EOF
chmod 755 "${jobFile}"

done

#################################################
# Submit jobs

for INPUT in ${INPUTS[@]}; do
echo "${INPUT}"
jobFile="${INPUT}.rmdupMappedReads.job"
qsub "${jobFile}"
done 

#######
# Check

#######
# Check

ls -ltrh "${DATADIR}"/sortedBAMs

# Dup rates
cd "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/rmdupMappedReads

for INPUT in ${INPUTS[@]}; do
printf "${INPUT}\t" 
while read FILE; do 
grep "in library" "${FILE}" | sed 's/\[bam_rmdup_core\]//' | sed "s/ in library.*$//"
done < <( ls -1 "${INPUT}".rmdupMappedReads.job.o* )
done  > ${WORKDIR}_rmdup_rates.txt











#################################################
#################################################
#################################################
# Step 1. Subset mapped data to be chr and SNP specific
# Allows for parallel processing

mkdir "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/gensnpBAMsPerChrom
cd "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/gensnpBAMsPerChrom


#################################################
# Generate jobs

TALLY=0

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpBAMsPerChrom/genChr${CURRCHR}snpBAMs"
mkdir -p "${WDDIR}"
cd "${WDDIR}"

INPUTDIR="${DATADIR}/sortedBAMs"

OUTPUTDIR="${DATADIR}/snpBAMs/chr${CURRCHR}snpBAMs"
mkdir -p "${OUTPUTDIR}"

SNPFILE="${DATADIR}/chrSNPbeds/chr${CURRCHR}snpBEDs"/spretus.${cellType4SNPs}.snps.chr${CURRCHR}.bed.200bp.gz

for INPUT in ${INPUTS[@]}; do
BASEINPUT=$( basename "${INPUT}" ".fq.gz" )
echo "${INPUT}" "${BASEINPUT}"

TALLY=$(( TALLY + 1 ))
echo "${TALLY}"


jobFile="${BASEINPUT}.chr${CURRCHR}snpBAMs.job"
cat << EOF > "${jobFile}"
#!/bin/bash -x
#\$ -l h_rt=7:59:59
#\$ -l h_vmem=16G
#\$ -cwd
#\$ -j y

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
module load samtools/0.1.19

hostname

printf "\n\nstart: %s\n\n" "\$( date )"

samtools view -b -L "${SNPFILE}" "${INPUTDIR}"/${BASEINPUT}.filtered.sorted.rmdup.bam > "${OUTPUTDIR}"/${BASEINPUT}.chr${CURRCHR}snps.bam

samtools index "${OUTPUTDIR}"/${BASEINPUT}.chr${CURRCHR}snps.bam

printf "\n\nfinish: %s\n" "\$( date )"

EOF
chmod 755 "${jobFile}"
done

done # CURRCHR

cd ..


#################################################
# Submit jobs

for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpBAMsPerChrom/genChr${CURRCHR}snpBAMs"
#WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpKallistoBAMsPerChrom/genChr${CURRCHR}snpBAMs"
mkdir -p "${WDDIR}"
cd "${WDDIR}"

for INPUT in ${INPUTS[@]}; do
echo "${INPUT}"
BASEINPUT=$( basename "${INPUT/.fq.gz/}" )

jobFile="${BASEINPUT}.chr${CURRCHR}snpBAMs.job"
qsub "${jobFile}"
done 

done # CURRCHR

cd ..


#################################################
# Check

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpBAMsPerChrom/chr${CURRCHR}snpBAMs"
#WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpKallistoBAMsPerChrom/chr${CURRCHR}snpBAMs"

#ls -lrtrh "${WDDIR}"

OUTPUTDIR="${DATADIR}/snpBAMs/chr${CURRCHR}snpBAMs"
#OUTPUTDIR="${DATADIR}/snpKallistoBAMs/chr${CURRCHR}snpBAMs"

ls -lrtrh "${OUTPUTDIR}"

done | more  # CURRCHR 






















#################################################
#################################################
#################################################
# Step 2. Pileups at SNPs & annotate
# 150706
# *** WHERE THE WORK GETS DONE ***
# *** PARAMETERS: Q>=13 ***
#
# samtools mpileup for chr  SNPs reads using region lists
# *** NOTE: instead of full BAM file, use subselected BAM with reads in +/- 100bp around SNP
# *** ie. *.chr19.${INPUT}_SNPs_100bp.accepted_hits.sorted.bam) -> MUCH FASTER
#
# samtools version 1.2
# http://biobits.org/samtools_primer.html
# http://massgenomics.org/2012/03/5-things-to-know-about-samtools-mpileup.html
# https://www.biostars.org/p/8566/

GENOME="${HOME}/refdata/iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"

PERBASEQSCORETHRESH=13

MAXDEPTH=10000
# NOT SURE THIS IS IMPORTANT: 
# https://github.com/samtools/samtools/issues/29 
# http://sourceforge.net/p/samtools/mailman/message/26912287/

#################################################
# Generate jobs

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 

echo "${CURRCHR}"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
mkdir -p "${WDDIR}"
cd "${WDDIR}"

OUTDIR="${DATADIR}/snpPileups/chr${CURRCHR}snpPileups"
mkdir -p "${OUTDIR}"

PILEUPFILEPREFIX="Patski.chr${CURRCHR}.SNPpileup"

SNPFILE="${DATADIR}/chrSNPbeds/chr${CURRCHR}snpBEDs"/spretus.${cellType4SNPs}.snps.chr${CURRCHR}.bed.gz

# BAM file list
INPUTDIR="${DATADIR}/snpBAMs/chr${CURRCHR}snpBAMs"
ls -1 "${INPUTDIR}"/*.bam > chr${CURRCHR}snpBAM_list.txt

jobFile="chr${CURRCHR}snpPileups.job"

cat << EOF > "${jobFile}"
#!/bin/bash -x
#\$ -l h_rt=7:59:59
#\$ -l h_vmem=16G
#\$ -cwd
#\$ -j y

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
# Need newer version here
module load samtools/1.3

hostname

printf "\n\nstart: %s\n\n" "\$( date )"

#//# ###
#//# # VCF format
#//# samtools mpileup -v -I --max-depth "${MAXDEPTH}" --min-BQ "${PERBASEQSCORETHRESH}"  \
#//# -f "${GENOME}" \
#//# -l "${SNPFILE}" \
#//# -b chr${CURRCHR}snpBAM_list.txt \
#//# > ${PILEUPFILEPREFIX}.vcf
#//# 
#//# bcftools view -Sc "${PILEUPFILEPREFIX}.vcf" > "${PILEUPFILEPREFIX}.vcf.tsv"

###
# pileup format

cd "${OUTDIR}"

samtools mpileup --max-depth ${MAXDEPTH} --min-BQ ${PERBASEQSCORETHRESH} \
-f "${GENOME}" \
-l "${SNPFILE}" \
-b "${WDDIR}"/chr${CURRCHR}snpBAM_list.txt \
> "${PILEUPFILEPREFIX}.tsv"

gzip -f "${PILEUPFILEPREFIX}.tsv"

# DEBUG 
# zcat "${PILEUPFILEPREFIX}.tsv.gz" | awk '{print NF; print $0; }' | head 
# zcat "${PILEUPFILEPREFIX}.tsv.gz" | awk 'BEGIN{FS="\t"} {print NF; print $0; }' | head 

printf "\n\nfinish: %s\n" "\$( date )"

EOF
chmod 755 "${jobFile}"

done # CURRCHR

cd ..



#################################################
# Submit jobs

for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
#mkdir "${WDDIR}"
cd "${WDDIR}"

#OUTDIR="${DATADIR}/snpPileups/chr${CURRCHR}snpPileups"
#OUTDIR="${DATADIR}/snpKallistoPileups/chr${CURRCHR}snpPileups"
# mkdir "${OUTDIR}"
#cd "${OUTDIR}"

jobFile="chr${CURRCHR}snpPileups.job"
qsub "${jobFile}"

done # CURRCHR

cd ..



#################################################
# Check

for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

PILEUPFILEPREFIX="Patski.chr${CURRCHR}.SNPpileup"
#ls -ltrh "${DATADIR}/snpPileups/chr${CURRCHR}snpPileups/${PILEUPFILEPREFIX}.tsv.gz"
du -sh "${DATADIR}/snpPileups/chr${CURRCHR}snpPileups/${PILEUPFILEPREFIX}.tsv.gz"

done # CURRCHR




#################################################
# Annotate & cleanup

DATADIR=$( cat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR )
cd "${DATADIR}/snpPileups"

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"

OUTDIR="${DATADIR}/snpPileups/chr${CURRCHR}snpPileups"
cd "${OUTDIR}"

PILEUPFILEPREFIX="Patski.chr${CURRCHR}.SNPpileup"

SNPFILE="${DATADIR}/chrSNPbeds/chr${CURRCHR}snpBEDs"/spretus.${cellType4SNPs}.snps.chr${CURRCHR}.bed.gz

TMPDIR=$(mktemp -dt "GBDIR.XXXXXXX")

###
# Add SNP info: Not the ID, but the SNP (ref/alt)
awk 'BEGIN{FS="\t"} NR==FNR {locus=sprintf("%s_%s", $1, $3); meHash[locus]=$4} \
FNR<NR { locus=sprintf("%s_%s", $1, $2); if(locus in meHash){ printf("%s\t%s\n", meHash[locus], $0) } else
{ printf("%s\t%s\n", "NA", $0 ) }}' \
<( zcat "${SNPFILE}" ) <( zcat "${PILEUPFILEPREFIX}.tsv.gz" ) \
> "${TMPDIR}/${PILEUPFILEPREFIX}_summary.tsv.tmp2"

# cat "${TMPDIR}/${PILEUPFILEPREFIX}_summary.tsv.tmp2" | awk 'BEGIN{FS="\t"} {print NF; print $0; }' | head 

###
# Adjusted coverage (Ignore '>' '<' & '*')
# https://www.biostars.org/p/8566/
# Ref to alt ratio
# Alt to ref ratio

# *** 160505 Updated to consider only ref and alt SNPs ***

# cat "${TMPDIR}/${PILEUPFILEPREFIX}_summary.tsv.tmp2" | \
# awk 'BEGIN{FS="\t"} {line=sprintf("%s\t%s\t%s\t%s", $1, $2, $3, $4); \
# for (i=5; i<=NF; i++){ if((i-5) % 3 == 0){ 
# adjCov=$i; refHits=0; altHits=0; refProp=0; seq=$(i+1); cleanSeq=seq; \
# if(adjCov>0 && (seq~/>/ || seq~/</ || seq~/*/) ) { adjCov-=split(seq,a,"<")-1; adjCov-=split(seq,a,">")-1; adjCov-=split(seq,a,"*")-1; \
# gsub(">", "", cleanSeq); gsub("<", "", cleanSeq); gsub("\*", "", cleanSeq); } \
# if(adjCov>0) { refHits+=split(seq,a,".")-1; refHits+=split(seq,a,",")-1; altHits+=gsub("[AaCcGgTt]", "x",seq); refProp=refHits/adjCov; altProp=altHits/adjCov } \
# line=sprintf("%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%s\t%s\t%s", line, $i, adjCov, refHits, altHits, refProp, altProp, cleanSeq, $(i+1), $(i+2)); \
# }} \
# printf("%s\n", line)}' \
# > "${TMPDIR}/${PILEUPFILEPREFIX}_summary.tsv.tmp3"

# \044 = $
# \052 = *
# \056 = .
# \136 = ^

# For testing:
#cat "${TMPDIR}/${PILEUPFILEPREFIX}_summary.tsv.tmp2" | head -n 1000 | 

cat "${TMPDIR}/${PILEUPFILEPREFIX}_summary.tsv.tmp2" | \
awk 'BEGIN{FS="\t"} { \
altSNP=substr($1, 2, 1); altSNPlower=tolower(altSNP); \
line=sprintf("%s\t%s\t%s\t%s\t%s\t%s", $1, altSNP, altSNPlower, $2, $3, $4); \
\
for (i=5; i<=NF; i++){ if((i-5) % 3 == 0){ 
cov=$i; refHits=0; altHits=0; refProp=0; altProp=0; seq=$(i+1); cleanSeq=seq; \
\
if(cov>0) { \
gsub(">", "", cleanSeq); \
gsub("<", "", cleanSeq); \
gsub("\\\052", "", cleanSeq); \
gsub("\\\136.", "", cleanSeq); \
gsub("\\\044", "", cleanSeq); \
cleanSeqAR = cleanSeq; \
refHits+=gsub("\\\056", "R", cleanSeqAR); \
refHits+=gsub(",", "r", cleanSeqAR); \
altHits+=gsub(altSNP, "A", cleanSeqAR); \
altHits+=gsub(altSNPlower, "a", cleanSeqAR); \
adjCov=refHits+altHits } \
\
if(adjCov>0) { refProp=refHits/adjCov; altProp=altHits/adjCov } \
else { cleanSeq="?"; cleanSeqAR="?" }\
\
line=sprintf("%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%s\t%s\t%s\t%s", line, cov, adjCov, refHits, altHits, refProp, altProp, cleanSeq, cleanSeqAR, $(i+1), $(i+2)); \
}} \
printf("%s\n", line)}' \
> "${TMPDIR}/${PILEUPFILEPREFIX}_summary.tsv.tmp3"

# more "${TMPDIR}/${PILEUPFILEPREFIX}_summary.tsv.tmp3"


###
# pileup header
# See 'chr${CURRCHR}snpBAM_list.txt' for DATASETs
HEADER="SNPid\tRefSNP\tAltSNP\tchr\tpos\trefbase"
while read LINE; do
DATASET=$( echo "${LINE}" | sed 's#^.*/\(.*\)\.chr.*#\1#' )
HEADER="${HEADER}\t${DATASET}_cov\t${DATASET}_trueCov\t${DATASET}_refCov\t${DATASET}_altCov\t${DATASET}_refProp\t${DATASET}_altProp\t${DATASET}_cleancalls\t${DATASET}_cleancallsAR\t${DATASET}_basecalls\t${DATASET}_quals"
done < <( cat "${WDDIR}"/chr${CURRCHR}snpBAM_list.txt )

###
# Add header, drop basecalls & quals column
printf "${HEADER}\n" | cat - <( cat "${TMPDIR}/${PILEUPFILEPREFIX}_summary.tsv.tmp3" ) | \
awk 'BEGIN{FS="\t"} \
\
NR==1{ for (i=1; i<=NF; i++){ \
if($i!~/quals$/ && $i!~/_basecalls$/ && $i!~/_cov$/){ colHash[i]=$i } \
} }; \
\
{line=$1; \
for (i=2; i<=NF; i++){ \
if(i in colHash){ line=sprintf("%s\t%s", line, $i) } \
}; \
printf("%s\n", line)}' \
| gzip > "${PILEUPFILEPREFIX}.annotated.tsv.gz" 

rm -rf "${TMPDIR}"

done # CURRCHR



#########
# Check # of SNPs!
# DON'T KNOW CAUSE BUT JUST FILTER THOSE PROBLEMATIC SNPS OUT FOR NOW

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo ""
echo "${CURRCHR}"

#//# WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/"snpPileups/chr${CURRCHR}snpPileups"
#//# mkdir "${WDDIR}"
#//# cd "${WDDIR}"
OUTDIR="${DATADIR}/snpPileups/chr${CURRCHR}snpPileups"
#OUTDIR="${DATADIR}/snpKallistoPileups/chr${CURRCHR}snpPileups"
# mkdir "${OUTDIR}"
cd "${OUTDIR}"

SNPPDATAFILE="Patski.chr${CURRCHR}.SNPpileup.tsv.gz"
zcat "${SNPPDATAFILE}" | awk 'BEGIN{IFS="\t"} {print NF}' 2>&1 | head
echo "---"
zcat "${SNPPDATAFILE}" | awk 'BEGIN{IFS="\t"} {print NF}' | grep 6 | wc -l
zcat "${SNPPDATAFILE}" | awk 'BEGIN{IFS="\t"} {print NF}' | grep -v 6 | wc -l
echo "---"

SNPPDATAFILE="Patski.chr${CURRCHR}.SNPpileup.annotated.tsv.gz"
zcat "${SNPPDATAFILE}" | awk 'BEGIN{IFS="\t"} {print NF}' 2>&1 | head
echo "---"
zcat "${SNPPDATAFILE}" | awk 'BEGIN{IFS="\t"} {print NF}' | grep 13 | wc -l
zcat "${SNPPDATAFILE}" | awk 'BEGIN{IFS="\t"} {print NF}' | grep -v 13 | wc -l
echo "---"

done 


###
# Clean up

NFIELDS=13

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

#//# WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/"snpPileups/chr${CURRCHR}snpPileups"
#//# mkdir "${WDDIR}"
#//# cd "${WDDIR}"
OUTDIR="${DATADIR}/snpPileups/chr${CURRCHR}snpPileups"
#OUTDIR="${DATADIR}/snpKallistoPileups/chr${CURRCHR}snpPileups"
# mkdir "${OUTDIR}"
cd "${OUTDIR}"

SNPPDATAFILE="Patski.chr${CURRCHR}.SNPpileup.annotated.tsv.gz"

cp -pr "${SNPPDATAFILE}" "${SNPPDATAFILE}.bak"
zcat "${SNPPDATAFILE}.bak" | awk -v numFields=${NFIELDS} 'NF==numFields{print $0}' | gzip > "${SNPPDATAFILE}"

done 


###
# Check counts

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

OUTDIR="${DATADIR}/snpPileups/chr${CURRCHR}snpPileups"
cd "${OUTDIR}"

SNPPDATAFILE="Patski.chr${CURRCHR}.SNPpileup.annotated.tsv.gz"

zcat "${SNPPDATAFILE}.bak" | wc -l
zcat "${SNPPDATAFILE}" | wc -l

done 





###
# Restore backups, if necessary.

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

OUTDIR="${DATADIR}/snpPileups/chr${CURRCHR}snpPileups"
cd "${OUTDIR}"

SNPPDATAFILE="Patski.chr${CURRCHR}.SNPpileup.annotated.tsv.gz"

cp -pr "${SNPPDATAFILE}.bak" "${SNPPDATAFILE}"

done 









#################################################
#################################################
#################################################
# 160510
# Step 3. Alt and Ref SNP coverage for each chrom: bedgraphs from mpileup tsv

cd "${DATADIR}/snpPileups"

################
# Per chromosome bedgraphs

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"

OUTDIR="${DATADIR}/snpPileups/chr${CURRCHR}snpPileups"
cd "${OUTDIR}"

SNPPDATAFILE="Patski.chr${CURRCHR}.SNPpileup.annotated.tsv.gz"

DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#' ) )
TALLY=-1
for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"

TALLY=$(( TALLY + 1 ))
ALTDATCOL=$(( (TALLY * 7) + 9 ))
REFDATCOL=$(( (TALLY * 7) + 8 ))
echo "${DATASET} ${TALLY} ${ALTDATCOL} ${REFDATCOL}"

###
# Bedgraphs for alt coverage for each line
zcat "${SNPPDATAFILE}" | awk -v datcol=${ALTDATCOL} 'NR>1 {printf("%s\t%s\t%s\t%d\n", $4, $5-1, $5, $datcol)}' | gzip > "${DATASET}".chr${CURRCHR}.altCov.bedgraph.gz

###
# Bedgraphs for ref coverage for each line
zcat "${SNPPDATAFILE}" | awk -v datcol=${REFDATCOL} 'NR>1 {printf("%s\t%s\t%s\t%d\n", $4, $5-1, $5, $datcol)}' | gzip > "${DATASET}".chr${CURRCHR}.refCov.bedgraph.gz

done # DATASET

done # CURRCHR









#################################################
# Aggregation into genome-wide maps

OUTDIR="${DATADIR}/snpPileups"
cd "${OUTDIR}"

############
# Initialize files

cat /dev/null > Patski.SNPpileup.annotated.tsv

CURRCHR="1"
WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#' ) )
for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"
echo "${DATASET}"

cat /dev/null > "${DATASET}".altCov.bedgraph
cat /dev/null > "${DATASET}".refCov.bedgraph

done # DATASET



############
# Populate files

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

CHROMDIR="${OUTDIR}"/"chr${CURRCHR}snpPileups"

SNPPDATAFILE="Patski.chr${CURRCHR}.SNPpileup.annotated.tsv.gz"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"

zcat "${CHROMDIR}"/Patski.chr${CURRCHR}.SNPpileup.annotated.tsv.gz >> Patski.SNPpileup.annotated.tsv

DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#' ) )
for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"
echo "${DATASET}"

###
# Bedgraphs for alt coverage for each line
zcat "${CHROMDIR}"/"${DATASET}".chr${CURRCHR}.altCov.bedgraph.gz >> "${DATASET}".altCov.bedgraph

###
# Bedgraphs for ref coverage for each line
zcat "${CHROMDIR}"/"${DATASET}".chr${CURRCHR}.refCov.bedgraph.gz >> "${DATASET}".refCov.bedgraph

done # DATASET

done # CURRCHR



############
# Compress files

gzip -f Patski.SNPpileup.annotated.tsv

CURRCHR="1"
WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#' ) )
for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"
echo "${DATASET}"

gzip -f "${DATASET}".altCov.bedgraph
gzip -f "${DATASET}".refCov.bedgraph

done # DATASET

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#################################################
#################################################
#################################################
# 160505
# Step 4: SNP coverage totals over diploid peaks

#################################################
# BEDMAP: genes overlapped by SNPs:
# Get average ratio/proportion of SNPs that fall within genes
# https://bedops.readthedocs.io/en/latest/content/reference/statistics/bedmap.html#bedmap

module load bedops/2.4.14

DATADIR=$( cat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR )

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
echo "${CURRCHR}"

OUTDIR="${DATADIR}/snpPileups/chr${CURRCHR}snpPileups/peakLevelAnalysis.chr${CURRCHR}.AllReads"
mkdir "${OUTDIR}"
cd "${OUTDIR}"

PILEUPFILEPREFIX="Patski.chr${CURRCHR}.SNPpileup"

#REFGENEFILE="$HOME/refdata/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/refFlat.chr${CURRCHR}.bed.gz"
PEAKDIR="peaksPE"


#########
# Prep input data

TMPDIR=$(mktemp -dt "GBDIR.XXXXXXX")

# Must unzip data & insert 4th colum
# zcat "${REFGENEFILE}" > "${TMPDIR}"/refFlat.chr${CURRCHR}.bed

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#' | grep -v input ) )  # No input 
for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"
echo "${DATASET}"

CHRPEAKFILE="${DATADIR}"/"${PEAKDIR}"/"${DATASET}".black6/"${DATASET}"_peaks.chr${CURRCHR}.bed.gz 
zcat "${CHRPEAKFILE}" > "${TMPDIR}"/"${DATASET}"_peaks.chr${CURRCHR}.bed

zcat ../"${DATASET}".chr${CURRCHR}.altCov.bedgraph.gz | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "X", $4}' > "${TMPDIR}"/"${DATASET}".chr${CURRCHR}.altCov.bed
zcat ../"${DATASET}".chr${CURRCHR}.refCov.bedgraph.gz | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "X", $4}' > "${TMPDIR}"/"${DATASET}".chr${CURRCHR}.refCov.bed

# sum, mean & std dev
bedmap --echo --sum --mean --stdev --count --delim "\t" "${TMPDIR}"/"${DATASET}"_peaks.chr${CURRCHR}.bed "${TMPDIR}"/"${DATASET}".chr${CURRCHR}.altCov.bed > "${DATASET}".chr${CURRCHR}.altCov.sumOverPeaks.bed
bedmap --echo --sum --mean --stdev --count --delim "\t" "${TMPDIR}"/"${DATASET}"_peaks.chr${CURRCHR}.bed "${TMPDIR}"/"${DATASET}".chr${CURRCHR}.refCov.bed > "${DATASET}".chr${CURRCHR}.refCov.sumOverPeaks.bed

# See '20170220_RNAseq_reprocessAllLibs_ASEstreamlined_workflow.sh' for bedops testing

done # DATASET 

rm -rf "${TMPDIR}"


#########
# Generate BEDs of proportions per peak

#MINCOV=10
MINCOV=5
 

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#' | grep -v input ) )
for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"
echo "${DATASET}"

# Average proportions - BED FORMAT
paste "${DATASET}".chr${CURRCHR}.altCov.sumOverPeaks.bed "${DATASET}".chr${CURRCHR}.refCov.sumOverPeaks.bed | awk -v minCov=${MINCOV} 'BEGIN{OFS="\t"}; ($7+$17) >= minCov {print $1, $2, $3, $4, $5, $6, $7/($7+$17), ($7+$17)}' > "${DATASET}".chr${CURRCHR}.minCov${MINCOV}.altProp.averageOverPeaks.bed

paste "${DATASET}".chr${CURRCHR}.altCov.sumOverPeaks.bed "${DATASET}".chr${CURRCHR}.refCov.sumOverPeaks.bed | awk -v minCov=${MINCOV} 'BEGIN{OFS="\t"}; ($7+$17) >= minCov {print $1, $2, $3, $7/($7+$17)}' > "${DATASET}".chr${CURRCHR}.minCov${MINCOV}.altProp.averageOverPeaks.bedgraph

done # DATASET 
 
done # CURRCHR
 
cd ../..
 
 












#################################################
# Aggregation into genome-wide coverage counts
#
# ***************************************************
# *** But need a separate table for each data set ***
# *** b/c number of peaks is different for each   ***
# ***************************************************

OUTDIR="${DATADIR}/snpPileups"
cd "${OUTDIR}"


############
# Initialize files

CURRCHR="1"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#' | grep -v input ) )
for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"
echo "${DATASET}"

printf "chr\tstart\tend" > "${DATASET}".SNPcoverage.sumOverPeaks.tsv
printf "\t${DATASET}.altCov.sumOverPeaks\t${DATASET}.altCov.meanOverPeaks\t${DATASET}.altCov.sdOverPeaks\t${DATASET}.altCov.snpCountOverPeaks" >> "${DATASET}".SNPcoverage.sumOverPeaks.tsv
printf "\t${DATASET}.refCov.sumOverPeaks\t${DATASET}.refCov.meanOverPeaks\t${DATASET}.refCov.sdOverPeaks\t${DATASET}.refCov.snpCountOverPeaks" >> "${DATASET}".SNPcoverage.sumOverPeaks.tsv

cat /dev/null > "${DATASET}".altCov.sumOverPeaks.bed
cat /dev/null > "${DATASET}".refCov.sumOverPeaks.bed

printf "\n" >> "${DATASET}".SNPcoverage.sumOverPeaks.tsv

done # DATASET







############
# Populate files

######
# SNPcoverage table

TMPDIR=$(mktemp -dt "GBDIR.XXXXXXX")

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
# CURRCHR="10"

echo "${CURRCHR}"

CHROMDIR="${OUTDIR}"/"chr${CURRCHR}snpPileups/peakLevelAnalysis.chr${CURRCHR}.AllReads"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#' | grep -v input ) )

for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"
echo "${DATASET}"

# cat "${CHROMDIR}"/PatskiDel1.chr"${CURRCHR}".altCov.sumOverPeaks.bed | cut -f 1-6 > "${TMPDIR}"/"${DATASET}".SNPcoverage.sumOverPeaks.tmp.tsv
cat "${CHROMDIR}"/"${DATASET}".chr"${CURRCHR}".altCov.sumOverPeaks.bed | cut -f 1-3 > "${TMPDIR}"/"${DATASET}".SNPcoverage.sumOverPeaks.tmp.tsv

paste \
<( cat "${TMPDIR}"/"${DATASET}".SNPcoverage.sumOverPeaks.tmp.tsv ) \
<( cat "${CHROMDIR}"/"${DATASET}".chr"${CURRCHR}".altCov.sumOverPeaks.bed | cut -f 4-7 ) \
<( cat "${CHROMDIR}"/"${DATASET}".chr"${CURRCHR}".refCov.sumOverPeaks.bed | cut -f 4-7 ) \
> "${TMPDIR}"/Patski.SNPcoverage.sumOverPeaks.tmp2.tsv

mv "${TMPDIR}"/Patski.SNPcoverage.sumOverPeaks.tmp2.tsv "${TMPDIR}"/"${DATASET}".SNPcoverage.sumOverPeaks.tmp.tsv

cat "${TMPDIR}"/"${DATASET}".SNPcoverage.sumOverPeaks.tmp.tsv >> "${DATASET}".SNPcoverage.sumOverPeaks.tsv

done # DATASET

done # CURRCHR 

rm -rf "${TMPDIR}"



######
# SNPcoverage bed files

# For each chrom
for CURRCHR in ${CHROMS[@]}; do 
# CURRCHR="10"

echo "${CURRCHR}"

CHROMDIR="${OUTDIR}"/"chr${CURRCHR}snpPileups/peakLevelAnalysis.chr${CURRCHR}.AllReads"

WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#'  | grep -v input ) )

for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"
echo "${DATASET}"

###
# Bedgraphs for alt coverage for each line
cat "${CHROMDIR}"/"${DATASET}".chr${CURRCHR}.altCov.sumOverPeaks.bed >> "${DATASET}".altCov.sumOverPeaks.bed

###
# Bedgraphs for ref coverage for each line
cat "${CHROMDIR}"/"${DATASET}".chr${CURRCHR}.refCov.sumOverPeaks.bed >> "${DATASET}".refCov.sumOverPeaks.bed

done # DATASET

done # CURRCHR



############
# Compress files

CURRCHR="1"
WDDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}/gensnpPileups/genchr${CURRCHR}snpPileups"
DATASETLIST=( $( cat "${WDDIR}/chr${CURRCHR}snpBAM_list.txt" | sed 's#^.*/\(.*\)\.chr.*#\1#' | grep -v input ) )

for DATASET in ${DATASETLIST[@]}; do 
# DATASET="PatskiWT"
echo "${DATASET}"

gzip -f "${DATASET}".SNPcoverage.sumOverPeaks.tsv
gzip -f "${DATASET}".altCov.sumOverPeaks.bed
gzip -f "${DATASET}".refCov.sumOverPeaks.bed

done # DATASET





############
# Check

OUTDIR="${DATADIR}/snpPileups"
#OUTDIR="${DATADIR}/snpKallistoPileups"
# mkdir "${OUTDIR}"
cd "${OUTDIR}"

find . -maxdepth 1 -name "*.bedgraph.gz" -print0 | xargs -0 ls -ltrh
