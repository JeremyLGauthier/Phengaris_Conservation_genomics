################################################
#title           :Phengaris population genomics scripts
#description     :This script combine all scripts for the analyses undertaken for the paper entitled: 
#"Conservation strategy insights for three protected Phengaris butterflies combining population genomics and landscape analyses"
#author		 : Jeremy Gauthier
#email	: jeremy.gauthier@vd.ch
#date            :31.03.2025
################################################


#Demultiplexing and cleaning


awk '{print $1}' list_barcodes | sort | uniq > list_lib

while read a
	do
	grep "$a" list_barcodes | awk '{print ">"$2"#^"$3}' | tr "#" "\n" > temp_barcodes.fasta
	cutadapt -e 0.17 --max-n 1 --minimum-length 30 -q 10 --no-indels -g file:temp_barcodes.fasta -o demux-{name}_R1_.fastq.gz -p demux-{name}_R2_.fastq.gz "$a"_L1_R1_001_*.fastq.gz "$a"_L1_R2_001_*.fastq.gz
	done < list_lib

for i in `ls demux-*.fastq.gz` ; do fastqc "$i" ; done


#Loci reconstruction and SNP calling using IpyRAD
ipyrad -p ./params_ipyrad90 -f -s 1234567 

------- ipyrad params file (v.0.6.5)--------------------------------------------
species_90			## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
./                              ## [1] [project_dir]: Project dir (made in curdir if not present)
                               ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                               ## [3] [barcodes_path]: Location of barcodes file
./*.fastq.gz                      ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)
                               ## [6] [reference_sequence]: Location of reference sequence file
pairddrad                      ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
                               ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
1000000                        ## [13] [maxdepth]: Max cluster depth within samples
0.90                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
0                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
0.05                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)
0.05                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)
2                              ## [21] [min_samples_locus]: Min # samples per locus for output
0.2                         ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)
20                           ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)
0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)
0, 0, 0, 0                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
p, s, v                        ## [27] [output_formats]: Output formats (see docs)
                               ## [28] [pop_assign_file]: Path to population assignment file


#Population structure
vcftools --vcf species_raw.vcf --min-alleles 2 --max-alleles 2 --remove-indels --max-missing 0.8 --recode --out species_snp_bi_miss08
random_select_snp_by_loci.sh species_snp_bi_miss08.recode.vcf 1snp_species_snp_bi_miss08.recode.vcf
vcf2structure_gn.sh 1snp_species_snp_bi_miss08.recode.vcf

##PCA using adegenet

data<-read.structure("1snp_species_snp_bi_miss08.recode.str",n.ind=xx,n.loc=xx,onerowperind=FALSE,col.lab=2,col.pop=1,col.others=NULL,row.marknames=NULL,NA.char=-9,ask=F) 


sansna<-scaleGen(data,scale=F, NA.method="mean")
pca1 <- dudi.pca(sansna,cent=F,scale=F,scannf=FALSE,nf=4)
barplot(pca1$eig[1:20],main="PCA eigenvalues", col=heat.colors(20))
s.class(pca1$li,xax=1,yax=2,pch=19)
(pca1$eig/sum(pca1$eig))*100
add.scatter.eig(pca1$eig[1:10], 3,1,2,posi = "topleft")
s.label(pca1$li)
add.scatter.eig(pca1$eig[1:20], 3,1,2,,posi = "bottomright")
s.class(pca1$li)
s.label(pca1$li,clabel = 0,)
s.label(pca1$li,boxes=FALSE,clabel=1)

##STRUCTURE
for i in `seq 1 3`
	do
#tested K
	for j in `seq 1 10`
		do
		mkdir st"$i"_"$j"_"$1".dir
		nb_sample=`grep -c "." $1 | awk '{print $1/2}'`
		nb_marker=`head -n 1 $1 | tr '\t' '\n' | grep "." -c | awk '{print $1-1}'`
		sed -e 's/KN/'"$j"'/g' -e 's/name_output/st'"$i"'_'"$j"'_results/g' -e 's/name_sample/'"$1"'/g' -e 's/nb_sample/'"$nb_sample"'/g' -e 's/nb_marker/'"$nb_marker"'/g' mainparams > ./st"$i"_"$j"_"$1".dir/mainparams
		cp extraparams ./st"$i"_"$j"_"$1".dir/
		cp structure ./st"$i"_"$j"_"$1".dir/
		r1=`shuf -i1-9 -n1`
		r2=`shuf -i1-9 -n1`
		r3=`shuf -i1-9 -n1`
		cp run.sh ./st"$i"_"$j"_"$1".dir/st"$i"_"$j".sh
		sed -e 's/structure/structure\ \-D\ '"$r1"''"$r2"''"$r3"'/g' ./st"$i"_"$j"_"$1".dir/st"$i"_"$j".sh > ./st"$i"_"$j"_"$1".dir/st"$i"_"$j"_2.sh
		cd ./st"$i"_"$j"_"$1".dir/
		chmod 755 ./st"$i"_"$j"_2.sh
		nohup ./st"$i"_"$j"_2.sh &
		cd ..
		done
	done


#Descriptive statistics
data<-read.structure("species_snp_bi_miss08.recode.str",n.ind=xx,n.loc=xx,onerowperind=FALSE,col.lab=1,col.pop=2,col.others=NULL,row.marknames=NULL,NA.char=-9,ask=F) 
data2 <- genind2hierfstat(data)
basicstat <- basic.stats(data2, diploid = TRUE, digits = 3)

Ho<-basicstat$Ho
Ho2<-data.frame(Ho)
results<-data.frame(colMeans(Ho2,na.rm = TRUE),sapply(Ho2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),table(data2$pop))

Fis<-basicstat$Fis
Fis2<-data.frame(Fis)
results<-data.frame(colMeans(Fis2,na.rm = TRUE),sapply(Fis2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),table(data2$pop))

Hs<-basicstat$Hs
Hs2<-data.frame(Hs)
results<-data.frame(colMeans(Hs2,na.rm = TRUE),sapply(Hs2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),table(data2$pop))

pairwise.WCfst(data2,diploid=TRUE)

#Migration rate

perl ./vcf2genepop.pl vcf=species_snp_bi_miss08.recode.vcf pops="$list_loc" > species_snp_bi_miss08.recode.genepop
divMigrate(infile = "species_snp_bi_miss08.recode.genepop", outfile = "species_out", boots = 100, stat = "all")

#Ne estimation
##position on genome
blastn -db GCA_963565745.1_ilPheArio1.1_genomic.fna -query species_raw.fasta -max_target_seqs 1 -outfmt 6 -out blast_species_on_genome
sort -k1,1 -k12,12gr -k11,11g -k3,3gr blast_species_on_genome | sort -u -k1,1 --merge > blast_species_on_genome_best
while read a
	do
	locus=`echo $a | awk '{print $1}'`
	pos=`echo $a | awk '{print $2}'`
	genot=`echo $a | cut -d " " -f 3- `
	scaff=`awk '{if ($1 == "'$locus'") print $3}' blast_species_on_genome_best`
	pos1=`awk '{if ($1 == "'$locus'") print $10}' blast_species_on_genome_best`
	pos2=`awk '{if ($1 == "'$locus'") print $11}' blast_species_on_genome_best`
	if [ "$pos1" -gt "$pos2" ]; then
	    newpos=$(($pos1+$pos))
	else
	    newpos=$(($pos1-$pos))
	fi
	echo $scaff $newpos $genot >> temp_vcf_good_pos
	done < temp_vcf_snp

## Ne estimations
library(dartR)
library(vcfR)
vcf <- vcfR::read.vcfR("species_good_pos_final.vcf")
pop <- read.table("species_pop_file.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
### Transform vcf to genind
genind <- vcfR::vcfR2genind(vcf)
### Define population map 
genind@pop <- pop$STRATA
### Transform the genind to genlight
genlight <- gi2gl(genind)
### Estimate Ne based on LD
LD.ne <- gl.LDNe(genlight)
write.table(LD.ne, "species_pop_Ne_LD.txt",quote=FALSE)

