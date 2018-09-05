# Great Apes project with all samples

Vilhelmiina Haavisto will continue with this project, and extract data for all the sequenced Great Apes. 

## BLAST for bonobos (Pan paniscus)
https://link.springer.com/article/10.1007%2Fs00251-017-0990-x#Sec2

All individuals can be found at https://www.ncbi.nlm.nih.gov/sra (filter for P.paniscus)

Used 
KY613168.1 
KY613128.1
KY613075.1
as query accessions (Papa-A, -B and -C)
https://www.ncbi.nlm.nih.gov/nuccore/KY613075.1
https://www.ncbi.nlm.nih.gov/nuccore/KY613128.1
https://www.ncbi.nlm.nih.gov/nuccore/KY613168.1

Blast version: megablast, set max target sequences to 5000. All other settings at defaults.

Downloaded all hits in fasta format, and HIT table for information on start and end position of match with query sequence (from the top of the page!). HIT table only needs to be downloaded once per blast (has all 3 queries in).

```
#Upload to server
scp *filename/file path* alexjvr@fgcz-c-047.uzh.ch:/srv/kenlab/alexjvr_p1795/GreatApes
```
In R (examples only as I did this for each hit table separately)
```
dze241291.ht <- read.table("dze241291_hittable.txt", header=F) #read in the table
dze241291.ht$INDV <- "Dzeeta" #add a column with the individual name

#After reading all the tables in, I merged them first by individual to make sure I got all of them
dzemerged <- rbind(dze241290.ht, dze241291.ht)

#Then merge them all
hittables.merged <- rbind(bonmerged, catmerged, chimerged, desmerged, dzemerged, hermerged, hormerged, kommerged, kosmerged, kummerged, lbmerged, natmerged, salmerged)

#And draw a histogram of the alignments lengths (located in the column aptly-named V4)
hist(hittables.merged$V4, xlab = "Alignment length", main = NULL, breaks = 50, col = "purple", ylim=c(0, 20000))
```
### Column name key
V1 = Query acc. version
V2 = Subject acc. version
V3 = % identity
V4 = Alignment length
V5 = Mismatches
V6 = Gap opens
V7 = Q start
V8 = Q end
V9 = S start
V10 = S end
V11 = E value
V12 = Bit score

NCBI insert size: 205 


# SRA-toolkit

We re-evaluated our methods. It might be better to use an optimised algorithm like bwa-mem to map all the SRA reads to the reference sequences and then call variants on these. Based on the BLAST results we have a distribution of the start site frequency of all the SRA reads as mapped to each one of the reference sequences. We expect that the regions of higher depth would be conserved between genes, and so sequences from all three loci are found with BLAST. 

So we've decided to use BWA mem to map the reads to a linked sequence of all three MHC genes + a reference gene (in this case [Interleuken 6](https://www.ncbi.nlm.nih.gov/nuccore/NC_000007.14?report=genbank&from=22725889&to=22732002)). 

First we'll download the data directly as paired end fastq data from NCBI using SRA-toolkit. 

This was downloaded directly in the working directory on fgcz47 server: 

```
alexjvr@fgcz-c-047:/srv/kenlab/alexjvr_p1795/GreatApes

fastq-dump --split-files SRR747648
```

## Mapping to our file containing MHC-A, -B, -C consensus sequences as well as IL6 sequence

```
/usr/local/ngseq/packages/Aligner/BWA/0.7.15/bin/bwa mem -t 10 /srv/kenlab/alexjvr_p1795/GreatApes/Mapping/MHCrefseq.fasta /srv/kenlab/alexjvr_p1795/GreatApes/SRR748147_1.fastq /srv/kenlab/alexjvr_p1795/GreatApes/SRR748147_2.fastq > aln-pe.SRR748147.sam

/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools import /srv/kenlab/alexjvr_p1795/GreatApes/Mapping/MHCrefseq.fasta aln-pe.SRR748147.sam aln-pe.SRR748147.bam

/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools sort aln-pe.SRR748147.bam -o SRR748147.bam.sorted

/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools flagstat SRR748147.bam.sorted
```

## 05.09.2018

```
alexjvr@fgcz-c-047:/srv/kenlab/alexjvr_p1795/GreatApes$

#Realigning files with new NoGaps MHCrefseq file
/usr/local/ngseq/packages/Aligner/BWA/0.7.15/bin/bwa mem -t 10 /srv/kenlab/alexjvr_p1795/GreatApes/Mapping/MHCrefseqNoGaps.fasta /srv/kenlab/alexjvr_p1795/GreatApes/SRR748147_1.fastq /srv/kenlab/alexjvr_p1795/GreatApes/SRR748147_2.fastq > aln-pe.SRR748147.ng.sam

/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools import /srv/kenlab/alexjvr_p1795/GreatApes/MHCrefseqNoGaps.fasta.fai aln-pe.SRR748147.ng.sam aln-pe.SRR748147.ng.bam

#Removing unmapped reads from .bam file (cuts down on size and comp. time)
/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools view -b -f 2 aln-pe.SRR748147.ng.bam > SRR748147.mo.ng.bam

#Sorting the mapped reads-only .bam file
samtools sort SRR748147.mo.ng.bam -o SRR748147.sorted.mo.ng.bam

#Adding and checking for read groups
/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools view -H SRR747648.sorted.bam | grep '^@RG'/
java -jar /usr/local/ngseq/packages/Tools/Picard/2.18.0/picard.jar AddOrReplaceReadGroups \
      I=SRR748147.sorted.mo.ng.bam \
      O=SRR748147.rg.sorted.mo.ng.bam \
      RGID=1 \
      RGLB=SRR748147 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=SRR748147


#Creating .dict and index files for GATK
java -jar /usr/local/ngseq/packages/Tools/Picard/2.18.0/picard.jar CreateSequenceDictionary R= MHCrefseqNoGaps.fasta O= MHCrefseqNoGaps.dict
/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools index SRR748147.rg.sorted.mo.ng.bam
samtools faidx MHCrefseqNoGaps.fasta

#SNP Calling
java -jar /usr/local/ngseq/packages/Variants/GATK/3.8.1.0/GenomeAnalysisTK.jar -R MHCrefseqNoGaps.fasta -T HaplotypeCaller -I SRR748147.rg.sorted.mo.ng.bam --emitRefConfidence GVCF -o SRR748147.g.vcf

#Data aggregation step - combining the gvcf files from the previous SNP calling step
java -jar /usr/local/ngseq/packages/Variants/GATK/3.8.1.0/GenomeAnalysisTK.jar\
   -T CombineGVCFs \
   -R MHCrefseqNoGaps.fasta \
   --variant SRR747648.g.vcf \
   --variant SRR748147.g.vcf \
   -o cohort.g.vcf

#Joint genotyping
java -jar /usr/local/ngseq/packages/Variants/GATK/3.8.1.0/GenomeAnalysisTK.jar \ -T GenotypeGVCFs \ -R MHCrefseqNoGaps.fasta \ --variant cohort.g.vcf \ -o cohort.output.vcf
```
## Some questions...

**1.	How does this final step call the variants? What are the default filters that resulted in the loss of so many loci?
Read Filters automatically applied to the data in GenotypeGVCFs processing (italics indicate those which are most likely influencing our dataset)**
* MalformedReadFilter: Filters out ‘malformed’ reads, including those with invalid alignment starts/ends, disagreeing headers, missing read groups, or empty/not supported CIGAR strings - seems to have most criteria so maybe this is what is doing the majority of the filtering?
* BadCigarFilter: filters out reads where the CIGAR strings do not match up, where there are clips, or where there are consecutive indels in the cigar (does our dataset have these?? I couldn’t find them, but maybe I’m not looking in the right place…)
* UnmappedReadFilter: filters out reads that are SAM flagged as unmapped – we already did this in an earlier step so I don’t think this filter is having much of an effect
* NotPrimaryAlignmentFilter: recognized SAM flag identifying secondary alignment (flag 256), so that records used for analysis are likely to be mapped in the right place and therefore relevant. Does not filter out supplementary alignments (flag 2048)
* FailsVendorQualityCheckFilter: filters out sam flag 512, which is reads that fail the vendor quality check
* DuplicateReadFilter: filters out duplicate reads

**2.	How are these loci distributed across our 4 reference genes? (both in terms of how many loci mapped to each gene and whether we have a representation of the whole gene in each case)**
27 loci mapped to IL6, 94 to MHC-A, 1 to MHC-B, 32 to MHC-C (-1 from grep output as they also appear in the header of the file, add up to 154 which was our number of sites)
``` grep -o 'IL6' cohort.output.vcf | wc -l ```
 
**3.	How was the mapping (bwa) affected by the removal of the insert (when we removed the -'s)? Can you compare how many loci we get when mapping to each of the references?**
```
samtools flagstat SRR748147.rg.sorted.mo.ng.bam
samtools flagstat SRR748147.sorted.rg.mappedonly.bam
```
SRR748147.rg.sorted.mo.ng.bam (mapped to no gaps reference) returns 141724 reads, second (maps to gapped reference) returns 69464 reads


**3 cont) And I think we need a literature/forum search to work out what the best option is for mapping loci when we expect there to be inserts.**

 

