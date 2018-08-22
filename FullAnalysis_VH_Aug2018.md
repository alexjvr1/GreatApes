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
