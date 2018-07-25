# Great Apes Sequence Prep  

Sequences will be extracted from NCBI using BLAST. 

The BioProject SRA data for the Great Apes project: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=189439

Great Apes Project website: http://biologiaevolutiva.org/greatape/data.html

### BLAST conditions: 

What was the query sequence? 

Are we using other sequences as well? 

Blast version?  

Reference? 


Download all hits in fasta format. 

Download HIT table for information on start and end position of match with query sequence. 

This is saved independently for each run. There are multiple runs per individual. 


## Sequence alignments and choosing loci

The first thing we need to figure out is where we're making cut-offs for different loci. 

Look at the distribution of start and end points for all the data. 


```
R

hist(sort(Dian.data$StartQuery), breaks=20)

```






