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

I'm going to import all the data into R, and then draw a histogram of all the sequence start points. 

Open R in the folder with all the data. On my computer this is here: 

/Users/alexjvr/2018.postdoc/GorillaMSc_201807
```
R version 3.5.0

Dian.data <- read.csv("DianSRX243537.csv", header=T) # import the csv file into R. I've added a row with headers at the top of this HIT file. If this isn't there, you can say header=F. 

hist(Dian.data$StartQuery, breaks=20)  # hist is the function to draw a histogram. Breaks is a way of binning the data. The more breaks, the smaller the bins. Call the column you want to plot by referring to the data and then the specifying the column ($). 

# to save any graphs or tables to pdf you can do the following
pdf("Figurename.pdf") ## give the pdf a name. Remember the extension

hist(Dian.data$StartQuery, breaks=20) ## all the code for the figure. This is just some simple code, but often these can have a lot more lines of code. 

dev.off() ## signals the end of the pdf file, and that it should be written to file. It will be saved in the folder where you're working, unless you change the path. 
```

Google any questions. Usually the top hits will be helpful, particularly on [R-bloggers](https://www.r-bloggers.com/basics-of-histograms/) and [Stackoverflow](https://stackoverflow.com/questions/16931895/exact-number-of-bins-in-histogram-in-r). 

Similarly, you can google anything to do with github and coding in github markdown. You can't break anything, so don't be afraid to try things! 







