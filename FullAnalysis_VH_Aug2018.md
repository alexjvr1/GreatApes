# Great Apes project with all samples

Vilhelmiina Haavisto will continue with this project, and extract data for all the sequenced Great Apes. 

##BLAST for bonobos (Pan paniscus)
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

Blast version: blastn, set max target sequences to 5000

Downloaded all hits in fasta format, and HIT table for information on start and end position of match with query sequence (from the top of the page!).

#Upload to server
scp *filename/file path* alexjvr@fgcz-c-047.uzh.ch:/srv/kenlab/alexjvr_p1795/GreatApes


