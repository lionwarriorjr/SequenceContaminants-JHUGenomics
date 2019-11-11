##How to download raw fastq files from NCIB Sequence Reads Archive (SRA)
1. Download SRAtool kit from NCBI (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software). We need this tool
to unpack the data downloaded from SRA
2. Unzip the SRAtool kit
3. You can move the SRAtool kit directory to 'Application' directory if you want
4. Find the sequence you want to download from SRA (https://www.ncbi.nlm.nih.gov/sra). Note that
most of the raw fastq files are really big because they are usually sequences of the 
whole genome. Try to look for files that are not as big.
5. Once you find the sequence that you want to download, click on the Run ID. You should be directed
to the 'Metadata' tab. Click on the 'Data access' tab and download the data.
6. The file you just downloaded is in SRA format, so we need to use
SRAtool kit to convert it into fastq file.
7. cd into the bin directory under your SRAtool kit directory and use the following command in your terminal. If the
the file you downloaded is paired end (you can find out from metadata tab), use
this command in your terminal
```shell script
./fastq-dump --split-files -O outputPath SRAfilePath
```
For single end fastq, use this instead
```shell script
./fastq-dump -O outputPath SRAfilePath
```


