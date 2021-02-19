# MeChIP2

# Reference Genome


Notes below are stil scratch notes


```download_genome``` rule downloads generates the June 2020 GRCm39/mm39 (mouse) reference genome from the UCSC Genome Browser. The full code has been incuded below 




```
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz
tar zvfx chromFA.tar.gz 
cat *.fa > mm10.fa
rm chr*.fa
```
