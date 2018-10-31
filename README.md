
### This code is associated with the paper from Di Rienzi et al., "Resilience of small intestinal beneficial bacteria to the toxicity of soybean oil fatty acids". eLife, 2018. http://dx.doi.org/10.7554/eLife.32581


#Lactobacillus_soybeanoil

#The files in this repository are the qiime, R, perl, and other bash commands utilized to process the data in the project, "Resilience of small intestinal lactobacilli to the toxicity of soybean oil fatty acids"



#### mappingfile.txt
###### provides the metadata for all mice microbiome samples

#### biomedit.pl
###### perl script used in the creation of Figure 4 – Figure supplement 2D

#### QiimeandRcommands.txt
###### Qiime and R commands for Figures 4-2A,B,D, 4-3, 4C,D, 4-2, 4B  

#### weighteddistancemetricpairs.pl
###### perl script used in the creation of Figure 4B

#### linoleicacid_threetimepointratio.R
###### R commands for Figures 4E, 4–4, 2, 2-1A,B

#### linoleicacidgrowthdata.txt
###### raw data for Figures 2 and 
###### note that missing metadata for the JW strains are filled in during the creation of the Figure using the information in JWstraininfo.txt

#### JWstraininfo.txt
###### provides the information on the L. reuteri strains isolated from a variety of animals
###### used in the creation of Figure 2

#### LactoGenomeAssembly.txt
###### pipeline using spades, Bandage, and nucmer to assemble the parental strains (LR0 and LJO)

#### GATK.txt
###### BWA-MEM, Picard, and GATK pipelines for identifying mutations in the strains and populations for evolved Lactobacilli

