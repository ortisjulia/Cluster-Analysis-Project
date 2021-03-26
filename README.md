# Cluster-Analysis-Project
Project belonging to the BINP28 course of the Bioinformatics program at Lund University. The question to be answered was the following: How many clusters are there in data and what is the probability of an individual belonging to each one? Use Structure-like approach.

**Introduction**
-----------

Cluster analysis is the task of grouping a set of objects in such a way that the ones that belong to the same group are more similar to each other than to those in other groups. It is a common technique for statistical data analysis, which is used in many fields, including bioinformatics (1).

Cluster analysis can be achieved by various algorithms that differ significantly in their understanding of what constitutes a cluster and how to identify them. The appropriate clustering algorithm and parameter settings depend on the individual data set and intended use of the results (1).

In this case, I am going to use cluster analysis for genetic clustering, which consists on using the similarity of genetic data to infer population structures (1).

There are six steps in cluster analysis (2):

1. Noise/outlier filtering

  Outliers are observations that are distant from the rest of the data due to errors or rare cases. They can negatively influence the analysis so it is better to filter them out. In my case, I should filter out Naxos2, which the instructions of the exercise already tell me that it is an outgroup to the other Taxa present in the file.

2. Choose distance measure

  The distance measure defines the dissimilarity between the objects that are being clustered.

3. Choose clustering criterion

  The clustering criterion helps us decide how the objects should be clustered. You can for example choose to minimize a cost function, like in optimization-based clustering, or to use a set of rules, such as in hierarchical clustering.

4. Choose clustering algorithm

  An algorithm to perform the clustering needs to be chosen based on the chosen distance measure and clustering criterion. Some popular cluster algorithms are K-means, hierarchical clustering, self-organization maps (SOM) and fuzzy clustering.

  When a K-means approach is used, K points are placed into the space represented by the object that are being clustered. These points represent initial group centroids, which are the cluster centers. Then, each object is assigned to the group that has the closest centroid. When all objects have been assigned, the positions of the K centroids are recalculated. Then the previous steps are repeated until the centroids no longer move. This produces a separation of the objects into groups from which the metric to be minimized can be calculated. The advantages of this algorithm is that it's fast and easy to interpret. However, it is easily affected by outliers, it has to be told how many groups (K) to find, no measure is provided of how well a data point fits in a cluster, and there's no guarantee to find the global optimum.

  There are different types of hierarchical clustering:
  - Agglomerative hierarchical clustering: starts with n clusters (one for each object) and at every step, two clusters are merged into one.
  - Divisive hierarchical clustering: starts with one cluster containing all objects and at every step a cluster is split into two.

  Moreover, there are also different types of linkages between the clusters:
  - Single linkage: the distance between two clusters is the minimum distance between the members of two clusters. It tends to generate long chains.
  - Complete linkage: the distance between two clusters is the maximum distance between the members of two clusters. It generates compact clusters.
  - Average linkage: the distance between two clusters is the mean distance between all pairs of objects of the two clusters. It is robust to outliers.
  - Average group linkage: the distance between two clusters is the distance between the means of the two clusters.

5. Validation of the results

  Check the resulting clusters for unexpected groups. For example, if a membrane protein is clustered with cytoplasmic proteins there might be something wrong.

6. Interpretation of the results

  It is important to know that the results depend on the choices made. Choosing different distance measures, clustering criterions or clustering algorithms will lead to different results.

In this project it is also important to know what genetic admixture is. According to Rius and Darling, intraspecific genetic admixture occurs when multiple divergent genetic linages come into gene flow contact and interbreed (3). In more plain words, it is the presence of DNA from a distantly-related species in an individual as a result of interbreeding between species that had been genetically isolated and had developed unique gene pools (4).

Admixture mapping is a method of gene mapping based on the hypothesis that differences in disease rates between populations are partly because of frequency differences in disease-causing genetic variants. In admixed populations these variants occur more often on chromosome segments inherited from the ancestral population with the higher disease variant frequency. In other words, gene flow between reproductively isolated populations results in chromosomal admixture with contributions from each contributing ancestral population (5).

**Software**
------------

In order to find out which clusters are in the data I tried out several software, including `VarClust` (6) and `variant-spark` (7). However, neither of them seemed to be adequate and didn't result in good results. For this reason, I decided to use `ADMIXTURE`, which is a high performance tool for estimating individual population allele frequencies from single nucleotide polymorphism (SNP) data (8). In addition to that software, there is a similar one called `STRUCTURE`, developed at Stanford University, that uses multi-locus genotype data to investigate population structure (9). Nevertheless, `ADMIXTURE` is easier to use, so even though I installed `STRUCTURE` in my local machine using the installation protocol found in their [website](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/install.html), I decided against using it because it took a long time to process the file and most of the time the software stopped working. All the information about how `STRUCTURE` is used can be found in their [documentation](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html).

For `ADMIXTURE` there is a quite good manual in their website: [ADMIXTURE documentation](http://dalexander.github.io/admixture/download.html). According to the documentation, `ADMIXTURE` uses a block relaxation approach to alternately update allele frequency and ancestry fraction parameters. Each block update is handled by solving a large number of independent convex optimization problems, which are tackled using a fast sequential quadratic programming algorithm. Convergence algorithm is accelerated using a novel quasi-Newton acceleration method (8, 9).

After using `ADMIXTURE` to obtain the result files, you need to use `R` or `Python` to plot the data. For some of the plots, I used [RStudio](https://rstudio.com/), which is an IDE to run `R`. The free version can be downloaded for Windows machine on [RStudio download](https://rstudio.com/products/rstudio/download/). Before installing RStudio, however, you need to install base R for Windows by going to the following [website](https://cran.rstudio.com/bin/windows/) and following the instructions. After installing base R, you can go back to the RStudio website and download it. After downloading it, you just need to follow the installation wizard instructions for a proper installation of the software.

For other plots I used `Python`, as it is easier to add the error bars in the plot. To run `Python` I used `Spyder`, which is an IDE that allows me to run `Python` and create `Python` scripts. To install `Spyder` on a Windows machine, you can download the `.exe` file from the [Spyder website](https://www.spyder-ide.org/) and follow the installation wizard instructions.

**ADMIXTURE**

**Installation**

In order to install `ADMIXTURE` you should first install `conda`. You can follow the instructions for its installation in the following website: [Anaconda installation](https://docs.anaconda.com/anaconda/install/linux/).

```
#To download Anaconda (conda) on Linux:
  #1 Download the Anaconda installer for Linux in your browser
  #2 Open a terminal and run the following:
     sha256sum /path/filename
  #3 Enter the following to install Anaconda for Python 3.7:
     bash ~/Downloads/Anaconda3-2020.02-Linux-x86_64.sh #The bash command has to be included whether or not you are using Bash shell
  #4 The installer makes you agree to the license agreement. Press Yes to agree after reading and reviewing it
  #5 Click enter to accept the default install location, CTRL-C to cancel the installation or specify an alternate installation directory. If you accept the default install location, the installer displays "PREFIX=/home/<user>/anaconda<3>" and continues to install
  #6 The installer asks you if you want to initialize Anaconda3 by running conda init. Say Yes
  #7 The installer finishes and thanks you for installing Anaconda3
  #8 The installer provides a link to PyCharm for Anaconda. No need to install it
  #9 Close and open your terminal window for the installation to take effect
  #10 Verify the installation by entering conda list on your terminal. If Anaconda is installed and working, this will display a list of installed packages and their versions.
```

Install `ADMIXTURE` using `conda` and the following command:

```
#To install ADMIXTURE
conda install admixture
```

Proceed by pressing `y` when asked.

**Usage**

`ADMIXTURE` input is binary `PLINK (.bed)`, ordinary `PLINK (.ped)` or `EIGENSTRAT (.geno)` formatted files and its output is simple space-delimited files containing the parameter estimates. In my case, I was given a `.vcf` file. For this reason, I need to get one of those formatted files. First, though, I need to take out the outgroup as it is not relevant in my analysis (it is an outlier) and can cause the algorithm to not work as good as it should. In order to do so I use the following command:

```
#To delete Naxos2 (the outgroup) from the VCF file
#First create a file with the samples you want (all minus Naxos2)
query -l ProjTaxa.vcf | head -n -1 > samples.txt

#Then use bcftools to select those samples from the ProjTaxa.vcf and put them in a new file called filtered.vcf
bcftools view -S samples.txt ProjTaxa.vcf > filtered.vcf

#The output will be a vcf file without the outgroup (just 15 samples instead of 16)

#To check if Naxos2 is present in the new file you can use:
bcftools query -l filtered.vcf

#The output should show all the other taxa minus Naxos2
```

After taking out the output, I can proceed to convert my vcf file to a `.bed` file using `PLINK`. This software can be installed using conda by following the command:

```
#Install PLINK
conda install plink
#You should press y when asked if you want to proceed
```

To get the `.ped` and `.map` files you can use the following command:

```
#To get .ped and .map from .vcf
vcftools --vcf filtered.vcf --plink --out filtered_taxa

#The software converts the unrecognized values used for CHROM: chrZ - Replacing with 0. Those values were giving me problems in the .bed file

#The output are a filtered_taxa.ped and filtered_taxa.map files.
```

You can also get the `.bed` files using the following command:

```
#To get .bed files
plink --file filtered_taxa --make-bed --out filtered_taxa
```

To use `ADMIXTURE` we need an input file and an idea of K, your belief of the number of ancestral populations. You should also have the associated support files alongside your main input file, in the same directory. For example, if your primary input file is a `.bed` file, you should have the associated `.bim` and `.fam` files in the same directory. Same goes for the `.ped` and `.map` files.

In order to choose the best K-value you can use `ADMIXTURE`'s cross-validation procedure, which means that a good value of K will exhibit a low cross-validation error compared to other K values. Cross-validation is enabled by just adding the `--cv` flag. In the default setting the cross-validation procedure will perform a 5-fold cross-validation that will be reported in the output. We can run the following command:

```
#To choose the best K for our analysis
for K in 1 2 3 4 5; do admixture --cv filtered_taxa.bed $K | tee log${K}.out; done
#The output will be five files with a .out extension with the cross-validation errors. To get them we can use grep:
grep -h CV log*.out
```

When I ran the analysis for all the 16 samples found in the vcf file I got no problem. However, when I ran it for the new vcf file with just 15 samples I got the following error: detected that all genotypes are missing from a SNP locus. The solution I found to make it work is to remove all loci where more than 99.9% of genotypes are missing. I used the following command:

```
#To remove all loci where more than 99.9% of genotypes are missing
plink --bfile filtered_taxa --geno 0.999 --make-bed --out filtered_taxa_geno

#The output is a new .bed file with only the non-missing loci
```

After removing the genotypes that where missing I could also run the cross-validation analysis on `ADMIXTURE` for the file with just 15 samples.

After finding out which is the most optimal K-value for the analysis, I can run `ADMIXTURE` like this:

```
#To run ADMIXTURE
admixture ProjTaxa.bed K #K is the number obtained for the best K value analysis
```

You can also get the standard errors by using `-B` flag:

```
#To get the standard errors
admixture -B ProjTaxa.ped K #K is the number obtained for the best K value analysis

#The output will be a .Q file, a .P file and a .Q_se, which is in the same unadorned file format as the point estimates
```

However, you need to be aware that if you calculate the standard errors the running time will be much longer. In the `ADMIXTURE` manual, it is explained the `-B` flag used to calculate the standard errors tells `ADMIXTURE` to perform point estimation and then also use a bootstrapping procedure to calculate the standard errors. This is why it takes much longer.

After obtaining the result files using `ADMIXTURE`, we need to use `R` or `Python` in order to obtain the plots. In the result section you can see all the plots for the data and how to obtain them.

**Results**

_Results when the outgroup is not removed_

For the file with all 16 samples (ProjTaxa.bed), the cross-validation results are the following:

```
#Cross-validation results for 16 samples
grep -h CV log*.out #grep -h finds the pattern CV in the files log*.out. The * refers to any number or character found after log and before the dot

#The output are the lines in which CV is found
CV error (K=1): 0.98901
CV error (K=2): 1.01794         
CV error (K=3): 1.12115        
CV error (K=4): 1.08782       
CV error (K=5): 1.32805
```

According to the `ADMIXTURE` manual, the most sensible decision is to choose the K value which gives the least cross-validation error. In this case, the one that gives the least amount of cross-validation error is K=1. So, for the data with all 16 samples the following code should be the one to run:

```
#Run ADMIXTURE
admixture ProjTaxa.bed 1
```

The output is a file for each parameter set: Q (the ancestry fractions) and P (the allele frequencies of the inferred ancestral populations). Note that the output files have a 1 in them, which indicates the number of populations (K) that was assumed for the analysis.

The Q estimates are output as a simple matrix, so we can plot it using the `read.table` and `plot` commands in R. You can also get the stacked bar-charts using the `barplot` command.

```
#To plot the results in R
tbl=read.table("ProjTaxa.1.Q")
par(mar=c(6,4,4,4))
barplot(t(as.matrix(tbl)), names.arg=c("8N05240", "8N05890", "8N06612", "8N73248", "8N73604", "K006", "K010", "K011", "K015", "K019", "Lesina_280", "Lesina_281", "Lesina_282", "Lesina_285", "Lesina_286", "Naxos2"), las=2, main="Structure-like cluster analysis", col="chartreuse3", xlab='Taxa', ylab="Ancestry", border=NA)

#par(mar()) allows the user to make the borders of the plot bigger so that the description of each column fits
#t(as.matrix) is used to transform the table into a matrix and transpose it
#names.arg is used to name the barplot's columns
#las is used to display the column names vertically instead of horizontally so that they all fit.
#main for the title of the plot
#col is used to choose the color of the plot
#xlab used to name the x-axis
#ylab used to name the y axis
#border is used to add a border line to the columns in the plot. NA means that none is added.
```

The resulting figure will be a barplot with 16 columns showing the ancestry of each individual studied. Nevertheless, in this case the plot is not relevant because they are all assigned to the same group, therefore, the plot is the same color for all the individuals.

![ProjTaxaPlot](https://user-images.githubusercontent.com/70640998/112621196-c1228380-8e29-11eb-810c-68900b06ed28.png)

**Figure 1.** Structure-like cluster analysis of ProjTaxa.vcf file including the outgroup (Naxos2). The best K value obtained was 1, which groups all the taxa found in the file in 1 cluster. This is probably due to the presence of outliers (the outgroup).

As I previously said, the reason for this results (figure 1) might be that I used an outgroup in the analysis. This group had obviously become a problem in the analysis as the cross-validation analysis showed that the best was to use a K of 1, which groups all the individuals together. This is clearly not what I want. Nevertheless, I performed this analysis with the outgroup to visually show that in this case it is better to remove it.

_Results when the outgroup is removed_

For the file with all 15 samples (filtered_taxa_geno.bed), the cross-validation results are the following:

```
#Cross-validation results for 16 samples
grep -h CV log*.out #grep -h finds the pattern CV in the files log*.out. The * refers to any number or character found after log and before the dot

#The output are the lines in which CV is found
CV error (K=1): 0.60277
CV error (K=2): 0.58891
CV error (K=3): 0.74078
CV error (K=4): 1.06903
CV error (K=5): 1.08233
```

As a good value of K will exhibit a low cross-validation error compared to the other K values, `ADMIXTURE` will be ran with a K of 2 because in the cross-validation results it was the best lowest, therefore the best one.

```
#Run ADMIXTURE
admixture filtered_taxa_geno.bed 2
```

The Q estimates are output as a simple matrix, so we can plot it using the `read.table` and `plot` commands in R. You can also get the stacked bar-charts using the `barplot` command.

```
#To plot the results in R
tbl=read.table("filtered_taxa_geno.2.Q")
par(mar=c(6,4,4,4))
barplot(t(as.matrix(tbl)), names.arg=c("8N05240", "8N05890", "8N06612", "8N73248", "8N73604", "K006", "K010", "K011", "K015", "K019", "Lesina_280", "Lesina_281", "Lesina_282", "Lesina_285", "Lesina_286"), las=2, main="Cluster analysis", col=c("chartreuse3", "darkblue"), xlab='Taxa', ylab="Ancestry", border=NA)

#par(mar()) allows the user to make the borders of the plot bigger so that the description of each column fits
#t(as.matrix) is used to transform the table into a matrix and transpose it
#names.arg is used to name the barplot's columns
#las is used to display the column names vertically instead of horizontally so that they all fit.
#main for the title of the plot
#col is used to choose the color of the plot
#xlab used to name the x-axis
#ylab used to name the y axis
#border is used to add a border line to the columns in the plot. NA means that none is added.
```

![filtered_taxa_geno_plot](https://user-images.githubusercontent.com/70640998/112621352-ed3e0480-8e29-11eb-8146-181243f19778.png)

**Figure 2.** Structure-like cluster analysis of ProjTaxa.vcf file excluding the outgroup (Naxos2). The best K value obtained was 2, which groups all the taxa found in the file in 2 clusters, seen in blue and green. We can observe that 8N taxa are all grouped together in cluster blue, while Lesina taxa are all grouped together in cluster green. Nevertheless, in K0 taxa both clusters can be found, seen by the presence of both blue and green in those taxa.

As you can observe in figure 2, this plot yields more interesting results. We can clearly see that a K of 2 separates the taxa in different groups. The 8N taxa, which are all grouped together, can be seen in blue. The Lesina taxa, also grouped together, are seen in green. However, we see something quite interesting, the K0 taxa have a mix of both blue and green, which suggests that they might be a descendant taxa groups from the other two taxa groups studied here. They also seem more closely related to the 8N group taxa than the Lesina taxa group, as they share higher ancestry percentage.

_Results when the outgroup is removed and with bootstrapping_

I rerun the `ADMIXTURE` analysis with the `-B` flag in order to obtain the standard errors of the analysis. As I had already cross-validated to find the best value for K, I did not repeat the cross-validation analysis. The K I used was 2, the same as previously. When `-B` flag is used without any other specification, the default bootstrapping analysis is performed. This is that 200 replicates are done. If you want a different amount of replicates you just need to specify it after the flag, such as: `-B500`, which would perform 500 replicates. You need to be aware that the more replicates you do, the more time it will take to run.

```
#Run ADMIXTURE and obtain standard error using bootstrapping
admixture -B filtered_taxa_geno.bed 2

#The analysis will take about an hour and it will output a .Q file, a .P file and a .Q_se file. This last one will contain the standard error.
```

In order to plot the results, I use the following `Python` code:

```
#Plot results with standard error bars:
#Import the packages I will need
import numpy as np
import matplotlib.pyplot as plt

#The amount of columns I will have (depends on the amount of taxa)
N = 15

#C1 is the cluster 1 data
C1 = (0, 0, 0, 0, 0, 0.28, 0.30, 0.17, 0.22, 0.17, 1, 1, 1, 1, 1)
#C1 standard error
C1Std = (0, 0, 0, 0, 0, 0.050, 0.063, 0.050, 0.050, 0.055, 0, 0, 0, 0, 0)

#C2 is the cluster 2 data
C2 = (1, 1, 1, 1, 1, 0.72, 0.70, 0.83, 0.78, 0.83, 0, 0, 0, 0, 0)
#C2 standard error
C2Std = (0, 0, 0, 0, 0, 0.050, 0.063, 0.050, 0.050, 0.055, 0, 0, 0, 0, 0)

# the x locations for the taxa
ind = np.arange(N)

# the width of the bars
width = 0.35

#To create the plot bars for C1
p1 = plt.bar(ind, C1, width, yerr=C1Std, color='green')

#To create the plot bars for C2
p2 = plt.bar(ind, C2, width, bottom=C1, color='blue')

#To create the actual plot
plt.ylabel('Ancestry') #Y-axis label
plt.xlabel('Taxa') #X-axis label
plt.title('Structure-like cluster analysis') #Title of the plot
plt.xticks(ind, ('8N05240', '8N05890', '8N06612', '8N73248', '8N73604', 'K006', 'K010', 'K011', 'K015', 'K019', 'Lesina_280', 'Lesina_281', 'Lesina_282', 'Lesina_285', 'Lesina_286'), rotation='vertical') #Name of the bars in the X-axis
plt.yticks(np.arange(0, 1.5, 0.2)) #Numbers in the Y-axis
plt.legend((p1[0], p2[0]), ('Cluster 1', 'Cluster 2')) #Legend of the plot
plt.show() #Plot
```

![filtered_taxa_geno_plot_errorbars](https://user-images.githubusercontent.com/70640998/112621424-01820180-8e2a-11eb-93d3-445626eea9ef.png)

**Figure 3.** Structure-like cluster analysis of ProjTaxa.vcf file excluding the outgroup (Naxos2). The best K value obtained was 2, which groups all the taxa found in the file in 2 clusters, seen in blue and green. We can see that this plot is exactly the same as the one in figure 2 but with the error bars.

As you can see in figure 3, the plot I obtained is the same as in the analysis without bootstrapping. Nevertheless, the bootstrapping gave me an extra file with the standard errors. In order to properly show the results, I use `Python` to add the error bars in the plot. You can clearly see that 8N taxa group and Lesina taxa group have a standard error of 0, which suggests a perfect placement in their respective clusters. Regarding K0 taxa, we can see that there is some standard error, which suggests that there are some differences in the cluster placement of K0 in the bootstrapping replicates. Nevertheless, the standard error is not too big, as seen in table 1.

**Table 1.** Standard error of each taxa group estimated by using bootstrapping. C1 and C2 stand for cluster 1 and cluster 2, respectively.

![filtered_taxa_geno_error_8N](https://user-images.githubusercontent.com/70640998/112621475-15c5fe80-8e2a-11eb-862d-88dd3f571c3d.jpg)

![filtered_taxa_geno_error_K0](https://user-images.githubusercontent.com/70640998/112621485-18c0ef00-8e2a-11eb-8a47-b10503862826.jpg)

![filtered_taxa_geno_error_Lesina](https://user-images.githubusercontent.com/70640998/112621497-1ced0c80-8e2a-11eb-8d7e-8abda3d373ec.jpg)

This table was obtained by using `R` and running the following script:

```
#Libraries needed
library(gridExtra)
library(grid)

#Get the data from the standard error file produced in admixture
se <- read.table("filtered_taxa_geno.2.Q_se")

#Transform it into a transposed matrix
bse <- t(as.matrix(se))

#Add column names to the matrix
colnames(bse) <- c("8N05240", "8N05890", "8N06612", "8N73248", "8N73604", "K006", "K010", "K011", "K015", "K019", "Lesina_280", "Lesina_281", "Lesina_282", "Lesina_285", "Lesina_286")

#Divide the data into 3 parts so it creates 3 smaller tables (it is done like this because one table did not fit)
bse_8N <- bse[1:2, 1:5] #8N taxa table
bse_K0 <- bse[1:2, 6:10] #K0 taxa table
bse_Lesina <- bse[1:2, 11:15] #Lesina taxa table

#Finally you can create the 3 grid tables
grid.table(bse_8N)
grid.table(bse_K0)
grid.table(bse_Lesina)
```

**Conclusions**

As we have seen in this structure-like cluster analysis of the data found in the `ProjTaxa.vcf` file, it is a good idea to remove outgroups from the data before starting. Moreover, we have seen that `ADMIXTURE` is a good software to get structure-like cluster plots in combination with `R`.

My results show that the taxa K006, K010, K011, K015 and K019 taxa (from now called K0 taxa group) have both clusters present, the one in which the Lesina taxa group (Lesina_280, Lesina_281, Lesina_282, Lesina_285 and Lesina_286 taxa) belongs to and the one in which the 8N taxa group (8N05240, 8N05890, 8N06612, 8N73248 and 8N73604 taxa) belongs to. This allows me to conclude that Lesina and 8N taxa groups were reproductively isolated populations that had gene flow resulting in the third taxa group that has a chromosomal admixture with contributions from each contributing ancestral population. Therefore, the data seems to point to the fact that Lesina and 8N taxa are ancestor taxa to K0 taxa groups. Moreover, the data seems to show that K0 taxa groups are more closely related to 8N taxa groups, as a higher percentage of ancestry is found between them than with the Lesina taxa groups. Nevertheless, phylogenetic analysis needs to be done in order to properly conclude that.

**References**
------------

1. Wikipedia contributors. (2021, February 8). Cluster analysis. Wikipedia. https://en.wikipedia.org/wiki/Cluster_analysis

2. Center for integrative bioinformatics VU. (2008). Clustering algorithms. https://www.ibi.vu.nl/teaching/masters/bi_tools/2008/tools_lec4_2008_handout.pdf

3. Rius, M., & Darling, J. A. (2014). How important is intraspecific genetic admixture to the success of colonising populations? Trends in Ecology & Evolution, 29(4), 233–242. https://doi.org/10.1016/j.tree.2014.02.003

4. Wikipedia contributors. (2020, December 18). Genetic admixture. Wikipedia. https://en.wikipedia.org/wiki/Genetic_admixture

5. Winkler, C. A., Nelson, G. W., & Smith, M. W. (2010). Admixture Mapping Comes of Age. Annual Review of Genomics and Human Genetics, 11(1), 65–89. https://doi.org/10.1146/annurev-genom-082509-141523

6. F. (2019). fasterius/VarClust. GitHub. https://github.com/fasterius/VarClust

7. A. (2021). aehrc/VariantSpark. GitHub. https://github.com/aehrc/VariantSpark

8. Alexander, D.H., Lange, K. Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. BMC Bioinformatics 12, 246 (2011). https://doi.org/10.1186/1471-2105-12-246

9. Structure Software for Population Genetics Inference. (2021). Structure Software. https://web.stanford.edu/group/pritchardlab/structure.html
