## Logbook
### Questions
Recommended IDE (integrated development environment) for Python?
~~Should data preparation for log2 transformation be done by removing rows with 0 activity or should a constant (1) be added to each value?~~
Regarding significance: is q-value not an option with the p-value distribution this data has?
Regarding significance: what alternative method do you suggest should be used? Bonferroni correction?
Regarding significance: should significant results be reported, or just aquired p-values? Worry of ad hoc since the full experiment has been run (without log2 transformation) prior to setting a threshold.

How would using porch instead of what I've already done affect the results? How does it operate differently?

### Ideas/plans/thoughts
Perform PCA on organs individually
Try to perform the analysis with porch
Try performing the analysis with weighted PCA
~~Save figures of activity of TFs found to be most significantly differentially expressed over time and across organs~~
Look further into the TFs found to be most significantly differentially expressed over time and across organs, their relevance, function etc.
Look into which genes explain most of the variation for interesting TFs, how much of the first principal component they make up
Test including a second principal component.
Look into issues with the method, how principal components depend on linear relationships

Consider the effect of TF regulation being both activating and repressing for when doing the PCA. Is this part of the power of PCA, that it is taken into consideration?

~~Compare results for TF activity aqcuired with PCA and the measured level of mRNA for the gene~~

### 16/1 2020
* Read *A Quick Guide to Organizing Computational Biology Projects* (William Stafford Noble, 2009)
* Followed the guide *Version Control with Git* provided by Lukas
* Created local directory “translong” and connected it with GitHub directory “translong” 
* Added README file to “translong” directory and added directories for data, experiments and source code with readme files
* Created Google document “Project plan” and set up structure
* Started looking for sources for project plan background
* Started looking into Miriam’s project directory to understand the code
* Created Google document “Notebook (temporary)” for tracking activities

### 17/1 2020
* Started working on a WBS for the project plan
* Started working on a time plan for the project plan
* Started working on code for fetching transcription factor data

### 20/1 2020
* Continued working on WBS for the project plan
* Continued working on time plan for the project plan
* Started writing code for performing PCA, reading and copying Miriam’s code to understand
* Started reading about transcription factors

### 21/1 2020
* Set up logbook in markdown format
* Continued reading up on transcription factors and ehancer regulatory effects
* Continued writing on background for project plan
* Defined project objectives for project plan

### 22/1 2020
* Continued writing on background for project plan
* Looked into using LaTeX
* Read *Transcription factors: from enhancer binding to developmental control* (Spitz, F., Furlong, E. E. M., 2012)
* Partly read *Princiapl component analysus* (Abdi, H., Williams, L. J., 2010)

### 23/1 2020
* Read *Determinants of enhancer and promoter activities of regulatory elements* (Andersson, R., Sandelin, A., 2020)
* Continued working on background, WBS and gantt chart for project plan

### 24/1 2020
* Read *Chapter Twenty-Four - Cell Fate Determination by Transcription Factors* (Gurdon, J. B., 2016)
* Read *Transcription factors read epigenetics* (Hughes, T. R., Lambert, S. A., 2017)
* Continued working on the project plan
* Tested parts of Miriam's code to get a better understanding of the project

### 27/1 2020
* Made refinements to the project plan (updated background, WBS and timeplan)
* Looked further into the ANOVA performed by Miriam
* Started writing code for creating a csv file with all TFs and their associated genes listed

### 28/1 2020
* Updated timeplan in project plan
* Updated "Background" for project plan in accordance with Luka's comments
* Continued writing code for creating a csv file with all TFs and their associated genes listed

### 29/1 2020
* Finished code for creating a csv file with all TFs and their associated genes listed, though it required some manual editing of dowloaded csv with urls for fetching csvs for target genes for all TFs
* Started writing code for creating a dataframe with expression values that has genes as row index, and organ and time point as column index

### 30/1 2020
* Attempted to write on code for finding gene symbols for each gene ensembl ID and inserting it into dataframe with expression values

### 31/1 2020
* Made functional code for replacing ensembl gene IDs with gene symbols
* Performed PCA for a single TF

### 3/2 2020
* Struggled with collecting gene targets for all TFs in one file, discovered some issues.
* Created a separate file for perfroming a test PCA on a single TF and attemped to gather its results in a dataframe where results from all TFs could be collected
* Looked into alternative methods for finding gene symbols for each gene ensembl ID and inserting it into dataframe with expression values
* Finished writing code for finding gene symbols for each gene ensembl ID and inserting it into dataframe with expression values
* Attempted to solve the issue that the gene list for each TF is read as a string and not a list

### 4/2 2020
* Handed in the first project plan
* Created separate file for performing ANOVA on results from analysis of a single TF
* Attempted to perform ANOVA (must look further into how it should be performed)
* Looked into ratio of variance explained by the first principal component for the single TF

### 5/2 2020
* Discussed with Lukas about doing the PCA for the full dataset or separately for organs. Separately makes sense, but must do on full dataset to be able to compare. Do separately to look for genes with large impact
* Discussed with Lukas about potentially using their "porch" package to perform the PCA
* Finished code for performing ANOVA
* Created a file for getting q-values (will need to perform PCA on multiple TFs)
* Created a file for plotting results of PCA of one TF and ran the program
* Added so that the a dataframe with the ratio of variance explained is created when PCA is performed

### 6/2 2020
* Finished code for performing PCA on all TFs, calculating p-values and calculating q-values for variance between organs, developmental stage and their combination
* Ran an experiment for all TFs

### 7/2 2020
* Looked further into understanding the code for q-value calculation
* Looked into using LaTeX, with the KTH template (struggled a lot)

### 10/2 2020
* Created a LaTeX document on Overleaf from a KTH template, moved over background text from the project plan to thesis introduction, and added sources and acronyms to function with LaTeX.
* Log transformed (log2) the data and performed another experiment, showing more promising results

### 11/2 2020
* Optimized function for replacing ensembl IDs with gene symbols
* Made an alternative function for replacing ensembl IDs with gene symbols using MyGene
* Updated code for plotting results
* Considered if preparation of data for log2 transformation should be done by adding constant (1) or removing rows with zero values. Removing rows loses a lot of data. ASK!
* Worked on an outline for the thesis in Overleaf.

### 12/2 2020
* Talked to Lukas and Gustavo about log2 transformation. Lukas suggested weighted PCA. Gustavo suggested adding constant or sklearn transformation.
* Read into weighted PCA, what it is and what it is for
* Looked into porch to see how it can be applied
* Remodeled parts of code to make more of it into functions

### 13/2 2020
* Looked into extracting the percentage each gene contributes to the first pricipal component of each TF
* Plotted the expression level of TFs found to have significant differences over organ and developmental stage (Hnf1a, Rtf1) to compare to the expression of their PCs.
* Looked further into th relavance of Hnf1a. Which genes it has been reported to affect and how it compares to ChIP-atlas data.

### 14/2 2020
* Added so that created figures are saved in the /exp folder in the specific subfolder for the day
* Refined the histogram of the p-value distribution to include all three groups (organ, dev_stage, organ:dev stage)
* Wrote some thoughts into the Thesis as outline for discussion
* Wrote experimental code that shows histogram of the distribution of how much genes defining a TF contribute to the first principal component for the TF. Applied it to the two TFs found most interesting
* Improved annotation of code

