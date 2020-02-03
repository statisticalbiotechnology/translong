## Logbook
### Questions
Recommended text editor for Python?
Which ANOVAs are appropriate?
Is it appropriate to perform the PCA for the full dataset at once, or only brain and liver data separately?
How should I handle the replicates for the same time point? Keep them separate or combine in some way?

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

### 3/2 2019
* Struggled with collecting gene targets for all TFs in one file, discovered some issues.
* Created a separate file for perfroming a test PCA on a single TF and attemped to gather its results in a dataframe where results from all TFs could be collected

* Looked into alternative methods for finding gene symbols for each gene ensembl ID and inserting it into dataframe with expression values
* Finished writing code for finding gene symbols for each gene ensembl ID and inserting it into dataframe with expression values
* Attempted to solve the issue that the gene list for each TF is read as a string and not a listS

