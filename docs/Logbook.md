## Logbook
### Questions

* Regarding significance: is q-value not an option with the p-value distribution this data has?

* Regarding significance: what alternative method do you suggest should be used? Bonferroni correction?

* Regarding significance: should significant results be reported, or just aquired p-values? Worry of ad hoc since the full experiment has been run (without log2 transformation) prior to setting a threshold.

* Why is the standardization before PCA done on each individual gene set and not on the full gene set once? We want to center the data to look at the variation specifically among the selected genes, and don't want it to be affected by other genes.

* If I am to plot the change in expression level of the single gene coding for a TF, should it be done with data that has been standardized similarly to the data used for PCA (but then for the full dataset once)? Important to note that 

* Ideas for interpretation of the fact that PCA on the full dataset leads to the first PC explaining >35% of the variation? General increase in expression over time? General difference between organs? Sampling error? Difference between organs is likely what is moslty captured, since it is a clear separation. May be cause for the p-values for C(dev_stage) being high.

* How does the sklean PCA function deal with replicates? How should I actually deal with replicates? No problem if trend is the same, but information is lost if they are different and we simply classify them as separate datapoints.

### Ideas/plans/thoughts
* ~~Plot PC1 vs PC2 to see if different organs cluster along PC1 axis~~

* ~~Perform PCA on organs individually~~

* ~~Save figures of activity of TFs found to be most significantly differentially expressed over time and across organs~~

* ~~Compare results for TF activity aqcuired with PCA and the measured level of mRNA for the gene~~

* ~~Try to reduce the size of the TF target gene sets~~

* ~~Use MACS2 and STRING binidng score for limiting TF target gene sets~~

* ~~Remove TFs with more than eg 500 genes~~

* ~~Plot number of genes in TF gene set against variance expained by the first principal component~~

* ~~Do PCA on full gene set~~


* In case much time is available: make a figure where each gene's activity for a TF is shown as colour coded squares for each sample, and compare to the overall actvity change for TF for each sample

* Look for significant TFs defined by their 2nd principal component

* Look for singificant TFs defined by their 1st and 2nd principal components

* Use porch for q-value calc (if the method is applicable)

* Look further into the TFs found to be most significantly differentially expressed over time and across organs, their relevance, function etc.

* Look into which genes explain most of the variation for interesting TFs, how much of the first principal component they make up

* Look into issues with the method, how principal components depend on linear relationships

* Consider the effect of TF regulation being both activating and repressing for when doing the PCA. Is this part of the power of PCA, that it is taken into consideration?

* Try performing the analysis with weighted PCA

* Try to perform the analysis with porch

* Make sunburst plots?

* Filter TF gene sets by only picking the n genes with the highest MACS2 score for any sample, with no hard MACS2 cutoff. n=500 may give too many genes, making it difficult to find which are important?

* Skip doing an ANOVA. Instead look at activity variation over time for predetermined, well studied TFs that are known to have varying activity and try to validate the method by seeing if the patterns are similar.

* TF gene expression could potentially be used as an indication of correct results for validation.

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

### 17/2 2020
* Ran PCA with genes potentially regulated by TF was defined as genes +-1kbp from TF binding site

### 18/2 2020
* Worked on transferreing desirable features of Overleaf thesis template to the KTH template
* Transferred what had previously been written in Overleaf template to KTH template
* Worked on code for creating stricter TF target gene sets (by setting a MACS2 score threshold, that is based on q-value from ChIP-seq peak calling)

### 19/2 2020
* Started writing on background for thesis (about principal component analysis)
* Looked further into creating smaller TF target gene sets (tested higher MACS2 score thresholds)

### 20/2 2020
* Continued writing on background for thesis (PCA)
* Tested higher MACS2 score thresholds (700, 800) for TF target gene set determination

### 21/2 2020
* Wrote code for plotting gene set size vs ratio of variance explained to see correlation and ran said code for datasets from different MACS2 score thresholds
* Started looking into plotting a heatmap comparing a change over  time of a TF's PC and each individual gene in its gene set

### 24/2 2020
* Continued writing on thesis (background - ANOVA)

### 25/2 2020
* Continued writing on thesis (background - ANOVA)
* Tried running PCA on full dataset (all genes) to see how much of the variance the first PC would explain
* Tested TF gene set generation with very high MACS2 score threshold.
* Improved code for fetching TF gene sets to be able to handle temporary loss of internet connection

### 26/2 2020
* Continued writing on thesis (background - ANOVA, materials and methods)
* Compared gene sets from different MACS2 score thesholds and ran a full experiment with a MACS2 threshold of 700 and a max TF gene set size of 500 genes (TFs with more genes were removed).

### 27/2 2020
* Continued writing on thesis (materials and methods)
* Made improvements to figures
* Optimized code a bit
* Looked into using porch for performing PCA

### 28/2 2020
* Wasted time
* Made efforts to run PCA with porch, but encountered many problems (structure of dataframes, handling of zero values in expression_df, getting the ratio of variance explained)
* Attempted to improve/change produced figures of TF 'activity' over time
* Reflected on issues with the project

### 2/3 2020
* Discussed my questions with Gus (handling of replicates is a question worth considering, but will not be handled differently in this project)
* Adapted code to be able to give results for more than the first PC (must still fix the cases when a TF has fewer genes than the requested number of PCs)
* Wrote code for performing PCA on organs separately

### 3/3 2020
* Fixed so code can give results for more than the first PC even in the cases when a TF has fewer genes than the requested number of PCs (Do PCA for fewer PCs and fill emptry results with NaN)
* Wrote code to plot results from PC1 against results from PC2, colouring points based on developmental stage and with shapes dependent on organ.

### 4/3 2020
* Reflected on what I've done in the project so far
* Thought about how to continue with the project
* Improved the plot of PC1 vs PC2 by separating datapoints from different organs by colour, which is gradiated to show at which timepoint they were sampled

### 5/3 2020
* Forgot to write in logbook

### 6/3 2020
* Reflected on how to continue the project, what can be reported as results, what should be done for validation...
* Discussed how to continue the project with Patrik and Gus.

### 9/3 2020
* Created figure with 4 subplots showing PC1 vs PC2 with two different coloring schemes for filtered (so no TFs with more than 500 genes) or unfiltered TF gene sets, to be presented for the weekly meeting.
* Discussed about how to continue the project. Came up with new workflow for testing how well the method works: First test for obvious case TFs (high and varying expression), check if "activity" pattern resembles expression, test with non-obvious TFs (always present but require activation), check literature if results seem accurate/probable.
* Wrote code for finding TFs whose gene expression is high and with high variance to be selected for continued analysis.
* Fixed bugs for function for finding gene sets of the N genes with highest MACS2 score

