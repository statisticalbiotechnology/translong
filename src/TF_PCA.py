import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
import time
#import mygene		#for changing ensembl names to "normal" gene names

#Function for performing the PCA. Takes dataframe with expression values as input
def my_pca(df, n_pc=1, normalize=True):
	df = df.dropna(axis = 0, how = 'all') 		#Remove rows with only NA values 
	x = df.values.T 							#Set x as transpose of only the numerical values of the dataframe
	if normalize:	
		x2 = preprocessing.scale(x)				#Standardize the data (center to mean and scale to unit variance)
	else:
		x2 = x
	pca = PCA(n_components = n_pc)				#Set PCA parameters
	pca.fit(x2)									#Fit data to model
#	my_pca.pca = pca							#Used for Patrik's dictionary method for components
	x3 = pca.fit_transform(x2)					#Transform the data and set x3 as principal components 
	out_df = pd.DataFrame(x3.transpose(), index=list(range(1,n_pc+1)), columns=df.columns) #Create dataframe of traspose of output with index set as nr of PC+1 (?)
	out_df = out_df.transpose()
	return out_df

#Function for creating csv with all TFs and their target genes. For each TF, identifies urls for target gene data and fetches the info and puts it into a list that is then converted to a dataframe and given as output
def TF_targets(TFrange):   
	chip = pd.read_csv('../data/Transfactors/chip_atlas_analysis_list_CORRECTED.csv', quotechar='|').loc[1255:1953,['Antigen','Target Genes (TSS }1k)','Target Genes (TSS }5k)', 'Target Genes (TSS }10k)']]
	chip = chip.set_index('Antigen')
	# ~ print(chip.loc[:,TFrange])	
	TF_gene_list = []
	dex = 0
	for url in chip.loc[:,TFrange]:							#OBS! This takes a long time. The original chip_analysis_list.csv that was downloaded did not work, as an extra '"' was added to the start and end of each row with '"' in it and to each '"' to double to '""'
		try:
			TF_gene_set = pd.read_csv(url, sep='\t')		#For a specific TF, read csv from url as a dataframe
			genes = TF_gene_set['Target_genes'].tolist()	#Take the contents of column 'Target_genes' and puts it into a list
			TF_gene_list.append(genes)						#Append the list for a specific TF to list with all TFs
			print('Genes for '+chip.index[dex]+' found')
		except HTTPError:									#If the url can't be reached, insert 'Not found' in the list and continue (to get correct index)
			genes = ['Not found']
			print('Genes for '+chip.index[dex]+' NOT found')
		dex=dex+1	
	TF_gene_sets = pd.DataFrame({'Genes':TF_gene_list}, index=chip.index)	#Create a dataframe from the list of TFs and their target genes
	return TF_gene_sets

###################### Import data ###################### 

targetrange = 'Target Genes (TSS }10k)'  #Set if want to find genes 1, 5 or 10 kb from TF binding site

try:	#Read csv with TFs and their target genes. If not available, fetch the data and create a csv
	TF_gene_sets = pd.read_csv('../data/Transfactors/TF_gene_sets.csv', index_col = 'Antigen')
except IOError:
	TF_gene_sets = TF_targets(targetrange)
	TF_gene_sets.to_csv('../data/Transfactors/TF_gene_sets.csv')

dataset = pd.read_csv('../data/genes.raw.htseq2.tsv', sep='\t', index_col = 'Gene') #Read csv with mRNA expression data

dataset2 = pd.read_csv('../data/E-MTAB-2328.sdrf.tsv', sep='\t')					#Read csv with specifications of assays


chars = pd.DataFrame()
chars['assay'] = dataset2.loc[:,'Assay Name'].str.slice(stop=6)
chars['dev_stage'] = dataset2.loc[:,'Characteristics[developmental stage]']
chars['organ'] = dataset2.loc[:,'Characteristics[organism part]']
chars = chars.drop_duplicates()
chars = chars.set_index('assay')

dev_stage = chars.loc[:,'dev_stage']
organ = chars.loc[:,'organ']

datasetT = dataset.T
datasetT['dev_stage'] = dev_stage
datasetT['organ'] = organ
datasetT.set_index(['organ','dev_stage'], inplace=True)
expdata = datasetT.T

print(expdata)

TFdata = expdata.loc[['ENSMUSG00000093778','ENSMUSG00000093768'],'liver']
print(TFdata)
# ~ for TF in TF_gene_sets.index:				#Performing the PCA
	# ~ genes = TF_gene_sets.loc[TF, 'Genes']
	# ~ TFdata = dataset.loc[genes]
	# ~ res = my_pca(TFdata)
	# ~ PCA_per_TF[TF] = res

TF = TF_gene_sets.index[1]
PCA_per_TF = pd.DataFrame(index=dataset.columns)
PCA_per_TF[TF] = 100
# ~ print(PCA_per_TF)


# ~ TF = TF_gene_sets.index[1]
# ~ print(TF)
# ~ genes = TF_gene_sets.loc[TF, 'Genes']
# ~ print(genes)


# ~ print(test)

# ~ test2 = pd.read_csv('../data/Transfactors/TF_gene_sets.csv', index_col = 'Antigen')
# ~ print(test2)
# ~ print(test.loc['aha-1'])
# ~ real_gene_names =pd.read_csv('../data/real_gene_names.txt', sep='\t')
# ~ real_gene_names = real_gene_names.set_index('Gene stable ID')

# ~ genes_components_per_TF = {}





# ~ #chip = read_chip_atlas('../data/Transfactors/Acaa2.10.tsv')
# ~ exp_data = pd.read_csv('../data/genes.raw.htseq2.tsv', sep='\t') #load in expression data, seperated by tab
# ~ dataset = exp_data.set_index('Gene')


# ~ pathwaydata = dataset.loc[genes]
# ~ pathwaydata = pathwaydata.dropna(axis = 0, how = 'all') #has to be done so the lists match, this makes the dropna in my_pca function obsolete
# ~ presentgenes = pathwaydata.index.values.tolist()
# ~ res, components = my_pca(pathwaydata)
# ~ pca_per_pathway[pathway] = res
# ~ components = components.tolist()[0]
# ~ pathwayname = reactome.loc[pathway, 'pathway_name']
# ~ innerdict = {}
# ~ for i in range(0, len(presentgenes)):
	# ~ component = components[i]
	# ~ gene = genes[i]
        # ~ if gene in real_gene_names.index:
            # ~ real_name = real_gene_names.loc[gene, "Gene name"]
            # ~ innerdict[real_name] = component
        # ~ else:
			# ~ innerdict[gene] = component
	# ~ sorted_innerdict = sorted(innerdict.items(), key = operator.itemgetter(1), reverse = True)
	# ~ genes_components_per_pathway[pathwayname] = sorted_innerdict



# ~ print(chip)
# ~ pca_per_TF = pd.DataFrame
# ~ for  TF in chip.index:
	# ~ genes = chip.loc[TF, 'genes']
	# ~ TFdata = dataset.loc[genes]
# ~ pc = my_pca(dataset)
# ~ print(pc)
