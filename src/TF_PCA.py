import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
import time
#import mygene		#for changing ensembl names to "normal" gene names

def my_pca(df, n_pc=1, normalize=True):	#
	df = df.dropna(axis = 0, how = 'all') # remove rows without any values
	x = df.values.T #set x as transpose of only the numerical values
	if normalize:	
		x2 = preprocessing.scale(x)	#Standardize the data (center to mean and scale to unit variance)
	else:
		x2 = x
	pca = PCA(n_components = n_pc)	#Set PCA parameters
	pca.fit(x2)	
#	my_pca.pca = pca	#Used for Patrik's dictionary method for components
	x3 = pca.fit_transform(x2)
	out_df = pd.DataFrame(x3.transpose(), index=list(range(1,n_pc+1)), columns=df.columns) #Create dataframe of traspose of output with index set as nr of PC+1 (?)
	out_df = out_df.transpose()
	return out_df

def TF_targets(TFrange): #Function for creating csv with all TFs and their target genes. For each TF, identifies urls for target gene data and fetches the info and puts it into a list that is then converted to a dataframe and given as output  
	chip = pd.read_csv('../data/Transfactors/chip_atlas_analysis_list.csv').loc[:,['Antigen','Target Genes (TSS }1k)','Target Genes (TSS }5k)', 'Target Genes (TSS }10k)']]
	chip = chip.dropna()
	chip = chip.set_index('Antigen')
	print(chip.loc[0:10,2])	#To go through full dataset, change loc[0:10,2] to loc[:,TFrange]
	TF_gene_list = []
	for url in chip.loc[0:10,2]:	#To go through full dataset, change loc[0:10,2] to loc[:,TFrange]	
		TF_gene_set = pd.read_csv(url, sep='\t')
		genes = TF_gene_set['Target_genes'].tolist()
		TF_gene_list.append(genes)
	TF_gene_sets = pd.DataFrame({'Genes':TF_gene_list}, index=chip.index[0:10])
	return TF_gene_sets

start = time.time()		#testing the runtime, can be removed
targetrange = 'Target Genes (TSS }10k)'  #Set if want to find genes 1, 5 or 10 kb from TF binding site

try:	#Try to read data w. TFs and their target genes. If not available, fetch the data and create a csv
	TF_gene_sets = pd.read_csv('../data/Transfactors/TF_gene_sets.csv', index_col = 'Antigen')
except IOError:
	TF_gene_sets = TF_targets(targetrange)
	TF_gene_sets.to_csv('../data/Transfactors/TF_gene_sets.csv')

dataset = pd.read_csv('../data/genes.raw.htseq2.tsv', sep='\t', index_col = 'Gene') #Read data w. mRNA expression values

end = time.time()		#testing the runtime, can be removed
print(end - start)

# ~ for TF in TF_gene_sets.index:
	# ~ genes = TF_gene_sets.loc[TF, 'Genes']
	# ~ TFdata = dataset.loc[genes]
	


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
