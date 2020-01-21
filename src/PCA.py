import pandas as pd
import numpy as np
import mygene
from sklearn.decomposition import PCA
from sklearn import preprocessing


def pca_copy(df, n_pc=1, normalize=True):
	df = df.dropna(axis = 0, how = 'all') # remove rows without any values
	x = df.values.T #set x as transpose of only the numerical values
	
	if normalize:	
		x2 = preprocessing.scale(x)	#Standardize the data (center to mean and scale to unit variance
	else:
		x2 = x
	pca = PCA(n_components = n_pc)	#Set PCA parameters
	pca.fit(x2)	#Perform PCA
	x3 = pca.fit_transform(x2)	#fit and transform data (why?)
	out_df = pd.DataFrame(x3.transpose(), index=list(range(1,n_pc+1)), columns=df.columns) #Create dataframe of traspose of output with index set as nr of PC+1 (?)
	out_df = out_df.transpose() #Transpose data (why?)
	return out_df


exp_data = pd.read_csv('../data/genes.raw.htseq2.tsv', sep='\t') #load in expression data, seperated by tab
dataset = exp_data.set_index('Gene')

pc = pca_copy(dataset)
print(pc)
