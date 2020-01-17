import requests
import pandas as pd

#url='http://dbarchive.biosciencedbc.jp/kyushu-u/hg19/target/CLOCK.1.tsv'


df = pd.read_csv('../data/Transfactors/chip_atlas_analysis_list.csv')
urls1 = df.loc[:,'Target Genes (TSS }1k)']
urls5 = df.loc[:,'Target Genes (TSS }5k)']
urls10 = df.loc[:,'Target Genes (TSS }10k)']
urls1.dropna(axis=0, how='any', inplace=True)


#for url in urls1:
	#	response = requests.get(url) 
#except:
#	print('Fail')
#else:
#	print('Success')
#	print(response.headers)
#def get_data(path, url)
#	response = requests.get(url, stream=True)
#	print(r.headers.get('content-type')
