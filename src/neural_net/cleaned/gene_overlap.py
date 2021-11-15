import pandas as pd

blist = pd.read_csv("boruta.csv", names=['genes'])
ref_list = pd.read_csv("cancer_gene_census.csv", names=['genes'])

blist = set(blist['genes'])
ref_list = set(ref_list['genes'])

print(*blist.intersection(ref_list), sep='\n')
print(len(blist.intersection(ref_list)))

# for idx, item in enumerate(ref_list):
# 	if item in full_list['Gene Symbol']:
# 		print(full_list['Gene Symbol'][idx])
