import pandas as pd

blist = pd.read_csv("boruta.csv", names=['blist'])
ref_list = pd.read_csv("cancer_gene_census.csv", names=['ref_list'])

blist = set(blist['blist'])
ref_list = set(ref_list['ref_list'])

print(blist.intersection(ref_list))
print(len(blist.intersection(ref_list)))
