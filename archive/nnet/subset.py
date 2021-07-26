def make_set(fname):
	genes = set()
	with open(fname, 'r') as fin:
		for i in fin:
			genes.add(i.strip())
	return genes

if __name__ == '__main__':
	b1 = 'cleaned/boruta-99-99-0.05.csv'
	b2 = 'cleaned/boruta-99-25-0.05.csv'
	b3 = 'cleaned/boruta-99-10-0.05.csv'
	first = make_set(b1)
	second = make_set(b2)
	third = make_set(b3)
	intersection = set.intersection(first, second, third)
	print(intersection)
