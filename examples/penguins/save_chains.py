import reticulate
import numpy as np
import csv

def save_chains(m,n,p, rep):
	for i in range(rep):
		if i%10 == 0:
			print(i)
		for j in range(m):
			X = reticulate.get_samples(n)
			filename = 'Mchains/m_'+str(m)+"_n_"+str(n)+'_/chain_m_'+"_"+str(i)+"_"+str(j)+'.csv'
			with open(filename,"w+") as my_csv:
				csvWriter = csv.writer(my_csv,delimiter=',')
				csvWriter.writerows(X)

m = [5]
n = 100000
p = 29
rep = 10
for i in m:
	save_chains(i,n,p,rep)