import os
import glob
import pandas as pd
import numpy as np

def smartAppend(table,name,value):
    """ helper function for appending in a dictionary """
    if name not in table.keys():
        table[name] = []
    table[name].append(value)

path_results = "/share/ScratchGeneral/anncuo/OneK1K/CRM_association/CD4pos_Tcells/"

if __name__ == '__main__':

    fname = os.path.join(path_results,"*.tsv")
    files = glob.glob(fname)
    # print (files)
    x = 0
    table = {}

    for file in files:
        #print(file)
        x += 1
        if x%500 == 0: 
            print(x)
        df = pd.read_csv(file, index_col=0)
        nsnps = int(len(df))
        if nsnps==0:
            continue
        line = str(file).split("/")
        # print(line)
        gene = str(line[-1]).split(".")[0]
        # print(gene)
        chrom = df['chrom'].values[0]
        # print(chrom)
        for i in range(nsnps):
            temp = {}
            temp['gene'] = gene
            temp['n_snps'] = nsnps
            temp['snp_id'] = df['variant'].values[i]
            temp['pv_raw'] = df['pv'].values[i]
            temp['pv_Bonf'] = nsnps * temp['pv_raw']
            if temp['pv_Bonf']>1: temp['pv_Bonf'] = 1
            if temp['pv_Bonf']<0: temp['pv_Bonf'] = 0

        for key in temp.keys():
            smartAppend(table, key, temp[key])


    for key in table.keys():
        table[key] = np.array(table[key])

    df = pd.DataFrame.from_dict(table)
    outfile = "summary.csv" 
    myp = os.path.join(path_results, outfile)
    df.to_csv(myp)
