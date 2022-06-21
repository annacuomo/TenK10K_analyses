import os
import pandas as pd

if __name__ == '__main__':

    revision_folder = "/hps/nobackup2/stegle/users/acuomo/all_scripts/struct_LMM2/sc_endodiff/debug_May2021/REVISION/"
    fvf_filename = revision_folder+"/CRM_interaction_chr20/fvf.csv"
    fvf = pd.read_csv(fvf_filename, index_col = 0)

    genes = fvf['feature'].unique()

    print("prepare job to submit")
    bsub = "bsub -R \"rusage[mem=80000]\" -M 80000"
    flags = "MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE"

    for i in range(len(genes)):

        gene = genes[i]
        n = fvf[fvf['feature']==gene].shape[0]

        for j in range(0,n,10):

            py = f"python interaction_test_for_10_snp_gene_pairs.py {i} {j}"
            cmd = f"{bsub} \"{flags} {py}\""
            print(cmd)
            os.system(cmd)

