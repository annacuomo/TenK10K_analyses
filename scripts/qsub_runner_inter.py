import os
import pandas as pd

if __name__ == '__main__':

    inputs_folder = "/share/ScratchGeneral/anncuo/OneK1K/input_files_CellRegMap/"
    fvf_filename = inputs_folder+"fvf_Bcell_eqtls.csv"
    #fvf_filename = inputs_folder+"fvf_Monocyte_eqtls.csv"
    fvf = pd.read_csv(fvf_filename, index_col = 0)

    print("prepare job to submit")
    qsub = "qsub -cwd -l mem_requested=600G -q short.q -r yes -N run_crm_inter_py -o stdout_run_crm_mono -e stderr_run_crm_mono -m ae -M a.cuomo@garvan.org.au -b y"
    
    
    for j in range(22):

        chrom = j+1
        if chrom not in [1, 2, 6, 7, 11, 12, 17, 19]: continue # ad hoc in case some chromosomes are already finished
        genes = fvf[fvf['chrom']==int(chrom)]['feature'].unique()

        for i in range(len(genes)):

            py = f"/share/ScratchGeneral/anncuo/jupyter/conda_notebooks/envs/cellregmap_notebook/bin/python run_interaction_Bcells_one_gene.py {chrom} {i}"
            # py = f"/share/ScratchGeneral/anncuo/jupyter/conda_notebooks/envs/cellregmap_notebook/bin/python run_interaction_Monocytes_sex_one_gene.py {chrom} {i}"
            cmd = f"{qsub} \"{py}\" "
            print(cmd)
            os.system(cmd)

