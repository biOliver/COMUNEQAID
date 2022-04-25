import os

scopID = snakemake.params[0]
comID = snakemake.params[1]
comID_list = snakemake.params[1].split(',')
userID = snakemake.params[2]

for comID in comID_list:
    uniqueID = f'{scopID}_{comID}'
    path_share = f'/projects/{userID}/COMUNEQAID/outs/{scopID}'
    path_log = f'{path_share}/scRNAseq/04_Log/{comID}'

    os.system(f'cp /projects/SCOP/resources/non-species-specific/workflow-description/current/* /projects/{userID}/COMUNEQAID/outs/{scopID}/scRNAseq/04_Log/{comID}/.')
    os.system(f'rm -r /projects/{userID}/COMUNEQAID/manage-dir/tmp-data/{uniqueID}')
    os.system(f'cp logs/* {path_log}/.')

os.system(f'rm logs/0*')
os.system(f'touch {snakemake.output}')
