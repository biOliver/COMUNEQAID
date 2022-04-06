import os

scopID = snakemake.params[0]
comID = snakemake.params[1]
userID = snakemake.params[2]

uniqueID = f'{scopID}_{comID}'

path_share = f'/projects/{userID}/COMUNEQAID/outs/{scopID}'
path_log = f'{path_share}/scRNAseq/04_Log/{comID}'

os.system(f'cp /projects/SCOP/resources/non-species-specific/workflowDesc/current/* /projects/{userID}/COMUNEQAID/outs/{scopID}/scRNAseq/04_Log/{comID}/.')
os.system(f'rm -r /projects/{userID}/COMUNEQAID/manage-dir/tmp-data/{uniqueID}')
os.system(f'mv logs/* {path_log}/.')
os.system(f'touch {snakemake.output}')
