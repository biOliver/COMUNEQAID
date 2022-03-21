import os

scopID = snakemake.params[0]
comID = snakemake.params[1]
userID = snakemake.params[2]

uniqueID = f'{scopID}_{comID}'

path_share = f'/projects/{userID}/COMUNEQAID/outs/{scopID}' # <- Change this to your own workflow directory
path_log = f'{path_share}/scRNAseq/04_Log/{comID}'

#os.system(f'cp /projects/oliver/SCOaaP/resources/workflowDesc/snRNAseq_workflow_v1.0.pdf /projects/oliver/SCOaaP/outs/{scopID}/scRNAseq/04_Log/{runID}/.')
os.system(f'rm -r /projects/SCOP/pipelines/COMUNEQAID/COMUNEQAID-app/app-data/tmp-data/{uniqueID}')
os.system(f'mv logs/* {path_log}/.')
os.system(f'touch {snakemake.output}')
