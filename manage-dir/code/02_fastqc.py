# Imports
import os

# Load run-specific data
scopID = snakemake.params[0]
comID_list = snakemake.params[1].split(',')
userID = snakemake.params[2]

# Functions
def unpackInfo(fname):
    tmp = open(path_tmp + fname, 'r').read().rstrip('\n')
    return(tmp)

# Directory vars
path_share = f'/projects/{userID}/COMUNEQAID/outs/{scopID}'
path_fastq = f'{path_share}/scRNAseq/02_FASTQ'
path_qc = f'{path_share}/scRNAseq/00_QC'

for comID in comID_list:

    # Updated directory vars

    uniqueID = f'{scopID}_{comID}'
    path_tmp = f'/projects/{userID}/COMUNEQAID/manage-dir/tmp-data/{uniqueID}' #dirData

    #path_sample_sheets = f'{path_tmp}/sample-sheets'

    # Run vars
    #seqID   = unpackInfo('/tmp_seqID.txt')
    #flowID   = unpackInfo('/tmp_flowID.txt')
    #protocol   = unpackInfo('/tmp_protocol.txt')
    workflow = unpackInfo('/tmp_workflow.txt')

    bcl_dict_path = open(f'{path_tmp}/pin-data-10x.pkl', 'rb')
    pinDict_10x = pickle.load(bcl_dict_path)
    bcl_dict_path.close()

    if (workflow == '10x + HTO'):
        bcl_dict_path = open(f'{path_tmp}/pin-data-hto.pkl', 'rb')
        pinDict_hto = pickle.load(bcl_dict_path)
        bcl_dict_path.close()

    seq_names = []

    for pin, (indices,seq_name) in pinDict_10x.items():
        seq_names.append(seq_name)

    for pin, (indices,seq_name) in pinDict_hto.items():
        seq_names.append(seq_name)

    seq_names = list(set(seq_names))

    for seq_name in seq_names:
        print(f'{seq_name}')

        flowID = seq_name.split('_')[3][1:]

        tmp_path_fastqs = f'{path_fastq}/{seq_name}/fastq-path/{flowID}'
        tmp_path_qc = f'{path_qc}/fastQC/{seq_name}'

        #print(tmp_path_fastqs)
        #print(tmp_path_qc)


        if os.path.isdir(f'{tmp_path_qc}'):
            print('##\t-\tpipestance exists - skipping')

        else:
            os.system(f'mkdir -p {tmp_path_qc}')
            os.system(f'time /tools/fastqc/0.11.5/fastqc {tmp_path_fastqs}/*.fastq.gz --threads 1 --outdir {tmp_path_qc}')

os.system(f'touch {snakemake.output[0]}')
