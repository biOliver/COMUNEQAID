# Imports
import os
import pandas as pd
import logging

# Init
logging.basicConfig(level=logging.INFO,
                    filename=snakemake.log[0],
                    format='%(message)s')
                    
logging.getLogger().addHandler(logging.StreamHandler())

# Load run-specific data
scopID = snakemake.params[0]
comID = snakemake.params[1]
userID = snakemake.params[2]

uniqueID = f'{scopID}_{comID}'

# Functions
def unpackInfo(fname):
    tmp = open(path_tmp + fname, 'r').read().rstrip('\n')
    return(tmp)

# Directory vars
path_share = f'/projects/{userID}/COMUNEQAID/outs/{scopID}'
path_bcls  = f'{path_share}/scRNAseq/01_BCL'
path_fastq = f'{path_share}/scRNAseq/02_FASTQ'
path_tmp = f'/projects/{userID}/COMUNEQAID/manage-dir/tmp-data/{uniqueID}'
path_sample_sheets = f'{path_tmp}/sample-sheets'

# Run vars
workflow    = unpackInfo('/tmp_workflow.txt')

bcl_dict_path = open(f'{path_tmp}/pin-data-10x.pkl', 'rb')
pinDict_10x = pickle.load(bcl_dict_path)
bcl_dict_path.close()

if (workflow == '10x + HTO'):
    bcl_dict_path = open(f'{path_tmp}/pin-data-hto.pkl', 'rb')
    pinDict_hto = pickle.load(bcl_dict_path)
    bcl_dict_path.close()

################################################################################
################################################################################
################################################################################

print('#'*80)
print('#####                                                                      #####')
print('#####                      Demultiplexing FASTQ files                      #####')
print('#####                                                       ** / *****     #####')
print('#'*80)
print('#####\n##')
print('##\tDemultiplexing libraries:')

for pin, (indices,seq_name) in pinDict_10x.items():

    print(f'##\t-\t{seq_name}')
    
    flowID = seq_name.split('_')[3][1:]

    dir_bcl = f'{path_bcls}/{seq_name}'
    sample_sheet = f'{path_sample_sheets}/samplesheet_{pin}.csv'
    dir_out = f'{path_fastq}/{seq_name}'

    if os.path.isfile(f'{dir_out}/interop_path/IndexMetricsOut.bin'):
        print('##\t-\tpipestance exists - skipping')

    else:
        year = int(seq_name[0:2])

        if year <= 20:
            versionRNA = 'RNA v3.0'
            prefix = 'SI-GA-'
        elif year >= 21:
            versionRNA = 'RNA v3.1'
            prefix = 'SI-TT-'

        os.system(f'mkdir -p {dir_out}')
        os.chdir(dir_out)

        if (versionRNA == 'RNA v3.1'):
            os.system(f'bcl2fastq --use-bases-mask=Y28,I10,I10,Y90 \
            --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=10 \
            --mask-short-adapter-reads=10 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r {snakemake.threads} -p {snakemake.threads} -w {snakemake.threads} -R {dir_bcl} \
            --output-dir=fastq-path \
            --interop-dir=interop_path \
            --sample-sheet={sample_sheet} \
            --barcode-mismatches=1 \
            --no-lane-splitting 2>&1 | tee logs/bcl2fastq.log -a')

os.chdir(f'/projects/{userID}/COMUNEQAID/manage-dir')
os.system(f'touch {snakemake.output[0]}')

print('##\n###')
print('#'*80)
