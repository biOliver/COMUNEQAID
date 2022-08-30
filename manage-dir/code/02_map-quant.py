# Imports
import os
import glob
import pandas as pd
import logging

# Init
logging.basicConfig(level=logging.INFO,
                    filename=snakemake.log[0],
                    format='%(message)s')
                    
logging.getLogger().addHandler(logging.StreamHandler())

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
path_res = f'/projects/SCOP/resources'

for comID in comID_list:
    uniqueID = f'{scopID}_{comID}'
    
    path_tmp = f'/projects/{userID}/COMUNEQAID/manage-dir/tmp-data/{uniqueID}'
    path_outs = f'{path_share}/scRNAseq/03_PipelineOut/{comID}'

    # Run vars
    organism    = unpackInfo('/tmp_organism.txt')
    workflow    = unpackInfo('/tmp_workflow.txt')

    bcl_dict_path = open(f'{path_tmp}/pin-data-10x.pkl', 'rb')
    pinDict_10x = pickle.load(bcl_dict_path)
    bcl_dict_path.close()

    indexDict_10x = dict()

    for pin, (indices,seq_name) in pinDict_10x.items():
        seq_year = int(seq_name[0:2])

        if seq_year <= 20:
            prefix = 'SI-GA-'
        elif seq_year >= 21:
            prefix = 'SI-TT-'

        for index in indices:
            indexDict_10x[prefix+index] = set()

    for pin, (indices,seq_name) in pinDict_10x.items():
        seq_year = int(seq_name[0:2])

        if seq_year <= 20:
            prefix = 'SI-GA-'
        elif seq_year >= 21:
            prefix = 'SI-TT-'

        for index in indices:
            indexDict_10x[prefix+index] = indexDict_10x[prefix+index].union({seq_name})

    if (workflow == '10x + HTO'):
        bcl_dict_path = open(f'{path_tmp}/pin-data-hto.pkl', 'rb')
        pinDict_hto = pickle.load(bcl_dict_path)
        bcl_dict_path.close()

        indexDict_hto = dict()

        for pin, (indices,seq_name) in pinDict_hto.items():
            for index in indices:
                indexDict_hto[index] = set()

        for pin, (indices,seq_name) in pinDict_hto.items():
            for index in indices:
                indexDict_hto[index] = indexDict_hto[index].union({seq_name})

    whitelist = f'{path_res}/non-species-specific/whitelist/3M-february-2018.txt'

    if organism == 'Mouse':
        path_ref     = f'{path_res}/mus-musculus/sc-ref/mouse_splici'
        salmon_index = f'{path_ref}/af_tutorial_splici/mm10_splici_idx'
        T2G          = f'{path_ref}/refdata-gex-mm10-2020-A/t2g.tsv'
        T3G          = f'{path_ref}/refdata-gex-mm10-2020-A/t2g_3col.tsv'

    elif organism == 'Human':
        path_ref     = f'{path_res}/homo-sapiens/sc-ref/human_splici'
        salmon_index = f'{path_ref}/af_tutorial_splici/grch38_splici_idx'
        T2G          = f'{path_ref}/refdata-gex-GRCh38-2020-A/t2g.tsv'
        T3G          = f'{path_ref}/refdata-gex-GRCh38-2020-A/t2g_3col.tsv'

    elif organism == 'Rat':
        path_ref     = f'{path_res}/rattus-norvegicus/sc-ref/rat_splici'
        salmon_index = f'{path_ref}/af_tutorial_splici/rnor_splici_idx'
        T2G          = f'{path_ref}/Rnor_6.0/t2g.tsv'
        T3G          = f'{path_ref}/Rnor_6.0/t2g_3col.tsv'
    
    elif organism == 'Rhesus':
        path_ref     = f'{path_res}/macaca-mulatta/sc-ref/rhesus_splici'
        salmon_index = f'{path_ref}/af_tutorial_splici/mmul_10_splici_idx'
        T2G          = f'{path_ref}/Mmul_10/t2g.tsv'
        T3G          = f'{path_ref}/Mmul_10/t2g_3col.tsv'

    elif organism == 'Mouse - optimized':
        path_ref     = f'{path_res}/mus-musculus/sc-ref/mouse_pool'
        salmon_index = f'{path_ref}/af_tutorial_splici/mouse-optimized_idx'
        T2G          = f'{path_ref}/mouse_mm10_optimized_v1/t2g.tsv'
        T3G          = f'{path_ref}/mouse_mm10_optimized_v1/t2g_3col.tsv'

    ################################################################################
    ################################################################################
    ################################################################################

    logging.info('#'*80)
    logging.info('#####                                                                      #####')
    logging.info('#####                      Mapping and quantification                      #####')
    logging.info('#####                                                      *** / *****     #####')
    logging.info('#'*80)
    logging.info('#####\n##')
    logging.info('##\tMapping and quantifying 10x libraries:')

    for index, seq_names in indexDict_10x.items():

        index_list = index.split(',')
    
        if len(index_list) == 1:
            logging.info(f'##\t-\tproccessing index:\t{index}')
        
            path_quant_out = f'{path_outs}/10x/{index}'

        elif len(index_list) > 1:
            indices_names_join = '-'.join(index_list)
            logging.info(f'##\t-\tproccessing indices:\t{indices_names_join}')

            path_quant_out = f'{path_outs}/10x/{indices_names_join}'
    
        seq_names_join = ', '.join(seq_names)
        logging.info(f'##\t\t(from {seq_names_join})')
        
    
        if os.path.isdir(f'{path_quant_out}'):
            logging.info('##\t-\tpipestance exists - skipping')
    
        else:
            files_read_1 = ''
            files_read_2 = ''

            for seq_name in seq_names:
        
                flowID = seq_name.split('_')[3][1:]
       
                for index_i in index_list:

                    files_read_1 += f'{path_fastq}/{seq_name}/fastq-path/{flowID}/*{index_i}*R1*.fastq.gz '
                    files_read_2 += f'{path_fastq}/{seq_name}/fastq-path/{flowID}/*{index_i}*R2*.fastq.gz '



            os.system(f'mkdir -p {path_quant_out}')

            logging.info(f'##\tmapping with alevin')
            os.system(f'/tools/anaconda/envs/blp841/salmon-new/bin/salmon alevin \
                        -i {salmon_index} -p {snakemake.threads} -l ISR --chromiumV3 --sketch \
                        -1 {files_read_1} -2 {files_read_2} \
                        -o {path_quant_out}/map --tgMap {T2G} 2>&1 | tee {path_quant_out}/alevin-fry.log -a')

            logging.info(f'##\tgenerating permit list')
            os.system(f'/tools/anaconda/envs/blp841/salmon-new/bin/alevin-fry generate-permit-list -d both -u {whitelist} \
                        -i {path_quant_out}/map -o {path_quant_out}/quant 2>&1 | tee {path_quant_out}/alevin-fry.log -a')

            logging.info(f'##\tcollating')
            os.system(f'/tools/anaconda/envs/blp841/salmon-new/bin/alevin-fry collate -t 1 -i {path_quant_out}/quant -r {path_quant_out}/map 2>&1 | tee {path_quant_out}/alevin-fry.log -a')

            logging.info(f'##\tquantifying')
            os.system(f'/tools/anaconda/envs/blp841/salmon-new/bin/alevin-fry quant -t {snakemake.threads} -i {path_quant_out}/quant \
                        -o {path_quant_out}/res --tg-map {T3G} --resolution cr-like --use-mtx 2>&1 | tee {path_quant_out}/alevin-fry.log -a')

            count_genes = 0
            count_cells = 0

            with open(f'{path_quant_out}/res/alevin/quants_mat_cols.txt', 'r') as f:
                for count_genes, line in enumerate(f):
                    pass

            with open(f'{path_quant_out}/res/alevin/quants_mat_rows.txt', 'r') as f:
                for count_cells, line in enumerate(f):
                    pass

            count_genes += 1
            count_cells += 1

            with open(f'{path_quant_out}/res/meta_info.json', 'w') as f:
                f.write('{')
                f.write('  "alt_resolved_cell_numbers": [],')
                f.write('  "dump_eq": false,')
                f.write(f'  "num_genes": {count_genes},')
                f.write(f'  "num_quantified_cells": {count_cells},')
                f.write('  "resolution_strategy": "CellRangerLike",')
                f.write('  "usa_mode": true')
                f.write('}')

    if (workflow == '10x + HTO'):
        logging.info('##\tMapping and quantifying HTO libraries:')

        for index, seq_names in indexDict_hto.items():

            index_list = index.split(',')

            if len(index_list) == 1:
                logging.info(f'##\t-\tproccessing index:\t{index}')

                path_quant_out = f'{path_outs}/hto/{index}'
                features_tsv_name = f'{index}-features.tsv'

            elif len(index_list) > 1:
                indices_names_join = '-'.join(index_list)
                logging.info(f'##\t-\tproccessing indices:\t{indices_names_join}')

                path_quant_out = f'{path_outs}/hto/{indices_names_join}'
                features_tsv_name = f'{indices_names_join}-features.tsv'

            seq_names_join = ', '.join(seq_names)
            logging.info(f'##\t\t(from {seq_names_join})')
        
            if os.path.isdir(f'{path_quant_out}'):
                logging.info('##\t-\tpipestance exists - skipping')

            else:
                files_read_1 = ''
                files_read_2 = ''

                for seq_name in seq_names:

                    flowID = seq_name.split('_')[3][1:]

                    for index_i in index_list:

                        files_read_1 += f'{path_fastq}/{seq_name}/fastq-path/{flowID}/*{index_i}*R1*.fastq.gz '
                        files_read_2 += f'{path_fastq}/{seq_name}/fastq-path/{flowID}/*{index_i}*R2*.fastq.gz '

                os.system(f'mkdir -p {path_quant_out}')
        
                os.system(f'cp {path_tmp}/hto-features/{features_tsv_name} {path_quant_out}')

                os.system('awk {0} {1}/{2} > {1}/t2g_hto.tsv'.format('\'{print $1\"\\t\"$1;}\'', path_quant_out, features_tsv_name))

                logging.info('##\tbuilding HTO index')
                os.system(f'/tools/anaconda/envs/blp841/salmon-new/bin/salmon index \
                            -t {path_quant_out}/{features_tsv_name} \
                            -i {path_quant_out}/hto_index --features -k 7 2>&1 | tee {path_quant_out}/alevin-fry.log -a')
    
                logging.info(f'##\tmapping with alevin')
                os.system(f'/tools/anaconda/envs/blp841/salmon-new/bin/salmon alevin -l ISR \
                            -i {path_quant_out}/hto_index -1 {files_read_1} -2 {files_read_2} \
                            --read-geometry 2[1-15] --bc-geometry 1[1-16] --umi-geometry 1[17-26] \
                            -o {path_quant_out}/map -p {snakemake.threads} --sketch 2>&1 | tee {path_quant_out}/alevin-fry.log -a')

                logging.info(f'##\tgenerating permit list')
                os.system(f'/tools/anaconda/envs/blp841/salmon-new/bin/alevin-fry generate-permit-list -d fw -i {path_quant_out}/map -o {path_quant_out}/quant -u {whitelist} 2>&1 | tee {path_quant_out}/alevin-fry.log -a')
            
                logging.info(f'##\tcollating')
                os.system(f'/tools/anaconda/envs/blp841/salmon-new/bin/alevin-fry collate -r {path_quant_out}/map -i {path_quant_out}/quant -t 1 2>&1 | tee {path_quant_out}/alevin-fry.log -a')
 
                logging.info(f'##\tquantifying')
                os.system(f'/tools/anaconda/envs/blp841/salmon-new/bin/alevin-fry quant -m {path_quant_out}/t2g_hto.tsv -i {path_quant_out}/quant -o {path_quant_out}/res \
                            -r cr-like -t {snakemake.threads} --use-mtx 2>&1 | tee {path_quant_out}/alevin-fry.log -a')

                with open(f'{path_quant_out}/res/alevin/quants_mat_cols.txt', 'r') as f:
                    for count_genes, line in enumerate(f):
                        pass

                with open(f'{path_quant_out}/res/alevin/quants_mat_rows.txt', 'r') as f:
                    for count_cells, line in enumerate(f):
                        pass

                count_genes += 1
                count_cells += 1

                with open(f'{path_quant_out}/res/meta_info.json', 'w') as f:
                    f.write('{\n')
                    f.write('  "alt_resolved_cell_numbers": [],\n')
                    f.write('  "dump_eq": false,\n')
                    f.write(f'  "num_genes": {count_genes},\n')
                    f.write(f'  "num_quantified_cells": {count_cells},\n')
                    f.write('  "resolution_strategy": "CellRangerLike",\n')
                    f.write('  "usa_mode": false\n')
                    f.write('}')

    logging.info('##\n###')
    logging.info('#'*80)

os.system(f'touch {snakemake.output}')
