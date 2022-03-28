# Importis
import os
import re
import pandas as pd
import numpy as np
import logging
from bs4 import BeautifulSoup

# Init
logging.basicConfig(level=logging.INFO,
                    filename=snakemake.log[0],
                    format='%(message)s')
                    
logging.getLogger().addHandler(logging.StreamHandler())


outputDict = dict()

translate_index_10x_singleIndex = {
    'A1': ['GGTTTACT','CTAAACGG','TCGGCGTC','AACCGTAA'],
    'A2': ['TTTCATGA','ACGTCCCT','CGCATGTG','GAAGGAAC'],
    'A3': ['CAGTACTG','AGTAGTCT','GCAGTAGA','TTCCCGAC'],
    'A4': ['TATGATTC','CCCACAGT','ATGCTGAA','GGATGCCG'],
    'A5': ['CTAGGTGA','TCGTTCAG','AGCCAATT','GATACGCC'],
    'A6': ['CGCTATGT','GCTGTCCA','TTGAGATC','AAACCGAG'],
    'A7': ['ACAGAGGT','TATAGTTG','CGGTCCCA','GTCCTAAC'],
    'A8': ['GCATCTCC','TGTAAGGT','CTGCGATG','AACGTCAA'],
    'A9': ['TCTTAAAG','CGAGGCTC','GTCCTTCT','AAGACGGA'],
    'A10': ['GAAACCCT','TTTCTGTC','CCGTGTGA','AGCGAAAG'],
    'A11': ['GTCCGGTC','AAGATCAT','CCTGAAGG','TGATCTCA'],
    'A12': ['AGTGGAAC','GTCTCCTT','TCACATCA','CAGATGGG'],
    'B1': ['GTAATCTT','TCCGGAAG','AGTTCGGC','CAGCATCA'],
    'B2': ['TACTCTTC','CCTGTGCG','GGACACGT','ATGAGAAA'],
    'B3': ['GTGTATTA','TGTGCGGG','ACCATAAC','CAACGCCT'],
    'B4': ['ACTTCATA','GAGATGAC','TGCCGTGG','CTAGACCT'],
    'B5': ['AATAATGG','CCAGGGCA','TGCCTCAT','GTGTCATC'],
    'B6': ['CGTTAATC','GCCACGCT','TTACTCAG','AAGGGTGA'],
    'B7': ['AAACCTCA','GCCTTGGT','CTGGACTC','TGTAGAAG'],
    'B8': ['AAAGTGCT','GCTACCTG','TGCTGTAA','CTGCAAGC'],
    'B9': ['CTGTAACT','TCTAGCGA','AGAGTGTG','GACCCTAC'],
    'B10': ['ACCGTATG','GATTAGAT','CTGACTGA','TGACGCCC'],
    'B11': ['GTTCCTCA','AGGTACGC','TAAGTATG','CCCAGGAT'],
    'B12': ['TACCACCA','CTAAGTTT','GGGTCAAG','ACTGTGGC'],
    'C1': ['CCACTTAT','AACTGGCG','TTGGCATA','GGTAACGC'],
    'C2': ['CCTAGACC','ATCTCTGT','TAGCTCTA','GGAGAGAG'],
    'C3': ['TCAGCCGT','CAGAGGCC','GGTCAATA','ATCTTTAG'],
    'C4': ['ACAATTCA','TGCGCAGC','CATCACTT','GTGTGGAG'],
    'C5': ['CGACTTGA','TACAGACT','ATTGCGTG','GCGTACAC'],
    'C6': ['ATTACTTC','TGCGAACT','GCATTCGG','CAGCGGAA'],
    'C7': ['GTCTCTCG','AATCTCTC','CGGAGGGA','TCAGAAAT'],
    'C8': ['GTTGAGAA','AGATCTGG','TCGATACT','CACCGCTC'],
    'C9': ['GCGCAGAA','ATCTTACC','TATGGTGT','CGAACCTG'],
    'C10': ['TCTCAGTG','GAGACTAT','CGCTTAGC','ATAGGCCA'],
    'C11': ['GAGGATCT','AGACCATA','TCCTGCGC','CTTATGAG'],
    'C12': ['TCTCGTTT','GGCTAGCG','ATGACCGC','CAAGTAAA'],
    'D1': ['CACTCGGA','GCTGAATT','TGAAGTAC','ATGCTCCG'],
    'D2': ['TAACAAGG','GGTTCCTC','ATCATGCA','CCGGGTAT'],
    'D3': ['ACATTACT','TTTGGGTA','CAGCCCAC','GGCAATGG'],
    'D4': ['CCCTAACA','ATTCCGAT','TGGATTGC','GAAGGCTG'],
    'D5': ['CTCGTCAC','GATCAGCA','ACAACAGG','TGGTGTTT'],
    'D6': ['CATGCGAT','TGATATTC','GTGATCGA','ACCCGACG'],
    'D7': ['ATTTGCTA','TAGACACC','CCACAGGG','GGCGTTAT'],
    'D8': ['GCAACAAA','TAGTTGTC','CGCCATCG','ATTGGCGT'],
    'D9': ['AGGAGATG','GATGTGGT','CTACATCC','TCCTCCAA'],
    'D10': ['CAATACCC','TGTCTATG','ACCACGAA','GTGGGTGT'],
    'D11': ['CTTTGCGG','TGCACAAA','AAGCAGTC','GCAGTTCT'],
    'D12': ['GCACAATG','CTTGGTAC','TGCACCGT','AAGTTGCA'],
    'E1': ['TGGTAAAC','GAAAGGGT','ACTGCTCG','CTCCTCTA'],
    'E2': ['GTGGTACC','TACTATAG','ACAAGGTA','CGTCCCGT'],
    'E3': ['AGGTATTG','CTCCTAGT','TCAAGGCC','GATGCCAA'],
    'E4': ['TTCGCCCT','GGATGGGC','AATCAATG','CCGATTAA'],
    'E5': ['CATTAGCG','TTCGCTGA','ACAAGAAT','GGGCTCTC'],
    'E6': ['CTGCGGCT','GACTCAAA','AGAAACTC','TCTGTTGG'],
    'E7': ['CACGCCTT','GTATATAG','TCTCGGGC','AGGATACA'],
    'E8': ['ATAGTTAC','TGCTGAGT','CCTACGTA','GAGCACCG'],
    'E9': ['TTGTTTCC','GGAGGAGG','CCTAACAA','AACCCGTT'],
    'E10': ['AAATGTGC','GGGCAAAT','TCTATCCG','CTCGCGTA'],
    'E11': ['AAGCGCTG','CGTTTGAT','GTAGCACA','TCCAATGC'],
    'E12': ['ACCGGCTC','GAGTTAGT','CGTCCTAG','TTAAAGCA'],
    'F1': ['GTTGCAGC','TGGAATTA','CAATGGAG','ACCCTCCT'],
    'F2': ['TTTACATG','CGCGATAC','ACGCGGGT','GAATTCCA'],
    'F3': ['TTCAGGTG','ACGGACAT','GATCTTGA','CGATCACC'],
    'F4': ['CCCAATAG','GTGTCGCT','AGAGTCGC','TATCGATA'],
    'F5': ['GACTACGT','CTAGCGAG','TCTATATC','AGGCGTCA'],
    'F6': ['CGGAGCAC','GACCTATT','ACTTAGGA','TTAGCTCG'],
    'F7': ['CGTGCAGA','AACAAGAT','TCGCTTCG','GTATGCTC'],
    'F8': ['CATGAACA','TCACTCGC','AGCTGGAT','GTGACTTG'],
    'F9': ['CAAGCTCC','GTTCACTG','TCGTGAAA','AGCATGGT'],
    'F10': ['GCTTGGCT','AAACAAAC','CGGGCTTA','TTCATCGG'],
    'F11': ['GCGAGAGT','TACGTTCA','AGTCCCAC','CTATAGTG'],
    'F12': ['TGATGCAT','GCTACTGA','CACCTGCC','ATGGAATG'],
    'G1': ['ATGAATCT','GATCTCAG','CCAGGAGC','TGCTCGTA'],
    'G2': ['TGATTCTA','ACTAGGAG','CAGCCACT','GTCGATGC'],
    'G3': ['CCTCATTC','AGCATCCG','GTGGCAAT','TAATGGGA'],
    'G4': ['GCGATGTG','AGATACAA','TTTCCACT','CACGGTGC'],
    'G5': ['GAGCAAGA','TCTGTGAT','CGCAGTTC','ATATCCCG'],
    'G6': ['CTGACGCG','GGTCGTAC','TCCTTCTT','AAAGAAGA'],
    'G7': ['GGTATGCA','CTCGAAAT','ACACCTTC','TAGTGCGG'],
    'G8': ['TATGAGCT','CCGATAGC','ATACCCAA','GGCTGTTG'],
    'G9': ['TAGGACGT','ATCCCACA','GGAATGTC','CCTTGTAG'],
    'G10': ['TCGCCAGC','AATGTTAG','CGATAGCT','GTCAGCTA'],
    'G11': ['TTATCGTT','AGCAGAGC','CATCTCCA','GCGGATAG'],
    'G12': ['ATTCTAAG','CCCGATTA','TGGAGGCT','GAATCCGC'],
    'H1': ['GTATGTCA','TGTCAGAC','CACGTCGG','ACGACATT'],
    'H2': ['TAATGACC','ATGCCTTA','GCCGAGAT','CGTATCGG'],
    'H3': ['CCAAGATG','AGGCCCGA','TACGTGAC','GTTTATCT'],
    'H4': ['GCCATTCC','CAAGAATT','TTGCCGGA','AGTTGCAG'],
    'H5': ['CCACTACA','GATTCTGG','TGCGGCTT','ATGAAGAC'],
    'H6': ['TAGGATAA','CCTTTGTC','GTACGCGG','AGCACACT'],
    'H7': ['AGCTATCA','CATATAAC','TCAGGGTG','GTGCCCGT'],
    'H8': ['TTGTTGAT','GCTCAACC','CAAAGTGG','AGCGCCTA'],
    'H9': ['ACACTGTT','CAGGATGG','GGCTGAAC','TTTACCCA'],
    'H10': ['GTAATTGC','AGTCGCTT','CACGAGAA','TCGTCACG'],
    'H11': ['GGCGAGTA','ACTTCTAT','CAAATACG','TTGCGCGC'],
    'H12': ['GACAGCAT','TTTGTACA','AGGCCGTG','CCATATGC']
}

translate_index_10x_dualIndex = {
    'A1': ['GTAACATGCG','AGGTAACACT'], 'A2': ['GTGGATCAAA','CAGGGTTGGC'],
    'A3': ['CACTACGAAA','ATCAGTCTAA'], 'A4': ['CTCTAGCGAG','GATGAAGATA'],
    'A5': ['GTAGCCCTGT','ATAGATGCTC'], 'A6': ['TAACGCGTGA','GAAGTTAGGG'],
    'A7': ['TCCCAAGGGT','AAAGGTAGTA'], 'A8': ['CGAAGTATAC','CTCCAAGTTC'],
    'A9': ['AAGTGGAGAG','GTAACAGGAA'], 'A10': ['CGTGACATGC','TTTAGACCAT'],
    'A11': ['CGGAACCCAA','TCCTCGAATC'], 'A12': ['CACCGCACCA','ATTGACAGTC'],
    'B1': ['ACAGTAACTA','AACGAACTGT'], 'B2': ['TCTACCATTT','GACTCTCCCG'],
    'B3': ['CACGGTGAAT','TGTGACGAAC'], 'B4': ['GTAGACGAAA','ACCACACTAG'],
    'B5': ['TCGGCTCTAC','AGACCATCGG'], 'B6': ['AATGCCATGA','GCATTACGTA'],
    'B7': ['GCCTTCGGTA','AAATCGTTGG'], 'B8': ['GCACTGAGAA','TTCACGCATA'],
    'B9': ['TATTGAGGCA','CACTTACCTG'], 'B10': ['GCCCGATGGA','CTAGACGATT'],
    'B11': ['TCTTACTTGC','CTAGAGGTCA'], 'B12': ['CGTCAAGGGC','GAGTGACCTA'],
    'C1': ['TGCGCGGTTT','TTTATCCTTG'], 'C2': ['CAATCCCGAC','TACTACTCGG'],
    'C3': ['ATGGCTTGTG','CACAACATTC'], 'C4': ['TTCTCGATGA','GTGCCCGACA'],
    'C5': ['TCCGTTGGAT','GCGAGAACGT'], 'C6': ['ACGACTACCA','TTAGGGTCGT'],
    'C7': ['CGCGCACTTA','AGAATACAGG'], 'C8': ['GCTACAAAGC','AGGGCACGTG'],
    'C9': ['TATCAGCCTA','AGGACGAAAC'], 'C10': ['AGAATGGTTT','TCCCACCCTC'],
    'C11': ['ATGGGTGAAA','AATTCCCAAG'], 'C12': ['TCGTCAAGAT','CCTGAGTTGC'],
    'D1': ['TGCAATGTTC','TTCGACAAGC'], 'D2': ['TTAATACGCG','ACCCGAGGTG'],
    'D3': ['CCTTCTAGAG','TCGTTGTATT'], 'D4': ['GCAGTATAGG','GTGCACGGAA'],
    'D5': ['TGGTTCGGGT','CTCCTGCCAC'], 'D6': ['CCCAGCTTCT','GTTTGGTGTC'],
    'D7': ['CCTGTCAGGG','GTTACGGGCT'], 'D8': ['CGCTGAAATC','GCAGACACCT'],
    'D9': ['TGGTCCCAAG','ACGCCAGAGG'], 'D10': ['ATGCGAATGG','CGACACTTGT'],
    'D11': ['CGAATATTCG','TTGCTTCCAG'], 'D12': ['GAATTGGTTA','CTACTAGAGT'],
    'E1': ['TTATTCGAGG','AGCAGGACAG'], 'E2': ['ATGGAGGGAG','AATGGGTTAT'],
    'E3': ['ACCAGACAAC','CCTAGTTCCT'], 'E4': ['AACCACGCAT','TAACCTGAAT'],
    'E5': ['CGCGGTAGGT','CAACATCCTG'], 'E6': ['TTGAGAGTCA','CTACCAGGTT'],
    'E7': ['GTCCTTCGGC','CTGTGCATGA'], 'E8': ['GAGCAAGGGC','CCAAGTCAAT'],
    'E9': ['TGTCCCAACG','TGGACATCGA'], 'E10': ['CACAATCCCA','TTGTGGATAT'],
    'E11': ['TCCGGGACAA','TGGCATTCAC'], 'E12': ['CGTCCACCTG','GTCATGAATG'],
    'F1': ['AAGATTGGAT','AAATCCCGCT'], 'F2': ['AAGGGCCGCA','GAGGAATCAG'],
    'F3': ['GAGAGGATAT','CCCATTTCAA'], 'F4': ['CCCACCACAA','AAGCGGAGGT'],
    'F5': ['CGGCTGGATG','GTGCTTATCA'], 'F6': ['TTGCCCGTGC','AATCTCACGC'],
    'F7': ['AATGTATCCA','TAAGCTCATT'], 'F8': ['CTCCTTTAGA','GAGCTATGTC'],
    'F9': ['GTCCCATCAA','GTCACGTTCG'], 'F10': ['CCGGCAACTG','TGTTAAACCG'],
    'F11': ['TTCACACCTT','GTGTACACTA'], 'F12': ['GAGACGCACG','ATGTTCATAG'],
    'G1': ['TGTAGTCATT','TACGATCAAG'], 'G2': ['CATGTGGGTT','TAAAGGAATC'],
    'G3': ['ATGACGTCGC','ATCCTGACCT'], 'G4': ['GCGCTTATGG','CTAGCCAGGC'],
    'G5': ['ATAGGGCGAG','ACTCGATGCA'], 'G6': ['GCGGGTAAGT','CTTAGTGCTA'],
    'G7': ['GTTTCACGAT','TTTGGCCGAA'], 'G8': ['TAAGCAACTG','TTGAGTATAG'],
    'G9': ['CCGGAGGAAG','AACATCCGCA'], 'G10': ['ACTTTACGTG','AGGGCGTTCA'],
    'G11': ['GATAACCTGC','GTTTCTAATG'], 'G12': ['CTTGCATAAA','AAGCCCTGAT'],
    'H1': ['ACAATGTGAA','TAACGGTACG'], 'H2': ['TAGCATAGTG','GACAGAGCCG'],
    'H3': ['CCCGTTCTCG','CCAATCCGTC'], 'H4': ['AGTTTCCTGG','CTGTGTGGCA'],
    'H5': ['AGCAAGAAGC','AGAAACACAA'], 'H6': ['CCTATCCTCG','GTTAGTATTC'],
    'H7': ['ACCTCGAGCT','ATCGAACACA'], 'H8': ['ATAAGGATAC','CCCTATCTAT'],
    'H9': ['AGAACTTAGA','AAAGGACTCG'], 'H10': ['TTATCTAGGG','TAGAGCCTTT'],
    'H11': ['ACAATCGATC','CATTCCGTCA'], 'H12': ['TGATGATTCA','CGACTCCTAC']}

translate_index_hto = {
    'D701': ['ATTACTCGAT','AGATCTCGGT'], 'D702': ['TCCGGAGAAT','AGATCTCGGT'],
    'D703': ['CGCTCATTAT','AGATCTCGGT'], 'D704': ['GAGATTCCAT','AGATCTCGGT'],
    'D705': ['ATTCAGAAAT','AGATCTCGGT'], 'D706': ['GAATTCGTAT','AGATCTCGGT'],
    'D707': ['CTGAAGCTAT','AGATCTCGGT'], 'D708': ['TAATGCGCAT','AGATCTCGGT'],
    'D709': ['CGGCTATGAT','AGATCTCGGT'], 'D710': ['TCCGCGAAAT','AGATCTCGGT'],
    'D711': ['TCTCGCGCAT','AGATCTCGGT'], 'D712': ['AGCGATAGAT','AGATCTCGGT']}

translate_bc_hto = {
    'A0451': 'TTCCTGCCATTACTA','A0452': 'CCGTACCTCATTGTT','A0453': 'GGTAGATGTCCTCAG',
    'A0454': 'TGGTGTCATTCTTGA','A0455': 'ATGATGAACAGCCAG','A0456': 'CTCGAACGCTTATCG',
    'A0457': 'CTTATCACCGCTCAA','A0458': 'TGACGCCGTTGTTGT','A0459': 'GCCTAGTATGATCCA',
    'A0460': 'AGTCACAGTATTCCA','A0461': 'TTGGCCTTTGTATCG','A0462': 'AACGCCAGTATGAAC',
    'A0463': 'GCTTCCGTATATCTG','A0464': 'AGCAACTCACTCTTT','A0465': 'CACTAGAGCTTGGAT',
    'A0301': 'ACCCACCAGTAAGAC','A0302': 'GGTCGAGAGCATTCA','A0303': 'CTTGCCGCATGTCAT',
    'A0304': 'AAAGCATTCTTCACG','A0305': 'CTTTGTCTTTGTGAG','A0306': 'TATGCTGCCACGGTA',
    'A0251': 'GTCAACTCTTTAGCG','A0252': 'TGATGGCCTATTGGG','A0253': 'TTCCGCCTCTCTTTG',
    'A0254': 'AGTAAGTTCAGCGTA'}

colnames_samp = ['Name','Tissue/Region','Index (10x)','Index (HTO)','BCL PIN (10x)','BCL PIN (HTO)','HTO','Loaded Cells (Sample)','Lane','Various Info']

# Load run-specific data
scopID = snakemake.params[0]
comID = snakemake.params[1]
userID = snakemake.params[2]

uniqueID = f'{scopID}_{comID}'

path_share = f'/projects/{userID}/COMUNEQAID/outs/{scopID}'
path_bcls = f'{path_share}/scRNAseq/01_BCL/'
path_tmp = f'/projects/{userID}/COMUNEQAID/manage-dir/tmp-data/{uniqueID}'
path_sample_sheets = f'{path_tmp}/sample-sheets'
path_hto_features = f'{path_tmp}/hto-features'
path_templates = f'/projects/SCOP/resources/non-species-specific/templates'
path_app_data = f'/projects/SCOP/pipelines/COMUNEQAID/COMUNEQAID-app/app-data/'

tmp_nume = pd.read_csv(f'{path_app_data}/settings/{uniqueID}_nume.csv')
tmp_sele = pd.read_csv(f'{path_app_data}/settings/{uniqueID}_sele.csv')
tmp_text = pd.read_csv(f'{path_app_data}/settings/{uniqueID}_text.csv')

outputDict['tmp_samples'] = tmp_nume['value'][0]

outputDict['tmp_workflow'] = tmp_sele['value'][0]
outputDict['tmp_seqType'] = tmp_sele['value'][1]
outputDict['tmp_organism'] = tmp_sele['value'][2]

outputDict['tmp_scopID'] = tmp_text['value'][0]
outputDict['tmp_comID'] = tmp_text['value'][1]

tmp_samp = pd.read_csv(f'{path_app_data}/tables/{uniqueID}_samp.csv',
                       names = colnames_samp, header = 0,
                       dtype = {'BCL PIN (10x)': str,'BCL PIN (HTO)': str})

outputDict['tmp_samp'] = tmp_samp.to_string()

if not os.path.exists(path_share):
    os.makedirs(f'{path_share}/scRNAseq')
    os.system(f'cp -r /projects/SCOP/resources/non-species-specific/templates/template_SCOP/* {path_share}/scRNAseq/.')

files_bcls = os.listdir(path_bcls)

os.system(f'mkdir -p {path_sample_sheets}')

if (outputDict['tmp_workflow'] == '10x + HTO'):
    os.system(f'mkdir -p {path_hto_features}')

################################################################################
################################################################################
################################################################################

logging.info('#'*80)
logging.info('####################                                        ####################')
logging.info('##########                         COMUNEQAID                         ##########')
logging.info('##########   ( COmbine MUltiple single Nuclei sEQuencing runs AID )   ##########')
logging.info('##########                                                            ##########')
logging.info('####################                                        ####################')
logging.info('#'*80)
logging.info('#')
logging.info('#')
logging.info('#')
logging.info('#'*80)
logging.info('#####                                                                      #####')
logging.info('#####                         Compiling user input                         #####')
logging.info('#####                                                        * / *****     #####')
logging.info('#'*80)
logging.info('#####\n##')
logging.info(f'##\tSCOP ID:\t{scopID}')
logging.info(f'##\tCOMUNEQA ID:\t{comID}')
logging.info('##')

tmp_samples = outputDict['tmp_samples']
tmp_workflow = outputDict['tmp_workflow']
tmp_seqType = outputDict['tmp_seqType']
tmp_organism = outputDict['tmp_organism']

logging.info(f'##\t\tnumber of samples:\t{tmp_samples}')
logging.info(f'##\t\tworkflow:\t\t{tmp_workflow}')
logging.info(f'##\t\ttype:\t\t\t{tmp_seqType}')
logging.info(f'##\t\torganism:\t\t{tmp_organism}')
logging.info(f'##\t\tavailable sequencing runs:')
for file_bcl in files_bcls:
    logging.info(f'##\t\t-\t{file_bcl}')
logging.info('##')
logging.info('#####\n##')
logging.info(f'##\tCreating pool table...')
if outputDict['tmp_workflow'] == '10x':
    agg_samp = tmp_samp.groupby('Index (10x)',as_index=False).agg(lambda x: '_'.join(x.unique()))
    sum_samp = tmp_samp.groupby('Index (10x)',as_index=False).sum()
    max_samp = tmp_samp.groupby('Index (10x)',as_index=False).max()

    poo_samp = pd.DataFrame(columns=['Index (10x)', 'Lane', 'Loaded Cells'])

    if '*' in max_samp['Lane'].tolist():
        poo_samp['Index (10x)'] = agg_samp['Index (10x)']
    else:
        poo_samp['Index (10x)'] = agg_samp['Index (10x)'] + '_L' + max_samp['Lane'].astype(str)

    poo_samp['Lane'] = max_samp['Lane']
    poo_samp['Loaded Cells'] = sum_samp['Loaded Cells (Sample)']
    poo_samp['BCL PIN (10x)'] = agg_samp['BCL PIN (10x)']

    poo_samp.to_csv(f'{path_tmp}/poolTable.csv', index=False)

elif outputDict['tmp_workflow'] == '10x + HTO':
    agg_samp = tmp_samp.groupby('Index (10x)',as_index=False).agg(lambda x: '_'.join(x.unique()))
    sum_samp = tmp_samp.groupby('Index (10x)',as_index=False).sum()
    max_samp = tmp_samp.groupby('Index (10x)',as_index=False).max()

    poo_samp = pd.DataFrame(columns=['Index (10x)', 'Index (HTO)', 'Lane', 'Loaded Cells'])

    if '*' in max_samp['Lane'].tolist():
        poo_samp['Index (10x)'] = agg_samp['Index (10x)']
        poo_samp['Index (HTO)'] = agg_samp['Index (HTO)']
    else:
        poo_samp['Index (10x)'] = agg_samp['Index (10x)'] + '_L' + max_samp['Lane'].astype(str)
        poo_samp['Index (HTO)'] = agg_samp['Index (HTO)'] + '_L' + max_samp['Lane'].astype(str)

    poo_samp['Lane'] = max_samp['Lane']
    poo_samp['Loaded Cells'] = sum_samp['Loaded Cells (Sample)']
    poo_samp['BCL PIN (10x)'] = agg_samp['BCL PIN (10x)']
    poo_samp['BCL PIN (HTO)'] = agg_samp['BCL PIN (HTO)']

    poo_samp.to_csv(f'{path_tmp}/poolTable.csv', index=False)
logging.info('##')

logging.info('##\t10x libraries:')
logging.info('##\t-\tinitiating PIN dictionary...')
pinDict_10x = dict()
aggregate_pins_10x_df = tmp_samp.groupby('BCL PIN (10x)',as_index=False).agg(lambda x: '_'.join(x.unique()))

for index, row in aggregate_pins_10x_df.iterrows():
    pins = row['BCL PIN (10x)'].split(',')
    
    for pin in pins:
        pinDict_10x[pin] = (set(),'')

logging.info('##\t-\tupdating PIN dictionary with info (indices and full seq name)...')
for index, row in aggregate_pins_10x_df.iterrows():
    pins = row['BCL PIN (10x)'].split(',')
    index = row['Index (10x)'].split('_')

    for pin in pins:
        an_iterator = filter(lambda x:f'_{pin}' in x, files_bcls)
        tmp_seqname = list(an_iterator)[0]

        pinDict_10x[pin] = (pinDict_10x[pin][0].union(index),tmp_seqname)

logging.info('##\t-\tdictionary complete:')
for pin, (indices,seq_name) in pinDict_10x.items():
    logging.info(f'##\t\tPIN:\t\t{pin}')
    logging.info(f'##\t\tSeq name:\t{seq_name}')
    logging.info(f'##\t\tIndices:\t{indices}')

logging.info('##\t-\tcreating samplesheet...')
for pin, (indices,seq_name) in pinDict_10x.items():

    flowID = seq_name.split('_')[3][1:]

    with open(f'{path_bcls}{seq_name}/RunParameters.xml', 'r') as f:
        tmp_xml = f.read()
        tmp_bea = BeautifulSoup(tmp_xml, "xml")
        
        tmp_R1 = tmp_bea.find('Read1NumberOfCycles').text
        tmp_R2 = tmp_bea.find('Read2NumberOfCycles').text
        tmp_I1 = tmp_bea.find('IndexRead1NumberOfCycles').text
        tmp_I2 = tmp_bea.find('IndexRead2NumberOfCycles').text
            
        if int(tmp_I2) == 0:
            versionRNA = 'RNA v3.0'
            prefix = 'SI-GA-'
        elif int(tmp_I2) > 0:
            versionRNA = 'RNA v3.1'
            prefix = 'SI-TT-'

        outputDict['tmp_versionRNA'] = versionRNA
   
    if versionRNA == 'RNA v3.0':
        os.system(f'cp {path_templates}/template_samplesheet/samplesheet_singleIndex.csv {path_sample_sheets}/samplesheet_{pin}.csv') 
    if versionRNA == 'RNA v3.1':
        os.system(f'cp {path_templates}/template_samplesheet/samplesheet_dualIndex.csv {path_sample_sheets}/samplesheet_{pin}.csv')

    for index in indices:
        for ind in index.split(','):
            # 3.0
            if versionRNA == 'RNA v3.0':
                ind_code_1 = translate_index_10x_singleIndex[ind][0]
                ind_code_2 = translate_index_10x_singleIndex[ind][1]
                ind_code_3 = translate_index_10x_singleIndex[ind][2]
                ind_code_4 = translate_index_10x_singleIndex[ind][3]

                ind_1 = prefix + ind + '_1'
                ind_2 = prefix + ind + '_2'
                ind_3 = prefix + ind + '_3'
                ind_4 = prefix + ind + '_4'

                with open(f'{path_sample_sheets}/samplesheet_{pin}.csv', 'a') as w:
                    w.write(f',{ind_1},{ind_1},{ind_code_1},{flowID}\n')
                    w.write(f',{ind_2},{ind_2},{ind_code_2},{flowID}\n')
                    w.write(f',{ind_3},{ind_3},{ind_code_3},{flowID}\n')
                    w.write(f',{ind_4},{ind_4},{ind_code_4},{flowID}\n')

            # 3.1
            if versionRNA == 'RNA v3.1':
                ind_code_1 = translate_index_10x_dualIndex[ind][0]
                ind_code_2 = translate_index_10x_dualIndex[ind][1]
                ind = prefix + ind

                with open(f'{path_sample_sheets}/samplesheet_{pin}.csv', 'a') as w:
                    w.write(f',{ind},{ind},{ind_code_1},{ind_code_2},{flowID},{scopID}\n')

bcl_dict_path = open(f'{path_tmp}/pin-data-10x.pkl', 'wb')
pickle.dump(pinDict_10x, bcl_dict_path)
bcl_dict_path.close()
logging.info('##')

if outputDict['tmp_workflow'] == '10x + HTO':
    logging.info('##\tHTO libraries:')
    logging.info('##\t-\tinitiating PIN dictionary...')

    pinDict_hto = dict()
    aggregate_pins_hto_df = tmp_samp.groupby('BCL PIN (HTO)',as_index=False).agg(lambda x: '_'.join(x.unique()))
    
    for index, row in aggregate_pins_hto_df.iterrows():
        pins = row['BCL PIN (HTO)'].split(',')

        for pin in pins:
            pinDict_hto[pin] = (set(),'')
    
    logging.info('##\t-\tupdating PIN dictionary with info (indices and full seq name)...')
    for index, row in aggregate_pins_hto_df.iterrows():
        pins = row['BCL PIN (HTO)'].split(',')
        index = row['Index (HTO)'].split('_')

        for pin in pins:
            an_iterator = filter(lambda x:f'_{pin}' in x, files_bcls)
            tmp_seqname = list(an_iterator)[0]

            pinDict_hto[pin] = (pinDict_hto[pin][0].union(index),tmp_seqname)
    
    logging.info('##\t-\tdictionary complete:')
    for pin, (indices,seq_name) in pinDict_hto.items():
        logging.info(f'##\t\tPIN:\t\t{pin}')
        logging.info(f'##\t\tSeq name:\t{seq_name}')
        logging.info(f'##\t\tIndices:\t{indices}')

    logging.info('##\t-\tupdating samplesheet...')
    for pin, (indices,seq_name) in pinDict_hto.items():
       
        flowID = seq_name.split('_')[3][1:]
####

        for index in indices:
            for ind in index.split(','):
                # 3.0
                if versionRNA == 'RNA v3.0':
                    ind_code = translate_index_hto[ind][0][0:8]

                    with open(f'{path_sample_sheets}/samplesheet_{pin}.csv', 'a') as w:
                        w.write(f',{ind},{ind},{ind_code},{flowID}\n')

                # 3.1
                if versionRNA == 'RNA v3.1':
                    ind_code_1 = translate_index_hto[ind][0]
                    ind_code_2 = translate_index_hto[ind][1]

                    with open(f'{path_sample_sheets}/samplesheet_{pin}.csv', 'a') as w:
                        w.write(f',{ind},{ind},{ind_code_1},{ind_code_2},{flowID},{scopID}\n')

    bcl_dict_path = open(f'{path_tmp}/pin-data-hto.pkl', 'wb')
    pickle.dump(pinDict_hto, bcl_dict_path)
    bcl_dict_path.close()

    logging.info(f'##\t-\tcreating features.tsv(s)...')
    for index, row in poo_samp.iterrows():

        tmpList = list()
        tmpPool_list = row['Index (HTO)'].split(',')
        
        if len(tmpPool_list) == 1:
            tmpPool = tmpPool_list[0]
            
        elif len(tmpPool_list) > 1:
            tmpPool = '-'.join(tmpPool_list)

        tmpLane = row['Lane']

        tmpDf = tmp_samp[np.logical_and(tmp_samp['Index (HTO)'] == row['Index (HTO)'],tmp_samp['Lane'] == row['Lane'])]

        for jndex, jow in tmpDf.iterrows():
            tmpList.append([jow['Name'],translate_bc_hto[jow['HTO']]])

        features = pd.DataFrame(tmpList)

        dname = f'{path_hto_features}/{tmpPool}-features.tsv'

        features.to_csv(f'{dname}', index=False, sep = '\t', header = False)
logging.info('##')

# Output variables to tmpFolder
logging.info(f'##\tSaving temporary variables to temporary folder...')
for key, val in outputDict.items():
    if key not in ('tmp_samp'):
        tmpFile = open(f'{path_tmp}/{key}.txt','w')
        tmpFile.write(str(val))
        tmpFile.close()
    
os.system('touch {}'.format(snakemake.output[0]))

logging.info('##')
logging.info('###')
logging.info('#'*80)
