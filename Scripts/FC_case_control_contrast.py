import sys
sys.path.append('../')
import cnvfc
import numpy as np
import pandas as pd
import pathlib as pal

root_p = pal.Path('../data/')
pheno_p = root_p / 'pheno/Pheno.csv'
connectome_p = root_p / 'preprocessed/connectome/sample_connectome/python/'
connectome_t = 'connectome_s{}_mist64.npy'
label_p = root_p / 'parcellation/Parcel_Information/MIST_64.csv'
out_p = root_p / 'processed/fc_profiles/'
if not out_p.is_dir():
    out_p.mkdir()

conn_mask = np.tril(np.ones((64, 64))).astype(bool)
pheno = pd.read_csv(pheno_p)
pheno.rename(columns={'Unnamed: 0': 'niak_id'}, inplace=True)
labels = pd.read_csv(label_p, sep=';')
roi_labels = labels.label.values

# Read in the paths to the pre-computed individual seed by seed FC matrices and put them in an array
paths = [(connectome_p / connectome_t.format(row.Subject)).resolve() for rid, row in pheno.iterrows()]
conn_stack = np.array([np.load(p)[conn_mask] for p in paths])

group = 'DX_GROUP'
regressors = ' + '.join(['SITE_ID', 'FD_scrubbed', 'AGE_AT_SCAN'])

# Define case-control contrast
glm = cnvfc.stats.glm_wrap_cc(conn_stack, pheno, group, case='Diagnosed', control='Control',
                                  regressors=regressors, report=True)
table, table_stand_beta, table_qval = cnvfc.tools.summarize_glm(glm, conn_mask, roi_labels)

# Store the results
table.to_csv(out_p / 'icc_case_vs_con.tsv', sep='\t')
table_stand_beta.to_csv(out_p / 'icc_case_vs_con_standardized_betas.tsv', sep='\t')
table_qval.to_csv(out_p / 'icc_case_vs_con_fdr_corrected_pvalues.tsv', sep='\t')
