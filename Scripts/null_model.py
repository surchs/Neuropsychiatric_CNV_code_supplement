import sys
sys.path.append('../')
import cnvfc
import numpy as np
import pandas as pd
import pathlib as pal

n_iter = 5000
root_p = pal.Path('../data/')
pheno_p = root_p / 'pheno/phenotypic_information.csv'
connectome_p = root_p / 'preprocessed/connectome/sample_connectome/python/'
connectome_t = 'connectome_s{}_cambridge64.npy'
out_p = root_p / 'processed/null_model/'
if not out_p.is_dir():
    out_p.mkdir()

case = 'Diagnosed'
control = 'Control'
group = 'DX_GROUP'
regressors_str = ' + '.join(['SITE_ID', 'FD_scrubbed', 'AGE_AT_SCAN'])

conn_mask = np.tril(np.ones((64, 64))).astype(bool)
pheno = pd.read_csv(pheno_p)

paths = [(connectome_p / connectome_t.format(row.Subject)).resolve() for rid, row in pheno.iterrows()]
conn_stack = np.array([np.load(p)[conn_mask] for p in paths])

# Run Null Model for deletion cases
permuted_betas = cnvfc.stats.permutation_glm(pheno, conn_stack, group, case, control,
                                             regressors=regressors_str, n_iter=5000, stand=False)
# Save the betas
np.save(out_p / f'icc_sample_null_model_{group}_case_vs_control.npy', permuted_betas)
