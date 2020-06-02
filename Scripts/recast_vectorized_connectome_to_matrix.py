import sys
sys.path.append('../')
import cnvfc
import numpy as np
import pandas as pd
import pathlib as pal

root_p = pal.Path('../')
pheno_p = root_p / 'data/pheno/phenotypic_information.csv'
connectome_p = root_p / 'data/preprocessed/connectome/sample_connectome/niak/'
connectome_t = 'connectome_s{}_cambridge64.mat'
connectome_n = 'sample_con'
out_p = root_p / 'data/preprocessed/connectome/sample_connectome/python/'
if not out_p.is_dir():
    out_p.mkdir()

conn_mask = np.tril(np.ones((64,64))).astype(bool)
pheno = pd.read_csv(pheno_p)

# Load the pre-computed individual vectorized connectomes and recast them into a seed by seed FC matrix
for rid, row in pheno.iterrows():
    sub = row.SUB_ID
    p = connectome_p / connectome_t.format(sub)
    new_name = connectome_t.format(sub).split('.')[0]
    # Check if output already exists
    if not p.is_file():
        print(f'{sub} does not have {p.resolve()}')
    conn_mat = cnvfc.tools.niak_conn2mat(p.resolve(), connectome_n, conn_mask)
    
    np.save((out_p / f'{new_name}.npy').resolve(), conn_mat)
    if rid%50==0:
        print(rid)
