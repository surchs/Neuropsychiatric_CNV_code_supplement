import numpy as np
import pandas as pd
from scipy import io as sio
from statsmodels.sandbox.stats.multicomp import multipletests as stm


def octvec2mat(vec, mask):
    """
    Takes a vector in Fortran style (Octave/Matlab)
    and maps it back into a mask correctly in
    python
    """
    vec_mat = mask.flatten(order='F').astype(bool)
    tmp = np.zeros_like(vec_mat, dtype=float)
    tmp[vec_mat] = vec
    vol = np.reshape(tmp, mask.shape, order='F')
    return vol


def oct2mat(vec, mask):
    temp = octvec2mat(vec, mask)
    temp += temp.T
    # Reset diagonal
    temp[np.eye(mask.shape[0]).astype(bool)] = temp[np.eye(mask.shape[0]).astype(bool)]/2
    return temp


def niak_conn2mat(connectome_p, contrast_name, mask):
    mat = sio.loadmat(connectome_p)
    conn = mat[contrast_name]['connectome'][0][0].squeeze()
    return oct2mat(conn, mask)


def conn2mat(conn, mask):
    conn_mat = np.zeros(shape=mask.shape)
    conn_mat[mask] = conn
    conn_mat += conn_mat.T
    conn_mat[np.eye(mask.shape[0]).astype(bool)] /= 2

    return conn_mat


def report_connectivity_alterations(table):
    n_conn = table.shape[0]
    mask_sig = table.qvals.values < 0.05
    if not any(mask_sig):
        return f'None of the {n_conn} connections are significant.'
    effects = table.loc[mask_sig].stand_betas.values

    n_eff = effects[effects < 0]
    p_eff = effects[effects > 0]
    n_pos = len(p_eff)
    n_neg = len(n_eff)
    n_hits = len(effects)
    report = f'There are {n_hits} significant connections out of {n_conn}.\n' \
        f'Of those, {n_pos} show an increase and {n_neg} show a decrease in connectivity.\n' \
        f'The range of effects is:\n' \
        f'Positive effects (z-scores): {np.round(np.min(p_eff), 3) if n_pos > 0 else "n.s."} to ' \
        f'{np.round(np.max(p_eff), 3) if n_pos > 0 else "n.s."}\n' \
        f'Negative effects (z-scores): {np.round(np.min(n_eff), 3) if n_neg > 0 else "n.s."} to ' \
        f'{np.round(np.max(n_eff), 3) if n_neg > 0 else "n.s."}'
    return report


def report_enrichment(table, corrected=True):
    # Determine the kind of enrichment
    if 'U' in table.columns:
        kind = 'categorical'
    elif 'pearson_r' in table.columns:
        kind = 'continuous'
    else:
        raise Exception('I am not sure if this is a continuous or categorical enrichment')
    n_region = table.shape[0]
    if corrected:
        test='q_value'
    if not corrected:
        test='p_value'
    if not any(table[test] < 0.05):
        report = f'{"Corrected" if corrected else "Uncorrected"} Statistics:\n' \
                 f'None of the {n_region} regions shows significant enrichment.'
        return report
    n_significant = np.sum(table[test]<0.05)
    rep_head = f'{"Corrected" if corrected else "Uncorrected"} Statistics:\n' \
               f'There are {n_significant} out of {n_region} regions that show significant enrichment.'
    if kind == 'categorical':
        rep_body = [f'FC in {row.node_name} is {"positively" if row.rank_biserial_correlation > 0 else "negatively"} ' \
                    f'enriched:\n    U: {row.U:.0f} {test}: {row[test]:.2E} ' \
                    f'rank_corr: {row.rank_biserial_correlation:.3f}.'
                    for rid, row in table.loc[table[test] < 0.05].iterrows()]
    else:
        rep_body = [f'FC in {row.node_name} is {"positively" if row.pearson_r > 0 else "negatively"} ' \
                        f'enriched:\n    {test}: {row[test]:.2E} ' \
                        f'pearson_r: {row.pearson_r:.3f}.'
                    for rid, row in table.loc[table[test] < 0.05].iterrows()]
    report = '\n'.join([rep_head] + rep_body)
    return report


def summarize_glm(glm_table, conn_mask, roi_labels):
    out_table = glm_table.copy()
    (fdr_pass, qval, _, _) = stm(glm_table.pvals, alpha=0.05, method='fdr_bh')
    out_table['qval'] = qval
    # Return to matrix form
    stand_beta_table = pd.DataFrame(conn2mat(out_table.stand_betas.values, conn_mask) , index=roi_labels, columns=roi_labels)
    qval_table = pd.DataFrame(conn2mat(out_table.pvals.values, conn_mask), index=roi_labels, columns=roi_labels)
    return out_table, stand_beta_table, qval_table
