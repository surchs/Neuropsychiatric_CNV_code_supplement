import sys
import time
import warnings
import scipy as sp
import numpy as np
import patsy as pat
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from .tools import conn2mat
from sklearn import linear_model as sln
from sklearn import preprocessing as skp
from sklearn import model_selection as skm
from statsmodels.sandbox.stats.multicomp import multipletests as stm


def make_weights(data, pattern):
    if not data.shape[1:] == pattern.shape:
        raise Exception(f'data and pattern shape mismatch: {data.shape[1:]}, {pattern.shape}')
    n_nodes = pattern.shape[0]
    n_data = data.shape[0]
    weights = np.array([[np.corrcoef(data[data_id, node_id, :],
                                     pattern[node_id, :])[0, 1]
                         for node_id in range(n_nodes)]
                        for data_id in range(n_data)])
    return weights


def categorical_enrichment(pheno, data, group, case=None, control=None, node_labels='node_name_not_specificed'):
    sub_mask, cases_mask = find_subset(pheno, group, [case, control])
    sub_data = data[sub_mask, ...]
    n_data = sub_data.shape[1]
    n_case = np.sum(cases_mask[case])
    n_control = np.sum(cases_mask[control])

    results = {key: list() for key in ['node', 'U', 'p_value',
                                       f'median_{group}_{case}',
                                       f'median_{group}_{control}',
                                       'rank_biserial_correlation']}
    # Conduct Mann-Whitney-U for each node region
    for node_id in range(n_data):
        u_right, p = sp.stats.mannwhitneyu(sub_data[cases_mask[case], node_id], sub_data[cases_mask[control], node_id],
                                           alternative='two-sided')
        u_left = n_case * n_control - u_right
        u_min = np.min([u_left, u_right])
        # Compute the median for the case and control groups
        median_case = np.median(sub_data[cases_mask[case], node_id])
        median_control = np.median(sub_data[cases_mask[control], node_id])
        # Compute rank biserial correlation
        r_b = 1 - (2 * u_min) / (n_case * n_control)
        # Determine if cases > controls or reverse
        case_gt_con = u_right > u_min
        if not case_gt_con:
            r_b = -r_b
        # Store the results
        results['node'].append(node_id + 1)
        results['U'].append(u_min)
        results['p_value'].append(p)
        results[f'median_{group}_{case}'].append(median_case)
        results[f'median_{group}_{control}'].append(median_control)
        results['rank_biserial_correlation'].append(r_b)

    results_table = pd.DataFrame(data=results)
    # Correct for multiple comparisons
    (fdr_pass, qval, _, _) = stm(results_table.p_value, alpha=0.05, method='fdr_bh')
    results_table['q_value'] = qval
    # Add the node_names
    results_table['node_name'] = node_labels

    return results_table


def continuous_enrichment(pheno, data, covariate, node_labels='node_name_not_specificed'):
    sub_mask = find_subset(pheno, covariate)
    sub_pheno = pheno.loc[sub_mask]
    sub_data = data[sub_mask, :]
    n_data = sub_data.shape[1]

    results = {key: list() for key in ['node', 'pearson_r', 'p_value']}
    for node_id in range(n_data):
        r, p = stats.pearsonr(sub_pheno[covariate], sub_data[:, node_id])
        # Store the results
        results['node'].append(node_id + 1)
        results['pearson_r'].append(r)
        results['p_value'].append(p)

    results_table = pd.DataFrame(data=results)
    # Correct for multiple comparisons
    (fdr_pass, qval, _, _) = stm(results_table.p_value, alpha=0.05, method='fdr_bh')
    results_table['q_value'] = qval
    # Add the node_names
    results_table['node_name'] = node_labels

    return results_table


def find_subset(pheno, column, cases=None):
    # TODO check pandas type input
    subset_mask = np.array(~pheno[column].isnull())
    if cases is not None and not not cases:
        all_cases = pheno.loc[subset_mask][column].unique()
        try:
            case_available = np.array([True if case in all_cases else False for case in cases])
        except TypeError as e:
            raise Exception(f'the attribute "cases" needs to be iterable but is: {type(cases)}') from e
        if not all(case_available):
            if not any(case_available):
                raise Exception(f'none of the requested cases of "{column}" are available')
            else:
                warnings.warn(
                    f'\nnot all requested cases of "{column}" are available: {list(zip(cases, case_available))}',
                    RuntimeWarning)
        case_masks = np.array([pheno[column] == case for case in cases])
        subset_mask = np.any(case_masks, 0)
        # Return the masked instances of the requested cases
        cases_dict = {case: case_masks[idx][subset_mask] for idx, case in enumerate(cases)}
        return subset_mask, cases_dict
    else:
        return subset_mask


def standardize(data, mask):
    scaler = skp.StandardScaler(with_mean=False, with_std=True)
    scaler.fit(data[mask, :])
    standardized_data = scaler.transform(data)
    return standardized_data


def find_contrast(design_matrix, contrast):
    # Find the contrast column
    contrast_columns = [(col_id, col) for col_id, col in enumerate(design_matrix.columns) if f'{contrast}' in col]
    if not len(contrast_columns) == 1:
        raise Exception(f'There is no single factor that matches {contrast}: {(list(design_matrix.columns))}')
    return contrast_columns


def fast_glm(data, design_matrix, contrast):
    # Does not compute p-values but operates in parallel
    contrast_id, contrast_name = find_contrast(design_matrix, contrast)[0]
    glm = sln.LinearRegression(fit_intercept=False, normalize=False, n_jobs=-2)
    res = glm.fit(design_matrix, data)
    betas = res.coef_[:, contrast_id]
    return betas


def glm(data, design_matrix, contrast):
    contrast_id, contrast_name = find_contrast(design_matrix, contrast)[0]
    n_data = data.shape[1]

    # Conduct the GLM
    betas = np.zeros(shape=n_data)
    pvals = np.zeros(shape=n_data)
    for conn_id in range(n_data):
        model = sm.OLS(data[:, conn_id], design_matrix)
        results = model.fit()
        betas[conn_id] = results.params[contrast_id]
        pvals[conn_id] = results.pvalues[contrast_id]

    return betas, pvals


def glm_wrap_cc(conn, pheno, group, case, control, regressors='', report=False, fast=False):
    # Make sure pheno and conn have the same number of cases
    if not conn.shape[0] == pheno.shape[0]:
        print(f'Conn ({conn.shape[0]}) and pheno ({pheno.shape[0]}) must be same number of cases')

    # Define the subset of the sample
    sub_mask, case_masks = find_subset(pheno, group, [case, control])
    sub_conn = conn[sub_mask, :]
    sub_pheno = pheno.loc[sub_mask]
    n_sub = np.sum(sub_mask)
    n_case = np.sum(case_masks[case])
    n_control = np.sum(case_masks[control])
    n_data = sub_conn.shape[1]
    if report:
        print(f'Selected sample based on group variable {group}.\n'
              f'cases: {case} (n={n_case})\n'
              f'controls: {control} (n={n_control})\n'
              f'original sample: n={pheno.shape[0]}; new sample: n={n_sub}\n'
              f'{n_data} data points available\n'
              f'standardized estimators are based on {group}=={control}')

    stand_conn = standardize(sub_conn, case_masks[control])

    # Construct design matrix
    if type(control) == str:
        contrast = f'C({group}, Treatment("{control}"))'
    else:
        contrast = f'C({group}, Treatment({control}))'
    formula = ' + '.join((regressors, contrast))
    dmat = pat.dmatrix(formula, sub_pheno, return_type='dataframe')

    if fast:
        betas = fast_glm(sub_conn, dmat, group)
        table = pd.DataFrame(data={'betas': betas})
    else:
        betas, pvals = glm(sub_conn, dmat, group)
        stand_betas, _ = glm(stand_conn, dmat, group)
        table = pd.DataFrame(data={'betas': betas, 'stand_betas': stand_betas, 'pvals': pvals})

    return table


def glm_wrap_continuous(conn, pheno, contrast, regressors, report=False, fast=False):
    # Make sure pheno and conn have the same number of cases
    if not conn.shape[0] == pheno.shape[0]:
        print(f'Conn ({conn.shape[0]}) and pheno ({pheno.shape[0]}) must be same number of cases')

    # Define the subset of the sample
    sub_mask = find_subset(pheno, contrast)
    sub_conn = conn[sub_mask, :]
    sub_pheno = pheno.loc[sub_mask]
    n_sub = sub_pheno.shape[0]
    n_data = sub_conn.shape[1]
    sub_conn_stand = standardize(sub_conn, np.ones(n_sub).astype(bool))

    if report:
        print(f'Selected sample based on contrast variable {contrast}.\n'
              f'Found {n_sub} subjects with no missing data for {contrast}\n'
              f'original sample: n={pheno.shape[0]}; new sample: n={n_sub}\n'
              f'{n_data} data points available\n'
              f'standardized estimators are based on all subjects with no missing data for {contrast}')

    formula = ' + '.join((regressors, contrast))
    design_matrix = pat.dmatrix(formula, sub_pheno, return_type='dataframe')

    if fast:
        betas = fast_glm(sub_conn, design_matrix, contrast)
        table = pd.DataFrame(data={'betas': betas})
    else:
        betas, pvals = glm(sub_conn, design_matrix, contrast)
        stand_betas, _ = glm(sub_conn_stand, design_matrix, contrast)
        table = pd.DataFrame(data={'betas': betas, 'stand_betas': stand_betas, 'pvals': pvals})

    return table


def permutation_glm_contrast(conn, pheno, n_a, n_a_1, n_b, n_b_1, n_iter=5000):
    # conn: connectome (n_sub, n_conn)
    # pheno: pheno pandas table, same order as conn
    # n_a / n_b: number of subs per group a / b
    # n_a_1 / n_b_1: number of subs in the first contrast of group a / b
    # n_iter: number of permutations

    # Some sanity checks because people can be stupid
    if n_a_1 > n_a:
        raise ValueError(f'n_a_1 ({n_a_1}) cannot be larger than n_a ({n_a})')
    if n_b_1 > n_b:
        raise ValueError(f'n_b_1 ({n_b_1}) cannot be larger than n_b ({n_b})')

    shsp = skm.ShuffleSplit(n_splits=n_iter, test_size=0.5)
    pheno.reset_index(inplace=True, drop=True)
    pheno.loc[:, 'group'] = 0
    n_subjects = pheno.shape[0]
    beta = np.zeros((conn.shape[1], 2, n_iter))

    start = time.time()
    for fold_id, (group_a, group_b) in enumerate(shsp.split(np.ones(n_subjects))):
        pheno_a = pheno.loc[group_a[:n_a]].copy()
        pheno_b = pheno.loc[group_b[:n_b]].copy()
        # Assign the groups
        pheno_a.loc[group_a[:n_a_1], 'group'] = 1
        pheno_b.loc[group_b[:n_b_1], 'group'] = 1
        mat1 = pat.dmatrix('sex_dummy + Age_Bin + FD_scrubbed_s1_s2 + group', pheno_a, return_type='dataframe')
        mat2 = pat.dmatrix('sex_dummy + Age_Bin + FD_scrubbed_s1_s2 + group', pheno_b, return_type='dataframe')
        beta[:, 0, fold_id] = fast_glm(conn[group_a[:n_a], :], mat1, 'group')
        beta[:, 1, fold_id] = fast_glm(conn[group_b[:n_b], :], mat2, 'group')

        elapsed = time.time() - start
        done = fold_id + 1
        remaining = n_iter - done
        time_left = (elapsed / done) * remaining

        sys.stdout.write('\r {}/{}. {:.2f}s left ({:.3f}s)'.format(done, n_iter, time_left, elapsed / done))
        sys.stdout.flush()
    sys.stdout.write('\r Done. Took {:.2f}s'.format(elapsed))
    sys.stdout.flush()

    return beta


def contrast_correlation(beta, mask, mode='region'):
    # mode: 'region' or 'global'
    if mode == 'region':
        n_iter = beta.shape[2]
        n_seed = mask.shape[0]
        corr = np.zeros((n_seed, n_iter))
        for fold_id in range(n_iter):
            # Map the betas back to a matrix
            beta1_mat = conn2mat(beta[:, 0, fold_id], mask)
            beta2_mat = conn2mat(beta[:, 1, fold_id], mask)
            # Compute the correlation
            corr[:, fold_id] = np.array([np.corrcoef(beta1_mat[nid, :], beta2_mat[nid, :])[0, 1] for nid in range(64)])
        return corr
    elif mode == 'global':
        n_iter = beta.shape[2]
        corr = np.array([np.corrcoef(beta[:, 0, fold_id], beta[:, 1, fold_id])[0, 1] for fold_id in range(n_iter)])
        return corr
    else:
        raise ValueError(f'mode {mode} is not configured!')


def permutation_glm(pheno, conn, group, case, control, regressors='', n_iter=1000, stand=False):
    # Make sure pheno and conn have the same number of cases
    if not conn.shape[0] == pheno.shape[0]:
        print(f'Conn ({conn.shape[0]}) and pheno ({pheno.shape[0]}) must be same number of cases')
    sub_mask, case_masks = find_subset(pheno, group, [case, control])
    sub_pheno = pheno.loc[sub_mask, :]
    sub_pheno.loc[:, 'group'] = ''
    sub_conn = conn[sub_mask, :]
    # Standardize the input data
    if stand:
        sub_conn = standardize(sub_conn, case_masks[control])
    n_conn = sub_conn.shape[1]

    n1 = np.sum(case_masks[case])
    n2 = np.sum(case_masks[control])
    r2 = n2 / (n1 + n2)

    splitter = skm.ShuffleSplit(n_splits=n_iter, test_size=r2, random_state=1)

    i_pheno = sub_pheno.copy()
    group_id = list(i_pheno.columns).index('group')

    betas = np.zeros(shape=(n_iter, n_conn))
    start = time.time()
    for fold_id, (train_id, test_id) in enumerate(splitter.split(sub_conn)):
        n1_ids = train_id[:n1]
        n2_ids = test_id[:n2]
        i_pheno = sub_pheno.copy()
        i_pheno['group'] = ''
        i_pheno.iloc[n1_ids, group_id] = 'case'
        i_pheno.iloc[n2_ids, group_id] = 'control'
        table = glm_wrap_cc(sub_conn, i_pheno, 'group', 'case', 'control', regressors=regressors,
                            report=False, fast=True)
        betas[fold_id, :] = table['betas'].values

        elapsed = time.time() - start
        done = fold_id + 1
        remaining = n_iter - done
        time_left = (elapsed / done) * remaining

        sys.stdout.write('\r {}/{}. {:.2f}s left ({:.3f}s per permutation)'.format(done, n_iter, time_left, elapsed / done))
        sys.stdout.flush()
    sys.stdout.write('\r Done. Took {:.2f}s'.format(elapsed))
    sys.stdout.flush()

    return betas
