import pandas as pd
import numpy as np


def cov_matrix(constellation, labels, cov_dfs, key='mean gap dark'):
    
    ids = constellation.sats_id
    
    df = pd.DataFrame(index=ids, columns=labels)
    
    for i in range(len(ids)):
        for j in range(len(labels)):
            tmp_mean_total = np.mean(cov_dfs[j].collapse_sats().to_df()[key].to_numpy()/3600)
            tmp_mean_partial = np.mean(cov_dfs[j].collapse_sats([ids[k] for k in range(len(ids)) if k != i ]).to_df()[key].to_numpy()/3600)
            df.iloc[i][j] = tmp_mean_partial - tmp_mean_total
    
    return df


