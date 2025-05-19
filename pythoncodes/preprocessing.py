import pandas as pd
import numpy as np
import itertools
import dhg
import multiprocessing as mp
import random

from sklearn.metrics import mutual_info_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import kendalltau, entropy, spearmanr

def calculate_corr_single(data, pair, type = 'Spearman',abss = False):
    # calculate the correlation between variables
    # data: a pandas df, cols = genes(vars), rows = cells(samples)
    # type: MI for mutual information, Corr for Kendall's correlation
    variable1 = pair[0]
    variable2 = pair[1]
    type = type.lower()
    if type == 'spearman':
        corr, pv = spearmanr(data[variable1], data[variable2])
    elif type == 'mi':
        MI = mutual_info_score(data[variable1], data[variable2])
        entropy_max = np.max([entropy(data[variable1]),entropy(data[variable2])])
        corr = MI/entropy_max
    elif type == 'pearson':
        return np.corrcoef(data[variable1], data[variable2])[0,1]
    else:
        corr, pv = kendalltau(data[variable1], data[variable2])
        if pv >= 0.2:
            corr = 0
    if abss:
        return abs(corr)
    else:
        return max(corr, 0)


def calculate_corr(Exp_mat, DB, type = 'Spearman',abss = True):

    corr_df = DB.iloc[:,range(2)]
    corr_df.columns = ['From','To']

    cor_tuples = list(corr_df.to_records(index=False))

    corr_args = [(Exp_mat, pair, type, abss) for pair in cor_tuples]

    with mp.Pool(processes=int(mp.cpu_count()/2)) as pool:
        # Use starmap to pass multiple arguments to the compute function
        results = pool.starmap(calculate_corr_single, corr_args)
    pool.close()
    pool.join()
    
    corr_df['Correlation'] = results

    return corr_df


def Generate_negative_3(df, args, sample_n = 0):

    # Check the current random seed
    curr_seed = args.seed
    # df: all positive samples
    # Generate all possible combinations of X, Y, and Z values
    N1 = set(df['Receptor'])
    N2 = set(df['TF'])
    N3 = set(df['TG'])

    all_combinations = list(itertools.product(N1, N2, N3))
    df_all = pd.DataFrame(all_combinations, columns=['Receptor','TF','TG'])
    df = df[['Receptor','TF','TG']]
    random.seed(curr_seed)
    # Filter out combinations that are not already in the DataFram
    merged_df = df_all.merge(df, on=['Receptor','TF','TG'], how='left', indicator=True)
    neg_samples = merged_df[merged_df['_merge'] == 'left_only']
    neg_samples = neg_samples.drop(columns='_merge')
    if sample_n == 0:
        neg_samples_filtered = neg_samples.sample(n = len(df))
    else:
        neg_samples_filtered = neg_samples.sample(n = sample_n)

    return neg_samples_filtered


def Filter_RTTDB(df1, df2, thres1, thres2, genes, args, first = 'score', sample_scale = 0.5, ood_frac = 0.5):

    df1.columns = ['Receptor','TF','score1']
    df2.columns = ['TF','TG','score2']
    df1 = df1.sort_values(by = 'score1', ascending = False)
    df2 = df2.sort_values(by = 'score2', ascending = False)

    df_all = pd.merge(df1, df2, on = 'TF', how = 'inner')
    df_all = df_all[['Receptor','TF','TG','score1','score2']]
    df_all['score'] = df_all['score1']*df_all['score2']
    df_all = df_all.sort_values(by = 'score', ascending = False)

    if first == 'RecTF' or first == 'TFTG':
        n_selected = int(len(df1)*thres1)
        df1_topk = df1.head(n_selected)
        df1_bottomk = df1.tail(n_selected)
        n_selected = int(len(df2)*thres1)
        df2_topk = df2.head(n_selected)
        df2_bottomk = df2.tail(n_selected)
        df_pos = pd.merge(df1_topk, df2_topk, on = 'TF', how = 'inner')
        df_neg = pd.merge(df1_bottomk, df2_bottomk, on = 'TF', how = 'inner')

    else: # first == 'score'        
        n_selected = int(len(df_all)*thres1*thres2)
        df_pos = df_all.head(n_selected)
        df_neg = df_all.tail(n_selected)

    # Generate the hypergraph framework
    df_pos = df_pos[['Receptor','TF','TG']]
    df_neg = df_neg[['Receptor','TF','TG']]
    df_hg = df_pos.copy()

    # Generate the negative samples
    df_pos['label'] = 1
    df_neg1 = Generate_negative_3(df_all[['Receptor','TF','TG']], args, sample_n = round(len(df_pos)*ood_frac))
    df_neg1['label'] = 0
    df_neg2 = df_neg.sample(n = (len(df_pos)-len(df_neg1)))
    df_neg2['label'] = 0

    df_pos = df_pos.sample(n = round(len(df_pos)*sample_scale))
    df_neg1 = df_neg1.sample(n = round(len(df_neg1)*sample_scale))
    df_neg2 = df_neg2.sample(n = round(len(df_neg2)*sample_scale))
    df_samples = pd.concat([df_pos,df_neg1,df_neg2],axis = 0)

    # Map gene names to indexes
    string_to_index = {string: index for index, string in enumerate(genes)}
    def map_strings_to_indexes(value):
        return string_to_index.get(value, value)
    
    # df_all['Receptor'] = df_all['Receptor'].applymap(map_strings_to_indexes)
    # df_all['TF'] = df_all['TF'].applymap(map_strings_to_indexes)
    # df_all['TG'] = df_all['TG'].applymap(map_strings_to_indexes)
    df_hg = df_hg.applymap(map_strings_to_indexes)
    df_all[['Receptor','TF','TG']] = df_all[['Receptor','TF','TG']].applymap(map_strings_to_indexes)
    df_samples[['Receptor','TF','TG']] = df_samples[['Receptor','TF','TG']].applymap(map_strings_to_indexes)

    edge_list = [tuple(row) for row in df_hg.values]
    hg = dhg.Hypergraph(len(genes), edge_list)

    return hg, df_samples, df_all


def data_scale(df, pca = False,n_comp = 500):
    # Extract numerical columns
    data = df.values
    # Create StandardScaler instance
    scaler = StandardScaler()

    # Fit and transform on numerical columns
    scaled_data = scaler.fit_transform(data)

    if pca:
        nc = min(n_comp,int(0.9*min(data.shape)))
        pca = PCA(n_components=nc)
        pca.fit(scaled_data)
        df_embeddings = pca.transform(scaled_data)
        df_embeddings = pd.DataFrame(df_embeddings, columns=[f'PC{i+1}' for i in range(nc)],index = df.index)           
    else:
        df_embeddings = pd.DataFrame(scaled_data, columns=df.columns, index = df.index)

    return df_embeddings.astype(np.float32)

# Do the normalization based on cell counts
def NormalizeData(df, log2 = False): # rows: cells, cols: genes
    exp = df.values
    ff = np.median(np.sum(exp,axis=1))/np.sum(exp,axis = 1)
    exp_nor = np.dot(np.diag(ff),exp)
    if log2:
        exp_lognor = np.log2(1+exp_nor)
        exp_df = pd.DataFrame(exp_lognor, index = df.index, columns = df.columns)
    else:
        exp_df = pd.DataFrame(exp_nor, index = df.index, columns = df.columns)

    return exp_df
