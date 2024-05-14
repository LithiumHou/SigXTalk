import pandas as pd
import numpy as np
import itertools
import dhg
import multiprocessing as mp

from sklearn.metrics import mutual_info_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import kendalltau, entropy, spearmanr




# filter the databases
def filter_DB(RecTFDB, TFTGDB, gene_all, recs = None):

    # RecTFDB = RecTFDB[RecTFDB.apply(lambda row: all(value in gene_all for value in row), axis=1)]
    # TFTGDB = TFTGDB[TFTGDB.apply(lambda row: all(value in gene_all for value in row), axis=1)]
    RecTFDB.columns = ['From','To']
    TFTGDB.columns = ['From','To']

    if recs is None:
        recs = list(set(RecTFDB['From']) - set(RecTFDB['To']))

    recs = set(recs).intersection(set(gene_all))
    recs = recs.intersection(set(RecTFDB['From']))
    recs = list(recs)
    tfs_1 = list(set(TFTGDB['From']))
    tfs_2 = list(set(RecTFDB['To']))
    tfs_ori = set(tfs_1).intersection(set(tfs_2))
    tfs_ori = list(set(tfs_ori).intersection(set(gene_all)))
    tgs_ori = list(set(TFTGDB['To']).intersection(set(gene_all)))

    tgs_new = list(set(tgs_ori) - set(recs+tfs_1))
    tfs_new = list(set(tfs_ori) - set(tgs_new+recs))
    recs_new = recs

    RecTFDB = RecTFDB[RecTFDB['From'].isin(recs_new)]
    RecTFDB = RecTFDB[RecTFDB['To'].isin(tfs_new)]
    TFTGDB = TFTGDB[TFTGDB['From'].isin(tfs_new)]
    TFTGDB = TFTGDB[TFTGDB['To'].isin(tgs_new)]

    return RecTFDB,TFTGDB

# Do the standardscaler pipeline on a dataframe
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


# calculate the correlation between variables
def calculate_corr_single(data, pair, type = 'MI'):
    # data: a pandas df, cols = genes(vars), rows = cells(samples)
    # type: MI for mutual information, Corr for Kendall's correlation
    variable1 = pair[0]
    variable2 = pair[1]
    if type == 'Spearman':
        corr, pv = spearmanr(data[variable1], data[variable2])
        return corr
    elif type == 'MI':
        MI = mutual_info_score(data[variable1], data[variable2])
        entropy_max = np.max([entropy(data[variable1]),entropy(data[variable2])])
        return MI/entropy_max
    else:
        corr, pv = kendalltau(data[variable1], data[variable2])
        if pv >= 0.2:
            corr = 0
        return corr


def calculate_corr(Exp_mat, DB, type = 'Spearman'):

    corr_df = DB.copy()
    corr_df.columns = ['From','To']

    cor_tuples = list(corr_df.to_records(index=False))

    corr_args = [(Exp_mat, pair, type) for pair in cor_tuples]

    with mp.Pool(processes=int(mp.cpu_count()/2)) as pool:
        # Use starmap to pass multiple arguments to the compute function
        results = pool.starmap(calculate_corr_single, corr_args)
    pool.close()
    pool.join()
    
    corr_df['Correlation'] = results
    return corr_df
    
def construct_graph(cor_mat,g_thres):
    
    cor_mat = cor_mat - np.eye(len(cor_mat))

    top_percent_threshold = np.percentile(cor_mat.values, (1-g_thres)*100)
    filtered_mat = (cor_mat > top_percent_threshold).astype(int)
    filtered_mat = filtered_mat.values
    # Convert the adj matrix to edge lists
    rows, cols = np.triu_indices_from(filtered_mat, k=1)
    edges_indices = filtered_mat[rows, cols] > 0
    edges = list(zip(rows[edges_indices], cols[edges_indices]))
    
    genes = cor_mat.columns
    G = dhg.Graph(len(genes),edges)
    
    return G


# construct the hypergraph using the correlation/mi matrix
def construct_hypergraph(df1, df2, hg_thres):

    df1.columns = ['Receptor','TF','Weight1']
    df2.columns = ['TF','TG','Weight2']

    df_all = pd.merge(df1,df2,on='TF')
    df_all['Weight'] = df_all['Weight1']*df_all['Weight2']
    df_all['Weight'] = df_all['Weight'].abs()
    del df_all['Weight1']
    del df_all['Weight2']

    df_sorted = df_all.sort_values(by='Weight',ascending=False)
    df_sorted.reindex(columns=['Receptor','TF','TG','Weight'])

    all_len = len(df_sorted)
    hg_index = int(hg_thres*all_len)
    df_hg = df_sorted.iloc[:hg_index]

    df_hg_filtered = df_hg[['Receptor','TF','TG']]

    return df_hg_filtered, df_sorted

# generate positive training samples for the model
def Generate_positive(df,sample_thres):

    df_sorted = df.sort_values(by='Weight',ascending=False)
    df_sorted.reindex(columns=['Receptor','TF','TG','Weight'])

    all_len = len(df_sorted)
    sample_index = int(sample_thres*all_len)    
    df_samples = df_sorted.iloc[:sample_index]
    df_samples['label'] = 1
    del df_samples['Weight']

    return df_samples

def Generate_negative_2(df):
    # Generate all possible combinations of X, Y values

    df2 = df.copy()
    df2.columns = ['From','To','label']
    N1 = set(df2.iloc[:,0])
    N2 = set(df2.iloc[:,1])

    all_combinations = list(itertools.product(N1, N2))
    df_all = pd.DataFrame(all_combinations, columns=['From','To'])

    # Filter out combinations that are not already in the DataFram
    merged_df = df_all.merge(df2, on=['From','To'], how='left', indicator=True)
    neg_samples = merged_df[merged_df['_merge'] == 'left_only']
    neg_samples = neg_samples.drop(columns='_merge')
    neg_samples['label'] = 0
    
    del all_combinations, merged_df
    return neg_samples

def Generate_negative_3(df):
    # df: all positive samples
    # Generate all possible combinations of X, Y, and Z values
    N1 = set(df['Receptor'])
    N2 = set(df['TF'])
    N3 = set(df['TG'])

    all_combinations = list(itertools.product(N1, N2, N3))
    df_all = pd.DataFrame(all_combinations, columns=['Receptor','TF','TG'])
    df = df[['Receptor','TF','TG']]

    # Filter out combinations that are not already in the DataFram
    merged_df = df_all.merge(df, on=['Receptor','TF','TG'], how='left', indicator=True)
    neg_samples = merged_df[merged_df['_merge'] == 'left_only']
    neg_samples = neg_samples.drop(columns='_merge')
    neg_samples_filtered = neg_samples.sample(n = len(df))  
    neg_samples_filtered['label'] = 0

    return neg_samples_filtered

def Generate_tail(df,sample_thres):

    df_sorted = df.sort_values(by='Weight',ascending=True)
    df_sorted.reindex(columns=['Receptor','TF','TG','Weight'])

    all_len = len(df_sorted)
    sample_index = int(sample_thres*all_len)    
    df_samples = df_sorted.iloc[:sample_index]
    df_samples['label'] = 0
    del df_samples['Weight']

    return df_samples
