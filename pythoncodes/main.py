import torch
import torch.nn.functional as F
import random
import numpy as np
import pandas as pd
import sys
import argparse
import dhg
import xgboost as xgb
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
# from sklearn.decomposition import PCA
# from sklearn.linear_model import LinearRegression

work_dir = '/pythoncodes'
sys.path.append(work_dir)

from preprocessing import *
from predictor import *
from training import *

# Parse the args for hyperlink prediction

parser = argparse.ArgumentParser(description="Process some inputs.")

parser.add_argument('--TG', type=int, default=1, help='Whether consider certain target genes only')
parser.add_argument('--thres', type=float, default=[0.1,0.05,0.1,0.75], help='The thresholds for expression, training set, hypergraph construction and XGBoost')
parser.add_argument('--corr_type', type=str, default='Kendall', help='The method to calculate correlations')
parser.add_argument('--train_size', type=float, default=0.8, help='Fraction of training samples.')
parser.add_argument('--hgnn_dims', type=int, default=[256,128], help='The dimension of hidden layer')
parser.add_argument('--linear_dims', type=int, default=[64, 32], help='The dimensions of the MLP layer')
parser.add_argument('--lr', type=float, default= 0.01, help='Initial learning rate.')
parser.add_argument('--epochs', type=int, default= 50, help='Number of epoch.')
parser.add_argument('--batch_size', type=int, default=256, help='The size of each batch')
parser.add_argument('--model_output', type=bool, default=False, help='Whether output the auroc and aupr values')
parser.add_argument('--seed', type=int, default=2024, help='Random seed')

args = parser.parse_args()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

random.seed = args.seed
torch.manual_seed(args.seed)
np.random.seed(args.seed)

print(f"Dataset: {args.project}, target type: {args.target_type}")


# Load and filter the input expression matrix
input_exp = pd.read_csv('/inputs/ExpressionCount.csv', header = 0, index_col=0, sep=' ') # row: genes, col: cells
input_exp = input_exp.astype(np.float32)
exp_thres = args.thres[0]
gene_all = input_exp.index[input_exp.astype(bool).sum(axis=1) > exp_thres*len(input_exp.columns)].tolist()
input_exp = input_exp.loc[gene_all]
input_exp_nor = data_scale(input_exp,pca=True,n_comp=250)

input_used = input_exp_nor.copy()
input_tensor = torch.tensor(input_used.values, dtype=torch.float32)
input_tensor = input_tensor.to(device)

# load the gene-gene interaction databases
LR = pd.read_csv('/inputs/LigRec.csv',header=0, sep = ' ')
Recs = list(set(LR['To']))
RecTFDB = pd.read_csv('/inputs/RecTFDB.txt',header=0, sep = ' ')
TFTGDB = pd.read_csv('/inputs/TFTGDB.txt',header=0, sep = ' ')
RecTFDB, TFTGDB = filter_DB(RecTFDB,TFTGDB,gene_all,recs=Recs)

# map gene names to indexes
string_to_index = {string: index for index, string in enumerate(gene_all)}
def map_strings_to_indexes(value):
    return string_to_index.get(value, value)

# construct the hypergraph and generate the training samples
df_rt = calculate_corr(input_exp.T, RecTFDB, type = args.corr_type)
df_tftg = calculate_corr(input_exp.T, TFTGDB, type = args.corr_type)
df_rt['Correlation'] = df_rt['Correlation'].abs()
df_tftg['Correlation'] = df_tftg['Correlation'].abs()

input_hg, input_all = construct_hypergraph(df_rt, df_tftg, args.thres[2])
input_hg = input_hg.applymap(map_strings_to_indexes)
input_all = input_all.applymap(map_strings_to_indexes)
num_genes = len(gene_all)
edge_list = [tuple(row) for row in input_hg.values]
hg = dhg.Hypergraph(num_genes, edge_list)
hg = hg.to(device)

pos_sample = Generate_positive(input_all, sample_thres=args.thres[1])
neg_sample1 = Generate_negative_3(input_all)
neg_sample1 = neg_sample1.sample(n = round(len(pos_sample)*0.5))
neg_sample2 = Generate_tail(input_all, sample_thres=args.thres[1])
neg_sample2 = neg_sample2.sample(n = len(pos_sample)-len(neg_sample1))
all_samples = pd.concat([pos_sample,neg_sample1,neg_sample2],axis = 0)
all_samples[['Receptor','TF','TG']] = all_samples[['Receptor','TF','TG']].applymap(map_strings_to_indexes)
training_data, val_data = train_test_split(all_samples, test_size=1-args.train_size-0.001, stratify=all_samples['label'],random_state=args.seed)

# Construct the model
mymodel = HGNNPredictor(
    in_channels = input_tensor.size()[1],
    hgnn_channels = args.hgnn_dims,
    linear_channels = args.linear_dims
).to(device)
        
mymodel = train(args, input_tensor, mymodel, hg, training_data, val_data, device)

project = args.project
tt = args.target_type
fname = '/home/jiawen/myMLnet/results/'+project+'/model_'+tt+'.pth'
torch.save(mymodel.state_dict(), fname)

mymodel.eval()
pred_results = input_all[['Receptor','TF','TG']]
pred_tensor = torch.tensor(pred_results.values).to(device)
predictions = mymodel(input_tensor,hg,pred_tensor)
predictions = F.relu(predictions)
pred_results['pred_label'] = predictions.cpu().detach().numpy()

mapping = {i: gene_all[i] for i in range(len(gene_all))}
for col in pred_results.columns[:-1]:
    pred_results[col] = pred_results[col].map(mapping)

pred_threshold = args.thres[3]
pred_filtered = pred_results.loc[pred_results['pred_label']>pred_threshold,:]

if args.TG == 1:
    TGs = pd.read_csv('/inputs/TG.csv',header=0, sep = ' ')
    TGs = TGs.values.flatten().tolist()
    pred_filtered = pred_filtered[pred_filtered['TG'].isin(TGs)]

pred_filtered = pred_filtered.reset_index(drop=True)


def XGBmodel_single(cur_row,exp):

    cur_TF = cur_row.iloc[0,1]
    cur_TG = cur_row.iloc[0,2]
    Recs_target = cur_row['Receptor'].tolist()
    exp_TF = np.array(exp[cur_TF])
    exp_TG = np.array(exp[cur_TG])

    exp_TFTG = exp_TF*exp_TG
    # exp_TFTG = np.log2(1+exp_TFTG)
    exp_recs = exp[Recs_target]
    # exp_recs = np.log2(1+exp_recs)

    xgb_model=xgb.XGBRegressor(n_estimators=1000, learning_rate=0.1, max_depth=5, random_state=2024,device = 'cuda')
    xgb_model.fit(exp_recs,exp_TFTG)

    importance = xgb_model.get_booster().get_score(importance_type='gain')
    df_imp = pd.DataFrame(list(importance.items()), columns=['Receptor', 'importance'])
    df_imp['TF'] = cur_TF
    df_imp['TG'] = cur_TG
    df_imp = df_imp[['Receptor','TF','TG','importance']]

    return df_imp

def PCR_single(cur_row,exp,n_comp = 5):

    cur_TF = cur_row.iloc[0,1]
    cur_TG = cur_row.iloc[0,2]
    Recs_target = cur_row['Receptor'].tolist()
    exp_TF = np.array(exp[cur_TF])
    exp_TG = np.array(exp[cur_TG])

    exp_TFTG = exp_TF*exp_TG
    exp_TFTG = np.log2(1+exp_TFTG)
    exp_recs = exp[Recs_target].values
    exp_recs = np.log2(1+exp_recs)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(exp_recs)
    # Perform PCA
    pca = PCA()
    X_train_pca = pca.fit_transform(X_scaled)
    X_train_pca_reduced = X_train_pca[:, :n_comp]
    model = LinearRegression()
    model.fit(X_train_pca_reduced, exp_TFTG)

    # Compute the original variable coefficients
    pc_coefficients = model.coef_
    loadings = pca.components_[:n_comp].T
    importance = loadings @ pc_coefficients

    return tuple(importance)

def XGBmodel_loop(df, exp):
    # Ensure the DataFrame has a unique index
    df = df.reset_index(drop=True)
    exp_nor = NormalizeData(exp,log2 = False)

    # Grouping DataFrame by 'y' and 'z'
    grouped = df.groupby(['TF', 'TG'])
    
    print("Start XGBoost regression ...")
    
    # with mp.Pool(24) as pool:
    #     results = pool.starmap(RFmodel_single, tasks)
    count = 0
    res_list = [None for _ in range(len(grouped))]
    with tqdm(total=len(grouped)) as pbar:
        pbar.set_description('Processing:')
        for _,group in grouped:
            imp = XGBmodel_single(group,exp_nor)
            res_list[count] = imp
            count+=1
            pbar.update(1)

    # Flatten the list of results and sort by the original index
    flat_results = pd.concat(res_list)

    return flat_results

results = XGBmodel_loop(pred_filtered.copy(), input_exp.copy().T)

fname = '/results/result_demo.csv'
results.to_csv(fname, index = False)
