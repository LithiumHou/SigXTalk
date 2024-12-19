import torch
import torch.nn.functional as F
import random
import sys
import argparse
import dhg

work_dir = '/home/jiawen/myMLnet/pythoncodes'
sys.path.append(work_dir)

from preprocessing import *
from predictor import *
from training import *

# Parse the args for hyperlink prediction

parser = argparse.ArgumentParser(description="Process some inputs.")

parser.add_argument('--project', type=str, help='The name of the output folder.')
parser.add_argument('--target_type', type=str,help='The name of the target cell type.')
parser.add_argument('--thres', type=float, default=[0.1,0.1], help='The thresholds for training set and hypergraph construction.')
parser.add_argument('--corr_type', type=str, default='Spearman', help='The method to calculate correlations, upper/lower case compatible.')
parser.add_argument('--train_size', type=float, default=0.66, help='Fraction of training samples.')
parser.add_argument('--hgnn_dims', type=int, default=[256, 128], help='The dimension of hidden layer.')
parser.add_argument('--linear_dims', type=int, default=[64, 32], help='The dimensions of the MLP layer.')
parser.add_argument('--lr', type=float, default= 0.005, help='Initial learning rate.')
parser.add_argument('--epochs', type=int, default= 100, help='Number of epochs.')
parser.add_argument('--batch_size', type=int, default=128, help='The size of each batch.')
parser.add_argument('--seed', type=int, default=1024, help='Random seed.')
parser.add_argument('--device', type=str, default='cuda', help='The device for training the neural network.')
parser.add_argument('--save_model', type=bool, default=False, help='Whether to save the model.')

args = parser.parse_args([])
device = torch.device(args.device if torch.cuda.is_available() else "cpu")

random.seed = args.seed
torch.manual_seed(args.seed)
np.random.seed(args.seed)

print(f"Dataset: {args.project}, target type: {args.target_type}")

# Load the input expression matrix
input_exp = pd.read_csv('/home/jiawen/myMLnet/pythoncodes/inputs/ExpressionCount.csv', header = 0, index_col=0, sep=' ') # row: genes, col: cells
input_exp = input_exp.astype(np.float32)
gene_all = input_exp.index.tolist()

# load the gene-gene interaction databases
RecTFDB = pd.read_csv('/home/jiawen/myMLnet/pythoncodes/inputs/RecTFDB.csv',header=0, sep = ' ')
TFTGDB = pd.read_csv('/home/jiawen/myMLnet/pythoncodes/inputs/TFTGDB.csv',header=0, sep = ' ')

# construct the hypergraph and generate the training samples
df_rt = calculate_corr(input_exp.T, RecTFDB, type = args.corr_type, abss = False)
df_tftg = calculate_corr(input_exp.T, TFTGDB, type = args.corr_type, abss = False)


hg, all_samples, input_all = Filter_RTTDB(RecTFDB, df_tftg, thres1 = args.thres[0], thres2 = args.thres[1], genes = gene_all,
                                            args = args, first = 'score', sample_scale = 0.75, ood_frac = 0.5)

# Construct the model
mymodel = HGNNPredictor(
    in_channels = input_exp.shape[1],
    hgnn_channels = args.hgnn_dims,
    linear_channels = args.linear_dims
)
        
mymodel = Train(args, input_exp, mymodel, hg, all_samples, device)

if args.save_model:
    fname = '/home/jiawen/myMLnet/results/'+args.project+'/model_'+args.target_type+'.pth'
    torch.save(mymodel.state_dict(), fname)

pred_results = Predict(Exp = input_exp, model = mymodel, input_samples = input_all,
                       hypergraph = hg, genes = gene_all, device = device)

# Save the results!
fname = '/home/jiawen/myMLnet/results/'+args.project+'/pathways_'+args.target_type+'.csv'
pred_results.to_csv(fname, index = False)

