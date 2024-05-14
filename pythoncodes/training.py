import pandas as pd
import torch
import numpy as np
import time
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import random

from torch.utils.data import Dataset, DataLoader, TensorDataset
from sklearn.metrics import roc_auc_score,average_precision_score


def train(args, Exp_tensor, model, hypergraph, training_data,val_data, device):

    random.seed = args.seed
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # Form dataloaders
    training_dataset = DataFrameDataset(training_data)
    training_load = DataLoader(training_dataset, batch_size=args.batch_size, shuffle=True)
    # val_dataset = DataFrameDataset(val_data)
    val_tensor = torch.tensor(val_data.values)

    epochs = args.epochs
    optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
    # optimizer = torch.optim.SGD(net.parameters(), lr=0.1)
    scheduler = optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.975)

    print('\n ----- Start training ... -----\n')

    for epoch_i in range(epochs):
        
        print('Starting [ Epoch', epoch_i+1, 'of', epochs, '] ...')
        start = time.time()
        running_loss = 0.0
        cur_lr = optimizer.param_groups[0]['lr']
        
        
        for train_x, train_y in training_load:
            
            model.train()
            optimizer.zero_grad()

            train_x = train_x.to(device)
            train_y = train_y.to(torch.float).to(device).view(-1,1)

            pred = model(Exp_tensor,hypergraph,train_x)
            pred = F.relu(pred)

            loss_value = F.binary_cross_entropy(pred, train_y)

            loss_value.backward()
            optimizer.step()
            # scheduler.step()

            running_loss += loss_value.item()

        scheduler.step()
        train_time = time.time()-start
        
        # Evaluation
        start = time.time()
        with torch.no_grad():
            model.eval()

            val_inputs = val_tensor[:,:-1].to(device)
            val_outputs = model(Exp_tensor, hypergraph, val_inputs)
            val_outputs = F.relu(val_outputs)
            val_y = val_tensor[:,-1]
            val_y = val_y.to(torch.float).to(device).view(-1,1)
            val_loss = F.binary_cross_entropy(val_outputs, val_y)
            AUC, AUPR, AUPR_norm = Evaluation(y_pred = val_outputs, y_true = val_y)

        eval_time = time.time()-start

        print(' - Training BCE loss: {:.4f}'.format(running_loss),', Validation BCE loss:{:.3F}'.format(val_loss))
        print(' - AUPR:{:.3F}'.format(AUPR),', AUROC: {:.3F}'.format(AUC), ', AUPR-norm: {:.3F}'.format(AUPR_norm))
        print(' - Elapsed training time:{:.3f}'.format(train_time),', validation time:{:.3f}'.format(eval_time))
        print(' - Current learning rate:{:.5f}'.format(cur_lr))
 
    print('\n ----- Training finished !! -----\n')

    if args.model_output:
        return model,AUC,AUPR
    else:
        return model


def Evaluation(y_pred, y_true):

    y_p = y_pred.cpu().detach().numpy()
    y_p = y_p.flatten()

    y_t = y_true.cpu().numpy().flatten().astype(int)

    AUC = roc_auc_score(y_true=y_t, y_score=y_p)
    AUPR = average_precision_score(y_true=y_t,y_score=y_p)
    AUPR_norm = AUPR/np.mean(y_t)


    return AUC, AUPR, AUPR_norm

class DataFrameDataset(Dataset):
    def __init__(self, dataframe):
        self.dataframe = dataframe

    def __len__(self):
        return len(self.dataframe)

    def __getitem__(self, idx):
        sample = self.dataframe.iloc[idx]
        # Process the sample if needed
        features = sample[['Receptor', 'TF', 'TG']].values
        labels = sample['label']
        return torch.tensor(features), torch.tensor(labels)

