import pandas as pd
import torch
import numpy as np
import time
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import random

from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score,average_precision_score
from sklearn.model_selection import train_test_split


def Train(args, Exp, model, hypergraph, samples, device):

    # Convert dataframe to torch tensor
    Exp_tensor = torch.tensor(Exp.values, dtype=torch.float32)
    Exp_tensor = Exp_tensor.to(device)

    hypergraph = hypergraph.to(device)
    model = model.to(device)

    # Form dataloaders
    training_data, val_data = train_test_split(samples, test_size=1-args.train_size-0.001, 
                                               stratify=samples['label'], random_state=args.seed)
    training_dataset = DataFrameDataset(training_data)
    training_load = DataLoader(training_dataset, batch_size=args.batch_size, shuffle=True)
    val_dataset = DataFrameDataset(val_data)
    val_tensor = torch.tensor(val_data.values)

    optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
    # optimizer = torch.optim.SGD(net.parameters(), lr=0.1)
    scheduler = optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.975)

    print('\n ----- Start training ... -----\n')


    for epoch_i in range(args.epochs):
        
        print('Starting [ Epoch', epoch_i+1, 'of', args.epochs, '] ...')
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

def Predict(Exp, model, input_samples, hypergraph, genes, device):

    input_tensor = torch.tensor(Exp.values, dtype=torch.float32)
    input_tensor = input_tensor.to(device)
    hypergraph = hypergraph.to(device)
    model.eval()
    pred_results = input_samples[['Receptor','TF','TG']]
    pred_tensor = torch.tensor(pred_results.values).to(device)
    predictions = model(input_tensor,hypergraph,pred_tensor)
    predictions = F.relu(predictions)
    pred_results['pred_label'] = predictions.cpu().detach().numpy()

    mapping = {i: genes[i] for i in range(len(genes))}
    for col in pred_results.columns[:-1]:
        pred_results[col] = pred_results[col].map(mapping)

    return pred_results

'''
# This function has been deprecated
def train_old(args, Exp_tensor, model, hypergraph, training_data,val_data, device):

    # Form dataloaders
    training_dataset = DataFrameDataset(training_data)

    def seed_worker(worker_id):
        worker_seed = torch.initial_seed() % 2**32
        np.random.seed(worker_seed)
        random.seed(worker_seed)
        
    g = torch.Generator()
    g.manual_seed(0)

    training_load = DataLoader(
        training_dataset, 
        batch_size=args.batch_size,
        num_workers=4,
        worker_init_fn=seed_worker,
        generator=g,
        shuffle=True
    )
    # val_dataset = DataFrameDataset(val_data)
    val_tensor = torch.tensor(val_data.values)

    epochs = args.epochs
    optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
    # optimizer = torch.optim.SGD(net.parameters(), lr=0.1)
    scheduler = optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.975)

    print('\n ----- Start training ... -----\n')

    # Record the metrics
    t_loss = []
    v_loss = []
    auc = []
    aupr = []
    times = []

    for epoch_i in range(epochs):

        if args.epoch_output:
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

        t_loss.append(running_loss)
        v_loss.append(val_loss.cpu().item())
        auc.append(AUC)
        aupr.append(AUPR)
        times.append(train_time + eval_time)

        if args.epoch_output:
            print(' - Training BCE loss: {:.4f}'.format(running_loss),', Validation BCE loss:{:.3F}'.format(val_loss))
            print(' - AUPR:{:.3F}'.format(AUPR),', AUROC: {:.3F}'.format(AUC), ', AUPR-norm: {:.3F}'.format(AUPR_norm))
            print(' - Elapsed training time:{:.3f}'.format(train_time),', validation time:{:.3f}'.format(eval_time))
            print(' - Current learning rate:{:.5f}'.format(cur_lr))
 
    print('\n ----- Training finished !! -----\n')


    if args.metrics_output:
        df = pd.DataFrame({
            'Epoch': list(range(1, epochs + 1)),
            'train_loss': t_loss,
            'validation_loss': v_loss,
            'AUROC': auc,
            'AUPR': aupr,
            'Running_Time': times
        })
        return model,df
    else:
        return model, AUC, AUPR
    
'''