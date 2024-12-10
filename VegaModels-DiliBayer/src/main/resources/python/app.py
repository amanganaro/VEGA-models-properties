import sys
import os
import pandas as pd
import numpy as np
import json
import torch
import tempfile
from torch.autograd import Variable
from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem.MolStandardize.standardize import Standardizer
from models.model_mtnn import MultiTaskNN


# Endpoint list
ASSAYS_TASKS = ['BSEPi', 'BSEPs', 'PGPi', 'PGPs', 'MRP4i', 'MRP3i', 'MRP3s', 'MRP2i', 'MRP2s', 'BCRPi', 'BCRPs', 'OATP1B1i', 'OATP1B3i', 'NRF2', 'LXR', 'AHR', 'PPARa', 'PPARg', 'PXR', 'FXR', 'MTX_MP', 'MTX_RC', 'MTX_FOM', 'PLD', 'PLD_HTS', 'HTX', 'ERS', 'ARE']
DILI_TASKS = ['DILI_majority', 'DILI_sensitive', 'DILI_secure']

standardizer = Standardizer(max_tautomers=10)
include_stereoinfo = False

def standardizeMol(mol):
    """
    Standardizer inspired by MELLODDY consortium
    :param mol: rdkit molecule object
    :return: cleaned rdkit molecule object
    """
    if mol is None:
        return mol
    try:
        mol = standardizer.charge_parent(mol)
        mol = standardizer.isotope_parent(mol)
        if include_stereoinfo is False:
            mol = standardizer.stereo_parent(mol)
        mol = standardizer.tautomer_parent(mol)
        return standardizer.standardize(mol)
    except Chem.rdchem.AtomValenceException:
        return None


def prepare_structures(input_df, smiles_col='smiles'):
    '''
    Standardize compounds, get canonical smiles and deduplicate structures based on canonical smiles
    '''
    df = input_df.copy()
    PandasTools.AddMoleculeColumnToFrame(df, smiles_col,'molecule', includeFingerprints=False)
    len_1 = len(df)
    df.dropna(subset=['molecule'], axis=0, inplace=True)
    len_2 = len(df)
    print(f'No. of missing molecules: {len_1-len_2}')

    # Standardize mols and get canonical smiles
    for idx in df.index:
        try:
            stand_mol = standardizeMol(df.loc[idx, 'molecule'])
            df.loc[idx, 'canonical_smiles'] = Chem.MolToSmiles(stand_mol, canonical=True)
        except:
            print(f'Failed to standardize {df.loc[idx, smiles_col]}')
            df.loc[idx, 'canonical_smiles'] = None
    
    df.drop(['molecule'], axis=1, inplace=True)
    df.dropna(subset=['canonical_smiles'], inplace=True)
    return df

def get_mean_prediction_ensemble(pred_task, tasks):
    ### Average predictions and get STD for uncertainty estimation
    mean_pred = pd.DataFrame()
    for task in tasks:
        pred_task_mean = np.mean(pred_task[task], axis=0)
        pred_task_std = np.std(pred_task[task], axis=0)
        mean_pred[f'{task}_mean'] = pred_task_mean
        mean_pred[f'{task}_std'] = pred_task_std
        # binary class label
        mean_pred[f'{task}_class'] = (pred_task_mean >= 0.5).astype(int)
    return mean_pred

def make_predictions_with_ensemble(Xtest):
    ### apply ensemble model on ONTOX tasks to test set
    tasks = DILI_TASKS + ASSAYS_TASKS
    # read HPs from json file
    script_dir = os.path.dirname(os.path.realpath(__file__))
    hp_path = os.path.join(script_dir, './models/hyperparameters.json')
    with open(hp_path, 'r') as openfile:
        hparam_dict = json.load(openfile)

        
    pred_task = {}
    for split_id in range(5): # ensemble with 5 models (one per split)
        print(f'Split: {split_id}')
        # load models and make predictions
        script_dir = os.path.dirname(os.path.realpath(__file__))
        model_file = f'models/model_mtnn_all_tasks_cddd_{split_id}.pth'
        model_file_path = os.path.join(script_dir, model_file)
        hparam = hparam_dict[str(split_id)]

        device='cpu'
        model = MultiTaskNN(input_size=512, params=hparam, n_tasks=len(tasks))
        with open(model_file_path, 'rb') as f:
            model.load_state_dict(torch.load(f, map_location=device))
        #model.load_state_dict(torch.load(model_file, map_location=device))
        model.to(device)
        model.eval()

        # get predictions from single model
        test_data = Variable(torch.from_numpy(Xtest).float())
        preds_list = model(test_data)
        predictions = [torch.Tensor.cpu(p).detach().numpy() for p in preds_list]
        predictions = np.array(predictions).squeeze().T

        # Handle cases where predictions might be 1D
        if len(predictions.shape) == 1:
            predictions = predictions.reshape(-1, len(tasks))

        # calculate predictions for all tasks
        for i, task in enumerate(tasks):
            if split_id == 0:
                pred_task[task] = predictions[:,i]
            else:
                pred_task[task] = np.vstack((pred_task[task], predictions[:,i]))
        
    # get mean prediction and STD from ensemble
    mean_pred = get_mean_prediction_ensemble(pred_task, tasks)
    return mean_pred, pred_task

def main():
    file_location = sys.argv[1]
    file_output_location = sys.argv[2]
    test_df = pd.read_csv(file_location)
    test_df = prepare_structures(test_df, smiles_col='smiles')
    x_test = test_df[[c for c in test_df if 'cddd' in c]].values
    mean_pred, pred_task = make_predictions_with_ensemble(x_test)
    test_pred = test_df.merge(mean_pred, right_index=True, left_index=True)
    test_pred.to_csv(file_output_location, index=False)


if __name__ == "__main__":
    main()
