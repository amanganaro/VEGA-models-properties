import pandas as pd

import pickle
import os
import argparse
import subprocess
import sys

# rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem


def check_smiles(smiles):
    try: Chem.MolFromSmiles(smiles)
    except: 
        print("invalid smiles.")
        sys.exit(1)

def ad_evaluation(smiles:str, radius=3, num_bits=1024, singlesmiles=False):
    """
    check similarity between compounds in dataset and target using morgan fingerprint as descriptors
    """
    fp1 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), radius, nBits=num_bits)
    ad_results = {}
    path_folder = os.path.join(file_dir, "data")
    path_folder_endpoint = os.path.join(path_folder, f'train_apical.csv')
    # name = folder.split("_")[0]
    data = pd.read_csv(path_folder_endpoint)
    smiles_list = list(data['smiles'])
    counter = 0
    for smi in smiles_list:
        # print(f"[DEBUGGING]: {smi}")
        try: 
            fp2 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), radius, nBits=num_bits)
            if counter < 3: 
                if DataStructs.DiceSimilarity(fp1, fp2) > 0.3: 
                    counter += 1
        except:
            print(f'error molecule {smi}')
            continue
        
        if counter < 3:
            ad_results["apical_data"] = 'out AD'
        else:
            ad_results["apical_data"] = 'in AD'

    if singlesmiles: return pd.DataFrame([ad_results], index=[smiles])
    else: return pd.DataFrame([ad_results])

def import_models():
    "import models find in models_apical folder"
    path = os.path.join(file_dir, "models_apical")

    models_path = [path]

    # import ML models
    models = {}

    for model_path in models_path:
        files = os.listdir(model_path)
        print("models_apical")
        for file in files:
            if "pkl" in file:
                path = os.path.join(model_path, file)
                with open(path, 'rb') as file:
                    models["models_apical"] = pickle.load(file)
            """
            elif "train" in file:
                path = os.path.join(model_path, file)
                columns = pd.read_excel(path)
                descriptors["models_apical"] = list(columns.keys()[2:])
            """
    
    return models

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='cardiotoxicity assessment')

    # Add an argument for input
    parser.add_argument('--output', required=True, help='CSV output file')
    parser.add_argument('--input', required=True, help='CSV cddd descriptors file input')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Access the input value
    input_file = args.input
    output_file = args.output
    # Your script logic here
    return {
        "input_file": input_file,
        "output_file": output_file
    }


def transform_data(target_CDDD:pd.DataFrame):
    """
    take as input the dictionary with descriptors for model.
    """
    target_CDDD.rename(columns=lambda x:
        str(int(x.split('_')[1]) - 1) if x.startswith('cddd_') and x.split('_')[1].isdigit()
        else x, inplace=True)
    # adjust columns name as str
    data = target_CDDD.loc[:, list(models['models_apical'].feature_names_in_)]

    return data

def inference(smiles, data:pd.DataFrame, singlemol=True):

    # Prediction
    pred = models['models_apical'].predict(data)
    if singlemol:
        results = {"models_apical": pred}

        results_final = {}
        for name, result in results.items():
            if result == 0: activity = 'Inactive'
            else: activity = 'Active'
            results_final[name] = activity   

        return pd.DataFrame(pd.DataFrame(results_final, index=[smiles]))
    else:
        results = pd.DataFrame([pred]).T.rename(columns={0:"ApicalModel"})
        # print(results)
        results['smiles'] = smiles
        return results

if __name__ == "__main__":

    file_dir=os.path.dirname(os.path.abspath(__file__))

    # ask smiles of chemicals to evaluate
    files = main()
    target_CDDD = pd.read_csv(files['input_file'])

    for smi in target_CDDD['smiles']:
        check_smiles(smi)

    # import models pipeline
    print("[INFO]: Models import...")
    models = import_models()
    print("done")

    print("[INFO]: AD_evaluation...")
    for n, smi in enumerate(target_CDDD['smiles']):
        ad = ad_evaluation(smi)
        if n == 0:
            ad_results = ad.copy()
        else:
            ad_results = pd.concat([ad_results, ad], axis=0)
    # debug
    ad_results.reset_index(inplace=True, drop=True)
    print("done")

    target_dict = transform_data(target_CDDD=target_CDDD)
    results = inference(smiles=target_CDDD['smiles'].tolist(),
                        data=target_dict, singlemol=False)

    results.reset_index(drop=True, inplace=True)
    results['smiles'] = target_CDDD['smiles'].tolist()
    merged_df = pd.concat([results, ad_results], axis=1)
    merged_df.to_csv(files["output_file"])
    print("Complete")
    


    
    
   