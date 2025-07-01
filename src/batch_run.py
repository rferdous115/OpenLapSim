import os
import itertools
import json 
import pandas as pd

from AccEnvCalc import AccEnvCalc
from AccEnvCalcBM import AccEnvCalcBM
from LapTimeSimCalc import LapTimeSimCalc
from PostProc import PostProc
from SetupFileLoader import SetupFileLoader
from BatchSetupFileLoader import BatchSetupFileLoader

def make_param_grid(grid_dict):
    """
    Given a dict of paramter: list_of_values, returns a list of parameter dicts for every combination. 
    Example input: 
    {'mu': [1.0, 1.2], 'clt': [0.5, 1.0], 'cx':[0.3, 0.4]}
    Output: [
        {'mu':1.0, 'clt':0.5, 'cx':0.3},
        {'mu':1.0, 'clt':0.5, 'cx':0.4},
        {'mu':1.0, 'clt':1.0, 'cx':0.3},
        ...
    ]
    """
    keys = list(grid_dict.keys())
    values_product = itertools.product(*(grid_dict[key] for key in keys))
    return [dict(zip(keys,vals)) for vals in values_product]

def run_one_case(params, output_folder):
    return

# def loadJSON(batchParams):
#     batchParamArray = np.arange(batchParams[0], batchParams[1]+batchParams[2], batchParams[2])
#     batchParamArray=batchParamArray.tolist()
#     print(f"batchParamArray before BatchSetupFileLoader: {batchParamArray}")
#     s = BatchSetupFileLoader(self.setupFilesPath + self.setupFileName, batchParamArray, batchWhich)
#     s.loadBatchJSON()

def main():
    # s = BatchSetupFileLoader(self.setupFilesPath + self.setupFileName, batchParamArray, batchWhich)
    # s.loadBatchJSON()

    # Need to import setupdict. or bring this into runopenlapsimBM.

    grid_dict = {
        'mu':   [1.0, 1.1, 1.2],
        'clt':  [0.5, 1.0],
        'cx':   [0.30, 0.35]
    }

    # Above needs to be userinput with a start stop end

    param_list = make_param_grid(grid_dict)
    print(param_list)

    # Load individual json files. 
    s = BatchSetupFileLoader("setupFiles", param_list)

    

if __name__ == "__main__":
    main()
