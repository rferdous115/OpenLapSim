"""
---------------------------
Setup File - OLS
---------------------------

This class loads the setupFile.json and creates a dictionary.
Below an example of setupFile.json (REMOUVE comments "#..." on real JSON file)

{
    "setupName" : "Gp2Dummy",
    "mcar"		: 728,        #[Kg]; total car mass
    "clt"		: 3.1,        #[100 pt.]; Lift coefficient (-)
    "cx"		: 1.0,        #[100 pt.]; Drag coefficient
    "afrcar"	: 1.0,        #[m2]; Frontal Area
    "mbrk"		: 7000,       #[Nm]; Max Braking Torque
    "gripx"		: 1.15,       #tyre friction coeff long
    "gripy"		: 1.40,       #tyre friction coeff lat
    "loadEff"   : 0.10,       #grip Load Effect % / 1KN of Fz
    "rtyre"		: 0.32,       #[m]; tyre radius
    "rGearRat"	: [10.0,7.8,6.1,7.8,5.2,4.5,4.0],  #Final Gear Ratio
    "reff"		: 0.95,       # drive line efficiency
    "EngNm"     : [200,300,430,380], # [Nm]; Engine Torque
    "EngRpm"    : [0,3000,7000,10000],  # [rpm]Engine rmp
    "rho"		: 1.22        #[Kg/m3]; air density
}

---------------------------
@author: Davide Strassera
@first release: 2019-12-21
by Python 3.7
---------------------------

"""

# Import Packages
import json

class BatchSetupFileLoader:

    def __init__(self, setupFileName, batchParamArray, batchWhich):
        self.setupFileName = setupFileName
        self.setupDict = {}
        self.batchParamArray = batchParamArray
        self.batchWhich = batchWhich
        self.setupDictList = [0] * len(self.batchParamArray)

    def batchWriteJSON(self, data, filename):
        with open(filename, "w") as f:
            json.dump(data, f, indent=4)

    def loadBatchJSON(self):
        #load json
        # with open(self.setupFileName) as f:
        #     data = json.load(f)

        for key, val in enumerate(self.batchParamArray):
            with open(self.setupFileName) as f:
                data = json.load(f)
            data[self.batchWhich] = val
            batch_file_name_list = list(str(self.setupFileName))
            batch_file_name_list.insert(-5, str(key))
            batch_file_name = ''.join(batch_file_name_list)
            print(batch_file_name)
            self.batchWriteJSON(data, batch_file_name)
            self.setupDictList[key] = data
            
        print(self.setupDictList)
        # set setupDict
        # self.setupDict = data

class BatchSetupFileLoader2:

    def __init__(self, setupFileName, param_list):
        self.setupFileName = setupFileName
        self.setupDict = {}
        self.param_list = param_list
        self.setupDictList = [0] * len(self.param_list[0])

    def batchWriteJSON(self, data, filename):
        with open(filename, "w") as f:
            json.dump(data, f, indent=4)

    def loadBatchJSON(self):
        #load json
        # with open(self.setupFileName) as f:
        #     data = json.load(f)
        i = 0
        for param in self.param_list:
            print(f"param {param}")
            for key, val in param.items():
                # for v in val:
                with open(self.setupFileName) as f:
                    data = json.load(f)
                print(f"Inputting {key} with vals {val}")
                data[key] = val
                batch_file_name_list = list(str(self.setupFileName))
                batch_file_name_list.insert(-5, str(i))
                batch_file_name = ''.join(batch_file_name_list)
                # print(batch_file_name)
                self.batchWriteJSON(data, batch_file_name)
                self.setupDictList[key] = data
                i += 1
                
                    # print(self.setupDictList)
                    # set setupDict
                    # self.setupDict = data
