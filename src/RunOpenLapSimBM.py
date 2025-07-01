"""
---------------------------
OpenLapSim - OLS
---------------------------

This is a steady state Lap Time Simulator for a simple point mass vehicle
with aero forces, constant tyre grip(x and y), engine torque map and gears.

Steps:
    1 - Select Files: TrackFile.txt and SetupFile.py
    2 - Calculate the Acceleration Envelope
    3 - Calculate the Lap Time Simulation (vcar)
    4 - Plot Results

---------------------------
@autor: Davide Strassera
@first release: 2019-12-21
by Python 3.7
---------------------------

"""

# ----------------------------------------------------------------------------

# import packages generic
import datetime
import matplotlib.pyplot as plt
import time
import numpy as np
# import dearpygui.dearpygui as dpg
import math
from math import sin, cos, log10
import os
import itertools
import json 
import pandas as pd

# import packages (OLP)
from AccEnvCalc import AccEnvCalc
from AccEnvCalcBM import AccEnvCalcBM
from LapTimeSimCalc import LapTimeSimCalc
from PostProc import PostProc
from SetupFileLoader import SetupFileLoader
from BatchSetupFileLoader import BatchSetupFileLoader
from BatchSetupFileLoader import BatchSetupFileLoader2

class RunOpenLapSim:

    def __init__(self, setupFileName, trackFileName,
                 bExport, bPlot, bPlotExtra, doBatch, batchWhich, batchParams):
        # inputs
        self.setupFileName = setupFileName
        self.trackFileName = trackFileName
        self.doBatch = doBatch
        self.batchWhich = batchWhich
        self.batchParams = batchParams
        self.bExport = bExport
        self.bPlot = bPlot
        self.bPlotExtra = bPlotExtra
        self.trackFilesPath = "trackFiles/"
        self.exportFilesPath = "exportFiles/"
        self.setupFilesPath = "setupFiles/"
        # outputs
        self.laptime = None
        self.vcarmax = None
        self.tcomp = None  # computational time

    @staticmethod
    def createExportSimFile(vcar, dist, exportFilesPath):
        time = datetime.datetime.now()
        timestrf = time.strftime("%b-%d-%Y")
        NewExportFileName = (exportFilesPath + "SimExport_"
                             + str(timestrf) + ".txt")
        newFile = open(NewExportFileName, "w")

        for i in range(len(dist)):
            lineToWrite = (str(dist[i]) + "\t" + str(vcar[i]) + "\n")
            newFile.write(lineToWrite)
        newFile.close()
        return NewExportFileName

    def run(self):
        print("---------------------------")
        print("OpenLapSim")
        print("---------------------------")

        # Computation time start
        tstart = time.time()

        # If Batch Run is chosen
        if self.doBatch == 1:
            batchParamArray = np.arange(self.batchParams[0], self.batchParams[1]+self.batchParams[2], self.batchParams[2])
            batchParamArray=batchParamArray.tolist()
            print(f"batchParamArray before BatchSetupFileLoader: {batchParamArray}")
            s = BatchSetupFileLoader(self.setupFilesPath + self.setupFileName, batchParamArray, batchWhich)
            s.loadBatchJSON()

            lapTimeArray = [0] * len(batchParamArray)
            
            # Run Acceleration Envelope
            for elem, val in enumerate(batchParamArray):
                print(f"elem {elem} -> {s.setupDictList[elem]}")
                aE = AccEnvCalcBM(s.setupDictList[elem])
                aE.Run()
            
                # Run Lap time Simulation
                trackFile = (self.trackFilesPath+self.trackFileName)
                l1 = LapTimeSimCalc(trackFile, aE.accEnvDict, 10)
                l1.Run()
                l2 = LapTimeSimCalc(trackFile, aE.accEnvDict,
                                    l1.lapTimeSimDict["vxaccEnd"])
                l2.Run()

                # set output channels from simulation for Export
                vcar = l2.lapTimeSimDict["vcar"]  # car speed [m/s]
                dist = l2.lapTimeSimDict["dist"]  # circuit dist [m]

                # export
                if self.bExport == 1:
                    RunOpenLapSim.createExportSimFile(vcar, dist, self.exportFilesPath)

                # Computation time end
                tend = time.time()
                tcomp = round(tend - tstart, 1)
                print("Computational time: ", tcomp)

                # Post Processing
        
                pP = PostProc(aE.accEnvDict, l2.lapTimeSimDict)
                lapTimeArray[elem] = pP.laptime
                print(f"LapTime {elem}: {lapTimeArray[elem]}")
                pP.printData()
                if self.bPlot == 1:
                    # pP.plotAccEnv()
                    pP.plotGGV()
                    pP.plotLapTimeSim()
                if self.bPlotExtra == 1:
                    pP.plotLapTimeSimExtra()
                    pP.plotAccEnvExtra()
                plt.show()  # plot all figure once at the end

                plt.plot(batchParamArray, lapTimeArray, 'bo')
                plt.xlabel(f"{batchWhich} (kg)")
                plt.ylabel("laptime (s)")
                plt.title(f"{batchWhich} vs laptime")

                # output values
                self.laptime = l2.lapTimeSimDict["laptime"]
                self.vcarmax = l2.lapTimeSimDict["vcarmax"]
                self.tcomp = tcomp

                print(f"Lap Times: {lapTimeArray}")

        else:
            # SetupFile obj instantiation
            s = SetupFileLoader(self.setupFilesPath + self.setupFileName)
            s.loadJSON()

            # Run Acceleration Envelope
            aE = AccEnvCalcBM(s.setupDict)
            aE.Run()

            # Run Lap time Simulation
            trackFile = (self.trackFilesPath+self.trackFileName)
            l1 = LapTimeSimCalc(trackFile, aE.accEnvDict, 10)
            l1.Run()
            l2 = LapTimeSimCalc(trackFile, aE.accEnvDict,
                                l1.lapTimeSimDict["vxaccEnd"])
            l2.Run()

            # set output channels from simulation for Export
            vcar = l2.lapTimeSimDict["vcar"]  # car speed [m/s]
            dist = l2.lapTimeSimDict["dist"]  # circuit dist [m]

            # export
            if self.bExport == 1:
                RunOpenLapSim.createExportSimFile(vcar, dist, self.exportFilesPath)

            # Computation time end
            tend = time.time()
            tcomp = round(tend - tstart, 1)
            print("Computational time: ", tcomp)

            # Post Processing

            pP = PostProc(aE.accEnvDict, l2.lapTimeSimDict)
            pP.printData()
            if self.bPlot == 1:
                # pP.plotAccEnv()
                pP.plotGGV()
                pP.plotLapTimeSim()
                pP.plotLapTimeSimAxAcc()
                pP.plotTrackMap()
            if self.bPlotExtra == 1:
                pP.plotLapTimeSimExtra()
                pP.plotAccEnvExtra()
                pP.AV()
            plt.show()  # plot all figure once at the end
            plt.ion()

            # output values
            self.laptime = l2.lapTimeSimDict["laptime"]
            self.vcarmax = l2.lapTimeSimDict["vcarmax"]
            self.tcomp = tcomp

# ----------------------------------------------------------------------------

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

if __name__ == '__main__':

    # SetupFile.json 
    setupFileName = "SetupFile.json"
    # TrackFile.txt  
    trackFileName = "TrackFile.txt"
    # Batch Run (Redundant)
    doBatch = 1
    batchWhich = "mcar" #choose out of mcar, clt, cx, afrcar, mbrk, gripx, gripy, loadEff, rtyre, reff, rho. No array parameters for now. 
    batchParams = [720, 740, 2] #min, max, step

    # Batch Run (New)
    grid_dict = {
        'mu':   [1.0, 1.1, 1.2],
        'clt':  [0.5, 1.0],
        'cx':   [0.30, 0.35]
    }

    param_list = make_param_grid(grid_dict) # list of dictionaries
    print(param_list) # 

    # Populate new JSONS
    s = BatchSetupFileLoader2("self.setupFilesPath" + setupFileName, param_list)
    



    # Additional Options
    bExport = 1  
    bPlot = 1 
    bPlotExtra = 1 

    # object instantiation
    runOpenLapSim = RunOpenLapSim(setupFileName, trackFileName,
                                  bExport, bPlot, bPlotExtra, doBatch, batchWhich, batchParams)
    runOpenLapSim.run()