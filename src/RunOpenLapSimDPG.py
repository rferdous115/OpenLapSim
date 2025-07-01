"""
---------------------------
OpenLapSimDPG - OLS
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
import dearpygui.dearpygui as dpg

# import packages (OLP)
from AccEnvCalc import AccEnvCalc
from LapTimeSimCalc import LapTimeSimCalc
from LapTimeSimCalcQSS import LapTimeSimCalcQSS
from PostProcDPG import PostProcDPG
from SetupFileLoader import SetupFileLoader
from BatchSetupFileLoader import BatchSetupFileLoader

class RunOpenLapSimDPG:

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
    
    def create_dpg_ui(self, vcar, dist, lapSimTimeDict):
        print("Creating DPG UI...")
        dpg.create_context()
        dpg.create_viewport(title='Custom Title', width=600, height=300)

        with dpg.window(tag="Primary Window"):
            dpg.add_text("Car Parameters")
            dpg.add_input_text(label="Mass (kg)", default_value="750", scientific=True)
            dpg.add_input_text(label="Lift Coefficient (-)", default_value="3.1", scientific=True)
            dpg.add_input_text(label="Drag Coefficient", default_value="1.0", scientific=True)
            dpg.add_input_text(label="Frontal Area (m^2)", default_value="1.0", scientific=True)
            dpg.add_input_text(label="Tyre Friction Coefficient Longitudinal", default_value="1.15", scientific=True)
            dpg.add_input_text(label="Grip Load Effective (% / kN of Fz)", default_value="0.10", scientific=True)
            dpg.add_input_text(label="Tire Radius", default_value="0.32", scientific=True)
            dpg.add_input_text(label="Gear Ratio", default_value="10.0, 7.8, 6.1, 7.8, 5.2, 4.5, 4.0", scientific=True)
            dpg.add_input_text(label="Drive Line Efficiency", default_value="0.95", scientific=True)
            dpg.add_input_text(label="Engine Torque [Nm]", default_value="200,300,430,380")
            dpg.add_input_text(label="Engine RPM", default_value="0,3000,7000,10000")
            dpg.add_input_text(label="Air Density (Kg/m3)", default_value="1.22", scientific=True)
            dpg.add_button(label="Select Track File", callback="Click")
            # dpg.add_slider_float(label="float", default_value=0.273, max_value=1)

            dist = np.arange(0, dist[len(dist)-1], 10)
            # vcar = [0, 40, 60, 80, 120]

            # print(len(dist))
            # print(len(vcar))

            dpg.add_button(label="Save") # NoF
            dpg.add_button(label="Simulate") # NoF

            vxacc = lapSimTimeDict["vxacc"]
            vxdec = lapSimTimeDict["vxdec"]
            vxcor = lapSimTimeDict["vxcor"]

            with dpg.plot(label="Velocity Profile", height=400, width=-1):
                                    
                # optionally create legend
                dpg.add_plot_legend()

                
                # FIX ME (DISTANCE) Distance Vs. Velocity Plot 
                # REQUIRED: create x and y axes
                dpg.add_plot_axis(dpg.mvXAxis, label="Distance (m)")
                dpg.add_plot_axis(dpg.mvYAxis, label="Velocity (km/h)", tag="y_axis_1")
                dpg.add_line_series(dist, vcar, label="Velocity Profile", parent="y_axis_1")


            with dpg.plot(label="Laptiem Simulation - Extra", height=400, width=-1):
                                    
                # optionally create legend
                dpg.add_plot_legend()
            
                #  Lap Time Simulation Extra Stuff 

                dpg.add_plot_axis(dpg.mvXAxis, label="Distance (m)")
                dpg.add_plot_axis(dpg.mvYAxis, label="Velocity (km/h)", tag="y_axis_2")
                dpg.add_line_series(dist, vcar, label="Velocity Profile", parent="y_axis_2")
                dpg.add_line_series(dist, vxcor, label="vxcor", parent="y_axis_2")
                dpg.add_line_series(dist, vxacc, label="vxacc", parent="y_axis_2")
                dpg.add_line_series(dist, vxdec, label="vxdec", parent="y_axis_2")

            
            dpg.add_button(label="Save") # NoF

        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.set_primary_window("Primary Window", True)
        dpg.start_dearpygui()
        dpg.destroy_context()

    def run(self):
        print("---------------------------")
        print("OpenLapSimDPG")
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
                aE = AccEnvCalc(s.setupDictList[elem])
                aE.Run()
            
                # Run Lap time Simulation
                trackFile = (self.trackFilesPath+self.trackFileName)
                l1 = LapTimeSimCalcQSS(trackFile, aE.accEnvDict, 10)
                l1.Run()
                l2 = LapTimeSimCalcQSS(trackFile, aE.accEnvDict,
                                    l1.lapTimeSimDict["vxaccEnd"])
                l2.Run()

                # set output channels from simulation for Export
                vcar = l2.lapTimeSimDict["vcar"]  # car speed [m/s]
                dist = l2.lapTimeSimDict["dist"]  # circuit dist [m]

                # export
                if self.bExport == 1:
                    RunOpenLapSimDPG.createExportSimFile(vcar, dist, self.exportFilesPath)

                # Computation time end
                tend = time.time()
                tcomp = round(tend - tstart, 1)
                print("Computational time: ", tcomp)

                # Post Processing
        
                pP = PostProcDPG(aE.accEnvDict, l2.lapTimeSimDict)
                pP.DPG()
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
            aE = AccEnvCalc(s.setupDict)
            aE.Run()

            # Run Lap time Simulation
            trackFile = (self.trackFilesPath+self.trackFileName)
            l1 = LapTimeSimCalcQSS(trackFile, aE.accEnvDict, 10)
            l1.Run()
            l2 = LapTimeSimCalcQSS(trackFile, aE.accEnvDict,
                                l1.lapTimeSimDict["vxaccEnd"])
            l2.Run()

            # set output channels from simulation for Export
            vcar = l2.lapTimeSimDict["vcar"]  # car speed [m/s]
            dist = l2.lapTimeSimDict["dist"]  # circuit dist [m]

            # export
            if self.bExport == 1:
                RunOpenLapSimDPG.createExportSimFile(vcar, dist, self.exportFilesPath)

            # Computation time end
            tend = time.time()
            tcomp = round(tend - tstart, 1)
            print("Computational time: ", tcomp)

            # Post Processing

            # pP = PostProcDPG(aE.accEnvDict, l2.lapTimeSimDict)
            # pP.DPG()
            # pP.printData()
            # if self.bPlot == 1:
            #     # pP.plotAccEnv()
            #     pP.plotGGV()
            #     pP.plotLapTimeSim()
            # if self.bPlotExtra == 1:
            #     pP.plotLapTimeSimExtra()
            #     pP.plotAccEnvExtra()
            # plt.show()  # plot all figure once at the end

            # output values
            self.laptime = l2.lapTimeSimDict["laptime"]
            self.vcarmax = l2.lapTimeSimDict["vcarmax"]
            self.tcomp = tcomp

            self.create_dpg_ui(vcar, dist, l2.lapTimeSimDict)

# ----------------------------------------------------------------------------


if __name__ == '__main__':

    # SetupFile.json 
    setupFileName = "SetupFile.json"
    # TrackFile.txt  
    trackFileName = "TrackFile.txt"
    # Batch Run 
    doBatch = 0
    batchWhich = "mcar" #choose out of mcar, clt, cx, afrcar, mbrk, gripx, gripy, loadEff, rtyre, reff, rho. No array parameters for now. 
    batchParams = [720, 740, 2] #min, max, step
    # Additional Options
    bExport = 1  
    bPlot = 1 
    bPlotExtra = 1 

    # object instantiation
    runOpenLapSimDPG = RunOpenLapSimDPG(setupFileName, trackFileName,
                                  bExport, bPlot, bPlotExtra, doBatch, batchWhich, batchParams)
    runOpenLapSimDPG.run()