import dearpygui.dearpygui as dpg
import math
from math import sin, cos, log10
# Import packages
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
from PostProc import PostProc

class PostProcDPG:
    def __init__(self, accEnvDict, lapSimTimeDict):
        self.size = 10
        # inputs
        self.GGVfull = accEnvDict["GGVfull"]
        self.GGVacc = lapSimTimeDict["GGVacc"]
        self.GGVdec = lapSimTimeDict["GGVdec"]

        self.vxvect = accEnvDict["vxvect"]
        self.ay = accEnvDict["ay"]
        self.axacc = accEnvDict["axacc"]
        self.axdec = accEnvDict["axdec"]
        self.dist = lapSimTimeDict["dist"]
        self.vcar = lapSimTimeDict["vcar"]
        self.laptime = lapSimTimeDict["laptime"]
        self.vcarmax = lapSimTimeDict["vcarmax"]
        self.vxacc = lapSimTimeDict["vxacc"]
        self.vxdec = lapSimTimeDict["vxdec"]
        self.vxcor = lapSimTimeDict["vxcor"]
        # extra channels (AccEnv)
        self.nGear = accEnvDict["nGear"]
        self.EngNm = accEnvDict["EngNm"]
        self.EngRpm = accEnvDict["EngRpm"]
        self.Fzaero = accEnvDict["Fzaero"]
        self.Fxaero = accEnvDict["Fxaero"]
        self.Fxgrip = accEnvDict["Fxgrip"]
        self.Fxdrive = accEnvDict["Fxdrive"]

    # def plotLapTimeSim(self):
    #     plt.figure(2, figsize=(self.size, self.size/2))
    #     plt.title("OpenLapSim - Lap Time Simulation")
    #     plt.plot(self.dist, self.vcar, 'b-', linewidth=2, label="vcar")
    #     plt.xlabel('distance [m]')
    #     plt.ylabel('velocity [m/s]')
    #     plt.legend()
    #     plt.grid(b=True, which='major', linestyle=':')
    #     plt.ylim(0, self.vcarmax*1.2)
    #     plt.xlim(0, max(self.dist))

    def DPG(self):
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

            with dpg.plot(label="Velocity Profile", height=400, width=-1):
                                    
                # optionally create legend
                dpg.add_plot_legend()

                # REQUIRED: create x and y axes
                dpg.add_plot_axis(dpg.mvXAxis, label="Distance (m)")
                
                with dpg.plot_axis(dpg.mvYAxis, label="Velocity (km/h)"):

                    # series belong to a y axis
                    dpg.add_line_series(self.dist, self.vcar, label="Velocity Profile")

            

            dpg.add_button(label="Save") # NoF

        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.set_primary_window("Primary Window", True)
        dpg.start_dearpygui()
        dpg.destroy_context()
