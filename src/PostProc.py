"""
---------------------------
Post Processing - OLS
---------------------------

This class computes the post processing of the simulation.

---------------------------
@autor: Davide Strassera
@first release: 2019-12-21
by Python 3.7
---------------------------

"""
# Import packages
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
from matplotlib.collections import LineCollection

class PostProc:

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


    def plotGGV(self, bPlotGGVfull=0):
        GGVacc = self.GGVacc
        GGVdec = self.GGVdec
        GGVfull = self.GGVfull

        xyz1 = GGVacc
        X1 = xyz1[:, 0]
        Y1 = xyz1[:, 1]
        Z1 = xyz1[:, 2]
        # Griddata
        ploty1, plotz1, = np.meshgrid(np.linspace(np.min(Y1), np.max(Y1), 30),
                                      np.linspace(np.min(Z1), np.max(Z1), 30))
        plotx1 = interp.griddata((Y1, Z1), X1, (ploty1, plotz1),
                                 method='linear', fill_value=0.0)

        xyz2 = GGVdec
        X2 = xyz2[:, 0]
        Y2 = xyz2[:, 1]
        Z2 = xyz2[:, 2]
        # Griddata
        ploty2, plotz2, = np.meshgrid(np.linspace(np.min(Y2), np.max(Y2), 30),
                                      np.linspace(np.min(Z2), np.max(Z2), 30))
        plotx2 = interp.griddata((Y2, Z2), X2, (ploty2, plotz2),
                                 method='linear', fill_value=0.0)

        xyz3 = GGVfull
        X3 = xyz3[:, 0]
        Y3 = xyz3[:, 1]
        Z3 = xyz3[:, 2]

        fig = plt.figure(5)
        ax = fig.add_subplot(111, projection='3d')
        surf1 = ax.plot_surface(plotx1, ploty1, plotz1,
                                cstride=1, rstride=1, cmap='coolwarm',
                                edgecolor='black', linewidth=0.2,
                                antialiased=True)
        surf2 = ax.plot_surface(plotx2, ploty2, plotz2,
                                cstride=1, rstride=1, cmap='coolwarm',
                                edgecolor='black', linewidth=0.2,
                                antialiased=True)
        if bPlotGGVfull == 1:
            ax.scatter(X3, Y3, Z3, color="black", label="GGVfull sparse")
            ax.legend()

        # Add a color bar which maps values to colors.
        fig.colorbar(surf1, shrink=0.5, aspect=5)

        plt.title("OpenLapSim - Performance Envelope")
        plt.xlabel('ax [m/s^2]')
        plt.ylabel('ay [m/s^2]')
        ax.set_zlabel('velocity [m/s]')

    def AV(self):
        GGVacc = self.GGVacc
        GGVdec = self.GGVdec
        GGVfull = self.GGVfull
        axacc = self.axacc 
        vxvect = self.vxvect
        
        plt.plot(vxvect, axacc)

    def plotAccEnvExtra(self):
        f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4,
                                               figsize=(self.size*1.5,
                                                        self.size/2))
        ax1.set_title("Forces[N] (vcar[m/s])")
        ax1.plot(self.Fzaero, self.vxvect, 'c-', label="Fzaero")
        ax1.plot(self.Fxaero, self.vxvect, 'm-', label="Fxaero")
        ax1.plot(self.Fxgrip, self.vxvect, 'r-', label="FxGrip")
        ax1.plot(self.Fxdrive, self.vxvect, 'g-', label="FxDrive")
        ax1.legend()
        ax1.grid(b=True, which='major', linestyle=':')

        ax2.set_title("Gear (vcar[m/s])")
        ax2.plot(self.nGear, self.vxvect, 'c-', label="nGear")
        ax2.legend()
        ax2.grid(b=True, which='major', linestyle=':')

        ax3.set_title("EngNm (vcar[m/s])")
        ax3.plot(self.EngNm, self.vxvect, 'c-', label="EngNm")
        ax3.legend()
        ax3.grid(b=True, which='major', linestyle=':')

        ax4.set_title("EngRpm (vcar[m/s])")
        ax4.plot(self.EngRpm, self.vxvect, 'c-', label="EngRpm")
        ax4.legend()
        ax4.grid(b=True, which='major', linestyle=':')

    def plotLapTimeSim(self):
        plt.figure(2, figsize=(self.size, self.size/2))
        plt.title("OpenLapSim - Lap Time Simulation")
        plt.plot(self.dist, self.vcar, 'b-', linewidth=2, label="vcar")
        plt.xlabel('distance [m]')
        plt.ylabel('velocity [m/s]')
        plt.legend()
        plt.grid(b=True, which='major', linestyle=':')
        # plt.ylim(0, self.vcarmax*1.2)
        plt.ylim(0, 100)
        plt.xlim(0, max(self.dist))
    
    def plotLapTimeSimAxAcc(self):
        plt.figure(4, figsize=(self.size, self.size/2))
        plt.title("vxvect vs axacc and axdec")
        plt.plot(self.vxvect, self.axacc, 'r-', linewidth=2, label="axacc",)
        plt.plot(self.vxvect, self.axdec, 'g-', linewidth=2, label="axdec",)
        plt.ylim(-50, 20)
        # plt.xlabel('distance [m]')
        # plt.ylabel('vxacc [m/s]')
        plt.legend()
        plt.grid(b=True, which='major', linestyle=':')
        # plt.ylim(0, self.vcarmax*1.2)
        # plt.xlim(0, max(self.dist))

        plt.figure(6, figsize=(self.size, self.size/2))
        plt.title("vxvect vs ay")
        plt.plot(self.vxvect, self.ay, 'r-', linewidth=2, label="ay",)
        # plt.plot(self.vxvect, self.ay, 'g-', linewidth=2, label="axdec",)
        # plt.xlabel('distance [m]')
        # plt.ylabel('vxacc [m/s]')
        plt.legend()
        plt.grid(b=True, which='major', linestyle=':')
        # plt.ylim(0, self.vcarmax*1.2)
        # plt.xlim(0, max(self.dist))

    def plotLapTimeSimExtra(self):
        plt.figure(3, (self.size, self.size/2))
        plt.title("Lap Time Simulation - Extra")
        plt.plot(self.dist, self.vxcor, 'c-', label="vxcor") 
        plt.plot(self.dist, self.vxacc, 'm-', label="vxacc")
        plt.plot(self.dist, self.vxdec, 'r-', label="vxdec")
        plt.plot(self.dist, self.vcar, 'b-', linewidth=2, label="vcar")
        plt.xlabel('distance [m]')
        plt.ylabel('velocity [m/s]')
        plt.legend()
        plt.grid(b=True, which='major', linestyle=':')
        plt.ylim(0, self.vcarmax*1.2)
        plt.xlim(0, max(self.dist))

    def printData(self):
        print("PostProc completed")
        print("---------------------------")
        print("LapTime: ", self.laptime, "[s]")
        print("TopSpeed: ", np.round(self.vcarmax*3.6, 1), "[Km/h]")
        print("---------------------------")

    def plotTrackMap(self, delta_s=10, cmap='viridis', linewidth=2):
            """
            Plot the track layout with color blending based on vehicle velocity.
            Requires self.curv (curvature in 1/m), self.dist, and self.vcar arrays.
            """

            track = np.loadtxt("trackFiles/trackFile.txt")
            dist = track[:, 0]
            curv = track[:, 1]

            # if self.curv is None:
            #     raise ValueError("Curvature data (self.curv) not available to plot track map.")

            # Initialize position arrays
            n = len(self.dist)
            x = np.zeros(n)
            y = np.zeros(n)
            theta = 0.0

            # Integrate heading and position
            for i in range(1, n):
                k = curv[i]
                ds = delta_s
                dtheta = k * ds
                theta_prev = theta
                theta += dtheta
                if k != 0:
                    R = 1.0 / k
                    x[i] = x[i-1] + R*(np.sin(theta) - np.sin(theta_prev))
                    y[i] = y[i-1] + R*(np.cos(theta_prev) - np.cos(theta))
                else:
                    # Treat as straight
                    x[i] = x[i-1] + ds * np.cos(theta_prev)
                    y[i] = y[i-1] + ds * np.sin(theta_prev)

            # Create line segments between points
            points = np.vstack([x, y]).T
            segments = np.stack([points[:-1], points[1:]], axis=1)

            # Compute midpoint velocities and normalize for color
            v_mid = 0.5 * (self.vcar[:-1] + self.vcar[1:])
            norm = (v_mid - v_mid.min()) / (v_mid.max() - v_mid.min())

            # Create a LineCollection with the normalized velocity map
            lc =  LineCollection(segments, cmap=cmap, norm=plt.Normalize(0, 1), linewidths=linewidth)
            lc.set_array(norm)

            # Plot
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.add_collection(lc)
            ax.set_aspect('equal', 'datalim')
            ax.autoscale()
            cbar = plt.colorbar(lc, ax=ax)
            cbar.set_label('Velocity [m/s]')
            plt.title('Track Map - Colored by Vehicle Velocity')
            plt.xlabel('X [m]')
            plt.ylabel('Y [m]')
            plt.show()
    
  

