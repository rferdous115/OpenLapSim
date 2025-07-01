"""
---------------------------
Acceleration Envelope Calculator - OLS
---------------------------

This  class computes the Performance Envelope (Ax, Ay, vcar).

---------------------------
@autor: Davide Strassera
@first release: 2019-12-21
by Python 3.7
---------------------------

"""

# Import Packages
import numpy as np
import scipy.constants as sc
import pandas as pd

class AccEnvCalc:

    def __init__(self, setupDict):
        # inputs
        self.mcar = setupDict["mcar"]
        self.clt = setupDict["clt"]
        self.cx = setupDict["cx"]
        self.afrcar = setupDict["afrcar"]
        self.mbrk = setupDict["mbrk"]
        self.gripx = setupDict["gripx"]
        self.gripy = setupDict["gripy"]
        self.loadEff = setupDict["loadEff"]
        self.rtyre = setupDict["rtyre"]
        self.rGearRat = setupDict["rGearRat"]
        self.reff = setupDict["reff"]
        self.EngNm = setupDict["EngNm"]
        self.EngRpm = setupDict["EngRpm"]
        self.rho = setupDict["rho"]

        # constants
        self.g = sc.g     # 9.80665
        self.pi = sc.pi    # 3.14159
        # parameters
        self.nSteps = 10
        self.nAx = 20
        self.LOAD_EFF_SCALE = 10000  # [N]
        # outputs
        self.accEnvDict = {
            "vxvect": None,  # longitudinal velocity vector? 
            "axacc": None, # longitudinal acceleration vector? 
            "axdec": None, # longitudinal deceleration vector? 
            "ay": None, # lateral acceleration 
            # extra Channels
            "nGear": None, # 
            "EngNm": None,
            "EngRpm": None,
            "Fzaero": None,
            "Fxaero": None,
            "gripx": None,
            "gripy": None,
            "Fxgrip": None,
            "Fxdrive": None,
            # GGV
            "GGVacc": None,
            "GGVdec": None,
            "GGVfull": None,
            }

    def Run(self):

        # Functions Definitions
        def Mfinaldrive(vx, EngNm, EngRpm, rGear):
            neng = np.zeros(len(self.rGearRat))
            meng = np.zeros(len(self.rGearRat))
            Mfinaldrive = np.zeros(len(self.rGearRat))
            for i in range(len(self.rGearRat)): 
                ntyre = vx/(2*self.pi*self.rtyre)*60 # calculate wheel RPM
                neng[i] = ntyre*self.rGearRat[i] # calculate engine RPM at given gear ratio.
                meng[i] = np.interp(neng[i], self.EngRpm, self.EngNm) # calculate engine torque at given engine RPM. 
                # check nengine is in range of rpm max(revlimit) and min(stall)
                if min(EngRpm) < neng[i] < max(EngRpm):
                    meng[i] = meng[i]
                else:
                    meng[i] = 0 # engine torque at given gear is zero. (doesn't run)
                Mfinaldrive[i] = meng[i]*self.rGearRat[i] # final drive torque (torque at the wheels)
            # from the Mfinaldrive[i] array select index which has max Torque
            indexArray = np.where(Mfinaldrive == np.amax(Mfinaldrive))
            # outMfinaldrive = max(Mfinaldrive)
            index = indexArray[0][0]  # WARN: can have 2 pos with same result!
            outMfinaldrive = Mfinaldrive[index]
            outmeng = meng[index] # output the max engine torque 
            outneng = neng[index] # output the max engine rpm 
            outnGear = index+1 # output the the gear ratio which outputs the max engine torque and rpm
            return outMfinaldrive, outmeng, outneng, outnGear

        def Fxaero(vx):
            outFxaero = 0.5*self.rho*pow(vx, 2)*self.afrcar*self.cx  # [N]
            return outFxaero

        def Fzaero(vx):
            outFzaero = 0.5*self.rho*pow(vx, 2)*self.afrcar*self.clt  # [N]
            return outFzaero

        def gripLoadEff(grip, fz):
            deltaGripLoadEff = grip*(self.loadEff*(fz/self.LOAD_EFF_SCALE))
            newGrip = grip - deltaGripLoadEff
            return newGrip

        # VxMax Calculation (forces equilibrium); F_tractive - F_xdrag approaches 0 at what vmax? 
        vxmax = 1
        e = 0.1 
        while e > 0.05: 
            e = ((Mfinaldrive(vxmax, self.EngNm, self.EngRpm, self.rGearRat)[0]
                 / self.rtyre)-Fxaero(vxmax))/self.mcar 
            vxmax += 0.1

        # Ax & Ay Calculation
        vxmax = np.round(vxmax, 1)
        nSteps = self.nSteps
        small = 0.00000001  # to avoid division by zero
        vxvect = np.linspace(small, vxmax, nSteps)

        with open('vxvect_list.txt', 'w') as f:
            for vx in vxvect:
                f.write("%s\n" % vx)

        ay = [0]*len(vxvect)
        for i in range(len(vxvect)):
            Fzaero_ = Fzaero(vxvect[i]) # calculate downforce at that speed
            gripYcurrent = gripLoadEff(self.gripy, Fzaero_+self.mcar*self.g) # calculate lateral grip coefficient based on normal force plus downforce
            ay[i] = (Fzaero_ + self.mcar*self.g) * gripYcurrent / self.mcar # lateral acc ay = (downforce + mg)*gripY/m 

        Fxbrk = self.mbrk/self.rtyre # max braking force
        axacc = [0]*len(vxvect)
        axdec = [0]*len(vxvect)
        Fxgrip = [0]*len(vxvect)
        Fxaero_ = [0]*len(vxvect)
        Fzaero_ = [0]*len(vxvect)
        Fxdrive = [0]*len(vxvect)
        outmeng = [0]*len(vxvect)
        outneng = [0]*len(vxvect)
        outnGear = [0]*len(vxvect)
        # print(f"vxvect len = {len(vxvect)}")s
        for i in range(len(vxvect)):
            Fzaero_[i] = Fzaero(vxvect[i]) # calcualate downforce at that speed. 
            gripXcurrent = gripLoadEff(self.gripx, Fzaero_[i]+self.mcar*self.g) # calculate long. grip coefficient based on normal force plus downforce
            Fxgrip[i] = (Fzaero_[i]+self.mcar*self.g) * gripXcurrent  # calculate tractive grip limit Fx
            Fxaero_[i] = Fxaero(vxvect[i]) # calculate long drag force at that speed. 
            outMfd, outmeng[i], outneng[i], outnGear[i] = Mfinaldrive(vxvect[i],
                                                                      self.EngNm,
                                                                      self.EngRpm,
                                                                      self.rGearRat) # calculate final torque, engine torque, engine RPM, and gear ratio at that speed. 
            Fxdrive[i] = outMfd*self.reff/self.rtyre # tractive force at the wheels
            # print(f"Numerator = {(min(Fxdrive[i], Fxgrip[i])-Fxaero_[i])}")
            axacc[i] = max(0, (min(Fxdrive[i], Fxgrip[i])-Fxaero_[i])/self.mcar) # Calculates net accelerative force acting on the car. F_drive/F_grip - F_aero. min prevents wheel slip. Makes sure actual traction force is less than drivetrain force. 
            axdec[i] = -(min(Fxbrk, Fxgrip[i])+Fxaero_[i])/self.mcar # Calculates net braking force. It's either the max braking force or present tire grip force + drag 
            print(f"{i}: {axdec}")

        self.accEnvDict = {
            "vxvect": vxvect,
            "axacc": axacc,
            "axdec": axdec,
            "ay": ay,
            # extra Channels
            "nGear": outnGear,
            "EngNm": outmeng,
            "EngRpm": outneng,
            "Fzaero": Fzaero_,
            "Fxaero": Fxaero_,
            "gripx": None,
            "gripy": None,
            "Fxgrip": Fxgrip,
            "Fxdrive": Fxdrive,
        }

        def generateGGV(axacc, axdec, ay, vxvect):
            """ This method generates a GGVacc and GGVdec surface using ellipse
            equation for combine, given the vectors (axacc,axdec,ay,vxvect)."""
            nAx = self.nAx #nAx = 20
            nVx = len(vxvect)
            size = nVx*nAx
            # GGV ACCELERATION
            GGVacc = np.zeros((size, 3))  # ay,ax,vx
            for i in range(nVx):
                ayStep = np.absolute(ay[i])/(nAx)
                for j in range(nAx):
                    ayreal = ay[i] - ayStep*j
                    axcombine = np.sqrt(np.absolute(np.power(axacc[i], 2) *
                                (1-(np.power(ayreal, 2)/np.power(ay[i], 2)))))
                    index = (i*nAx)+j
                    GGVacc[index, 0] = np.round(axcombine, 2)
                    GGVacc[index, 1] = np.round(ayreal, 2)
                    GGVacc[index, 2] = np.round(vxvect[i], 2)
            # GGV DECELERATION
            GGVdec = np.zeros((size, 3))  # ax,ay,vx
            for i in range(nVx):
                ayStep = np.absolute(ay[i])/(nAx)
                for j in range(nAx):
                    ayreal = ay[i] - ayStep*j
                    axcombine = - np.sqrt(np.absolute(np.power(axdec[i], 2) *
                                  (1-(np.power(ayreal, 2)/np.power(ay[i], 2)))))
                    index = ((i*nAx)+j)
                    GGVdec[index, 0] = np.round(axcombine, 2)
                    GGVdec[index, 1] = np.round(ayreal, 2)
                    GGVdec[index, 2] = np.round(vxvect[i], 2)
            return GGVacc, GGVdec

        ay = self.accEnvDict["ay"]
        axdec = self.accEnvDict["axdec"]
        axacc = self.accEnvDict["axacc"]
        vxvect = self.accEnvDict["vxvect"]

        GGVacc, GGVdec = generateGGV(axacc, axdec, ay, vxvect)

        # Mirror the GGV to left and concat GGVacc and GGVdec
        GGVaccLeft = GGVacc*[1, -1, 1]
        GGVacc = np.concatenate((GGVacc, GGVaccLeft))
        GGVdecLeft = GGVdec*[1, -1, 1]
        GGVdec = np.concatenate((GGVdec, GGVdecLeft))
        GGVfull = np.concatenate((GGVacc, GGVdec))

        self.accEnvDict["GGVacc"] = GGVacc
        self.accEnvDict["GGVdec"] = GGVdec
        self.accEnvDict["GGVfull"] = GGVfull

        print("AccEnvCalc completed")
