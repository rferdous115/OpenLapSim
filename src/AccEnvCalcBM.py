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
from casadi import *


class AccEnvCalcBM:

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

        # Bicycle Model Variables
        # self.wheelbase = setupDict["wheelbase"]
        # self.cogx = setupDict["cogx"]
        # self.cogz = setupDict["cogz"]

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

        LOAD_EFF_SCALE = 10000

        # aE.mcar = 740
        cogx = 1.78
        cogz = 0.31
        # aE.g = 9.81
        wheelbase = 3.12

        W_R = self.mcar * self.g * cogx / wheelbase # long. weight transfer rear 
        W_F = self.mcar * self.g - W_R # long. weight transfer front


        Fxbrk = self.mbrk/self.rtyre # max braking force
        axacc = [0]*len(vxvect)
        axacc_bm = [0]*len(vxvect)
        axdec = [0]*len(vxvect)
        Fxgrip = [0]*len(vxvect)
        Fxgrip_R = [0]*len(vxvect)
        Fxgrip_F = [0]*len(vxvect)
        Fxaero_ = [0]*len(vxvect)
        Fzaero_ = [0]*len(vxvect)
        Fxdrive = [0]*len(vxvect)
        loadTransfer_ = [0]*len(vxvect)
        outmeng = [0]*len(vxvect)
        outneng = [0]*len(vxvect)
        outnGear = [0]*len(vxvect)

        # Ay lateral acceleration
        ay = [0]*len(vxvect)
        for i in range(len(vxvect)):
            Fzaero_[i] = Fzaero(vxvect[i]) # calculate downforce at that speed
            gripYcurrent = gripLoadEff(self.gripy, Fzaero_[i]+self.mcar*self.g) # calculate lateral grip coefficient based on normal force plus downforce
            ay[i] = (Fzaero_[i] + self.mcar*self.g) * gripYcurrent / self.mcar # lateral acc ay = (downforce + mg)*gripY/m

            # front and rear down-loads
            Fz_F = Fzaero_[i] + W_F
            Fz_R = Fzaero_[i] + W_R

            # front/rear grip limits
            Fx_grip_F = Fz_F * self.gripx
            Fx_grip_R = Fz_R * self.gripx

            ay[i] = (Fz_F + Fz_R) * gripYcurrent / self.mcar # lateral acc ay = (downforce + mg)*gripY/m

        #Ax longitudinal acceleration

        for i in range(len(vxvect)):

            LOAD_EFF_SCALE = 10000

            # aE.mcar = 740
            cogx = 1.78
            cogz = 0.31
            # aE.g = 9.81
            wheelbase = 3.12
            # aE.gripx = 1.15

            # deltaN = self.mcar * axacc[i] * cogz / (self.g * wheelbase)
            Fzaero_[i] = Fzaero(vxvect[i]) # calcualate downforce at that speed. 
            # loadTransfer_[i] = self.mcar

            # Bicycle Model
            # gripXcurrent_F = gripLoadEff(self.gripx, Fzaero_[i]+ W_F - deltaN)
            # gripXcurrent_R = gripLoadEff(self.gripx, Fzaero_[i]+ W_F + deltaN)
            # Fxgrip_F[i] = (Fzaero_[i] + W_F - deltaN) * gripXcurrent_F
            # Fxgrip_R[i] = (Fzaero_[i] + W_R + deltaN) * gripXcurrent_R


            gripXcurrent = gripLoadEff(self.gripx, Fzaero_[i]+self.mcar*self.g) #!!! # calculate long. grip coefficient based on normal force plus downforce
            # gripXcurrent = self.gripx
            Fxgrip[i] = (Fzaero_[i]+self.mcar*self.g) * gripXcurrent  # calculate tractive grip limit Fx
            Fxaero_[i] = Fxaero(vxvect[i]) # calculate long drag force at that speed. 
            outMfd, outmeng[i], outneng[i], outnGear[i] = Mfinaldrive(vxvect[i],
                                                                      self.EngNm,
                                                                      self.EngRpm,
                                                                      self.rGearRat) # calculate final torque, engine torque, engine RPM, and gear ratio at that speed. 
            Fxdrive[i] = outMfd*self.reff/self.rtyre # tractive force at the wheels

            ### START

            #Bicycle Model 

            # FORWARD LONG ACCELERATION

            W_R = self.mcar * self.g * cogx / wheelbase # long. weight transfer rear 
            W_F = self.mcar * self.g - W_R # long. weight transfer front

            # [old] i = 5; Fz_selfro[i] = 4108.378429795827; Fxgrip[i] = 11584.638986282649; Fxselfro_[i] = 1325.2833644502668; Fxdrive[i] = 7338.889980531738; axacc[i] = 8.126495427137122
            # [old] vxvect[25] = 42.806122453877556; Fz_aero[i] = 3465.0005500432576; Fxgrip[i] = 11008.174362709; Fxaero_[i] = 1117.7421129171798; Fxdrive[i] = 8813.373435419537; axacc[i] = 10.399501787165347
            # [trying now] vxvect[47] = 80.47551020448981; Fz_aero[i] = 12246.697941394536; Fxgrip[i] = 18054.66353677579; Fxaero_[i] = 3950.5477230304955; Fxdrive[i] = 4590.467997206904; axacc[i] = 0.8647571272654173
            # [trying now] vxvect[19] = 32.53265306734695; Fz_aero[i] = 2001.3843180002784; Fxgrip[i] = 9661.314616045227; Fxaero_[i] = 645.6078445162188; Fxdrive[i] = 11425.606197737097; axacc[i] = 12.183387529093254
            C_0 = Fxdrive[i] # Fxdrive[i]
            C_1 = Fzaero_[i]*0.5 + W_F # Fz_aero[i] + weight transfer front
            # C_2 = self.mcar * cogz / (self.g * wheelbase)
            C_2 = self.mcar * cogz / (self.g * wheelbase)
            C_3 = Fxaero_[i] # Fxaero[i]
            C_4 = Fzaero_[i]*0.5 + W_R # Fz_aero[i] + weight transfer rear
            C_5 = Fxbrk
            C_6 = 0


            a = MX.sym('a')
            b = MX.sym('b')
            a2 = MX.sym('a2') # acceleration
            g = MX.sym('g') # grip
            F = MX.sym('F')
            x = MX.sym('x')
            x2 = MX.sym('x2')
            y = MX.sym('y')

            # dynamic weight transfer
            deltaN = (self.mcar * a * cogz)/wheelbase

            factor_one = g - ((g * F * self.loadEff) / LOAD_EFF_SCALE)
            # factor_two = min(1.0, max(0.0, g - g * F / LOAD_EFF_SCALE))

            # SKELETON

            # def gripLoadEff(grip, fz):
            # deltaGripLoadEff = grip*(self.loadEff*(fz/self.LOAD_EFF_SCALE))
            # newGrip = grip - deltaGripLoadEff
            # return newGrip


            gripLoadEffFunc = Function('gripLoadEffFunc', [g, F],\
                                    [F, factor_one], \
                                        ['g', 'F'], ['F_in', 'F_loadEff'])

            # gripLoadEffFuncSX = gripLoadEffFunc.expand()

            maxFunc = Function('maxFunc', [x,y],\
                                [y, 0.5*(x + y + sqrt((x-y)**2))], \
                                    ['x', 'y'], ['x_in', 'm'])

            minFunc = Function('minFunc', [x, y], \
                        [y, 0.5*(x + y - sqrt((x-y)**2))], \
                            ['x', 'y'], ['x_in', 'm'])

            # max0Func = Function('max0Func', [x],\
            #             (x, x + sqrt(x**2)) / 2)

            # arg = [gripx, x - y * a]
            # res = gripLoadEffFunc(arg)

            # res = gripLoadEffFunc(self.gripx, x - y * a)
            # F_long = ('F_long', [x, y, a], 
            #             [(x - y * a) * res[1] - C_3])

            resC1_C2 = gripLoadEffFunc(self.gripx, C_1 - C_2 * a)
            resC2_C4 = gripLoadEffFunc(self.gripx, C_4 + C_2 * a)

            # resR = gripLoadEffFunc(self.gripx, C_4 + C_2*a)
            # rear_trac = (C_4 + C_2*a) * resR[1]    # rear grip force    

            # resC1_C2 = self.gripx
            # resC2_C4 = self.gripx

            # resMin01 = minFunc(C_0, (C_1 - C_2 * a) * resC1_C2[1] - C_3)
            # resMin02 = minFunc(C_0, (C_4 + C_2 * a) * resC2_C4[1] - C_3)

            # resMaxFunc01 = maxFunc(0, 0.5*(C_0 + (C_1 - C_2 * a) * resC1_C2[1] - C_3 - sqrt((C_0-(C_1 - C_2 * a)*resC1_C2[1]+C_3)**2+0.0001)))
            resMaxFunc01 = [0,0]
            # resMaxFunc02 = maxFunc(0, 0.5*(C_0 + (C_4 + C_2 * a) * resC2_C4[1] - C_3 - sqrt((C_0-(C_4 + C_2 * a)*resC2_C4[1]+C_3)**2+0.0001)))
            resMaxFunc02 = maxFunc(0, 0.5*(C_0 + (C_4 + C_2 * a) * resC2_C4[1]- sqrt((C_0-(C_4 + C_2 * a)*resC2_C4[1])**2+0.0001)) - C_3)

            # resMaxFunc01 = maxFunc(0, 0.5*(C_0 + (C_1 - C_2 * a) * resC1_C2 - C_3 - sqrt((C_0-(C_1 - C_2 * a)*resC1_C2+C_3)**2))) 
            # resMaxFunc02 = maxFunc(0, 0.5*(C_0 + (C_4 + C_2 * a) * resC2_C4 - C_3 - sqrt((C_0-(C_4 + C_2 * a)*resC2_C4+C_3)**2)))

            # print(f"resMaxFunc01 + resMaxFunc02= {resMaxFunc01[1] + resMaxFunc02[1]}")
            # resMaxFunc01 = maxFunc(0, resMin01)
            # resMaxFunc02 = maxFunc(0, resMin02)

            constraint =  a - ((resMaxFunc01[1] + resMaxFunc02[1]) / self.mcar)
            # print(type(constraint))
                
            nlp = {'x': a, 'f': ((resMaxFunc01[1] + resMaxFunc02[1]) / self.mcar)
                                            , 'g': constraint}

            solver = nlpsol('S', 'ipopt', nlp)
            # print(S)

            sol = solver(x0=0.1, lbg=0, ubg=0)
            axacc[i] = sol['x'].__float__()

            # axacc[i] = max(0, (min(Fxdrive[i], Fxgrip[i])-Fxaero_[i])/self.mcar)

            # LONG DECELERATION

            # resC1a_C2a = gripLoadEffFunc(self.gripx, C_1 - C_2 * a2)
            # resC2a_C4a = gripLoadEffFunc(self.gripx, C_4 + C_2 * a2)

            # resMaxFunc03 = minFunc(0, 0.5*(C_5 + (C_1 - C_2 * a2) * resC1a_C2a[1] - C_3 - sqrt((C_5-(C_1 - C_2 * a2)*resC1a_C2a[1]+C_3)**2)))
            # resMaxFunc04 = minFunc(0, 0.5*(C_5 + (C_4 + C_2 * a2) * resC2a_C4a[1] - C_3 - sqrt((C_5-(C_4 + C_2 * a2)*resC2a_C4a[1]+C_3)**2)))

            # # resMaxFunc03 = 0.5*(C_5 + (C_1 - C_2 * a2) * resC1a_C2a[1] - C_3 - sqrt((C_5-(C_1 - C_2 * a2)*resC1a_C2a[1]+C_3)**2))
            # # resMaxFunc04 = 0.5*(C_5 + (C_4 + C_2 * a2) * resC2a_C4a[1] - C_3 - sqrt((C_5-(C_4 + C_2 * a2)*resC2a_C4a[1]+C_3)**2))

            # constraint2 =  a2 - ((resMaxFunc03[1] + resMaxFunc04[1])/ self.mcar)
            
            # nlp2 = {'x': a2, 'f': (-(resMaxFunc03[1] + resMaxFunc04[1]) / self.mcar)
            #                                 , 'g': constraint2}

            # solver2 = nlpsol('S2', 'ipopt', nlp2)
            # # print(S)

            # sol2 = solver2(x0=-0.1, lbg=0, ubg=0)
            # axdec[i] = sol2['x'].__float__()
            # print(f"{i}: {axdec}")

            
            # front and rear down-loads
            Fz_F = Fzaero_[i]/2 + W_F
            Fz_R = Fzaero_[i]/2 + W_R

            # front/rear grip limits
            Fx_grip_F = Fz_F * self.gripx
            Fx_grip_R = Fz_R * self.gripx

            # distribute the brake force front/rear (simplet: split equally or use ration of normal loads)
            ratio_F = Fz_F /(Fz_F + Fz_R)
            Fx_brake_F = min(Fxbrk*ratio_F, Fx_grip_F)
            Fx_brake_R = min(Fxbrk*(1-ratio_F), Fx_grip_R)

            # total brake 
            Fx_brake = Fx_brake_F + Fx_brake_R
            axdec[i] = -(min(Fx_brake, Fxgrip[i])+Fxaero_[i])/self.mcar

            # axdec[i] = -(min(Fxbrk, Fxgrip[i])+Fxaero_[i])/self.mcar
            
            ### END

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

        print("AccEnvCalcBM completed")
