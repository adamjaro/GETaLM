
#_____________________________________________________________________________
# Electron beam-gas using Lifshitz QED, Eq. 93.16
#
#_____________________________________________________________________________

from ctypes import c_double
import atexit
import numpy as np

import ROOT as rt
from ROOT import TRandom3, gROOT, addressof, TDatabasePDG, TLorentzVector
from ROOT import TMath, TF2, TF1, TRandom3, gRandom, TH1D
from ROOT import TFoam, std

from particle import particle
from beam import beam

#_____________________________________________________________________________
class gen_beam_gas:
    #_____________________________________________________________________________
    def __init__(self, parse, tree=None):
        self.parse = parse

        #electron energy, GeV
        self.Ee = parse.getfloat("main", "Ee")

        print("Ee, GeV =", self.Ee)

        #Z of the nucleus
        self.Z = 1
        if parse.has_option("main", "Z"):
            self.Z = parse.getint("main", "Z")
        print("Z:", self.Z)

        #minimal photon energy, GeV
        self.emin = parse.getfloat("main", "emin")
        print("emin, GeV =", self.emin)

        #alpha r_e^2
        self.ar2 = 7.297*2.818*2.818*1e-2 # m barn

        #electron mass
        self.me = TDatabasePDG.Instance().GetParticle(11).Mass()

        #maximal delta
        self.dmax = 200.
        if parse.has_option("main", "dmax"):
            self.dmax = parse.getfloat("main", "dmax")

        print("dmax:", self.dmax)

        #cross section formula
        #self.eqpar = self.eq_93p16(self)
        #self.dSigDwDt = TF2("dSigDwDt", self.eqpar, self.emin, self.Ee, 0, self.dmax)
        #self.dSigDwDt.SetNpx(2000)
        #self.dSigDwDt.SetNpy(2000)
        #self.dSigDwDt.SetNpx(100)
        #self.dSigDwDt.SetNpy(100)
        #gRandom.SetSeed(5572323)

        #total integrated cross section over all delta (to 1e5)
        #dSigInt = TF2("dSigInt", self.eqpar, self.emin, self.Ee, 0, 1e5)
        #sigma_tot = dSigInt.Integral(self.emin, self.Ee, 0, 1e5)
        #print("Total cross section, mb:", sigma_tot)

        #FOAM integration
        self.eq_foam = self.eq_93p16_foam(self)
        self.rand_foam = TRandom3()
        self.rand_foam.SetSeed(123)
        self.foam = TFoam("foam")
        self.foam.SetkDim(2)
        self.foam.SetnCells(2000)
        self.foam.SetRhoInt(self.eq_foam)
        self.foam.SetPseRan(self.rand_foam)
        self.foam.Initialize()

        #to print the total cross section from FOAM at the end
        atexit.register(self.finish)

        #uniform generator for azimuthal angles
        self.rand = TRandom3()
        self.rand.SetSeed(5572323)

        #chamber pressure for z-vertex
        self.pressure_par = self.eq_pressure(self)
        zmin = self.pressure_par.zmin
        zmax = self.pressure_par.zmax
        #custom range in z, only within the default limits given by pressure data
        if parse.has_option("main", "zmin") and ( parse.getfloat("main", "zmin") > zmin ):
            zmin = parse.getfloat("main", "zmin")
        if parse.has_option("main", "zmax") and ( parse.getfloat("main", "zmax") < zmax ):
            zmax = parse.getfloat("main", "zmax")
        print("zmin, zmax (mm):", zmin, zmax)
        self.pressure_func = TF1("pressure", self.pressure_par, zmin, zmax)

        #electron lattice
        self.lat = self.load_lattice()

        #beam transverte shape for vertex in x and y
        self.beam_par = self.eq_beam_sigma(self.lat, self.pressure_par.hz, parse)

        #tree output from the generator
        tlist = ["true_phot_w", "true_phot_delta"]
        tlist += ["true_phot_theta", "true_phot_phi", "true_phot_E"]
        tlist += ["vtx_x", "vtx_y", "vtx_z", "divx", "divy"]
        self.tree_out = self.set_tree(tree, tlist)

        print("Beam-gas generator initialized")

    #_____________________________________________________________________________
    def generate(self, add_particle):

        #return

        #photon energy and delta in nucleus rest frame
        #w = c_double(0)
        #d = c_double(0)
        #self.dSigDwDt.GetRandom2(w, d)

        #w = w.value
        #d = d.value

        #photon energy and delta in nucleus rest frame by FOAM
        wdvec = std.vector("double")(2)
        self.foam.MakeEvent()
        self.foam.GetMCvect(wdvec.data())

        w = wdvec[0]*self.Ee
        d = wdvec[1]*self.dmax

        #set the tree output
        self.tree_out.true_phot_w = w
        self.tree_out.true_phot_delta = d

        #polar angle theta
        theta = d*self.me/self.Ee

        #uniform azimuthal angle
        phi = 2. * TMath.Pi() * self.rand.Rndm()

        #photon
        phot = add_particle( particle(22) )
        phot.stat = 1
        phot.pxyze_prec = 9

        #photon vector, negative pz
        px = w*TMath.Sin(theta)*TMath.Cos(phi)
        py = w*TMath.Sin(theta)*TMath.Sin(phi)
        pz = w*TMath.Cos(theta)

        phot.vec.SetPxPyPzE(px, py, -pz, w)

        #photon kinematics in generator output
        self.tree_out.true_phot_theta = phot.vec.Theta()
        self.tree_out.true_phot_phi = phot.vec.Phi()
        self.tree_out.true_phot_E = phot.vec.E()

        #scattered electron, initialize as beam and constrain with the photon
        electron = add_particle( beam(self.Ee, 11, -1) )
        electron.stat = 1
        electron.pxyze_prec = 9

        electron.vec -= phot.vec

        #z-vertex from the pressure
        zpos = self.pressure_func.GetRandom()
        phot.vz = zpos
        electron.vz = zpos

        #xy vertex and divergence from lattice
        xpos, ypos, divx, divy = self.beam_par(zpos)
        phot.vx = xpos
        phot.vy = ypos
        electron.vx = xpos
        electron.vy = ypos

        #vertex position and divergence in output tree
        self.tree_out.vtx_x = xpos
        self.tree_out.vtx_y = ypos
        self.tree_out.vtx_z = zpos
        self.tree_out.divx = divx
        self.tree_out.divy = divy

        #divergence in x by rotation along y
        phot.vec.RotateY(divx)
        electron.vec.RotateY(divx)

        #divergence in y by rotation along x
        phot.vec.RotateX(divy)
        electron.vec.RotateX(divy)

    #_____________________________________________________________________________
    class eq_93p16:
        def __init__(self, gen):
            self.gen = gen
        def __call__(self, x, par):

            #Eq. 93.16 from Lifshitz QED for bremsstrahlung

            #photon energy and angle
            w = x[0]
            d = x[1]

            #initial and final electron E and E'
            E = self.gen.Ee
            Efin = E - w

            t1 = 8.*self.gen.Z*self.gen.Z*self.gen.ar2*(1./w)*(Efin/E)*d/((1+d**2)**2)

            t2 = ( (E/Efin) + (Efin/E) - 4*d*d/((1+d**2)**2) )*TMath.Log(2.*E*Efin/(self.gen.me*w))

            t3 = 0.5*( (E/Efin) + (Efin/E) + 2 - 16*d*d/((1+d**2)**2) )

            sig = t1*(t2 - t3)

            return sig

    #_____________________________________________________________________________
    class eq_93p16_foam:
        def __init__(self, gen):
            self.gen = gen
        def __call__(self, ndim, x):

            #Eq. 93.16 from Lifshitz QED for bremsstrahlung, version for TFoam

            #photon energy and angle, FOAM input ranges from 0 to 1
            w = x[0]*self.gen.Ee
            d = x[1]*self.gen.dmax

            if w < self.gen.emin: return 0.

            #initial and final electron E and E'
            E = self.gen.Ee
            Efin = E - w

            t1 = 8.*self.gen.Z*self.gen.Z*self.gen.ar2*(1./w)*(Efin/E)*d/((1+d**2)**2)

            t2 = ( (E/Efin) + (Efin/E) - 4*d*d/((1+d**2)**2) )*TMath.Log(2.*E*Efin/(self.gen.me*w))

            t3 = 0.5*( (E/Efin) + (Efin/E) + 2 - 16*d*d/((1+d**2)**2) )

            sig = t1*(t2 - t3)

            return sig

    #_____________________________________________________________________________
    def finish(self):

        foam_int = c_double(0)
        foam_int_err = c_double(0)
        self.foam.Finalize(foam_int, foam_int_err)

        sig = foam_int.value*self.dmax*(self.Ee-self.emin)
        sig_err = foam_int_err.value*self.dmax*(self.Ee-self.emin)

        print("Total cross section (mb):", sig, "+/-", sig_err)

    #_____________________________________________________________________________
    class eq_pressure:
        def __init__(self, gen):

            #pressure from input spreadsheet for z-vertex

            import pandas as pd
            from scipy.interpolate import interp1d

            #open the input
            self.xls = pd.read_excel(gen.parse.get("main", "pressure_xlsx").strip("\"'"),\
                sheet_name=gen.parse.get("main", "pressure_sheet").strip("\"'"),\
                usecols=gen.parse.get("main", "pressure_usecols").strip("\"'"),\
                skiprows=gen.parse.getint("main", "pressure_skiprows"),\
                nrows=gen.parse.getint("main", "pressure_nrows"),\
                index_col=None, header=None, engine="openpyxl")

            #z position in mm and inverted sign for detector convention
            self.xls[1] = -1e3*self.xls[1]

            #additional flat vacuum level at lower z below the xlsx if requested
            if gen.parse.has_option("main", "pressure_lowz_start"):
                lowz_start = gen.parse.getfloat("main", "pressure_lowz_start")
                lowz_end = gen.parse.getfloat("main", "pressure_lowz_end")
                lowz_mbar = gen.parse.getfloat("main", "pressure_lowz_mbar")

                lowz_data = [{1:lowz_start, 2:lowz_mbar}, {1:lowz_end, 2:lowz_mbar}]

                #add the additional pressure level to the DataFrame
                self.xls = pd.concat( [pd.DataFrame(lowz_data), self.xls], ignore_index=True )

            #range in z from the input
            self.zmin = self.xls[1][ self.xls[1].index[0] ]
            self.zmax = self.xls[1][ self.xls[1].index[-1] ]

            #print("zmin, zmax:", self.zmin, self.zmax)

            #linear interpolation from the input data
            self.interp = interp1d(self.xls[1], self.xls[2], kind="linear")

            #spacing in z for vertices in x and y
            beam_z_spacing = 400
            if gen.parse.has_option("main", "beam_z_spacing"):
                beam_z_spacing = gen.parse.getint("main", "beam_z_spacing")

            print("beam_z_spacing:", beam_z_spacing)
            self.hz = TH1D("pressure_hz", "pressure_hz", beam_z_spacing, self.zmin, self.zmax)

        def __call__(self, x, par):

            return self.interp(x[0])

    #_____________________________________________________________________________
    class eq_beam_sigma:
        def __init__(self, lat, hz, parse):

            #beam sigma from input lattice for vertex in x and y

            #beam emittance
            eps_x = parse.getfloat("main", "eps_x")
            eps_y = parse.getfloat("main", "eps_y")

            #spacing in z
            self.hz = hz

            #interpolation for beta(s) given in lattice convention
            from scipy.interpolate import CubicHermiteSpline, interp1d
            beta_x = CubicHermiteSpline(lat["s"], lat["beta_x"], -2*lat["alpha_x"])
            beta_y = CubicHermiteSpline(lat["s"], lat["beta_y"], -2*lat["alpha_y"])

            #print("bins:", self.hz.GetNbinsX(), self.hz.GetBinWidth(0))

            #interpolation for angular divergence as a function of s in lattice convention
            idiv_x = interp1d(lat["s"], self.beam_get_divergence(eps_x, lat["alpha_x"], lat["beta_x"]), kind="linear")
            idiv_y = interp1d(lat["s"], self.beam_get_divergence(eps_y, lat["alpha_y"], lat["beta_y"]), kind="linear")

            #interpolation for xpos
            f=interpolate.interp1d(lat["orbit_z"], lat["orbit_x"]*1e3)

            #horizontal and vertical Gaussians
            self.gx = {}
            self.gy = {}

            #angular divergence
            self.divx = {}
            self.divy = {}

            for ibin in range(self.hz.GetNbinsX()+1):

                #beam sigma in mm at a given z, z for beta(s) in meters in lattice convention
                spos = -1e-3*self.hz.GetBinCenter(ibin)
                sigma_x = 1e3*np.sqrt(eps_x*beta_x(spos))
                sigma_y = 1e3*np.sqrt(eps_y*beta_y(spos))
                
                #value of x(zpos) after interpolation. Minus sign in order to be consistence with detector convention of z and x axis. 
                xmean=-f(spos)

                #Gaussians for a given sigma at a given z
                self.gx[ibin] = TF1("beam_sigma_x_"+str(ibin), "gaus", xmean-5.*sigma_x, xmean+5.*sigma_x)
                self.gx[ibin].SetParameters(1, xmean, sigma_x)
                self.gy[ibin] = TF1("beam_sigma_y_"+str(ibin), "gaus", -5.*sigma_y, 5.*sigma_y)
                self.gy[ibin].SetParameters(1, 0, sigma_y)

                #print(ibin, self.hz.GetBinCenter(ibin), sigma_x, sigma_y)

                #divergence Gaussians
                dx = idiv_x(spos)
                dy = idiv_y(spos)
                self.divx[ibin] = TF1("beam_div_x_"+str(ibin), "gaus", -5.*dx, 5.*dx)
                self.divx[ibin].SetParameters(1, 0, dx)
                self.divy[ibin] = TF1("beam_div_y_"+str(ibin), "gaus", -5.*dy, 5.*dy)
                self.divy[ibin].SetParameters(1, 0, dy)

        def __call__(self, zpos):

            #Gaussians at a given z
            ibin = self.hz.FindBin(zpos)
            if ibin < 0:
                ibin = 0
            if ibin > self.hz.GetNbinsX():
                ibin = self.hz.GetNbinsX()

            #generate horizontal and vertical vertex position
            xpos = self.gx[ibin].GetRandom()
            ypos = self.gy[ibin].GetRandom()

            #angular divergence
            dx = self.divx[ibin].GetRandom()
            dy = self.divy[ibin].GetRandom()

            return xpos, ypos, dx, dy

        def beam_get_divergence(self, eps, alpha, beta):

            #angular divergence in rad

            return np.sqrt( eps*( (1. + alpha**2)/beta ) )

    #_____________________________________________________________________________
    def load_lattice(self):

        from pandas import DataFrame

        #input lattice
        inp = open(self.parse.get("main", "lattice_txt").strip("\"'"), "r")

        #lattice dataframe
        col = ["name", "key", "s", "length", "angle", "x_pitch", "magnet_x", "magnet_z", "magnet_theta",\
            "orbit_x", "orbit_z", "orbit_theta", "dispersion", "dispersion_derivative", "beta_x", "alpha_x",\
            "beta_y", "alpha_y", "field", "gradient"]
        val = []

        #skip file header
        for i in range(3): inp.readline()

        #load the lattice
        while(True):
            line = inp.readline()
            if(len(line) == 0): break

            lin = line.split()

            #IP6 marker is identical to the drift before
            if lin[0] == "IP6": continue

            #name and key as string, values as float
            for i in range(len(lin)):
                if i < 2:
                    lin[i] = str(lin[i])
                else:
                    lin[i] = float(lin[i])

            val.append( lin )

        df = DataFrame(val, columns=col)

        #range with pressure data, in lattice convention
        #df = df.query("s>-33 and s<10")
        df = df.query("s>-33 and s<40")

        return df

    #_____________________________________________________________________________
    def set_tree(self, tree, tlist):

        #set output to the tree

        #tree variables
        struct = "struct gen_beam_gas { Double_t "
        for i in tlist:
            struct += i + ", "
        struct = struct[:-2] + ";};"
        gROOT.ProcessLine( struct )
        tree_out = rt.gen_beam_gas()

        #put zero to all variables
        for i in tlist:
            exec("tree_out."+i+"=0")

        #add variables to the tree
        if tree is not None:
            for i in tlist:
                tree.Branch(i, addressof(tree_out, i), i+"/D")

        return tree_out


















