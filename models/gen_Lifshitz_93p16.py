
#_____________________________________________________________________________
# Bremsstrahlung photons by Lifshitz QED, Eq. 93.16
#
#_____________________________________________________________________________

from ctypes import c_double

import ROOT as rt
from ROOT import TRandom3, gROOT, AddressOf, TDatabasePDG, TLorentzVector
from ROOT import TMath, TF2, Double, TRandom3, gRandom

from particle import particle
from beam import beam

#_____________________________________________________________________________
class gen_Lifshitz_93p16:
    #_____________________________________________________________________________
    def __init__(self, parse, tree):

        #electron and proton energy, GeV
        self.Ee = parse.getfloat("main", "Ee")
        self.Ep = parse.getfloat("main", "Ep")

        print "Ee, GeV =", self.Ee
        print "Ep, GeV =", self.Ep

        #A and Z of the nucleus
        self.A = 1
        self.Z = 1
        if parse.has_option("main", "A"):
            self.A = parse.getint("main", "A")
        if parse.has_option("main", "Z"):
            self.Z = parse.getint("main", "Z")
        print "A:", self.A
        print "Z:", self.Z

        #minimal photon energy, GeV
        self.emin = parse.getfloat("main", "emin")
        print "emin, GeV =", self.emin

        #alpha r_e^2
        self.ar2 = 7.297*2.818*2.818*1e-2 # m barn

        #electron and nucleus mass
        self.me = TDatabasePDG.Instance().GetParticle(11).Mass()
        self.mp = TDatabasePDG.Instance().GetParticle(2212).Mass()
        self.mn = self.A*self.mp

        #nucleus beam vector
        nvec = TLorentzVector()
        pz_a = TMath.Sqrt(self.Ep**2-self.mp**2)*self.Z
        en_a = TMath.Sqrt(pz_a**2 + self.mn**2)
        nvec.SetPxPyPzE(0, 0, pz_a, en_a)

        #boost vector of nucleus beam
        self.nbvec = nvec.BoostVector()

        #electron beam vector
        evec = TLorentzVector()
        evec.SetPxPyPzE(0, 0, -TMath.Sqrt(self.Ee**2-self.me**2), self.Ee)

        #electron beam energy in nucleus beam rest frame
        evec.Boost(-self.nbvec.x(), -self.nbvec.y(), -self.nbvec.z())
        self.Ee_n = evec.E()

        print "Ee_n, GeV:", self.Ee_n

        #minimal photon energy in nucleus rest frame
        eminv = TLorentzVector()
        eminv.SetPxPyPzE(0, 0, -self.emin, self.emin)
        eminv.Boost(-self.nbvec.x(), -self.nbvec.y(), -self.nbvec.z())
        emin_n = eminv.E()
        print "emin_n, GeV:", emin_n

        #maximal delta in nucleus frame
        dmax_n = 100.
        if parse.has_option("main", "dmax_n"):
            dmax_n = parse.getfloat("main", "dmax_n")

        print "dmax_n:", dmax_n

        #cross section formula
        self.dSigDwDt = TF2("dSigDwDt", self.eq, emin_n, self.Ee_n, 0, dmax_n)
        self.dSigDwDt.SetNpx(2000)
        self.dSigDwDt.SetNpy(2000)
        gRandom.SetSeed(5572323)

        #total integrated cross section over all delta (to 1e5)
        dSigInt = TF2("dSigInt", self.eq, emin_n, self.Ee_n, 0, 1e5)
        sigma_tot = dSigInt.Integral(emin_n, self.Ee_n, 0, 1e5)

        print "Total cross section, mb:", sigma_tot

        #uniform generator for azimuthal angles
        self.rand = TRandom3()
        self.rand.SetSeed(5572323)

        #tree output from the generator
        tlist = ["true_phot_w", "true_phot_delta", "true_phot_theta_n"]
        tlist += ["true_phot_theta", "true_phot_phi", "true_phot_E"]
        tlist += ["true_el_theta", "true_el_phi", "true_el_E"]
        self.tree_out = self.set_tree(tree, tlist)

        print "Lifshitz_93p16 parametrization initialized"

    #_____________________________________________________________________________
    def generate(self, add_particle):

        #photon energy and delta in nucleus rest frame
        #w = Double(0)
        #d = Double(0)
        w = c_double(0)
        d = c_double(0)
        self.dSigDwDt.GetRandom2(w, d)

        w = w.value
        d = d.value

        #polar angle theta
        theta_n = d*self.me/self.Ee_n

        #set the tree output
        self.tree_out.true_phot_w = w
        self.tree_out.true_phot_delta = d
        self.tree_out.true_phot_theta_n = theta_n

        #uniform azimuthal angle
        phi_n = 2. * TMath.Pi() * self.rand.Rndm() #uniform azimuthal angle

        #photon
        phot = add_particle( particle(22) )
        phot.stat = 1
        phot.pxyze_prec = 9

        #set the photon vector in nucleus rest frame
        px_n = w*TMath.Sin(theta_n)*TMath.Cos(phi_n)
        py_n = w*TMath.Sin(theta_n)*TMath.Sin(phi_n)
        pz_n = w*TMath.Cos(theta_n)

        phot.vec.SetPxPyPzE(px_n, py_n, pz_n, w)

        #transform the photon vector to laboratory frame
        phot.vec.Boost(-self.nbvec.x(), -self.nbvec.y(), -self.nbvec.z())

        #rotate the photon for pz < 0
        phot.vec.SetPxPyPzE(phot.vec.Px(), phot.vec.Py(), -phot.vec.Pz(), phot.vec.E())

        #photon kinematics in generator output
        self.tree_out.true_phot_theta = phot.vec.Theta()
        self.tree_out.true_phot_phi = phot.vec.Phi()
        self.tree_out.true_phot_E = phot.vec.E()

        #scattered electron, initialize as beam
        electron = add_particle( beam(self.Ee, 11, -1) )
        electron.stat = 1
        electron.pxyze_prec = 9

        #constrain the scattered electron with the photon
        electron.vec -= phot.vec

        #electron kinematics in generator output
        self.tree_out.true_el_theta = electron.vec.Theta()
        self.tree_out.true_el_phi = electron.vec.Phi()
        self.tree_out.true_el_E = electron.vec.E()

    #_____________________________________________________________________________
    def eq(self, x):

        #Eq. 93.16

        #photon energy and angle
        w = x[0]
        d = x[1]

        #initial and final electron E and E'
        E = self.Ee_n
        Efin = E - w

        t1 = 8.*self.Z*self.Z*self.ar2*(1./w)*(Efin/E)*d/((1+d**2)**2)

        t2 = ( (E/Efin) + (Efin/E) - 4*d*d/((1+d**2)**2) )*TMath.Log(2.*E*Efin/(self.me*w))

        t3 = 0.5*( (E/Efin) + (Efin/E) + 2 - 16*d*d/((1+d**2)**2) )

        sig = t1*(t2 - t3)

        return sig

    #_____________________________________________________________________________
    def set_tree(self, tree, tlist):

        #set output to the tree

        #tree variables
        struct = "struct gen_Lifshitz_93p16 { Double_t "
        for i in tlist:
            struct += i + ", "
        struct = struct[:-2] + ";};"
        gROOT.ProcessLine( struct )
        tree_out = rt.gen_Lifshitz_93p16()

        #put zero to all variables
        for i in tlist:
            exec("tree_out."+i+"=0")

        #add variables to the tree
        if tree is not None:
            for i in tlist:
                tree.Branch(i, AddressOf(tree_out, i), i+"/D")

        return tree_out








