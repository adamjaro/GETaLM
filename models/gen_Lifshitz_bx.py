
#_____________________________________________________________________________
# Bremsstrahlung photons and electrons generated for individual bunch
# crossings, based on Lifshitz QED, Eq. 93.16
#
#_____________________________________________________________________________

from ctypes import c_double
import atexit

import ROOT
from ROOT import TDatabasePDG, TLorentzVector, TMath
from ROOT import TRandom3, TFoam, std, TF1
from ROOT import gROOT, addressof

from particle import particle
from beam import beam
from beam_effects import beam_effects

#_____________________________________________________________________________
class gen_Lifshitz_bx:
    #_____________________________________________________________________________
    def __init__(self, parse, tree, hepmc_attrib):

        #electron and proton beam energy, GeV
        self.Ee = parse.getfloat("main", "Ee")
        Ep = parse.getfloat("main", "Ep")

        print("Ee (GeV):", self.Ee)
        print("Ep (GeV):", Ep)

        #A and Z of the nucleus
        nA = 1
        self.Z = 1
        if parse.has_option("main", "A"):
            nA = parse.getint("main", "A")
        if parse.has_option("main", "Z"):
            self.Z = parse.getint("main", "Z")
        print("A:", nA)
        print("Z:", self.Z)

        #alpha r_e^2
        self.ar2 = 7.297*2.818*2.818*1e-2 # m barn

        #electron and nucleus mass
        self.me = TDatabasePDG.Instance().GetParticle(11).Mass()
        mp = TDatabasePDG.Instance().GetParticle(2212).Mass()
        mn = nA*mp

        #nucleus boost vector
        nvec = TLorentzVector()
        pz_a = TMath.Sqrt(Ep**2-mp**2)*self.Z
        en_a = TMath.Sqrt(pz_a**2 + mn**2)
        nvec.SetPxPyPzE(0, 0, pz_a, en_a)
        self.nboost = nvec.BoostVector()

        #electron beam energy in nucleus beam rest frame
        evec = TLorentzVector()
        evec.SetPxPyPzE(0, 0, -TMath.Sqrt(self.Ee**2-self.me**2), self.Ee)
        evec.Boost(-self.nboost.x(), -self.nboost.y(), -self.nboost.z())
        self.Ee_n = evec.E()
        print("Ee_n (GeV):", self.Ee_n)

        #minimal photon energy in nucleus rest frame
        emin = 0.5 # GeV
        if parse.has_option("main", "emin"):
            emin = parse.getfloat("main", "emin")
        if parse.has_option("Lifshitz_bx", "emin"):
            emin = parse.getfloat("Lifshitz_bx", "emin")
        print("emin (GeV):", emin)
        eminv = TLorentzVector()
        eminv.SetPxPyPzE(0, 0, -emin, emin)
        eminv.Boost(-self.nboost.x(), -self.nboost.y(), -self.nboost.z())
        self.emin_n = eminv.E()
        print("emin_n (GeV):", self.emin_n)

        #maximal delta in nucleus frame
        self.dmax_n = 100.
        if parse.has_option("main", "dmax_n"):
            self.dmax_n = parse.getfloat("main", "dmax_n")
        if parse.has_option("Lifshitz_bx", "dmax_n"):
            self.dmax_n = parse.getfloat("Lifshitz_bx", "dmax_n")
        print("dmax_n:", self.dmax_n)

        #FOAM integration
        self.eq_foam = self.eq_93p16_foam(self)
        self.rand_foam = TRandom3()
        self.rand_foam.SetSeed(123)
        self.foam = TFoam("foam")
        self.foam.SetkDim(2)
        self.foam.SetnCells(4000)
        self.foam.SetRhoInt(self.eq_foam)
        self.foam.SetPseRan(self.rand_foam)
        self.foam.Initialize()

        #make samples and integrate for total cross section
        for i in range(1000000): self.foam.MakeEvent()
        foam_int = (c_double(0), c_double(0))
        self.foam.GetIntegMC(foam_int[0], foam_int[1])

        #total cross section (mb)
        sigma_tot = foam_int[0].value*self.dmax_n*(self.Ee_n-self.emin_n)
        sigma_tot_err = foam_int[1].value*self.dmax_n*(self.Ee_n-self.emin_n)

        print("Total cross section (mb):", sigma_tot, "+/-", sigma_tot_err)

        #uniform generator for azimuthal angles
        self.rand = TRandom3()
        self.rand.SetSeed(5572323)

        #Poisson distribution for number of interactions per bunch crossing
        self.nbunch = 0
        if parse.has_option("Lifshitz_bx", "nbunch"):
            self.nbunch = parse.getint("Lifshitz_bx", "nbunch")
        print("nbunch:", self.nbunch)
        #nbunch = 0 sets for single interaction simulated per event
        if self.nbunch > 1:
            self.pois = self.make_pois(parse, self.nbunch, sigma_tot)

        #tree output from the generator
        if self.nbunch == 0:
            tlist = ["true_phot_w", "true_phot_delta", "true_phot_theta_n"]
            tlist += ["true_phot_theta", "true_phot_phi", "true_phot_en"]
            tlist += ["true_phot_px", "true_phot_py", "true_phot_pz"]
            tlist += ["true_el_theta", "true_el_phi", "true_el_en"]
            tlist += ["true_el_px", "true_el_py", "true_el_pz"]
        else:
            tlist = ["num_interactions"]

        self.tree_out = self.set_tree(tree, tlist)

        #event attributes for hepmc
        self.hepmc_attrib = hepmc_attrib
        self.hepmc_vtx_start = 0 # first particle vertex id

        #local beam effects in case the model is used as background
        self.beff = None
        if parse.has_section("Lifshitz_bx"):
            self.beff = beam_effects(parse, None, "Lifshitz_bx")

        #to print the total cross section from FOAM at the end
        atexit.register(self.finish)

        print("Lifshitz_bx parametrization initialized")

    #_____________________________________________________________________________
    def generate(self, add_particle):

        #number of interactions in bunch crossing
        ni = 1
        if self.nbunch > 1:
            ni = int(TMath.Floor(self.pois.GetRandom()))
            self.tree_out.num_interactions = ni
            self.hepmc_attrib["num_interactions"] = ni

        for i in range(ni):
            self.make_phot_el(add_particle, i+self.hepmc_vtx_start)

    #_____________________________________________________________________________
    def make_phot_el(self, add_particle, ivtx):

        #photon energy and delta in nucleus rest frame by FOAM
        self.foam.MakeEvent()
        wdvec = std.vector("double")(2)
        self.foam.GetMCvect(wdvec.data())
        w = wdvec[0]*self.Ee_n
        d = wdvec[1]*self.dmax_n

        #polar (theta) and azimuthal (phi) angles in nucleus rest frame
        theta_n = d*self.me/self.Ee_n
        phi_n = 2. * TMath.Pi() * self.rand.Rndm() #uniform azimuthal angle

        #photon particle
        phot = add_particle( particle(22) )
        phot.stat = 1
        phot.vtx_id = ivtx
        phot.pxyze_prec = 9

        #set the photon vector in nucleus rest frame
        px_n = w*TMath.Sin(theta_n)*TMath.Cos(phi_n)
        py_n = w*TMath.Sin(theta_n)*TMath.Sin(phi_n)
        pz_n = w*TMath.Cos(theta_n)
        phot.vec.SetPxPyPzE(px_n, py_n, pz_n, w)

        #transform the photon vector to laboratory frame
        phot.vec.Boost(-self.nboost.x(), -self.nboost.y(), -self.nboost.z())

        #rotate the photon for pz < 0
        phot.vec.SetPxPyPzE(phot.vec.Px(), phot.vec.Py(), -phot.vec.Pz(), phot.vec.E())

        #scattered electron, initialize as beam
        electron = add_particle( beam(self.Ee, 11, -1) )
        electron.stat = 1
        electron.vtx_id = ivtx
        electron.pxyze_prec = 9

        #constrain the scattered electron with the photon
        electron.vec -= phot.vec

        #set the tree output and hepmc attributes for the case of one interaction
        if self.nbunch == 0:

            #cross section formula
            self.tree_out.true_phot_w = w
            self.tree_out.true_phot_delta = d
            self.tree_out.true_phot_theta_n = theta_n
            self.hepmc_attrib["true_phot_w"] = w
            self.hepmc_attrib["true_phot_delta"] = d
            self.hepmc_attrib["true_phot_theta_n"] = theta_n

            #photon and electron kinematics
            self.tree_out.true_phot_theta = phot.vec.Theta()
            self.tree_out.true_phot_phi = phot.vec.Phi()
            self.tree_out.true_phot_en = phot.vec.E()
            self.tree_out.true_phot_px = phot.vec.Px()
            self.tree_out.true_phot_py = phot.vec.Py()
            self.tree_out.true_phot_pz = phot.vec.Pz()
            self.tree_out.true_el_theta = electron.vec.Theta()
            self.tree_out.true_el_phi = electron.vec.Phi()
            self.tree_out.true_el_en = electron.vec.E()
            self.tree_out.true_el_px = electron.vec.Px()
            self.tree_out.true_el_py = electron.vec.Py()
            self.tree_out.true_el_pz = electron.vec.Pz()

            #kinematics in hepmc attributes
            self.hepmc_attrib["true_phot_theta"] = phot.vec.Theta()
            self.hepmc_attrib["true_phot_phi"] = phot.vec.Phi()
            self.hepmc_attrib["true_phot_en"] = phot.vec.E()
            self.hepmc_attrib["true_el_theta"] = electron.vec.Theta()
            self.hepmc_attrib["true_el_phi"] = electron.vec.Phi()
            self.hepmc_attrib["true_el_en"] = electron.vec.E()

        #beam effects for the photon and electron (local when as background)
        if self.beff is not None:
            self.beff.apply((phot, electron))

    #_____________________________________________________________________________
    class eq_93p16_foam:
        def __init__(self, gen):
            self.gen = gen
        def __call__(self, ndim, x):

            #Eq. 93.16 from Lifshitz QED for bremsstrahlung, version for TFoam

            #photon energy and angle, FOAM input ranges from 0 to 1
            w = x[0]*self.gen.Ee_n
            d = x[1]*self.gen.dmax_n

            if w < self.gen.emin_n: return 0.

            #initial and final electron E and E'
            E = self.gen.Ee_n
            Efin = E - w

            t1 = 8.*self.gen.Z*self.gen.Z*self.gen.ar2*(1./w)*(Efin/E)*d/((1+d**2)**2)

            t2 = ( (E/Efin) + (Efin/E) - 4*d*d/((1+d**2)**2) )*TMath.Log(2.*E*Efin/(self.gen.me*w))

            t3 = 0.5*( (E/Efin) + (Efin/E) + 2 - 16*d*d/((1+d**2)**2) )

            sig = t1*(t2 - t3)

            return sig

    #_____________________________________________________________________________
    def make_pois(self, parse, nbunch, sigma_tot):

        #Poisson distribution for number of interactions per bunch crossing

        #beam velocity (units of c)
        beta = TMath.Sqrt(self.Ee**2-self.me**2)/self.Ee
        print("Electron beam beta:", beta)

        #bunch spacing, sec
        Tb = parse.getfloat("Lifshitz_bx", "circ")/(beta*TMath.C()*nbunch)
        print("Bunch spacing (micro sec):", 1e6*Tb)
        print("Bunch frequency (MHz):", 1e-6/Tb)

        #luminosity per bunch crossing, mb^-1
        Lb = parse.getfloat("Lifshitz_bx", "L_inst")*1e-27*Tb
        print("Luminosity per bunch crossing (mb^-1):", Lb)

        #mean number of interactions per bunch crossing
        lam = sigma_tot*Lb
        print("Mean number of interactions per bunch crossing:", lam)

        #corresponding Poisson distribution
        pois = TF1("Pois", "TMath::Power([0], Int_t(TMath::Floor(x)) )\
            *TMath::Exp(-[0])/TMath::Factorial( Int_t(TMath::Floor(x)) )", 0, 12.*lam)
        pois.SetParameter(0, lam)
        pois.GetRandom() # some initialization in TF1 happens here

        return pois

    #_____________________________________________________________________________
    def finish(self):

        foam_int = c_double(0)
        foam_int_err = c_double(0)
        self.foam.Finalize(foam_int, foam_int_err)

        sig = foam_int.value*self.dmax_n*(self.Ee_n-self.emin_n)
        sig_err = foam_int_err.value*self.dmax_n*(self.Ee_n-self.emin_n)

        print("Total cross section (mb):", sig, "+/-", sig_err)

    #_____________________________________________________________________________
    def set_tree(self, tree, tlist):

        #set output to the tree

        #tree variables
        struct = "struct gen_Lifshitz_bx { Double_t "
        for i in tlist:
            struct += i + ", "
        struct = struct[:-2] + ";};"
        gROOT.ProcessLine( struct )
        tree_out = ROOT.gen_Lifshitz_bx()

        #put zero to all variables
        for i in tlist:
            exec("tree_out."+i+"=0")

        #add variables to the tree
        if tree is not None:
            for i in tlist:
                tree.Branch(i, addressof(tree_out, i), i+"/D")

        return tree_out














