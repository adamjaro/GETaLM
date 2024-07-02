
import ROOT as rt
from ROOT import TDatabasePDG, TF1, TRandom3, gRandom, TMath
from ROOT import gROOT, addressof

from photon import photon
from beam import beam

#_____________________________________________________________________________
class gen_zeus:
    #Bethe-Heitler bremsstrahlung photon according to ZEUS
    #Eur. Phys. J. C (2011) 71: 1574
    #_____________________________________________________________________________
    def __init__(self, parse, tree=None):

        #electron beam, GeV
        self.Ee = parse.getfloat("main", "Ee")
        #proton beam, GeV
        self.Ep = parse.getfloat("main", "Ep")

        print("Ee =", self.Ee, "GeV")
        print("Ep =", self.Ep, "GeV")

        #minimal photon energy, GeV
        self.emin = parse.getfloat("main", "emin")
        print("emin =", self.emin)

        #maximal photon angle
        self.tmax = 1.5e-3
        if parse.has_option("main", "tmax"):
            self.tmax = parse.getfloat("main", "tmax")

        print("tmax =", self.tmax)

        #electron and proton mass
        self.me = TDatabasePDG.Instance().GetParticle(11).Mass()
        self.mp = TDatabasePDG.Instance().GetParticle(2212).Mass()
        self.mep = self.me * self.mp

        #normalization,  4 alpha r_e^2
        self.ar2 = 4*7.297*2.818*2.818*1e-2 # m barn

        #parametrizations for dSigma/dE_gamma and dSigma/dtheta
        gRandom.SetSeed(5572323)
        self.eq1par = self.eq1(self)
        self.dSigDe = TF1("dSigDe", self.eq1par, self.emin, self.Ee)

        self.theta_const = 1e-11 # constant term in theta formula
        self.eq2par = self.eq2(self)
        self.dSigDtheta = TF1("dSigDtheta", self.eq2par, 0, self.tmax)

        #uniform generator for azimuthal angles
        self.rand = TRandom3()
        self.rand.SetSeed(5572323)

        #configure direct TTree output from the generator, used when 'tree'
        #argument is provided
        tlist = ["en", "theta", "phi"]
        tlist += ["true_phot_en", "true_phot_theta", "true_phot_phi"]
        tlist += ["true_phot_px", "true_phot_py", "true_phot_pz"]
        tlist += ["true_el_en", "true_el_theta", "true_el_phi"]
        tlist += ["true_el_px", "true_el_py", "true_el_pz"]
        self.tree_out = self.set_tree(tree, tlist)

        print("ZEUS parametrization initialized")
        print("Total cross section: "+str(self.dSigDe.Integral(self.emin, self.Ee))+" mb")

    #_____________________________________________________________________________
    class eq1:
        def __init__(self, gen):
            self.gen = gen
        def __call__(self, x, par):

            #E_gamma
            Eg = x[0]
            #electron and proton energy
            Ee = self.gen.Ee
            Ep = self.gen.Ep

            #scattered electron Ee'
            Escat = Ee - Eg
            #if Escat < 1e-5: return 0.

            t1 = Escat/(Eg*Ee)
            t2 = (Ee/Escat) + (Escat/Ee) - 2./3
            t3 = TMath.Log(4*Ep*Ee*Escat/(self.gen.mep*Eg)) - 1./2

            return self.gen.ar2*t1*t2*t3

    #_____________________________________________________________________________
    class eq2:
        def __init__(self, gen):
            self.gen = gen
        def __call__(self, x, par):

            #photon angular distribution
            t = x[0]

            return self.gen.theta_const * t/(( (self.gen.me/self.gen.Ee)**2 + t**2 )**2)

    #_____________________________________________________________________________
    def generate(self, add_particle):

        #initialize the scattered electron as a beam electron
        electron = add_particle( beam(self.Ee, 11, -1) )
        electron.stat = 1
        electron.pxyze_prec = 9

        #kinematics for the Bethe-Heitler bremsstrahlung photon
        #energy and polar angle
        en = self.dSigDe.GetRandom()
        theta = self.dSigDtheta.GetRandom()

        #azimuthal angle
        phi = 2. * TMath.Pi() * self.rand.Rndm()

        #set direct tree outputs for energy, theta and phi
        self.tree_out.en = en
        self.tree_out.theta = theta
        self.tree_out.phi = phi

        #put the Bethe-Heitler bremsstrahlung photon to the event
        phot = add_particle( photon(en, theta, phi) )
        phot.pxyze_prec = 9 # increase kinematics precision for the photon

        #constrain the scattered electron with the photon
        electron.vec -= phot.vec

        #photon and electron kinematics in tree output
        self.tree_out.true_phot_en = phot.vec.Energy()
        self.tree_out.true_phot_theta = phot.vec.Theta()
        self.tree_out.true_phot_phi = phot.vec.Phi()
        self.tree_out.true_phot_px = phot.vec.Px()
        self.tree_out.true_phot_py = phot.vec.Py()
        self.tree_out.true_phot_pz = phot.vec.Pz()
        self.tree_out.true_el_en = electron.vec.Energy()
        self.tree_out.true_el_theta = electron.vec.Theta()
        self.tree_out.true_el_phi = electron.vec.Phi()
        self.tree_out.true_el_px = electron.vec.Px()
        self.tree_out.true_el_py = electron.vec.Py()
        self.tree_out.true_el_pz = electron.vec.Pz()

    #_____________________________________________________________________________
    def set_tree(self, tree, tlist):

        #direct output to TTree if provided, Double_t values

        #tree variables
        struct = "struct gen_zeus { Double_t "
        for i in tlist:
            struct += i + ", "
        struct = struct[:-2] + ";};"
        gROOT.ProcessLine( struct )
        tree_out = rt.gen_zeus()

        #put zero to all variables
        for i in tlist:
            exec("tree_out."+i+"=0")

        #add variables to the tree
        if tree is not None:
            for i in tlist:
                tree.Branch(i, addressof(tree_out, i), i+"/D")

        return tree_out























