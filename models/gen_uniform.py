
#_____________________________________________________________________________
# Photons with uniform energy for simulations testing
#
#_____________________________________________________________________________

import ROOT as rt
from ROOT import TRandom3, gROOT, addressof, TDatabasePDG, TMath

from particle import particle
from beam import beam

#_____________________________________________________________________________
class gen_uniform:
    #_____________________________________________________________________________
    def __init__(self, parse, tree, hepmc_attrib):

        print("Uniform configuration:")
        
        #electron and proton beam energy, GeV
        self.Ee = parse.getfloat("main", "Ee")
        self.Ep = parse.getfloat("main", "Ep")
        print("Ee =", self.Ee, "GeV")
        print("Ep =", self.Ep, "GeV")
        
        #minumum and maximum energy, GeV
        self.emin = parse.getfloat("main", "emin")
        self.emax = parse.getfloat("main", "emax")
        print("emin =", self.emin)
        print("emax =", self.emax)

        #pdg for generated particle, electron or photon
        self.pdg = 22
        if parse.has_option("main", "pdg"):
            self.pdg = parse.getint("main", "pdg")

        print("pdg =", self.pdg)

        #angular range
        self.theta_min = 0.
        self.theta_max = 0.

        #maximal polar angle as pi - theta_max in rad
        if parse.has_option("main", "pitheta_max"):
            self.theta_min = TMath.Pi() - parse.getfloat("main", "pitheta_max")
            self.theta_max = TMath.Pi()

        #angles as mlt = -log_10(pi - theta)
        if parse.has_option("main", "mlt_min"):
            mlt_min = parse.getfloat("main", "mlt_min")
            mlt_max = parse.getfloat("main", "mlt_max")
            print("mlt_min =", mlt_min)
            print("mlt_max =", mlt_max)
            self.theta_min = TMath.Pi() - 10.**(-mlt_min)
            self.theta_max = TMath.Pi() - 10.**(-mlt_max)

        print("theta_min =", self.theta_min)
        print("theta_max =", self.theta_max)

        #generator functions for photons and electrons
        self.gen_func = {}
        self.gen_func[22] = self.gen_phot
        self.gen_func[11] = self.gen_el

        #test for pdg
        if self.gen_func.get(self.pdg) is None:
            print("Fatal: pdg", self.pdg, "is not supported")
            raise KeyError

        #uniform generator
        self.rand = TRandom3()
        self.rand.SetSeed(5572323)

        #electron mass
        self.me = TDatabasePDG.Instance().GetParticle(11).Mass()

        #set the output tree
        if self.pdg == 22:
            tnam = ["true_phot_theta", "true_phot_phi", "true_phot_E"]
        if self.pdg == 11:
            tnam = ["true_el_pT", "true_el_theta", "true_el_phi", "true_el_E"]
        self.out = self.make_tree(tree, tnam)

        #event attributes for hepmc
        self.hepmc_attrib = hepmc_attrib

        print("Uniform generator initialized")

    #_____________________________________________________________________________
    def gen_phot(self, add_particle):

        #photon generator

        #energy for the current event
        en = self.emin + (self.emax - self.emin) * self.rand.Rndm()

        #polar and azimuthal angles
        theta = self.theta_min + (self.theta_max - self.theta_min) * self.rand.Rndm()
        phi = 2. * TMath.Pi() * self.rand.Rndm()

        #photon momentum
        px = en * TMath.Sin(theta) * TMath.Cos(phi)
        py = en * TMath.Sin(theta) * TMath.Sin(phi)
        pz = en * TMath.Cos(theta)

        #make the photon
        phot = particle(22)
        phot.vec.SetPxPyPzE(px, py, pz, en)
        phot.stat = 1
        phot.pxyze_prec = 9

        #photon kinematics in the output
        self.out.true_phot_theta = phot.vec.Theta()
        self.out.true_phot_phi = phot.vec.Phi()
        self.out.true_phot_E = phot.vec.E()
        self.hepmc_attrib["true_phot_theta"] = phot.vec.Theta()
        self.hepmc_attrib["true_phot_phi"] = phot.vec.Phi()
        self.hepmc_attrib["true_phot_E"] = phot.vec.E()

        #put the photon to the event
        add_particle( phot )

    #_____________________________________________________________________________
    def gen_el(self, add_particle):

        #electron generator

        #energy
        en = self.emin + (self.emax - self.emin) * self.rand.Rndm()

        #polar angle
        theta = self.theta_min + (self.theta_max - self.theta_min) * self.rand.Rndm()

        #azimuthal angle
        phi = 2. * TMath.Pi() * self.rand.Rndm()

        #electron momentum
        ptot = TMath.Sqrt(en**2 - self.me**2)
        px = ptot * TMath.Sin(theta) * TMath.Cos(phi)
        py = ptot * TMath.Sin(theta) * TMath.Sin(phi)
        pz = ptot * TMath.Cos(theta)

        #make the electron
        el = particle(11)
        el.vec.SetPxPyPzE(px, py, pz, en)
        el.stat = 1
        el.pxyze_prec = 9

        #electron kinematics in output tree
        self.out.true_el_pT = el.vec.Pt()
        self.out.true_el_theta = el.vec.Theta()
        self.out.true_el_phi = el.vec.Phi()
        self.out.true_el_E = el.vec.E()

        #electron kinematics in hepmc
        self.hepmc_attrib["true_el_pT"] = self.out.true_el_pT
        self.hepmc_attrib["true_el_theta"] = self.out.true_el_theta
        self.hepmc_attrib["true_el_phi"] = self.out.true_el_phi
        self.hepmc_attrib["true_el_E"] = self.out.true_el_E

        #put the electron to the event
        add_particle( el )

    #_____________________________________________________________________________
    def generate(self, add_particle):

        #beam electron
        ebeam = add_particle( beam(self.Ee, 11, -1) )
        ebeam.stat = 4
        ebeam.pxyze_prec = 9
        
        #beam proton
        pbeam = add_particle( beam(self.Ep, 2212, 1) )
        pbeam.stat = 4
        pbeam.pxyze_prec = 9

        #generator for a given pdg

        self.gen_func[self.pdg](add_particle)
        


    #_____________________________________________________________________________
    def make_tree(self, tree, tnam):

        #create the tree variables
        tcmd = "struct gen_out { Double_t "
        for i in tnam:
            tcmd += i + ", "
        tcmd = tcmd[:-2] + ";};"
        gROOT.ProcessLine( tcmd )
        self.out = rt.gen_out()

        #put zero to all variables
        for i in tnam:
            exec("self.out."+i+"=0")

        #set the variables in the tree
        if tree is not None:
            for i in tnam:
                tree.Branch(i, addressof(self.out, i), i+"/D")

        return self.out









