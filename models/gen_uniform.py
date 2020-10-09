
#_____________________________________________________________________________
# Photons with uniform energy for simulations testing
#
#_____________________________________________________________________________

import ROOT as rt
from ROOT import TRandom3, gROOT, AddressOf

from particle import particle

#_____________________________________________________________________________
class gen_uniform:
    #_____________________________________________________________________________
    def __init__(self, parse, tree):

        #minumum and maximum photon energy, GeV
        self.emin = parse.getfloat("main", "emin")
        self.emax = parse.getfloat("main", "emax")
        print "emin =", self.emin
        print "emax =", self.emax

        #uniform generator for photon energies
        self.rand = TRandom3()
        self.rand.SetSeed(5572323)

        #set the output tree
        tnam = ["gen_E"]
        self.out = self.make_tree(tree, tnam)

        print "Uniform generator initialized"

    #_____________________________________________________________________________
    def generate(self, add_particle):

        #energy for the current event
        en = self.emin + (self.emax - self.emin) * self.rand.Rndm()
        self.out.gen_E = en

        #print en

        #make the photon
        phot = particle(22)
        phot.vec.SetPxPyPzE(0, 0, -en, en)
        phot.stat = 1
        phot.pxyze_prec = 9

        #put the photon to the event
        add_particle( phot )


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
                tree.Branch(i, AddressOf(self.out, i), i+"/D")

        return self.out









