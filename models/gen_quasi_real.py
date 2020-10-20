
#_____________________________________________________________________________
# Quasi-real photoproduction at low Q^2
#
# d^2 sigma / dxdy is Eq. II.6 in Amaldi, ed., Study on ep facility, 1979, page 320,
# Conf.Proc. C790402 (1979) 1-474, http://inspirehep.net/record/151135
#
# Version 2 improvement: the cross section is transformed to log_10(x) and log_10(y)
# for arbitrary kinematics reach
#
# The total gamma-proton cross section is from Donnachie, Landshoff, Total cross sections,
# Phys.Lett. B296 (1992) 227-232, http://inspirehep.net/record/337839
#
# Similar procedure was followed in S. Levonian, H1LUMI, H1-04/93-287 (1993),
# https://www-h1.desy.de/~levonian/papers/h1lumi.ps.gz
# 
#_____________________________________________________________________________

import math
from math import pi
import atexit

import ROOT as rt
from ROOT import TF2, Double, TMath, TRandom3, gROOT, AddressOf, TDatabasePDG
from ROOT import TLorentzVector

from particle import particle

#_____________________________________________________________________________
class gen_quasi_real:
    #_____________________________________________________________________________
    def __init__(self, parse, tree):

        print "Quasi-real configuration:"

        #electron and proton beam energy, GeV
        self.Ee = parse.getfloat("main", "Ee")
        self.Ep = parse.getfloat("main", "Ep")
        print "Ee =", self.Ee, "GeV"
        print "Ep =", self.Ep, "GeV"

        #electron and proton mass
        self.me = TDatabasePDG.Instance().GetParticle(11).Mass()
        mp = TDatabasePDG.Instance().GetParticle(2212).Mass()

        #boost vector pbvec of proton beam
        pbeam = TLorentzVector()
        pbeam.SetPxPyPzE(0, 0, TMath.Sqrt(self.Ep**2-mp**2), self.Ep)
        self.pbvec = pbeam.BoostVector()

        #electron beam energy Ee_p in proton beam rest frame
        ebeam = TLorentzVector()
        ebeam.SetPxPyPzE(0, 0, -TMath.Sqrt(self.Ee**2-self.me**2), self.Ee)
        ebeam.Boost(-self.pbvec.x(), -self.pbvec.y(), -self.pbvec.z()) # transform to proton beam frame
        self.Ee_p = ebeam.E()

        #center-of-mass squared s, GeV^2
        self.s = self.get_s(self.Ee, self.Ep)
        print "s =", self.s, "GeV^2"
        print "sqrt(s) =", TMath.Sqrt(self.s), "GeV"

        #range in x
        xmin = parse.getfloat("main", "xmin")
        xmax = parse.getfloat("main", "xmax")
        print "xmin =", xmin
        print "xmax =", xmax

        #range in u = log_10(x)
        umin = TMath.Log10(xmin)
        umax = TMath.Log10(xmax)
        print "umin =", umin
        print "umax =", umax

        #range in y
        ymin = parse.getfloat("main", "ymin")
        ymax = parse.getfloat("main", "ymax")

        #range in W
        wmin = -1.
        wmax = -1.
        if parse.has_option("main", "Wmin"):
            wmin = parse.getfloat("main", "Wmin")
            print "Wmin =", wmin
        if parse.has_option("main", "Wmax"):
            wmax = parse.getfloat("main", "Wmax")
            print "Wmax =", wmax

        #adjust range in y according to W
        if wmin > 0 and ymin < wmin**2/self.s:
            ymin = wmin**2/self.s
        if wmax > 0 and ymax > wmax**2/self.s:
            ymax = wmax**2/self.s

        print "ymin =", ymin
        print "ymax =", ymax

        #range in v = log_10(y)
        vmin = TMath.Log10(ymin)
        vmax = TMath.Log10(ymax)
        print "vmin =", vmin
        print "vmax =", vmax

        #range in Q2
        self.Q2min = parse.getfloat("main", "Q2min")
        self.Q2max = parse.getfloat("main", "Q2max")
        print "Q2min =", self.Q2min
        print "Q2max =", self.Q2max

        #constant term in the cross section
        self.const = TMath.Log(10)*TMath.Log(10)*(1./137)/(2.*math.pi)

        #cross section formula for d^2 sigma / dxdy, Eq. II.6
        #transformed as x -> u = log_10(x) and y -> v = log_10(y)
        self.eq = TF2("d2SigDuDvII6", self.eq_II6_uv, umin, umax, vmin, vmax)

        self.eq.SetNpx(1000)
        self.eq.SetNpy(1000)

        #uniform generator for azimuthal angles
        self.rand = TRandom3()
        self.rand.SetSeed(5572323)

        #generator event variables in output tree
        tnam = ["gen_u", "gen_v", "true_x", "true_y", "true_Q2", "true_W2"]
        tnam += ["true_el_Q2"]
        tnam += ["true_el_pT", "true_el_theta", "true_el_phi", "true_el_E"]

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

        #counters for all generated and selected events
        self.nall = 0
        self.nsel = 0

        #print generator statistics at the end
        atexit.register(self.show_stat)

        #total integrated cross section
        self.sigma_tot = self.eq.Integral(umin, umax, vmin, vmax)
        print "Total integrated cross section for a given x and y range:", self.sigma_tot, "mb"

        print "Quasi-real photoproduction initialized"

    #_____________________________________________________________________________
    def generate(self, add_particle):

        #get the x and y for a given range in Q^2
        while True:

            #values of u = log_10(x) and v = log_10(y) from the cross section
            u = Double(0)
            v = Double(0)
            self.eq.GetRandom2(u, v)

            #x and y from the transformation
            x = 10.**u
            y = 10.**v

            #scattered electron energy and polar angle in proton beam rest frame
            en_p = self.Ee_p*(1 - y)
            theta_p = 2 * TMath.ASin( 0.5*TMath.Sqrt(x*y*self.s/((1-y)*self.Ee_p**2)) )

            #prevent unphysical energy and angle
            if theta_p < 0. or theta_p > pi: continue
            if en_p**2 < self.me**2: continue

            #event Q^2
            Q2 = x*y*self.s

            self.nall += 1 # increment all events counter

            #select the range in Q^2
            if Q2 < self.Q2min or Q2 > self.Q2max: continue

            self.nsel += 1 # increment selected counter

            break

        #scattered electron in the event
        el = add_particle( particle(11) )
        el.stat = 1
        el.pxyze_prec = 9

        #kinematics for scattered electron in proton beam rest frame
        phi_p = 2. * TMath.Pi() * self.rand.Rndm() #uniform azimuthal angle
        elp_p = TMath.Sqrt(en_p**2 - self.me**2) #total momentum

        #set the scattered electron vector in proton beam rest frame
        px_p = elp_p*TMath.Sin(theta_p)*TMath.Cos(phi_p)
        py_p = elp_p*TMath.Sin(theta_p)*TMath.Sin(phi_p)
        pz_p = elp_p*TMath.Cos(theta_p)

        el.vec.SetPxPyPzE(px_p, py_p, pz_p, en_p)

        #transform the scattered electron vector to the laboratory frame, negative z direction
        el.vec.Boost(-self.pbvec.x(), -self.pbvec.y(), -self.pbvec.z()) # transform to lab
        el.vec.SetPxPyPzE(el.vec.Px(), el.vec.Py(), -el.vec.Pz(), el.vec.E()) # rotate for pz<0

        #tree output with generator kinematics
        self.out.gen_u = u
        self.out.gen_v = v
        self.out.true_x = x
        self.out.true_y = y
        self.out.true_Q2 = Q2
        self.out.true_W2 = self.s*y

        #electron kinematics
        self.out.true_el_pT = el.vec.Pt()
        self.out.true_el_theta = el.vec.Theta()
        self.out.true_el_phi = el.vec.Phi()
        self.out.true_el_E = el.vec.E()

        #Q^2 from electron energy and angle
        self.out.true_el_Q2 = 2.*self.Ee*el.vec.E()*(1.-TMath.Cos(TMath.Pi()-el.vec.Theta()))

        #print Q2, self.out.gen_el_Q2

    #_____________________________________________________________________________
    def eq_II6_uv(self, val):

        #II.6 transformed as x -> u = log_10(x) and y -> v = log_10(y)
        u = val[0]
        v = val[1]

        sig = self.const*( 1. + (1.-10**v)**2 )*(1.-10.**u)

        #term of the total gamma-p cross section by Donnachie, Landshoff, 1992
        sig *= 0.0677*((10.**v)*self.s)**0.0808 + 0.129*((10.**v)*self.s)**-0.4525 # mb

        return sig

    #_____________________________________________________________________________
    def get_s(self, Ee, Ep):

        #calculate the CMS squared s

        #proton mass
        mp = TDatabasePDG.Instance().GetParticle(2212).Mass()

        #CMS energy squared s, GeV^2
        s = 2.*Ee*Ep + self.me**2 + mp**2
        s += 2*TMath.Sqrt(Ee**2 - self.me**2) * TMath.Sqrt(Ep**2 - mp**2)

        #print "sqrt(s):", TMath.Sqrt(s)

        return s

    #_____________________________________________________________________________
    def show_stat(self):

        #print generator statistics at the end
        print "Total generated events:", self.nall
        print "Selected events:", self.nsel

        if self.nall <= 0: return

        print "Cross section of generated sample:", self.sigma_tot*float(self.nsel)/float(self.nall), "mb"












