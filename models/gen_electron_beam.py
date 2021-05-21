
#_____________________________________________________________________________
# Beam electrons with a given energy for simulations testing
#
#_____________________________________________________________________________

from ROOT import TMath, TDatabasePDG, TF1

from particle import particle

#_____________________________________________________________________________
class gen_electron_beam:
    #_____________________________________________________________________________
    def __init__(self, parse):

        #electron beam energy, GeV
        self.Ee = parse.getfloat("main", "Ee")
        print("Ee =", self.Ee)

        #energy spread from input relative energy spread
        self.espread = None
        if parse.has_option("main", "espread"):
            sig_e = self.Ee * parse.getfloat("main", "espread")
            self.espread = TF1("espread", "gaus", -12*sig_e, 12*sig_e)
            self.espread.SetParameters(1, 0, sig_e)

        #electron mass
        self.me = TDatabasePDG.Instance().GetParticle(11).Mass()

        print("Electron beam initialized")

    #_____________________________________________________________________________
    def generate(self, add_particle):

        #energy spread
        de = 0
        if self.espread is not None:
            de = self.espread.GetRandom()

        #electron energy and momentum along z
        en = self.Ee + de
        pz = -TMath.Sqrt(en**2 - self.me**2)

        #beam Lorentz vector
        beam = particle(11)
        beam.vec.SetPxPyPzE(0, 0, pz, en)
        beam.stat = 1
        beam.pxyze_prec = 9

        #put the beam electron to the event
        add_particle( beam )

















