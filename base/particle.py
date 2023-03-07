
from ROOT import TLorentzVector, TDatabasePDG

#generic particle

#_____________________________________________________________________________
class particle:
    #_____________________________________________________________________________
    def __init__(self, pdg):
        #particle Lorentz vector
        self.vec = TLorentzVector()
        #index in particle list
        self.idx = 0
        #status code 1=decay to be tracked, 3=beam particle
        self.stat = 1
        #pdg code
        self.pdg = pdg
        #particle database for pass and codes
        self.pdgdat = TDatabasePDG.Instance()
        #mass, GeV
        self.mass = self.pdgdat.GetParticle(self.pdg).Mass()
        #parent particle id
        self.parent_id = 0
        #vertex coordinates, mm
        self.vx = 0.
        self.vy = 0.
        self.vz = 0.
        self.vtx_id = 0
        #precision for momentum and energy
        self.pxyze_prec = 6

    #_____________________________________________________________________________
    def write_tx(self, track_list):

        #output line in TX format

        #Geant code and momentum
        lin = "TRACK:  "+str(self.pdgdat.ConvertPdgToGeant3(self.pdg))
        pxyz_form = " {0:."+str(self.pxyze_prec)+"f}"
        lin += pxyz_form.format( self.vec.Px() )
        lin += pxyz_form.format( self.vec.Py() )
        lin += pxyz_form.format( self.vec.Pz() )

        #track id
        lin += " " + str(len(track_list))

        #start and stop vertex and pdg
        lin += " 1 0 " + str(self.pdg)

        track_list.append(lin)

    #_____________________________________________________________________________
    def write_tparticle(self, particles, ipos):

        #write to TParticle clones array

        p = particles.ConstructedAt(ipos)

        p.SetMomentum(self.vec)
        p.SetPdgCode(self.pdg)
        p.SetProductionVertex(self.vx, self.vy, self.vz, 0)

    #_____________________________________________________________________________
    def make_hepmc_particle(self, hepmc):

        #create HepMC3 particle

        p = hepmc.GenParticle(hepmc.FourVector(self.vec.Px(), self.vec.Py(), self.vec.Pz(), self.vec.E()), self.pdg, self.stat)

        return p










