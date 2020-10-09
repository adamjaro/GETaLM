
from math import sqrt

from ROOT import TDatabasePDG

from particle import particle

#_____________________________________________________________________________
class beam(particle):
    #beam particle
    #_____________________________________________________________________________
    def __init__(self, en, pdg, zsign):
        #energy, pdg, sign of z-momentum

        particle.__init__(self, pdg)
        self.pdg = pdg
        #status code for beam particles
        self.stat = 201
        #set kinematics for beam particle
        m = TDatabasePDG.Instance().GetParticle(pdg).Mass()
        pz = zsign*sqrt( en**2 - m**2 )
        self.vec.SetPxPyPzE(0, 0, pz, en)



