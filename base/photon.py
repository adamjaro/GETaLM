
from math import sqrt, cos, log, tan, pi

from particle import particle

#bremsstrahlung photon

#_____________________________________________________________________________
class photon(particle):
    #_____________________________________________________________________________
    def __init__(self, en, theta, phi):
        particle.__init__(self, 22)
        #set the photon given its energy en and angles theta and phi
        pt = en*sqrt( (1.-cos(2.*theta))/2. )
        #for pz negative
        theta = pi - theta
        eta = -log( tan(theta/2.) )
        #set Lorentz vector
        self.vec.SetPtEtaPhiE(pt, eta, phi, en)
        #status for final particle
        self.stat = 1



