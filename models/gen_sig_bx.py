
#_____________________________________________________________________________
# Signal from selected model appended by bunch crossing background
# from gen_Lifshitz_bx
#
#_____________________________________________________________________________

from gen_quasi_real import gen_quasi_real
from gen_Lifshitz_bx import gen_Lifshitz_bx
from beam_effects import beam_effects

#_____________________________________________________________________________
class gen_sig_bx:
    #_____________________________________________________________________________
    def __init__(self, parse, tree, hepmc_attrib):

        #signal model
        sig_name = parse.get("main", "signal_model").strip("\"'")
        if sig_name == "quasi-real":
            self.sig_model = gen_quasi_real(parse, tree, hepmc_attrib)
        else:
            print("In gen_sig_bx: invalid signal model specified")
            exit()

        #background by gen_Lifshitz_bx
        self.bkg_model = gen_Lifshitz_bx(parse, tree, hepmc_attrib)

    #_____________________________________________________________________________
    def generate(self, add_particle):

        #signal particles in event
        self.sig_particles = []
        self.sig_model.generate(self.add_sig_particle)

        #beam effects via gen_Lifshitz_bx for signal particles
        self.bkg_model.beff.apply( self.sig_particles )

        #put signal particles to the event
        for i in self.sig_particles:
            add_particle(i)

        #generate background particles
        self.bkg_model.generate(add_particle)

    #_____________________________________________________________________________
    def add_sig_particle(self, particle):

        #add signal particle
        self.sig_particles.append(particle)

        return particle




















