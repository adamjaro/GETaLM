
import atexit

import ROOT as rt
from ROOT import TFile, TTree, gROOT, addressof, TClonesArray

from beam import beam

#_____________________________________________________________________________
class file_output:
    #output from the generator
    #_____________________________________________________________________________
    def __init__(self, parse):

        #create the individual outputs
        self.set_write_tx = False
        if parse.has_option("main", "write_tx"):
            self.set_write_tx = parse.getboolean("main", "write_tx")

        self.set_write_root = True
        if parse.has_option("main", "write_root"):
            self.set_write_root = parse.getboolean("main", "write_root")

        self.set_write_hepmc = False
        if parse.has_option("main", "write_hepmc"):
            self.set_write_hepmc = parse.getboolean("main", "write_hepmc")

        if self.set_write_tx: self.make_tx(parse)

        self.ltree = None
        if self.set_write_root: self.make_root(parse)

        self.hepmc_attrib = {}
        if self.set_write_hepmc: self.make_hepmc(parse)

    #_____________________________________________________________________________
    def make_tx(self, parse):

        #TX output

        nam = parse.get("main", "nam").strip("\"'") + ".tx"
        print("TX output name:", nam)

        self.tx_out = open(nam, "w")
        self.tx_ievt = 1

    #_____________________________________________________________________________
    def make_root(self, parse):

        #ROOT output

        nam = parse.get("main", "nam").strip("\"'") + ".root"
        print("ROOT output name:", nam)

        self.out_root = TFile(nam, "recreate")
        #tree variables, all Double_t
        tlist = ["phot_en", "phot_theta", "phot_phi"]
        tlist += ["el_en", "el_theta", "el_phi"]
        #C structure holding the variables
        struct = "struct tree_out { Double_t "
        for i in tlist: struct += i + ", "
        struct = struct[:-2] + ";};"
        gROOT.ProcessLine( struct )
        #create the output tree
        self.tree_out = rt.tree_out() # instance of the C structure
        self.ltree = TTree("ltree", "ltree")
        for i in tlist:
            exec("self.tree_out."+i+"=0")
            self.ltree.Branch(i, addressof(self.tree_out, i), i+"/D")

        #particles array
        self.particles_out = TClonesArray("TParticle")
        self.particles_out.SetOwner(True)
        self.ltree.Branch("particles", self.particles_out)

        atexit.register(self.close_root)

    #_____________________________________________________________________________
    def make_hepmc(self, parse):

        #HepMC3 output

        global hepmc
        from pyHepMC3 import HepMC3 as hepmc

        nam = parse.get("main", "nam").strip("\"'") + ".hepmc"
        print("HepMC3 output name:", nam)

        self.hepmc_out = hepmc.WriterAscii(nam, hepmc.GenRunInfo())
        self.hepmc_ievt = 0

        #electron and proton beam energy to create primary vertex
        self.hepmc_Ee = parse.getfloat("main", "Ee")
        self.hepmc_Ep = parse.getfloat("main", "Ep")

    #_____________________________________________________________________________
    def write_tx(self, tracks):

        #TX Starlight format

        if not self.set_write_tx: return

        #tracks and vertex position in cm
        tracks_tx = []
        vx = 0.
        vy = 0.
        vz = 0.
        #tracks loop
        for t in tracks:
            #only final particles
            if t.stat != 1: continue

            vx = t.vx/10.
            vy = t.vy/10.
            vz = t.vz/10.

            t.write_tx(tracks_tx)

        #number of tracks for event and vertex lines
        ntrk = str(len(tracks_tx))

        #event line
        evtlin = "EVENT: "+str(self.tx_ievt)+" "+ntrk+" 1"
        self.tx_out.write(evtlin+"\n")

        #vertex line
        vtxlin = "VERTEX:"
        vtx_prec = 9
        if abs(vx)<1e-9 and abs(vy)<1e-9 and abs(vz)<1e-9:
            vtx_prec = 0
        vtx_form = " {0:."+str(vtx_prec)+"f}"
        vtxlin += vtx_form.format(vx)
        vtxlin += vtx_form.format(vy)
        vtxlin += vtx_form.format(vz)
        vtxlin += " 0 1 0 0 "+ntrk
        self.tx_out.write(vtxlin+"\n")

        #track lines
        for tlin in tracks_tx:
            self.tx_out.write(tlin+"\n")

        self.tx_ievt += 1

    #_____________________________________________________________________________
    def write_root(self, tracks):

        #ROOT output

        if not self.set_write_root: return

        #initialize the particles array
        ipos = 0
        self.particles_out.Clear("C")

        t = self.tree_out

        for i in tracks:
            #select the final photon and electron
            if i.stat != 1: continue

            #put the particles to TParticles clones array
            i.write_tparticle(self.particles_out, ipos)
            ipos += 1

            #final photon
            if i.pdg == 22:

                t.phot_en    = i.vec.Energy()
                t.phot_theta = i.vec.Theta()
                t.phot_phi   = i.vec.Phi()

            #final electron
            if i.pdg == 11:

                t.el_en     = i.vec.Energy()
                t.el_theta  = i.vec.Theta()
                t.el_phi    = i.vec.Phi()

        #fill the tree
        self.ltree.Fill()

    #_____________________________________________________________________________
    def write_hepmc(self, tracks):

        #HepMC3 format

        if not self.set_write_hepmc: return

        #hepmc event
        evt = hepmc.GenEvent(hepmc.Units.GEV, hepmc.Units.MM)
        evt.set_event_number(self.hepmc_ievt);

        #event attributes
        for i in self.hepmc_attrib:
            attr = hepmc.DoubleAttribute(self.hepmc_attrib[i])
            evt.add_attribute(i, attr)

        #primary vertex
        vtx = hepmc.GenVertex( hepmc.FourVector(tracks[0].vx, tracks[0].vy, tracks[0].vz, 0) )

        #beam electron and proton as incoming particles to primary vertex
        beam_el = beam(self.hepmc_Ee, 11, -1)
        beam_prot = beam(self.hepmc_Ep, 2212, 1)
        vtx.add_particle_in( hepmc.GenParticle(hepmc.FourVector(beam_el.vec.Px(), beam_el.vec.Py(), beam_el.vec.Pz(),\
            beam_el.vec.E()), beam_el.pdg, 4) )
        vtx.add_particle_in( hepmc.GenParticle(hepmc.FourVector(beam_prot.vec.Px(), beam_prot.vec.Py(),\
            beam_prot.vec.Pz(), beam_prot.vec.E()), beam_prot.pdg, 4) )

        #tracks loop
        for t in tracks:
            #only final particles
            if t.stat != 1: continue

            vtx.add_particle_out( t.make_hepmc_particle(hepmc) )

        evt.add_vertex(vtx)

        self.hepmc_out.write_event(evt)

        self.hepmc_ievt += 1

    #_____________________________________________________________________________
    def close_root(self):

        self.out_root.Write()
        self.out_root.Close()

















