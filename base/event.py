
import configparser
from sys import stdout

from ROOT import gROOT

from beam_effects import beam_effects
from file_output import file_output

from gen_h1 import gen_h1
from gen_zeus import gen_zeus
from gen_quasi_real import gen_quasi_real
from gen_electron_beam import gen_electron_beam
from gen_read_py import gen_read_py
from gen_uniform import gen_uniform
from gen_Lifshitz_93p16 import gen_Lifshitz_93p16
from gen_beam_gas import gen_beam_gas
from gen_Lifshitz_bx import gen_Lifshitz_bx
from gen_sig_bx import gen_sig_bx

#_____________________________________________________________________________
class event:
    #generated event
    #_____________________________________________________________________________
    def __init__(self, config):

        #input configuration
        print("Generator configuration:", config)
        parse = configparser.RawConfigParser(inline_comment_prefixes=(";","#"))
        parse.read(config)

        if not parse.has_section("main"):
            print("Failed to load the configuration")
            return

        print("ROOT version:", gROOT.GetVersion())

        #output
        self.out = file_output(parse)

        #physics model
        model = parse.get("main", "model").strip("\"'")
        if model == "h1":
            self.gen = gen_h1(parse)
        elif model == "zeus":
            self.gen = gen_zeus(parse, self.out.ltree, self.out.hepmc_attrib)
        elif model == "quasi-real":
            self.gen = gen_quasi_real(parse, self.out.ltree, self.out.hepmc_attrib)
        elif model == "electron-beam":
            self.gen = gen_electron_beam(parse)
        elif model == "read-py":
            self.gen = gen_read_py(parse, self.out.ltree, self.out.hepmc_attrib)
        elif model == "uniform":
            self.gen = gen_uniform(parse, self.out.ltree, self.out.hepmc_attrib)
        elif model == "Lifshitz_93p16":
            self.gen = gen_Lifshitz_93p16(parse, self.out.ltree, self.out.hepmc_attrib)
        elif model == "beam-gas":
            self.gen = gen_beam_gas(parse, self.out.ltree)
        elif model == "Lifshitz_bx":
            self.gen = gen_Lifshitz_bx(parse, self.out.ltree, self.out.hepmc_attrib)
        elif model == "sig_bx":
            self.gen = gen_sig_bx(parse, self.out.ltree, self.out.hepmc_attrib)
        else:
            print("Invalid generator specified")
            exit()

        #beam effects
        self.beff = None
        if parse.has_section("beam_effects") == True:
            self.beff = beam_effects(parse, self.out.ltree, "beam_effects", self.out.hepmc_attrib)

        #tracks in the event
        self.tracks = []

        #run
        nev = parse.getint("main", "nev") # number of events to generate
        print("Number of events:", nev)
        self.event_loop(nev)

    #_____________________________________________________________________________
    def event_loop(self, nev):

        iprint = int(nev/12)
        if iprint == 0:
            iprint = 12

        for i in range(nev):
            if i%iprint == 0 and i>0:
                print("{0:.1f} %".format(100.*i/nev))
                stdout.flush()
            self.generate()

        print("All done")

    #_____________________________________________________________________________
    def generate(self):

        #generate the event

        #clear tracks list
        self.tracks = []

        #run the physics model
        self.gen.generate(self.add_particle)

        #apply the beam effects
        if self.beff is not None:
            self.beff.apply(self.tracks)

        #TX, ROOT and HepMC3 outputs
        self.out.write_tx(self.tracks)
        self.out.write_root(self.tracks)
        self.out.write_hepmc(self.tracks)

    #_____________________________________________________________________________
    def add_particle(self, part):

        #keep track of particle index in event
        part.idx = len(self.tracks)+1
        self.tracks.append(part)

        return part





















