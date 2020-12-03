
import ConfigParser
from sys import stdout

from beam_effects import beam_effects
from file_output import file_output

from gen_h1 import gen_h1
from gen_zeus import gen_zeus
from gen_quasi_real import gen_quasi_real
from gen_electron_beam import gen_electron_beam
from gen_read_py import gen_read_py
from gen_uniform import gen_uniform
from gen_Lifshitz_93p16 import gen_Lifshitz_93p16

#_____________________________________________________________________________
class event:
    #generated event
    #_____________________________________________________________________________
    def __init__(self, config):

        #input configuration
        print "Generator configuration:", config
        parse = ConfigParser.RawConfigParser()
        parse.read(config)

        #output
        self.out = file_output(parse)

        #physics model
        model = parse.get("main", "model").strip("\"'")
        if model == "h1":
            self.gen = gen_h1(parse)
        elif model == "zeus":
            self.gen = gen_zeus(parse)
        elif model == "quasi-real":
            self.gen = gen_quasi_real(parse, self.out.ltree)
        elif model == "electron-beam":
            self.gen = gen_electron_beam(parse)
        elif model == "read-py":
            self.gen = gen_read_py(parse, self.out.ltree)
        elif model == "uniform":
            self.gen = gen_uniform(parse, self.out.ltree)
        elif model == "Lifshitz_93p16":
            self.gen = gen_Lifshitz_93p16(parse, self.out.ltree)
        else:
            print "Invalid generator specified"
            exit()

        #beam effects
        self.beff = None
        if parse.has_section("beam_effects") == True:
            self.beff = beam_effects(parse)

        #tracks in the event
        self.tracks = []

        #run
        nev = parse.getint("main", "nev") # number of events to generate
        print "Number of events:", nev
        self.event_loop(nev)

    #_____________________________________________________________________________
    def event_loop(self, nev):

        iprint = nev/12

        for i in xrange(nev):
            if i%iprint == 0:
                print 100*i/nev, "%"
                stdout.flush()
            self.generate()

        print "All done"

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

        #TX and ROOT outputs
        self.out.write_tx(self.tracks)
        self.out.write_root(self.tracks)

    #_____________________________________________________________________________
    def add_particle(self, part):

        #keep track of particle index in event
        part.idx = len(self.tracks)+1
        self.tracks.append(part)

        return part

 



















