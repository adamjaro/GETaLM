
#_____________________________________________________________________________
#
# Reads scattered electrons and event variables from Pythia6 ascii outputs
#
#_____________________________________________________________________________

import ROOT as rt
from ROOT import gROOT, addressof

from particle import particle

#_____________________________________________________________________________
class gen_read_py:
    #_____________________________________________________________________________
    def __init__(self, parse, tree):

        #open the input
        nam = parse.get("main", "input").strip("\"'")
        print("Input:", nam)
        self.inp = open(nam, "r")

        #skip the header, 6 lines
        for i in range(6): self.inp.readline()

        #tree output for true x, y and Q^2
        tlist = ["true_x", "true_y", "true_Q2", "true_W2", "true_Nu"]
        tlist += ["true_el_pT", "true_el_theta", "true_el_phi", "true_el_E"]
        self.tree_out = self.set_tree(tree, tlist)

        print("gen_read_py initialized")

    #_____________________________________________________________________________
    def generate(self, add_particle):

        lin = ""
        while lin.find("Event finished") < 0:
            lin = self.inp.readline()
            if lin == "": raise IOError("No more events")

            line = lin.split()

            #kinematics from event header
            if line[0] != "I" and line[0] != "I," and len(line) == 30:

                self.tree_out.true_y = float(line[10])
                self.tree_out.true_Q2 = float(line[11])
                self.tree_out.true_x = float(line[12])
                self.tree_out.true_W2 = float(line[13])
                self.tree_out.true_Nu = float(line[14])

            #skip non-particle lines and header
            if len(line) != 14: continue
            if line[0] == "I": continue

            #status, pdg and mother particle
            stat = int(line[1])
            pdg = int(line[2])
            mot = int(line[3])

            #select the scattered electron
            if stat != 1 or pdg != 11 or mot != 3:
                continue

            #print line

            #electron momentum and energy
            px = float(line[6])
            py = float(line[7])
            pz = float(line[8])
            en = float(line[9])

            #generate the electron
            el = add_particle( particle(11) )
            el.vec.SetPxPyPzE(px, py, pz, en)
            el.stat = 1
            el.pxyze_prec = 9

            #electron kinematics directly in output tree
            self.tree_out.true_el_pT = el.vec.Pt()
            self.tree_out.true_el_theta = el.vec.Theta()
            self.tree_out.true_el_phi = el.vec.Phi()
            self.tree_out.true_el_E = el.vec.E()

        #print self.tree_out.true_y, self.tree_out.true_Q2, self.tree_out.true_x, self.tree_out.true_W2, self.tree_out.true_Nu

    #_____________________________________________________________________________
    def set_tree(self, tree, tlist):

        #set output to the tree

        #tree variables
        struct = "struct gen_read_py_out { Double_t "
        for i in tlist:
            struct += i + ", "
        struct = struct[:-2] + ";};"
        gROOT.ProcessLine( struct )
        tree_out = rt.gen_read_py_out()

        #put zero to all variables
        for i in tlist:
            exec("tree_out."+i+"=0")

        #add variables to the tree
        if tree is not None:
            for i in tlist:
                tree.Branch(i, addressof(tree_out, i), i+"/D")

        return tree_out












