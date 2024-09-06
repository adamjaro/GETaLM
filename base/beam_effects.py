
import configparser

import ROOT as rt
from ROOT import TF1, gROOT, addressof

#_____________________________________________________________________________
class beam_effects:
    #angular divergence and emittance
    #_____________________________________________________________________________
    def __init__(self, parse, tree=None, section="beam_effects", hepmc_attrib=None):

        # flag to use or not use the beam effects
        self.use_beam_effects = False
        if parse.has_section(section) == True:
            self.use_beam_effects = parse.getboolean(section, "use_beam_effects")

        print("Beam effects configuration in:", section)
        print("use_beam_effects =", self.use_beam_effects)

        if self.use_beam_effects == False: return

        #beam size at IP in x, sigma_x in mm
        sig_x = parse.getfloat(section, "sig_x")
        print("sig_x =", sig_x)
        self.vtx_x = self.make_gaus("vtx_x", sig_x)

        #beam size at IP in y, sigma_y in mm
        sig_y = parse.getfloat(section, "sig_y")
        print("sig_y =", sig_y)
        self.vtx_y = self.make_gaus("vtx_y", sig_y)

        #bunch length along z
        sig_z = parse.getfloat(section, "sig_z")
        print("sig_z =", sig_z)
        self.vtx_z = self.make_gaus("vtx_z", sig_z)

        #angular divergence in x, horizontal, rad
        theta_x = parse.getfloat(section, "theta_x")
        print("theta_x =", theta_x)
        self.div_x = self.make_gaus("div_x", theta_x)

        #angular divergence in y, vertical, rad
        theta_y = parse.getfloat(section, "theta_y")
        print("theta_y =", theta_y)
        self.div_y = self.make_gaus("div_y", theta_y)

        #tree output from beam effects
        tlist = ["beff_vx", "beff_vy", "beff_vz", "beff_tx", "beff_ty"]
        tlist += ["beff_phot_en", "beff_phot_theta", "beff_phot_phi"]
        tlist += ["beff_phot_px", "beff_phot_py", "beff_phot_pz"]
        tlist += ["beff_el_en", "beff_el_theta", "beff_el_phi"]
        tlist += ["beff_el_px", "beff_el_py", "beff_el_pz"]
        self.tree_out = self.set_tree(tree, tlist)

        #event attributes for hepmc
        self.hepmc_attrib = hepmc_attrib

    #_____________________________________________________________________________
    def apply(self, tracks):
        #apply beam effects

        if not self.use_beam_effects: return

        #beam size in x, y and z
        xpos = self.vtx_x.GetRandom()
        ypos = self.vtx_y.GetRandom()
        zpos = self.vtx_z.GetRandom()

        #angular divergence in x and y
        tx = self.div_x.GetRandom()
        ty = self.div_y.GetRandom()

        #beam size and divergence in output tree
        if self.tree_out is not None:

            self.tree_out.beff_vx = xpos
            self.tree_out.beff_vy = ypos
            self.tree_out.beff_vz = zpos

            self.tree_out.beff_tx = tx
            self.tree_out.beff_ty = ty

        #apply to the final particles
        for i in tracks:
            #select only final particles
            if i.stat != 1: continue

            #vertex position
            i.vx = xpos
            i.vy = ypos
            i.vz = zpos

            #divergence in x by rotation along y
            i.vec.RotateY(tx)

            #divergence in y by rotation along x
            i.vec.RotateX(ty)

            #particle kinematics in output tree
            if self.tree_out is not None:

                #photon
                if i.pdg == 22:

                    self.tree_out.beff_phot_en    = i.vec.Energy()
                    self.tree_out.beff_phot_theta = i.vec.Theta()
                    self.tree_out.beff_phot_phi   = i.vec.Phi()
                    self.tree_out.beff_phot_px = i.vec.Px()
                    self.tree_out.beff_phot_py = i.vec.Py()
                    self.tree_out.beff_phot_pz = i.vec.Pz()

                #electron
                if i.pdg == 11:

                    self.tree_out.beff_el_en    = i.vec.Energy()
                    self.tree_out.beff_el_theta = i.vec.Theta()
                    self.tree_out.beff_el_phi   = i.vec.Phi()
                    self.tree_out.beff_el_px = i.vec.Px()
                    self.tree_out.beff_el_py = i.vec.Py()
                    self.tree_out.beff_el_pz = i.vec.Pz()

            #particle kinematics in hepmc output
            if self.hepmc_attrib is not None:

                #photon
                if i.pdg == 22:

                    self.hepmc_attrib["beff_phot_en"] = i.vec.Energy()
                    self.hepmc_attrib["beff_phot_theta"] = i.vec.Theta()
                    self.hepmc_attrib["beff_phot_phi"] = i.vec.Phi()
                    self.hepmc_attrib["beff_phot_px"] = i.vec.Px()
                    self.hepmc_attrib["beff_phot_py"] = i.vec.Py()
                    self.hepmc_attrib["beff_phot_pz"] = i.vec.Pz()

                #electron
                if i.pdg == 11:

                    self.hepmc_attrib["beff_el_en"] = i.vec.Energy()
                    self.hepmc_attrib["beff_el_theta"] = i.vec.Theta()
                    self.hepmc_attrib["beff_el_phi"] = i.vec.Phi()
                    self.hepmc_attrib["beff_el_px"] = i.vec.Px()
                    self.hepmc_attrib["beff_el_py"] = i.vec.Py()
                    self.hepmc_attrib["beff_el_pz"] = i.vec.Pz()

    #_____________________________________________________________________________
    def make_gaus(self, name, sig):

        gx = TF1(name, "gaus", -12*sig, 12*sig)
        gx.SetParameters(1, 0, sig)

        return gx

    #_____________________________________________________________________________
    def set_tree(self, tree, tlist):

        #set output to the tree if provided
        if tree is None:
            return None

        #tree variables
        struct = "struct beff2_out { Double_t "
        for i in tlist:
            struct += i + ", "
        struct = struct[:-2] + ";};"
        gROOT.ProcessLine( struct )
        tree_out = rt.beff2_out()

        #put zero to all variables
        for i in tlist:
            exec("tree_out."+i+"=0")

        #add variables to the tree
        for i in tlist:
            tree.Branch(i, addressof(tree_out, i), i+"/D")

        return tree_out











