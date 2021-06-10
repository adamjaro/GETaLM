# GETaLM - Generator for Electron Tagger and Luminosity Monitor

GETaLM is a set of event generators and related tools for simulations
of luminosity monitor and electron tagger at the EIC.

### Prerequisites

The generator is based on ROOT6 package and Python3. No installation is necessary.

### Running the generator

The executable is followed by the name of the steering card:

<pre><code> ./GETaLM name_of_file </pre></code>

e.g. to generate the bremsstrahlung process the command is

<pre><code> ./GETaLM lumi_18x275.ini </pre></code>

### Directory structure

The GETaLM executable is placed in the top directory. Individual directories are as follows:

- base: base classes for event generation, independent of the physics model
- models: implementations of the individual physics models
- cards: example steering cards

### List of program files

- GETaLM: the generator executable
- init-gen.sh: PATH initialization for bash-type shell
- init-gen.csh: PATH initialization for c-type shell
- base/event.py: manager class for event generation
- base/file_output.py: implementation for generator output
- base/beam_effects.py: effects of beam angular divergence and vertex spread
- base/particle.py: base class for generated particle
- base/beam.py: implementation of particle base for beam particle
- base/photon.py: implementation of particle base for a photon
- models/gen_h1.py: bremsstrahlung generator following H1 parametrization
- models/gen_zeus.py: bremsstrahlung generator following ZEUS parametrization
- models/gen_Lifshitz_93p16.py: bremsstrahlung generator following the QED calculation
- models/gen_quasi_real.py: generator for quasi-real photoproduction
- models/gen_electron_beam.py: utility function to generate beam electrons with effects of beam divergence
- models/gen_read_py.py: reader utility for scattered electrons in Pythia6 events
- models/gen_uniform.py: utility function for electrons or photons with uniform energy and angular distribution
- cards/lumi_18x275.ini: steering card for bremsstrahlung event generation
- cards/quasireal_18x275.ini: steering card for quasi-real photoproduction
- cards/read_py.ini: steering card for the reader for Pythia6 events
- cards/gen_uni.ini: steering card for uniform electrons
- cards/el_beam_18.ini: steering card for beam electrons














