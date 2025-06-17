# Bremsstrahlung example configurations for Phase-1 EIC

## Scope

Generator configurations for bremsstrahlung events in ep and eRu for Phase-1 EIC assumptions.
Beam energies and luminosity values are taken from slides by Elke Aschenauer:

https://agenda.infn.it/event/43344/contributions/250126/attachments/130534/194297/Early.Science.ECA.v2.pptx

All information is given on pages 20 and 21 in the above slides.

## List of files

### brems_ep_10x250_hidiv.ini

ep at 10x250 GeV, total cross section for E_gamma > 0.1 GeV: 226 mb

### brems_eRu_10x115_single.ini

eRu at 10x115 GeV, single bremsstrahlung interaction per simulated event, total cross section for E_gamma > 0.1 GeV: 419 barn

### brems_eRu_10x115_multiple_bx.ini

eRu at 10x115 GeV with multiple bressstrahlung interactions per simulated event, total cross section for E_gamma > 0.1 GeV: 419 barn
(same as with single interactions above). In-bunch pileup as mean number of interactions in bunch crossing is 3.5 on average during
the fill and 9.6 at the beginning of the fill.

