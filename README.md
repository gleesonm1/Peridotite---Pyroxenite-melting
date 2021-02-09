# Peridotite---Pyroxenite-melting
Here I present the codes used in Gleeson et al. (2020; EPSL) and Gleeson and Gibson (2021; Geochem. Geophys. Geosyst.).

All code in this repository is written for MATLAB R2018b.

Presented within this repository are 3 functions:
  1. MeltPX.m - MATLAB version of the MeltPX function presented in Lambart et al. (2016).
  2. dirichletRnd.m - function to generate samples from a Dirichlet distribution. This code was writen by Mo Chen and is critical to the codes presented in this repository.
  3. MeltPXtraceDirichletFe.m - function to calculate the trace element and Fe isotope composition of melts produced by a two component mantle.
  
Alongside these three functions, I present 2 example codes of how they may be used.
  1. xxx - Uses a Dirichlet mixing function to determine the contribution of pyroxenitic melt to a basalt (based on trace element composition) and predicts the Fe isotope composition of this basalt.
  2. LongPlots.m - Uses the 3 functions outlined above to predict how trace element chemistry and crustal thickness will changes along an oceanic spreading centre.
