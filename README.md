# Orbital-Collisions
Code to simulate collisions of objects in orbit, based on the NASA Breakup Model (EVOLVE 4.0)

This code was produced as part of an investigation submitted as a paper to Advances in Space Science, by Wynand Joubert and Steven Tingay, from the International Centre for Radio Astronomy Research, Curtin University.  The abstract of the paper, as submitted is:

"In this paper we consider the use of the Murchison Widefield Array (MWA) for rapid response observations of the debris clouds produced by collisions between objects in Earth orbit.  With an increasing density of objects in Low Earth Orbit, the risk of new debris clouds is also increasing.  The MWA constitutes a wide field, rapid response passive radar system and we examine its likely performance in the detection and characterisation of debris clouds.  In order to undertake this work, we adapt the NASA EVOLVE 4.0 breakup model for our purposes, utilising the EVOLVE outputs to produce representative dynamic debris clouds.  We find that the MWA is likely to detect a large fraction (>70%) of debris cloud fragments for collision masses between 100 kg and 1000 kg for orbits in the lower part of LEO, if the MWA can achieve close to optimal detection sensitivity.  Useful detection fractions are still achieved for more conservative assumptions, however.  The detection fraction of fragments decreases as a function of altitude and inversely with collision mass.  Encouragingly, we find that the wide field nature of the MWA allows the full evolving debris clouds to be observed in a single observation, with only $\sim2\%$ of the debris fragments escaping the sensitive portion of the field of view after 100 seconds, for all collision masses and altitudes.  These results show that the MWA is an intrinsically useful facility for rapid characterisation of debris clouds, but work is required to achieve the data processing within a timeframe that is useful to provide rapid alerts."

As the paper is accepted, we will update this README with links to a public pre-print as well as the online journal.

The code base consists of an implementation (in C) of the EVOLVE 4.0 NASA Breakup Model, by Johnson, N.L., Krisko, P.H., Liou, J.-C., and Anz-Meador, P.D., "NASA's new breakup model of evolve 4.0", Advances in Space Research, 2001, #9, pp 1377-1384}: DOI: 10.1016/S0273-1177(01)00423-9, ADS: https://ui.adsabs.harvard.edu/abs/2001AdSpR..28.1377J.

The code base also includes a python wrapper to the C code that produces the plots and analysis that are presented in our paper.  To set up the code:

Download main.c and main.h and compile e.g.: gcc -o collison main.c 

The resultant executable can be run via command line arguments as: collision <charlen> <mass> <output> <numloops>, where:
  charlen: is the lower limit on characteristic length (m)
  mass: is the collision mass (kg)
  output: indicates output to terminal or not (0=terminal; 1=not)
  numloops: number of times to run the simulation
  
If the executable is run without command line arguments, it will prompt the user for the first three arguments above and run a single simulation.

Run the phython script collision.py, after setting the parameters at the top of the script.  collision.py runs the C code and produces all the plots used in the paper with the default settings (however note that the individual realisation of the simulation will be different to that presented in the paper).
