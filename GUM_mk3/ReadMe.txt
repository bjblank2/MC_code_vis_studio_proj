Debugging GUM_mk3:
This code runs a Metropolis MC algorithm to find the ground state magnetic configuration of NiMnIn alloys.
The code works for one implementation of the algorithm and does not for another.
The implemention that I need help with is runMetropolis3().
I have set up main.cpp so that you shouldn't need to change anything. It is set up for a system that SHOULD result in a total
energy of -7.373979425049739 eV/atom... but it doesn't when I use runMetropolis3(). That's the problem.
For debugging, there is one main source files that you will want to look at:
- monte_carlo.cpp
The rest are supporting files that define usefull objects and classes that are used for setting up the simultaion 
and occasionally in it as well.
In monte_carlo.cpp there are three functions that will be relavent:
1) runMetropolis3() ... because that's the one that isn't working
2) evalSiteEnergySpin() which evaluates the energy of a single atom in terms of spin ECI only
3) evalSiteEnergyAll() which evaluates the energy of the simulation cell in terms of all ECI
There is till a LOT of code in these three functions though...