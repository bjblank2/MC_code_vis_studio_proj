Debugging GUM_mk3:
This code runs a Metropolis MC algorithm to find the ground state magnetic configuration of NiMnIn alloys.
The code includes four main files that you will want to look at:
main.cpp - the main function that runs everything
monte_carlo.cpp - the collection of functinos needed to run the monte carlo simulation
rule.cpp - the object container for storage and use instructions for ECI
sim_cell.cpp - the object file that contains the state of the simulation cell
All other files are ether input-output txt files or header files (or file_io.cpp which you can ignore)
 
There are several MC inplementations in monte_carlo.cpp. The implemention that I need help with is runMetropolis3().
I have set up main.cpp so that you shouldn't need to change anything. 
The one source files that you will want to look at the most is:
- monte_carlo.cpp
The rest are supporting files that define usefull objects and classes that are used for setting up the simultaion 
and occasionally in it as well.
In monte_carlo.cpp there are three functions that will be most relavent:
1) runMetropolis3() 
2) evalSiteEnergySpin() which evaluates the energy of a single atom in terms of spin ECI only
3) evalSiteEnergyAll() which evaluates the energy of the simulation cell in terms of all ECI
