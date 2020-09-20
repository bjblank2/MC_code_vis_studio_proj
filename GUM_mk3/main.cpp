#include <iostream>
#include <string>
#include "rule.h"
#include "monte_carlo.h"
using namespace std;

SimCell sim_cell;
int main(void) {
	bool use_poscar = true; // if true, read POSCAR file to fill unit cell. The unit cell will be used as the basic building block of the simulation cell
	string input_file = "POSCAR";
	vector<Rule> mc_rules; // Initalize the "Rules" that will be used to properly apply the ECI to each site in the simulation cell
	fillRuleList(mc_rules, "Mart_Rules.txt", 0); // Populate mc_rules. Rules are defined as follows, [species (0,1,2 for Ni,Mn,In)]:[NA if monomer,dist of dymer, AB,AC,BC if trimer]:no-spin(0) or spin(1): value
	int shape[3] = { 4, 4, 4 }; // The dementions of the simulation cell in terms of the unit cell
	vector<int> species{ 512, 512, 0 }; // Number of Ni, Mn, In atoms respectivley 
	cout << shape[0] << ',' << shape[1] << ',' << shape[2] << '\n';
	cout << species[0] << ',' << species[1] << ',' << species[2] << '\n';
	vector<float> dist_list{ 0.0 }; // Initalize list of discances used in mc_rules
	fillDistList(dist_list, mc_rules); // Fill the dist_list
	// Create the simulation cell object. Arguments (POSCAR_file, dist_list, shape, species numbs, cutoff (currently unused), sim_type (also unused), phase_init (aust/mart), spin_init (AFM/FM/RAND), species_init (Ordered/Random), bool use_poscar) 
	// If use_poscar = true, the poscar should reflect the prop of atom types in vector<int> species
	sim_cell.initSimCell(input_file, dist_list, shape, species, 1, string("DEFAULT"), string("MART"), string("RAND"), string("ORDERED"),use_poscar); // Create and initalize the simulation cell
  	cout << "Testing" << '\n';
	cout << "begining MC" << '\n';
	// All of the following MC functinos have the same argument format: (Number of MC moves on each site in the simulation cell per temperature increment, Initial temp, final temp, temp inc, sim_cell object, mc_rules object)
	//runMetropolis1(10, 1, 50, 1, sim_cell, mc_rules); // Hold over code from BEG implementation
	//runMetropolis2(10, 1, 50, 1, mc_rules, mc_rules); // Another BEG implementation
	runMetropolis3(100, 2000, 10, -10, sim_cell, mc_rules); // Run MC using implementation number 3 (The fast one that assumes only spin flips)
	//runMetropolis4(300, 1000, 1, -5, sim_cell, mc_rules); // Run MC using implementation number 4 (Relax with species and spin flips, then re-run using only spin flips)
	//runMetropolis5(200, 1500, 1, -10, sim_cell, mc_rules); // Run MC using implementaion number 5 (Always use spin and species filps)
	//runMetropolis6(300, 1000, 1, -10, sim_cell, mc_rules); // Only for dibugging. Not real MC
	int exit;
	std::cin >> exit;
}