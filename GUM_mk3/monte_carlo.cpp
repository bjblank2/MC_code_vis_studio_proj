#include "monte_carlo.h"
#include <cmath>
#include <random>
extern SimCell sim_cell;

// Fill list of distances from mc_rules
void fillDistList(vector<float> &dist_list, vector<Rule> mc_rules) {
	// Loop through all mc_rules
	for (int i = 0; i < mc_rules.size(); i++) {
		// for two atom terms
		if (mc_rules[i].GetLength() == 2) {
			if (find(dist_list.begin(), dist_list.end(), mc_rules[i].GetDists()[0]) == dist_list.end()) {
				dist_list.push_back(mc_rules[i].GetDists()[0]);
			}
		}
		// for three atom terms
		else if (mc_rules[i].GetLength() == 3) {
			if (find(dist_list.begin(), dist_list.end(), mc_rules[i].GetDists()[0]) == dist_list.end()) {
				dist_list.push_back(mc_rules[i].GetDists()[0]);
			}
			if (find(dist_list.begin(), dist_list.end(), mc_rules[i].GetDists()[1]) == dist_list.end()) {
				dist_list.push_back(mc_rules[i].GetDists()[1]);
			}
			if (find(dist_list.begin(), dist_list.end(), mc_rules[i].GetDists()[2]) == dist_list.end()) {
				dist_list.push_back(mc_rules[i].GetDists()[2]);
			}
		}
	}
}

// Fill list of mc_rules objects from rule file
void fillRuleList(vector<Rule> &list, const char * rule_file, int offset) {
	
	vector<float> distances;
	vector<int> spins;
	vector<int> species;
	float energy_contribution = 0;
	int rule_type = 0;
	int rule_length = 0;
	string phase = "";
	string rule_line;
	vector<string> rule_lines;
	ifstream rule_list;
	rule_list.open(rule_file);
	// Parce rule list txt file
	if (rule_list.is_open()) {
		while (getline(rule_list, rule_line))
		{
			rule_lines.push_back(rule_line);
		}
		rule_list.close();
		for (int i = 0; i < rule_lines.size(); i++) {
			list.push_back(Rule(rule_lines[i]));
		}
	}
	else cout << "Unable to open file";
}

// ONLY USED IN MC1 and MC2 functions! DON'T NEED TO LOOK AT THESE
float evalSiteEnergy1(float temp, int site, SimCell &sim_cell, vector<Rule> &mc_rules) {
	float Kb = 0.000086173324;
	float uB = .000057883818012;
	float H = 0;
	float site_energy = 0;
	int site_species = sim_cell.atom_list[site].getSpecies();
	int site_spin = sim_cell.atom_list[site].getSpin();
	int site_phase = sim_cell.atom_list[site].getPhase();
//	site_energy -= Kb * temp * log(8)*(1 - pow(site_phase, 2));
	site_energy -= 3 * uB*H*site_spin;
	for (int rule = 0; rule < mc_rules.size(); rule++) {
		if (mc_rules[rule].GetLength() == 1) {
			if (site_species == mc_rules[rule].GetSpecies()[0]) {
				site_energy += mc_rules[rule].GetEnrgCont();
				//sim_cell.clust_count[rule] += 1;
			}
		}
		else if (mc_rules[rule].GetLength() == 2) {
			int neighbor_spin;
			int neighbor_species;
			for (int neighbor = 0; neighbor < sim_cell.atom_list[site].getNumbNeighbors(site, sim_cell); neighbor++) {
				neighbor_spin = sim_cell.atom_list[site].getNeighborSpin(neighbor, sim_cell);
				neighbor_species = sim_cell.atom_list[site].getNeighborSpecies(neighbor, sim_cell);
				if (compf(sim_cell.atom_list[site].getNeighborDist(neighbor, sim_cell), mc_rules[rule].GetDists()[0])) {
					vector<int> pair1{ site_species,neighbor_species };
					vector<int> pair2{ neighbor_species,site_species };
					if (pair1 == mc_rules[rule].GetSpecies() || pair2 == mc_rules[rule].GetSpecies()) {
						if (mc_rules[rule].GetType() == 0) {
							site_energy += 0.5 * mc_rules[rule].GetEnrgCont();
							//sim_cell.clust_count[rule] += 1;
						}
						else if (mc_rules[rule].GetType() == 1) {
							site_energy += 0.5 * mc_rules[rule].GetEnrgCont() * site_spin * neighbor_spin;
							//sim_cell.clust_count[rule] += site_spin * neighbor_spin;

						}
					}
				}
			}

		}
		//else if (mc_rules[rule].GetLength() == 3) {
		//	int neighbor1_spin;
		//	int neighbor1_species;
		//	int neighbor2_spin;
		//	int neighbor2_species;
		//	for (int neighbor1 = 0; neighbor1 < sim_cell.atom_list[site].getNumbNeighbors(site, sim_cell); neighbor1++) {
		//		neighbor1_spin = sim_cell.atom_list[site].getNeighborSpin(neighbor1, sim_cell);
		//		neighbor1_species = sim_cell.atom_list[site].getNeighborSpecies(neighbor1, sim_cell);
		//		for (int neighbor2 = 0; neighbor2 < sim_cell.atom_list[site].getNumbNeighbors(site, sim_cell); neighbor2++) {
		//			neighbor2_spin = sim_cell.atom_list[site].getNeighborSpin(neighbor2, sim_cell);
		//			neighbor2_species = sim_cell.atom_list[site].getNeighborSpecies(neighbor2, sim_cell);
		//			float dist_ab = sim_cell.atom_list[site].getNeighborDist(neighbor1, sim_cell);
		//			float dist_ac = sim_cell.atom_list[site].getNeighborDist(neighbor2, sim_cell);
		//			float dist_bc = sim_cell.findAtomDists(sim_cell.atom_list[site].getNeighborIndex(neighbor1, sim_cell), sim_cell.atom_list[site].getNeighborIndex(neighbor2, sim_cell));
		//			int atom_a = site_species;
		//			int atom_b = neighbor1_species;
		//			int atom_c = neighbor2_species;
		//			bool cluster_match = false;
		//			vector<int> trip{ atom_a,atom_b,atom_c };
		//			vector<float> dists{ dist_ab,dist_bc,dist_ac };
		//			vector<float> rule_dists = mc_rules[rule].GetDists();
		//			if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
		//				cluster_match = true;
		//			}
		//			trip = { atom_a,atom_c,atom_b };
		//			dists = { dist_ac,dist_bc,dist_ab };
		//			if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
		//				cluster_match = true;
		//			}
		//			trip = { atom_b,atom_a,atom_c };
		//			dists = { dist_ab,dist_ac,dist_bc };
		//			if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
		//				cluster_match = true;
		//			}
		//			trip = { atom_b,atom_c,atom_a };
		//			dists = { dist_bc,dist_ac,dist_ab };
		//			if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
		//				cluster_match = true;
		//			}
		//			trip = { atom_c,atom_a,atom_b };
		//			dists = { dist_ac,dist_ab,dist_bc };
		//			if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
		//				cluster_match = true;
		//			}
		//			trip = { atom_c,atom_b,atom_a };
		//			dists = { dist_bc,dist_ab,dist_ac };
		//			if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
		//				cluster_match = true;
		//			}
		//			if (cluster_match == true) {
		//				if (mc_rules[rule].GetType() == 0) {
		//					site_energy += 0.5 * mc_rules[rule].GetEnrgCont();
		//					//sim_cell.clust_count[rule] += 1;
		//				}
		//				else if (mc_rules[rule].GetType() == 1) {
		//					site_energy += 0.5 * mc_rules[rule].GetEnrgCont() * site_spin * neighbor1_spin * neighbor2_spin;
		//					//sim_cell.clust_count[rule] += site_spin * neighbor1_spin * neighbor2_spin;
		//				}
		//			}
		//		}
		//	}
		//}
	}
	return site_energy;
}
// ONLY USED IN MC1 and MC2 functions! DON'T NEED TO LOOK AT THESE
float evalSiteEnergy2(float temp, int site, SimCell &sim_cell, vector<Rule> &mc_rules) {
	float Kb = 0.000086173324;
	float uB = .000057883818012;
	float H = 0;
	float site_energy = 0;
	int numb_neighbors;
	int site_species = sim_cell.atom_list[site].getSpecies();
	int site_spin = sim_cell.atom_list[site].getSpin();
	int site_phase = sim_cell.atom_list[site].getPhase();
		site_energy -= Kb * temp * log(8)*(1 - pow(site_phase, 2));
	site_energy -= 3 * uB*H*site_spin;
	
	for (int rule = 0; rule < mc_rules.size(); rule++) {
		if (find(mc_rules[rule].GetSpecies().begin(), mc_rules[rule].GetSpecies().end(), site_species) != mc_rules[rule].GetSpecies().end()) {
			if (mc_rules[rule].GetLength() == 1) {
					site_energy += mc_rules[rule].GetEnrgCont();
					//sim_cell.clust_count[rule] += 1;
			}
			else if (mc_rules[rule].GetLength() == 2) {
				int neighbor_spin;
				int neighbor_species;
				numb_neighbors = sim_cell.atom_list[site].getNumbNeighbors(site, sim_cell);
				for (int neighbor = 0; neighbor < numb_neighbors; neighbor++) {
					neighbor_species = sim_cell.atom_list[site].getNeighborSpecies(neighbor, sim_cell);
					if (find(mc_rules[rule].GetSpecies().begin(), mc_rules[rule].GetSpecies().end(), neighbor_species) != mc_rules[rule].GetSpecies().end()) {
						if (compf(sim_cell.atom_list[site].getNeighborDist(neighbor, sim_cell), mc_rules[rule].GetDists()[0])) {
							if (mc_rules[rule].GetType() == 0) {
								site_energy += 0.5 * mc_rules[rule].GetEnrgCont();
								//sim_cell.clust_count[rule] += 1 
							}
							else if (mc_rules[rule].GetType() == 1) {
								neighbor_spin = sim_cell.atom_list[site].getNeighborSpin(neighbor, sim_cell);
								site_energy += 0.5 * mc_rules[rule].GetEnrgCont() * site_spin * neighbor_spin;
								//sim_cell.clust_count[rule] += site_spin * neighbor_spin;

							}
						}
					}
				}

			}
			else if (mc_rules[rule].GetLength() == 3) {
				int neighbor1_spin;
				int neighbor1_species;
				int neighbor2_spin;
				int neighbor2_species;
				for (int neighbor1 = 0; neighbor1 < sim_cell.atom_list[site].getNumbNeighbors(site, sim_cell); neighbor1++) {
					neighbor1_spin = sim_cell.atom_list[site].getNeighborSpin(neighbor1, sim_cell);
					neighbor1_species = sim_cell.atom_list[site].getNeighborSpecies(neighbor1, sim_cell);
					for (int neighbor2 = 0; neighbor2 < sim_cell.atom_list[site].getNumbNeighbors(site, sim_cell); neighbor2++) {
						neighbor2_spin = sim_cell.atom_list[site].getNeighborSpin(neighbor2, sim_cell);
						neighbor2_species = sim_cell.atom_list[site].getNeighborSpecies(neighbor2, sim_cell);
						float dist_ab = sim_cell.atom_list[site].getNeighborDist(neighbor1, sim_cell);
						float dist_ac = sim_cell.atom_list[site].getNeighborDist(neighbor2, sim_cell);
						float dist_bc = sim_cell.findAtomDists(sim_cell.atom_list[site].getNeighborIndex(neighbor1, sim_cell), sim_cell.atom_list[site].getNeighborIndex(neighbor2, sim_cell));
						int atom_a = site_species;
						int atom_b = neighbor1_species;
						int atom_c = neighbor2_species;
						bool cluster_match = false;
						vector<int> trip{ atom_a,atom_b,atom_c };
						vector<float> dists{ dist_ab,dist_bc,dist_ac };
						vector<float> rule_dists = mc_rules[rule].GetDists();
						if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
							cluster_match = true;
						}
						trip = { atom_a,atom_c,atom_b };
						dists = { dist_ac,dist_bc,dist_ab };
						if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
							cluster_match = true;
						}
						trip = { atom_b,atom_a,atom_c };
						dists = { dist_ab,dist_ac,dist_bc };
						if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
							cluster_match = true;
						}
						trip = { atom_b,atom_c,atom_a };
						dists = { dist_bc,dist_ac,dist_ab };
						if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
							cluster_match = true;
						}
						trip = { atom_c,atom_a,atom_b };
						dists = { dist_ac,dist_ab,dist_bc };
						if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
							cluster_match = true;
						}
						trip = { atom_c,atom_b,atom_a };
						dists = { dist_bc,dist_ab,dist_ac };
						if (compv(dists, rule_dists) && trip == mc_rules[rule].GetSpecies()) {
							cluster_match = true;
						}
						if (cluster_match == true) {
							if (mc_rules[rule].GetType() == 0) {
								site_energy += 0.5 * mc_rules[rule].GetEnrgCont();
								//sim_cell.clust_count[rule] += 1;
							}
							else if (mc_rules[rule].GetType() == 1) {
								site_energy += 0.5 * mc_rules[rule].GetEnrgCont() * site_spin * neighbor1_spin * neighbor2_spin;
								//sim_cell.clust_count[rule] += site_spin * neighbor1_spin * neighbor2_spin;
							}
						}
					}
				}
			}
		}
	}
	return site_energy;
}

// Evaluate the energy controbution from a singal atom, includeing spin chem rules.
float evalSiteEnergyAll(float temp, int site, map<string, float> &rule_map_spin, map<string, float> &rule_map_chem, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	// Initalize variables
	float Kb = 0.00008617333262;
	float uB = .000057883818012;
	float H = 0;
	string key;
	float enrg = 0;
	map<string, float>::iterator rule_itr;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;
	int neighbor_site2;
	int neighbor_species2;
	int neighbor_spin2;
	int dist_2_ind;
	float neighbor_dist2;
	float neighbor_dist3;
	enrg -= 3 * uB*H*site_spin; // add energy from external magnetic field (if there is one... right now its hard coded to 0)
	key = "_" + to_string(site_species) + ",0,"; // make key for one atom term
	rule_itr = rule_map_chem.find(key);
	enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // apply rule for one atom term
	// for 2 atom terms
	// loop through all neighbors in the neighbor list (list only includes atoms that have distances listed in the dist_list (the list of distances used in a rule))
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		neighbor_site1 = neighbor_index_list[site][i];
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ","; // make key string
		rule_itr = rule_map_chem.find(key);
		enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // if the key was in the chem rule map, apply the rule, otherwise do nothing
		rule_itr = rule_map_spin.find(key);
		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * site_spin * neighbor_spin1 : 0.0; // if the key was in the spin rule map, apply the rule, otherwise do nothing
		// for 2 atom terms
		// loop through all neighbors in the neighbor list
		for (int j = 0; j < neighbor_index_list[site].size(); j++) {
			if (i != j) {
				neighbor_site2 = neighbor_index_list[site][j];
				neighbor_species2 = atom_species[neighbor_site2];
				neighbor_spin2 = atom_spin[neighbor_site2];
				vector<int>::iterator dist_2_itr = find(neighbor_index_list[neighbor_site1].begin(), neighbor_index_list[neighbor_site1].end(), neighbor_site2);
				if (dist_2_itr != neighbor_index_list[neighbor_site1].end()) {
					dist_2_ind = distance(neighbor_index_list[neighbor_site1].begin(), dist_2_itr);
					neighbor_dist2 = neighbor_dist_list[neighbor_site1][dist_2_ind];
					neighbor_dist3 = neighbor_dist_list[site][j];
					// make key for three atom term
					key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_species2);
					key +="," + to_string(neighbor_dist1) + "," + to_string(neighbor_dist2) + "," + to_string(neighbor_dist3) + ",";
					rule_itr = rule_map_chem.find(key);
					enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second / 6 : 0.0; // apply rule (if possible). Devide by six because of tripple count degeneracy of three atom terms
				}
			}
		}
	}
	return enrg;
}

// Evaluate the energy contribution from a singal atom using only spin rules. (Assuming the spin is all that changed)
float evalSiteEnergySpin(float temp, int site, map<string, float> &rule_map_spin, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	// Initalize variables and constants
	float Kb = 0.00008617333262;
	float uB = .000057883818012;
	float H = 0;
	string key;
	float enrg = 0;
	map<string, float>::iterator rule_itr;
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	enrg -= 3 * uB*H*site_spin; // apply magnetic field (if any, right now its hard coded to zero)
	// for 2 atom terms
	// Loop through all neighbors in the neighbor list (neighbor list only includes distances included in the dist_list (the list of distances used in a rule))
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		// gather info to build key
		neighbor_site1 = neighbor_index_list[site][i];
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ","; // create key string
		rule_itr = rule_map_spin.find(key);
  		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * site_spin * neighbor_spin1 : 0.0; // if key is found in rule list, apply the rule, else do nothing.
	}
	return enrg;
}

// calculate the contribution to the spin-product order paramiter from a single site
float calcMag2(int site, vector<int> &atom_spin, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list){
	int site_spin = atom_spin[site];
	float mag_product = 0;
	// loop through all neighbors in neighbor list
	for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
		if (neighbor_dist_list[site][i] == .5  or neighbor_dist_list[site][i] == 0.75) {
			mag_product += site_spin * atom_spin[neighbor_index_list[site][i]]; // if neighbor is directly above/below or directly adjacent... do the product
		}
	}
	return mag_product / 6; // devide by 6 for normalization
}

// claculate the stag mag order paramiter contribution from a single atom
float stagMag(int site, int spin, SimCell& sim_cell) {
	int mag = 0;
	int check = fmod(sim_cell.atom_list[site].pos[0] * sim_cell.atom_list[site].pos[1] * 2,2); // check to see what sublattice the site is on
	if (check == 0) { mag += spin; }
	else { mag += 0.0; }
	return mag;
}

// Curently unused!!
float calcSpecies(int site, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	int site_species = 0;
	int neighbor_species = 0;
	if (atom_species[site] == 1) { site_species = 1; }
	else if (atom_species[site] == 2) { site_species = -1; }
	float species_product = 0;
	for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
		if (neighbor_dist_list[site][i] == .5 or neighbor_dist_list[site][i] == 0.75) {
			if (atom_species[neighbor_index_list[site][i]] == 1) { neighbor_species = 1; }
			else if (atom_species[neighbor_index_list[site][i]] == 2) { neighbor_species = -1; }
			species_product += site_species * neighbor_species;
		}
	}
	return species_product / 6;
}
// Curently unused!!
float delSiteEnergySpin(float temp, int site, int old_spin, map<string, float> &rule_map_spin, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	float Kb = 0.000086173324;
	float uB = .000057883818012;
	float H = 0;
	string key;
	float enrg = 0;
	map<string, float>::iterator rule_itr;
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	//	site_energy -= Kb * temp * log(8)*(1 - pow(site_phase, 2));
	enrg -= 3 * uB*H*site_spin;
	// for 2 atom terms
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		neighbor_site1 = neighbor_index_list[site][i];
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ",";
		rule_itr = rule_map_spin.find(key);
		if (site_spin == 0) {
			enrg -= (rule_itr != rule_map_spin.end()) ? rule_itr->second * old_spin * neighbor_spin1 : 0.0;
		}
		else {
			enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * 2 * site_spin * neighbor_spin1 : 0.0;
		}
	}
	return enrg;
}
// Overloaded def of evalLattice that is only used for BEG. Don't need to look at it
float evalLattice(float temp, SimCell &sim_cell, vector<Rule> &mc_rules) {
	float e_total = 0;
	for (int site = 0; site < sim_cell.numb_atoms; site++) {
		e_total += evalSiteEnergy1(temp, site, sim_cell, mc_rules);
	}
	return e_total;// / sim_cell.numb_atoms * 16;// - 66.8295069760671;
}

// Evaluate the total energy of the simulation cell in it's current state
float evalLattice(float temp, map<string, float> &rule_map_spin, map<string, float> &rule_map_chem, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	float enrg = 0;
	// loop thruogh simulation cell and add everything up
	for (int site = 0; site < atom_species.size(); site++) {
		enrg += evalSiteEnergyAll(temp, site, rule_map_spin, rule_map_chem, atom_spin, atom_species, neighbor_index_list, neighbor_dist_list);
	}
	return enrg;
}

// Evaluate the spin contribution to the total energy of the simulation cell
float evalLatticeSpin(float temp, map<string, float> &rule_map_spin, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	float enrg = 0;
	// loop through everything and add it up
	for (int site = 0; site < atom_species.size(); site++) {
		enrg += evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin, atom_species, neighbor_index_list, neighbor_dist_list);
	}
	return enrg;// / atom_species.size() * 16 - 66.8295069760671;
}

// MC for BEG only. Don't need to look at it (and it's prety ugly so you probably shouldn't)
void runMetropolis1(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules) {
	float Kb = .0000861733035;
	float e_total = 0;
	float e_site_old = 0;
	float e_site_new = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int numb_atoms = sim_cell.numb_atoms;
	int flip_count = 0;
	int flip_count2 = 0;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	float initial_enrg = evalLattice(temp1, sim_cell, mc_rules);
	cout << initial_enrg;
	cout << "\n";
	for (float temp = temp1; temp < temp2; temp += temp_inc) {
		e_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Spin
				e_site_old = evalSiteEnergy1(temp, site, sim_cell, mc_rules);
				old_spin = sim_cell.atom_list[site].getSpin();
				spin_same = true;
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) {
						new_spin = -1;
					}
					else if (spin_rand <= 0.6666666666666666) {
						new_spin = 0;
					}
					else {
						new_spin = 1;
					}
					if (new_spin != old_spin) { spin_same = false; }
				}
				sim_cell.atom_list[site].setSpin(new_spin);
				e_site_new = evalSiteEnergy1(temp, site, sim_cell, mc_rules);
				if (e_site_new <= e_site_old) {
					e_total += e_site_new;
					flip_count2 += 1;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand < keep_prob) {
						e_total += e_site_new;
						flip_count += 1;
					}
					else {
						sim_cell.atom_list[site].setSpin(old_spin);
						e_total += e_site_old;
					}
				}
				current_spin = sim_cell.atom_list[site].getSpin();
				spin_total2 += current_spin;
				if (sim_cell.atom_list[site].getSpecies() != 0) {
					for (int neighbors = 0; neighbors < 6; neighbors++) {
						spin_total += sim_cell.atom_list[site].getSpin() * sim_cell.atom_list[site].getNeighborSpin(neighbors, sim_cell);
					}
				}
			}
			spin_avg += spin_total;
			spin_avg2 += spin_total2;
			e_avg += e_total;
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / passes / numb_atoms * 16;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms / 6 * 2;
		cout << " , ";
		cout << spin_avg2 / passes / numb_atoms;
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}
// MC for BEG only. Don't need to look at it (ugly here too. save yourself the trouble and don't look!)
void runMetropolis2(float passes, float temp1, float temp2, float temp_inc, vector<Rule> &mc_rules) {
	float Kb = .0000861733035;
	float e_total = 0;
	float e_site_old = 0;
	float e_site_new = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int numb_atoms = sim_cell.numb_atoms;
	vector<SimCell::Atom> &atom_list = sim_cell.atom_list;
	int flip_count = 0;
	int flip_count2 = 0;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	float initial_enrg = evalLattice(temp1, sim_cell, mc_rules);
	cout << initial_enrg;
	cout << "\n";
	for (float temp = temp1; temp < temp2; temp += temp_inc) {
		e_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Spin
				e_site_old = evalSiteEnergy1(temp, site, sim_cell, mc_rules);
				old_spin = atom_list[site].getSpin();
				spin_same = true;
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) {
						new_spin = -1;
					}
					else if (spin_rand <= 0.6666666666666666) {
						new_spin = 0;
					}
					else {
						new_spin = 1;
					}
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_list[site].setSpin(new_spin);
				e_site_new = evalSiteEnergy1(temp, site, sim_cell, mc_rules);
				if (e_site_new <= e_site_old) {
					e_total += e_site_new;
					flip_count2 += 1;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand < keep_prob) {
						e_total += e_site_new;
						flip_count += 1;
					}
					else {
						atom_list[site].setSpin(old_spin);
						e_total += e_site_old;
					}
				}
				current_spin = atom_list[site].getSpin();
				spin_total2 += current_spin;
				if (atom_list[site].getSpecies() != 0) {
					for (int neighbors = 0; neighbors < 6; neighbors++) {
						spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(neighbors, sim_cell);
					}
				}
			}
			spin_avg += spin_total;
			spin_avg2 += spin_total2;
			e_avg += e_total;
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / passes / numb_atoms * 16;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms / 6 * 2;
		cout << " , ";
		cout << spin_avg2 / passes / numb_atoms;
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}

// Run MC with ONLY spin flips
void runMetropolis3(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules) {
	// Pull some basic information out of the sim_cell object (an object mostly used for setup that contains inforamtion on the
	// position, spin, atomic species, and neighbor list for each atom in the simulation cell)
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, float> rule_map_spin;
	map <string, float> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into a map for spin and a map for chemestry rules to make lookup faster.
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	// Loop through all mc_rules
	for (int i = 0; i < mc_rules.size(); i++) {
		string rule_key = "_"; // add starting charicter to key string for parcing purposes.
		// If rule is single atom term...
		if (mc_rules[i].GetLength() == 1) {
			rule_key = "_" + to_string(mc_rules[i].GetSpecies()[0]) + ",0,"; // make key string
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));// assign key atom pair 
		}
		// If rule is two atom term...
		else if (mc_rules[i].GetLength() == 2) {
			vector<int> species = mc_rules[i].GetSpecies();
			float dist = mc_rules[i].GetDists()[0];
			// make degenerate keys for [0,1]->[1,0] type rules
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont())); // assign key atom pair if chem type
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont())); // assign key atom pair if spin type
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont())); // assign key atom pair if chem type
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont())); // assign key atom pair if spin type
			}
		}
		// if rule is three atom term...
		else if (mc_rules[i].GetLength() == 3) {
			vector<int> trip = mc_rules[i].GetSpecies();
			vector<float> dists = mc_rules[i].GetDists();
			// make keys and assign key value pairs for all degenerate three atom ruels, (always chem rules)
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
	}
	// Extract species list, spin list, neighbor index list and neighbor dist list from sim_cell object to make them more accesable latter.
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	// initalize variables 
	float Kb = 0.00008617333262;
	float e_total = 0;
	float e_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	int spin_avg_Mn = 0;
	int spin_total_Mn = 0;
	float spin_avg_Mn2 = 0;
	float spin_total_Mn2 = 0;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	// Initalize random number generator. This part might give your computer some problems if it is a Mac FYI !!!!!!!!!!
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// Calculate the total energy of the initalized simulation cell
	float initial_enrg = evalLattice(temp1, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	// Calculate the total energy contribution form magnetic ECI
	float initial_spin_cont = evalLatticeSpin(temp1, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms << "\n";
	cout << "Temp, energy, spin avg, spin avg Mn only, spin product, favorable flips, unfavorable flips \n";
	// Enter main MC loop
	//Begin decreasing temperature
	for (int temp = temp1; temp >= temp2; temp+=temp_inc) {
		e_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		spin_avg_Mn = 0;
		spin_avg_Mn2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		// loop through the simulation cell specified number of times
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			spin_total_Mn = 0;
			spin_total_Mn2 = 0;
			// loop through every atom in the cell
			for (int site = 0; site < numb_atoms; site++) {
				// Flip a spin
				old_spin = atom_spin_list[site]; // get pre-flip spin
				e_site = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list); // get pre-flip energy
				// pick a new spin at random (don't want to waist time picking same spin that site already has)
				spin_same = true;
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) {
						new_spin = -1;
					}
					else if (spin_rand <= 0.6666666666666666) {
						new_spin = 0;
					}
					else {
						new_spin = 1;
					}
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin; // flip the spin
				e_site_new = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list); // find energy after the flip
				e_flip = (e_site_new - e_site); // find net change in energy
				// apply chose keep/throw out using Metropolis
				if (e_flip < 0) {
					flip_count2 += 1;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_flip));
					if (keep_rand < keep_prob) {
						flip_count += 1;
					}
					else {
						atom_spin_list[site] = old_spin;
						e_flip = 0;
					}
				}
				initial_enrg += e_flip; // update the total energy
				e_total += initial_enrg / numb_atoms; // normalize
				current_spin = atom_spin_list[site]; 
				spin_total += current_spin; // get spin total
				if (atom_species_list[site] == 1) {
					spin_total_Mn += current_spin; // Mn only spin totla
					spin_total_Mn2 += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list); // calc spin product order param
					//spin_total_Mn2 += stagMag(site, atom_spin_list[site], sim_cell); // calc stag mag order param
				}
			}
			spin_avg += spin_total;
			spin_avg_Mn += spin_total_Mn;
			spin_avg_Mn2 += spin_total_Mn2;
			e_avg += e_total / passes;
		}
		// Output normaized avarages of tracked quantities
		cout << temp;
		cout << " , ";
		cout << e_avg / numb_atoms;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg_Mn / passes / sim_cell.species_numbs[1];
		cout << " , ";  
		cout << spin_avg_Mn2 / passes / sim_cell.species_numbs[1];
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
	writeSuperCell(atom_species_list, atom_spin_list, sim_cell); // write the filnal state of the simulation cell to txt file (for debugging and sanity checking mostly)
}

// This is not the function I am testing with. You can look through it if you want but its not as important 
void runMetropolis4(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules) {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, float> rule_map_spin;
	map <string, float> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < mc_rules.size(); i++) {
		string rule_key = "_";
		if (mc_rules[i].GetLength() == 1) {
			rule_key = "_" + to_string(mc_rules[i].GetSpecies()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
		else if (mc_rules[i].GetLength() == 2) {
			vector<int> species = mc_rules[i].GetSpecies();
			float dist = mc_rules[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
		}
		else if (mc_rules[i].GetLength() == 3) {
			vector<int> trip = mc_rules[i].GetSpecies();
			vector<float> dists = mc_rules[i].GetDists();
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	float Kb = .0000861733035;
	float e_total = 0;
	float e_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	int spin_avg_Mn = 0;
	int spin_total_Mn = 0;
	float spin_avg_Mn2 = 0;
	float spin_total_Mn2 = 0;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// Find ground state configuration
	int temp_inc_relax = -10;
	for (int temp = temp1; temp >= temp2; temp += temp_inc_relax) {
		cout << temp << "\n";
		if (temp <= 20) {
			temp_inc_relax = -1;
		}
		for (int i = 0; i < passes; i++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						e_site = evalSiteEnergyAll(temp, site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_site += evalSiteEnergyAll(temp, rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						e_site_new = evalSiteEnergyAll(temp, site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_site_new += evalSiteEnergyAll(temp, rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_flip = e_site_new - e_site;
						if (e_flip < 0) {}
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb*temp)*(e_flip));
							if (keep_rand < keep_prob) {}
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
				// Flip Spin
				//e_site_old = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				if (e_flip < 0) {}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_flip));
					if (keep_rand < keep_prob) {}
					else { atom_spin_list[site] = old_spin; }
				}
			}
		}
	}
	// begin actual MC run //
	writeSuperCell(atom_species_list, atom_spin_list, sim_cell);
	float initial_enrg = evalLattice(temp1, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(temp1, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + -4.176806541934177 << ", " << initial_spin_cont / numb_atoms << "\n";
	for (int temp = temp1; temp >= temp2; temp += temp_inc) {
		e_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		spin_avg_Mn = 0;
		spin_avg_Mn2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			spin_total_Mn = 0;
			spin_total_Mn2 = 0;
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Spin
				//e_site_old = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				if (e_flip < 0) { flip_count += 1; }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; }
					else {
						atom_spin_list[site] = old_spin;
						e_flip = 0;
					}
				}
				initial_enrg += e_flip;
				e_total += initial_enrg / numb_atoms;
				current_spin = atom_spin_list[site];
				spin_total += current_spin;
				if (atom_species_list[site] == 1) {
					spin_total_Mn += current_spin;
					spin_total_Mn2 += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list);
				}
			}
			spin_avg += spin_total;
			spin_avg_Mn += spin_total_Mn;
			spin_avg_Mn2 += spin_total_Mn2;
			e_avg += evalLattice(temp1, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list) / (passes);        
			//e_avg += e_total / (passes);
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / numb_atoms;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg_Mn / passes / sim_cell.species_numbs[1];
		cout << " , ";
		cout << spin_avg_Mn2 / passes / sim_cell.species_numbs[1];
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}

// This is not the function I am testing with. You can look through it if you want but its not as important 
void runMetropolis5(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules) {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, float> rule_map_spin;
	map <string, float> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < mc_rules.size(); i++) {
		string rule_key = "_";
		if (mc_rules[i].GetLength() == 1) {
			rule_key = "_" + to_string(mc_rules[i].GetSpecies()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
		else if (mc_rules[i].GetLength() == 2) {
			vector<int> species = mc_rules[i].GetSpecies();
			float dist = mc_rules[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
		}
		else if (mc_rules[i].GetLength() == 3) {
			vector<int> trip = mc_rules[i].GetSpecies();
			vector<float> dists = mc_rules[i].GetDists();
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	float Kb = .0000861733035;
	float e_total = 0;
	float e_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	int spin_avg_Mn = 0;
	int spin_total_Mn = 0;
	float spin_avg_Mn2 = 0;
	float spin_total_Mn2 = 0;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// begin actual MC run //
	float initial_enrg = evalLattice(temp1, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(temp1, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + -4.918174088267291 << ", " << initial_spin_cont / numb_atoms << "\n";
	for (int temp = temp1; temp >= temp2; temp += temp_inc) {
		e_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		spin_avg_Mn = 0;
		spin_avg_Mn2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			spin_total_Mn = 0;
			spin_total_Mn2 = 0;
			float order_enrgy_running = 0;
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						e_site = evalSiteEnergyAll(temp, site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_site += evalSiteEnergyAll(temp, rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						e_site_new = evalSiteEnergyAll(temp, site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_site_new += evalSiteEnergyAll(temp, rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_flip = e_site_new - e_site;
						if (e_flip < 0) {}
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb*temp)*(e_flip));
							if (keep_rand < keep_prob) {}
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
				// Flip Spin
				//e_site_old = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				if (e_flip < 0) { flip_count += 1; }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; }
					else {
						atom_spin_list[site] = old_spin;
						e_flip = 0;
					}
				}
				current_spin = atom_spin_list[site];
				spin_total += current_spin;
				if (atom_species_list[site] == 1) {
					spin_total_Mn += current_spin;
					spin_total_Mn2 += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list);
				}
			}
			spin_avg += spin_total;
			spin_avg_Mn += spin_total_Mn;
			spin_avg_Mn2 += spin_total_Mn2;
			e_avg += evalLattice(temp1, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list) / (passes);
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / numb_atoms; // -4.918174088267291;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg_Mn / passes / sim_cell.species_numbs[1];
		cout << " , ";
		cout << spin_avg_Mn2 / passes / sim_cell.species_numbs[1];
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
	writeSuperCell(atom_species_list, atom_spin_list, sim_cell);
}

// Not actual MC. Only for debugging
void runMetropolis6(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules) {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, float> rule_map_spin;
	map <string, float> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < mc_rules.size(); i++) {
		string rule_key = "_";
		if (mc_rules[i].GetLength() == 1) {
			rule_key = "_" + to_string(mc_rules[i].GetSpecies()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
		else if (mc_rules[i].GetLength() == 2) {
			vector<int> species = mc_rules[i].GetSpecies();
			float dist = mc_rules[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
		}
		else if (mc_rules[i].GetLength() == 3) {
			vector<int> trip = mc_rules[i].GetSpecies();
			vector<float> dists = mc_rules[i].GetDists();
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	float Kb = .0000861733035;
	float e_total = 0;
	float e_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	int spin_avg_Mn = 0;
	int spin_total_Mn = 0;
	float spin_avg_Mn2 = 0;
	float spin_total_Mn2 = 0;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// begin actual MC run //
	float initial_enrg = evalLattice(temp1, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(temp1, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms << ", " << initial_spin_cont / numb_atoms << "\n";
}

// Write the current state of the simulation cell to a txt file. Mostly for debugging.
void writeSuperCell(vector<int> &atom_species, vector<int>& atom_spins, SimCell &sim_cell) {
	ofstream OUT_file;
	OUT_file.open("OUTCAR.txt");
	if (OUT_file.is_open()) {
		OUT_file << "Ni2Mn(1-X)In(X)\n" << " 5.76 \n";
		OUT_file << sim_cell.sup_cell[0] <<" 0 0\n" << "0 "<< sim_cell.sup_cell[1] << " 0\n" << "0 0 " << sim_cell.sup_cell[2] * 1.5 << "\n";
		OUT_file << "Ni Mn In\n";
		OUT_file << sim_cell.species_numbs[0] << " " << sim_cell.species_numbs[1] << " " << sim_cell.species_numbs[2] << "\n";
		for (int i = 0; i < sim_cell.numb_atoms; i++) {
			OUT_file << sim_cell.atom_list[i].pos[0] << " " << sim_cell.atom_list[i].pos[1] << " " << sim_cell.atom_list[i].pos[2] << " " << atom_species[i] << " " << atom_spins[i] << "\n";
		}
	}
}

// not used
bool compf(float x1, float x2, float eps) {
	return (fabs( x1- x2) < eps);
}
// not used
bool compv(vector<float> &x1, vector<float> &x2, float eps) {
	bool equal = true;
	for (int i = 0; i < x1.size(); i++) {
		if (fabs(x1[i] - x2[i]) > eps) {
			equal = false;
		}
	}
	return equal;
}

//
//
//
//void clacBEGParams(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
//	int home_spin = atom_list[site].getSpin();
//	int home_phase = atom_list[site].getPhase();
//	int home_species = atom_list[site].getSpecies();
//	int neighbor_spin;
//	int neighbor_phase;
//	int neighbor_species;
//	int neighbor_order;
//	string neighbor_plain;
//	J_K[0] = 0;
//	J_K[1] = 0;
//	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
//		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
//		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
//		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
//		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
//		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
//		for (int i = 0; i < cluster_rules.size(); i++) {
//			if (neighbor_order == cluster_rules[i].getOrder()) {
//				if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
//					if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
//						if (neighbor_plain == cluster_rules[i].getPlain() || cluster_rules[i].getPlain() == "ALL") {
//							if (cluster_rules[i].getNeighborArrangment() == "PERM") {
//								if (home_species != neighbor_species) {
//									if (cluster_rules[i].getPhase() == 1) {
//										J_K[0] += cluster_rules[i].getEnergyContribution();
//									}
//									if (cluster_rules[i].getPhase() == 0) {
//										J_K[1] += cluster_rules[i].getEnergyContribution();
//									}
//								}
//							}
//							if (cluster_rules[i].getNeighborArrangment() == "COMB") {
//								if (cluster_rules[i].getPhase() == 1) {
//									J_K[0] += cluster_rules[i].getEnergyContribution();
//								}
//								if (cluster_rules[i].getPhase() == 0) {
//									J_K[1] += cluster_rules[i].getEnergyContribution();
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		for (int i = 0; i < spin_rules.size(); i++) {
//			if (neighbor_order == spin_rules[i].getOrder()) {
//				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
//					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
//						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
//							if (spin_rules[i].getNeighborArrangment() == "PERM") {
//								if (home_species != neighbor_species) {
//									if (spin_rules[i].getPhase() == 1) {
//										J_K[0] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									}
//									if (spin_rules[i].getPhase() == 0) {
//										J_K[1] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									}
//								}
//							}
//							if (spin_rules[i].getNeighborArrangment() == "COMB") {
//								if (spin_rules[i].getPhase() == 1) {
//									J_K[0] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//								if (spin_rules[i].getPhase() == 0) {
//									J_K[1] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	J_K[0] -= .0;
//	J_K[1] -= .0;
//}

//
//void init_calcJK(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
//	int home_spin = atom_list[site].getSpin();
//	int home_phase = atom_list[site].getPhase();
//	int home_species = atom_list[site].getSpecies();
//	int neighbor_spin;
//	int neighbor_phase;
//	int neighbor_species;
//	int neighbor_order;
//	string neighbor_plain;
//	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
//		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
//		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
//		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
//		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
//		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
//		for (int i = 0; i < cluster_rules.size(); i++) {
//			if (neighbor_order == cluster_rules[i].getOrder()) {
//				if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
//					if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
//						if (neighbor_plain == cluster_rules[i].getPlain() || cluster_rules[i].getPlain() == "ALL") {
//							if (cluster_rules[i].getNeighborArrangment() == "PERM") {
//								if (home_species != neighbor_species) {
//									if (cluster_rules[i].getPhase() == 1) {
//										atom_list[site].J += cluster_rules[i].getEnergyContribution();
//									}
//									if (cluster_rules[i].getPhase() == 0) {
//										atom_list[site].K += cluster_rules[i].getEnergyContribution();
//									}
//								}
//							}
//							if (cluster_rules[i].getNeighborArrangment() == "COMB") {
//								if (cluster_rules[i].getPhase() == 1) {
//									atom_list[site].J += cluster_rules[i].getEnergyContribution();
//								}
//								if (cluster_rules[i].getPhase() == 0) {
//									atom_list[site].K += cluster_rules[i].getEnergyContribution();
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		for (int i = 0; i < spin_rules.size(); i++) {
//			if (neighbor_order == spin_rules[i].getOrder()) {
//				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
//					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
//						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
//							if (spin_rules[i].getNeighborArrangment() == "PERM") {
//								if (home_species != neighbor_species) {
//									if (spin_rules[i].getPhase() == 1) {
//										atom_list[site].J += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									}
//									if (spin_rules[i].getPhase() == 0) {
//										atom_list[site].K += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									}
//								}
//							}
//							if (spin_rules[i].getNeighborArrangment() == "COMB") {
//								if (spin_rules[i].getPhase() == 1) {
//									atom_list[site].J += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//								if (spin_rules[i].getPhase() == 0) {
//									atom_list[site].K += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//void pair_calcJK(int site, int neighbor, vector<Atom> &atom_list, vector<float> J_K) {
//	J_K[0] = (atom_list[site].J + atom_list[neighbor].J) / 2;
//	J_K[1] = (atom_list[site].K + atom_list[neighbor].K) / 2;
//}
//
//void re_calcJK(int site, int old_home_spin, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
//	atom_list[site].J = 0;
//	atom_list[site].K = 0;
//	int home_spin = atom_list[site].getSpin();
//	int home_phase = atom_list[site].getPhase();
//	int home_species = atom_list[site].getSpecies();
//	int neighbor_index;
//	int neighbor_spin;
//	int neighbor_phase;
//	int neighbor_species;
//	int neighbor_order;
//	string neighbor_plain;
//	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
//		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
//		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
//		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
//		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
//		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
//		neighbor_index = atom_list[site].neighbors[neighbor];
//		for (int i = 0; i < cluster_rules.size(); i++) {
//			if (neighbor_order == cluster_rules[i].getOrder()) {
//				if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
//					if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
//						if (neighbor_plain == cluster_rules[i].getPlain() || cluster_rules[i].getPlain() == "ALL") {
//							if (cluster_rules[i].getNeighborArrangment() == "PERM") {
//								if (home_species != neighbor_species) {
//									if (cluster_rules[i].getPhase() == 1) {
//										atom_list[site].J += cluster_rules[i].getEnergyContribution();
//									}
//									if (cluster_rules[i].getPhase() == 0) {
//										atom_list[site].K += cluster_rules[i].getEnergyContribution();
//									}
//								}
//							}
//							if (cluster_rules[i].getNeighborArrangment() == "COMB") {
//								if (cluster_rules[i].getPhase() == 1) {
//									atom_list[site].J += cluster_rules[i].getEnergyContribution();
//								}
//								if (cluster_rules[i].getPhase() == 0) {
//									atom_list[site].K += cluster_rules[i].getEnergyContribution();
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		for (int i = 0; i < spin_rules.size(); i++) {
//			if (neighbor_order == spin_rules[i].getOrder()) {
//				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
//					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
//						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
//							if (spin_rules[i].getNeighborArrangment() == "PERM") {
//								if (home_species != neighbor_species) {
//									if (spin_rules[i].getPhase() == 1) {
//										atom_list[site].J += -1*spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//										atom_list[neighbor_index].J -= spin_rules[i].getEnergyContribution()*old_home_spin*neighbor_spin;
//										atom_list[neighbor_index].J += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									}
//									if (spin_rules[i].getPhase() == 0) {
//										atom_list[site].K += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//										atom_list[neighbor_index].K -= spin_rules[i].getEnergyContribution()*old_home_spin*neighbor_spin;
//										atom_list[neighbor_index].K += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									}
//								}
//							}
//							if (spin_rules[i].getNeighborArrangment() == "COMB") {
//								if (spin_rules[i].getPhase() == 1) {
//									atom_list[site].J += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									atom_list[neighbor_index].J -= spin_rules[i].getEnergyContribution()*old_home_spin*neighbor_spin;
//									atom_list[neighbor_index].J += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//								if (spin_rules[i].getPhase() == 0) {
//									atom_list[site].K += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									atom_list[neighbor_index].K -= spin_rules[i].getEnergyContribution()*old_home_spin*neighbor_spin;
//									atom_list[neighbor_index].K += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//void clacBEGParams(vector<float> &J_K) {
//	J_K[0] = -0.002; //-0.822868 , -0.771559
//	J_K[1] = -0.026;
//}
//
//
//void clacBEGParamsNEW(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
//	int home_phase = atom_list[site].getPhase();
//	int home_species = atom_list[site].getSpecies();
//	int neighbor_spin;
//	int neighbor_phase;
//	int neighbor_species;
//	int neighbor_order;
//	string neighbor_plain;
//	J_K[0] = 0;
//	J_K[1] = 0;
//	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
//		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
//		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
//		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
//		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
//		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
//		for (int i = 0; i < cluster_rules.size(); i++) {
//			if (neighbor_order == cluster_rules[i].getOrder()) {
//				if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
//					if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
//						if (neighbor_plain == cluster_rules[i].getPlain() || cluster_rules[i].getPlain() == "ALL") {
//							if (cluster_rules[i].getNeighborArrangment() == "PERM") {
//								if (home_species != neighbor_species) {
//									if (cluster_rules[i].getPhase() == 1) {
//										J_K[0] += cluster_rules[i].getEnergyContribution();
//									}
//									if (cluster_rules[i].getPhase() == 0) {
//										J_K[1] += cluster_rules[i].getEnergyContribution();
//									}
//								}
//							}
//							if (cluster_rules[i].getNeighborArrangment() == "COMB") {
//								if (cluster_rules[i].getPhase() == 1) {
//									J_K[0] += cluster_rules[i].getEnergyContribution();
//								}
//								if (cluster_rules[i].getPhase() == 0) {
//									J_K[1] += cluster_rules[i].getEnergyContribution();
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	J_K[0] -= .0;
//	J_K[1] -= .0;
//}
//
//float evalSiteEnergySINGLE(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
//	float Kb = 0.000086173324;
//	float uB = .000057883818012;
//	float H = 0;
//	float site_energy = 0;
//	int site_phase = atom_list[site].getPhase();
//	int neighbor_phase;
//	int site_spin = atom_list[site].getSpin();
//	int sig1;
//	int sig2;
//	// select wether to use fixed or on the fly J-K calcuations 
//	//calcBEGParams(J_K); // Fixed J-K
//	clacBEGParamsNEW(site, atom_list, cluster_rules, spin_rules, J_K);  // on the fly J-K
//	for (int neighbor = 0; neighbor < 8; neighbor++) {
//		neighbor_phase = atom_list[site].getNeighborPhase(1, neighbor, atom_list);
//		sig1 = 1 - pow(site_phase, 2);
//		sig2 = 1 - pow(neighbor_phase, 2);
//		site_energy += J_K[0] * site_phase*neighbor_phase + J_K[1] * sig1*sig2;
//	}
//	site_energy /= 8; ////////////////////////////////////////////////////////////////////////// AAAAAAAAAAAAAAAHHHHHHHHH !!!!!!!!!! ////////////
//	site_energy -= Kb * temp * log(8)*(1 - pow(site_phase, 2));
//	site_energy -= 3 * uB*H*site_spin;
//	// add mag contribution
//	int home_spin = atom_list[site].getSpin();
//	int home_phase = atom_list[site].getPhase();
//	int home_species = atom_list[site].getSpecies();
//	int neighbor_spin;
//	int neighbor_species;
//	int neighbor_order;
//	string neighbor_plain;
//	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
//		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
//		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
//		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
//		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
//		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
//		for (int i = 0; i < spin_rules.size(); i++) {
//			if (neighbor_order == spin_rules[i].getOrder()) {
//				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
//					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
//						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
//							if (spin_rules[i].getNeighborArrangment() == "PERM") {
//								if (home_species != neighbor_species) {
//									if (abs(home_phase) == 1) {
//										if (abs(spin_rules[i].getPhase()) == 1){
//											site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//										}
//									}
//									else if (home_phase == 0) {
//										if (spin_rules[i].getPhase() == 0) {
//											site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//										}
//									}
//								}
//							}
//							else if (spin_rules[i].getNeighborArrangment() == "COMB") {
//								if (spin_rules[i].getPhase() == 1 == abs(home_phase)) {
//									site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//								else if (spin_rules[i].getPhase() == 0 == home_phase) {
//									site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	return site_energy;
//}
//
//float evalSiteEnergyTOTAL(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
//	float site_energy = evalSiteEnergySINGLE(temp, site, atom_list, cluster_rules, spin_rules, J_K);
//	for (int neighbor = 0; neighbor < 26; neighbor++) {
//		int n_site = atom_list[site].neighbors[neighbor];
//		site_energy += evalSiteEnergySINGLE(temp, n_site, atom_list, cluster_rules, spin_rules, J_K);
//	}
//	int home_spin = atom_list[site].getSpin();
//	int home_phase = atom_list[site].getPhase();
//	int home_species = atom_list[site].getSpecies();
//	int neighbor_spin;
//	int neighbor_species;
//	int neighbor_order;
//	int neighbor_phase;
//	string neighbor_plain;
//	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
//		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
//		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
//		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
//		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
//		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
//		for (int i = 0; i < spin_rules.size(); i++) {
//			if (neighbor_order == spin_rules[i].getOrder()) {
//				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
//					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
//						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
//							if (spin_rules[i].getNeighborArrangment() == "PERM") {
//								if (home_species != neighbor_species) {
//									if (spin_rules[i].getPhase() == 1 == abs(home_phase)) {
//										site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									}
//									else if (spin_rules[i].getPhase() == 0 == home_phase) {
//										site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//									}
//								}
//							}
//							if (spin_rules[i].getNeighborArrangment() == "COMB") {
//								if (spin_rules[i].getPhase() == 1 == abs(home_phase)) {
//									site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//								else if (spin_rules[i].getPhase() == 0 == home_phase) {
//									site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	return site_energy;
//}
//
//float evalSiteEnergy3(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
//	float Kb = 0.000086173324;
//	float uB = .000057883818012;
//	float H = 0;
//	float site_energy = 0;
//	int site_phase = atom_list[site].getPhase();
//	int neighbor_phase;
//	int site_spin = atom_list[site].getSpin();
//	int sig1;
//	int sig2;
//	// select wether to use fixed or on the fly J-K calcuations 
//	//clacBEGParams(J_K); // Fixed J-K
//	clacBEGParams(site, atom_list, cluster_rules, spin_rules, J_K);  // on the fly J-K
//	for (int neighbor = 0; neighbor < 8; neighbor++) {
//		neighbor_phase = atom_list[site].getNeighborPhase(1, neighbor, atom_list);
//		sig1 = 1 - pow(site_phase, 2);
//		sig2 = 1 - pow(neighbor_phase, 2);
//		site_energy += J_K[0] * site_phase*neighbor_phase + J_K[1] * sig1*sig2;
//	}
//	site_energy /= 8; ////////////////////////////////////////////////////////////////////////// AAAAAAAAAAAAAAAHHHHHHHHH !!!!!!!!!! ////////////
//	site_energy -= Kb * temp * log(2)*(1 - pow(site_phase, 2));
//	site_energy -= 3 * uB*H*site_spin;
//	// add mag contribution
//	return site_energy;
//}
//
//float evalSiteEnergy4(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
//	float site_energy = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
//	for (int neighbor = 0; neighbor < 26; neighbor++) {
//		int n_site = atom_list[site].neighbors[neighbor];
//		site_energy += evalSiteEnergy3(temp, n_site, atom_list, cluster_rules, spin_rules, J_K);
//	}
//	return site_energy / 27;
//}
//
//float evalSiteEnergy5(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
//	float Kb = .0000861733035;
//	int site_phase = atom_list[site].getPhase();
//	int site_spin = atom_list[site].getSpin();
//	int site_species = atom_list[site].getSpecies();
//	int neighbor_phase;
//	int neighbor_spin;
//	int neighbor_species;
//	int neighbor_order;
//	string neighbor_plain;
//
//	float current_enrg = 0;
//	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
//		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
//		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
//		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
//		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
//		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
//		for (int clust_rule = 0; clust_rule < cluster_rules.size(); clust_rule++) {
//			if (neighbor_order == cluster_rules[clust_rule].getOrder()) {
//				if (find(cluster_rules[clust_rule].home_species.begin(), cluster_rules[clust_rule].home_species.end(), site_species) != cluster_rules[clust_rule].home_species.end()) {
//					if (find(cluster_rules[clust_rule].neighbor_species.begin(), cluster_rules[clust_rule].neighbor_species.end(), neighbor_species) != cluster_rules[clust_rule].neighbor_species.end()) {
//						if (neighbor_plain == cluster_rules[clust_rule].getPlain() || cluster_rules[clust_rule].getPlain() == "ALL") {
//							if (cluster_rules[clust_rule].getNeighborArrangment() == "PERM") {
//								if (site_species != neighbor_species) {
//									if (cluster_rules[clust_rule].getPhase() == 1) {
//										current_enrg += cluster_rules[clust_rule].getEnergyContribution()*site_phase*neighbor_phase;
//									}
//									else if (cluster_rules[clust_rule].getPhase() == 0) {
//										current_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
//									}
//								}
//							}
//							else if (cluster_rules[clust_rule].getNeighborArrangment() == "COMB") {
//								if (cluster_rules[clust_rule].getPhase() == 1) {
//									current_enrg += cluster_rules[clust_rule].getEnergyContribution()*site_phase*neighbor_phase;
//								}
//								else if (cluster_rules[clust_rule].getPhase() == 0) {
//									current_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		for (int spin_rule = 0; spin_rule < spin_rules.size(); spin_rule++) {
//			if (neighbor_order == spin_rules[spin_rule].getOrder()) {
//				if (find(spin_rules[spin_rule].home_species.begin(), spin_rules[spin_rule].home_species.end(), site_species) != spin_rules[spin_rule].home_species.end()) {
//					if (find(spin_rules[spin_rule].neighbor_species.begin(), spin_rules[spin_rule].neighbor_species.end(), neighbor_species) != spin_rules[spin_rule].neighbor_species.end()) {
//						if (neighbor_plain == spin_rules[spin_rule].getPlain() || spin_rules[spin_rule].getPlain() == "ALL") {
//							if (spin_rules[spin_rule].getNeighborArrangment() == "PERM") {
//								if (site_species != neighbor_species) {
//									if (spin_rules[spin_rule].getPhase() == 1 and abs(site_phase) == 1) {
//										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//									}
//									if (spin_rules[spin_rule].getPhase() == 0 and site_phase == 0) {
//										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//									}
//									if (spin_rules[spin_rule].getPhase() == 1 and abs(neighbor_phase) == 1) {
//										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//									}
//									if (spin_rules[spin_rule].getPhase() == 0 and neighbor_phase == 0) {
//										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//									}
//								}
//							}
//							else if (spin_rules[spin_rule].getNeighborArrangment() == "COMB") {
//								if (spin_rules[spin_rule].getPhase() == 1 and abs(site_phase) == 1) {
//									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//								}
//								if (spin_rules[spin_rule].getPhase() == 0 and site_phase == 0) {
//									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//								}
//								if (spin_rules[spin_rule].getPhase() == 1 and abs(neighbor_phase) == 1) {
//									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//								}
//								if (spin_rules[spin_rule].getPhase() == 0 and neighbor_phase == 0) {
//									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	current_enrg += Kb * temp*log(8)*(1 - pow(site_phase, 2));
//	return current_enrg;
//}
//
//float evalSiteEnergy6(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
//	float Kb = 0.000086173324;
//	float uB = .000057883818012;
//	float H = 0;
//	float site_energy = 0;
//	int site_phase = atom_list[site].getPhase();
//	int neighbor_phase;
//	int neighbor_site;
//	int site_spin = atom_list[site].getSpin();
//	float avg_J;
//	float avg_K;
//	int sig1;
//	int sig2;
//	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
//		//neighbor_site = atom_list[site].getNeighbor(1, neighbor, atom_list);
//		//neighbor_phase = atom_list[site].getNeighborPhase(1, neighbor, atom_list);
//		neighbor_site = atom_list[site].neighbors[neighbor];
//		neighbor_phase = atom_list[neighbor_site].getPhase();
//		avg_J = (atom_list[site].J + atom_list[neighbor_site].J) / 2;
//		avg_K = (atom_list[site].K + atom_list[neighbor_site].K) / 2;
//		sig1 = 1 - pow(site_phase, 2);
//		sig2 = 1 - pow(neighbor_phase, 2);
//		site_energy += avg_J * site_phase*neighbor_phase + avg_K * sig1*sig2;
//	}
//	site_energy /=  atom_list[site].getNumbNeighbors();
//	site_energy -= Kb * temp * log(2)*(1 - pow(site_phase, 2));
//	site_energy -= 3 * uB*H*site_spin;
//	// add mag contribution
//	return site_energy;
//}
//
//void eval_flip(float temp, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, int site, int new_state[3], vector<float> &flip_enrgs) {
//	float Kb = .0000861733035;
//	int site_phase = atom_list[site].getPhase();
//	int site_spin = atom_list[site].getSpin();
//	int site_species = atom_list[site].getSpecies();
//	int new_phase = new_state[0];
//	int new_spin = new_state[1];
//	int new_species = new_state[2];
//	int neighbor_phase;
//	int neighbor_spin;
//	int neighbor_species;
//	int neighbor_order;
//	string neighbor_plain;
//
//	float current_enrg = 0;
//	float new_enrg = 0;
//	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
//		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
//		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
//		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
//		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
//		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
//		for (int clust_rule = 0; clust_rule < cluster_rules.size(); clust_rule++) {
//			if (neighbor_order == cluster_rules[clust_rule].getOrder()) {
//				if (find(cluster_rules[clust_rule].home_species.begin(), cluster_rules[clust_rule].home_species.end(), site_species) != cluster_rules[clust_rule].home_species.end()) {
//					if (find(cluster_rules[clust_rule].neighbor_species.begin(), cluster_rules[clust_rule].neighbor_species.end(), neighbor_species) != cluster_rules[clust_rule].neighbor_species.end()) {
//						if (neighbor_plain == cluster_rules[clust_rule].getPlain() || cluster_rules[clust_rule].getPlain() == "ALL") {
//							if (cluster_rules[clust_rule].getNeighborArrangment() == "PERM") {
//								if (site_species != neighbor_species) {
//									if (cluster_rules[clust_rule].getPhase() == 1) {
//										current_enrg += cluster_rules[clust_rule].getEnergyContribution()*site_phase*neighbor_phase;
//										new_enrg += cluster_rules[clust_rule].getEnergyContribution()*new_phase*neighbor_phase;
//									}
//									else if (cluster_rules[clust_rule].getPhase() == 0) {
//										current_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
//										new_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(new_phase, 2))*(1 - pow(neighbor_phase, 2));
//									}
//								}
//							}
//							else if (cluster_rules[clust_rule].getNeighborArrangment() == "COMB") {
//								if (cluster_rules[clust_rule].getPhase() == 1) {
//									current_enrg += cluster_rules[clust_rule].getEnergyContribution()*site_phase*neighbor_phase;
//									new_enrg += cluster_rules[clust_rule].getEnergyContribution()*new_phase*neighbor_phase;
//								}
//								else if (cluster_rules[clust_rule].getPhase() == 0) {
//									current_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
//									new_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(new_phase, 2))*(1 - pow(neighbor_phase, 2));
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		for (int spin_rule = 0; spin_rule < spin_rules.size(); spin_rule++) {
//			if (neighbor_order == spin_rules[spin_rule].getOrder()) {
//				if (find(spin_rules[spin_rule].home_species.begin(), spin_rules[spin_rule].home_species.end(), site_species) != spin_rules[spin_rule].home_species.end()) {
//					if (find(spin_rules[spin_rule].neighbor_species.begin(), spin_rules[spin_rule].neighbor_species.end(), neighbor_species) != spin_rules[spin_rule].neighbor_species.end()) {
//						if (neighbor_plain == spin_rules[spin_rule].getPlain() || spin_rules[spin_rule].getPlain() == "ALL") {
//							if (spin_rules[spin_rule].getNeighborArrangment() == "PERM") {
//								if (site_species != neighbor_species) {
//									if (spin_rules[spin_rule].getPhase() == 1 and abs(site_phase) == 1) {
//										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//									}
//									if (spin_rules[spin_rule].getPhase() == 1 and abs(new_phase) == 1) {
//										new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin / 2;
//									}
//									if (spin_rules[spin_rule].getPhase() == 0 and site_phase == 0) {
//										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//									}
//									if (spin_rules[spin_rule].getPhase() == 0 and new_phase == 0) {
//										new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin / 2;
//									}
//									if (spin_rules[spin_rule].getPhase() == 1 and abs(neighbor_phase) == 1) {
//										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//										new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin / 2;
//									}
//									if (spin_rules[spin_rule].getPhase() == 0 and neighbor_phase == 0) {
//										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//										new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin / 2;
//									}
//								}
//							}
//							else if (spin_rules[spin_rule].getNeighborArrangment() == "COMB") {
//								if (spin_rules[spin_rule].getPhase() == 1 and abs(site_phase) == 1) {
//									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//								}
//								if (spin_rules[spin_rule].getPhase() == 1 and abs(new_phase) == 1) {
//									new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin / 2;
//								}
//								if (spin_rules[spin_rule].getPhase() == 0 and site_phase == 0) {
//									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//								}
//								if (spin_rules[spin_rule].getPhase() == 0 and new_phase == 0) {
//									new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin / 2;
//								}
//								if (spin_rules[spin_rule].getPhase() == 1 and abs(neighbor_phase) == 1) {
//									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//									new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin / 2;
//								}
//								if (spin_rules[spin_rule].getPhase() == 0 and neighbor_phase == 0) {
//									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
//									new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin / 2;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	current_enrg += Kb * temp*log(8)*(1 - pow(site_phase, 2));
//	new_enrg += Kb * temp*log(8)*(1 - pow(site_phase, 2));
//	flip_enrgs[0] = current_enrg;
//	flip_enrgs[1] = new_enrg;
//}
//
//
//float evalCluster(vector<Atom> &atom_list, Cluster &cluster, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K, float temp) {
//	float Kb = .000086173324;
//	float total_H = 0;
//	int site;
//	float total_H_inc;
//	float inc_count = 0;
//	int site_phase;
//	int neighbor_phase;
//	int proposed_link = 0;
//	if (cluster.clusterSize() == 1) {
//		site = cluster.cluster_list[0];
//		total_H = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
//	}
//	else {
//		for (int i = 0; i < cluster.clusterSize(); i++) {
//			site = cluster.cluster_list[i];
//			site_phase = atom_list[site].getPhase();
//			clacBEGParams(site, atom_list, cluster_rules, spin_rules, J_K);
//			total_H_inc = 0;
//			inc_count = 0;
//			for (int neighbor = 0; neighbor < 8; neighbor++) {
//				proposed_link = atom_list[site].getNeighbor(1, neighbor, atom_list);
//				if (cluster.inList(proposed_link)) {
//					neighbor_phase = atom_list[site].getNeighborPhase(1, neighbor, atom_list);
//					total_H_inc += J_K[0] * site_phase*neighbor_phase + J_K[1] * (1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
//					inc_count += 1;
//				}
//			}
//			total_H += (total_H_inc / inc_count + Kb * temp*log(2)*pow(site_phase,2));
//		}
//	}
//	return total_H;
//}
//
//void flipCluster(int seed_phase, int new_phase, vector<Atom> &atom_list, Cluster &cluster, bool reset) {
//	int site;
//	int old_phase;
//	if (reset == false) {
//		if (seed_phase*new_phase == -1) {
//			for (int cluster_site = 0; cluster_site < cluster.clusterSize(); cluster_site++) {
//				site = cluster.cluster_list[cluster_site];
//				atom_list[site].setPhase(new_phase);
//			}
//		}
//		else {
//			for (int cluster_site = 0; cluster_site < cluster.clusterSize(); cluster_site++) {
//				site = cluster.cluster_list[cluster_site];
//				if ((seed_phase == 1 and new_phase == 0) or (seed_phase == 0 and new_phase == -1)) {
//					old_phase = atom_list[site].getPhase();
//					if (old_phase == 1) {
//						atom_list[site].setPhase(0);
//					}
//					else if (old_phase == 0) {
//						atom_list[site].setPhase(-1);
//					}
//				}
//				if ((seed_phase == -1 and new_phase == 0) or (seed_phase == 0 and new_phase == 1)) {
//					old_phase = atom_list[site].getPhase();
//					if (old_phase == -1) {
//						atom_list[site].setPhase(0);
//					}
//					else if (old_phase == 0) {
//						atom_list[site].setPhase(1);
//					}
//				}
//			}
//		}
//	}
//	else {
//		for (int cluster_site = 0; cluster_site < cluster.clusterSize(); cluster_site++) {
//			site = cluster.cluster_list[cluster_site];
//			if ((seed_phase == 1 and new_phase == 0) or (seed_phase == 0 and new_phase == -1)) {
//				old_phase = atom_list[site].getPhase();
//				if (old_phase == 0) {
//					atom_list[site].setPhase(1);
//				}
//				else if (old_phase == -1) {
//					atom_list[site].setPhase(0);
//				}
//			}
//			if ((seed_phase == -1 and new_phase == 0) or (seed_phase == 0 and new_phase == 1)) {
//				old_phase = atom_list[site].getPhase();
//				if (old_phase == 0) {
//					atom_list[site].setPhase(-1);
//				}
//				else if (old_phase == 1) {
//					atom_list[site].setPhase(0);
//				}
//			}
//		}
//	}
//}
////////////////////////////////////////////////////////////////////////////////////////////////
//// Befor you run a MC algorithm be sure to check if it is configured for heating or cooling //
////////////////////////////////////////////////////////////////////////////////////////////////
//void runMetropolis1(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
//	float Kb = .0000861733035;
//	float e_total = 0;
//	float e_site_old = 0;
//	float e_site_new = 0;
//	float spin_rand = 0;
//	float phase_rand = 0;
//	float keep_rand = 0;
//	int old_phase = 0;
//	int new_phase = 0;
//	int old_spin = 0;
//	int new_spin = 0;
//	int current_spin = 0;
//	int current_phase = 0;
//	bool phase_same;
//	bool spin_same;
//	float e_avg = 0;
//	float spin_avg = 0;
//	float spin_total = 0;
//	float spin_avg2 = 0;
//	float spin_total2 = 0;
//	float phase_total = 0;
//	float phase_avg = 0;
//	float phase_total_abs = 0;
//	float phase_avg_abs = 0;
//	float keep_prob = 0;
//	int numb_atoms = size(atom_list);
//	float current_J;
//	float current_K;
//	float new_J;
//	float new_K;
//	float atom_avg_J;
//	float atom_avg_K;
//	float pass_avg_J;
//	float pass_avg_K;
//	int flip_count = 0;
//	int flip_count2 = 0;
//	vector<float> J_K = { 0,0 };
//	std::mt19937_64 rng;
//	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
//	rng.seed(ss);
//	std::uniform_real_distribution<double> unif(0, 1);
//	cout << evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K);
//	cout << "\n";
//	for (float temp = temp1; temp < temp2; temp += temp_inc) {
//		e_avg = 0;
//		phase_avg = 0;
//		phase_avg_abs = 0;
//		spin_avg = 0;
//		spin_avg2 = 0;
//		pass_avg_J = 0;
//		pass_avg_K = 0;
//		flip_count = 0;
//		flip_count2 = 0;
//		for (int i = 0; i < passes; i++) {
//			e_total = 0;
//			phase_total = 0;
//			phase_total_abs = 0;
//
//			spin_total = 0;
//			spin_total2 = 0;
//			atom_avg_J = 0;
//			atom_avg_K = 0;
//			for (int site = 0; site < atom_list.size(); site++) {
//				// Flip Phase
//				bool keep = false;
//				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				//e_site_old = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);////////////////////////          //////
//				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				//old_phase = atom_list[site].getPhase();
//				//phase_same = true;
//				//while (phase_same == true) {
//				//	phase_rand = unif(rng);
//				//	if (phase_rand <= 0.3333333333333333) {
//				//		new_phase = -1;
//				//	}
//				//	else if (phase_rand <= 0.6666666666666666) {
//				//		new_phase = 0;
//				//	}
//				//	else {
//				//		new_phase = 1;
//				//	}
//				//	if (new_phase != old_phase) { phase_same = false; }
//				//}
//				//atom_list[site].setPhase(new_phase);
//				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				//e_site_new = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);////////////////////////         ///////
//				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				//if (e_site_new <= e_site_old) {
//				//	flip_count2 += 1;
//				//	keep = true;
//				//	e_total += e_site_new;
//				//}
//				//else {
//				//	keep_rand = unif(rng);
//				//	keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
//				//	if (keep_rand <= keep_prob) {
//				//		keep = true;
//				//		e_total += e_site_new;
//				//		flip_count += 1;
//				//	}
//				//	else {
//				//		atom_list[site].setPhase(old_phase);
//				//		keep = false;
//				//		e_total += e_site_old;
//				//	}
//				//}
//				//current_phase = atom_list[site].getPhase();
//				//phase_total += current_phase;
//				//phase_total_abs += abs(current_phase);
//
//				//atom_avg_J += J_K[0];
//				//atom_avg_K += J_K[1];
//
//				// Flip Spin
//				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				e_site_old = evalSiteISING(temp, site, atom_list, cluster_rules, spin_rules, J_K);                   ////////////////
//				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				current_J = J_K[0];
//				current_K = J_K[1];
//				old_spin = atom_list[site].getSpin();
//				spin_same = true;
//				while (spin_same == true) {
//					spin_rand = unif(rng);
//					if (spin_rand <= 0.3333333333333333) {
//						new_spin = -1;
//					}
//					else if (spin_rand <= 0.6666666666666666) {
//						new_spin = 0;
//					}
//					else {
//						new_spin = 1;
//					}
//					if (new_spin != old_spin) { spin_same = false; }
//				}
//				atom_list[site].setSpin(new_spin);
//				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				e_site_new = evalSiteISING(temp, site, atom_list, cluster_rules, spin_rules, J_K);                     //////////////
//				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				new_J = J_K[0];
//				new_K = J_K[1];
//				if (e_site_new <= e_site_old) {
//					e_total += e_site_new;
//					flip_count2 += 1;
//					keep = true;
//				}
//				else {
//					keep_rand = unif(rng);
//					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
//					if (keep_rand < keep_prob) {
//						e_total += e_site_new;
//						flip_count += 1;
//						keep = true;
//					}
//					else {
//						atom_list[site].setSpin(old_spin);
//						e_total += e_site_old;
//						keep = false;
//					}
//				}
//				current_spin = atom_list[site].getSpin();
//				spin_total2 += current_spin;
//				if (atom_list[site].getSpecies() != 0) {
//					for (int neighbors = 0; neighbors < 6; neighbors++) {
//						spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
//					}
//				}
//			}
//			phase_avg += phase_total;
//			phase_avg_abs += phase_total_abs;
//			spin_avg += spin_total;
//			spin_avg2 += spin_total2;
//			e_avg += e_total;
//			pass_avg_J += atom_avg_J;
//			pass_avg_K += atom_avg_K;
//		}
//		cout << temp;
//		cout << " , ";
//		cout << e_avg / passes / numb_atoms * 16;
//		cout << " , ";
//		cout << phase_avg / passes / numb_atoms;
//		cout << " , ";
//		cout << phase_avg_abs / passes / numb_atoms;
//		cout << " , ";
//		cout << spin_avg / passes / numb_atoms / 6 * 2;
//		cout << " , ";
//		cout << spin_avg2 / passes / numb_atoms;
//		cout << " , ";
//		cout << pass_avg_J / passes / numb_atoms;
//		cout << " , ";
//		cout << pass_avg_K / passes / numb_atoms;
//		cout << " , ";
//		cout << pass_avg_K / pass_avg_J;
//		cout << " , ";
//		cout << flip_count;
//		cout << " , ";
//		cout << flip_count2;
//		cout << "\n";
//	}
//}
//
//void runMetropolis2(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
//	float Kb = .0000861733035;
//	float e_total = 0;
//	float e_site_old = 0;
//	float e_site_new = 0;
//	float spin_rand = 0;
//	float phase_rand = 0;
//	float keep_rand = 0;
//	int old_phase = 0;
//	int new_phase = 0;
//	int old_spin = 0;
//	int new_spin = 0;
//	int current_spin = 0;
//	int current_phase = 0;
//	bool both_same;
//	float e_avg = 0;
//	float spin_avg = 0;
//	float spin_total = 0;
//	float spin_avg2 = 0;
//	float spin_total2 = 0;
//	float phase_total = 0;
//	float phase_avg = 0;
//	float keep_prob = 0;
//	int numb_atoms = size(atom_list);
//	float current_J;
//	float current_K;
//	float new_J;
//	float new_K;
//	float atom_avg_J;
//	float atom_avg_K;
//	float pass_avg_J;
//	float pass_avg_K;
//	int flip_count = 0;
//	int flip_count2 = 0;
//	vector<float> J_K = { 0,0 };
//	std::mt19937_64 rng;
//	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
//	rng.seed(ss);
//	std::uniform_real_distribution<double> unif(0, 1);
//	float e_lattice = evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K); //only for debug
//	cout << evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K);
//	cout << "\n";
//	for (float temp = temp1; temp < temp2; temp += temp_inc) {
//		e_avg = 0;
//		phase_avg = 0;
//		spin_avg = 0;
//		spin_avg2 = 0;
//		pass_avg_J = 0;
//		pass_avg_K = 0;
//		flip_count = 0;
//		flip_count2 = 0;
//		//e_total = evalLattice(temp, atom_list, cluster_rules, spin_rules);
//		for (int i = 0; i < passes; i++) {
//			//e_total = evalLattice(temp, atom_list, cluster_rules, spin_rules);
//			e_total = 0;
//			phase_total = 0;
//			spin_total = 0;
//			spin_total2 = 0;
//			atom_avg_J = 0;
//			atom_avg_K = 0;
//			for (int site = 0; site < atom_list.size(); site++) {
//				//cout << atom_list[site].getSpecies();
//				// Flip Phase and Spin
//				bool keep = false;
//				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				e_site_old = evalSiteEnergy4(temp, site, atom_list, cluster_rules, spin_rules, J_K);/////////////////////////////MADE IT 4!!!!!!!
//				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				current_J = J_K[0];
//				current_K = J_K[1];
//				old_phase = atom_list[site].getPhase();
//				old_spin = atom_list[site].getSpin();
//				both_same = true;
//				while (both_same == true) {
//					phase_rand = unif(rng);
//					spin_rand = unif(rng);
//					if (phase_rand <= 0.3333333333333333) {
//						new_phase = -1;
//					}
//					else if (phase_rand <= 0.6666666666666666) {
//						new_phase = 0;
//					}
//					else {
//						new_phase = 1;
//					}
//					if (spin_rand <= 0.3333333333333333) {
//						new_spin = -1;
//					}
//					else if (spin_rand <= 0.6666666666666666) {
//						new_spin = 0;
//					}
//					else {
//						new_spin = 1;
//					}
//					if ((old_phase == new_phase) and (old_spin == new_spin)) { both_same = true; }
//					else { both_same = false; }
//				}
//				atom_list[site].setSpin(new_spin);
//				atom_list[site].setPhase(new_phase);
//				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				e_site_new = evalSiteEnergy4(temp, site, atom_list, cluster_rules, spin_rules, J_K);//MADE IT 4!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				new_J = J_K[0];
//				new_K = J_K[1];
//				if (e_site_new <= e_site_old) {
//					flip_count2 += 1;
//					keep = true;
//					e_total += e_site_new;
//					atom_avg_J += new_J;
//					atom_avg_K += new_K;
//				}
//				else {
//					keep_rand = unif(rng);
//					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
//					if (keep_rand <= keep_prob) {
//						keep = true;
//						flip_count += 1;
//						e_total += e_site_new;
//						atom_avg_J += new_J;
//						atom_avg_K += new_K;
//					}
//					else {
//						atom_list[site].setPhase(old_phase);
//						atom_list[site].setSpin(old_spin);
//						keep = false;
//						e_total += e_site_old;
//						atom_avg_J += current_J;
//						atom_avg_K += current_K;
//					}
//				}
//				current_phase = atom_list[site].getPhase();
//				phase_total += abs(current_phase);
//				current_spin = atom_list[site].getSpin();
//				if (atom_list[site].getSpecies() != 0) {
//					for (int neighbors = 0; neighbors < 6; neighbors++) {
//						spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
//					}
//				}
//			}
//			phase_avg += phase_total;
//			spin_avg += spin_total;
//			spin_avg2 += spin_total2;
//			e_avg += e_total;
//			pass_avg_J += atom_avg_J;
//			pass_avg_K += atom_avg_K;
//		}
//		cout << temp;
//		cout << " , ";
//		cout << e_avg / passes / numb_atoms * 16;
//		cout << " , ";
//		cout << phase_avg / passes / numb_atoms;
//		cout << " , ";
//		cout << spin_avg / passes / numb_atoms / 6 * 2;
//		cout << " , ";
//		cout << spin_avg2 / passes / numb_atoms;
//		cout << " , ";
//		cout << pass_avg_J / passes / numb_atoms;
//		cout << " , ";
//		cout << pass_avg_K / passes / numb_atoms;
//		cout << " , ";
//		cout << pass_avg_K / pass_avg_J;
//		cout << " , ";
//		cout << flip_count;
//		cout << " , ";
//		cout << flip_count2;
//		cout << "\n";
//	}
//}
//// Mixed Cluster and Wolf Algorithm from Bouadci and Crneiro
//void runMetropolis3(int sub_passes, int total_passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
//	float Kb = .0000861733035;
//	float e_total = 0;
//	float H_site_old = 0;
//	float H_site_new = 0;
//	float spin_rand = 0;
//	float phase_rand = 0;
//	float keep_rand = 0;
//	int old_phase = 0;
//	int new_phase = 0;
//	int old_spin = 0;
//	int new_spin = 0;
//	int current_spin = 0;
//	int current_phase = 0;
//	bool spin_same;
//	float e_avg = 0;
//	float spin_avg = 0;
//	float spin_total = 0;
//	float spin_avg2 = 0;
//	float spin_total2 = 0;
//	float phase_total = 0;
//	float phase_avg = 0;
//	float keep_prob = 0;
//	int numb_atoms = size(atom_list);
//	float atom_avg_J;
//	float atom_avg_K;
//	float pass_avg_J;
//	float pass_avg_K;
//	int flip_count = 0;
//	int flip_count2 = 0;
//	float phase_avg_abs = 0;
//	float phase_total_abs = 0;
//	Cluster cluster;
//	float rand;
//	float H_cluster_old;
//	float H_cluster_new;
//	float prob;
//	int seed_site;
//	int seed_phase;
//	bool phase_same;
//	vector<float> J_K = { 0,0 };
//	std::mt19937_64 rng;
//	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
//	rng.seed(ss);
//	std::uniform_real_distribution<double> unif(0, 1);
//	std::random_device rd;     // only used once to initialise (seed) engine
//	std::mt19937 rng_i(rd());    // random-number engine used (Mersenne-Twister in this case)
//	std::uniform_int_distribution<int> uni(0, atom_list.size() - 1); // guaranteed unbiased
//	cout << evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K);
//	cout << "\n";
//	for (int site = 0; site < atom_list.size(); site++) {
//		init_calcJK(site, atom_list, cluster_rules, spin_rules);
//	}
//	for (float temp = temp1; temp < temp2; temp += temp_inc) {
//		//cout << "\n" << temp << "\n";
//		e_avg = 0;
//		phase_avg = 0;
//		phase_avg_abs = 0;
//		spin_avg = 0;
//		spin_avg2 = 0;
//		flip_count = 0;
//		flip_count2 = 0;
//		//cout << "starting spin flips \n";
//		for (int i = 0; i < total_passes; i++) {
//			for (int j = 0; j < sub_passes; j++) {
//				spin_total = 0;
//				spin_total2 = 0; 
//				for (int site = 0; site < atom_list.size(); site++) {
//					// Flip Spin
//					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//					H_site_old = evalSiteEnergy6(temp, site, atom_list, cluster_rules, spin_rules);                   ////////////////
//					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//					old_spin = atom_list[site].getSpin();
//					spin_same = true;
//					while (spin_same == true) {
//						spin_rand = unif(rng);
//						if (spin_rand <= 0.3333333333333333) {
//							new_spin = -1;
//						}
//						else if (spin_rand <= 0.6666666666666666) {
//							new_spin = 0;
//						}
//						else {
//							new_spin = 1;
//						}
//						if (new_spin != old_spin) { spin_same = false; }
//					}
//					atom_list[site].setSpin(new_spin);
//					re_calcJK(site, old_spin, atom_list, cluster_rules, spin_rules);
//					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//					H_site_new = evalSiteEnergy6(temp, site, atom_list, cluster_rules, spin_rules);                   ////////////////
//					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//					if (H_site_new <= H_site_old) {
//						flip_count2 += 1;
//					}
//					else {
//						keep_rand = unif(rng);
//						keep_prob = exp(-1 / (Kb*temp)*(H_site_new - H_site_old));
//						if (keep_rand < keep_prob) {
//							flip_count += 1;
//						}
//						else {
//							atom_list[site].setSpin(old_spin);
//							re_calcJK(site, new_spin, atom_list, cluster_rules, spin_rules);
//						}
//					}
//					current_spin = atom_list[site].getSpin();
//					spin_total2 += current_spin;
//					if (atom_list[site].getSpecies() != 0) {
//						for (int neighbors = 0; neighbors < 6; neighbors++) {
//							spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
//						}
//					}
//				}
//				spin_avg += spin_total/numb_atoms;
//				spin_avg2 += spin_total2/numb_atoms;
//			}
//		    // Starting Mixed Cluster/Wolff moves
//			seed_site = uni(rng_i);
//			seed_phase = atom_list[seed_site].getPhase();
//			phase_same = true;
//			while (phase_same == true)
//			{
//				phase_rand = unif(rng);
//				if (phase_rand <= 0.3333333333333333) {
//					new_phase = -1;
//				}
//				else if (phase_rand <= 0.6666666666666666) {
//					new_phase = 0;
//				}
//				else {
//					new_phase = 1;
//				}
//				if (new_phase != seed_phase) { phase_same = false; }
//			}
//			//cout << "\entering cluster algorythm \n";
//			if (seed_phase*new_phase == -1) {
//				//cout << "running wolff algorythm \n cluster will be accepted \n";
//				cluster.plant_cluster(seed_site, atom_list);
//				cluster.growClusterWolff(temp, atom_list);
//				//cout << "cluster grown, size = " << cluster.clusterSize() << "\n";
//				flipCluster(seed_phase, new_phase, atom_list, cluster);
//			}
//			else {
//				//cout << "running mixed cluster algorythm \n";
//				cluster.plant_cluster(seed_site, atom_list);
//				cluster.growClusterMixed(temp, new_phase, atom_list);
//				//cout << "cluster grown, size = " << cluster.clusterSize() << "\n";
//				H_cluster_old = evalCluster(atom_list, cluster, cluster_rules, spin_rules, J_K, temp);
//				flipCluster(seed_phase, new_phase, atom_list, cluster);
//				H_cluster_new = evalCluster(atom_list, cluster, cluster_rules, spin_rules, J_K, temp);
//				//cout << "H_cluster_new " << H_cluster_new << " H_cluster_old " << H_cluster_old << "\n";
//				if (H_cluster_new <= H_cluster_old) {
//					//cout << "accepting MC cluster flip: new energy < old energy \n";
//				}
//				else {
//					rand = unif(rng);
//					prob = exp(-1 / (Kb*temp)*(H_cluster_new - H_cluster_old));
//					//cout << "Prob " << prob << "\n";
//					if (rand < prob) {
//						//cout << "accepting MC cluster flip \n";
//					}
//					else {
//						//cout << "rejecting MC cluster flip \n";
//						flipCluster(new_phase, seed_phase, atom_list, cluster, true);
//					}
//				}
//			}
//			e_avg += evalLattice(temp, atom_list, cluster_rules, spin_rules, J_K);
//			for (int k = 0; k < numb_atoms; k++) { phase_avg += abs(atom_list[k].getPhase()); }
//		}
//		phase_total = phase_avg/total_passes;
//		e_total = e_avg/total_passes;//evalLattice(temp, atom_list, cluster_rules, spin_rules, J_K);
//		cout << temp;
//		cout << " , ";
//		cout << e_total; // e_avg / passes / numb_atoms * 16;
//		cout << " , ";
//		cout << phase_total / numb_atoms;// / passes / numb_atoms;
//		cout << " , ";
//		cout << spin_avg / (sub_passes*total_passes) / 6 * 2;
//		cout << " , ";
//		cout << spin_avg2 / (sub_passes*total_passes);
//		cout << " , ";
//		cout << flip_count;
//		cout << " , ";
//		cout << flip_count2;
//		cout << "\n";
//	}
//}