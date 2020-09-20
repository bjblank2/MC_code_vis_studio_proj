#include "rule.h"
 

// This is the rule class (mc_rules object). It servs as a container for the MC rules and is used to create the rule map
Rule::Rule(void) {
}    

Rule::Rule(float _energy_contribution, int _rule_type, string _phase, vector<int> _species, vector<float> _distances, vector<int> _spins) {
	energy_contribution = _energy_contribution;
	rule_type = _rule_type;
	phase = _phase;
	species = _species;
	distances = _distances;
	spins = _spins;
}

Rule::Rule(string  inputline) {
	inputline.erase(std::remove(inputline.begin(), inputline.end(), ']'), inputline.end());
	inputline.erase(std::remove(inputline.begin(), inputline.end(), '['), inputline.end());
	vector<string> split_string = split(inputline, ":");
	if (split_string.size() != 4) cout << "error";
	vector<string>s_species = split(split_string[0], ",");
	for (int i = 0; i < s_species.size(); i++) {
		species.push_back(stoi(s_species[i]));
	}
	vector<string>s_dists = split(split_string[1], ",");
	for (int i = 0; i < s_dists.size(); i++) {
		distances.push_back(stof(s_dists[i]));
	}
	rule_type = stoi(split_string[2]);
	energy_contribution = stof(split_string[3]);
}

string Rule::GetPhase() {
	return phase;
}
float Rule::GetEnrgCont() {
	return energy_contribution;
}
int Rule::GetType() {
	return rule_type;
}
int Rule::GetLength() {
	return species.size();
}
vector<int> Rule::GetSpecies() {
	return species;
}
vector<float> Rule::GetDists() {
	return distances;
}
vector<int> Rule::GetSpins() {
	return spins;
}