#pragma once
#ifndef rule_h
#define rule_h
#include "file_io.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;


class Rule {
private: 
	float energy_contribution;
	int rule_type;
	string phase;
	vector<int> species;
	vector<float> distances;
	vector<int> spins;
public:

	Rule(void);
	Rule(float _energy_contribution, int _rule_type, string _phase, vector<int> _species, vector<float> _distances, vector<int> _spins);
	Rule(string  inputline);
	string GetPhase();
	float GetEnrgCont();
	int GetType();
	int GetLength();
	vector<int> GetSpecies();
	vector<float> GetDists();
	vector<int> GetSpins();

};

#endif