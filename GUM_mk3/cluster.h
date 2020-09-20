#pragma once
#ifndef cluster_h
#define cluster_h
#include "atom.h"
#include <vector>
#include <chrono>
#include <string>
#include <iostream>
#include <cmath>
#include <random>
using namespace std;

class Cluster {
private:
	vector<int> branch;
	vector<int> root;
	int seed;
	int new_phase;
public:
	void setNewPhase(int phase);
	int getNewPhase();
	//void plant_cluster(int _seed, vector<Atom> &atom_list);
	vector<int> cluster_list;
	bool continue_growth;
	Cluster(void);
	Cluster(int _next_seed, vector<int> &_cluster_list);
//	bool inList(int site);
//	int clusterSize();
//	void growClusterWolff(float temp, vector<Atom> &atom_list);
//	void growClusterMixed(float temp, int _new_phase, vector<Atom> &atom_list);
//	void setSiteStatus(int site, string status, vector<Atom> &atom_list);
//	string getSiteStatus(int site, vector<Atom> &atom_list);
};
#endif // !cluster_h
