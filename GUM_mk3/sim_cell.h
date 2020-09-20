#pragma once
#pragma once
#ifndef sim_cell_h
#define sim_cell_h
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "file_io.h"
using namespace std;


class SimCell {
public:
	float cutoff;
	string sim_type;
	vector<int> species_types;
	vector<int> species_numbs;
	//vector<int> clust_count{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	int sup_cell[3];
	float cell_dim[3];
	int numb_atoms;
	int numb_cells[3];
	float unit_LC[3];
	int X_num;
	class Atom {
	private:
		int species;
		int spin;
		int phase;
		string cluster_status;
	public:
		vector<float> neighbor_dists;
		vector<int> neighbors;
		int index;
		float pos[3];
		Atom(void);
		Atom(int _index, int _species, int _spin, int _phase, vector<float> _pos);
		void setSpin(int _spin);
		void setSpecies(int _species);
		void setPhase(int _phase);
		int getSpin(void);
		int getSpecies(void);
		int getPhase(void);
		int getNeighborSpin(int _neighbor, SimCell &sim_cell);
		int getNeighborSpecies(int _neighbor, SimCell &sim_cell);
		int getNeighborPhase(int _neighbor, SimCell &sim_cell);
		int getNeighborIndex(int _neighbor, SimCell &sim_cell);
		float getNeighborDist(int _neighbor, SimCell &sim_cell);
		int getNumbNeighbors(int _site, SimCell &sim_cell);
	};
	vector<Atom> atom_list;
	vector<Atom> unit_cell;

	SimCell(void);
	SimCell(string POSCAR_file, int _sup_cell[3], vector<int> &_species_numbs, float _cutoff, string _sim_type, string phase_init, string spin_init, string species_init);
	void initSimCell(string POSCAR_file, vector<float> &dist_list, int _sup_cell[3], vector<int> &_species_numbs, float _cutoff, string _sim_type, string phase_init, string spin_init, string species_init, bool use_poscar);
	void fillUnitCell(string POSCAR_file, bool use_poscar);
	void fillAtomList(vector<vector<float>> &_pos_list, vector<int> &_species_list, vector<float> dist_list, string phase_init, string spin_init, string species_init);
	void make_supercell(vector<vector<float>> &_pos_list, vector<int> &_species_list, string phase_init);
	void setNeighborDists(vector<float> &dist_list);
	float findAtomDists(int atom1, int atom2);
	};
#endif // !sim_cell_h
