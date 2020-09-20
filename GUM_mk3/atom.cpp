//#include "atom.h"
//
//Atom::Atom(void) {
//	species = 10;
//	spin = 10;
//	phase = 10;
//	index = 10;
//	cluster_status = "unknown";
//}
//Atom::Atom(int _index, int _species, int _spin, int _phase, float _pos[3]) {
//	index = _index;
//	species = _species;
//	spin = _spin;
//	phase = _phase;
//	pos[0] = _pos[0];
//	pos[1] = _pos[1];
//	pos[2] = _pos[2];
//	cluster_status = "unknown";
//}
//void Atom::setNeighbor(int _order, string _plain, int _index) {
//	if (_order == 1) {
//		neighbors_1.push_back(_index);
//		neighbor_plain_1.push_back(_plain);
//	}
//	if (_order == 2) {
//		neighbors_2.push_back(_index);
//		neighbor_plain_2.push_back(_plain);
//	}
//	if (_order == 3) {
//		neighbors_3.push_back(_index);
//		neighbor_plain_3.push_back(_plain);
//	}
//	neighbors.push_back(_index);
//	neighbor_plains.push_back(_plain);
//	neighbor_orders.push_back(_order);
//}
//int Atom::getNeighbor(int _order, int _neighbor, vector<Atom> &atom_list) {
//	int neighbor_index;
//	if (_order == 1) {
//		neighbor_index = neighbors_1[_neighbor];
//	}
//	else if (_order == 2) {
//		neighbor_index = neighbors_2[_neighbor];
//	}
//	else if (_order == 3) {
//		neighbor_index = neighbors_3[_neighbor];
//	}
//	else {
//		cout << "error";
//		cout << '\n';
//		neighbor_index= 1000000000;
//	}
//	return neighbor_index;
//}
//int Atom::getNeighborSpin(int _neighbor, vector<Atom> &atom_list) {
//	int neighbor_index = neighbors[_neighbor];
//	int neighbor_spin = atom_list[neighbor_index].getSpin();
//	return neighbor_spin;
//}
//int Atom::getNeighborSpecies(int _neighbor, vector<Atom> &atom_list) {
//	int neighbor_index = neighbors[_neighbor];
//	int neighbor_species = atom_list[neighbor_index].getSpecies();
//	return neighbor_species;
//}
//int Atom::getNeighborPhase(int _neighbor, vector<Atom> &atom_list) {
//	int neighbor_index = neighbors[_neighbor];
//	int neighbor_phase = atom_list[neighbor_index].getPhase();
//	return neighbor_phase;
//}
//void Atom::setSpin(int _spin) {
//	spin = _spin;
//}
//void Atom::setSpecies(int _species) {
//	species = _species;
//}
//void Atom::setPhase(int _phase) {
//	phase = _phase;
//}
//int Atom::getSpin(void) {
//	return spin;
//}
//int Atom::getSpecies(void) {
//	return species;
//}
//int Atom::getPhase(void) {
//	return phase;
//}