#include "cluster.h"

// This class is used only for BEG mixed cluster MC 
Cluster::Cluster(void) {
	continue_growth = true;
	seed = 0;
	cluster_list = {};
}

Cluster::Cluster(int _seed, vector<int> &_cluster_list) {
	seed = _seed;
	continue_growth = true;
	cluster_list = _cluster_list;
}

void Cluster::setNewPhase(int phase) {
	new_phase = phase;
}

int Cluster::getNewPhase() {
	return new_phase;
}

//void Cluster::plant_cluster(int _seed, vector<Atom> &atom_list) {
//	for (int i = 0; i < atom_list.size(); i++) {
//		setSiteStatus(i, "unknown", atom_list);
//	}
//	seed = _seed;
//	setSiteStatus(seed, "seed", atom_list);
//	branch = {};
//	root = {};
//	cluster_list = {};
//	cluster_list.push_back(seed);
//	root.push_back(seed);
//	continue_growth = true;
//}
//
//void Cluster::setSiteStatus(int site, string status, vector<Atom> &atom_list) {
//	atom_list[site].setClusterStatus(status);
//}
//
//string Cluster::getSiteStatus(int site, vector<Atom> &atom_list) {
//	return atom_list[site].getClusterStatus();
//}
//
//bool Cluster::inList(int site) {
//	bool status;
//	if (std::find(cluster_list.begin(), cluster_list.end(), site) != cluster_list.end()) {
//		status = true;
//	}
//	else {
//		status = false;
//	}
//	return status;
//}
//
//int Cluster::clusterSize() {
//	return cluster_list.size();
//}
//
//void Cluster::growClusterWolff(float temp, vector<Atom> &atom_list) {
//	float Kb = .000086173324;
//	float B = 1 / (Kb*temp);
//	float rand;
//	float prob;
//	int neighbor_site;
//	int neighbor = 0;
//	int numb_neighbors = atom_list[root.back()].getNumbNeighbors();
//	float root_J;
//	float link_J;
//	std::mt19937_64 rng;
//	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
//	rng.seed(ss);
//	std::uniform_real_distribution<double> unif(0, 1);
//	while (continue_growth == true) {
//		root_J = atom_list[root.back()].J;
//		neighbor_site = atom_list[root.back()].getNeighbor(1, neighbor, atom_list);
//		if (getSiteStatus(neighbor_site, atom_list) == "unknown") {
//			if (atom_list[neighbor].getPhase() != atom_list[seed].getPhase()) {
//				setSiteStatus(neighbor_site, "cap", atom_list);
//				neighbor += 1;
//				if (neighbor >= 8) {
//					root.pop_back();
//					neighbor = 0;
//					if (root.size() == 0) {
//						continue_growth = false;
//						root.push_back(seed);
//					}
//				}
//			}
//			else {
//				rand = unif(rng);
//				link_J = -(root_J + atom_list[neighbor_site].J) / (2);
//				prob = 1 - exp(-2 * B * link_J);
//				if (rand <= prob) {
//					setSiteStatus(neighbor_site, "root", atom_list);
//					root.push_back(neighbor_site);
//					cluster_list.push_back(neighbor_site);
//					neighbor = 0;
//					if (cluster_list.size() >= atom_list.size()) { continue_growth = false; }
//				}
//				else {
//					setSiteStatus(neighbor_site, "cap", atom_list);
//					neighbor += 1;
//					if (neighbor >= 8) {
//						root.pop_back();
//						neighbor = 0;
//						if (root.size() == 0) { 
//							continue_growth = false;
//							root.push_back(seed);
//						}
//					}
//				}
//			}
//		}
//		else {
//			neighbor++;
//			if (neighbor >= 8) {
//				neighbor = 0;
//				setSiteStatus(neighbor_site, "cap", atom_list);
//				root.pop_back();
//				if (root.size() == 0) {
//					continue_growth = false;
//					root.push_back(seed);
//				}
//			}
//		}
//	}
//}
//
//void Cluster::growClusterMixed(float temp, int _new_phase, vector<Atom> &atom_list) {
//	float Kb = .000086173324;
//	float B = 1 / (Kb*temp);
//	float rand;
//	float prob;
//	int neighbor_site;
//	int neighbor = 0;
//	int numb_neighbors = atom_list[root.back()].getNumbNeighbors();
//	float root_J;
//	float root_K;
//	float link_J;
//	float link_K;
//	float seed_phase = atom_list[seed].getPhase();
//	setNewPhase(_new_phase);
//	std::mt19937_64 rng;
//	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
//	rng.seed(ss);
//	std::uniform_real_distribution<double> unif(0, 1);
//	if ((seed_phase == 1 and getNewPhase() == 0) or (seed_phase == 0 and getNewPhase() == -1)) {
//		while (continue_growth == true) {
//			root_J = atom_list[root.back()].J;
//			root_K = atom_list[root.back()].K;
//			neighbor_site = atom_list[root.back()].getNeighbor(1, neighbor, atom_list);
//			if (getSiteStatus(neighbor_site, atom_list) == "unknown") {
//				if (atom_list[neighbor_site].getPhase() == 1 or atom_list[neighbor_site].getPhase() == 0) {
//					if (atom_list[neighbor_site].getPhase() == atom_list[root.back()].getPhase()) {
//						rand = unif(rng);
//						link_J = -(root_J + atom_list[neighbor_site].J) / (2);
//						link_K = -(root_K + atom_list[neighbor_site].K) / (2);
//						prob = 1 - exp(B * (-link_J - link_K / 3));
//						if (rand <= prob) {
//							setSiteStatus(neighbor_site, "root", atom_list);
//							root.push_back(neighbor_site);
//							cluster_list.push_back(neighbor_site);
//							neighbor = 0;
//							if (cluster_list.size() >= atom_list.size()) { continue_growth = false; }
//						}
//						else {
//							setSiteStatus(neighbor_site, "cap", atom_list);
//							neighbor += 1;
//							if (neighbor >= 8) {
//								root.pop_back();
//								neighbor = 0;
//								if (root.size() == 0) {
//									continue_growth = false;
//									root.push_back(seed);
//								}
//							}
//						}
//					}
//					else {
//						rand = unif(rng);
//						link_J = -(root_J + atom_list[neighbor_site].J) / 2;
//						link_K = -(root_K + atom_list[neighbor_site].K) / 2;
//						prob = 1 - exp(B * (-link_J + link_K / 3));
//						if (rand <= prob) {
//							setSiteStatus(neighbor_site, "root", atom_list);
//							root.push_back(neighbor_site);
//							cluster_list.push_back(neighbor_site);
//							neighbor = 0;
//							if (cluster_list.size() >= atom_list.size()) { continue_growth = false; }
//						}
//						else {
//							setSiteStatus(neighbor_site, "cap", atom_list);
//							neighbor += 1;
//							if (neighbor >= 8) {
//								root.pop_back();
//								neighbor = 0;
//								if (root.size() == 0) {
//									continue_growth = false;
//									root.push_back(seed);
//								}
//							}
//						}
//					}
//				}
//				else {
//					setSiteStatus(neighbor_site, "cap", atom_list);
//					neighbor += 1;
//					if (neighbor >= 8) {
//						root.pop_back();
//						neighbor = 0;
//						if (root.size() == 0) {
//							continue_growth = false;
//							root.push_back(seed);
//						}
//					}
//				}
//			}
//			else { 
//				neighbor++;
//				if (neighbor >= 8) {
//					neighbor = 0;
//					setSiteStatus(neighbor_site, "cap", atom_list);
//					root.pop_back();
//					if (root.size() == 0) {
//						continue_growth = false;
//						root.push_back(seed);
//					}
//				}
//			}
//		}
//	}
//	else if ((seed_phase == -1 and getNewPhase() == 0) or (seed_phase == 0 and getNewPhase() == 1)) {
//		while (continue_growth == true) {
//			root_J = atom_list[root.back()].J;
//			root_K = atom_list[root.back()].K;
//			neighbor_site = atom_list[root.back()].getNeighbor(1, neighbor, atom_list);
//			if (getSiteStatus(neighbor_site, atom_list) == "unknown") {
//				if (atom_list[neighbor_site].getPhase() == -1 or atom_list[neighbor_site].getPhase() == 0) {
//					if (atom_list[neighbor_site].getPhase() == atom_list[root.back()].getPhase()) {
//						rand = unif(rng);
//						link_J = -(root_J + atom_list[neighbor_site].J) / (2);
//						link_K = -(root_K + atom_list[neighbor_site].K) / (2);
//						prob = 1 - exp(B * (-link_J - link_K / 3));
//						if (rand <= prob) {
//							setSiteStatus(neighbor_site, "root", atom_list);
//							root.push_back(neighbor_site);
//							cluster_list.push_back(neighbor_site);
//							neighbor = 0;
//							if (cluster_list.size() >= atom_list.size()) { continue_growth = false; }
//						}
//						else {
//							setSiteStatus(neighbor_site, "cap", atom_list);
//							neighbor += 1;
//							if (neighbor >= 8) {
//								root.pop_back();
//								neighbor = 0;
//								if (root.size() == 0) {
//									continue_growth = false;
//									root.push_back(seed);
//								}
//							}
//						}
//					}
//					else {
//						rand = unif(rng);
//						link_J = -(root_J + atom_list[neighbor_site].J) / 2;
//						link_K = -(root_K + atom_list[neighbor_site].K) / 2;
//						prob = 1 - exp(B * (-link_J + link_K / 3));
//						if (rand <= prob) {
//							setSiteStatus(neighbor_site, "root", atom_list);
//							root.push_back(neighbor_site);
//							cluster_list.push_back(neighbor_site);
//							neighbor = 0;
//							if (cluster_list.size() >= atom_list.size()) { continue_growth = false; }
//						}
//						else {
//							setSiteStatus(neighbor_site, "cap", atom_list);
//							neighbor += 1;
//							if (neighbor >= 8) {
//								root.pop_back();
//								neighbor = 0;
//								if (root.size() == 0) {
//									continue_growth = false;
//									root.push_back(seed);
//								}
//							}
//						}
//					}
//				}
//				else {
//					setSiteStatus(neighbor_site, "cap", atom_list);
//					neighbor += 1;
//					if (neighbor >= 8) {
//						root.pop_back();
//						neighbor = 0;
//						if (root.size() == 0) {
//							continue_growth = false;
//							root.push_back(seed);
//						}
//					}
//				}
//			}
//			else {
//				neighbor++;
//				if (neighbor >= 8) {
//					neighbor = 0;
//					setSiteStatus(neighbor_site, "cap", atom_list);
//					root.pop_back();
//					if (root.size() == 0) {
//						continue_growth = false;
//						root.push_back(seed);
//					}
//				}
//			}
//		}
//	}
//}
