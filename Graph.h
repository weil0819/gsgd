/*
 * Preprocess.cpp
 *
 *  Created on: 02 Sept, 2019
 *      Author: Wesley
 */

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"

using namespace std;

class Graph {
private:
	string dir; //input graph directory
	ui n, m; //number of nodes and edges of the graph

	int eps_a2, eps_b2, miu; // eps_a2/eps_b2 = eps^2

	ui *pstart; //offset of neighbors of nodes
	int *edges; //adjacent ids of edges
	int *min_cn; //minimum common neighbor: -2 means not similar; -1 means similar; 0 means not sure; > 0 means the minimum common neighbor

	int *pa;
	int *rank; //pa and rank are used for the disjoint-set data structure

	int *cid; //cluster id 

	int *degree;
	int *similar_degree; //number of adjacent edges with similarity no less than epsilon

	double *ploc;	// location of nodes

	vector<pair<int,int> > noncore_cluster;

public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph() ;

	void baseline(const char *eps_s, int miu) ;

	void cluster_noncore_vertices(int eps_a2, int eps_b2, int mu) ;

	void output(const char *eps_s, const char *miu) ; 

	void cluster_count(const char *eps_s, const char *miu) ;

private:
	int naive_similar_check(int u, int v, int eps_a2, int eps_b2) ;

	int find_root(int u) ;
	void my_union(int u,int v) ;

	void get_eps(const char *eps_s, int &eps_a, int &eps_b) ; 


};

#endif
