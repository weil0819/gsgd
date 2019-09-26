/*
 * Graph.h
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
	int dist;	// threshold of distance

	ui *pstart; //offset of neighbors of nodes
	int *edges; //adjacent ids of edges 
	ui *reverse; //the position of reverse edge in edges
	int *min_cn; //minimum common neighbor: -2 means not similar; -1 means similar; 0 means not sure; > 0 means the minimum common neighbor

	int *pa;
	int *rank; //pa and rank are used for the disjoint-set data structure

	int *cid; //cluster id 

	int *degree;
	int *similar_degree; //number of adjacent edges with similarity no less than epsilon

	float *ploc;		// location of nodes

	vector<pair<int,int> > noncore_cluster;	// for non-core nodes, <cluster ID, node ID>
	vector<int> invalid_cluster;			// invalid cluster IDs

public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph() ;

	void gsgd(const char *eps_s, int miu, int dist) ;

	void dgcd(const char *eps_s, int miu, int dist) ;

	void cluster_noncore_vertices(int eps_a2, int eps_b2, int miu) ;

	void output(const char *eps_s, const char *miu, const char *dist) ; 

	void cluster_count(const char *eps_s, const char *miu, const char *dist) ;

	

private:
	int naive_similar_check(int u, int v, int eps_a2, int eps_b2) ;

	double euclidean_dist2(int u, int v) ;
	double euclidean_dist2(float ux, float uy, float vx, float vy) ;

	double mcc_radius2(int v1, int v2, int v3) ;

	int find_root(int u) ;
	void my_union(int u,int v) ;

	void get_eps(const char *eps_s, int &eps_a, int &eps_b) ; 
	void prune_and_cross_link(int dist) ;
	ui binary_search(const int *array, ui b, ui e, int val) ;

	void eps_neighbor(ui v, int eps_a2, int eps_b2, int dist, vector<pair<ui, ui> > &H, vector<ui> &S) ;

};

#endif
