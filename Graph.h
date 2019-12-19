/*
 * Graph.h
 *
 *  Created on: 02 Sept, 2019
 *      Author: Wesley
 */

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "HashSet.h"

using namespace std;


struct Circle {
	static const Circle INVALID;

	double centerX;
	double centerY;
	double radius;	
};


class Graph {
private:
	string dir; 	//input graph directory
	ui n, m; 		//number of nodes and edges of the graph

	int eps_a2, eps_b2, miu; // eps_a2/eps_b2 = eps^2 
	int gamma;				 // threshold of distance

	ui *pstart; 	//offset of neighbors of nodes
	int *edges; 	//adjacent ids of edges 
	ui *reverse; 	//the position of reverse edge in edges
	int *min_cn; 	//minimum common neighbor: -2 means not similar; -1 means similar; 0 means not sure; > 0 means the minimum common neighbor

	int *pa;
	int *rank; 		//pa and rank are used for the disjoint-set data structure

	int *cid; 		//cluster id 

	int *degree;
	int *similar_degree; 	//number of adjacent edges with similarity no less than epsilon

	float *ploc;			// location of nodes

	vector<pair<int,int> > noncore_cluster;	// for non-core nodes, <cluster ID, node ID>
	vector<int> invalid_cluster;			// invalid cluster IDs


public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph() ;

	// AI3D
	void gsgd(const char *eps_s, int mu, int gamma) ;

	void dgcd(const char *eps_s, int mu, int gamma) ;


	// cluster + naive MCC  
	void baseline(const char *eps_s, int mu, int gamma) ; 

	// cluster + random MCC 
	void advance(const char *eps_s, int mu, int gamma) ;

	// greedy top-k 
	void greedy_topk(const char *eps_s, int mu, int gamma, int k) ; 

	// swap top-k 
	void swap_topk(const char *eps_s, int mu, int gamma, int k) ; 

	// naive top-k
	void topk(const char *eps_s, int mu, int gamma, int k) ;


	void cluster_noncore_vertices(int eps_a2, int eps_b2, int miu) ;
	void cluster_noncore_vertices(int eps_a2, int eps_b2, int miu, unordered_set<int> &US) ;

	void output(const char *eps_s, const char *miu, const char *gamma) ; 

	void cluster_count(const char *eps_s, const char *miu, const char *gamma) ;  

	void renew_cluster(int eps_a2, int eps_b2, int mu, vector<int> &curList, vector<vector<int> > &cluster) ; 

	void floyd_diameter(int eps_a2, int eps_b2, int mu, int gamma, vector<int> &cluster, vector<int> &output) ;

	
private:
	int naive_similar_check(int u, int v, int eps_a2, int eps_b2) ;
	int naive_similar_check(int u, int v, int eps_a2, int eps_b2, unordered_set<int> &US) ;

	double euclidean_dist2(int u, int v) ;
	double euclidean_dist2(float ux, float uy, float vx, float vy) ;

	double mcc_radius2(int v1, int v2, int v3) ;
	double mcc_radius2(int v1, int v2, int v3, double &centerX, double &centerY) ;

	int find_root(int u) ;
	int find_root(int *pa, int u) ;

	void my_union(int u,int v) ; 
	void my_union(int *pa, int *rank, int u, int v) ;
	

	void get_eps(const char *eps_s, int &eps_a, int &eps_b) ; 

	void prune_and_cross_link(int dist) ;

	ui binary_search(const int *array, ui b, ui e, int val) ;
	ui binary_search(const int *array, ui e, int val) ;

	void eps_neighbor(ui v, int eps_a2, int eps_b2, int dist, vector<pair<ui, ui> > &H, vector<ui> &S) ;

	bool is_subset(vector<int> &A, vector<int> &B) ; 

	// Circle Functions ... 
	Circle make_mcc(vector<int> &cluster) ;

	Circle make_circle_one_point(vector<int> &cluster, int end, ui u) ;

	Circle make_circle_two_point(vector<int> &cluster, int end, ui u, ui v) ;

	Circle make_diameter(ui u, ui v) ; 

	Circle make_circumcircle(ui a, ui b, ui c) ;

	double make_cross(double px, double py, double qx, double qy) ; 

	double make_distance(double px, double py, double qx, double qy) ;

	bool make_contains(double centerX, double centerY, double radius, ui p) ;

};

#endif
