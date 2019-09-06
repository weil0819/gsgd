/*
 * Graph.cpp
 *
 *  Created on: 02 Sept, 2019
 *      Author: Wesley
 */

#include "Utility.h"
#include "Graph.h"

Graph::Graph(const char *_dir) {
	dir = string(_dir);

	n = m = 0;

	eps_a2 = eps_b2 = miu = 0;
	dist = 0;

	pstart = NULL;
	edges = NULL;
	min_cn = NULL;

	cid = NULL;

	degree = NULL;
	similar_degree = NULL;

	pa = NULL;
	rank = NULL;

	ploc = NULL;
}

Graph::~Graph() {
	if(pstart != NULL) {
		delete[] pstart;
		pstart = NULL;
	}
	if(edges != NULL) {
		delete[] edges;
		edges = NULL;
	}
	if(min_cn != NULL) {
		delete[] min_cn;
		min_cn = NULL;
	}
	if(cid != NULL) {
		delete[] cid;
		cid = NULL;
	}
	if(degree != NULL) {
		delete[] degree;
		degree = NULL;
	}
	if(similar_degree != NULL) {
		delete[] similar_degree;
		similar_degree = NULL;
	}
	if(pa != NULL) {
		delete[] pa;
		pa = NULL;
	}
	if(rank != NULL) {
		delete[] rank;
		rank = NULL;
	}
	if(ploc != NULL) {
		delete[] ploc;
		ploc = NULL;
	}
}

void Graph::read_graph() {
	// Read-1: read degree file, including n, m.
	FILE *f = open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != (int)sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	printf("\tn = %u; m = %u\n", n, m/2);

	degree = new int[n];
	fread(degree, sizeof(int), n, f);

	fclose(f);

	// Read-2: read adjacent file.
	f = open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart == NULL) pstart = new ui[n+1];
	if(edges == NULL) edges = new int[m];
	if(min_cn == NULL) min_cn = new int[m];
	memset(min_cn, 0, sizeof(int)*m);

	int *buf = new int[n];

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) fread(buf, sizeof(int), degree[i], f);

		for(ui j = 0;j < degree[i];j ++) edges[pstart[i] + j] = buf[j];

		pstart[i+1] = pstart[i] + degree[i];

		++ degree[i];		// closed degree = open degree + 1
	}

	delete[] buf;

	fclose(f);

	// Read-3: read location file.
	f = open_file((dir + string("/loc.txt")).c_str(), "r");

	if(ploc == NULL) ploc = new float[2*n];
	ui ind = 0;
	while(!feof(f)){
		float a, b;
		fscanf(f, "%f %f", &a, &b);
		ploc[ind++] = a;
		ploc[ind++] = b;
	}

	fclose(f);
}


void Graph::cluster_noncore_vertices(int eps_a2, int eps_b2, int miu) {
	if(cid == NULL) cid = new int[n];		// assign each node a cluster ID
	for(ui u = 0; u < n; u++) cid[u] = n;	// initialize cluster ID

	for(ui u = 0; u < n; u++) {
		if(similar_degree[u] >= miu) {		// u is a core node
			int x = find_root(u);
			if(u < cid[x]) cid[x] = u;		// update cluster ID, which should be minimal ID in cluster
		}
	}

	noncore_cluster.clear();
	noncore_cluster.reserve(n);

	for(ui i = 0;i < n;i ++) {
		if(similar_degree[i] >= miu) {		// i should be core node
			for(ui j = pstart[i];j < pstart[i+1];j ++) {
				// when edges[j] is also a core node or i and edges[j] are in same cluster
				if(similar_degree[edges[j]] >= miu||pa[i] == pa[edges[j]]) continue;

				if(min_cn[j] >= 0) {
					min_cn[j] = naive_similar_check(i, edges[j], eps_a2, eps_b2);
				}

				if(min_cn[j] == -1) noncore_cluster.pb(mp(cid[pa[i]], edges[j]));
			}
		}
	}
}


void Graph::output(const char *eps_s, const char *miu, const char *dist) {
	printf("\t*** Start write result into disk!\n");

	string out_name = dir + "/result-" + string(eps_s) + "-" + string(miu) + "-" + string(dist) + ".txt";
	FILE *fout = open_file(out_name.c_str(), "w");

	fprintf(fout, "c/n vertex_id cluster_id\n");

	int mu = atoi(miu);
	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		fprintf(fout, "c %d %d\n", i, cid[pa[i]]);
	}

	sort(noncore_cluster.begin(), noncore_cluster.end());
	noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
	for(ui i = 0;i < noncore_cluster.size();i ++) {
		fprintf(fout, "n %d %d\n", noncore_cluster[i].second, noncore_cluster[i].first);
	}

	fclose(fout);
}


void Graph::cluster_count(const char *eps_s, const char *miu, const char *dist) {
	printf("\t*** Start write clusters into disk!\n");

	string out_name = dir + "/cluster-" + string(eps_s) + "-" + string(miu) + "-" + string(dist) + ".txt";
	FILE *fout = open_file(out_name.c_str(), "w");

	// fprintf(fout, "cluster_id \t #vertex \t vertex_id \n");

	vector<vector<int> > cluster_set(n);		// cluster_set[i] is the i-th cluster

	// Consider core nodes first.
	int mu = atoi(miu);
	for(ui u = 0; u < n; u++) {
		if(similar_degree[u] >= mu) {
			int cluster_id = cid[pa[u]];
			cluster_set[cluster_id].pb(u);
		}
	}

	// Consider non-core nodes.
	sort(noncore_cluster.begin(), noncore_cluster.end());
	noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end()); 

	for(ui i = 0;i < noncore_cluster.size();i ++) {
		cluster_set[noncore_cluster[i].first].pb(noncore_cluster[i].second);
	}

	// Iterate each cluster and write them into disk.
	for(int i = 0; i < n; i++) {
		if(cluster_set[i].size() > 0) {
			fprintf(fout, "cluster-%d, \t #cluster=%d \n", i, cluster_set[i].size());
			for(int j = 0; j < cluster_set[i].size(); j++) {
				fprintf(fout, "%d\t", cluster_set[i][j]);
			}
			fprintf(fout, "\n************************************************************************\n");
		}
	}

	fclose(fout);

}

void Graph::baseline(const char *eps_s, int miu, int dist) {
	// Pruning-1: prune all neighbors whose distance is larger than dist.
	for(ui u = 0; u < n; u++) {
		for(ui i = pstart[u]; i < pstart[u+1]; i++) {
			ui v = edges[i];
			float ux = ploc[2*u], uy = ploc[2*u+1];
			float vx = ploc[2*v], vy = ploc[2*v+1];
			double dis = (ux-vx)*(ux-vx) + (uy-vy)*(uy-vy);
			if(dis > (double)dist*dist) {
				// Mark v as empty and degree minus 1.
				edges[i] = u;	// mark u's neighbor same as u
				--degree[u];

				// cross-link to find (v,u) ???
			}
		}
	}

	// Step-1: Compute eps value (including eps_a2 and eps_b2).
	int eps_a2 = 0, eps_b2 = 0;
	get_eps(eps_s, eps_a2, eps_b2);
	eps_a2 *= eps_a2;
	eps_b2 *= eps_b2;

	// Step-2: Compute each node's similar_degree value.
	if(similar_degree == NULL) similar_degree = new int[n];
	memset(similar_degree, 0, sizeof(int)*n);

	for(ui u = 0; u < n; u++) {			// scan u 
		if(degree[u] < miu) {			// degree[u] is an upper bound of the size of its eps-neighorhood
			similar_degree[u] = -1;		 
			continue ;
		}
		for(ui i = pstart[u]; i < pstart[u+1]; i++) {
			ui v = edges[i];			// v is u's neighbor 
			if(min_cn[i] < 0) continue;	// min_cn[i] means the i-th position of edges 
			if(u == v) {				// for far neighbor, omit computing the similarity 
				min_cn[i] = -2;
				continue ;
			} 
			min_cn[i] = naive_similar_check(u, v, eps_a2, eps_b2);
			if(min_cn[i] == -1) ++similar_degree[u];	// compute N_eps[u]
		}
	}

	// Step-3: Initialize pa and rank for the disjoint-set data structure.
	if(pa == NULL) pa = new int[n];
	if(rank == NULL) rank = new int[n];

	for(ui i = 0; i < n; i++) {
		pa[i] = i;
		rank[i] = 0;
	}

	// Step-4: Cluster all core nodes (including merge clusters).
	for(ui u = 0; u < n; u++) {
		if(similar_degree[u] < miu) continue;	// u is not a core node
		for(ui i = pstart[u]; i < pstart[u+1]; i++) {
			ui v = edges[i]; 			
			if(u == v || similar_degree[v] < miu) continue;	// // v is marked node or v is not a core node 
			// If v is also a core node and sim(u,v)>=eps, merge them into one cluster.
			if(min_cn[i] == -1) my_union(u,v);
		}
	}

	// Step-5: Cluster non-core nodes.
	cluster_noncore_vertices(eps_a2, eps_b2, miu);


	// Step-6: Return qualified subgraph.



}


int Graph::naive_similar_check(int u, int v, int eps_a2, int eps_b2) {
	int du = degree[u], dv = degree[v];		// du and dv are 

	ui i = pstart[u], j = pstart[v];
	int cn = 2;	// number of common neighbor between u and v
	while(i < pstart[u+1]&&j < pstart[v+1]) { 
		if(edges[i] == u) {
			++ i;
			continue ;
		}
		if(edges[j] == v) {
			++ j;
			continue ;
		}

		if(edges[i] < edges[j]) ++ i;
		else if(edges[i] > edges[j]) ++ j;
		else {
			++ cn;
			++ i;
			++ j;
		}
	}

	// structural similarity(u,v) >= epsilon => cn^2 >= eps^2*(du*dv)
	if(((long long)cn)*((long long)cn)*eps_b2 >= ((long long)du)*((long long)dv)*eps_a2) return -1;

	return -2;
}


int Graph::find_root(int u) {
	int x = u;
	while(pa[x] != x) x = pa[x];

	while(pa[u] != x) {
		int tmp = pa[u];
		pa[u] = x;
		u = tmp;
	}

	return x;
}


void Graph::my_union(int u, int v) {
	int ru = find_root(u);
	int rv = find_root(v);

	if(ru == rv) return ;

	if(rank[ru] < rank[rv]) pa[ru] = rv;
	else if(rank[ru] > rank[rv]) pa[rv] = ru;
	else {
		pa[rv] = ru;
		++ rank[ru];
	}
}


void Graph::get_eps(const char *eps_s, int &eps_a, int &eps_b) {
	int i = 0;
	eps_a = 0; eps_b = 1;
	while(eps_s[i] != '\0'&&eps_s[i] != '.') {
		eps_a = eps_a*10 + (eps_s[i]-'0');
		++ i;
	}

	if(eps_s[i] == '.') {
		++ i;
		while(eps_s[i] != '\0') {
			eps_a = eps_a*10 + (eps_s[i]-'0');
			eps_b *= 10;
			++ i;
		}
	}
}
