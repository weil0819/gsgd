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
	reverse = NULL;
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
	if(reverse != NULL) {
		delete[] reverse;
		reverse = NULL;
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
	if(reverse == NULL) reverse = new ui[m];
	// if(cid == NULL) cid = new int[n];
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
	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu && invalid_cluster[cid[pa[i]]] == 0) {
		fprintf(fout, "c %d %d\n", i, cid[pa[i]]);
	}

	sort(noncore_cluster.begin(), noncore_cluster.end());
	noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
	for(ui i = 0;i < noncore_cluster.size();i ++) {
		if(invalid_cluster[noncore_cluster[i].first] == 0) {
			fprintf(fout, "n %d %d\n", noncore_cluster[i].second, noncore_cluster[i].first);
		}	
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
		if(similar_degree[u] >= mu && invalid_cluster[cid[pa[u]]] == 0) {
			int cluster_id = cid[pa[u]];
			cluster_set[cluster_id].pb(u);
		}
	}

	// Consider non-core nodes.
	sort(noncore_cluster.begin(), noncore_cluster.end());
	noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end()); 

	for(ui i = 0;i < noncore_cluster.size();i ++) {
		if(invalid_cluster[noncore_cluster[i].first] == 0) {
			cluster_set[noncore_cluster[i].first].pb(noncore_cluster[i].second);
		}		
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


void Graph::gsgd(const char *eps_s, int miu, int dist) {
#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif
	// Pruning-1: prune all neighbors whose distance is larger than dist.
	// Pruning-2: Find edge (u,v) in N[v] using binary search, and build cross link between (u,v) and (v,u).
	prune_and_cross_link(dist);

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
			ui v = edges[i];

			// v is neighbor of u, we only check smaller ID node, that is, u < v.
			if(min_cn[i] < 0 || v < u) continue;		// min_cn[i] means the i-th position of edges 

			if(u == v) {	
				min_cn[i] = -2;
				continue ;
			} 

			min_cn[i] = naive_similar_check(u, v, eps_a2, eps_b2);
			min_cn[reverse[i]] = min_cn[i];				// half computation

			if(min_cn[i] == -1) {
				++similar_degree[u];	// compute N_eps[u] 
				++similar_degree[v];	// compute N_eps[v] 
			}
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
			if(u == v || similar_degree[v] < miu) continue;	// v is marked node or v is not a core node 
			// If v is also a core node and sim(u,v)>=eps, merge them into one cluster.
			if(min_cn[i] == -1) my_union(u,v);
		}
	}

	// Pruning-2: prune all core clusters whose MCC is larger than dist.
	invalid_cluster.clear();
	invalid_cluster.reserve(n);


	// Step-5: Cluster non-core nodes.
	cluster_noncore_vertices(eps_a2, eps_b2, miu);

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#endif

	
	// Step-6: Compute MCC for each cluster and return qualified subgraph.
	vector<vector<int> > cluster_set(n);		// cluster_set[i] is the i-th cluster

	// Consider non-core nodes.
	sort(noncore_cluster.begin(), noncore_cluster.end());
	noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end()); 

	for(ui i = 0;i < noncore_cluster.size();i ++) {
		cluster_set[noncore_cluster[i].first].pb(noncore_cluster[i].second);
	}

	int CID = 0;
	for(int c = 0; c < n; c++) {
		bool valid = true;
		if(cluster_set[c].size() > 0) ++CID;
		if(cluster_set[c].size() > 2) {			// only non-core nodes can lay on the boundary of the MCC
			for(int i = 2; i < cluster_set[c].size(); i++) {	// v1
				ui v1 = cluster_set[c][i];
				for(int j = 0; j < i; j++) {					// v2
					ui v2 = cluster_set[c][j];
					for(int k = j + 1; k < i; k++) {			// v3
						ui v3 = cluster_set[c][k];
						double radius2 = mcc_radius2(v1, v2, v3);
						if(radius2 < (double)dist*dist) continue ;
						invalid_cluster[c] = 1;
						valid = false;
						break;
					}
					if(!valid) break;
				}
				if(!valid) break;
			}
			// if(valid) ++CID;
		}
	}
	printf("CID = %d\n", CID);

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("Clustering time: %lld, MMC time: %lld, Total time: %lld\n", mtime1, mtime-mtime1, mtime);
#endif
	
}


void Graph::dgcd(const char *eps_s, int miu, int dist) {

	if(cid == NULL) cid = new int[n];		// assign each node a cluster ID 
	for(ui u = 0; u < n; u++) cid[u] = n;	// initialize cluster ID

	// Step-1: Compute eps value (including eps_a2 and eps_b2).
	int eps_a2 = 0, eps_b2 = 0;
	get_eps(eps_s, eps_a2, eps_b2);
	eps_a2 *= eps_a2;
	eps_b2 *= eps_b2;

	// Step-2: Label all vertices as unprocessed -- false.
	bool *processed = new bool[n];

	// Step-3: Initialize queue Q, hash table H, cluster ID cid.
	queue<ui> Q;
	vector<pair<ui, ui> > H;
	int CID = 1;

	// Step-4: Iterate each unprocessed user v.
	for(ui v = 0; v < n; v++) {
		if(!processed[v]) {		// unprocessed user -- v

			// Step 4-1: Compute epsilon-neighbor of v.
			vector<ui> S;		// v's epsilon-neighbor list
			eps_neighbor(v, eps_a2, eps_b2, dist, H, S);

			// Step 4-2: Decide v is a core node or not.
			if(S.size() < miu) continue ;	// v is not a core mode.		
			
			// Ohterwise, v is a core node.
			cid[v] = CID;					// assign CID to v 
			processed[v] = true;
			
			// Step 4-3: Collect v's epsilon-neighbor.
			for(int i = 0; i < S.size(); i++) {
				ui u = S[i];
				cid[u] = CID;				
				if(!processed[u]) Q.push(u);
			}

			// Step 4-4: Scan each vertex in Q.
			while(!Q.empty()) {
				ui k = Q.front();
				Q.pop();

				if(!processed[k]) {
					S.clear();

					// Collect k's epsilon-neighbor list
					eps_neighbor(k, eps_a2, eps_b2, dist, H, S);

					if(S.size() < miu) continue ;  

					processed[k] = true;	// k is a core node

					// all k's epsilon-neighbor belong to same cluster
					for(int i = 0; i < S.size(); i++) {		
						ui h = S[i];
						cid[h] = CID;
						if(!processed[h]) Q.push(h);
					}
				}

			}

			// Step 4-5: Update cluster ID.			
			++CID;
		}
	}
	printf("CID = %d\n", CID);
}


void Graph::eps_neighbor(ui v, int eps_a2, int eps_b2, int dist, vector<pair<ui, ui> > &H, vector<ui> &S) {

	S.clear();		// S stores the epsilon-neighborhood of user -- v 

	// Iterate each neighbor of v and add it into N_gs(v).
	vector<ui> N_gs_v; 
	N_gs_v.pb(v); 
	for(ui i = pstart[v]; i < pstart[v+1]; i++) {
		ui u = edges[i];
		if(euclidean_dist2(v,u) <= dist*dist) {	// euclidean distance <= gamma
			N_gs_v.pb(u);
		}
	}

	// Iterate each N_gs(v) and decide which one is in Hash Table.
	for(int i = 0; i < N_gs_v.size(); i++) {
		ui v_i = v, v_j = N_gs_v[i];			// v_j is v_i's geo-social neighbor

		// (v_i,v_j) or (v_j,v_i) exists
		if(find(H.begin(), H.end(), mp(v_i, v_j)) != H.end() || find(H.begin(), H.end(), mp(v_j, v_i)) != H.end()) {
			S.pb(v_j);
		} 
		else { 
			// Compute N_gs(v_j).
			vector<ui> N_gs_j; 
			N_gs_j.pb(v_j); 

			for(ui j = pstart[v_j]; j < pstart[v_j+1]; j++) {
				ui u = edges[j];
				if(euclidean_dist2(v_j,u) <= dist*dist) {
					N_gs_j.pb(u);
				}
			}

			// Compute delta(v_i,v_j).
			vector<int> intersect;
			sort(N_gs_v.begin(),N_gs_v.end());   
			sort(N_gs_j.begin(),N_gs_j.end());   
			set_intersection(N_gs_v.begin(),N_gs_v.end(),N_gs_j.begin(),N_gs_j.end(),back_inserter(intersect));
			int cn = (int)intersect.size();

			// structural similarity(v_i,v_j) >= epsilon
			if(((long long)cn)*((long long)cn)*eps_b2 >= ((long long)N_gs_v.size())*((long long)N_gs_j.size())*eps_a2) {
				S.pb(v_j); 
				H.pb(mp(v_i,v_j)); // a pair of users -- (v_i,v_j), belong to each other's epsilon-neighborhood 
				H.pb(mp(v_j,v_i)); 
			}
		}
	}
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


double Graph::euclidean_dist2(int u, int v) {
	float ux = ploc[2*u], uy = ploc[2*u+1];
	float vx = ploc[2*v], vy = ploc[2*v+1];

	return (double)(ux-vx)*(ux-vx) + (uy-vy)*(uy-vy);
}


double Graph::euclidean_dist2(float ux, float uy, float vx, float vy) {
	return (double)(ux-vx)*(ux-vx) + (uy-vy)*(uy-vy);
}


double Graph::mcc_radius2(int v1, int v2, int v3) {
	double radius2 = 0.0;

	// Step-1: compute pairwise distance. 
	vector<double> dist2;
	dist2.reserve(3);
	dist2[0] = euclidean_dist2(v1, v2);
	dist2[1] = euclidean_dist2(v1, v3);
	dist2[2] = euclidean_dist2(v2, v3);

	// Step-2: judge the triangle (obtuse or right or acute) and find the center of the MCC.

	int indMax = 0;
	for(int i = 1; i < dist2.size(); i++) {
		if(dist2[i] > dist2[indMax]) indMax = i;
	}
	double part1 = dist2[indMax];	// longest_edges^2
	double part2 = 0.0;				// other_edges^2
	for(int i = 0; i < dist2.size(); i++) {
		if(i != indMax) part2 += dist2[i];
	}

	if(part1 >= part2) {			// obtuse or right triangle
		if(indMax == 0) {
			float centerX = (ploc[2*v1] - ploc[2*v2]) / 2;
			float centerY = (ploc[2*v1+1] - ploc[2*v2+1]) / 2;
			radius2 = euclidean_dist2(ploc[2*v1], ploc[2*v1+1], centerX, centerY);
		}else if(indMax == 1) {
			float centerX = (ploc[2*v2] - ploc[2*v3]) / 2;
			float centerY = (ploc[2*v2+1] - ploc[2*v3+1]) / 2;
			radius2 = euclidean_dist2(ploc[2*v2], ploc[2*v2+1], centerX, centerY);
		}else {
			float centerX = (ploc[2*v3] - ploc[2*v1]) / 2;
			float centerY = (ploc[2*v3+1] - ploc[2*v1+1]) / 2;
			radius2 = euclidean_dist2(ploc[2*v3], ploc[2*v3+1], centerX, centerY);
		}
	}else {							// acute triangle
		float a1 = ploc[2*v2] - ploc[2*v1], b1 = ploc[2*v2+1] - ploc[2*v1+1];
		float a2 = ploc[2*v3] - ploc[2*v1], b2 = ploc[2*v3+1] - ploc[2*v1+1];
		double c1 = (a1 * a1 + b1 * b1) / 2, c2 = (a2 * a2 + b2 * b2) / 2;
		double d = a1 * b2 - a2 * b1;
		float centerX = ploc[2*v1] + (c1 * b2 - c2 * b1) / d;
		float centerY = ploc[2*v1+1] + (a1 * c2 - a2 * c1) / d;
		radius2 = euclidean_dist2(ploc[2*v1], ploc[2*v1+1], centerX, centerY);
	}

	return radius2;
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


void Graph::prune_and_cross_link(int dist) {

	for(ui u = 0; u < n; u++) {
		for(ui i = pstart[u]; i < pstart[u+1]; i++) {
			ui v = edges[i];
			// v is neighbor of u, we only check smaller ID node, that is, u < v.
			if(v < u) {
				//if(min_cn[i] == 0) min_cn[i] = -2;
				continue; //this edge has already been checked
			}

			// u < v
			ui r_id = binary_search(edges, pstart[v], pstart[v+1], u);	// return u's position in v's adjacent array 
			reverse[i] = r_id;		// i is v's position in u's adjacent array
			reverse[r_id] = i;		// <position1, position2>
			
			double dist2 = euclidean_dist2(u,v);
			
			if(dist2 > (double)dist*dist) {
				// Mark v as empty and degree minus 1.
				edges[i] = u;		// mark u's neighbor same as u
				edges[r_id] = v;	// mark v's neighbor same as v
				--degree[u];
				--degree[v];
			}
		}
	}
}


ui Graph::binary_search(const int *array, ui b, ui e, int val) {
	-- e;

	if(array[e] < val) return e+1;
	while(b < e) {
		ui mid = b + (e-b)/2;
		if(array[mid] >= val) e = mid;
		else b = mid+1;
	}

	return e;
}



