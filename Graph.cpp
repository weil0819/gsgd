/*
 * Graph.cpp
 *
 *  Created on: 02 Sept, 2019
 *      Author: Wesley
 */

#include "Utility.h"
#include "Graph.h"

bool myComp(const vector<int> &A, const vector<int> &B) {
	return A.size() > B.size();
}

const Circle Circle::INVALID{0, 0, -1};
static default_random_engine randGen((std::random_device())());

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


void Graph::cluster_noncore_vertices(int eps_a2, int eps_b2, int miu, unordered_set<int> &US) {
	if(cid == NULL) cid = new int[n];		// assign each node a cluster ID
	for(ui u = 0; u < n; u++) cid[u] = n;	// initialize cluster ID

	for(auto u: US) {
		if(similar_degree[u] >= miu) {		// u is a core node
			int x = find_root(u);
			if(u < cid[x]) cid[x] = u;		// update cluster ID, which should be minimal ID in cluster
		}
	}

	noncore_cluster.clear();
	noncore_cluster.reserve(n);

	for(auto i: US) {
		if(similar_degree[i] >= miu) {		// i should be core node
			for(ui j = pstart[i];j < pstart[i+1];j ++) {
				// when edges[j] is also a core node or i and edges[j] are in same cluster
				if(US.find(edges[j]) == US.end() || similar_degree[edges[j]] >= miu||pa[i] == pa[edges[j]]) continue;

				if(min_cn[j] >= 0) {
					min_cn[j] = naive_similar_check(i, edges[j], eps_a2, eps_b2, US);
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
			fprintf(fout, "cluster-%d, \t #cluster=%d \n", i, (int)cluster_set[i].size());
			for(int j = 0; j < cluster_set[i].size(); j++) {
				fprintf(fout, "%d\t", cluster_set[i][j]);
			}
			fprintf(fout, "\n************************************************************************\n");
		}
	}

	fclose(fout);

}


void Graph::gsgd(const char *eps_s, int miu, int dist) {

	// Pruning-1: prune all neighbors whose distance is larger than dist.
	// Pruning-2: Find edge (u,v) in N[v] using binary search, and build cross link between (u,v) and (v,u).
	prune_and_cross_link(dist);

	// Step-1: Compute eps value (including eps_a2 and eps_b2).
	int eps_a2 = 0, eps_b2 = 0;
	get_eps(eps_s, eps_a2, eps_b2);
	eps_a2 *= eps_a2;
	eps_b2 *= eps_b2;

#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif

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
		}
	}

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
	// printf("CID = %d\n", CID);
}


void Graph::baseline(const char *eps_s, int mu, int gamma) {

	int eps_a2 = 0, eps_b2 = 0;
	get_eps(eps_s, eps_a2, eps_b2);
	eps_a2 *= eps_a2;
	eps_b2 *= eps_b2;

#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif

	// Cluster-Tric
	ui *deg = new ui[n];
	for(ui i = 0;i < n;i ++) deg[i] = pstart[i+1]-pstart[i];

	ui *adj = new ui[n];
	memset(adj, 0, sizeof(ui)*n);

	ui *similar = new ui[m];
	memset(similar, 0, sizeof(ui)*m);

	ui *pend = new ui[n];		

	for(ui i = 0;i < n;i ++) { 

		ui &end = pend[i] = pstart[i]; 

		ui j = pstart[i+1]; 

		while(true) {

			while(end < j&&(deg[edges[end]] < deg[i]||(deg[edges[end]]==deg[i]&&edges[end]<i))) ++ end;

			while(j > end&&(deg[edges[j-1]] > deg[i]||(deg[edges[j-1]]==deg[i]&&edges[j-1]>i))) -- j;

			if(end >= j) break;

			swap(edges[end], edges[j-1]);
		}
		sort(edges+pend[i], edges+pstart[i+1]);
	}
	
	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;

		for(ui j = pstart[u];j < pend[u];j ++) {
			ui v = edges[j];

			for(ui k = pstart[v];k < pend[v];k ++) if(adj[edges[k]]) {
				++ similar[j];		// edges[j]
				++ similar[k];		// edges[k]
				++ similar[adj[edges[k]] - 1];
			}
		}

		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}

	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pend[u];j ++) {
			ui v = edges[j];

			similar[j] += 2;

			if(((long long)similar[j])*((long long)similar[j])*eps_b2 >= ((long long)(deg[u]+1))*((long long)(deg[v]+1))*eps_a2) similar[j] = 1;
			else similar[j] = 0;

			ui r_id = binary_search(edges+pend[v], pstart[v+1]-pend[v], u) + pend[v];
			similar[r_id] = similar[j];
		}
	}

	delete[] pend; pend = NULL;
	delete[] deg; deg = NULL;
	delete[] adj; adj = NULL;

	if(similar_degree == NULL) similar_degree = new int[n];
	memset(similar_degree, 0, sizeof(int)*n);

	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pstart[i+1];j ++) {
		if(similar[j] == 1) ++ similar_degree[i];
	}

	if(pa == NULL) pa = new int[n];
	if(rank == NULL) rank = new int[n];
	for(ui i = 0;i < n;i ++) {
		pa[i] = i;
		rank[i] = 0;
	}

	for(ui i = 0;i < n;i ++) {
		if(similar_degree[i] < mu) continue;

		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(similar_degree[edges[j]] < mu) continue;		

			if(similar[j] == 1) my_union(pa, rank, i, edges[j]);
		}
	}

	if(cid == NULL) cid = new int[n];
	for(ui i = 0;i < n;i ++) cid[i] = n;

	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		int x = find_root(pa, i);
		if(i < cid[x]) cid[x] = i;
	}

	vector<pair<int,int> > noncore_cluster;
	noncore_cluster.reserve(n);

	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(similar_degree[edges[j]] >= mu) continue;

			if(similar[j] == 1) noncore_cluster.pb(mp(cid[pa[i]], edges[j]));
		}
	}

	delete[] similar; similar = NULL;

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#endif

	// Collect a set of clusters ... 
	vector<vector<int> > clusters(n);

	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		clusters[cid[find_root(pa, i)]].pb(i);
	}

	sort(noncore_cluster.begin(), noncore_cluster.end());
	noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
	for(ui i = 0;i < noncore_cluster.size();i ++) {
		clusters[noncore_cluster[i].first].pb(noncore_cluster[i].second);
	}

	// Iterate each cluster ... 
	for(ui c = 0; c < n; c++) {
		if(clusters[c].size() < mu) continue ;
		printf("============================ORIGINAL CLUSTER: %d has %d nodes==============================\n", c, (int)clusters[c].size());
		vector<vector<int> > final_clusters;
		vector<int> preList;
		for(int i = 2; i < clusters[c].size(); i++) {
			ui v_i = clusters[c][i];				// Vi
			for(int j = 0; j < i; j++) {
				ui v_j = clusters[c][j];			// Vj
				for(ui k = j+1; k < i; k++) {
					ui v_k = clusters[c][k];		// Vk
					double centerX = 0.0, centerY = 0.0;
					double r2 = mcc_radius2(v_i, v_j, v_k, centerX, centerY);	// r^2					

					if(r2 < (double)gamma*gamma) {
						vector<int> curList;		// a set of vertices in circle
						for(ui ii = 0; ii < clusters[c].size(); ii++) {
							double dist2 = euclidean_dist2(ploc[2*clusters[c][ii]], ploc[2*clusters[c][ii]+1], centerX, centerY);
							if(dist2 <= r2) curList.pb(clusters[c][ii]);
						}
						
						// Compute cluster based on "curList"
						if(curList.size() < mu) continue; 

						int isSub = 1;
						if(preList.empty()) preList = curList;
						if(!is_subset(curList, preList) && !is_subset(preList, curList)) isSub = 0; 

						if((isSub && curList.size() > preList.size()) || !isSub) {
							preList = curList;

							vector<vector<int> > cur_clusters;
							renew_cluster(eps_a2, eps_b2, mu, curList, cur_clusters);

							final_clusters.insert(final_clusters.end(), cur_clusters.begin(), cur_clusters.end());	
						} 
					}
				}
			}
		}
		// Printing ... 		
		if(!final_clusters.empty()) { 			
			sort(final_clusters.begin(), final_clusters.end(), myComp);
			for(int i = 0; i < final_clusters.size(); i++) {
				int j = i-1;
				while(j >= 0 && !is_subset(final_clusters[j], final_clusters[i])) --j;
				if(j != -1) continue ;

				printf("cluster-%d \t #cluster=%d \n", i, (int)final_clusters[i].size());
				for(auto u: final_clusters[i]) printf("%d\t", u);
					printf("\n************************************************************************\n");
				break;
			}			
		}
	}

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("Clustering time: %lld, MMC time: %lld, Total time: %lld\n", mtime1, mtime-mtime1, mtime);
#endif	
}


void Graph::advance(const char *eps_s, int mu, int gamma) {

	// Get eps value ... 
	int eps_a2 = 0, eps_b2 = 0;
	get_eps(eps_s, eps_a2, eps_b2);
	eps_a2 *= eps_a2;
	eps_b2 *= eps_b2;

#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif

	// Cluster-Tric ... 
	ui *deg = new ui[n];
	for(ui i = 0;i < n;i ++) deg[i] = pstart[i+1]-pstart[i];

	ui *adj = new ui[n];
	memset(adj, 0, sizeof(ui)*n);

	ui *similar = new ui[m];
	memset(similar, 0, sizeof(ui)*m);

	ui *pend = new ui[n];		

	for(ui i = 0;i < n;i ++) { 

		ui &end = pend[i] = pstart[i]; 

		ui j = pstart[i+1]; 

		while(true) {

			while(end < j&&(deg[edges[end]] < deg[i]||(deg[edges[end]]==deg[i]&&edges[end]<i))) ++ end;

			while(j > end&&(deg[edges[j-1]] > deg[i]||(deg[edges[j-1]]==deg[i]&&edges[j-1]>i))) -- j;

			if(end >= j) break;

			swap(edges[end], edges[j-1]);
		}
		sort(edges+pend[i], edges+pstart[i+1]);
	}
	
	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;

		for(ui j = pstart[u];j < pend[u];j ++) {
			ui v = edges[j];

			for(ui k = pstart[v];k < pend[v];k ++) if(adj[edges[k]]) {
				++ similar[j];		// edges[j]
				++ similar[k];		// edges[k]
				++ similar[adj[edges[k]] - 1];
			}
		}

		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}

	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pend[u];j ++) {
			ui v = edges[j];

			similar[j] += 2;

			if(((long long)similar[j])*((long long)similar[j])*eps_b2 >= ((long long)(deg[u]+1))*((long long)(deg[v]+1))*eps_a2) similar[j] = 1;
			else similar[j] = 0;

			ui r_id = binary_search(edges+pend[v], pstart[v+1]-pend[v], u) + pend[v];
			similar[r_id] = similar[j];
		}
	}

	delete[] pend; pend = NULL;
	delete[] deg; deg = NULL;
	delete[] adj; adj = NULL;

	if(similar_degree == NULL) similar_degree = new int[n];
	memset(similar_degree, 0, sizeof(int)*n);

	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pstart[i+1];j ++) {
		if(similar[j] == 1) ++ similar_degree[i];
	}

	// Cluter core nodes ... 
	if(pa == NULL) pa = new int[n];
	if(rank == NULL) rank = new int[n];
	for(ui i = 0;i < n;i ++) {
		pa[i] = i;
		rank[i] = 0;
	}

	for(ui i = 0;i < n;i ++) {
		if(similar_degree[i] < mu) continue;

		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(similar_degree[edges[j]] < mu) continue;		

			if(similar[j] == 1) my_union(pa, rank, i, edges[j]);
		}
	}

	if(cid == NULL) cid = new int[n];
	for(ui i = 0;i < n;i ++) cid[i] = n;

	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		int x = find_root(pa, i);
		if(i < cid[x]) cid[x] = i;
	}

	// Cluster non-core nodes ... 
	vector<pair<int,int> > noncore_cluster;
	noncore_cluster.reserve(n);

	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(similar_degree[edges[j]] >= mu) continue;

			if(similar[j] == 1) noncore_cluster.pb(mp(cid[pa[i]], edges[j]));
		}
	}

	delete[] similar; similar = NULL;

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#endif

	// Collect a set of clusters ... 
	vector<vector<int> > clusters(n);

	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		clusters[cid[find_root(pa, i)]].pb(i);
	}

	sort(noncore_cluster.begin(), noncore_cluster.end());
	noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
	for(ui i = 0;i < noncore_cluster.size();i ++) {
		clusters[noncore_cluster[i].first].pb(noncore_cluster[i].second);
	}

	// Iterate each cluster ... 
	for(ui c = 0; c < n; c++) {
		if(clusters[c].size() < mu) continue ;

		// printf("======================ORIGINAL CLUSTER: %d has %d nodes========================\n", c, (int)clusters[c].size());
		// for(int i = 0; i < clusters[c].size(); i++) {
		// 	printf("%d\t", clusters[c][i]);
		// }
		// printf("\n");

		// Is there a way to compute MCC for such a cluster?
		// then decide radius/diameter 
		// https://www.nayuki.io/page/smallest-enclosing-circle
		// https://www.codeproject.com/Articles/1165267/Coding-Challenge-Smallest-Circle-Problem 
		// https://blog.csdn.net/niiick/article/details/89153096 
		// https://blog.csdn.net/wu_tongtong/article/details/79362339 (*)
		
		Circle C = make_mcc(clusters[c]);
		// printf("==============ORIGINAL CLUSTER: %d has %d nodes===============\n", c, (int)clusters[c].size());
		// printf("C.centerX = %.2f, C.centerY= %.2f, C.radius = %.2f \n", C.centerX, C.centerY, C.radius);

		if(C.radius > gamma) {
			printf("============================ORIGINAL CLUSTER: %d has %d nodes==============================\n", c, (int)clusters[c].size());
			printf("\t*** centerX = %.2f, centerY= %.2f, radius = %.2f \n", C.centerX, C.centerY, C.radius);
			// for(int i = 0; i < clusters[c].size(); i++) {
			// 	printf("%d\t", clusters[c][i]);
			// }
			// printf("\n");

			// reduce clusters[c] to make it as a cluster in MCC
			// vector<int> curList;
			// for(int i = 0; i < clusters[c].size(); i++) {
			// 	ui u = clusters[c][i];

			// 	double dist = make_distance(C.centerX, C.centerY, ploc[2*u], ploc[2*u+1]);
			// 	if(dist <= gamma) {
			// 		curList.pb(u);
			// 	}
			// }
			// printf("\t*** There are %d vertices in circle \n", (int)curList.size());

			// vector<vector<int> > output;
			// renew_cluster(eps_a2, eps_b2, mu, curList, output);

			// if(output.empty()) {
			// 	printf("\t*** new cluster is empty \n");
			// 	continue ;
			// }
			// for(int i = 0; i < output.size(); i++) {
			// 	printf("\t*** new cluster has %d nodes \n", (int)output[i].size());
			// 	for(int j = 0; j < output[i].size(); j++) {
			// 		printf("%d\t", output[i][j]);
			// 	}
			// 	printf("\n");
			// }

		}

		// utilize Floyd-Warshall algorithm to calculate shortest-path between any two nodes in clusters[c].
		// calculate e(v) and d(G), nodes with e(v) > gamma can be considered as candidate to drop
		// vector<int> output;
		// floyd_diameter(eps_a2, eps_b2, mu, gamma, clusters[c], output);	

		// vector<vector<int> > output;
		// renew_cluster(eps_a2, eps_b2, mu, clusters[c], output);
		// if(output.empty()) continue ;
		// for(int i = 0; i < output.size(); i++) {
		// 	printf("new cluster has %d nodes \n", (int)output[i].size());
		// 	for(int j = 0; j < output[i].size(); j++) {
		// 		printf("%d\t", output[i][j]);
		// 	}
		// 	printf("\n");
		// }
		
	}

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("Clustering time: %lld, MMC time: %lld, Total time: %lld\n", mtime1, mtime-mtime1, mtime);
#endif
}


/*
// long-time cost renew implementation (keep it as backup)
void Graph::renew_cluster(int eps_a2, int eps_b2, int mu, vector<int> &curList) {
	unordered_set<int> US(curList.begin(), curList.end());
	
	if(degree == NULL) degree = new int[n];
	memset(degree, 0, sizeof(int)*n);

	for(auto i: curList) {
		for(ui j = pstart[i]; j < pstart[i+1]; j++) {
			if(US.find(edges[j]) != US.end()) ++degree[i];
		}
	}

	if(similar_degree == NULL) similar_degree = new int[n];
	memset(similar_degree, 0, sizeof(int)*n);

	if(min_cn == NULL) min_cn = new int[m];
	memset(min_cn, 0, sizeof(int)*m);

	int flag = 0;
	for(auto i: curList) {
		if(degree[i] < mu) {			
			similar_degree[i] = -1;		 
			continue ;
		}
		flag = 1;
		for(ui j = pstart[i]; j < pstart[i+1]; j++) {
			if(min_cn[j] < 0 || US.find(edges[j]) == US.end()) continue;
			min_cn[j] = naive_similar_check(i, edges[j], eps_a2, eps_b2, US);
			if(min_cn[j] == -1) ++similar_degree[i];
		}
	}
	if(!flag) return ;

	if(pa == NULL) pa = new int[n];
	if(rank == NULL) rank = new int[n];
	for(ui i = 0;i < n;i ++) {
		pa[i] = i;
		rank[i] = 0;
	}
	flag = 0;
	for(auto i: curList) {								
		if(similar_degree[i] < mu) continue;	
		flag = 1;			
		for(ui j = pstart[i];j < pstart[i+1];j ++) {		
			if(US.find(edges[j]) == US.end() || similar_degree[edges[j]] < mu) continue;	
			if(min_cn[j] == -1) my_union(i, edges[j]);		
		}
	}
	if(!flag) return ;

	cluster_noncore_vertices(eps_a2, eps_b2, mu, US);

}
*/

void Graph::renew_cluster(int eps_a2, int eps_b2, int mu, vector<int> &curList, vector<vector<int> > &cluster) {
	// Redefine "n,m,pstart,edges"
	unordered_map<int, int> UM;	
	unordered_map<int, int> UM_R;	

	ui nn = curList.size();
	ui mm = 0;
	for(ui i = 0; i < nn; i++) {
		UM[curList[i]] = i;
		UM_R[i] = curList[i];
	}

	ui *deg = new ui[nn];
	memset(deg, 0, sizeof(ui)*nn);

	for(ui i = 0; i < nn; i++) {
		ui u = curList[i];
		for(ui j = pstart[u]; j < pstart[u+1]; j++) {
			if(UM.find(edges[j]) != UM.end()) ++deg[i];
		}
		mm += deg[i];
	}

	ui *ps = new ui[nn+1];
	int *ed = new int[mm];

	ui x = 0;
	ps[0] = 0;
	for(ui i = 0; i < nn; i++) {
		if(deg[i] > 0) {
			ui u = curList[i];
			for(ui j = pstart[u]; j < pstart[u+1]; j++) {
				if(UM.find(edges[j]) != UM.end()) ed[x++] = UM[edges[j]];	
			}
		}
		ps[i+1] = ps[i] + deg[i];
	}

	// Cluster-Tric
	ui *adj = new ui[nn];
	memset(adj, 0, sizeof(ui)*nn);

	ui *similar = new ui[mm];
	memset(similar, 0, sizeof(ui)*mm);

	ui *pend = new ui[nn];

	for(ui i = 0; i < nn; i++) {
		ui &end = pend[i] = ps[i];
		ui j = ps[i+1];

		while(true) {
			while(end < j&&(deg[ed[end]] < deg[i]||(deg[ed[end]]==deg[i]&&ed[end]<i))) ++ end;

			while(j > end&&(deg[ed[j-1]] > deg[i]||(deg[ed[j-1]]==deg[i]&&ed[j-1]>i))) -- j; 

			if(end >= j) break;

			swap(ed[end], ed[j-1]);
		}
		sort(ed+pend[i], ed+ps[i+1]);
	}

	for(ui u = 0;u < nn;u ++) {
		for(ui j = ps[u];j < pend[u];j ++) adj[ed[j]] = j+1;

		for(ui j = ps[u];j < pend[u];j ++) {
			ui v = ed[j];

			for(ui k = ps[v];k < pend[v];k ++) if(adj[ed[k]]) {
				++ similar[j];		// ed[j]
				++ similar[k];		// ed[k]
				++ similar[adj[ed[k]] - 1];
			}
		}

		for(ui j = ps[u];j < pend[u];j ++) adj[ed[j]] = 0;
	}

	for(ui u = 0; u < nn; u ++) {
		for(ui j = ps[u];j < pend[u];j ++) {
			ui v = ed[j];

			similar[j] += 2;

			if(((long long)similar[j])*((long long)similar[j])*eps_b2 >= ((long long)(deg[u]+1))*((long long)(deg[v]+1))*eps_a2) similar[j] = 1;
			else similar[j] = 0;

			ui r_id = binary_search(ed+pend[v], ps[v+1]-pend[v], u) + pend[v];
			similar[r_id] = similar[j];
		}
	}

	delete[] pend; pend = NULL;
	delete[] adj; adj = NULL;
	delete[] deg; deg = NULL;

	int *similar_degree = new int[nn];
	memset(similar_degree, 0, sizeof(int)*nn);

	for(ui i = 0;i < nn;i ++) for(ui j = ps[i];j < ps[i+1];j ++) {
		if(similar[j] == 1) ++ similar_degree[i];
	}

	int *pa = new int[nn];
	int *rank = new int[nn];
	for(ui i = 0;i < nn;i ++) {
		pa[i] = i;
		rank[i] = 0;
	}

	int flag = 0;
	for(ui i = 0;i < nn;i ++) {
		if(similar_degree[i] < mu) continue;
		flag = 1;
		for(ui j = ps[i];j < ps[i+1];j ++) {
			if(similar_degree[ed[j]] < mu) continue;		
			if(similar[j] == 1) my_union(pa, rank, i, ed[j]);
		}
	}

	// Special Case -- no core nodes/clusters.
	if(!flag) {
		delete[] ps; ps = NULL;
		delete[] ed; ed = NULL;
		delete[] pa; pa = NULL;
		delete[] rank; rank = NULL;
		delete[] similar; similar = NULL;
		delete[] similar_degree; similar_degree = NULL;

		return ;
	}

	// General Case -- core nodes existing.
	int *cid = new int[nn];
	for(ui i = 0;i < nn;i ++) cid[i] = nn;

	for(ui i = 0;i < nn;i ++) if(similar_degree[i] >= mu) {
		int x = find_root(pa, i);
		if(i < cid[x]) cid[x] = i;
	}

	vector<pair<int,int> > noncore_cluster;
	noncore_cluster.reserve(nn);

	for(ui i = 0;i < nn;i ++) if(similar_degree[i] >= mu) {
		for(ui j = ps[i];j < ps[i+1];j ++) {
			if(similar_degree[ed[j]] >= mu) continue;
			if(similar[j] == 1) noncore_cluster.pb(mp(cid[pa[i]], ed[j]));
		}
	}

	// Printing ...	maybe more than one clusters in each triangle.
	vector<vector<int> > cluster_set(nn);

	for(ui u = 0; u < nn; u++) {
		if(similar_degree[u] >= mu) {
			int cluster_id = cid[pa[u]];
			cluster_set[cluster_id].pb(u);
		}
	}

	// Special Case -- no non-core nodes.
	if(!noncore_cluster.empty()) {
		sort(noncore_cluster.begin(), noncore_cluster.end());
		noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end()); 

		for(ui i = 0;i < noncore_cluster.size();i ++) {
			cluster_set[noncore_cluster[i].first].pb(noncore_cluster[i].second);		
		}
	}

	for(ui i = 0; i < nn; i++) {
		if(cluster_set[i].size() > 0) {
			vector<int> new_cluster;		

			for(auto u: cluster_set[i]) {
				new_cluster.pb(UM_R[u]);				
			}
			sort(new_cluster.begin(), new_cluster.end());
			cluster.pb(new_cluster);
		}
	}

	delete[] ps; ps = NULL;
	delete[] ed; ed = NULL;
	delete[] cid; cid = NULL;
	delete[] pa; pa = NULL;
	delete[] rank; rank = NULL;
	delete[] similar; similar = NULL;
	delete[] similar_degree; similar_degree = NULL;

}


void Graph::floyd_diameter(int eps_a2, int eps_b2, int mu, int gamma, vector<int> &cluster, vector<int> &output) {
	unordered_map<int, int> node_map;
	unordered_map<int, int> node_map_re;
	int N = cluster.size();
	for(int i = 0; i < N; i++) {
		node_map[i] = cluster[i];
		node_map_re[cluster[i]] = i;
	}

	vector<vector<double> > dist(N, vector<double>(N,INF));

	// Compute the weight for each edge
	for(ui i = 0; i < N; i++) {
		ui u = cluster[i];
		dist[i][i] = 0.0;
		for(ui j = pstart[u]; j < pstart[u+1]; j++) {
			ui v = edges[j];
			if(node_map_re.find(v) != node_map_re.end()) {
				int k = node_map_re[v];
				dist[i][k] = sqrt(euclidean_dist2(u,v));
				dist[k][i] = dist[i][k];
				printf("dist[%d][%d] = %.2f \n", i,k, dist[i][k]);
			}
		}
	}

	vector<double> ecc(N, 0.0);
	for(int k = 0; k < N; k++) {
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				if (dist[i][k] + dist[k][j] < dist[i][j]) dist[i][j] = dist[i][k] + dist[k][j];                    
			}			
		}
	}

	// ecc > 2gamma, ==> candidate nodes
	int cnt = 0;
	for(int i = 0; i < N; i++) {
		ecc[i] = *max_element(dist[i].begin(), dist[i].end());
	}
	sort(ecc.begin(), ecc.end());
}


void Graph::reduce_cluster(int eps_a2, int eps_b2, int mu, int gamma, vector<int> &cluster, vector<int> &output) {

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


int Graph::naive_similar_check(int u, int v, int eps_a2, int eps_b2,  unordered_set<int> &US) {
	int du = degree[u], dv = degree[v];		

	ui i = pstart[u], j = pstart[v];
	int cn = 2;	
	while(i < pstart[u+1]&&j < pstart[v+1]) { 

		if(US.find(edges[i]) == US.end()) ++i;
		if(US.find(edges[j]) == US.end()) ++j;

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


double Graph::mcc_radius2(int v1, int v2, int v3, double &centerX, double &centerY) {
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
			centerX = (ploc[2*v1] - ploc[2*v2]) / 2;
			centerY = (ploc[2*v1+1] - ploc[2*v2+1]) / 2;
			radius2 = euclidean_dist2(ploc[2*v1], ploc[2*v1+1], centerX, centerY);
		}else if(indMax == 1) {
			centerX = (ploc[2*v2] - ploc[2*v3]) / 2;
			centerY = (ploc[2*v2+1] - ploc[2*v3+1]) / 2;
			radius2 = euclidean_dist2(ploc[2*v2], ploc[2*v2+1], centerX, centerY);
		}else {
			centerX = (ploc[2*v3] - ploc[2*v1]) / 2;
			centerY = (ploc[2*v3+1] - ploc[2*v1+1]) / 2;
			radius2 = euclidean_dist2(ploc[2*v3], ploc[2*v3+1], centerX, centerY);
		}
	}else {							// acute triangle
		float a1 = ploc[2*v2] - ploc[2*v1], b1 = ploc[2*v2+1] - ploc[2*v1+1];
		float a2 = ploc[2*v3] - ploc[2*v1], b2 = ploc[2*v3+1] - ploc[2*v1+1];
		double c1 = (a1 * a1 + b1 * b1) / 2, c2 = (a2 * a2 + b2 * b2) / 2;
		double d = a1 * b2 - a2 * b1;
		centerX = ploc[2*v1] + (c1 * b2 - c2 * b1) / d;
		centerY = ploc[2*v1+1] + (a1 * c2 - a2 * c1) / d;
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


int Graph::find_root(int *pa, int u) {
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


void Graph::my_union(int *pa, int *rank, int u, int v) {
	int ru = find_root(pa, u);
	int rv = find_root(pa, v);

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


ui Graph::binary_search(const int *array, ui e, int val) {
	ui b = 0;
	-- e;
	while(b < e) {
		ui mid = b + (e-b)/2;
		if(array[mid] >= val) e = mid;
		else b = mid+1;
	}
	return e;
}


bool Graph::is_subset(vector<int> &A, vector<int> &B) {
	sort(A.begin(), A.end());
	sort(B.begin(), B.end());

	return includes(A.begin(), A.end(), B.begin(), B.end());
}


// ****** Circle Functions ******
Circle Graph::make_mcc(vector<int> &cluster) { 
	// shuffle original cluster
	// unsigned seed = chrono::system_clock::now().time_since_epoch().count();	
	// shuffle(cluster.begin(), cluster.end(), default_random_engine(seed));
	shuffle(cluster.begin(), cluster.end(), randGen);

	// progressively add vertices to MCC or recompute MCC
	Circle C = Circle::INVALID;
	for(ui i = 0; i < cluster.size(); i++) {
		ui p = cluster[i];
		if(C.radius < 0 || !make_contains(C.centerX, C.centerY, C.radius, p)) {
			C = make_circle_one_point(cluster, i+1, p);
		}
	}
	return C;
}

// One boundary vertex known
Circle Graph::make_circle_one_point(vector<int> &cluster, int end, ui p) {
	Circle C{ploc[2*p], ploc[2*p+1], 0};

	for(int i = 0; i < end; i++) {
		ui q = cluster[i];
		if(!make_contains(C.centerX, C.centerY, C.radius, q)) {
			if(C.radius == 0) C = make_diameter(p,q);				// circle with diameter = |p,q|
			else C = make_circle_two_point(cluster, i+1, p, q);		// circle with triangle (p,q,?)
		}
	}
	return C;
}

// Two boundary vertices known
Circle Graph::make_circle_two_point(vector<int> &cluster, int end, ui p, ui q) {
	Circle curC = make_diameter(p,q);
	Circle left  = Circle::INVALID;
	Circle right = Circle::INVALID;

	// For each vertex not in the two-point circle
	double qpX = ploc[2*q] - ploc[2*p];
	double qpY = ploc[2*q+1] - ploc[2*p+1];

	for(int i = 0; i < end; i++) {
		ui r = cluster[i];
		if(make_contains(curC.centerX, curC.centerY, curC.radius, r)) continue ;

		// Form a circumcircle and classify it on left or right side
		double cross = make_cross(qpX, qpY, ploc[2*r]-ploc[2*p], ploc[2*r+1]-ploc[2*p+1]);

		Circle C = make_circumcircle(p, q, r);

		if (C.radius < 0)
			continue;
		else if (cross > 0 && 
			(left.radius < 0 || make_cross(qpX, qpY, C.centerX-ploc[2*p], C.centerY-ploc[2*p+1]) > make_cross(qpX, qpY, left.centerX-ploc[2*p], left.centerY-ploc[2*p+1])))
			left = C;
		else if (cross < 0 && 
			(right.radius < 0 || make_cross(qpX, qpY, C.centerX-ploc[2*p], C.centerY-ploc[2*p+1]) < make_cross(qpX, qpY, right.centerX-ploc[2*p], right.centerY-ploc[2*p+1])))
			right = C;
	}

	// Select which circle to return
	if (left.radius < 0 && right.radius < 0)
		return curC;
	else if (left.radius < 0)
		return right;
	else if (right.radius < 0)
		return left;
	else
		return left.radius <= right.radius ? left : right;

}

// Two vertices circle
Circle Graph::make_diameter(ui u, ui v) {
	double centerX = (ploc[2*u] + ploc[2*v]) / 2;
	double centerY = (ploc[2*u+1] + ploc[2*v+1]) / 2;

	return Circle{centerX, centerY, max(make_distance(centerX, centerY, ploc[2*u], ploc[2*u+1]), make_distance(centerX, centerY, ploc[2*v], ploc[2*v+1]))};
}

// Three vertices circle
Circle Graph::make_circumcircle(ui a, ui b, ui c) {
	// Mathematical algorithm from Wikipedia: Circumscribed circle
	double ox = (min(min(ploc[2*a], ploc[2*b]), ploc[2*c]) + max(min(ploc[2*a], ploc[2*b]), ploc[2*c])) / 2;
	double oy = (min(min(ploc[2*a+1], ploc[2*b+1]), ploc[2*c+1]) + max(min(ploc[2*a+1], ploc[2*b+1]), ploc[2*c+1])) / 2;

	double ax = ploc[2*a] - ox;
	double ay = ploc[2*a+1] - oy;
	double bx = ploc[2*b] - ox; 
	double by = ploc[2*b+1] - oy;
	double cx = ploc[2*c] - ox; 
	double cy = ploc[2*c+1] - oy;

	double d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2;

	if (d == 0)
		return Circle::INVALID;

	double x = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d;
	double y = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d;

	double centerX = ox + x; 
	double centerY = oy + y;

	double radius = max(max(make_distance(centerX, centerY, ploc[2*a], ploc[2*a+1]), make_distance(centerX, centerY, ploc[2*b], ploc[2*b+1])), 
		make_distance(centerX, centerY, ploc[2*c], ploc[2*c+1]));

	return Circle{centerX, centerY, radius};
}

double Graph::make_cross(double px, double py, double qx, double qy) {
	return px * qy - py * qx;
}

double Graph::make_distance(double px, double py, double qx, double qy) {
	return hypot(px - qx, py - qy);
}

bool Graph::make_contains(double centerX, double centerY, double radius, ui p) {
	return make_distance(centerX, centerY, ploc[2*p], ploc[2*p+1]) <= radius * MULTIPLICATIVE_EPSILON;
}



