/*
 * Synprep.cpp
 * 
 * Obv-1: n=30,000, m=300,000
 * Obv-2: n=400,000, m=4,000,000
 * 
 *  Created on: 2 Dec, 2019
 *      Author: Wesley
 */


#include "../Utility.h"
#include <ext/hash_map>

#define MAXBUF 8000000
#define hash_map __gnu_cxx::hash_map

using namespace std;


ui n, m;

vector<ui> nodes;
vector<pair<ui, ui> > edges;
vector<pair<ui, pair<double, double> > > node_loc;
hash_map<ui, ui> pairNodes;

void input(const char *dir, const char *file_name);
void output(const char *dir);

int main(int argc, char *argv[]){

	if(argc < 3){
		printf("Usage: ./Synprep dir file_name \n");
		exit(1);
	}

	input(argv[1], argv[2]);
	printf("Finish reading the graph \n");

	output(argv[1]);
	printf("Finish writing the graph \n");

	return 0;

}


void input(const char *dir, const char *file_name) {
	// Read file 
	string name = string(dir) + "/" + string(file_name);
	FILE *fin = open_file(name.c_str(), "r");

	// fscanf(fin, "%d %d", &n, &m);

	// printf("\tn = %u; m = %u\n", n, m);

	while(!feof(fin)){
		int a, b, c;
		fscanf(fin, "%d %d %s", &a, &b, &c);

		if(a == b){
			nodes.pb(a);
		}else{
			nodes.pb(a);
			nodes.pb(b);
			edges.pb(mp(a, b));
			edges.pb(mp(b, a));
		}
	}

	// Sort nodes/edges and remove duplicated ones 
	sort(nodes.begin(), nodes.end());
	nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());

	sort(edges.begin(), edges.end());
	edges.erase(unique(edges.begin(), edges.end()), edges.end());

	printf("\tnodes.size() = %u; edges.size()/2 = %u, degree = %.2f\n", (int)nodes.size(), (int)edges.size()/2, (double)edges.size()/nodes.size());

	fclose(fin);

}


void output(const char *dir) {
	// map from old to new 
	vector<ui> nodesVec;
	vector<pair<ui, ui> > edgesVec;

	for(ui i = 0; i < nodes.size(); i++) {
		pairNodes.insert(pair<ui, ui>(nodes[i], i));	// <old, new>
	}

	for(vector<ui>::iterator it1 = nodes.begin(); it1 != nodes.end(); ++it1){ 
		if(pairNodes.find(*it1) != pairNodes.end()) {
			nodesVec.pb(pairNodes.find(*it1)->second);
		}		
	}

	for(vector<pair<ui, ui> >::iterator it2 = edges.begin(); it2 != edges.end(); ++it2){
		if(pairNodes.find((*it2).first) != pairNodes.end() && pairNodes.find((*it2).second) != pairNodes.end() ) {
			edgesVec.pb(mp(pairNodes.find((*it2).first)->second, pairNodes.find((*it2).second)->second));
		}
	}

	sort(nodesVec.begin(), nodesVec.end());
	nodesVec.erase(unique(nodesVec.begin(), nodesVec.end()), nodesVec.end());

	sort(edgesVec.begin(), edgesVec.end());
	edgesVec.erase(unique(edgesVec.begin(), edgesVec.end()), edgesVec.end());

	n = (ui)nodesVec.size();
	m = (ui)edgesVec.size();

	printf("In new file, n = %d, m = %d \n", n, m);

	// gaussian distribution
	double mean = 20.0, std = 40.0;
	for(ui i = 0; i < n; i++) {
		// http://c.biancheng.net/view/643.html
		random_device rd;
		default_random_engine generator{rd()};
		normal_distribution<double> distribution(mean,std);

		// generate x
		double x = distribution(generator);
		while(x > 200 || x < 0) {
			x = distribution(generator);
		}

		// generate y
		double y = distribution(generator);
		while(y > 200 || y < 0) {
			y = distribution(generator);
		}

		node_loc.pb(mp(i, mp(x,y)));
	}

	string name1 = string(dir) + "b_degree.bin";
	string name2 = string(dir) + "b_adj.bin"; 

	FILE *fout1 = open_file(name1.c_str(), "wb");
	FILE *fout2 = open_file(name2.c_str(), "wb");

	ui tt = sizeof(ui);
	fwrite(&tt, sizeof(ui), 1, fout1);
	fwrite(&n, sizeof(ui), 1, fout1);
	fwrite(&m, sizeof(ui), 1, fout1); 

	ui *buf = new ui[n];
	ui j = 0;
	for(ui i = 0; i < n; i++){
		ui cnt = 0;
		while((edgesVec[j].first == i) && (j < (ui)edgesVec.size())){
			buf[cnt] = edgesVec[j].second;
			cnt++;
			j++;
		}
		fwrite(&cnt, sizeof(ui), 1, fout1);
		fwrite(buf, sizeof(ui), cnt, fout2);
	}

	delete[] buf;

	fclose(fout1);
	fclose(fout2);

	string name3 = string(dir) + "loc.txt";
	FILE *fout3 = open_file(name3.c_str(), "w");

	for(vector<pair<ui, pair<double, double> > >::iterator it = node_loc.begin(); it != node_loc.end(); ++it){
		fprintf(fout3, "%f %f", (*it).second.first, (*it).second.second);
		fputs("\n", fout3);
	}

	fclose(fout3);

}

