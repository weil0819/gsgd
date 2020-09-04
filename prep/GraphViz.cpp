/*
 * GraphViz.cpp
 * 
 * Obv-1: n=58228, m=428156
 * Obv-2: not all nodes have locations
 * Obv-3: re-order node ID
 * 
 *  Created on: 31 Jan, 2020
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
vector<pair<ui, pair<float, float> > > node_loc;
hash_map<ui, ui> pairNodes;

void input(const char *dir, const char *edges_file, const char *checkins_file);
void output(const char *dir);

int main(int argc, char *argv[]){

	if(argc < 4){
		printf("Usage: ./GraphViz dir edges_file checkins_file \n");
		exit(1);
	}

	input(argv[1], argv[2], argv[3]);
	printf("Finish reading the graph \n");

	output(argv[1]);
	printf("Finish writing the graph \n");

	return 0;

}


void input(const char *dir, const char *edges_file, const char *checkins_file){

	// Read-1: read edges file first.
	string name1 = string(dir) + "/" + string(edges_file);
	FILE *fin = open_file(name1.c_str(), "r");

	while(!feof(fin)){
		ui a, b;
		fscanf(fin, "%d %d", &a, &b);

		if(a == b){
			nodes.pb(a);
		}else{
			nodes.pb(a);
			nodes.pb(b);
			edges.pb(mp(a, b));
			edges.pb(mp(b, a));
		}
	}

	// Sort nodes/edges and remove duplicated ones.
	sort(nodes.begin(), nodes.end());
	nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());

	sort(edges.begin(), edges.end());
	edges.erase(unique(edges.begin(), edges.end()), edges.end());

	n = (ui)nodes.size();
	m = (ui)edges.size();

	printf("In original file, n = %d, m = %d \n", n, m);

	fclose(fin);

	// Read-2: read checkins file.
	string name2 = string(dir) + "/" + string(checkins_file);
	fin = open_file(name2.c_str(), "r");

	char buf[MAXBUF];
	char *p;
	const char delim[] = " \t\n";

	ui node;					// current node ID
	float latitude, longitude;	// location
	int preNode = -1;			// previous node ID

	while(!feof(fin)){
		fgets(buf, MAXBUF, fin);
		if(buf[strlen(buf)-1] == '\n'){
			buf[strlen(buf)-1] = '\0';
		}
		int cnt = 0;			// index of fields
		p = strtok(buf, delim);							
		while(p){
			if(cnt == 0) {
				node = atoi(p);
				if(node == preNode) break;		// only acquire the first check-in location for each node
				preNode = node;
				++cnt;
			}else if(cnt == 2) {
				latitude = atof(p);
				++cnt;
			}
			else if(cnt== 3) {
				longitude = atof(p);
				++cnt;
			}
			else cnt++;
			p = strtok(NULL, delim);
		}
		if(cnt > 0) node_loc.pb(mp(node, mp(latitude,longitude)));
	}

	sort(node_loc.begin(), node_loc.end());
	node_loc.erase(unique(node_loc.begin(), node_loc.end()), node_loc.end());

	fclose(fin);

}


void output(const char *dir){
	// map
	vector<ui> nodesVec;
	vector<pair<ui, ui> > edgesVec;
	vector<pair<ui, pair<float, float> > > locationsVec;

	for(ui i = 0; i < node_loc.size(); i++) {
		pairNodes.insert(pair<ui, ui>(node_loc[i].first, i));	// <old, new>
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

	for(vector<pair<ui, pair<float, float> > >::iterator it3 =  node_loc.begin(); it3 != node_loc.end(); ++it3) {
		locationsVec.pb(mp(pairNodes.find((*it3).first)->second, mp((*it3).second.first,(*it3).second.second)));
	}

	sort(nodesVec.begin(), nodesVec.end());
	nodesVec.erase(unique(nodesVec.begin(), nodesVec.end()), nodesVec.end());

	sort(edgesVec.begin(), edgesVec.end());
	edgesVec.erase(unique(edgesVec.begin(), edgesVec.end()), edgesVec.end());

	sort(locationsVec.begin(), locationsVec.end());
	locationsVec.erase(unique(locationsVec.begin(), locationsVec.end()), locationsVec.end());

	n = (ui)nodesVec.size();
	m = (ui)edgesVec.size();

	printf("In new file, n = %d, m = %d \n", n, m);

	// Write to files -- Nodes and Edges.
	string name1 = string(dir) + "/Nodes.csv";
	string name2 = string(dir) + "/Edges.csv";

	ofstream fout1, fout2;
	fout1.open(name1,ios::out);
	fout2.open(name2,ios::out);

	fout1<<"Id"<<','<<"Latitude"<<','<<"Longitude"<<endl;
	fout2<<"Source"<<','<<"Target"<<endl;

	for(vector<pair<ui, pair<float, float> > >::iterator it = locationsVec.begin(); it != locationsVec.end(); ++it){
		fout1<<(*it).first<<','<<(*it).second.first<<','<<(*it).second.second<<endl;
	}

	for(vector<pair<ui, ui> >::iterator it = edgesVec.begin(); it != edgesVec.end(); ++it) {
		fout2<<(*it).first<<','<<(*it).second<<endl;
	}

	fout1.close();
	fout2.close();

}
