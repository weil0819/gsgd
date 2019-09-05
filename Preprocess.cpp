/*
 * Preprocess.cpp
 *
 *  Created on: 05 Sept, 2019
 *      Author: Wesley
 */

#include "Utility.h"
#include <ext/hash_map>

#define MAXBUF 8000000
#define hash_map __gnu_cxx::hash_map

using namespace std;

ui n, m;

vector<ui> nodes;
vector<pair<ui, ui> > edges;
vector<pair<ui, pair<double, double> > > node_loc;
hash_map<ui, ui> pairNodes;

void input(const char *dir, const char *edges_file, const char *checkins_file);
void output(const char *dir);

int main(int argc, char *argv[]){

	if(argc < 4){
		printf("Usage: ./Preprocess dir edges_file checkins_file \n");
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

	printf("n = %d, m = %d \n", n, m);

	// re-order node from 0 to n-1.
	for(ui i = 0; i < n; i++){
		pairNodes.insert(pair<ui, ui>(nodes[i], i));	// <old, new>
		// printf("(%d, %d) \n", nodes[i], i);
	}

	fclose(fin);

	// Read-2: read checkins file.
	string name2 = string(dir) + "/" + string(file_name2);
	fin = open_file(name2.c_str(), "r");

	// TODO ...
	
	char buf[MAXBUF];
	char *p;

	while(!feof(fin)){
		fgets(buf, MAXBUF, fin);
		if(buf[strlen(buf)-1] == '\n'){
			buf[strlen(buf)-1] = '\0';
		}
		p = strtok(buf, " ");
		int cnt = 0;
		ui tmpVex;
		set<ui> tmpSet;
		while(p){
			ui endP = (ui)atoi(p);
			if(cnt == 0){
				tmpVex = endP;
			}else{
				if(endP == 1){
					tmpSet.insert(cnt);
				}
			}
			++cnt;
			p = strtok(NULL, " ");
		}

		if(find(nodes.begin(), nodes.end(), tmpVex) != nodes.end()){
			nodeAttr.pb(mp(tmpVex, tmpSet));
		}

	}

	nodeAttr.pop_back();

	fclose(fin);

	string name3 = string(dir) + "/" + string(file_name3);

	fin = open_file(name3.c_str(), "r");

	fgets(buf, MAXBUF, fin);
	if(buf[strlen(buf)-1] == '\n'){
		buf[strlen(buf)-1] = '\0';
	}
	p = strtok(buf, " ");
	int cnt = 0;
	ui tmpVex;
	set<ui> tmpSet;
	while(p){
		++cnt;
		ui endP = (ui)atoi(p);
		if(endP == 1){
			tmpSet.insert(cnt);
		}
		p = strtok(NULL, " ");
	}

	fclose(fin);

	nodeAttr.pb(mp(ID, tmpSet));

}


void output(const char *dir){
	// map
	vector<ui> nodesVec;
	vector<pair<ui, ui> > edgesVec;
	vector<pair<ui, set<ui> > > nodeAttrVec;

	for(vector<ui>::iterator it1 = nodes.begin(); it1 != nodes.end(); ++it1){
		nodesVec.pb(pairNodes.find(*it1)->second);
	}

	for(vector<pair<ui, ui> >::iterator it2 = edges.begin(); it2 != edges.end(); ++it2){
		edgesVec.pb(mp(pairNodes.find((*it2).first)->second, pairNodes.find((*it2).second)->second));
	}

	for(vector<pair<ui, set<ui> > >::iterator it3 = nodeAttr.begin(); it3 != nodeAttr.end(); ++it3){
		nodeAttrVec.pb(mp(pairNodes.find((*it3).first)->second, (*it3).second));
	}

	sort(nodesVec.begin(), nodesVec.end());
	nodesVec.erase(unique(nodesVec.begin(), nodesVec.end()), nodesVec.end());
	sort(edgesVec.begin(), edgesVec.end());
	edgesVec.erase(unique(edgesVec.begin(), edgesVec.end()), edgesVec.end());
	sort(nodeAttrVec.begin(), nodeAttrVec.end());
	nodeAttrVec.erase(unique(nodeAttrVec.begin(), nodeAttrVec.end()), nodeAttrVec.end());


	string name1 = string(dir) + "/b_degree.bin";
	string name2 = string(dir) + "/b_adj.bin";

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

	string name3 = string(dir) + "/attr_adj.txt";

	FILE *fout3 = open_file(name3.c_str(), "w");

	for(vector<pair<ui, set<ui> > >::iterator it1 = nodeAttrVec.begin(); it1 != nodeAttrVec.end(); ++it1){
		fprintf(fout3, "%d ", (*it1).first);
		for(set<ui>::iterator it2 = (*it1).second.begin(); it2 != (*it1).second.end(); ++it2){
			fprintf(fout3, "%d ", *it2);
		}
		fputs("\n", fout3);
	}

	fclose(fout3);

}
