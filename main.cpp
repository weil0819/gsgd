/*
 * main.cpp
 *
 *  Created on: 02 Sept, 2019
 *      Author: Wesley
 */

#include "Utility.h"
#include "Graph.h"

void usage() {
	printf("Usage: [0]exe [1]alg [2]graph-dir [3]similarity-threshold [4]density-threshold [5]distance-threshold [6]k value [7 optional]output [8 optional] cluster\n");
}

int main(int argc, char *argv[]) {
	if(argc < 6) {
		usage() ;
		return 0;
	}

	printf("**** Geo-Social Group Detection (Release): %s, %s, %s, %s, %s, %s *** ", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
	printf("\n");

#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#else
	int start, end1, end;
	start = clock();
#endif

	Graph *graph = new Graph(argv[2]);
	graph->read_graph();
	printf("\t*** Finished loading graph!\n");

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#else
	end1 = clock();
#endif

	if(strcmp(argv[1], "gsgd") == 0) graph->gsgd(argv[3], atoi(argv[4]), atoi(argv[5]));
	else if(strcmp(argv[1], "gdcd") == 0) graph->dgcd(argv[3], atoi(argv[4]), atoi(argv[5]));
	else if(strcmp(argv[1], "naive") == 0) graph->baseline(argv[3], atoi(argv[4]), atoi(argv[5]));
	else if(strcmp(argv[1], "random") == 0) graph->advance(argv[3], atoi(argv[4]), atoi(argv[5]));
	else if(strcmp(argv[1], "greedy") == 0) graph->greedy_topk(argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6])); 
	else if(strcmp(argv[1], "swap") == 0) graph->swap_topk(argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
	else if(strcmp(argv[1], "topk") == 0) graph->topk(argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
	else if(strcmp(argv[1], "evalVII") == 0) graph->evalVII(argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
	else usage();

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("Total time without IO: %lld microsecond, %lf millisecond, %lf second\n", mtime-mtime1, (float)(mtime-mtime1)/1000, (float)(mtime-mtime1)/1000000);
#endif

	if(argc > 7 && strcmp(argv[7], "output") == 0) graph->output(argv[3], argv[4], argv[5]);

	if(argc > 8 && strcmp(argv[8], "cluster") == 0) graph->cluster_count(argv[3], argv[4], argv[5]);

	return 0;
}

