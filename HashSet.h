/*
 * HashSet.h
 *
 *  Created on: 13 Nov, 2019
 *      Author: Wesley
 */

#ifndef _HASH_SET_
#define _HASH_SET_

class HashSet {
private:
	ui *head;
	ui *val;
	ui *next;
	int cur, capacity, P;

public:
	HashSet(int size) ;
	~HashSet() ;

	void insert(int v) ;
	int find(int v) ;
};

#endif