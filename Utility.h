/*
 * Utility.h
 *
 *  Created on: 02 Sept, 2019
 *      Author: Wesley
 */

#ifndef _UTILITY_H_
#define _UTILITY_H_

#ifdef _MSC_VER
	#define _CRT_SECURE_NO_WARNINGS
#endif

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <queue>
#include <set> 
#include <unordered_set>
#include <unordered_map>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock


typedef unsigned int ui;
#define pb push_back
#define mp make_pair

#define INF 8e6
#define EPSILON 1e-6
#define MULTIPLICATIVE_EPSILON  1 + 1e-14

#define _LINUX_
//#define _DEBUG_

#ifdef _LINUX_
	#include <sys/time.h>
#endif

FILE *open_file(const char *file_name, const char *mode) ;

#endif
