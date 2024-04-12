/*
* Copyright(c) 2019 Sameer Wagh.
* This file is part of the secure-nn project, is was taken from the file Functionalities.h written by Sameer Wagh.
* Many changes and additions have been made and only part of the file written by Sameer Wagh has been copied 
* only for the use of this project.
*/

#pragma once
#include "tools.h"
#include "connect.h"
#include "globals.h"
using namespace std;

extern void start_time();
extern void start_communication();
extern void end_time(string str);
extern void end_communication(string str);

/******************************** preliminaries ********************************/
void aggregateCommunication();
void funcTruncate2PC(vector<myType> &a, size_t power, size_t size, size_t party_1, size_t party_2);
void funcReconstruct2PC(const vector<myType> &a, size_t size, string str);
void parallelPC(smallType* c, size_t start, size_t end, int t, 
				const smallType* share_m, const myType* r, 
				const smallType* beta, const smallType* betaPrime, size_t dim);
void funcPrivateCompareMPC(const vector<smallType> &share_m, const vector<myType> &r, 
							const vector<smallType> &beta, vector<smallType> &betaPrime, 
							size_t size, size_t dim);
void funcXORModuloOdd2PC(vector<smallType> &bit, vector<myType> &shares, vector<myType> &output, size_t size);
void funcShareConvertMPC(vector<myType> &a, size_t size);
void funcRELUPrime3PC(const vector<myType> &a, vector<myType> &b, size_t size);
void funcDivisionMPC(const vector<myType> &a, const vector<myType> &b, vector<myType> &quotient, 
						size_t size);

/******************************** Secure Computing Protocols ********************************/
void funcDotProductMPC(const vector<myType> &a, const vector<myType> &b, vector<myType> &c, size_t size);
void funcComputeMSB3PC(const vector<myType> &a, vector<myType> &b, size_t size);
void funcSign3PC(const vector<myType> &a, vector<myType> &b, size_t size);
void subSigmoid(vector<myType> g, vector<myType> &h);

void funcIsEqual(const vector<myType> &a, const vector<myType> &b, vector<myType> &beta, size_t size);//judge the equality of two numbers
void funcMax(vector<myType> &a, vector<myType> &max, size_t size1, size_t size2);//find the max value

/******************************** Debug ********************************/
void debugDotProd();
void debugComputeMSB();
void debugSign();