#include <string>
#include <memory>
#include <vector>
#include "globals.h"
using namespace std;

// 双精度

// 定点数截断
void initializeLR(vector<myType> &weights);
void trainLR(vector<myType> trainX, vector<myType> trainY, vector<myType> &weights);
void testLR(vector<myType> testX, vector<myType> testY, vector<myType> weights, vector<size_t> &counter);
void normalizeLR(vector<myType> &x);
void getPredictY(vector<myType> h, vector<myType> &predictY);
void getAccuracy(vector<myType> predictY, vector<myType> testY, vector<size_t> &counter);