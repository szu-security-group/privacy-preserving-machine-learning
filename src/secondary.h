#ifndef SECONDARY_H
#define SECONDARY_H

#pragma once
#include <sstream>
#include <fstream>
#include <iomanip>
#include <thread>
#include "secCompMultiParty.h"
#include "basicSockets.h"
#include "../utils/TedKrovetzAesNiWrapperC.h"
#include "tools.h"
#include "globals.h"

void parseInputs(int argc, char* argv[]);
void loadData(char* filename_train_data, char* filename_train_labels, char* filename_test_data, char* filename_test_labels);
void readMiniBatch(string phase, vector<myType> &X, vector<myType> &Y);

void logisticRegression();

void decisionTree();

void initializeMPC();
void deleteObjects();
#endif


#pragma once