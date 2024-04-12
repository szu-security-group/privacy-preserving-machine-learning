#include "secondary.h"
#include "basicSockets.h"
#include "tools.h"
#include "LogisticRegression.h"
#include "DecisionTree.h"
#include <iostream>
using namespace std;

//this player number
int partyNum;

//For faster DGK computation
smallType additionModPrime[PRIME_NUMBER][PRIME_NUMBER];
smallType multiplicationModPrime[PRIME_NUMBER][PRIME_NUMBER];

//communication
extern string * addrs;
extern BmrNet ** communicationSenders;
extern BmrNet ** communicationReceivers;

// dataset preprocessing
vector<myType> trainData, testData;
vector<myType> trainLabels, testLabels;
size_t trainDataBatchCounter = 0;
size_t trainLabelsBatchCounter = 0;
size_t testDataBatchCounter = 0;
size_t testLabelsBatchCounter = 0;

// 解析输入参数
void parseInputs(int argc, char* argv[]){
    assert((sizeof(double) == sizeof(myType)) && "sizeof(double) != sizeof(myType)");

	if (argc < 10) 
		print_usage(argv[0]);

	if (strcmp(argv[1], "STANDALONE") == 0)
		NUM_OF_PARTIES = 1;
	else if (strcmp(argv[1], "3PC") == 0)
		NUM_OF_PARTIES = 3;
	else if (strcmp(argv[1], "4PC") == 0)
		NUM_OF_PARTIES = 4;

	partyNum = atoi(argv[2]);
	
	if (partyNum < 0 or partyNum > 4) 
		print_usage(argv[0]);

	loadData(argv[6], argv[7], argv[8], argv[9]);
}

// 初始化三方参数
void initializeMPC()
{
	//populate offline module prime addition and multiplication tables
	for (int i = 0; i < PRIME_NUMBER; ++i)
		for (int j = 0; j < PRIME_NUMBER; ++j)
		{
			additionModPrime[i][j] = (i + j) % PRIME_NUMBER;
			multiplicationModPrime[i][j] = (i * j) % PRIME_NUMBER;//0-66
		}
}

// 加载训练和预测数据
void loadData(char* filename_train_data, char* filename_train_labels, char* filename_test_data, char* filename_test_labels){
	float temp;
	ifstream f(filename_train_data);
	for (int i = 0; i < TRAINING_DATA_SIZE * input_size; ++i)
	{
		f >> temp;
		trainData.push_back(floatToMyType(temp));
	}
	f.close();

	ifstream g(filename_train_labels);
	for (int i = 0; i < TRAINING_DATA_SIZE * output_size; ++i)
	{
		g >> temp;
		trainLabels.push_back(floatToMyType(temp));
	}
	g.close();

	ifstream h(filename_test_data);
	for (int i = 0; i < TEST_DATA_SIZE * input_size; ++i)
	{
		h >> temp;
		testData.push_back(floatToMyType(temp));
	}
	h.close();

	ifstream k(filename_test_labels);
	for (int i = 0; i < TEST_DATA_SIZE * output_size; ++i)
	{
		k >> temp;
		testLabels.push_back(floatToMyType(temp));
	}
	k.close();	
}

// 样本处理，随机批量
void readMiniBatch(string phase, vector<myType> &X, vector<myType> &Y){
	size_t s = trainData.size();//data size
	size_t t = trainLabels.size();//label size
	if (phase == "training")
	{
		for (int i = 0; i < input_size * MINI_BATCH_SIZE; ++i){
			X[i] = trainData[(trainDataBatchCounter + i)%s];
		}
		for (int i = 0; i < output_size * MINI_BATCH_SIZE; ++i){
			Y[i] = trainLabels[(trainLabelsBatchCounter + i)%t];
		}
		trainDataBatchCounter += input_size * MINI_BATCH_SIZE;
		trainLabelsBatchCounter += output_size * MINI_BATCH_SIZE;

		if (trainDataBatchCounter > s)
			trainDataBatchCounter -= s;
		if (trainLabelsBatchCounter > t)
			trainLabelsBatchCounter -= t;
	}

	size_t p = testData.size();//data size
	size_t q = testLabels.size();//label size
	if (phase == "testing")
	{
		for (int i = 0; i < input_size * MINI_BATCH_SIZE; ++i){
			X[i] = testData[(testDataBatchCounter + i)%p];
		}
		for (int i = 0; i < output_size * MINI_BATCH_SIZE; ++i){
			Y[i] = testLabels[(testLabelsBatchCounter + i)%q];
		}
		testDataBatchCounter += input_size * MINI_BATCH_SIZE;
		testLabelsBatchCounter += output_size * MINI_BATCH_SIZE;

		if (testDataBatchCounter > p)
			testDataBatchCounter -= p;
		if (testLabelsBatchCounter > q)
			testLabelsBatchCounter -= q;
	}
}

// 逻辑回归模型
void logisticRegression(){
	log_print("logistic regression");
	if(!STANDALONE){
		initializeMPC();
	}

	vector<myType> trainX(input_size*MINI_BATCH_SIZE);
	vector<myType> trainY(output_size*MINI_BATCH_SIZE);
	vector<myType> testX(input_size*MINI_BATCH_SIZE);
	vector<myType> testY(output_size*MINI_BATCH_SIZE);
	vector<myType> weights(input_size, 0);//initialize weights
	vector<size_t> counter(2, 0);//预测正确，总预测次数
	initializeLR(weights);

	// 训练
	start_m();
	for(int i = 0;i < NUM_ITERATIONS;i++){
		readMiniBatch("training", trainX, trainY);
		trainLR(trainX, trainY, weights);
	}
	end_m("training");

	// 测试
	start_m();
	for(int i = 0;i < TEST_DATA_SIZE;i++){
		readMiniBatch("testing", testX, testY);
		// 预测一条数据
		testLR(testX, testY, weights, counter);
		cout << "batch " << i << "' " << "accuracy:" << counter[0] << " out of " << counter[1] << " (" << (double(counter[0]*100)/counter[1]) << " %)" << endl;
	}
	end_m("testing");
}

// ID3决策树模型
void decisionTree(){
	log_print("decision tree");

	if(!STANDALONE)
		initializeMPC();

	start_m();
	DecisionTree dt = DT_train(trainData, trainLabels);
	end_m("training");

	start_m();
	DT_test(testData, testLabels, dt);
	end_m("testing");
}

void deleteObjects()
{
	//close connection
	for (int i = 0; i < NUM_OF_PARTIES; i++)
	{
		if (i != partyNum)
		{
			delete communicationReceivers[i];
			delete communicationSenders[i];
		}
	}
	delete[] communicationReceivers;
	delete[] communicationSenders;

	delete[] addrs;
}