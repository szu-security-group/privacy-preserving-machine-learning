#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include "LogisticRegression.h"
#include "Functionalities.h"
using namespace std;

// 初始化权重
void initializeLR(vector<myType> &weights){
    size_t lower = 30;
	size_t higher = 50;
	size_t decimation = 10000;
	size_t size = weights.size();
	vector<myType> temp(size);

	for (size_t i = 0; i < size; ++i)
		temp[i] = 0;
		// temp[i] = floatToMyType(1);
		// temp[i] = floatToMyType((float)(rand() % (higher - lower) + lower)/decimation);
	
	if (partyNum == PARTY_S)
		for (size_t i = 0; i < size; ++i)
			weights[i] = temp[i];
	else if (partyNum == PARTY_A or partyNum == PARTY_D)
		for (size_t i = 0; i < size; ++i)
			weights[i] = temp[i];
	else if (partyNum == PARTY_B or partyNum == PARTY_C)		
		for (size_t i = 0; i < size; ++i)
			weights[i] = 0;
}

void trainLR(vector<myType> trainX, vector<myType> trainY, vector<myType> &weights){
    if(STANDALONE){
        vector<myType> g(MINI_BATCH_SIZE, 0);
		for(int j = 0; j < MINI_BATCH_SIZE; j++){
			for(int k = 0; k < input_size; k++)
                // g(x^i)=w^T * x^i
				g[j] = addModuloOdd(g[j], multiplyMyTypesSA(trainX[j*input_size + k], weights[k], FLOAT_PRECISION));//第j个样本
		}
        // 激活函数：h = subSigmoid(g)
        vector<myType> h(MINI_BATCH_SIZE, 0);
		subSigmoid(g, h);
        // update weigths based on error
        vector<myType> error(input_size, 0);
        for(int j = 0; j < input_size; j++){
            for(int k = 0; k < MINI_BATCH_SIZE; k++){
                error[j] = addModuloOdd(error[j], multiplyMyTypesSA(subtractModuloOdd(trainY[k], h[k]), trainX[k*input_size+j], FLOAT_PRECISION));
            }
            weights[j] = addModuloOdd(weights[j], dividePlainSA(multiplyMyTypesSA(LEARNING_RATE, error[j], FLOAT_PRECISION), MINI_BATCH_SIZE));
        }
    }

	if(MPC){
		vector<myType> h(MINI_BATCH_SIZE, 0);
		vector<myType> g(MINI_BATCH_SIZE, 0);
		for(int j = 0; j < MINI_BATCH_SIZE; j++){
			vector<myType> x(input_size, 0);//样本点
			vector<myType> z(input_size, 0);
			for(int k = 0; k < input_size; k++)
				x[k] = trainX[j*input_size+k];//第j个样本点
			
			funcDotProductMPC(x, weights, z, input_size);//g=x*weights, Z_L

			for(int k = 0; k < input_size; k++)
				g[j] = addModuloOdd(z[k], g[j]);
		}
		subSigmoid(g, h);
		vector<myType> error(input_size, 0);
		subtractVectors(trainY, h, h, MINI_BATCH_SIZE);
		for(int j = 0; j < input_size; j++){
			vector<myType> x(MINI_BATCH_SIZE, 0);
			vector<myType> x_mul_h(MINI_BATCH_SIZE, 0);
			for(int k = 0; k < MINI_BATCH_SIZE; k++)
				x[k] = trainX[k*input_size+j];

			funcDotProductMPC(x, h, x_mul_h, MINI_BATCH_SIZE);

			for(int k = 0; k < MINI_BATCH_SIZE; k++){
				error[j] = addModuloOdd(error[j], x_mul_h[k]);
			}
		}
		if (PRIMARY)
			funcTruncate2PC(error, LOG_MINI_BATCH + LOG_LEARNING_RATE, input_size, PARTY_A, PARTY_B);

		addVectors(weights, error, weights, input_size);
	}
}

void testLR(vector<myType> testX, vector<myType> testY, vector<myType> weights, vector<size_t> &counter){
	vector<myType> predictY(MINI_BATCH_SIZE, 0);
	vector<myType> h(MINI_BATCH_SIZE, 0);

	if(STANDALONE){
		for (int j = 0; j < MINI_BATCH_SIZE; j++){
			for (int k = 0; k < input_size; k++)
				h[j] += multiplyMyTypesSA(testX[j*input_size + k], weights[k], FLOAT_PRECISION);
				// h[j] = addModuloOdd(h[j], multiplyMyTypesSA(testX[j*input_size + k], weights[k], FLOAT_PRECISION));

			if(h[j] < (MINUS_ONE >> 1))//正数
				if(h[j] < floatToMyType(0.5))
					predictY[j] = 0;
				else
					predictY[j] = floatToMyType(1);
			else
				predictY[j] = 0;				
				
			counter[1]++;
			if(predictY[j] == testY[j])
				counter[0]++;
		}
	}

	if(MPC){
		for (int j = 0; j < MINI_BATCH_SIZE; j++){
			vector<myType> x(input_size, 0);//x
			vector<myType> b(input_size, 0);

			for (int k = 0; k < input_size; k++)
				x[k] = testX[j*input_size + k];//第j个测试样本

			//归一化
			// LR_normalize(x);
				
			//2.sum(x*weights)=sum(h)
			funcDotProductMPC(x, weights, b, input_size);

			for(int k = 0;k < input_size;k++)
				h[j] += b[k];
		}
		getPredictY(h, predictY);
		getAccuracy(predictY, testY, counter);
	}
}

void getPredictY(vector<myType> h, vector<myType> &predictY){
	//h>0.5,y=1; h<0.5,y=0
	size_t size = h.size();
	vector<myType> f_share(size);//0.5
	vector<myType> h_sub_f(size);

	if(partyNum == PARTY_C){
		vector<myType> f(size, floatToMyType(0.5));//0.5
		vector<myType> f_1(size);
		vector<myType> f_2(size);
		
		splitIntoShares(f, f_1, f_2, size);//<f>1, <f>2
		sendVector<myType>(f_1, PARTY_A, size);
		sendVector<myType>(f_2, PARTY_B, size);
	}  

	if(PRIMARY){
		receiveVector<myType>(f_share, PARTY_C, size);		
		subtractVectors(h, f_share, h_sub_f, size);
	}
	funcSign3PC(h_sub_f, predictY, size);
}

void getAccuracy(vector<myType> predictY, vector<myType> testY, vector<size_t> &counter){
	vector<myType> temp_y(MINI_BATCH_SIZE), temp_ty(MINI_BATCH_SIZE);

	if(partyNum == PARTY_B)
		sendTwoVectors<myType>(predictY, testY, PARTY_A, MINI_BATCH_SIZE, MINI_BATCH_SIZE);

	if(partyNum == PARTY_A){
		receiveTwoVectors<myType>(temp_y, temp_ty, PARTY_B, MINI_BATCH_SIZE, MINI_BATCH_SIZE);
		addVectors<myType>(temp_y, predictY, temp_y, MINI_BATCH_SIZE);
		addVectors<myType>(temp_ty, testY, temp_ty, MINI_BATCH_SIZE);
	}

	for(size_t i = 0; i < MINI_BATCH_SIZE; i++){
		// cout << temp_y[i] << " " << temp_ty[i] << ";";
		counter[1]++;
		if(temp_y[i] == temp_ty[i])
			counter[0]++;
	}
	// cout << "Rolling accuracy: " << counter[0] << " out of " << counter[1] << " (" << (counter[0]*100/counter[1]) << " %)" << endl;
}
