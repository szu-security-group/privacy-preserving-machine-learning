#include <iostream>
#include <string>
#include "secondary.h"
#include "connect.h"
#include "AESObject.h"
#include "Functionalities.h"
using namespace std;

extern int partyNum;
int NUM_OF_PARTIES;

AESObject* aes_common;
AESObject* aes_indep;
AESObject* aes_a_1;
AESObject* aes_a_2;
AESObject* aes_b_1;
AESObject* aes_b_2;
AESObject* aes_c_1;
ParallelAESObject* aes_parallel;

int main(int argc, char** argv){

/******************************** Preprocess ********************************/
    // 参数预处理
    parseInputs(argc, argv);

    // AES加密——>随机数
    aes_indep = new AESObject(argv[4]);
	aes_common = new AESObject(argv[5]);
	aes_a_1 = new AESObject("files/keys/keyD");
	aes_a_2 = new AESObject("files/keys/keyD");
	aes_b_1 = new AESObject("files/keys/keyD");
	aes_b_2 = new AESObject("files/keys/keyD");
	aes_c_1 = new AESObject("files/keys/keyD");
	aes_parallel = new ParallelAESObject(argv[5]);

    if (!STANDALONE)
	{
		initializeCommunication(argv[3], partyNum);//初始化3方（或多方）通信
		synchronize(2000000);
	}

    if (PARALLEL)
		aes_parallel->precompute();

/******************************** run Models ********************************/
    logisticRegression();
	// decisionTree();

/******************************** debug ********************************/
	// initializeMPC();
	// // 乘法安全计算协议
	// debugDotProd();

/******************************** clean up ********************************/
    delete aes_common;
	delete aes_indep;
	delete aes_a_1;
	delete aes_a_2;
	delete aes_b_1;
	delete aes_b_2;
	delete aes_c_1;
	delete aes_parallel;
	if (partyNum != PARTY_S)
		deleteObjects();

    return 0;
}