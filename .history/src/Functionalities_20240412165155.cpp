#pragma once
#include "Functionalities.h"
#include <algorithm>    // std::rotate
#include <thread>
using namespace std;

/******************************** preliminaries ********************************/
// Share Truncation, truncate shares of a by power (in place) (power is logarithmic)
void funcTruncate2PC(vector<myType> &a, size_t power, size_t size, size_t party_1, size_t party_2)
{
	assert((partyNum == party_1 or partyNum == party_2) && "Truncate called by spurious parties");

	if (partyNum == party_1)
		for (size_t i = 0; i < size; ++i)
			a[i] = static_cast<uint64_t>(static_cast<int64_t>(a[i]) >> power);

	if (partyNum == party_2)
		for (size_t i = 0; i < size; ++i)
			a[i] = - static_cast<uint64_t>(static_cast<int64_t>(- a[i]) >> power);
}

void funcReconstruct2PC(const vector<myType> &a, size_t size, string str){
	assert((partyNum == PARTY_A or partyNum == PARTY_B) && "Reconstruct called by spurious parties");

	vector<myType> temp(size);
	if (partyNum == PARTY_B)
		sendVector<myType>(a, PARTY_A, size);

	if (partyNum == PARTY_A)
	{
		receiveVector<myType>(temp, PARTY_B, size);
		addVectors<myType>(temp, a, temp, size);
	
		cout << str << ": ";
		for (size_t i = 0; i < size; ++i)
			print_linear(temp[i], DEBUG_PRINT);
		cout << endl;
	}
}

//Thread function for parallel private compare
void parallelPC(smallType* c, size_t start, size_t end, int t, 
				const smallType* share_m, const myType* r, 
				const smallType* beta, const smallType* betaPrime, size_t dim)
{//share_m=x,dim=64
	size_t index3, index2;
	size_t PARTY;

	smallType bit_r, a, tempM;
	myType valueX;
	thread_local int shuffle_counter = 0;
	thread_local int nonZero_counter = 0;

	// for (size_t i = 0; i < 10; ++i){
	// 	for (size_t j = 0; j < dim; ++j)
	// 		cout << (int)share_m[i*dim+j] << " ";
	// 	cout << "--------------------------" << endl;
	// }


	//Check the security of the first if condition
	for (size_t index2 = start; index2 < end; ++index2)
	{//数
		//valueX=t=r+1 or valueX=r
		if (beta[index2] == 1 and r[index2] != MINUS_ONE)
			valueX = r[index2] + 1;
		else
			valueX = r[index2];
		// cout << "valueX:" << valueX << endl;//valueX:6 6 8 9 10 10 12 12 13 15

		if (beta[index2] == 1 and r[index2] == MINUS_ONE)
		{//r为最大值，对任意x都有x<r
			//One share of zero and other shares of 1
			//Then multiply and shuffle
			for (size_t k = 0; k < dim; ++k)
			{//位dim=64
				index3 = index2*dim + k;
				c[index3] = aes_common->randModPrime();
				// cout << c[index3] << " ";
				if (partyNum == PARTY_A)
					c[index3] = subtractModPrime((k!=0), c[index3]);

				c[index3] = multiplyModPrime(c[index3], aes_parallel->randNonZeroModPrime(t, nonZero_counter));
			}
			// cout << endl;
		}
		else
		{//beta=0 || beta=1,r[index2]!=MINUS_ONE
			//Single for loop
			a = 0;
			for (size_t k = 0; k < dim; ++k)
			{//位
				index3 = index2*dim + k;
				c[index3] = a;
				tempM = share_m[index3];//tempM=x[bit]
				bit_r = (smallType)((valueX >> (63-k)) & 1);//对valueX(t or r)从高到低取位（64位）
				// cout << (int)bit_r << " ";

				//w_bit=a=bit_r xor bit_x
				if (bit_r == 0)
					a = addModPrime(a, tempM);//a=tempM
				else
					a = addModPrime(a, subtractModPrime((partyNum == PARTY_A), tempM));//a=1-tempM

				if (!beta[index2])
				{//beta=0
					//c=jr-x+j+(x xor r)=jr-x+j+a=j(1+r)-x+a
					if (partyNum == PARTY_A)
						c[index3] = addModPrime(c[index3], 1+bit_r);//p0--j=1:1+r+a p1--j=0:a
					c[index3] = subtractModPrime(c[index3], tempM);//1+r+a-x or a-x
				}
				else
				{//beta=1
					//c=-jt+x+j+(x xor t)=j(1-t)+x+a
					if (partyNum == PARTY_A)
						c[index3] = addModPrime(c[index3], 1-bit_r);//p0--j=1:1-t+a
					c[index3] = addModPrime(c[index3], tempM);//1-t+a+x
				}
				// cout << (int)c[index3] << " ";
				c[index3] = multiplyModPrime(c[index3], aes_parallel->randNonZeroModPrime(t, nonZero_counter));//di=si*ci
			}
			// cout << endl;
		}
		aes_parallel->AES_random_shuffle(c, index2*dim, (index2+1)*dim, t, shuffle_counter);
	}
	aes_parallel->counterIncrement();
}

// Private Compare functionality
//比较share_m和r的大小
void funcPrivateCompareMPC(const vector<smallType> &share_m, const vector<myType> &r, 
							const vector<smallType> &beta, vector<smallType> &betaPrime, 
							size_t size, size_t dim)
{
	log_print("funcPrivateCompareMPC");

	assert(dim == BIT_SIZE && "Private Compare assert issue");
	size_t sizeLong = size*dim;
	size_t index3, index2;
	size_t PARTY;

	if (THREE_PC)
		PARTY = PARTY_C;//p2
	else if (FOUR_PC)
		PARTY = PARTY_D;


	if (PRIMARY)
	{
		smallType bit_r, a, tempM;
		vector<smallType> c(sizeLong);
		myType valueX;

		if (PARALLEL)
		{
			thread *threads = new thread[NO_CORES];
			int chunksize = size/NO_CORES;//chunksize=2

			for (int i = 0; i < NO_CORES; i++)
			{
				int start = i*chunksize;//start:0,2,4,6
				int end = (i+1)*chunksize;//end:2,4,6,10
				if (i == NO_CORES - 1)
					end = size;
				
				threads[i] = thread(parallelPC, c.data(), start, end, i, share_m.data(), 
									r.data(), beta.data(), betaPrime.data(), dim);
				// threads[i] = thread(parallelPC, ref(c.data()), start, end, i, ref(share_m.data()), 
				// 					ref(r.data()), ref(beta.data()), ref(betaPrime.data()), dim);
			}

			for (int i = 0; i < NO_CORES; i++)
				threads[i].join();

			delete[] threads;
		}
		else
		{
			//Check the security of the first if condition
			for (size_t index2 = 0; index2 < size; ++index2)
			{
				if (beta[index2] == 1 and r[index2] != MINUS_ONE)
					valueX = r[index2] + 1;//beta=1,(x<=r)=(x<r+1) x xor r+1
				else
					valueX = r[index2];//x xor r

				if (beta[index2] == 1 and r[index2] == MINUS_ONE)
				{
					//One share of zero and other shares of 1
					//Then multiply and shuffle
					for (size_t k = 0; k < dim; ++k)
					{
						index3 = index2*dim + k;
						c[index3] = aes_common->randModPrime();
						if (partyNum == PARTY_A)
							c[index3] = subtractModPrime((k!=0), c[index3]);

						c[index3] = multiplyModPrime(c[index3], aes_common->randNonZeroModPrime());

						// cout << c[index3] << " ";
					}
					// cout << endl;
				}
				else
				{
					//Single for loop
					a = 0;
					for (size_t k = 0; k < dim; ++k)
					{
						index3 = index2*dim + k;
						c[index3] = a;
						tempM = share_m[index3];

						bit_r = (smallType)((valueX >> (63-k)) & 1);

						if (bit_r == 0)
							a = addModPrime(a, tempM);
						else
							a = addModPrime(a, subtractModPrime((partyNum == PARTY_A), tempM));

						if (!beta[index2])
						{
							if (partyNum == PARTY_A)
								c[index3] = addModPrime(c[index3], 1+bit_r);
							c[index3] = subtractModPrime(c[index3], tempM);
						}
						else
						{
							if (partyNum == PARTY_A)
								c[index3] = addModPrime(c[index3], 1-bit_r);
							c[index3] = addModPrime(c[index3], tempM);
						}

						c[index3] = multiplyModPrime(c[index3], aes_common->randNonZeroModPrime());
					}
				}
				aes_common->AES_random_shuffle(c, index2*dim, (index2+1)*dim);
			}
		}
		sendVector<smallType>(c, PARTY, sizeLong);//send c to p2
	}

	if (partyNum == PARTY)
	{//3pc---p2
		vector<smallType> c1(sizeLong);
		vector<smallType> c2(sizeLong);

		receiveVector<smallType>(c1, PARTY_A, sizeLong);
		receiveVector<smallType>(c2, PARTY_B, sizeLong);

		for (size_t index2 = 0; index2 < size; ++index2)
		{
			betaPrime[index2] = 0;//betaPrime=0:if beta=0,x<=r;else x>r
			for (int k = 0; k < dim; ++k)
			{
				index3 = index2*dim + k;
				if (addModPrime(c1[index3], c2[index3]) == 0)
				{
					betaPrime[index2] = 1;//c=0,betaPrime=1:if beta=0, x>r,else x<=r 
					break;
				}	
			}
			// cout << (int)betaPrime[index2] << " ";
		}
		// cout << endl;
	}
}

// XOR shares with a public bit into output.
void funcXORModuloOdd2PC(vector<smallType> &bit, vector<myType> &shares, vector<myType> &output, size_t size)
{//(eta'', eta', theta, size)
	if (partyNum == PARTY_A)
	{//eta'+eta''-2*eta''eta'
		for (size_t i = 0; i < size; ++i)
		{//按位
			if (bit[i] == 1)//eta''=1, output=eta'+1-2*eta'=1-eta'
				output[i] = subtractModuloOdd<smallType, myType>(1, shares[i]);
			else//eat''=0, output=eta'
				output[i] = shares[i];
		}
	}

	if (partyNum == PARTY_B)
	{//eta'-2*eta''eta'
		for (size_t i = 0; i < size; ++i)
		{
			if (bit[i] == 1)//eta''=1, output=eta'-2*eta'=-eta'
				output[i] = subtractModuloOdd<smallType, myType>(0, shares[i]);
			else//eta''=0, output=eta'
				output[i] = shares[i];
		}
	}
}

// Convert shares of a in \Z_L to shares in \Z_{L-1} (in place)
// a \neq L-1
void funcShareConvertMPC(vector<myType> &a, size_t size)
{
	log_print("funcShareConvertMPC");

	vector<myType> r(size);
	vector<smallType> etaDP(size);
	vector<smallType> alpha(size);
	vector<smallType> betai(size);
	vector<smallType> bit_shares(size*BIT_SIZE);
	vector<myType> delta_shares(size);
	vector<smallType> etaP(size);
	vector<myType> eta_shares(size);
	vector<myType> theta_shares(size);
	size_t PARTY;

	if (THREE_PC)
		PARTY = PARTY_C;//p2
	else if (FOUR_PC)
		PARTY = PARTY_D;
	

	if (PRIMARY)
	{
		vector<myType> r1(size);
		vector<myType> r2(size);
		vector<myType> a_tilde(size);

		populateRandomVector<myType>(r1, size, "COMMON", "POSITIVE");
		populateRandomVector<myType>(r2, size, "COMMON", "POSITIVE");
		addVectors<myType>(r1, r2, r, size);//r=r1+r2

		if (partyNum == PARTY_A)
			wrapAround(r1, r2, alpha, size);//alpha=wrap(r1,r2,L)表示溢出, r1+r2>=L>MINUS_ONE--->alpha=1
		// cout << "r;r1;r2:" ;
		// for (size_t i=0;i<size;i++)
		// 	cout <<r[i] << ";" << r1[i] << ";" << r2[i] << " ";
		// cout << endl;

		if (partyNum == PARTY_A)
		{
			addVectors<myType>(a, r1, a_tilde, size);//a_tilde=a+r1
			wrapAround(a, r1, betai, size);//betai=wrap(a,r1,L)---betai=(a > MINUS_ONE - r1)
		}

		if (partyNum == PARTY_B)
		{
			addVectors<myType>(a, r2, a_tilde, size);//a_tilde=a+r2
			wrapAround(a, r2, betai, size);//betai=wrap(a,r2,L)---betai=(a > MINUS_ONE - r2)
		}

		// cout << "a_tilde:" ;
		// for (size_t i=0;i<size;i++)
		// 	cout << a_tilde[i] << " ";
		// cout << endl;

		populateBitsVector(etaDP, "COMMON", size);//随机位数:eta''
		sendVector<myType>(a_tilde, PARTY_C, size);//a_tilde送p2
	}


	if (partyNum == PARTY_C)
	{//p2
		vector<myType> x(size);
		vector<smallType> delta(size);
		vector<myType> a_tilde_1(size);	
		vector<myType> a_tilde_2(size);	
		vector<smallType> bit_shares_x_1(size*BIT_SIZE);
		vector<smallType> bit_shares_x_2(size*BIT_SIZE);
		vector<myType> delta_shares_1(size);
		vector<myType> delta_shares_2(size);

		receiveVector<myType>(a_tilde_1, PARTY_A, size);
		receiveVector<myType>(a_tilde_2, PARTY_B, size);
		addVectors<myType>(a_tilde_1, a_tilde_2, x, size);//x=a_tilde_1+a_tilde_2 重构a_tilde

		// cout << "a_tilde,a_tilde_1,a_tilde_2:" ;
		// for(size_t i=0;i<size;i++)
		// 	cout << x[i] << "," <<a_tilde_1[i] << "," << a_tilde_2[i] << " ";
		// cout << endl; 
		wrapAround(a_tilde_1, a_tilde_2, delta, size);//delta=wrap(a_tilde_1,a_tilde_2,L)

		// cout << "delta:" ;
		// for(size_t i=0;i<size;i++)
		// 	cout << (int)delta[i] << " ";
		// cout << endl;

		sharesOfBits(bit_shares_x_1, bit_shares_x_2, x, size, "INDEP");//p2 generates shares {x1}{x2},and sends to p0,p1
		sendVector<smallType>(bit_shares_x_1, PARTY_A, size*BIT_SIZE);
		sendVector<smallType>(bit_shares_x_2, PARTY_B, size*BIT_SIZE);
		
		sharesModuloOdd<smallType>(delta_shares_1, delta_shares_2, delta, size, "INDEP");//p2 generates shares {delta1}{delta2},and sends to p0,p1
		sendVector<myType>(delta_shares_1, PARTY_A, size);
		sendVector<myType>(delta_shares_2, PARTY_B, size);
	}

	if (PRIMARY)
	{//p0 or p1
		receiveVector<smallType>(bit_shares, PARTY_C, size*BIT_SIZE);
		receiveVector<myType>(delta_shares, PARTY_C, size);
	}

	// to fellow the definition of SecureNN, r is set to r - 1
  for (size_t i = 0; i < size; ++i)
    if(r[i] != 0)
      r[i] = r[i] - 1;
	
	//betaPrime=1:beta=0,x>r;beta=1,x<=r
	//batePrime=0:beta=0,x<=r;beta=1,x>r
	funcPrivateCompareMPC(bit_shares, r, etaDP, etaP, size, BIT_SIZE);//etaDP=beta=eta'';etaP=betaPrime=eta'

	if (partyNum == PARTY)
	{//3pc---p2
		vector<myType> eta_shares_1(size);
		vector<myType> eta_shares_2(size);

		sharesModuloOdd<smallType>(eta_shares_1, eta_shares_2, etaP, size, "INDEP");//p2 generates {etaP1}{etaP2}
		sendVector<myType>(eta_shares_1, PARTY_A, size);
		sendVector<myType>(eta_shares_2, PARTY_B, size);
	}

	if (PRIMARY)
	{//p0 or p1
		receiveVector<myType>(eta_shares, PARTY, size);
		funcXORModuloOdd2PC(etaDP, eta_shares, theta_shares, size);//etaDP=eta'', eta=theta=eta' xor eta''
		//theta=beta+(1-j)·(-alpha-1)+delta+eta
		addModuloOdd<myType, smallType>(theta_shares, betai, theta_shares, size);//add beta
		// to fellow the definition of SecureNN, to add delta
    	addModuloOdd<myType, myType>(theta_shares, delta_shares, theta_shares, size);//add delta
		
		if (partyNum == PARTY_A) {
      	//theta=beta-(alpha+1)+delta+eta
      		for (size_t i = 0; i < size; ++i)
        		alpha[i] += 1;
			subtractModuloOdd<myType, smallType>(theta_shares, alpha, theta_shares, size);//-(alpha+1)
		}
		
		subtractModuloOdd<myType, myType>(a, theta_shares, a, size);//output=a-theta
	}
}

// 3PC: PARTY_A, PARTY_B hold shares in a, want shares of RELU' in b.
void funcRELUPrime3PC(const vector<myType> &a, vector<myType> &b, size_t size)
{
	log_print("funcRELUPrime3PC");
	assert(THREE_PC && "funcRELUPrime3PC called in non-3PC mode");

	vector<myType> twoA(size, 0);
	myType j = 0;

	for (size_t i = 0; i < size; ++i)
		twoA[i] = (a[i] << 1);//左移一位，2倍

	// if(PRIMARY){
	// 	funcReconstruct2PC(a, MINI_BATCH_SIZE, "a");
	// 	funcReconstruct2PC(twoA, MINI_BATCH_SIZE, "twoA");
	// }

	funcShareConvertMPC(twoA, size);//L->L-1
	funcComputeMSB3PC(twoA, b, size);//MSB(2a)

	if (partyNum == PARTY_A)
		j = floatToMyType(1);

	if (PRIMARY){
		for (size_t i = 0; i < size; ++i)
			b[i] = j - b[i];//b=j-MSB(2a): 0,MSB(2a)=1,负; 1,MSB(2a)=0,正
		// cout << b[0] << " ";
	}
}

//All parties start with shares of a number in a and b and the quotient is in quotient.
//a is numerator(分子); b is denumerator(分母)
//a/b=quotient
void funcDivisionMPC(const vector<myType> &a, const vector<myType> &b, vector<myType> &quotient, 
							size_t size)
{
	log_print("funcDivisionMPC");
	if (THREE_PC)
	{
		vector<myType> varQ(size, 0);
		vector<myType> varP(size, 0);
		vector<myType> varD(size, 0);
		vector<myType> tempZeros(size, 0);
		vector<myType> varB(size, 0);
		vector<myType> input_1(size, 0), input_2(size, 0);

		for (size_t i = 0; i < size; ++i)
		{
			varP[i] = 0;
			quotient[i] = 0;//商
		}

		//分子a：1010<<1，分母b：101<<1
		for (size_t looper = 1; looper < FLOAT_PRECISION+1; ++looper)
		{
			if (PRIMARY)
			{//p0 or p1
				for (size_t i = 0; i < size; ++i)
					input_1[i] = -b[i];//input_1+b=L

				//z=x-2^i*y-v_i=x+2^i*(-y)-v_i
				funcTruncate2PC(input_1, looper, size, PARTY_A, PARTY_B);//input_1 >> looper(i)
				// cout << input_1[0] << "," << a[0] << " ";
				addVectors<myType>(input_1, a, input_1, size);//a+2^i*(-b)
				// cout << input_1[0] << endl;
				subtractVectors<myType>(input_1, varP, input_1, size);//z=a+2^i*(-b)-varP
			}
			funcRELUPrime3PC(input_1, varB, size);//varB判断正负:varB=0,input_1<0;varB=1,input_1>0

			//Get the required shares of y/2^i and 2^FLOAT_PRECISION/2^i in input_1 and input_2
			for (size_t i = 0; i < size; ++i)
					input_1[i] = b[i];

			if (PRIMARY)
				funcTruncate2PC(input_1, looper, size, PARTY_A, PARTY_B);//b/2^i

			if (partyNum == PARTY_A)
				for (size_t i = 0; i < size; ++i)
					input_2[i] = (1 << FLOAT_PRECISION);

			if (partyNum == PARTY_B)
				for (size_t i = 0; i < size; ++i)
					input_2[i] = 0;

			if (PRIMARY)
				funcTruncate2PC(input_2, looper, size, PARTY_A, PARTY_B);

			// funcSelectShares3PC(input_1, varB, varD, size);
			// funcSelectShares3PC(input_2, varB, varQ, size);

			vector<myType> A_one(size, 0), B_one(size, 0), C_one(size, 0);
			vector<myType> A_two(size, 0), B_two(size, 0), C_two(size, 0);

			if (HELPER)
			{
				vector<myType> A1_one(size, 0), A2_one(size, 0), 
							   B1_one(size, 0), B2_one(size, 0), 
							   C1_one(size, 0), C2_one(size, 0);

				vector<myType> A1_two(size, 0), A2_two(size, 0), 
							   B1_two(size, 0), B2_two(size, 0), 
							   C1_two(size, 0), C2_two(size, 0);

				populateRandomVector<myType>(A1_one, size, "INDEP", "POSITIVE");
				populateRandomVector<myType>(A2_one, size, "INDEP", "POSITIVE");
				populateRandomVector<myType>(B1_one, size, "INDEP", "POSITIVE");
				populateRandomVector<myType>(B2_one, size, "INDEP", "POSITIVE");
				populateRandomVector<myType>(A1_two, size, "INDEP", "POSITIVE");
				populateRandomVector<myType>(A2_two, size, "INDEP", "POSITIVE");
				populateRandomVector<myType>(B1_two, size, "INDEP", "POSITIVE");
				populateRandomVector<myType>(B2_two, size, "INDEP", "POSITIVE");


				addVectors<myType>(A1_one, A2_one, A_one, size);
				addVectors<myType>(B1_one, B2_one, B_one, size);
				addVectors<myType>(A1_two, A2_two, A_two, size);
				addVectors<myType>(B1_two, B2_two, B_two, size);

				for (size_t i = 0; i < size; ++i)
					C_one[i] = A_one[i] * B_one[i];

				for (size_t i = 0; i < size; ++i)
					C_two[i] = A_two[i] * B_two[i];

				splitIntoShares(C_one, C1_one, C2_one, size);
				splitIntoShares(C_two, C1_two, C2_two, size);

				sendSixVectors<myType>(A1_one, B1_one, C1_one, A1_two, B1_two, C1_two, PARTY_A, size, size, size, size, size, size);
				sendSixVectors<myType>(A2_one, B2_one, C2_one, A2_two, B2_two, C2_two, PARTY_B, size, size, size, size, size, size);
				// sendThreeVectors<myType>(A1_one, B1_one, C1_one, PARTY_A, size, size, size);
				// sendThreeVectors<myType>(A2_one, B2_one, C2_one, PARTY_B, size, size, size);
				// sendThreeVectors<myType>(A1_two, B1_two, C1_two, PARTY_A, size, size, size);
				// sendThreeVectors<myType>(A2_two, B2_two, C2_two, PARTY_B, size, size, size);

			}

			if (PRIMARY)
			{
				receiveSixVectors<myType>(A_one, B_one, C_one, A_two, B_two, C_two, PARTY_C, size, size, size, size, size, size);
				// receiveThreeVectors<myType>(A_one, B_one, C_one, PARTY_C, size, size, size);
				// receiveThreeVectors<myType>(A_two, B_two, C_two, PARTY_C, size, size, size);
				
				vector<myType> E_one(size), F_one(size), temp_E_one(size), temp_F_one(size);
				vector<myType> E_two(size), F_two(size), temp_E_two(size), temp_F_two(size);
				myType temp_one, temp_two;

				subtractVectors<myType>(input_1, A_one, E_one, size);
				subtractVectors<myType>(varB, B_one, F_one, size);
				subtractVectors<myType>(input_2, A_two, E_two, size);
				subtractVectors<myType>(varB, B_two, F_two, size);


				thread *threads = new thread[2];

				threads[0] = thread(sendFourVectors<myType>, ref(E_one), ref(F_one), ref(E_two), ref(F_two), adversary(partyNum), size, size, size, size);
				threads[1] = thread(receiveFourVectors<myType>, ref(temp_E_one), ref(temp_F_one), ref(temp_E_two), ref(temp_F_two), adversary(partyNum), size, size, size, size);

				for (int i = 0; i < 2; i++)
					threads[i].join();

				delete[] threads;

				//HEREEEEEEE
				// if (partyNum == PARTY_A)
				// 	sendFourVectors<myType>(E_one, F_one, E_two, F_two, adversary(partyNum), size, size, size, size);
				// else
				// 	receiveFourVectors<myType>(temp_E_one, temp_F_one, temp_E_two, temp_F_two, adversary(partyNum), size, size, size, size);	

				// if (partyNum == PARTY_B)
				// 	sendFourVectors<myType>(E_one, F_one, E_two, F_two, adversary(partyNum), size, size, size, size);
				// else
				// 	receiveFourVectors<myType>(temp_E_one, temp_F_one, temp_E_two, temp_F_two, adversary(partyNum), size, size, size, size);	


				// sendTwoVectors<myType>(E_one, F_one, adversary(partyNum), size, size);
				// receiveTwoVectors<myType>(temp_E_one, temp_F_one, adversary(partyNum), size, size);
				// sendTwoVectors<myType>(E_two, F_two, adversary(partyNum), size, size);
				// receiveTwoVectors<myType>(temp_E_two, temp_F_two, adversary(partyNum), size, size);


				addVectors<myType>(E_one, temp_E_one, E_one, size);
				addVectors<myType>(F_one, temp_F_one, F_one, size);
				addVectors<myType>(E_two, temp_E_two, E_two, size);
				addVectors<myType>(F_two, temp_F_two, F_two, size);

				for (size_t i = 0; i < size; ++i)
				{
					varD[i] = input_1[i] * F_one[i];
					temp_one = E_one[i] * varB[i];
					varD[i] = varD[i] + temp_one;

					if (partyNum == PARTY_A)
					{
						temp_one = E_one[i] * F_one[i];
						varD[i] = varD[i] - temp_one;
					}
				}
				
				for (size_t i = 0; i < size; ++i)
				{
					varQ[i] = input_2[i] * F_two[i];
					temp_two = E_two[i] * varB[i];
					varQ[i] = varQ[i] + temp_two;

					if (partyNum == PARTY_A)
					{
						temp_two = E_two[i] * F_two[i];
						varQ[i] = varQ[i] - temp_two;
					}
				}

				addVectors<myType>(varD, C_one, varD, size);
				funcTruncate2PC(varD, FLOAT_PRECISION, size, PARTY_A, PARTY_B);

				addVectors<myType>(varQ, C_two, varQ, size);
				funcTruncate2PC(varQ, FLOAT_PRECISION, size, PARTY_A, PARTY_B);
			}

			addVectors<myType>(varP, varD, varP, size);
			addVectors<myType>(quotient, varQ, quotient, size);
		}
	}

}



/******************************** Secure Computing Protocols ********************************/
// 三方乘法安全计算协议
void funcDotProductMPC(const vector<myType> &a, const vector<myType> &b, 
						   vector<myType> &c, size_t size) 
{
	log_print("funcDotProductMPC");

	if (THREE_PC)
	{
		vector<myType> A(size, 0), B(size, 0), C(size, 0);

		if (HELPER)
		{
			vector<myType> A1(size, 0), A2(size, 0), 
						   B1(size, 0), B2(size, 0), 
						   C1(size, 0), C2(size, 0);

			populateRandomVector<myType>(A1, size, "a_1", "POSITIVE");
			populateRandomVector<myType>(A2, size, "a_2", "POSITIVE");
			populateRandomVector<myType>(B1, size, "b_1", "POSITIVE");
			populateRandomVector<myType>(B2, size, "b_2", "POSITIVE");
			populateRandomVector<myType>(C1, size, "c_1", "POSITIVE");

			// populateRandomVector<myType>(A1, size, "INDEP", "POSITIVE");
			// populateRandomVector<myType>(A2, size, "INDEP", "POSITIVE");
			// populateRandomVector<myType>(B1, size, "INDEP", "POSITIVE");
			// populateRandomVector<myType>(B2, size, "INDEP", "POSITIVE");

			addVectors<myType>(A1, A2, A, size);//A=A1+A2
			addVectors<myType>(B1, B2, B, size);//B=B1+B2

			for (size_t i = 0; i < size; ++i)
				C[i] = A[i] * B[i];

			// splitIntoShares(C, C1, C2, size);
			subtractVectors<myType>(C, C1, C2, size);//C2=C-C1
			sendVector<myType>(C2, PARTY_B, size);

			// sendThreeVectors<myType>(A1, B1, C1, PARTY_A, size, size, size);
			// sendThreeVectors<myType>(A2, B2, C2, PARTY_B, size, size, size);
		}

		if (PRIMARY)
		{
			if (partyNum == PARTY_A)
			{
				populateRandomVector<myType>(A, size, "a_1", "POSITIVE");//A1
				populateRandomVector<myType>(B, size, "b_1", "POSITIVE");//B1
				populateRandomVector<myType>(C, size, "c_1", "POSITIVE");//C1
			}

			if (partyNum == PARTY_B)
			{
				populateRandomVector<myType>(A, size, "a_2", "POSITIVE");//A2
				populateRandomVector<myType>(B, size, "b_2", "POSITIVE");//B2
				receiveVector<myType>(C, PARTY_C, size);//receive C2
			}			

			// receiveThreeVectors<myType>(A, B, C, PARTY_C, size, size, size);
			vector<myType> E(size), F(size), temp_E(size), temp_F(size);
			myType temp;

			subtractVectors<myType>(a, A, E, size);//E=a-A
			subtractVectors<myType>(b, B, F, size);//F=b-B

			thread *threads = new thread[2];

			threads[0] = thread(sendTwoVectors<myType>, ref(E), ref(F), adversary(partyNum), size, size);
			threads[1] = thread(receiveTwoVectors<myType>, ref(temp_E), ref(temp_F), adversary(partyNum), size, size);
	
			for (int i = 0; i < 2; i++)
				threads[i].join();

			delete[] threads;

			
			addVectors<myType>(E, temp_E, E, size);
			addVectors<myType>(F, temp_F, F, size);

			

			for (size_t i = 0; i < size; ++i)
			{//output=X·F+E·Y-jE·F
				c[i] = a[i] * F[i];
				temp = E[i] * b[i];
				c[i] = c[i] + temp;

				if (partyNum == PARTY_A)
				{
					temp = E[i] * F[i];
					c[i] = c[i] - temp;
				}
			}
			addVectors<myType>(c, C, c, size);//add C
			funcTruncate2PC(c, FLOAT_PRECISION, size, PARTY_A, PARTY_B);
		}
	}
}

//Compute MSB of a and store it in b
//3PC: output is shares of MSB in \Z_L
void funcComputeMSB3PC(const vector<myType> &a, vector<myType> &b, size_t size)
{
	log_print("funcComputeMSB3PC");
	assert(THREE_PC && "funcComputeMSB3PC called in non-3PC mode");
	
	vector<myType> ri(size);
	vector<smallType> bit_shares(size*BIT_SIZE);
	vector<myType> LSB_shares(size);
	vector<smallType> beta(size);
	vector<myType> c(size);	
	vector<smallType> betaP(size);
	vector<smallType> gamma(size);
	vector<myType> theta_shares(size);

	if (partyNum == PARTY_C)
	{//p2 picks x from Z/L-1--r
		vector<myType> r1(size);
		vector<myType> r2(size);
		vector<myType> r(size);
		vector<smallType> bit_shares_r_1(size*BIT_SIZE);
		vector<smallType> bit_shares_r_2(size*BIT_SIZE);
		vector<myType> LSB_shares_1(size);
		vector<myType> LSB_shares_2(size);

		for (size_t i = 0; i < size; ++i)
		{
			r1[i] = aes_indep->randModuloOdd();
			r2[i] = aes_indep->randModuloOdd();
		}

		addModuloOdd<myType, myType>(r1, r2, r, size);
		sharesOfBits(bit_shares_r_1, bit_shares_r_2, r, size, "INDEP");//位<r[i]>0,<r[i]>1
		sharesOfLSB(LSB_shares_1, LSB_shares_2, r, size, "INDEP");//<r[0]>j---最低有效位
		sendVector<smallType>(bit_shares_r_1, PARTY_A, size*BIT_SIZE);
		sendVector<smallType>(bit_shares_r_2, PARTY_B, size*BIT_SIZE);
		sendTwoVectors<myType>(r1, LSB_shares_1, PARTY_A, size, size);
		sendTwoVectors<myType>(r2, LSB_shares_2, PARTY_B, size, size);
	}

	if (PRIMARY)
	{//p0 or p1
		vector<myType> temp(size);
		receiveVector<smallType>(bit_shares, PARTY_C, size*BIT_SIZE);
		receiveTwoVectors<myType>(ri, LSB_shares, PARTY_C, size, size);//receive ri and LSB_shares
		
		// for(size_t i=0;i<size;i++)
		// 	cout << a[i] << endl;

		addModuloOdd<myType, myType>(a, a, c, size);//c=2*a:-9 -7 -5 -3 0 0 2 4 6 8
		addModuloOdd<myType, myType>(c, ri, c, size);//r=2*a+x

		thread *threads = new thread[2];

		threads[0] = thread(sendVector<myType>, ref(c), adversary(partyNum), size);
		threads[1] = thread(receiveVector<myType>, ref(temp), adversary(partyNum), size);

		for (int i = 0; i < 2; i++)
			threads[i].join();

		delete[] threads;

		//HEREEEEEEE
		// if (partyNum == PARTY_A)
		// 	sendVector<myType>(c, adversary(partyNum), size);
		// else
		// 	receiveVector<myType>(temp, adversary(partyNum), size);

		// if (partyNum == PARTY_B)
		// 	sendVector<myType>(c, adversary(partyNum), size);
		// else
		// 	receiveVector<myType>(temp, adversary(partyNum), size);		


		// sendVector<myType>(c, adversary(partyNum), size);
		// receiveVector<myType>(temp, adversary(partyNum), size);

		addModuloOdd<myType, myType>(c, temp, c, size);//重构r
		populateBitsVector(beta, "COMMON", size);
	}

	funcPrivateCompareMPC(bit_shares, c, beta, betaP, size, BIT_SIZE);//betaP=beta'=beta xor (bit_shares>c)

	if (partyNum == PARTY_C)
	{//p2
		vector<myType> theta_shares_1(size);
		vector<myType> theta_shares_2(size);

		sharesOfBitVector(theta_shares_1, theta_shares_2, betaP, size, "INDEP");
		sendVector<myType>(theta_shares_1, PARTY_A, size);//beta'份额
		sendVector<myType>(theta_shares_2, PARTY_B, size);
	}

	vector<myType> prod(size), temp(size);
	if (PRIMARY)
	{//p0 or p1
		// theta_shares is the same as gamma (in older versions);
		// LSB_shares is the same as delta (in older versions);
		receiveVector<myType>(theta_shares, PARTY_C, size);
		
		myType j = 0;
		if (partyNum == PARTY_A)
			j = floatToMyType(1);

		for (size_t i = 0; i < size; ++i)
			theta_shares[i] = (1 - 2*beta[i])*theta_shares[i] + j*beta[i];//gamma=beta'+j*beta-2*beta*beta'

		for (size_t i = 0; i < size; ++i)
			LSB_shares[i] = (1 - 2*(c[i] & 1))*LSB_shares[i] + j*(c[i] & 1);//delta=x[0]+j*r[0]-2*r[0]*x[0]
	}

	funcDotProductMPC(theta_shares, LSB_shares, prod, size);

	if (PRIMARY)
	{
		populateRandomVector<myType>(temp, size, "COMMON", "NEGATIVE");
		for (size_t i = 0; i < size; ++i)
			b[i] = theta_shares[i] + LSB_shares[i] - 2*prod[i] + temp[i];//alpha=gamma+delta-2*theta
	}
}

// 判断正负函数安全计算协议
void funcSign3PC(const vector<myType> &a, vector<myType> &b, size_t size)
{
	log_print("funcSign3PC");
	assert(THREE_PC && "funcSign3PC called in non-3PC mode");

	vector<myType> twoA(size, 0);
	myType j = 0;

	for (size_t i = 0; i < size; ++i)
		twoA[i] = (a[i] << 1);//左移一位，2倍

	funcShareConvertMPC(twoA, size);//L->L-1
	funcComputeMSB3PC(twoA, b, size);//MSB(2a)

	if (partyNum == PARTY_A)
		j = floatToMyType(1);

	if (PRIMARY){
		for (size_t i = 0; i < size; ++i)
			b[i] = j - b[i];//b=j-MSB(2a): 0,MSB(2a)=1,负; 1,MSB(2a)=0,正
	}
}

// 替换后的激活函数安全计算协议
void subSigmoid(vector<myType> g, vector<myType> &h){
    if(STANDALONE){
        for(int i = 0; i < MINI_BATCH_SIZE; i++){
            // positive or negative
            if(g[i] < (MINUS_ONE >> 1)){
                if(g[i] > floatToMyType(6)){
                    h[i] = floatToMyType(1);
                }else{
                    h[i] = multiplyMyTypesSA(divideMyTypeSA(floatToMyType(1), floatToMyType(12)), g[i], FLOAT_PRECISION) + floatToMyType(0.5);
                }
            }else{
                if(g[i] < (MINUS_ONE - floatToMyType(6))){
                    h[i] = 0;
                }else{
                    h[i] = multiplyMyTypesSA(divideMyTypeSA(floatToMyType(1), floatToMyType(12)), g[i], FLOAT_PRECISION) + floatToMyType(0.5);
                }
            }
        }
    }

	if(MPC){
		size_t size = g.size();
		//B=u+6, C=u-6
		vector<myType> A(size, 0), B(size, 0), C(size, 0), b(size, 0), c(size, 0);
		//theta1=(u+6>0) xor (u-6>0)
		vector<myType> theta1(size, 0), theta2(size, 0);

		if(partyNum == PARTY_C){
			vector<myType> A0(size, 0), A1(size, 0), num_6(size, floatToMyType(6));
			splitIntoShares(num_6, A0, A1, size);
			sendVector<myType>(A0, PARTY_A, size);
			sendVector<myType>(A1, PARTY_B, size);
		}
		if(PRIMARY){
			receiveVector<myType>(A, PARTY_C, size);
			addVectors<myType>(g, A, B, size);
			subtractVectors<myType>(g, A, C, size);
		}
		//b=(u+6>0)
		funcSign3PC(B, b, size);
		//c=(u-6>0)
		funcSign3PC(C, c, size);
		//<b·c>_i
		funcDotProductMPC(b, c, theta1, size);

		if(PRIMARY){
			for(size_t i = 0;i < size; i++){
				theta1[i] = b[i] + c[i] - 2 * theta1[i];
			}	
		}
		funcDotProductMPC(B, theta1, theta2, size);
		if(PRIMARY){
			for(size_t i = 0;i < size; i++){
				theta2[i] = floatToMyType(1.0/12) * theta2[i] >> FLOAT_PRECISION;// >> FLOAT_PRECISION
			}
		}
		addVectors<myType>(theta2, c, h, size);
	}
}

// 判断两个数是否相等
void funcIsEqual(const vector<myType> &a, const vector<myType> &b, vector<myType> &beta, size_t size){
	vector<myType> c(size);
	vector<myType> w(size);
	vector<myType> alpha(size);

	funcDotProductMPC(a, b, c, size);//c=a*b

	for(size_t i = 0; i < size; i++)
		// w[i] = a[i] + b[i] - 2 * c[i];
		w[i] = a[i] - b[i];

	if(partyNum == PARTY_A){//u1 = alpha1
		vector<myType> u(size, 0);
		vector<myType> u1(size);
		vector<myType> u2(size);//-u1
		populateRandomVector<myType>(u1, size, "INDEP", "POSITIVE");

		addVectors<myType>(u1, u, alpha, size);

		subtractVectors<myType>(u, u1, u2, size);
		sendVector<myType>(u2, PARTY_B, size);
	}

	if(partyNum == PARTY_B){//u2 - u1 = alpha2
		vector<myType> u(size, 0);
		vector<myType> u1(size);
		vector<myType> u1_opp(size);
		vector<myType> u2(size);
		receiveVector<myType>(u1_opp, PARTY_A, size);
		populateRandomVector<myType>(u2, size, "INDEP", "POSITIVE");

		addVectors<myType>(u2, u1_opp, alpha, size);
		
		subtractVectors<myType>(u, u2, u1, size);
		sendVector<myType>(u1, PARTY_A, size);
	}

	if(partyNum == PARTY_A){//u1 - u2 = alpha1
		vector<myType> u2_opp(size);
		receiveVector<myType>(u2_opp, PARTY_B, size);

		addVectors<myType>(u2_opp, alpha, alpha, size);
	}

	if(PRIMARY){//beta = w - alpha
		subtractVectors<myType>(w, alpha, beta, size);
		sendVector<myType>(beta, PARTY_C, size);
	}

	if(partyNum == PARTY_C){
		vector<myType> beta1(size);
		vector<myType> beta2(size);
		receiveVector<myType>(beta1, PARTY_A, size);
		receiveVector<myType>(beta2, PARTY_B, size);
		addVectors<myType>(beta1, beta2, beta, size);
		sendVector<myType>(beta, PARTY_A, size);
		sendVector<myType>(beta, PARTY_B, size);
	}

	if(PRIMARY){
		receiveVector<myType>(beta, PARTY_C, size);
	}
}

void funcMax(vector<myType> &a, vector<myType> &max, size_t size1, size_t size2){
	vector<myType> diff(size2), beta(size2);
	max[0] = a[0];//initialize

	for (size_t i = 1; i < size1; i++){
		diff[0] = max[0] - a[i];//求差

		funcRELUPrime3PC(diff, beta, size2);//判断正负: beta=1,diff>=0;beta=0,diff<0
		funcDotProductMPC(diff, beta, max, size2);

		max[0] = max[0] + a[i];
	}
}


/******************************** Debug ********************************/
void debugDotProd(){
	size_t size = 10;
	vector<myType> a(size, 0), b(size, 0), c(size);
	vector<myType> temp(size);

	populateRandomVector<myType>(temp, size, "COMMON", "NEGATIVE");
	for (size_t i = 0; i < size; ++i)
	{
		if (partyNum == PARTY_A)
			a[i] = temp[i] + floatToMyType(i);
		else
			a[i] = temp[i];
	}

	populateRandomVector<myType>(temp, size, "COMMON", "NEGATIVE");
	for (size_t i = 0; i < size; ++i)
	{
		if (partyNum == PARTY_A)
			b[i] = temp[i] + floatToMyType(i);
		else
			b[i] = temp[i];
	}

	// if (PRIMARY)
	// 	for (size_t i = 0; i < size; ++i)
	// 		a[i] = aes_indep->get64Bits();


	// if (PRIMARY)
	// 	for (size_t i = 0; i < size; ++i)
	// 		b[i] = aes_indep->get64Bits();

	funcDotProductMPC(a, b, c, size);

	if (PRIMARY)
		funcReconstruct2PC(c, size, "c");
}

void debugComputeMSB()
{
	size_t size = 10;
	vector<myType> a(size, 0);

	if (partyNum == PARTY_A)
		for (size_t i = 0; i < size; ++i)
			a[i] = i - 5;

	if (THREE_PC)
	{
		vector<myType> c(size);
		funcComputeMSB3PC(a, c, size);

		if (PRIMARY)
			funcReconstruct2PC(c, size, "c");
	}
}

void debugSign()
{
	size_t size = 10;
	vector<myType> inputs(size, 0);

	if (partyNum == PARTY_A)
		for (size_t i = 0; i < size; ++i)
			inputs[i] = aes_indep->get8Bits() - aes_indep->get8Bits();

	if (THREE_PC)
	{
		vector<myType> outputs(size, 0);
		funcSign3PC(inputs, outputs, size);
		if (PRIMARY)
		{
			funcReconstruct2PC(inputs, size, "inputs");
			funcReconstruct2PC(outputs, size, "outputs");
		}
	}
}


extern CommunicationObject commObject;
void aggregateCommunication()
{
	vector<myType> vec(4, 0), temp(4, 0);
	vec[0] = commObject.getSent();
	vec[1] = commObject.getRecv();
	vec[2] = commObject.getRoundsSent();
	vec[3] = commObject.getRoundsRecv();

	if (THREE_PC)
	{
		if (partyNum == PARTY_B or partyNum == PARTY_C)
			sendVector<myType>(vec, PARTY_A, 4);

		if (partyNum == PARTY_A)
		{
			receiveVector<myType>(temp, PARTY_B, 4);
			addVectors<myType>(vec, temp, vec, 4);
			receiveVector<myType>(temp, PARTY_C, 4);
			addVectors<myType>(vec, temp, vec, 4);
		}
	}

	if (partyNum == PARTY_A)
	{
		cout << "------------------------------------" << endl;
		cout << "Total communication: " << (float)vec[0]/1000000 << "MB (sent) and " << (float)vec[1]/1000000 << "MB (recv)\n";
		cout << "Total calls: " << vec[2] << " (sends) and " << vec[3] << " (recvs)" << endl;
		cout << "------------------------------------" << endl;
	}
}