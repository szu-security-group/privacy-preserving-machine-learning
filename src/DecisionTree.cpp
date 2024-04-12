#include <iostream>
#include <math.h>
#include <algorithm>
#include "Functionalities.h"
#include "DecisionTree.h"
using namespace std;


//size_t整数, myType浮点数
DecisionTree DT_train(vector<myType> trainData, vector<myType> trainLabels){
    //[(x1,...),(x2,...),...,(xn,...)]; [y1,y2,...,yn]

    DecisionTree dt = DecisionTree(trainData, trainLabels);
    cout << "having done decision tree training" << endl;
    return dt;
}

//决策树初始化
DecisionTree::DecisionTree(vector<myType> &train_data, vector<myType> &train_label){
    this->train_data = train_data;
    this->train_label = train_label;

    //每个特征的类别/取值, input_size = features_num
    for (size_t i = 0; i < TRAINING_DATA_SIZE; i++)
        for (size_t j = 0; j < input_size; j++)
            featureValues[j].insert(train_data[i*input_size + j]);//每个特征索引下的取值，mytype类型，即p0和p1分别持有一个份额（加性）
    

    vector<size_t> dataset(TRAINING_DATA_SIZE);//样本索引
    vector<size_t> features(input_size);//特征索引
    
    for(size_t i = 0;i < TRAINING_DATA_SIZE;i++)
        dataset[i] = i;
    for (size_t i = 0; i < input_size; i++)
        features[i] = i; 

    DT_root = createDT(dataset, features);//创建决策树   
}

//样本数&特征数->创建决策树
TreeNode* DecisionTree::createDT(vector<size_t> &dataset, vector<size_t> &features){
    TreeNode *root = new TreeNode();
    map<myType, size_t> label_count = labelCount(dataset);//统计样本标签数(label, number)
    // cout << "features: " << features.size() << endl;

    //数据集中样本只有一类标签
    if (STANDALONE){
        if (label_count.size() == 1){
            root->result = label_count.begin()->first;//该类标签
            root->is_leaf = 1;
            return root;
        }
    }
    if(MPC){
        size_t stop_size = 1;
        vector<myType> stop(stop_size, 0);
        //problem：p0方标签{0，1}，p1方标签只有0    solve：在正常情况下可以计算，但是由于这里的P2方持有全零数据样本，为了满足这一特殊情况，计算时只用P0方的统计次数来计算
        if (partyNum == PARTY_A){
            if (label_count.size() == 1){
                stop[0] = 1;
            }
            sendVector<myType>(stop, PARTY_B, stop_size);
            sendVector<myType>(stop, PARTY_C, stop_size);
        }

        if(partyNum == PARTY_B){
            receiveVector<myType>(stop, PARTY_A, stop_size);
        }

        if(partyNum == PARTY_C){
            receiveVector<myType>(stop, PARTY_A, stop_size);
        }
        
        if(stop[0] != 0){//一类标签
            root->result = label_count.begin()->first;//该类标签
            root->is_leaf = 1;
            return root;
        }
    }
    
    // cout << "features: " << features.size() << endl;
    //如果特征集为空，则为单节点树，类别标记为数据集中样本数最多的类（概率）
    if (features.size() == 0){
        root->result = getMaxLabel(label_count);
        root->is_leaf = 1;//叶节点
        return root;
    }


    //其他情况：计算当前数据集每个特征的信息增益
    map<size_t, myType> gains;//信息增益 gains(feature, value)

    for (size_t i = 0; i < features.size(); i++){
        gains[features[i]] = caculateGain(dataset, features[i]);
    }

    //（最大信息增益，特征）
    size_t max_gain_feature = getMaxGainFeature(gains);

    vector<size_t> sub_features = features;
    sub_features.erase(find(sub_features.begin(), sub_features.end(), max_gain_feature));//将最大信息增益的特征从现有属性中去掉，已作为节点

    if(STANDALONE){
        for(myType value : featureValues[max_gain_feature]){
            TreeNode *branch = new TreeNode();//创建分支
            vector<size_t> sub_dataset = splitDataset(dataset, max_gain_feature, value);

            //如果子集为空，则分支节点为叶节点，类别为数据集中出现次数最多的标签
            if(sub_dataset.size() == 0){
                branch->is_leaf = 1;
                branch->result = getMaxLabel(label_count);
                branch->attribute = max_gain_feature;
                branch->attribute_value = value;
                root->branches.push_back(branch);
            }
            //递归创建
            else{
                branch = createDT(sub_dataset, sub_features);
                branch->attribute = max_gain_feature;
                branch->attribute_value = value;
                root->branches.push_back(branch);
            }
        }
    }
    
    if(MPC){
        size_t size = featureValues[max_gain_feature].size();
        if(partyNum == PARTY_A){
            vector<size_t> tran_size(1, featureValues[max_gain_feature].size());
            size = tran_size[0];
            sendVector<size_t>(tran_size, PARTY_B, 1);
        }
        if(partyNum == PARTY_B){
            vector<size_t> trans_size(1);
            receiveVector<size_t>(trans_size, PARTY_A, 1);
            size = trans_size[0];
        }

        vector<myType> values(size, 0);
        size_t j = 0;
        if(partyNum == PARTY_A){
            for(myType value : featureValues[max_gain_feature]){
                values[j] = value;//此情况下(dataset)：p0方持有两类，p1方只有0
                j++;
            }
        }


        for(size_t i = 0;i < size;i++){
            TreeNode *branch = new TreeNode();//创建分支节点
            vector<size_t> sub_dataset = splitDataset(dataset, max_gain_feature, values[i]);

            //如果子集为空，则分支节点为叶节点，类别为数据集中出现次数最多的标签
            if(sub_dataset.size() == 0){
                branch->is_leaf = 1;
                branch->result = getMaxLabel(label_count);
                branch->attribute = max_gain_feature;
                branch->attribute_value = values[i];
                root->branches.push_back(branch);
            }
            //递归创建
            else{
                branch = createDT(sub_dataset, sub_features);
                branch->attribute = max_gain_feature;
                branch->attribute_value = values[i];
                root->branches.push_back(branch);
            }
        }
    }
    
    return root;
}

//统计当前数据集中标签出现的次数，包括三方下统计标签出现的次数
map<myType, size_t> DecisionTree::labelCount(vector<size_t> &dataset){
    map<myType, size_t> result;
    vector<myType> label1(1);
    vector<myType> label2(1);
    vector<myType> label(1);
    vector<myType> beta(1);

    if (STANDALONE){
        for (auto i : dataset)
            result[train_label[i]]++; 
        return result;
    }

    if (MPC){
        // 一般情况下：
        //类别标签1的次数
        label1[0] = train_label[dataset[0]];
        for (auto i : dataset){
            label[0] = train_label[i];
            funcIsEqual(label1, label, beta, 1);//相邻两个标签值是否相等
            
            if (beta[0] == 0)
                result[label1[0]]++;
            else
                label2[0] = label[0];
        }

        //类别标签2的次数
        for (auto i : dataset){
            label[0] = train_label[i];
            funcIsEqual(label2, label, beta, 1);
            
            if(partyNum == PARTY_B){
                if(beta[0] == 0){
                    result[label2[0]]++;
                }
            }else{
                if (beta[0] == 0 && label1[0] != label2[0])
                    result[label2[0]]++;             
            }

        }
    }

    // cout << "label_" << label1[0] << ": " << result[label1[0]] << "; label_" << label2[0] << ": " << result[label2[0]] << endl;

    return result;
}

//数据集中出现次数最多的标签(result, max_time)
myType DecisionTree::getMaxLabel(map<myType, size_t> &label_count){
    size_t max_time = 0;
    myType result;

    for (auto label : label_count){
        if(max_time <= label.second){
            max_time = label.second;
            result = label.first;
        }
    }
    
    return result;
}

//计算信息增益
myType DecisionTree::caculateGain(vector<size_t> &dataset, size_t &feature){
    vector<size_t> sub_dataset;
    myType result = 0;
    myType EntD = caculateEntropy(dataset);

    if (STANDALONE){
        for (myType value : featureValues[feature]){
            sub_dataset = splitDataset(dataset, feature, value);//使用当前属性feature对数据集进行划分，得到多个子集：value1对应dataset1，value2对应dataset2
            result += multiplyMyTypesSA(divideMyTypeSA(floatToMyType(sub_dataset.size()), floatToMyType(dataset.size())), caculateEntropy(sub_dataset), FLOAT_PRECISION);//for v=1 to V：sum(Dv/D * Ent(Dv))
        }
    }

    if (MPC){
        vector<myType> values;
        for(myType value : featureValues[feature]){
            values.push_back(value);//此情况下(dataset)：p0方持有两类，p1方只有0
        }
         
        //保证双方的特征类别数相等
        size_t size = values.size();
        vector<size_t> trans_size(1, size);
        if(partyNum == PARTY_A){
            sendVector<size_t>(trans_size, PARTY_B, 1);
        }
        if(partyNum == PARTY_B){
            receiveVector<size_t>(trans_size, PARTY_A, 1);
            size = trans_size[0];
            if(size > 0)
                values.push_back(0);
        }

        //pi*logpi
        for(size_t i = 0; i < size;i++){
            sub_dataset = splitDataset(dataset, feature, values[i]);
            // cout << values[i] << "; " << sub_dataset.size() << endl;
        
            vector<myType> a(1, 0);
            vector<myType> b(1, 0);
            vector<myType> c(1, 0);//c=a/b
            vector<myType> e(1);//e=c*d
            vector<myType> beta(1);
            if(partyNum == PARTY_A){
                a[0] = floatToMyType(sub_dataset.size());
                b[0] = floatToMyType(dataset.size());
            }

            funcIsEqual(a, c, beta, 1);

            if(beta[0] == 0){
                e[0] = 0;
            }
            else{
                vector<myType> d(1, caculateEntropy(sub_dataset));
                funcDivisionMPC(a, b, c, 1);
                funcDotProductMPC(c, d, e, 1);
            }

            result += e[0];
        }
    }
    
    return  EntD - result;//Gain(D,a)=Ent(D)-Sum(Dv/D * Ent(Dv))
}

//计算信息熵
myType DecisionTree::caculateEntropy(vector<size_t> &dataset){
    map<myType, size_t> label_count = labelCount(dataset);
    size_t len_1 = dataset.size();
    size_t len_2;
    vector<size_t> l(1, label_count.size());//两类
    myType result = 0;

    if(STANDALONE){
        for (auto count : label_count){
            myType pi;
            pi = divideMyTypeSA(floatToMyType(count.second), floatToMyType(len_1));//当前样本集合中，第i类样本所占的比例pi（i=0，1）
            result -= multiplyMyTypesSA(pi, floatToMyType(log2(pi >> FLOAT_PRECISION)), FLOAT_PRECISION);//公式：(-p_0 * log2p_0)+(-p_1 * log2p_1)
        }
    }

    if(MPC){
        if (partyNum == PARTY_A)
            sendVector<size_t>(l, PARTY_B, 1);
    
        if(partyNum == PARTY_B)
            receiveVector<size_t>(l, PARTY_A, 1);
    
        len_2 = l[0];
        
        vector<myType> pi(len_2);//第i类样本的概率
        vector<myType> i_times(len_2, 0);//第i类样本的数量
        vector<myType> len(len_2, 0);//样本总数
        vector<myType> log_pi(len_2);
        vector<myType> pi_log_pi(len_2);
        
        if (partyNum == PARTY_A){
            size_t i = 0;
            for (auto count : label_count){
                i_times[i] = floatToMyType(count.second);
                i = i + 1;
            }
            for (size_t j = 0; j < len_2; j++){
                len[j] = floatToMyType(len_1);
            }
        }

        //pi*logpi
        funcDivisionMPC(i_times, len, pi, len_2);

        /*method1*/
        vector<myType> ui(len_2);
        if (partyNum == PARTY_C){
            vector<myType> u_1(len_2), u_2(len_2);
            splitIntoShares(ui, u_1, u_2, len_2);

            sendVector<myType>(u_1, PARTY_A, len_2);
            sendVector<myType>(u_2, PARTY_B, len_2);
        }
        
        if(PRIMARY){
            receiveVector<myType>(ui, PARTY_C, len_2);//u0+u1=0
        }

        if(partyNum == PARTY_A){//p0: u0
            vector<myType> pi_add_u0(len_2);
            for(size_t j = 0;j < len_2;j++){
                pi_add_u0[j] = pi[j] + ui[j];
                log_pi[j] = ui[j];
            }
            sendVector<myType>(pi_add_u0, PARTY_B, len_2);
        }

        if(partyNum == PARTY_B){//p1: log(x0+u0+u1+x1)+u1
            vector<myType> pi_add_u1(len_2);
            receiveVector<myType>(pi_add_u1, PARTY_A, len_2);
            for(size_t j = 0;j < len_2;j++){
                myType flag = floatToMyType(log2(pi_add_u1[j] + ui[j] + pi[j]) - log2(floatToMyType(1)));
                log_pi[j] = flag + ui[j];
            } 
        }

        funcDotProductMPC(pi, log_pi, pi_log_pi, len_2);

        for (size_t j = 0; j < len_2; j++)
            result = result - pi_log_pi[j];
    }

    return result;
}

//根据特征值划分数据集，feature：第几个特征，value：第几个特征的特征值
vector<size_t> DecisionTree::splitDataset(vector<size_t> &dataset, size_t &feature, myType &value){
    vector<size_t> result;

    if (STANDALONE){
        for (size_t i : dataset){
            if (train_data[i*input_size + feature] == value)
                result.push_back(i); 
        }
    }
    
    if (MPC){
        for (size_t i : dataset){
            size_t size = 1;
            vector<myType> x(size, train_data[i*input_size + feature]);//取出每个样本的某个feature
            vector<myType> y(size, value);
            vector<myType> beta(size);
            funcIsEqual(x, y, beta, size);
            
            if(beta[0] == 0)
                result.push_back(i);
        }
        
        return result;
    }
}

//获取最大信息增益对应的特征
size_t DecisionTree::getMaxGainFeature(map<size_t, myType> &gains){
    size_t max_gain_feature;
    
    //standalone
    if(STANDALONE){
        myType max_gain = 0;  
        for(auto gain : gains){
            if(max_gain <= gain.second){
                max_gain = gain.second;
                max_gain_feature = gain.first;
            }
        }
    }
    
    //MPC
    if(MPC){
        size_t size2 = 1;
        vector<myType> max(size2);
        vector<myType> max_gain(size2), beta(size2, 0), flagGain;
        for(auto gain : gains){
            flagGain.push_back(gain.second);
        }

        funcMax(flagGain, max, flagGain.size(), size2);

        for (auto gain : gains){
            max_gain[0] = gain.second;

            funcIsEqual(max_gain, max, beta, size2);

            if(beta[0] == 0){
                max_gain_feature = gain.first;
            }
        }
    }

    return max_gain_feature;
}

myType DecisionTree::classify(vector<myType> &test_data, TreeNode *root){
    //如果决策树节点是叶节点，直接返回结果
    if(root->is_leaf == 1){
        // cout << root->attribute << endl;
        return root->result;
    }
    
    if (STANDALONE){
        for (auto node : root->branches){
            if(test_data[node->attribute] == node->attribute_value)
                return classify(test_data, node);
        }
    }
    
    if (MPC){
        for (auto node : root->branches){
            //判断两个数是否相等
            size_t size = 1;
            vector<myType> x(size, test_data[node->attribute]);
            vector<myType> y(size, node->attribute_value);
            vector<myType> beta(size);

            funcIsEqual(x, y, beta, size);

            if (beta[0] == 0){
                return classify(test_data, node);                
            }
        }
    }
    
    return 0;
}

void DT_test(vector<myType> testData, vector<myType> testLabels, DecisionTree dt){
    vector<size_t> counter(2, 0);//(预测正确，预测总次数)

    TreeNode *root = dt.DT_root;
    vector<myType> test_data(input_size);
    vector<myType> predict_y(TEST_DATA_SIZE);
    for (size_t i = 0; i < TEST_DATA_SIZE; i++)
    {
        for (size_t j = 0; j < input_size; j++)
            test_data[j] = testData[i * input_size+j];
    
        predict_y[i] = dt.classify(test_data, root);
    }

    if (STANDALONE){
        for (size_t i = 0; i < TEST_DATA_SIZE; i++){
            if (predict_y[i] == testLabels[i])
                counter[0]++;
            counter[1]++;       
        }
    }

    if (MPC){
        vector<myType> temp_y(TEST_DATA_SIZE), temp_ty(TEST_DATA_SIZE);

        if(partyNum == PARTY_B)
		    sendTwoVectors<myType>(predict_y, testLabels, PARTY_A, TEST_DATA_SIZE, TEST_DATA_SIZE);

	    if(partyNum == PARTY_A){
		    receiveTwoVectors<myType>(temp_y, temp_ty, PARTY_B, TEST_DATA_SIZE, TEST_DATA_SIZE);
		    addVectors<myType>(temp_y, predict_y, temp_y, TEST_DATA_SIZE);
		    addVectors<myType>(temp_ty, testLabels, temp_ty, TEST_DATA_SIZE);
	    }

        for (size_t i = 0; i < TEST_DATA_SIZE; i++){
            if (temp_y[i] == temp_ty[i])
                counter[0]++;

            counter[1]++; 
        }
    }
    cout << "accuracy: " << counter[0] << " out of " << counter[1] << " (" << double(counter[0]*100)/counter[1] << " %)" << endl;
}
