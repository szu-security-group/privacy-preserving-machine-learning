#include <vector>
#include <map>
#include <set>
#include "globals.h"
using namespace std;

//定义决策树结点
struct TreeNode
{
    /* data */
    size_t is_leaf = 0;
    myType result = 0;
    vector<TreeNode*> branches;
    size_t attribute = 0;
    myType attribute_value = 0;
};

//决策树
class DecisionTree{
   public:
      vector<myType> train_data;
      vector<myType> train_label;
      map<size_t, set<myType>> featureValues;//每个特征值的取值(feature, value)
      TreeNode *DT_root;//决策树根节点

      //dataset数据样本索引，feature特征索引，value特征对应的不同特征取值
      DecisionTree(vector<myType> &train_data, vector<myType> &train_label);
      map<myType, size_t> labelCount(vector<size_t> &dataset);//统计当前数据集(记录数据样本的index)每个标签出现的次数
      vector<size_t> splitDataset(vector<size_t> &dataset, size_t &feature, myType &value);//分割数据集


      //gain
      myType caculateEntropy(vector<size_t> &dataset);//计算信息熵
      myType caculateGain(vector<size_t> &dataset, size_t &feature);//计算信息增益
      size_t getMaxGainFeature(map<size_t, myType> &gains);//获取最大信息增益的特征/属性 gains(feature, value)

      myType getMaxLabel(map<myType, size_t> &label_count);//获取出现次数最多的标签 label_count(feature, times)

      TreeNode* createDT(vector<size_t> &dataset, vector<size_t> &features);//创建决策树

      myType classify(vector<myType> &test_data, TreeNode *root);//测试集
};

DecisionTree DT_train(vector<myType> trainData, vector<myType> trainLabels);

void DT_test(vector<myType> testData, vector<myType> testLabels, DecisionTree dt);