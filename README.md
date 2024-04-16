# Privacy-Preserving Machine Learning
This repository lists the code for our work on privacy-preserving logistic regression (PPLR) and privacy-preserving decision tree (PPDT).
We focus on horizontally partitioned data, where medical institutions or data providers hold data samples with the same attributes.
With the help of the secondary server, these model enables two parties or medical institutions holding shares of input samples to train these models without learning any input data or model parameters.
This work is based on the work of [SecureNN](https://github.com/snwagh/securenn-public).

### Privacy-Preserving Logistic Regression
![model](https://github.com/RitaRun/PPML/assets/40885936/dfa87448-ba72-482a-8b95-7e6079dba4ff)
The PPLR model is based on a secure three-party computation of logistic regression using secure multiparty computation.
In this model, the non-linear activation function typically used in logistic regression is replaced with a computationally friendly linear activation function.

### Privacy-Preserving ID3 Decision Tree
PPDT uses a system framework similar to PPLR.
We have implemented multiple secure computing protocols such as logarithmic functions based on this framework.

# Build
Our implementation is based on the [SecureNN](https://github.com/snwagh/securenn-public) code.
We use C++ to further develop this project, and this work should work on Linux or WSL.
Here, we utilize WSL: Ubuntu-20.04.
### Installation
```
sudo apt-get install g++
sudo apt-get install make
sudo apt-get install libssl-dev

git clone [https://github.com/RitaRun/PPML.git](https://github.com/szu-security-group/privacy-preserving-machine-learning.git)
cd ${your_path}$
```
### Execution
See the makefile for details.
```
// plaintext(standalone), PPLR
make SPECTF_LRSA

// 3PC, PPLR
make SPECTF_LR3PC
```
