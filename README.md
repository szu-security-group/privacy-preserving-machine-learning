# Privacy-Preserving Logistic Regression
This repository lists the code for our work on privacy-preserving logistic regression (PPLR).
The PPLR model is based on a secure three-party computation of logistic regression using secure multiparty computation.
In this model, the non-linear activation function typically used in logistic regression is replaced with a computationally friendly linear activation function.

We focus on horizontally partitioned data, where medical institutions or data providers hold data samples with the same attributes.
With the help of the secondary server, this model enables two parties or medical institutions holding shares of input samples to train the Logistic Regression model without learning any input data or model parameters.

![model](https://github.com/RitaRun/PPML/assets/40885936/dfa87448-ba72-482a-8b95-7e6079dba4ff)

This work is based on the work of [SecureNN](https://github.com/snwagh/securenn-public).

# Build
Our implementation is based on the [SecureNN](https://github.com/snwagh/securenn-public) code.
We use C++ to further develop this project, and this work should work on Linux or WSL.
Here, we utilize WSL: Ubuntu-20.04.
### Installation
```
sudo apt-get install g++
sudo apt-get install make
sudo apt-get install libssl-dev

git clone https://github.com/RitaRun/PPML.git
cd ${your_path}$
```
### Execution
See the makefile for details.
```
// plaintext(standalone)
make SPECTFLRSA

// 3PC
make SPECTFLR3PC
```
