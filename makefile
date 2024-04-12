CXX=g++
SRC_CPP_FILES     := $(wildcard src/*.cpp)
OBJ_CPP_FILES     := $(wildcard utils/*.cpp)
OBJ_FILES    	  := $(patsubst src/%.cpp, src/%.o,$(SRC_CPP_FILES))
OBJ_FILES    	  += $(patsubst utils/%.cpp, utils/%.o,$(OBJ_CPP_FILES))
HEADER_FILES       = $(wildcard src/*.h)

FLAGS := -g -O0 -w -std=c++11 -pthread -msse4.1 -maes -msse2 -mpclmul -fpermissive -fpic
# FLAGS := -g -O3 -w -std=c++11 -pthread -msse4.1 -maes -msse2 -mpclmul -fpermissive -fpic
LIBS := -lcrypto -lssl
OBJ_INCLUDES := -I 'lib_eigen/' -I 'utils/Miracl/' -I 'utils/'
BMR_INCLUDES := $($(OBJ_INCLUDES), -L./)


all: BMRPassive.out

BMRPassive.out: $(OBJ_FILES)
	g++ $(FLAGS) -o $@ $(OBJ_FILES) $(BMR_INCLUDES) $(LIBS)
%.o: %.cpp $(HEADER_FILES)
	$(CXX) $(FLAGS) -c $< -o $@ $(OBJ_INCLUDES)




clean:
	rm -rf BMRPassive.out
	rm -rf src/*.o util/*.o

local: BMRPassive.out
	sh local_run


################################################### STANDALONE ###################################################
# algorithm: logistic regression, dataset: SPECTF
SPECTF_LRSA: BMRPassive.out
	./BMRPassive.out STANDALONE 4 files/parties_ip/parties_localhost files/keys/keyA files/keys/keyAB files/dataset/lr/SPECTF/train_data_p0 files/dataset/lr/SPECTF/train_label_p0 files/dataset/lr/SPECTF/test_data_p0 files/dataset/lr/SPECTF/test_label_p0

# algorithm: decision tree, dataset: Lung
LUNG_DTSA: BMRPassive.out
	./BMRPassive.out STANDALONE 4 files/parties_ip/parties_localhost files/keys/keyA files/keys/keyAB files/dataset/dt/lung_cancer/train_data_p0 files/dataset/dt/lung_cancer/train_label_p0 files/dataset/dt/lung_cancer/test_data_p0 files/dataset/dt/lung_cancer/test_label_p0

################################################### 3PC ###################################################
# Run all three parties with xterm terminal for A.
# algorithm: logistic regression, dataset: SPECTF
SPECTF_LR3PC: BMRPassive.out
	./BMRPassive.out 3PC 2 files/parties_ip/parties_localhost files/keys/keyC files/keys/keyCD files/dataset/lr/SPECTF/train_data_p0 files/dataset/lr/SPECTF/train_label_p0 files/dataset/lr/SPECTF/test_data_p0 files/dataset/lr/SPECTF/test_label_p0 >result2.txt &
	./BMRPassive.out 3PC 1 files/parties_ip/parties_localhost files/keys/keyB files/keys/keyAB files/dataset/lr/SPECTF/train_data_p1 files/dataset/lr/SPECTF/train_label_p1 files/dataset/lr/SPECTF/test_data_p1 files/dataset/lr/SPECTF/test_label_p1 >result1.txt &
	./BMRPassive.out 3PC 0 files/parties_ip/parties_localhost files/keys/keyA files/keys/keyAB files/dataset/lr/SPECTF/train_data_p0 files/dataset/lr/SPECTF/train_label_p0 files/dataset/lr/SPECTF/test_data_p0 files/dataset/lr/SPECTF/test_label_p0 >result0.txt &

# algorithm: decision tree, dataset: Lung
LUNG_DT3PC: BMRPassive.out
	./BMRPassive.out 3PC 2 files/parties_ip/parties_localhost files/keys/keyC files/keys/keyCD files/dataset/dt/lung_cancer/train_data_p0 files/dataset/dt/lung_cancer/train_label_p0 files/dataset/dt/lung_cancer/test_data_p0 files/dataset/dt/lung_cancer/test_label_p0 >result2.txt &
	./BMRPassive.out 3PC 1 files/parties_ip/parties_localhost files/keys/keyB files/keys/keyAB files/dataset/dt/lung_cancer/train_data_p1 files/dataset/dt/lung_cancer/train_label_p1 files/dataset/dt/lung_cancer/test_data_p1 files/dataset/dt/lung_cancer/test_label_p1 >result1.txt &
	./BMRPassive.out 3PC 0 files/parties_ip/parties_localhost files/keys/keyA files/keys/keyAB files/dataset/dt/lung_cancer/train_data_p0 files/dataset/dt/lung_cancer/train_label_p0 files/dataset/dt/lung_cancer/test_data_p0 files/dataset/dt/lung_cancer/test_label_p0 >result0.txt &