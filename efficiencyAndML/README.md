# Example of using hipe4ml

https://github.com/hipe4ml/hipe4ml 
documentation: https://hipe4ml.github.io

pip3 install hipe4ml

brew install libomp

pip3 install pylink3

pip3 install flake8

pip3 install uproot

pip3 install xxhash

pip3 install lz4

### check if properly installed

tests/run_tests.sh

### if having errors, then highly likely some packages are missing

pip3 install numpy  pandas  xgboost  matplotlib sklearn

# Example of binary training using XGBoost:

cd Binary/xgboost

python3 binary_jpsi.py
