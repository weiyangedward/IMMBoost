"""

This script find the best parameter through grid search
on a subset of training data. Then a linear svm model 
is trained using all training data and the best 
parameter and predict on test data. Note that liblinear 
library is needed.

"""

import sys
sys.path.append('../../src/liblinear-2.1/python')
from liblinearutil import *
sys.path.append('../../src/libsvm-3.21/tools')
from grid import *

# get training and test data
yTrain, xTrain = svm_read_problem(sys.argv[1])
yTest, xTest = svm_read_problem(sys.argv[2])

# grid search for best parameter on subset of training data
rate, grid = find_parameters(sys.argv[5], '-s 0 -t 0 -v 5')

# train svm
prob = problem(yTrain, xTrain)
param = parameter(" ".join(['-s 2 -B 1 -c', str(grid['c'])]))
m = train(prob, param)

# pred on test data
p_label, p_acc, p_val = predict(yTest, xTest, m)
save_model(sys.argv[3], m)

# output prediction to file
predout = open(sys.argv[4],'w')
predout.write('\n'.join(map(str,p_val)))
predout.close()

