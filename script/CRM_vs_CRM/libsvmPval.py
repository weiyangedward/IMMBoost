"""

This script finds the best parameter through grid search
on a subset of training data, and then train a full 
model using all training data along with the best 
parameters. Then the model is use to make prediction on 
test data.

"""
from __future__ import print_function
import sys, os

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# eprint("1..............\n")
# eprint(os.path.realpath(__file__))
# eprint(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), '../src/libsvm-3.21/python'))

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), '../src/libsvm-3.21/python'))
from svmutil import *

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), '../src/libsvm-3.21/tools'))
from grid import *

# import train and test data
yTrain, xTrain = svm_read_problem(sys.argv[1])
yTest, xTest = svm_read_problem(sys.argv[2])

# grid search for best parameter on a subset of training data
rate, grid = find_parameters(sys.argv[5], '-s 0 -t 2 -v 5 -out null')

# train linear SVM
prob = svm_problem(yTrain, xTrain)
param = svm_parameter(" ".join(['-s 0 -t 2 -c', str(grid['c']), '-g', str(grid['g'])]))
m = svm_train(prob, param)

# prediction on test data
p_label, p_acc, p_val = svm_predict(yTest, xTest, m)

# save model
svm_save_model(sys.argv[3], m)

# output confident score
predout = open(sys.argv[4],'w')
predout.write('\n'.join(map(str,p_val)))
predout.close()

