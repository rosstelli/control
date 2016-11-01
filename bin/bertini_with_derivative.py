from sympy import Symbol, cos, MutableDenseMatrix, Add, Mul, Integer, Matrix
from oct2py import octave
import pickle

#import scipy.io as sio
#result = sio.loadmat('test.mat')
#print result
#result = octave.reduced_determinant('test.txt', 1)
#print "REDUCED DETERMINANT"
#print result.flat
#print "'UNICODE'"
#print result.unicode
#print "'Pickle'"
#print result.pickle


print "DERIVATIVE VECTOR"
result2 = octave.derivative_vector('test.txt')
#print result2.pickle
representation = eval(result2.pickle) #pickle.loads(result2.pickle)
print representation

dimensions = representation.shape
print representation[1][0]


filecontents = "CONFIG\nTRACKTYPE: 1;\nEND;\nINPUT\nvariable_group x, y;\nfunction f;\nf = x^6 + y - 1;\nEND;"



#equation = Symbol(result.pickle)
#print "'EQUATION'"
#print equation

directory = '/var/tmp/'
filename = 'test3.txt'

targetfile = open(directory + filename, 'w')
targetfile.truncate()
targetfile.write(filecontents)
targetfile.close()

from subprocess import call
#call(["./bertini", directory + filename])
call(["rm", "-f", directory+filename])

print "DONE."

