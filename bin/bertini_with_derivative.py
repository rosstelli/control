from sympy import Symbol, MutableDenseMatrix, Add, Mul, Integer, Matrix
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
#cmpnds = octave.derivative_vector('/home/ross/test.txt')
compounds, l, r = octave.human_parser('/home/ross/test.txt')
#print "First compound is: " + compounds[0]
#print l
#print r
num_of_compounds = l.shape[0]
num_of_reactions = l.shape[1]

result2 = octave.derivative_vector2(l, r)
#print result2.pickle
representation = eval(result2.pickle) #pickle.loads(result2.pickle)
#print result2
print "\n\n\n"

# Check the number of elements



dimensions = representation.shape

#for x in range (0, dimensions[0]):
#	print representation[x]
#print representation[0]
#print dimensions[0]


filecontents = "CONFIG\nTRACKTYPE: 1;\nEND;\nINPUT\nvariable_group x_1"
 # variable_group x, y;\nfunction f;\nf = x^6 + y - 1;\nEND;
# Grab the variables
for x in range (2, num_of_compounds+1):
	filecontents += ", x_" + str(x)

# Grab the reaction constants
filecontents += ";\nconstant k_1"
for x in range (2, num_of_reactions+1):
	filecontents += ", k_" + str(x)
filecontents += ";\n"

# assign values to the reaction constants


for x in range (1, num_of_reactions + 1):
	filecontents += "k_" + str(x) + " = 1;\n"


# Grab the functions from the derivative vector
filecontents += "function f_0"
for x in range(1, dimensions[0]):
	filecontents += ", f_" + str(x)

filecontents += ";\n"

for x in range(0, dimensions[0]):
	filecontents += "f_" + str(x) + " = " + str(representation[x]) + ";\n"

filecontents += "\nEND;"

print filecontents


#equation = Symbol(result.pickle)
#print "'EQUATION'"
#print equation

directory = '/var/tmp/'
filename = 'test3.txt'

# Write to file.
targetfile = open(directory + filename, 'w')
targetfile.truncate()
targetfile.write(filecontents)
targetfile.close()

from subprocess import call
call(["./bertini", directory + filename])
call(["rm", "-f", directory+filename])

print "DONE."

