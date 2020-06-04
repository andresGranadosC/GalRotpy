

import sys

import os

# define the name of the directory to be created
path = sys.argv[1]
print(path)

try:
    os.mkdir(path)
except OSError:
    print ("Creation of the directory %s failed" % path)
else:
    print ("Successfully created the directory %s " % path)

#print(sys.argv)

#print ( (str(sys.argv[1]), "is delicious. Would you like to try some?\n") )

#print ("Or would you rather have the", str(sys.argv[2]) , "?\n")

#if (len(sys.argv)>3 ):
#    print ("Usage details: try_this.py <input1> <input2>")
#    sys.exit
