

import sys

print(len(sys.argv))

print(sys.argv)

print ( (str(sys.argv[1]), "is delicious. Would you like to try some?\n") )

print ("Or would you rather have the", str(sys.argv[2]) , "?\n")

if (len(sys.argv)>3 ):

    print ("Usage details: try_this.py <input1> <input2>")
    sys.exit
