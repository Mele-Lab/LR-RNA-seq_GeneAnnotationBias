#usage: python getMultiReads.py split_multiple_times.pkl multiReads.out

import sys
import pickle

pickleFile = open(sys.argv[1],"rb")

data = pickle.load(open(sys.argv[1], "rb"))

print (data, file=open(sys.argv[2],'a'))
