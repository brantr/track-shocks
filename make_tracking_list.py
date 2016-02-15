#!/usr/bin/env python
import sys

fname = sys.argv[1]

fp = open(fname,"r")
fl = fp.readlines()
fp.close()

for i in range(len(fl)):
	l = fl[i].strip('abcdefghijklmnopqrstuvwxyz./\n')
	l = l.lstrip('0')
	print(l)