#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
#import matplotlib

class coverageclass():
    
    def __init__(self,fname):
	try:
	    fp = open(fname)
	except:
	    raise IOError
	self.__coverage = np.zeros((1,3))
	self.__contignames = []
	self.__currenthistolength = 1000
	self.__histo = np.zeros(self.__currenthistolength)
	self.__s = 0
	self.__n = 0
	self.__s2 = 0
	self.__contigname = ""
	first = True
	for line in fp.readlines():
	    if line[0] == "#":
		self.__contigname = line.split()[1]
	    elif line.strip() == "":
		if first:
		    self.__coverage[0,0] = self.__n
		    self.__coverage[0,1] = self.__s/self.__n
		    self.__coverage[0,2] = np.sqrt(self.__n*self.__s2-self.__s*self.__s)/np.sqrt(self.__n*self.__n-self.__n)
		    first = False
		else:
		    self.__coverage = np.concatenate((self.__coverage,np.reshape(np.array((self.__n,self.__s/self.__n,np.sqrt(self.__n*self.__s2-self.__s*self.__s)/np.sqrt(self.__n*self.__n-self.__n))),(1,3))),axis=0)
		self.__contignames.append(self.__contigname)
		self.__n = 0
		self.__s = 0
		self.__s2 = 0
	    else:
		v = line.split()
		self.__n  += 1
		self.__s  += float(v[1])
		self.__s2 += float(v[1])**2
		self.add_histo(int(v[1]))
	fp.close()
	
		
    def add_histo(self,n):
	if n >= self.__currenthistolength:
	    self.__histo = np.concatenate((self.__histo,np.zeros(n - self.__currenthistolength + 1)))
	    self.__currenthistolength = len(self.__histo)
	self.__histo[n] += 1
	    
    def __iter__(self):
	for i in range(len(self.__coverage)):
	    yield self.__contigname[i],self.__coverage[i]


    def get_coverage(self):
	return self.__coverage
    def get_contignames(self):
	return self.__contignames
    def get_histo(self):
	return self.__histo

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--coveragefile")
    args = parser.parse_args()

    cov = coverageclass(args.coveragefile)

    c = cov.get_coverage()
    n = cov.get_contignames()
    h = cov.get_histo()

    for contig in c:
	print contig[0],contig[1],contig[2]
    
    for i in range(len(h)):
	print >> sys.stderr,i,h[i]
    
if __name__ == "__main__":
    main()

