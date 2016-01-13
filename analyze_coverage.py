#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import argparse
import sys,math
import matplotlib.pyplot as plt

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
		    self.__coverage[0,2] = self.__s2/self.__n
		    first = False
		else:
		    self.__coverage = np.concatenate((self.__coverage,np.reshape(np.array((self.__n,self.__s/self.__n,self.__s2/self.__n)),(1,3))),axis=0)
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
    x = np.arange(len(h))


    cov1h = np.dot(x,h)/np.sum(h)
    cov2h = np.dot(x*x,h)/np.sum(h)
    con_mean   = np.mean(c[:,1])
    con_stddev = np.sqrt(len(c[:,1]) * np.sum(c[:,2]) - np.sum(c[:,1])**2)/np.sqrt(len(c[:,1])**2-len(c[:,1]))
    print "mean coverage (average over all matched sequences): %8.2lf ± %8.2lf"%(cov1h,np.sqrt(cov2h - cov1h**2))
    print "mean coverage (average over contigs)              : %8.2lf ± %8.2lf"%(con_mean,con_stddev)

    fig = plt.figure()
    
    #ph1 = fig.add_subplot(2,2,1)
    #ph1.plot(x,h)
    #ph1.set_yscale('log')
    #ph1.set_xlabel('coverage')
    #ph1.set_ylabel('#(bp with coverage)')
    
    ph2 = fig.add_subplot(2,1,1)
    ph2.plot(x,h)
    ph2.set_yscale('log')
    ph2.set_xscale('log')
    ph2.set_xlabel('coverage')
    ph2.set_ylabel('#(bp with coverage)')
    
    nc1 = fig.add_subplot(2,1,2)
    nc1.plot(c[:,0],c[:,1],'r+')
    nc1.set_yscale('log')
    nc1.set_xscale('log')
    nc1.set_xlabel('contig length')
    nc1.set_ylabel('contig mean coverage')
    
    plt.show()

    #for contig in c:
	#print contig[0],contig[1],contig[2]
    
    #for i in range(len(h)):
	#print >> sys.stderr,i,h[i]
    
if __name__ == "__main__":
    main()

