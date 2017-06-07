#!/usr/bin/env python

# written by Saatvik Agarwal, Stanford University, 2011
# modified by Michael Hiller to output whole axt entries instead of bed windows

# given an axt file, it outputs all entire entries where at least one window passes the %id and entropy filter

import sys
import math
import tempfile
import subprocess
import shlex
import string
import os

#Declare global counter arrays
targetA = []
targetC = []
targetT = []
targetG = []
numMatches = []

#Declare filter parameters globally
minSeqIdent = 0.0
minEntropy = 0.0
lenAlignment = 0.0
queryLen = 0.0

#This will return the number of As, Cs, Ts, Gs or matches
#between index i (inclusive) and j (exclusive). j > i. It takes in an array arr that
#counts the number of As, Cs, Ts, etc. upto and including index i
def getNum(arr, i, j):
	if i == 0:
		return arr[j - 1]
	else:
		return arr[j - 1] - arr[i - 1]


#Returns the seq identity of the alignment between i and j
#j must be > i otherwise divide by 0 error
def getSeqIdent(i, j):
	return 100.0 * getNum(numMatches, i, j) / (j - i)


def getEntropy(i, j):
	a = float(getNum(targetA, i, j))
	t = float(getNum(targetT, i, j))
	c = float(getNum(targetC, i, j))
	g = float(getNum(targetG, i, j))
	
	if a + t + c + g == 0:
		return 0.0

	aPercent = a / (a + t + c + g)
	tPercent = t / (a + t + c + g)
	gPercent = g / (a + t + c + g)
	cPercent = c / (a + t + c + g)

	entropy = 0.0
	if (aPercent != 0):
		entropy -= aPercent * math.log(aPercent)
	if (cPercent != 0):
		entropy -= cPercent * math.log(cPercent)
	if (gPercent != 0):
		entropy -= gPercent * math.log(gPercent)
	if (tPercent != 0):
		entropy -= tPercent * math.log(tPercent)

	return entropy / math.log(2)



#Initialize global counter arrays
def initCounterArrays():
	global targetA
	global targetC
	global targetT
	global targetG
	global numMatches
	
	targetA = [0]*lenAlignment
	targetC = [0]*lenAlignment
	targetT = [0]*lenAlignment
	targetG = [0]*lenAlignment
	numMatches = [0]*lenAlignment

# output the entire axt entry if there is any sub-alignment of this alignment that passes our filters
def checkAxt(info, targetSeq, querySeq):
	start = 0
	#print targetSeq
	while True:
		if start >= lenAlignment:
			break
		if targetSeq[start] == '-':
			start += 1
			continue
		stoppedAt = -1
		for j in xrange(lenAlignment):
			end = start + j + windowSize
			if end > lenAlignment:
				break
			if getSeqIdent(start, end) >= minSeqIdent and getEntropy(start, end) >= minEntropy:
				stoppedAt = end
			else:
				break
		if stoppedAt == -1:
			start += 1
		else:
			fout.write(info + "\n" + targetSeq + "\n" + querySeq + "\n\n")		# axt requires an empty line after the info, targetseq, queryseq
			break

# processes one axt entry
def processHit(info, targetSeq, querySeq):
	
	#Access globals
	global targetA
	global targetC
	global targetT
	global targetG
	global numMatches
	global lenAlignment
	global queryLen
	
	#Get query and target sequence
	splitted = info.split()
	chrom = splitted[1]
	chromStart = int(splitted[2])
	name = splitted[4]
#	print querySeq
#	string = ''
#	for i in xrange(len(querySeq)):
#		if querySeq[i] == targetSeq[i]:
#			string += '|'
#		else:
#			string += ' '
#	print string
#	print targetSeq
	lenAlignment = len(querySeq)

	#Initialize global counter arrays	
	initCounterArrays()
	
	#Count As, Cs, Ts and Gs and the number of matches, mismatches and gaps in the alignment
	for i in xrange(lenAlignment):
		if (i > 0):
			targetA[i] = targetA[i - 1] 
			targetC[i] = targetC[i - 1] 
			targetT[i] = targetT[i - 1] 
			targetG[i] = targetG[i - 1] 
			numMatches[i] = numMatches[i - 1] 

		if (targetSeq[i].lower() == 'a'):
			targetA[i] += 1
		if (targetSeq[i].lower() == 'c'):
			targetC[i] += 1
		if (targetSeq[i].lower() == 't'):
			targetT[i] += 1
		if (targetSeq[i].lower() == 'g'):
			targetG[i] += 1
			
		if (querySeq[i].lower() == targetSeq[i].lower()):
			numMatches[i] += 1
			
	#Check to see if there is any sub-alignment of this alignment that passes our filters
	checkAxt(info, targetSeq, querySeq)

def printUsage():
	print "Usage is filterAxtIdentityEntropy.py input.axt minSeqIdent minEntropy windowSize output.axt"


if len(sys.argv) != 6:
	printUsage()
	sys.exit(1)

#Take in parameters
fin = open(sys.argv[1])
minSeqIdent = float(sys.argv[2])
minEntropy = float(sys.argv[3])
windowSize = int(sys.argv[4])
fout = open(sys.argv[5], 'w')

#Read hits line by line
lines = [line for line in fin.readlines() if not line.startswith('#')]
if len(lines) == 0:
	sys.exit(0)
fin.close()
i = 0
while i < len(lines):
	if lines[i].strip() == "":
		i += 1
		continue
	info = lines[i].strip()
	targetSeq = lines[i + 1].strip()
	querySeq = lines[i + 2].strip()
	processHit(info, targetSeq, querySeq)
	i += 3
fout.close()
