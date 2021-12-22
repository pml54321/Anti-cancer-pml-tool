
import re, os, sys,math
from collections import Counter
import pandas as pd
import pickle
import numpy as np
import pickle
import webbrowser
from threading import Timer
import pep as p
def GDPC(fastas, **kw):
	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
	}

	groupKey = group.keys()
	baseNum = len(groupKey)
	dipeptide = [g1+'.'+g2 for g1 in groupKey for g2 in groupKey]

	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key

	encodings = []
	header = ['#'] + dipeptide
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])

		code = [name]
		myDict = {}
		for t in dipeptide:
			myDict[t] = 0

		sum = 0
		for j in range(len(sequence) - 2 + 1):
			myDict[index[sequence[j]]+'.'+index[sequence[j+1]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]] + 1
			sum = sum +1

		if sum == 0:
			for t in dipeptide:
				code.append(0)
		else:
			for t in dipeptide:
				code.append(myDict[t]/sum)
		encodings.append(code)

	return encodings

def GAAC(fastas, **kw):
	group = {
		'alphatic': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharge': 'KRH',
		'negativecharge': 'DE',
		'uncharge': 'STCPNQ'
	}

	groupKey = group.keys()

	encodings = []
	header = ['#']
	for key in groupKey:
		header.append(key)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		count = Counter(sequence)
		myDict = {}
		for key in groupKey:
			for aa in group[key]:
				myDict[key] = myDict.get(key, 0) + count[aa]

		for key in groupKey:
			code.append(myDict[key]/len(sequence))
		encodings.append(code)

	return encodings

def AAC(fastas, **kw):
#     myOrder = myAAorder[args.order] if args.order != None else 'ACDEFGHIKLMNPQRSTVWY'
#     kw = {'order': myOrder}
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	#AA = 'ARNDCQEGHILKMFPSTWYV'
	encodings = []
	header = ['#']
	for i in AA:
		header.append(i)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		count = Counter(sequence)
		for key in count:
			count[key] = count[key]/len(sequence)
		code = [name]
		for aa in AA:
			code.append(count[aa])
		encodings.append(code)
	return encodings

def GTPC(fastas, **kw):
	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
	}

	groupKey = group.keys()
	baseNum = len(groupKey)
	triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key

	encodings = []
	header = ['#'] + triple
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])

		code = [name]
		myDict = {}
		for t in triple:
			myDict[t] = 0

		sum = 0
		for j in range(len(sequence) - 3 + 1):
			myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] + 1
			sum = sum +1

		if sum == 0:
			for t in triple:
				code.append(0)
		else:
			for t in triple:
				code.append(myDict[t]/sum)
		encodings.append(code)

	return encodings

def DPC(fastas, **kw):
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	encodings = []
	diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
	header = ['#'] + diPeptides
	encodings.append(header)

	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		tmpCode = [0] * 400
		for j in range(len(sequence) - 2 + 1):
			tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] +1
		if sum(tmpCode) != 0:
			tmpCode = [i/sum(tmpCode) for i in tmpCode]
		code = code + tmpCode
		encodings.append(code)
	return encodings

def TPC(fastas, **kw):
	AA ='ACDEFGHIKLMNPQRSTVWY'
	encodings = []
	triPeptides = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
	header = ['#'] + triPeptides
	encodings.append(header)

	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		tmpCode = [0] * 8000
		for j in range(len(sequence) - 3 + 1):
			tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j+1]]*20 + AADict[sequence[j+2]]] = tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j+1]]*20 + AADict[sequence[j+2]]] +1
		if sum(tmpCode) != 0:
			tmpCode = [i/sum(tmpCode) for i in tmpCode]
		code = code + tmpCode
		encodings.append(code)
	return encodings


def read_text(text):
    fileInput = text
    count=1
    print ("Converting to FASTA...")
    lst=[]
    for x in range(len(text)):
        lst1=[]

        lst1.append(str(count))
        lst1.append(text[x])
        lst.append(lst1)
        count = count + 1
    return(lst)
    print ("Done.")


def readFasta(file):
	if os.path.exists(file) == False:
		print('Error: "' + file + '" does not exist.')
		sys.exit(1)

	with open(file) as f:
		records = f.read()

	if re.search('>', records) == None:
		print('The input file seems not in fasta format.')
		sys.exit(1)

	records = records.split('>')[1:]
	myFasta = []
	for fasta in records:
		array = fasta.split('\n')
		name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
		myFasta.append([name, sequence])
	return myFasta



def Count_CTDC(seq1, seq2):
	sum = 0
	for aa in seq1:
		sum = sum + seq2.count(aa)
	return sum


def CTDC(fastas, **kw):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}

	groups = [group1, group2, group3]
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

	encodings = []
	header = ['#']
	for p in property:
		for g in range(1, len(groups) + 1):
			header.append(p + '.G' + str(g))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		for p in property:
			c1 = Count_CTDC(group1[p], sequence) / len(sequence)
			c2 = Count_CTDC(group2[p], sequence) / len(sequence)
			c3 = 1 - c1 - c2
			code = code + [c1, c2, c3]
		encodings.append(code)
	return encodings


def Count(aaSet, sequence):
	number = 0
	for aa in sequence:
		if aa in aaSet:
			number = number + 1
	cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
	cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

	code = []
	for cutoff in cutoffNums:
		myCount = 0
		for i in range(len(sequence)):
			if sequence[i] in aaSet:
				myCount += 1
				if myCount == cutoff:
					code.append((i + 1) / len(sequence) * 100)
					break
		if myCount == 0:
			code.append(0)
	return code


def CTDD(fastas, **kw):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}

	groups = [group1, group2, group3]
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

	encodings = []
	header = ['#']
	for p in property:
		for g in ('1', '2', '3'):
			for d in ['0', '25', '50', '75', '100']:
				header.append(p + '.' + g + '.residue' + d)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		for p in property:
			code = code + Count(group1[p], sequence) + Count(group2[p], sequence) + Count(group3[p], sequence)
		encodings.append(code)
	return encodings

def CTDT(fastas, **kw):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}

	groups = [group1, group2, group3]
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

	encodings = []
	header = ['#']
	for p in property:
		for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
			header.append(p + '.' + tr)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
		for p in property:
			c1221, c1331, c2332 = 0, 0, 0
			for pair in aaPair:
				if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
					c1221 = c1221 + 1
					continue
				if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
					c1331 = c1331 + 1
					continue
				if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
					c2332 = c2332 + 1
			code = code + [c1221/len(aaPair), c1331/len(aaPair), c2332/len(aaPair)]
		encodings.append(code)
	return encodings


