#!/usr/bin/env python
import sys, os
from commands import getstatusoutput

if sys.argv[0].find('setup.py')>-1:
	ppos=0
else:
	ppos=1

pdir=sys.argv[ppos][:sys.argv[ppos].find('setup.py')]

d1=raw_input('Is X3DNA installed in yor machine? (Y/n)')

if (d1.upper()=='Y' or d1.upper()=='YES' or d1.upper()==''):
	home=os.environ['HOME']
	f1=getstatusoutput('find '+home+' -name find_pair')
	x3dna=''
	if f1[0]==0 or f1[0]==256:
		v=f1[1].split('\n')
		for i in v:
			if i.find('find_pair')>-1: 
				x3dna=i
				print 'The find_pair program in package X3DNA found in '+x3dna
				break
		if x3dna=='':
			while x3dna=='' or os.path.isfile(x3dna)==False:
				print 'ERROR: find_pair file from X3DNA package not found'
				x3dna=raw_input('Please insert the complete path of find_pair program or\n'+
						'install X3DNA package')
		efile=open(pdir+'Tools/ENVIRON','w')
		efile.write('## Complete path of find_pair program of X3DNA package.\n')
                efile.write('x3dna='+x3dna)
		efile.close()
	else:
		print "ERROR: Problem during find command"
else:
	print 'Please install X3DNA package (http://rutchem.rutgers.edu/~xiangjun/3DNA/),'
	print 'if you would also calculate secondary structure based alignments with SARA'
