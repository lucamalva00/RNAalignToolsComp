import sys, re

try:
	lines=open(sys.argv[1],'r').readlines()
	print '# Parsing of the 3DNA file '+sys.argv[1]
	print '# FORMAT: <chain>:<residue>-<chain>:<residue>'
	for line in lines:
		mobj=re.search(' #[ ]*[0-9]* [\|x+] ([-A-Z0-9]):[.]*([0-9]+)_:[^:]*:[.]*([0-9]+)_:([-A-Z0-9])',line)
		if mobj:
			#print mobj.groups()
			data=mobj.groups()
			if data[0]=='-':
				ch1=' '
			else:
				ch1=data[0]
			if data[3]=='-':
                                ch2=' '
                        else:
                                ch2=data[3]
			if int(data[1])>int(data[2]):
				print ch1+':'+data[2]+'-'+ch2+':'+data[1]
			else:
				print ch1+':'+data[1]+'-'+ch2+':'+data[2] 
except:
	print 'File '+sys.argv[1]+' not found'	
