import sys,os
from commands import getstatusoutput
pdir=sys.argv[0][:sys.argv[0].find('rnass.py')]
parser3dna='python '+pdir+'out3dna.py '
pdbtools=pdir+'../Tools/'
fvar=open(pdbtools+'ENVIRON','r').readlines()
px3dna=''
x3dna=''

for line in fvar:
        if line.find('x3dna')==0:
                try:
                        px3dna=line[6:].split()[0].strip()
                        if os.path.isfile(px3dna)==False:
                                print 'Warning: find_pair program of X3DNA package not found\n'
                        if px3dna!='' :
                                x3dna=px3dna[:px3dna.index('bin/find_pair')]
                                px3dna=px3dna+' -t '
                except:
                        pass


def generate_ss(pdbfile,chain):
	fileoutput=None
        if chain!='' and  chain!=' ':
		output=pdbfile+'.'+chain
	else:
		output=pdbfile
  	npdb=checkModel(pdbfile)
        cpdb=extractChain(npdb,chain)
	if cpdb==None: cpdb=pdbfile
        o1=getstatusoutput(px3dna+cpdb+' '+output+'.ss')
        if o1[0]!=0:
		up=o1[1].find('unknown residue')
		if up>-1:
			print 'Unknown nucleotide:',o1[1][up:].split()[2]
			o1=getstatusoutput(px3dna[:-3]+cpdb+' '+output+'.ss')
		else:
			pass
        if o1[0]==0:
                o2=getstatusoutput(parser3dna+output+'.ss > '+output+'.pairs')
                if o2[0]==0: fileoutput=output+'.pairs'
                try:
			if npdb!=pdbfile: os.remove(npdb)
                        if cpdb!=pdbfile: os.remove(cpdb)
                        os.remove(output+'.ss')
                        os.remove('col_chains.scr')
                        os.remove('bestpairs.pdb')
                        os.remove('bp_order.dat')
                        os.remove('col_helices.scr')
                        os.remove('ref_frames.dat')
                        os.remove('hel_regions.pdb')
                except:
                        pass
        else:
                pass
        return fileoutput


def verify_ss(ssfile,chain,num=1):
        rss=False
        count=0
        #if chain=='_': chain=' '
        try:
                fss=open(ssfile,'r').readlines()
                for line in fss:
                        if line[0]!='#':
                                pair=line[:-1].split('-')
				if chain=='_': chain=pair[0][0]
				if pair[0][0]==chain.upper() and pair[1][0]==chain.upper(): count=count+1
                        else:
                                pass
                if count>num: rss=True
        except:
                pass
        return rss,chain


def checkModel(pdbfile):
	outfile=pdbfile
        pdbnmr=[]
	pdbtxt=open(pdbfile,'r').read()
	if pdbtxt.find('\nENDMDL'+74*' '+'\n')>-1:
        	lines=open(pdbfile,'r').readlines()
        	fm=0
        	for line in lines:
               		if line[:6]=='MODEL ': fm=fm+1
                	if fm==1: pdbnmr.append(line)
                	if line[:6]=='ENDMDL':
                        	pdbnmr.append(line)
                        	break
        	if len(pdbnmr)>3:
                	outfile=pdbfile+'.nmr'
                	open(outfile,'w').writelines(pdbnmr)
        return outfile


def extractChain(pdbfile,chain):
	npdb=[]
	filename=pdbfile+'.'+chain
	lines=open(pdbfile,'r').readlines()
	for line in  lines:
		if (line[:4]=='ATOM' or line[:6]=='HETATM'):
			if chain!='_':
				if line[21]==chain: npdb.append(line)
			else:
				if line[21]==' ': npdb.append(line)
	if npdb!=[]:
		filename=pdbfile+'.'+chain
		open(filename,'w').writelines(npdb)
	else:
		filename=None
	return filename
			


if __name__=="__main__":
	if len(sys.argv)>=2:
		try:
			ch=sys.argv[2]
			if ch==' ': ch='_'
		except:
			ch='_'
		try:
                        num=int(sys.argv[3])
                except:
                        num=1
		ssfile=generate_ss(sys.argv[1],ch)
		ssok,chain=verify_ss(ssfile,ch,num)
		print 'PDB:', sys.argv[1],'chain:',ch,'Secondary struture file saved in ',ssfile,ssok,chain
	else:
		print 'python ss.py pdbfile chain'


