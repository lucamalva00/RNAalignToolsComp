import sys
pdir=sys.argv[0][:sys.argv[0].find('rna2seq.py')]
pdbtools=pdir+'../Tools'
sys.path.append(pdbtools)
from NucleicAcidTools import *

def getRNASequence(filename,chain,atoms=['P'],het='Y'):
        if chain=='_': chain=' '
        try:
                opdb=readNucleicAcidChain(filename,chain,rna_at,het)
		seq=opdb.getSequence()
	except:
		seq=''
	return seq


if __name__=="__main__":
        if len(sys.argv)>1:
		filename=sys.argv[1]
                try:
                        ch=sys.argv[2]
                        if ch==' ': ch='_'
                except: 
                        ch='_'
		try:
			atom=sys.argv[3]
			if atom in ['C3\'','C3*','P']:
				atoms=[atom]
			else:
				atoms=['C3\'']
		except:
			atoms=['C3\'']
		try:
			head=sys.argv[4]
		except:
			head=''
		seq=getRNASequence(filename,ch,atoms=['P'],het='Y')
		if seq!='':
			if head!='': seq='>'+head+'\n'+seq
			print seq
		else:
			print 'ERROR: No nucleotides found in PDB file',filename,'chain',ch
		
	else:
		print 'usage: python rna2seq.py pdbfile [chain, atom, header] '
