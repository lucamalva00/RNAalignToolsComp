#!/usr/bin/env python
##  Copyright (C) 2008  Emidio Capriotti & Marc A. Marti-Renom
## 
##  This program and all program in this package are free software;
##  you can redistribute it and/or modify it under the terms of the
##  GNU General Public License as published by the Free Software 
##  Foundation; either version 2 of the License, or (at your option)
##  any later version.
## 
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
## 
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA
## 
##  contacts: Emidio Capriotti o Marc A. Marti-Renom
##            Structural Genomics Unit
##            Department of Bioinformatics
##            Centro de Investigacion Pricipe Felipe (CIPF)
##            Avenida Autopista del Saler 16
##            46012 Valencia
##            Spain
##            e-mail: {ecapriotti}{mmarti}@cipf.es, emidio@biocomp.unibo.it
##
## 
##  Last Update November 2008   
##
##
##  The SARA method developed by Capriotti and Marti-Renom aligns two RNA structures
##  using a unit-vector approach introduced in Kedem et al. 1999. The program is 
##  ispired by the MAMMOTH algorithm for protein structural alignment described in 
##  Ortiz et al. 20002. The optimization procedure is based on MaxSub algorithm 
##  developed by Siew et al. 2000. The evaluation of the returned RNA structural 
##  alignment is based on ProSup algorithm implemented by Lackner et al. 2000.
##  
##  REFERENCES
##
##  Capriotti E, Marti-Renom M. (2008)
##  RNA structure alignment by a unit-vector approach. 
##  Bioinformatics, 24; i112-i116. 
##  
##  Kedem K, Chew LP, Elber R. (1999)
##  Unit-vector RMS (URMS) as a tool to analyze molecular dynamics trajectories.
##  Proteins. 37 :554-564.
## 
##  Ortiz AR, Strauss CE, Olmea O. (2002)
##  MAMMOTH (matching molecular models obtained from theory): an automated method for
##  model comparison.
##  Protein Sci. 11 :2606-2621.
##
##  Siew N, Elofsson A, Rychlewski L, Fischer D. (2000)
##  MaxSub: an automated measure for the assessment of protein structure prediction
##  quality.
##  Bioinformatics. 2000 16 :776-785.
##
##  Lackner P, Koppensteiner WA, Sippl MJ, Domingues FS.
##  ProSup: a refined tool for protein structure alignment.
##  Protein Eng. 2000 13: 745-752.
##


import sys, os
from commands import getstatusoutput
from optparse import OptionParser
if sys.argv[0].find('runsara.py')>-1:
        ppos=0
else:
        ppos=1
pdir=sys.argv[ppos][:sys.argv[ppos].find('runsara.py')]

pdbtools=pdir+'Tools/'
sys.path.append(pdbtools)
sara=pdir+'sara.py '
parser3dna='python '+pdir+'Utils/out3dna.py '
fvar=open(pdbtools+'ENVIRON','r').readlines()
px3dna=''
x3dna=''
#SERVER OPTION
#os.environ['X3DNA'] = '/usr/local/apache/htdocs/services/SARA/program/Utils/X3DNA'

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
#px3dna='/home/emidio/X3DNA/bin/find_pair -t '
from NucleicAcidTools import *


def generate_ss(pdbfile,chain):
	fileoutput=None
        if chain!='' and  chain!=' ':
		output=pdbfile+'.'+chain
	else:
		output=pdbfile
  	npdb=checkModel(pdbfile)
        cpdb=extractChain(npdb,chain)
	if cpdb==None: cpdb=pdbfile
        o1=getstatusoutput('export X3DNA='+x3dna+';'+px3dna+cpdb+' '+output+'.ss')
	if o1[0]!=0: o1=getstatusoutput(px3dna[:-3]+cpdb+' '+output+'.ss')
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
			

def checkpdb(pdbfile,chain,atoms):
	ratoms=None
	altatoms=[]
	for atom in atoms:
                if atom[-1]=='*': altatoms.append(atom[:-1]+"'")
                if atom[-1]=="'": altatoms.append(atom[:-1]+"*")
	pdb=readPDBNucleicAcid(pdbfile,chain,atoms)
	tc=pdb.getNucleicAcidChain(0)
	if tc==None:
		pdb=readPDBNucleicAcid(pdbfile,chain,altatoms)
		tc=pdb.getNucleicAcidChain(0)
		if tc!=None: ratoms=altatoms
	else:
		ratoms=atoms
	return ratoms


def getatom(atoms):
	atom=atoms[0]
	if atom[-1]=="'": atom=atom[:-1]+"\\\'"
	return atom

def verify_ss(ssfile,chain,num=3):
	rss=False
	count=0
	if chain=='_': chain=' '
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
	return rss 

def build_ssfile(ssfile1,chain1,ssfile2,chain2,outfile):
	ofile=None
	try:
		out=open(outfile,'w')
		fss1=open(ssfile1,'r').readlines()
		fss2=open(ssfile2,'r').readlines()
		for line in fss1:
			vec=line.split('-')
			if vec[0][0]==chain1.upper() and vec[1][0]==chain1.upper():
				out.write('A'+vec[0][1:]+'-'+'A'+vec[1][1:]+'\n')
		for line in fss2:
			vec=line.split('-')
			if vec[0][0]==chain2.upper() and vec[1][0]==chain2.upper():
				out.write('B'+vec[0][1:]+'-'+'B'+vec[1][1:]+'\n')
		out.close()
		ofile=outfile
	except:
		'Warning: Secondary Structure File not created'
	return ofile



if __name__=="__main__":
        parser = OptionParser('\nUsage : python %prog  pdbfile1 chain1 pdbfile2 chain2 \n[-a atom -o output -p pdb_output_file -g gapopen -e gapext -b SSTrue -c cutoff -l unit-vector -s EvaluationTrue]')
        parser.add_option ("-g", "--gapopen", action="store",type="string",dest="go", help="Gap_Opening")
        parser.add_option ("-e", "--gapext",action="store",type="string",dest="ge", help="Gap_Extension")
	parser.add_option ("-a", "--atom",action="store",type="string",dest="atom", help="Atom Type")
        parser.add_option ("-o", "--output",action="store",type="string",dest="alout", help="Alignment Output File")
        parser.add_option ("-p", "--pdboutput",action="store",type="string",dest="pdbout", help="Structural Output File")
	parser.add_option ("-n", "--limitation",action="store",type="string",dest="slim", help="Length limit of Structure ")
	parser.add_option ("-l", "--unit-vector",action="store",type="string",dest="lv", help="Unit-Vector Elements")
	parser.add_option ("-b", "--non-local",action="store_true",dest="biuval", help="Non Local Unit Vector Alignment")
	parser.add_option ("-c", "--cutoff",action="store",type="string",dest="cutoff", help="Distance thereshold for aligned atoms")
	parser.add_option ("-s", "--evaluate",action="store_true",dest="beva", help="Evaluate Structural Alignment")
	parser.add_option ("--het",action="store_true",dest="het", help="Read HETATM")
	(options, args) = parser.parse_args()
	sopt=''
	if not options.go:
                go=7.0
        else:
		try:
                	go=float(options.go)
		except:
			go=7.0
		if go<=0.0: go=7.0
        if not options.ge:
                ge=0.45
        else:
		try:
                	ge=float(options.ge)
		except:
			ge=0.45
		if ge<=0.0: ge=0.45
	if not options.cutoff:
                dist=4.0
        else:
                try:
                        dist=float(options.cutoff)
                except:
                        dist=4.0
		if dist<=2.0: dist=2.0
	if not options.atom:
                atom=['C3*']
	else:
		if options.atom=='P':
			atom=['P']
		else:
			atom=['C3*']
	if not options.alout:
		rout='OSARA_'+str(os.getpid())
         	outfile=rout+'.txt'
        else:
                outfile=options.alout
	if not options.biuval:
                bial=False
        else:
                bial=options.biuval
	if not options.beva:
                eva=False
        else:
                eva=options.beva
	if not options.het:
                het='N'
        else:
                het=options.het
                if het==True:
                        het='Y'
                else:
                        het='N'
	if not options.lv:
		uvl=None
	else:
		try:
			uvl=int(options.lv)
		except:
			uvl=None
        if not options.pdbout:
                if not options.alout: 
			pdbout=rout
		else:
			pdbout=None
        else:
		pdbout=options.pdbout
	if not options.slim:
		slim=500
	else:
		try:
			slim=int(options.slim)
		except:
			slim=500
	if pdbout==None: 
		sopt=' -o '+outfile+' -v -c 3.5 -g '+str(go)+' -e '+str(ge)+' -c '+str(dist)+' '
	else:
		sopt=' -o '+outfile+' -p '+pdbout+' -v -c 3.5 -g '+str(go)+' -e '+str(ge)+' -c '+str(dist)+' '
	if eva==True: sopt=sopt+' -s '
	if het=='Y': sopt=sopt+' --het '

	
	if len(args)==4:
		pdbfile1=sys.argv[1]
		chain1=sys.argv[2]
		pdbfile2=sys.argv[3]
		chain2=sys.argv[4]
		if os.path.isfile(pdbfile1)==True and os.path.isfile(pdbfile2)==True:
			f1text=open(pdbfile1,'r').read()
                        f2text=open(pdbfile2,'r').read()
                        if f1text.find('\nATOM  ')>0 and f2text.find('\nATOM  ')>0:
				atoms1=checkpdb(pdbfile1,chain1,atom)
				atoms2=checkpdb(pdbfile2,chain2,atom)
				if atoms1==atoms2 and (atoms1!=None or atoms2!=None):
					atom=atoms1
					pdb1=readPDBNucleicAcid(pdbfile1,chain1,atom)
					if pdb1.getNucleicAcidChain(0)==None:
						if chain1.islower()==True:
							nchain1=chain1.upper()
						else:
							nchain1=chain1.lower()
						pdb1=readPDBNucleicAcid(pdbfile1,nchain1,atom)
						if pdb1.getNucleicAcidChain(0)!=None: chain1=nchain1
					pdb2=readPDBNucleicAcid(pdbfile2,chain2,atom)
					if pdb2.getNucleicAcidChain(0)==None:
						if chain2.islower()==True:
                                                        nchain2=chain2.upper()
                                                else:
                                                        nchain2=chain2.lower()
                                                pdb2=readPDBNucleicAcid(pdbfile2,nchain2,atom)
                                                if pdb2.getNucleicAcidChain(0)!=None: chain2=nchain2
					tc1=pdb1.getNucleicAcidChain(0)
					tc2=pdb2.getNucleicAcidChain(0)
					if atom[0]!='P' and (tc1==None or tc2==None):
						pdb1=readPDBNucleicAcid(pdbfile1,chain1,['P'])
						pdb2=readPDBNucleicAcid(pdbfile2,chain2,['P'])
						tc1=pdb1.getNucleicAcidChain(0)
						tc2=pdb2.getNucleicAcidChain(0)
						if tc1!=None and tc2!=None: atom=['P']
					if chain1=='_':
						och1=pdb1.getNucleicAcidChain(0)
						if och1!=None: chain1=pdb1.getChains()[0]
						if chain1==' ': chain1='_'
					else:
						och1=pdb1.getNucleicAcidChainByChain(chain1)
					if chain2=='_':
						och2=pdb2.getNucleicAcidChain(0)
						if och2!=None: chain2=pdb2.getChains()[0]
						if chain2==' ': chain2='_'
					else:
						och2=pdb2.getNucleicAcidChainByChain(chain2)
					if och1!=None and och2!=None:
						len1=och1.nuacLen
						len2=och2.nuacLen
						if len1>8 and len2>8 and len1<=slim and len2<=slim:
							if bial==True:
								ss1=generate_ss(pdbfile1,chain1)
								ss2=generate_ss(pdbfile2,chain2)
								if ss1!=None and ss2!=None:
									rss1=verify_ss(ss1,chain1)
									rss2=verify_ss(ss2,chain2)
									if rss1==True and rss2==True:
										if uvl==None: uvl=3	
										#print 'python '+sara+pdbfile1+' '+chain1+' '+pdbfile2+' '+chain2+sopt+' -l 3 -b '
										com='python '+sara+' '+pdbfile1+' '+chain1+' '+pdbfile2+' '+chain2+ \
                                                                                        sopt+' -l '+str(uvl)+' -b '+' -a '+getatom(atom)
										#print com
										saraout=getstatusoutput('python '+sara+' '+pdbfile1+' '+chain1+' '+pdbfile2+' '+chain2+ \
											sopt+' -l '+str(uvl)+' -b '+' -a '+getatom(atom))
										#print 'python '+sara+' '+pdbfile1+' '+chain1+sopt+' -l '+str(uvl)+' -b '+' -a '+getatom(atom)
										if saraout[0]==0:
											if pdbout==None: 
												print 'Output file: '+outfile
											else:
												print 'Output file: '+pdbout+'.pdb'
											print saraout[1]
											if pdbout!=None:
												if os.path.isfile(outfile)==True: os.remove(outfile)
										else:
											print 'Error:'
											print saraout[1]

									else:
										if rss1==False: print 'Secondary Structure elements in '+pdbfile1+' chain '+chain1+' too short or not found'
										if rss2==False: print 'Secondary Structure elements in '+pdbfile2+' chain '+chain2+' too short or not found'
										if uvl==None: uvl=6
										com='python '+sara+' '+pdbfile1+' '+chain1+' '+pdbfile2+' '+chain2+ \
                                                                                        sopt+' -l '+str(uvl)+' -a '+getatom(atom)
										#print com
										saraout=getstatusoutput('python '+sara+' '+pdbfile1+' '+chain1+' '+pdbfile2+' '+chain2+ \
											sopt+' -l '+str(uvl)+' -a '+getatom(atom))
										if saraout[0]==0:
											if pdbout==None:
												print 'Output file: '+outfile
											else:
												print 'Output file: '+pdbout+'.pdb'
											print saraout[1]
											if pdbout!=None:
												if os.path.isfile(outfile)==True: os.remove(outfile)
										else:
											print 'Error:'
											print saraout[1]
											if pdbout!=None:
												if os.path.isfile(outfile)==True: os.remove(outfile)
								else:
									if ss1==None: print 'Secondary Structure of '+pdbfile1+' not calculated'
									if ss2==None: print 'Secondary Structure of '+pdbfile2+' not calculated'
									if uvl==None: uvl=6
									com='python '+sara+' '+pdbfile1+' '+chain1+' '+pdbfile2+' '+chain2+ \
                                                                                        sopt+' -l '+str(uvl)+' -a '+getatom(atom)
									#print com
									saraout=getstatusoutput('python '+sara+' '+pdbfile1+' '+chain1+' '+pdbfile2+' '+chain2+ \
											sopt+' -l '+str(uvl)+' -a '+getatom(atom)) 
									if saraout[0]==0:
										if pdbout==None:
											print 'Output file: '+outfile
										else:
											print 'Output file: '+pdbout+'.pdb'
										print saraout[1]
										if pdbout!=None:
											if os.path.isfile(outfile)==True: os.remove(outfile)
									else:
										print 'Error:'
										print saraout[1]
							else:
								if uvl==None: uvl=6
								com='python '+sara+' '+pdbfile1+' '+chain1+' '+pdbfile2+' '+chain2+ \
                                                                                        sopt+' -l '+str(uvl)+' -a '+getatom(atom)
								#print com
								saraout=getstatusoutput('python '+sara+' '+pdbfile1+' '+chain1+' '+pdbfile2+' '+chain2+ \
											sopt+' -l '+str(uvl)+' -a '+getatom(atom)) 
								if saraout[0]==0:
									if pdbout==None:
										print 'Output file: '+outfile
									else:
										print 'Output file: '+pdbout+'.pdb'
									print saraout[1]
									if pdbout!=None:
										if os.path.isfile(outfile)==True: os.remove(outfile)
								else:
									print 'Error:'
									print saraout[1]
						else:
							if len1<=8: print 'Error: Chain '+chain1+' in '+pdbfile1+' shorter than 9 nucleotides'
							if len2<=8: print 'Error: Chain '+chain2+' in '+pdbfile2+' shorter than 9 nucleotides'
							if len1>slim: print 'Error: Structure '+pdbfile1+' longer than '+str(slim)+' nucleotides'
							if len2>slim: print 'Error: Structure '+pdbfile2+' longer than '+str(slim)+' nucleotides'
					else:
						if och1==None: print 'Error: Chain '+chain1+' not found in '+pdbfile1
						if och2==None: print 'Error: Chain '+chain2+' not found in '+pdbfile2
				elif atoms1==None:
					print 'Error: Chain '+chain1+' not found in '+pdbfile1
				elif atoms2==None:
					print 'Error: Chain '+chain2+' not found in '+pdbfile2
				else:
					print 'Error: Different atom numeclature in the pdb files due to PDB remediates'
				
			else:
				if f1text.find('\nATOM  ')==-1: print 'Error: '+pdbfile1+' Incorrect PDB file'
                                if f2text.find('\nATOM  ')==-1: print 'Error: '+pdbfile2+' Incorrect PDB file'

		else:
			if os.path.isfile(pdbfile1)==False: print 'Error: PDB file '+pdbfile1+' not found'
			if os.path.isfile(pdbfile2)==False: print 'Error: PDB file '+pdbfile2+' not found'
		fss1=pdbfile1+'.'+chain1+'.pairs'
		fss2=pdbfile2+'.'+chain2+'.pairs'
		if os.path.isfile(fss1)==True: os.remove(fss1)
		if os.path.isfile(fss2)==True: os.remove(fss2)
	else:
		print 'python runsara.py pdbfile1 chain1 pdbfile2 chain2 [-a atom -o output -p pdb_output_file -g gapopen -e gapext -b SSTrue -c cutoff -l unit-vector -s EvaluationTrue]'
