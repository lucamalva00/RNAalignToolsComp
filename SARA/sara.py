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
##  contact:  Emidio Capriotti or Marc A. Marti-Renom
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
##  -------------------------------------------------------------------------------------


import os, sys, math, time
import numpy as Numeric
from optparse import OptionParser
if sys.argv[0].find('sara.py')>-1:
        ppos=0
else:
        ppos=1
pdir=sys.argv[ppos][:sys.argv[ppos].find('sara.py')]

pdbtools=pdir+'Tools/'
sys.path.append(pdbtools)
from inputoutput import AtomNucleotideList, readSecondaryStructure, calculateDistances
from inputoutput import writeoutput, writeoutput_supair, writeoutput_onlypair, writePDBOut, selectDataBySecStru, getSuperimposed3D
from maxsub import maxsub, pair_maxsub, superfind
from align import alignv, alignvx
from evaluate import Superimpose3D
from NucleicAcidTools import *

aaname1=['A','C','D','E','F','G','H','I','K','L','M', \
		'N','P','Q','R','S','T','V','W','Y','X','Z',\
		'a','c','t','g','u','n']

dna_at=['C1*','C2','C2*','C3*','C4','C4*','C5','C5*','C6','C8','N1','N2','N3','N4','N6',
	'N7','N9','O1P','O2','O2P','O3*','O4','O4*','O5*','O6','P']

dna_base=['C2','C4','C5','C6','C8','N1','N2','N3','N4','N6','N7','N9','O2','O4','O6']

		
rna_at=['P','1H2','1H2*','1H4','1H5*','1H5M','1H6','2H2','2H2*','2H4','2H5*',
        '2H5M','2H6','3H5M','C1*','C2','C2*','C3*','C4','C4*','C5','C5*','C5M',
        'C6','C8','H1','H1*','H2','H3','H3*','H3T','H4*','H5','H5T','H6','H8','N1',
        'N2','N3','N4','N6','N7','N9','O1P','O2','O2*','O2P','O3*','O4','O4*','O5*','O6',
	'C1\'','C2\'','C3\'','C4\'','C5\'','H1\'','H2\'','H3\'','H4\'','H5\'','H5\'\'',
	'HO2\'','HO3\'','HO5\'','O2\'','O3\'','O4\'','O5\'']


rna_base=['C2','C4','C5','C5M','C6','C8','N1','N2','N3','N4','N6','N7','N9','O2','O4','O6']


try:
	import psyco
	psyco.full()
except:
	pass


	
	




if __name__=="__main__":
	t0=time.time()
	parser = OptionParser('\nUsage : python %prog  pdbfile1 chain1 pdbfile2 chain2 \n[-g gep_open -e exntension -l vector -a atom -o output -p pdb_output_file -c cutoff]')
    	parser.add_option ("-g", "--gapopen", action="store",type="string",dest="go", help="Gap_Opening - Default = 7.00")
	parser.add_option ("-e", "--gapext",action="store",type="string",dest="ge", help="Gap_Extension - Default = 0.45")
	parser.add_option ("-a", "--atom",action="store",type="string",dest="atom", help="Atom Type - Default = C3*")
	parser.add_option ("-l", "--unit-length",action="store",type="string",dest="lv", help="Unit Vector Dimension")
	parser.add_option ("-o", "--output",action="store",type="string",dest="alout", help="Alignment Output File")
	parser.add_option ("-p", "--pdboutput",action="store",type="string",dest="pdbout", help="Structural Output File")
	parser.add_option ("-b", "--non-local",action="store_true",dest="biuval", help="Secondary Structure Unit Vector Alignment")
	parser.add_option ("-d", "--distance",action="store",type="string",dest="ldist", help="Distance thereshold for maxsub - Default = 3.5")
	parser.add_option ("-c", "--cutoff",action="store",type="string",dest="cutoff", help="Distance thereshold algned atoms - Default = 4.0")
	parser.add_option ("-m", "--maxsub-cutoff",action="store",type="string",dest="maxoff", help="MaxSub maximum cutoff Distance")
	parser.add_option ("-v", "--verbose",action="store_true",dest="vpar", help="Verbose Output")
	parser.add_option ("-s", "--evaluate",action="store_true",dest="beva", help="Evaluate Structural Alignment")
	parser.add_option ("--het",action="store_true",dest="het", help="Read HETATM")
	(options, args) = parser.parse_args()
	if not options.go:
		go=7.0
	else:
		try:
			go=float(options.go)
		except:
			go=7.0
		if go<0.0: go=7.0
	if not options.ge:
		ge=0.45
	else:
		try:
			ge=float(options.ge)
		except:
			ge=0.45
		if ge<0.0: ge=0.45
	if not options.atom:
		atom=['C3*']
	else:
		if options.atom.upper()=='BASE':
			atom=rna_base
		else:
			if options.atom in rna_at: 
				atom=[options.atom]
			else:
				atom=['C3*']
	if not options.lv:
		lv=6
	else:
		try:
			lv=int(options.lv)
			if lv<=2 or lv>=11: lv=6
		except:
			lv=6
	if not options.alout:
		alout='output.txt'
	else:
		alout=options.alout
	if not options.pdbout:
		pdbout=''
	else:
		pdbout=options.pdbout
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
	if not options.ldist:
		rmstol=3.5
	else:
		try:
			rmstol=float(options.ldist)
		except:
			rmstol=3.5
		if rmstol<2.0: rmstol=2.0
	if not options.cutoff:
		dist=4.0
	else:
		try:
			dist=float(options.cutoff)
		except:
			dist=4.0
		if dist<=2.0: dist=2.0
	if not options.maxoff:
		maxoff=7.0
	else:
		try:
			maxoff=float(options.maxoff)
		except:
			maxoff=7.0
	if not options.vpar:
		alpar={}
	else:
		alpar={}
		alpar['GapOpen']=go
		alpar['GapExt']=ge
		alpar['LenUV']=lv
		alpar['CutOff']=dist
		alpar['SecStructure']=bial
		alpar['Atom']=atom[0]
	if pdbout==alout:
		if pdbout!='output.txt':
			print 'Warning: Alignment output renamed to output.txt' 
			alout='output.txt'
		else:
			print 'Warning: Alignment output renamed to alignment.txt'
			alout='alignment.txt'
	zth=-100.0
	nearwin=int(lv/2)
	pos1=[]
	pos2=[]
	dc1={}
	dc2={}
	if len(args)>=4:
		pdbfile1=sys.argv[1]
		chain1=sys.argv[2]
		pdbfile2=sys.argv[3]
		chain2=sys.argv[4]
		pos1,aa1,v1,npdb1,lipos1=AtomNucleotideList(pdbfile1,chain1,atom,het)
		pos2,aa2,v2,npdb2,lipos2=AtomNucleotideList(pdbfile2,chain2,atom,het)
		if len(pos1)>lv and len(pos2)>lv:
			uv1 = calculateDistances(v1)
			uv2 = calculateDistances(v2)
			if bial==True:
				fss1=pdbfile1+'.'+chain1+'.pairs'
				fss2=pdbfile2+'.'+chain2+'.pairs'
				dc1,rdc1,liss1,cont1=readSecondaryStructure(pdbfile1,fss1,chain1,atom,het)
				dc2,rdc2,liss2,cont2=readSecondaryStructure(pdbfile2,fss2,chain2,atom,het)
			t1=time.time()
			st1=t1-t0
			dt1=time.gmtime(st1)
			print '1) Read PDB Files: '+str(dt1[3])+'h '+str(dt1[4])+'m '+str(dt1[5])+ \
			    's '+str(int(10*(st1-int(st1))))
			print '   Len A: '+str(len(pos1))+' Len B: '+str(len(pos2))
			if bial==False or dc1=={} or dc2=={}:
				if (dc1=={} or dc2=={}) and bial==True: 
					if dc1=={}: print 'Warning: Secodary Structure in '+pdbfile1+'.pairs not found'
					if dc2=={}: print 'Warning: Secodary Structure in '+pdbfile2+'.pairs not found'
					print 'Alignment without Secondaty Structure information'
					alpar['SecStructure']=False
					bial=False
				sc,iv1,iv2,match=alignv(uv1,uv2,lv,go,ge)
			else:
				npos1,naa1,nv1,nuv1=selectDataBySecStru(pos1,aa1,v1,uv1,liss1,lv)
				npos2,naa2,nv2,nuv2=selectDataBySecStru(pos2,aa2,v2,uv2,liss2,lv)
				sc,iv1,iv2,match=alignvx(nv1,nv2,nuv1,nuv2,lv,go,ge)
				#sc,iv1,iv2,match,smap=alignvcont(uv1,uv2,dc1,dc2,lv,go,ge)
			#print iv1
			#print iv2
			#print match
			t2=time.time()
			st2=t2-t1
			dt2=time.gmtime(st2)
			print '2) Structural Alignment : '+str(dt2[3])+'h '+str(dt2[4])+'m '+str(dt2[5])+ \
			    's '+str(int(10*(st2-int(st2))))
			norm=min(len(v1),len(v2))
			#print liss1
			#print liss2
			if bial==False:
				psi,nali,norm,rms,evalue,zscore,znew,pid,lsup,CP,CE,r,ir1,ir2=maxsub(pos1,pos2,aa1,aa2, \
				    v1,v2,iv1,iv2,match,len(match),norm,lv,pdbout,rmstol,zth,atom[0],maxoff)
			else:
				psi,nali,norm,rms,evalue,zscore,znew,pid,lsup,CP,CE,r,ir1,ir2=pair_maxsub(npos1,npos2,naa1,naa2, \
				    nv1,nv2,iv1,iv2,match,len(match),norm,lv,pdbout,rmstol,zth,atom[0],maxoff)
			if (CE !=None and CP !=None and r !=None).all():
				t3=time.time()
				st3=t3-t2
				dt3=time.gmtime(st3)
				print '3) Refinement of the Alignment : '+str(dt3[3])+'h '+str(dt3[4])+'m '+str(dt3[5])+ \
				    's '+str(int(10*(st3-int(st3))))
				if bial==False:
					if eva==True:
						svc1,svc2=getSuperimposed3D(npdb1,npdb2,CP,CE,r,atom[0])
						sc,siv1,siv2,smatch=Superimpose3D(svc1,svc2,dist,go=0.0,ge=0.0)
						spsi,snali,norm,srms,spid,lsup,sir1,sir2=superfind(pos1,pos2,aa1,aa2,svc1,svc2,siv1,siv2,smatch,len(smatch),norm,atom[0])
						if snali>0:
							writeoutput(pdbfile1,chain1,siv1,aa1,sir1,lipos1,dc1,pdbfile2,chain2,siv2,aa2,sir2,lipos2,dc2, \
							    lsup,spsi,snali,norm,srms,sc,spid,alout,pdbout,t0,alpar)
						else:
							print 'Output not writed: Structure Alignment too short'
							pdbout=''
					else:
						if nali>0:
							writeoutput(pdbfile1,chain1,iv1,aa1,ir1,lipos1,dc1,pdbfile2,chain2,iv2,aa2,ir2,lipos2,dc2, \
								    lsup,psi,nali,norm,rms,sc,pid,alout,pdbout,t0,alpar)
						else:
							print 'Output not writed: Structure Alignment too short'
                                                        pdbout=''						
				else:
					if eva==True:
						svc1,svc2=getSuperimposed3D(npdb1,npdb2,CP,CE,r,atom[0])
						sc,siv1,siv2,smatch=Superimpose3D(svc1,svc2,dist,go=0.0,ge=0.0)
						spsi,snali,norm,srms,spid,lsup,sir1,sir2=superfind(pos1,pos2,aa1,aa2,svc1,svc2,siv1,siv2,smatch,len(smatch),norm,atom='C3*')
						if snali>0:
							writeoutput_supair(pdbfile1,chain1,siv1,aa1,sir1,lipos1,dc1,pdbfile2,chain2,siv2,aa2,sir2,lipos2,dc2, \
							    lsup,spsi,snali,norm,srms,sc,spid,alout,pdbout,t0,alpar)
						else:
							print 'Output not writed: Structure Alignment too short'
							pdbout=''
					else:
						if nali>0:
							writeoutput_onlypair(pdbfile1,chain1,iv1,naa1,ir1,lipos1,dc1,pdbfile2,chain2,iv2,naa2,ir2,lipos2,dc2, \
							    	lsup,psi,nali,norm,rms,sc,pid,alout,pdbout,t0,alpar)
						else:
							print 'Output not writed: Structure Alignment too short'
                                                        pdbout=''
				if pdbout!='':
					writePDBOut(pdbout,npdb1,npdb2,CP,CE,r)
					os.system('cat '+alout+' '+pdbout+' > '+pdbout+'.pdb')
				if os.path.isfile(pdbout)==True: os.remove(pdbout)
			else:
				print 'Output not writed: Structure Alignment too short'
			t4=time.time()
			st4=t4-t0
			dt4=time.gmtime(st4)
			print 'Global Time: '+str(dt4[3])+'h '+str(dt4[4])+'m '+str(dt4[5])+ \
			    's '+str(int(10*(st4-int(st4))))
		else:

			if pos1==[]:
				print 'ERROR: Atom '+atom[0], 'Not found in '+pdbfile1,chain1
			elif len(pos1)<=(lv):
				print 'ERROR: Atom '+atom[0], ' lower than UV number '+pdbfile1,chain1
			if pos2==[]:
				print 'ERROR: Atom '+atom[0], 'Not found in '+pdbfile2,chain2
			elif len(pos2)<=(lv):
				print 'ERROR: Atom '+atom[0], ' lower than UV number '+pdbfile1,chain1
	else:
		print 'Usage : python sara.py  pdbfile1 chain1 pdbfile2 chain2 \n[-g gep_open -e exntension -l vector -a atom -o output -p pdb_output_file]'
