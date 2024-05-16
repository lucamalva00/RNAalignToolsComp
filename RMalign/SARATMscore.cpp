#include "needle.h"
#include "RNAalign.h"
#include <iostream>
#include "Kabsch.h"
#include <string.h>
#include <time.h>
#include <math.h>
#include "debug_function.h"
#include "global_variable.h"
#include "SARATMscore.h"
using namespace  std;

double u_D0 = 0.0;
int main(int argc,char* argv[])
{
	clock_t t_begin,t_end;
	t_begin = clock();

	if (argc < 2)
		print_help();
	bool A_opt,B_opt,Ac_opt,Bc_opt,u_opt,a_opt,h_opt,d_opt,atom_opt,o_opt,s_opt;
	A_opt = B_opt = Ac_opt = Bc_opt = u_opt = a_opt = h_opt = d_opt = atom_opt = o_opt  = s_opt = false;
	char apdb[100],bpdb[100],ac[100],bc[100], atom[100],A0[100],superposed_file[100]; // char * apdb error will happen 
	char ss[100];
	for (int i=0; i < argc; i++)
		{
			//cout <<  i << "  "  <<  argv[i] << " " << argv[i+1]<< endl;
			if ( !strcmp(argv[i],"-B") && (i+1) < argc) {strcpy(apdb,argv[i+1]);A_opt = true;}
			if ( !strcmp(argv[i],"-A") && (i+1) < argc) {strcpy(bpdb,argv[i+1]);B_opt = true;}
			if ( !strcmp(argv[i],"-Bc") && (i+1) < argc) {strcpy(ac,argv[i+1]);Ac_opt = true;}
			if ( !strcmp(argv[i],"-Ac") && (i+1) < argc) {strcpy(bc,argv[i+1]);Bc_opt = true;}
			if ( !strcmp(argv[i],"-d") && (i+1) < argc) {u_D0 = atof(argv[i+1]);d_opt = true;}
			if ( !strcmp(argv[i],"-a") && (i+1) < argc) {strcpy(A0,argv[i+1]);a_opt = true;}
			if ( !strcmp(argv[i],"-h") && i < argc) {h_opt = true;}
			if ( !strcmp(argv[i],"-u") && (i+1) < argc) {d0_scale = atof(argv[i+1]);u_opt = true;}
			if ( !strcmp(argv[i],"-t") && (i+1) < argc) {strcpy(atom,argv[i+1]);atom_opt = true;}
			if ( !strcmp(argv[i],"-o") && (i+1) < argc) {strcpy(superposed_file,argv[i+1]);o_opt = true;}
			if ( !strcmp(argv[i],"-s") && (i+1) < argc) {strcpy(ss,argv[i+1]);s_opt = true;}
		}
	if (h_opt){
		print_help();
		exit(EXIT_FAILURE);
	}

	if (atom_opt)
		item = atom;

	if (A_opt && B_opt){
		aname = apdb;
		bname = bpdb;
	}
	if (!A_opt){
		cout << "Please provide structure B.(to be superposed)" << endl;
		exit( EXIT_FAILURE);
	}
	if (!B_opt){
		cout << "Please provide Structure A.(fixed)" << endl;
		exit(EXIT_FAILURE);
	}
	if (a_opt)
	{
		if (!strcmp(A0,"T"))
		{
			a_opt = true;
		}
		else if (!strcmp(A0,"F"))
		{
			a_opt = false;
		}
		else
		{
			cout << "Wrong value for -a , it should be T or F." << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (s_opt)
	{
		if (!strcmp(ss,"T"))
		{
			s_opt = true;
		}
		else if(!strcmp(ss,"F"))
		{
			s_opt = false;
		}
		else
		{
			cout << "Wrong value for -s , it should be T or F."  << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (u_opt)
	{
		if ( d0_scale <= 0){
			cout << "Wrong value for option -u ,it should be > 0" << endl;
			exit(EXIT_FAILURE);
		}
		else
			Lnorm_ass = d0_scale;
	}
	if (d_opt)
		if (u_D0 <= 0  ){
			cout << "Wrong value for option -d , it shoud be > 0" << endl;
			exit(EXIT_FAILURE);
		}
	if (Ac_opt)
		achain = ac;
	if (Bc_opt)
		bchain = bc;

	
	//cout << "TMScore between "<< aname << " : " << achain << "   "  << bname << " : "<< bchain << endl;
	// **********************************************************************
	// ********************load data*****************************************
	//***********************************************************************
	load_PDB(aname,bname,Ac_opt,Bc_opt,xa,ya,xresno,yresno,aseq,bseq,item,achain,bchain);
	alen = aseq.size(); 
	blen = bseq.size();
	xlen = alen;
	ylen = blen;
	minlen = getmin(xlen,ylen);
	if (locate_memory_temp(alen,blen)){
		cout << "[***]allocating memory fail!" << endl;
		exit(1);
	}
	//print_2Array(xa,alen,3);
	//print_2Array(ya,blen,3);
	//print_1Array(xresno,alen);
	//print_1Array(yresno,blen);

	//**********************************************************************
	//***************************set parameters*****************************
	//**********************************************************************
	int step = 1 , score_method = 8;
	int i ;
	int * alignment0 = new int[blen+1];
	//int * alignment  = new int[blen+1];
	double TM = -1;

	for ( i =0;i < blen; i++)
		alignment0[i] = -1;

	//double ddcc = 0.4 ;// for what ??
	//if (Lnorm <= 40) ddcc = 0.1;
	// *********************** end set *******************************

	// ***************************************************************
	//       get initial alignment with gapless threading
	//****************************************************************
	//output_fasta(aname,aseq,achain);
	//output_fasta(bname,bseq,bchain);
	
	//if (alen >= 500 || blen >= 500){cout << "Too bigger RAN absort " << endl ;exit(1);}
	string align_file = SARA_align(aname,bname,achain,bchain);
	float SARAscore1,SARAscore2;
	//float identity = read_SARA_align(alignment0,align_file,achain,bchain,xresno,yresno,SARAscore1,SARAscore2);
	float identity = read_SARA_align(alignment0,align_file,achain,bchain,SARAscore1,SARAscore2);
	cout << "The SARA normorlized score by " <<  aname << ":" << achain << " is " << SARAscore1 << endl;
	cout << "The SARA normorlized score by " << bname << ":" << bchain << " is " << SARAscore2 << endl;
	
	int n_ali = 0;
	for ( i =0;i < blen; i++)
		if (alignment0[i] >= 0)
			n_ali ++;
	
	step = 1 ; score_method = 0;
	//parameter_set4rowsearch(alen,blen);
	//TM=detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
	
	//parameter_set4final(blen);
	//TM = detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
	//   cout << "The Tmscore was normolized by length of structure A is " 
	//		<< blen << " " << bname.substr(0,4) << ":" << bchain << " is " << TM <<endl;
	double rms;
	int alignment_length;
	{
		int i, j, k;

		k=0;
		for(j=0; j<blen; j++) // loop for y structure 
		{
			i=alignment0[j];  // i means the order of the x structure
			if(i>=0) // aligned  copy coordinate
			{
				r1[k][0]=xa[i][0]; 
				r1[k][1]=xa[i][1];
				r1[k][2]=xa[i][2]; // orignal x 

				r2[k][0]=ya[j][0];
				r2[k][1]=ya[j][1];
				r2[k][2]=ya[j][2]; // orignal y

				xtm[k][0]=xa[i][0];
				xtm[k][1]=xa[i][1];
				xtm[k][2]=xa[i][2]; // x

				ytm[k][0]=ya[j][0];
				ytm[k][1]=ya[j][1];
				ytm[k][2]=ya[j][2];  // y             

				k++;
			}
			else if(i!=-1)
			{
				cout << "Wrong map!" << endl;
			}       
		}
		alignment_length = k;
		Kabsch(r1, r2, k, 1, &rms, t, u); // the import parameter is rms and t and u
	} // for the rmsd
	//for (double alpha = 0.1;alpha <= 0.5;)
	//{
	rms = sqrt(rms/alignment_length);
	double alpha = 0.5;
	parameter_set4final(blen,alpha);
	TM = detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
	    cout << "The Tmscore was normolized by length " 
			<< blen << " " << bname.substr(0,4) << ":" << bchain << " Score " << TM << " rmsd " << rms << 
			" alignment length: " << alignment_length  << 
			" iden: "<<  identity << endl;

	parameter_set4final(alen,alpha);
	TM = detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
	    cout << "The Tmscore was normolized by length " 
			<< alen << " " << aname.substr(0,4) << ":" << achain << " Score " << TM << " rmsd " << rms << 
			" alignment length: " << alignment_length << 
			" iden: " << identity << endl;
	//alpha = alpha + 0.05;
	//}
	if(a_opt)
	{
		parameter_set4final((alen+blen)*0.5,alpha);
		TM = detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
	    cout << "The Tmscore was normolized by average length " << (alen+blen)/2 <<" Score " << TM <<" rmsd " << 
			rms << " alignemnt length: " << alignment_length << endl; 
	}

	if(u_opt) // Lnorm_ass = user defined length
	{
		parameter_set4final(Lnorm_ass);
		TM = detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
		cout << "The Tmscore was normolized by user defined length is  " << TM << endl;
	}
	if(d_opt)
	{
		d0_scale = u_D0;
		parameter_set4scale(blen, d0_scale);
		TM = detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
		cout << "The Tmscore was normolized by user defined d0 and length  " << blen  << " "
			<<bname.substr(0,4)<<":"<<bchain <<" is  " << TM << " identity: " << identity <<endl;
		parameter_set4scale(alen, d0_scale);
		TM = detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
		cout << "The Tmscore was normolized by user defined d0 and length  " << alen  << " "
			<<aname.substr(0,4)<<":"<<achain <<" is  " << TM << " identity: " << identity <<endl;
		parameter_set4scale((blen+alen)*0.5, d0_scale);
		TM = detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
		cout << "The Tmscore was normolized by user defined d0 and length (A+B)/2 is  " << TM << endl;
	}

	//double * distance;
	//TMmax = TMSCORE(t,u,n_ali,distance);
	//cout << TM << " TMa " << TMmax/alen << " TMb " << TMmax/blen << " TMc " << TMmax *2/ (alen + blen)<< endl;
	t_end = clock();
	float tim = ((float)t_end - (float)t_begin)/CLOCKS_PER_SEC;
	cout << "Total running time is " << tim << endl;
	//system("rm -f *.fasta *.align");
	return 0;
}


string SARA_align(string aname,string bname,string achain,string bchain)
{
	stringstream ios;
	string outfile = aname + bname + achain + bchain;
	ios << "sara.py " << aname  << " "  << achain << " " << bname << " "  << bchain << " -a \"C3'\" -o " << 
		aname << bname << achain << bchain << " " <<  ">/dev/null 2>&1";
	//cout << ios.str() << endl;
	if (access(outfile.c_str(),0))
		system(ios.str().c_str());
	else
		return outfile;
	return outfile;
}


float read_SARA_align(int * alignment,string align_file,string achain,string bchain,int *areno,int* breno,
		float& SARAscore1,float& SARAscore2)
{
	fstream in(align_file.c_str());
	if (!in) {cout << "Can't open aling file " << align_file << endl;exit(1);}
	string AligTargetA("REMARK SARAEQVA"),AligTargetB("REMARK SARAEQVB");
	string AlenTag("REMARK SARALENA"),BlenTag("REMARK SARALENB"),SARAscore("REMARK SARASCORE");
	string temp;
	stringstream ios;
	float score1,score2,score;
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,AlenTag.size(),AlenTag) == 0)
			break;
	}
	alen = atoi(temp.substr(AlenTag.size(),temp.size() - AlenTag.size()).c_str());
	// reassign the alen beacuse SARA can deal with the NMR structure
	// but TMscore only all models in NMR, so it will decrease the TMscore
	// and, SARA align the sequence but ignore the last nt.
	// this cause the TMscore is smaller than the 1;
	// cout << alen << endl;
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,BlenTag.size(),BlenTag) == 0)
			break;
	}
	blen = atoi(temp.substr(BlenTag.size(),temp.size() - BlenTag.size()).c_str());
	//cout << blen << endl;
	
	// for the score
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,SARAscore.size(),SARAscore) == 0)
			break;
	}
	//cout << temp << endl;
	ios << temp ; ios >> temp >> temp >> score;
	//cout << score << endl;
	
	// for the alignment residue number
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,AligTargetA.size(),AligTargetA) == 0)
			break;
	}
	vector<string> AlignA,AlignB;
	split_whitespace(temp,AlignA);
	//cout << AlignA[0] << "  AAA " << AlignA[AlignA.size()-1] << endl;
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,AligTargetB.size(),AligTargetB) ==0)
			break;
	}
	split_whitespace(temp,AlignB);
	
	if (AlignA.size() != AlignB.size()){cout << "[***]Read alignment file error! Because Alignment size of A != B";exit(1);}
	//cout << "SARA alignment size: " << AlignB.size() << endl;
	for(int i=2;i<AlignB.size();i++)
	{
		int Aresindex = alen;
		int Bresindex = blen;
		for(int k = 0;k < alen;k++)
		{
			if (areno[k] == atoi(AlignA[i].c_str()))
			{
				Aresindex = k;
				break; // find the correspond residue number, break the loop
			}
		}
		for(int k = 0;k< blen; k++)
		{
			if (breno[k] == atoi(AlignB[i].c_str()))
			{
				Bresindex = k;
				break;
			}
		}
		if (Aresindex == alen || Bresindex == blen){cout << "Assign the alignment residue error! " << endl; exit(1);}
		else {alignment[Bresindex] = Aresindex;}
	}
	string AlignTag("REMARK SARAALI");
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,AlignTag.size(),AlignTag) ==0)
			break;
	}
	ios.str("");ios.clear();
	getline(in,temp);getline(in,temp);
	string SARAseqa,SARAseqb,SARAalign;
	ios << temp; ios >> temp >> temp >> SARAseqa;
	//cout << SARAseqa << endl;
	ios.str("");ios.clear();
	getline(in,temp);getline(in,temp);getline(in,temp);
	ios << temp ; ios >> temp >> temp >> SARAseqb;
	//cout << SARAseqb << endl;
	getline(in,temp);
	ios.str("");ios.clear();
	ios << temp; ios >> temp >> temp >> SARAalign;
	//cout << SARAalign << endl;
	float idencout(0);
	for(int k=0;k<SARAalign.size()-1;k++)
	{
		if (SARAseqa[k] == SARAseqb[k] && SARAseqa[k] != '-')
			idencout ++;
	}
	float identity = idencout / SARAalign.size();
	
	// for normalized SARAscore target
	string alignfile1(aname+achain+".normolized"),alignfile2(bname+bchain+".normolized");
	ios.str("");ios.clear();
	ios << "sara.py " << aname << " " << achain << " " << aname << " " << achain <<  " -a \"C3'\" -o " <<
		alignfile1 << " >/dev/null 2>&1";
	if (access(alignfile1.c_str(),0))
		system(ios.str().c_str());
	fstream in1(alignfile1.c_str());
	
	if (!in1.good()){cout << "Can't not read normolized file " << alignfile1 << endl; exit(1);}
	while(getline(in1,temp) && in1.good())
	{
		if (temp.compare(0,SARAscore.size(),SARAscore) == 0)
			break;
	}
	ios.str("");ios.clear();
	ios << temp ; ios >> temp >> temp >> score1;
	SARAscore1 = score/score1;
	
	//for normolized SARAscore  
	ios.str("");ios.clear();
	ios << "sara.py " << bname << " " << bchain << " " << bname << " " << bchain <<  " -a \"C3'\" -o " <<
		alignfile2 << " >/dev/null 2>&1";
	if (access(alignfile2.c_str(),0))
		system(ios.str().c_str());
	fstream in2(alignfile2.c_str());
	
	if (!in2.good()){cout << "Can't not read normolized file " << alignfile2 << endl; exit(1);}
	while(getline(in2,temp) && in2.good())
	{
		if (temp.compare(0,SARAscore.size(),SARAscore) == 0)
			break;
	}
	ios.str("");ios.clear();
	ios << temp ; ios >> temp >> temp >> score2;
	SARAscore2 = score/score2;
	// end normalized target
	//for(int k =0;k<blen;k++)
	//cout << k << ":" << alignment[k] << " ";
	return identity;
}

void split_whitespace(const string &str, vector<string> &result)
{
	string::size_type i, j, len = str.size();
	for (i = j = 0; i < len;)
	{
		while (i < len&&::isspace(str[i]))
			i++;
		j = i;
		while (i < len&&!::isspace(str[i]))
			i++;
		if (j < i)
		{
			result.push_back(str.substr(j, i - j));
		}
	}
}

float read_SARA_align(int * alignment,string align_file,string achain,string bchain,float& SARAscore1,float& SARAscore2)
{
	fstream in(align_file.c_str());
	if (!in) {cout << "Can't open aling file " << align_file << endl;exit(1);}
	string AlenTag("REMARK SARALENA"),BlenTag("REMARK SARALENB"),SARAscore("REMARK SARASCORE");
	string temp;
	stringstream ios;
	float score1,score2,score;
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,AlenTag.size(),AlenTag) == 0)
			break;
	}
	alen = atoi(temp.substr(AlenTag.size(),temp.size() - AlenTag.size()).c_str());
	// reassign the alen beacuse SARA can deal with the NMR structure
	// but TMscore only all models in NMR, so it will decrease the TMscore
	// and, SARA align the sequence but ignore the last nt.
	// this cause the TMscore is smaller than the 1;
	// cout << alen << endl;
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,BlenTag.size(),BlenTag) == 0)
			break;
	}
	blen = atoi(temp.substr(BlenTag.size(),temp.size() - BlenTag.size()).c_str());
	//cout << blen << endl;
	
	// for the score
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,SARAscore.size(),SARAscore) == 0)
			break;
	}
	//cout << temp << endl;
	ios << temp ; ios >> temp >> temp >> score;
	//cout << score << endl;
	
	// for the alignment residue number
	
	//cout << "SARA alignment size: " << AlignB.size() << endl;
	string AlignTag("REMARK SARAALI");
	while(getline(in,temp) && in.good())
	{
		if (temp.compare(0,AlignTag.size(),AlignTag) ==0)
			break;
	}
	ios.str("");ios.clear();
	getline(in,temp);getline(in,temp);
	string SARAseqa,SARAseqb,SARAalign;
	ios << temp; ios >> temp >> temp >> SARAseqa;
	//cout << SARAseqa << endl;
	ios.str("");ios.clear();
	getline(in,temp);getline(in,temp);getline(in,temp);
	ios << temp ; ios >> temp >> temp >> SARAseqb;
	//cout << SARAseqb << endl;
	getline(in,temp);
	ios.str("");ios.clear();
	ios << temp; ios >> temp >> temp >> SARAalign;
	//cout << SARAalign << endl;
	float idencout(0);
	for(int k=0;k<SARAalign.size()-1;k++)
	{
		if (SARAseqa[k] == SARAseqb[k] && SARAseqa[k] != '-')
			idencout ++;
	}
	float identity = idencout / SARAalign.size();
	
	// for assign the alignment
	//cout << alen << " " << blen << " " << SARAseqa << endl 
	//	<< SARAseqb << endl << SARAalign << endl;
	int aindex(0),bindex(0);
	for(int k=0;k<SARAalign.size()-1;k++)
	{
		if (SARAseqa[k] != '-' && SARAseqb[k] != '-')
		{
			alignment[bindex] = aindex;
			bindex ++ ; aindex ++;
		}
		else if(SARAseqa[k] != '-' && SARAseqb[k] == '-')
			aindex ++;
		else if(SARAseqa[k] == '-' && SARAseqb[k] == '-')
			continue;
		else if(SARAseqa[k] =='-' && SARAseqb[k] != '-')
			bindex ++;
	}

	// ending assign the alignment
	// for normalized SARAscore target
	string alignfile1(aname+achain+".normolized"),alignfile2(bname+bchain+".normolized");
	ios.str("");ios.clear();
	ios << "sara.py " << aname << " " << achain << " " << aname << " " << achain <<  " -a \"C3'\" -o " <<
		alignfile1 << " >/dev/null 2>&1";
	if (access(alignfile1.c_str(),0))
		system(ios.str().c_str());
	fstream in1(alignfile1.c_str());
	
	if (!in1.good()){cout << "Can't not read normolized file " << alignfile1 << endl; exit(1);}
	while(getline(in1,temp) && in1.good())
	{
		if (temp.compare(0,SARAscore.size(),SARAscore) == 0)
			break;
	}
	ios.str("");ios.clear();
	ios << temp ; ios >> temp >> temp >> score1;
	SARAscore1 = score/score1;
	
	//for normolized SARAscore  
	ios.str("");ios.clear();
	ios << "sara.py " << bname << " " << bchain << " " << bname << " " << bchain <<  " -a \"C3'\" -o " <<
		alignfile2 << " >/dev/null 2>&1";
	if (access(alignfile2.c_str(),0))
		system(ios.str().c_str());
	fstream in2(alignfile2.c_str());
	
	if (!in2.good()){cout << "Can't not read normolized file " << alignfile2 << endl; exit(1);}
	while(getline(in2,temp) && in2.good())
	{
		if (temp.compare(0,SARAscore.size(),SARAscore) == 0)
			break;
	}
	ios.str("");ios.clear();
	ios << temp ; ios >> temp >> temp >> score2;
	SARAscore2 = score/score2;
	// end normalized target
	//for(int k =0;k<blen;k++)
	//cout << k << ":" << alignment[k] << " ";
	return identity;
}
