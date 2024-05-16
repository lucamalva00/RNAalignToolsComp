#include "needle.h"
#include "RNAalign.h"
#include <iostream>
#include "Kabsch.h"
#include <string.h>
#include <time.h>
#include <math.h>
#include "debug_function.h"
#include "global_variable.h"
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
	output_fasta(aname,aseq,achain);
	output_fasta(bname,bseq,bchain);
	string align_file = needle_align(aname,bname,achain,bchain);
	float identity = read_needle_align_return(alignment0,align_file,aname,bname);
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
	double alpha = 0.25;
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
