#include "RNAalign.h"
#include <iostream>
#include "Kabsch.h"
#include <string.h>
#include <time.h>
#include <math.h>
#include "global_variable.h"
#include "debug_function.h"


using namespace std;


// argument variable  add with the pro



double u_D0 = 0.0;

int main(int argc, char* argv[])
{
	
	clock_t t_begin,t_end;
	t_begin = clock();

	if (argc < 2)
		print_help();

	bool A_opt,B_opt,Ac_opt,Bc_opt,u_opt,a_opt,h_opt,d_opt,atom_opt,o_opt,s_opt,aln_opt;
	A_opt = B_opt = Ac_opt = Bc_opt = u_opt = a_opt = h_opt = d_opt = atom_opt = o_opt = s_opt = aln_opt =false;
	char apdb[100],bpdb[100],ac[100],bc[100], atom[100],A0[100],superposed_file[100]; // char * apdb error will happen 
	char ss[100],aln_file[100];
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
			if ( !strcmp(argv[i],"-i") && (i+1) < argc) {strcpy(aln_file,argv[i+1]);aln_opt = true;}
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
		parameter_set4search(alen,blen);
		int step = 40 , score_method = 8;
		int i ;
		int * alignment0 = new int[blen+1];
		int * alignment  = new int[blen+1];
		double TM,TMmax = -1;

		for ( i =0;i < blen; i++)
			alignment0[i] = -1;

		double ddcc = 0.4 ;// for what ??
		if (Lnorm <= 40) ddcc = 0.1;
		// *********************** end set *******************************

		// ***************************************************************
		//       get initial alignment with gapless threading
		//****************************************************************

		get_initial(xa,ya,alen,blen,alignment0);
		TM=detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
		//cout << TM << endl;
		
		//print_1Array(alignment0,blen); cout << endl;
		// print out is import , it can help us kown how the get_initial function;
		if ( TM > TMmax)
		{
			TMmax = TM;
		}

		TM=DP_iter(xa, ya, alen, blen, t, u, alignment, 0, 2, 30);

		if ( TM > TMmax)
		{
			TMmax = TM;
			for (int i =0 ; i < blen; i++)
				alignment0[i] = alignment[i];
		}



		// ********************************************************************
		// ******************* end first alignment*****************************
		
		
		// ***************** so here actually retain for second structure *****
		// ***************** but we here hasnot decide how to add *************
		
		if (s_opt)
		{
			char * align_home;
			align_home = getenv("RNAalign_HOME");
			string home_string = align_home;
			string command1 = home_string + "/x3dna-dssr -i=" + aname + " -o=" + aname +".x3dna-dssr 2> /dev/null";
			string command2 = home_string + "/x3dna-dssr -i=" + bname + " -o=" + bname +".x3dna-dssr 2> /dev/null";
			
			system (command1.c_str());
			system (command2.c_str());
		}
		//  for initial secondary 
		if (s_opt)
		{
			//read_dssra(aname,achain);  // only 
			
			read_a_second_structual(aname,achain);
			read_b_second_structual(bname,bchain);
			
			//read_dssrb(bname,bchain);
			double gap_open = - 1.0;
			NWDP_TM(secx,secy,xlen,ylen,gap_open,alignment);
			TM=detailed_search(xa, ya, xlen, ylen,alignment, t, u, step, score_method);
			if(TM>TMmax)
			{
				TMmax=TM;
				for(int i=0; i<ylen; i++)
				{
					alignment0[i]=alignment[i];
				}
			} 
			if(TM > TMmax*0.2)
			{
				TM=DP_iter(xa, ya, xlen, ylen, t, u, alignment, 0, 2, 30);
				if(TM>TMmax)
				{
					TMmax=TM;
					for(int i=0; i<ylen; i++)
					{
						alignment0[i]=alignment[i];
					}
				}   
			}
		}
		//print_1Array(secy,blen);cout  << endl;
		//print_1Array(secx,alen);cout  << endl;
		//exit(0);
		//print_1Array(alignment0,blen);





		// ******************************************************************
		// ******************* for local superposed *************************
		// ******************************************************************
		if ( get_initial_local(xa, ya, alen, blen, alignment))
		{
			TM = detailed_search(xa,ya,alen,blen,alignment,t,u,step,score_method);
			if ( TM > TMmax)
			{
				TMmax = TM;
				for ( int i =0;i < blen; i++)
					alignment0[i] = alignment[i];
			}
			if ( TM > TMmax*ddcc)  // this is a very magic step, but why do this ?
			{
				TM = DP_iter(xa,ya,alen,ylen,t,u,alignment,0,2,2);
				if ( TM > TMmax)
				{
					TMmax = TM;
					for( int i =0;i< blen;i++)
						alignment0[i] = alignment[i];
				}
			}
		}
		else
		{
			cout << "Waring: initial alignment from local superposition fail!" << endl << endl;
		}
		//print_1Array(alignment0,blen); cout << endl;



		// *************************************************************************
		// ****************** end alignment from local similarity*******************
		// *************************************************************************
		

		// **************************************************************************
		// ************** get inintial alignment from interget local and secondary **
		// **************************************************************************
		

		if (s_opt)
		{
			get_initial_ssplus(xa, ya, alen, blen, alignment0, alignment);	
			TM=detailed_search(xa, ya, alen, blen, alignment, t, u, step, score_method);
		}
		else
		{
			get_initial_splus(xa, ya, alen, blen, alignment0, alignment);	
			TM=detailed_search(xa, ya, alen, blen, alignment, t, u, step, score_method);
		}
		//print_1Array(alignment0,blen); cout << endl;
		if ( TM > TMmax)
		{
			TMmax = TM;
			for ( int i =0;i < blen; i++)
				alignment0[i] = alignment[i];
		}
		if ( TM > TMmax*ddcc)  // this is a very magic step, but why do this ?
		{
			TM = DP_iter(xa,ya,alen,ylen,t,u,alignment,0,2,2);
			if ( TM > TMmax)
			{
				TMmax = TM;
				for( int i =0;i< blen;i++)
					alignment0[i] = alignment[i];
			}
		}
		//print_1Array(alignment0,blen); cout << endl;



		// ************************************************************************
		// **************** end initial bound two information *********************
		// ************************************************************************



		/*********************************************************************************/
		/*        get initial alignment based on fragment gapless threading              */ 
		/*********************************************************************************/   
		
		get_initial_fgt(xa, ya, alen, blen, xresno, yresno, alignment);
		TM=detailed_search(xa, ya, alen, blen, alignment, t, u, step, score_method);
		//print_1Array(alignment0,blen); cout << endl;
		if ( TM > TMmax)
		{
			TMmax = TM;
			for ( int i =0;i < blen; i++)
				alignment0[i] = alignment[i];
		}
		if ( TM > TMmax*ddcc)  // this is a very magic step, but why do this ?
		{
			TM = DP_iter(xa,ya,alen,ylen,t,u,alignment,0,2,2);
			if ( TM > TMmax)
			{
				TMmax = TM;
				for( int i =0;i< blen;i++)
					alignment0[i] = alignment[i];
			}
		}
		//print_1Array(alignment0,blen); cout << endl;

		// ********************************************************************************
		// *********** end alignment based on fragment gapless thresding ****************88
		// ********************************************************************************




		// ********************************************************************************
		//  end of all alignment , the following code will not  any more
		// ********************************************************************************

		bool flag = false;
		for (int i =0; i < blen; i++)
			if ( alignment0[i] >= 0)
			{
				flag = true;
				break;
			}
		if ( !flag)
		{
			cout << "There is no alignment the two RNA !! " << endl;
			exit(1);
		}

		//*********************************************************************************//
		//       Detailed TMscore search engine  --> prepare for final TMscore             //
		//*********************************************************************************// 

		step = 1;
		score_method = 8;
		TM=detailed_search(xa, ya, alen, blen, alignment0, t, u, step, score_method);
		
		
		int n_ali8 , k = 0 , n_ali = 0;
		int *m1, *m2;
		double d;
		m1 = new int[alen];
		m2 = new int[blen];
		if (!m1 || !m2){
			cout << "the last allocating memory fail!!!" << endl;
			exit(0);
		}
		
		do_rotation(xa, xt, alen, t, u);
		k = 0;
		for(int j=0; j < blen; j++)
		{
			i = alignment0[j];
			if ( i >= 0)
			{
				n_ali ++ ;
				d = sqrt (dist(&xt[i][0], &ya[j][0]));
				
				if ( d <= score_d8)
				{
					m1[k] = i;
					m2[k] = j;

					xtm[k][0] = xa[k][0];
					xtm[k][1] = xa[k][1];
					xtm[k][2] = xa[k][2];
					
					
					ytm[k][0] = ya[k][0];
					ytm[k][1] = ya[k][1];
					ytm[k][2] = ya[k][2];

					k ++ ;
				}
			
			}
		}
		n_ali8 = k;

    //*********************************************************************************//
    //                               Final TMscore                                     //
    //                     Please set parameters for output                            //
    //*********************************************************************************//
    double rmsd, TM1, TM2;
	double d0_out=5.0;  
    step=1;
    score_method=0;

	double t0[3], u0[3][3];
	double d0_0, TM_0;
	double Lnorm_0=blen;
	
	// for rmsd calculate
	int j = 0;
	for( i=0; i<blen; i++)
	{
		if (alignment0[i] >= 0)
		{
			for(int k=0;k<3;k++)
			{
				xtm[j][k] = xa[alignment0[i]][k];
				ytm[j][k] = ya[i][k];
			}
			j++;
		}
	}

	Kabsch(ytm, xtm, n_ali, 0, &rmsd, t0, u0);
	rmsd = sqrt( rmsd / n_ali);
	// end RMSD
	//cout << " rmsd: " << rmsd << endl;
	parameter_set4final(Lnorm_0);
	d0A=d0;
	d0_0=d0A;
	//cout << "Normolized by  length of A." << endl;
	TM1=TMscore8_search(xtm, ytm, n_ali, t0, u0, step, score_method, &rmsd);
	//TM1 = detailed_search(xa, ya, alen, blen, alignment, t0, u0, step, score_method);
	
	TM_0=TM1;

	//normalized by length of structure B
	parameter_set4final(alen+0.0);
	d0B=d0;
	//cout << "Normolized by  length of A." << endl;
	TM2=TMscore8_search(xtm, ytm, n_ali, t, u, step, score_method, &rmsd);
	//TM2 = detailed_search(xa, ya, alen, blen, alignment, t, u, step, score_method);





	if(a_opt)
	{
		//normalized by average length of structures A, B
		Lnorm_0=(alen+blen)*0.5;
		parameter_set4final(Lnorm_0);
		d0a=d0;
		d0_0=d0a;
		TM3=TMscore8_search(xtm, ytm, n_ali, t0, u0, step, score_method, &rmsd);
		//TM3 = detailed_search(xa, ya, alen, blen, alignment, t, u, step, score_method);
		TM_0=TM3;
	}
	if(u_opt)
	{	
		//normalized by user assigned length		
		parameter_set4final(Lnorm_ass);		
		d0u=d0;		
		d0_0=d0u;
		Lnorm_0=Lnorm_ass;
		TM4=TMscore8_search(xtm, ytm, n_ali, t, u, step, score_method, &rmsd);	
		TM_0=TM4;
	}
	if(d_opt)
	{
		//scaled by user assigned d0
		d0_scale = u_D0;
		parameter_set4scale(blen, d0_scale);
		d0_out=d0_scale;
		d0_0=d0_scale;
		//Lnorm_0=ylen;
		Lnorm_d0=Lnorm_0;
		TM5=TMscore8_search(xtm, ytm, n_ali, t0, u0, step, score_method, &rmsd);	
		TM_0=TM5;
	}
	//print_1Array(alignment0,blen); cout << endl;
	//print_1Array(m1,alen); cout << endl;
	//print_1Array(m2,blen); cout << endl;
	output_alignment(aname, bname, alen, blen, t0, u0, TM1, TM2, rmsd, d0_out, m1, m2, n_ali8, n_ali, 
			TM_0, Lnorm_0, d0_0,xresno,yresno,a_opt,u_opt,d_opt,o_opt,superposed_file);


	//cout << " rmsd: " << rmsd << endl;
	t_end = clock();
	float tim = ((float)t_end - (float)t_begin)/CLOCKS_PER_SEC;
	//cout << t_begin << " " <<  t_end << " " << CLOCKS_PER_SEC;
	cout << endl << "Total running time is " << tim << endl; 
	//// ---   free the memory ----

	free_memory();

	delete [] alignment;
	delete [] alignment0;
	delete [] m1;
	delete [] m2;
	return 0;
}
