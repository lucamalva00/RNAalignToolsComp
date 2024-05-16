#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include "RNAalign.h"
#include "function.h"
#include <string.h>
void PrintErrorAndQuit(char* sErrorString)
{
	cout << sErrorString << endl;
	exit(0);
}

template<class T>
void print_1Array(T * arrat, int size)
{
	for (int i=0; i<size;i++)
		cout << arrat[i] << " ";
}


void print_1Array(int *array, int size)
{	
	for (int i = 0; i < size ; i++)
		cout << array[i] << "  ";
	cout << endl;
}

void print_1Array(char *array, int size)
{	
	for (int i = 0; i < size ; i++)
		cout << array[i] << "  ";
	cout << endl;
}

bool locate_memory_temp (int alen, int blen){
	
	int minlen =  getmin(alen,blen);
	NewArray(&r1, minlen, 3);
	NewArray(&r2, minlen, 3);
    NewArray(&xtm, minlen, 3);
	NewArray(&ytm, minlen, 3);
	NewArray(&xt, alen, 3);
	NewArray(&score, alen+1, blen+1);
	NewArray(&path, alen+1, blen+1);
	NewArray(&val, alen+1, blen+1);
	secx = new int [alen];
	secy = new int [blen];
	if (r1 && r2 && xtm && ytm && xt && score && path && val && secx && secy)
		return false; 
	else
		return true;
}

void print_help(){
	cout << endl;
	cout << " ****************************************************** " << endl;
	cout << " RNAalign : A RNA structural alignment algorithm " << endl;
	cout << " any bugs report by email zhengjinfang1220@gmail.com " << endl;
	cout << " ****************************************************** " << endl ;
	cout << endl;
	cout << " Usage :: " << " RNAalign  [Opition]" << endl << endl <<
		" -A input filename of structure A, PDB format " << endl <<
		" -B input filename of structure B, PDB format " << endl <<
		" -Ac chain ids of structure A, string format, " <<
		" if not assigned, RNAalign will read all chains" << endl <<
		" -Bc chain ids of structure B, string format, "
		" if not assigned, RNAalign will read all chains" << endl <<
		" -u RMscore normalized by user defined length, " << 
		"!!!warning it shoud be >= min length of the two structure" << endl <<
		" -t RNA atom item to read; default atom C3'  " << endl << 
		" -a RMscore normalized by the average of two structure, T or F(default T)" << endl <<
		" -o the output file of structre B superposed on A " << endl <<
		" -d user define the d0 used in nomorlizing RMscore " << endl << 
		" -s add second structral information, T or F (default F)" << endl <<
		" -h print this help" << endl << endl << 
		" Example usage:" << endl <<
		" RNAalign -A PDB1.pdb -Ac A -B PDB2.pdb -Bc A " << endl  <<
		" RNAalign -A PDB1.pdb -B PDB2.pdb -t \"C4'\" "<< endl << 
		" RNAalign -A PDB1.pdb -B PDB2.pdb -u 50 -d 5.0 -o PDB1.sup " << endl ;
	cout << " RNAalign -A PDB1.pdb -Ac A -B PDB2.pdb -Bc A -s T  " << endl << endl;
	
	
	exit(EXIT_SUCCESS);
}


void load_PDB(string& aname,string& bname,bool Ac_opt,bool Bc_opt,double ** &acoor,
		double ** &bcoor,int* &aresno,int * &bresno,string& aseq,string& bseq,
		string &item,string &achain, string &bchain)
{
	
	ifstream ifile1(aname.c_str());
	ifstream ifile2(bname.c_str());
	
	if (!ifile1){
		cout << " cann't open " << aname << endl;
		exit(EXIT_FAILURE);
	}
	if (!ifile2){
		cout << " cann't open " << bname << endl;
		exit(EXIT_FAILURE);
	}

	int alen,blen;
	alen = get_APDB_length(ifile1,Ac_opt,item,achain);
	blen = get_BPDB_length(ifile2,Bc_opt,item,bchain); 


	///   allocate memory for coordrites
	NewArray(&acoor,alen,3);
	NewArray(&bcoor,blen,3);
	///  allocate memory for 
	aresno = new int[alen];
	bresno = new int[blen];

	chidx = new char[alen+1];
	chidy = new char[blen+1];
	chidx[alen]='\0';
	chidy[blen]='\0';
	ifstream ifile3(aname.c_str());
	ifstream ifile4(bname.c_str());
	

	read_PDB(ifile3,acoor,achain,aresno,aseq,chidx,item,Ac_opt);
	read_PDB(ifile4,bcoor,bchain,bresno,bseq,chidy,item,Bc_opt);

	//DeleteArray(&acoor,alen);
	//DeleteArray(&bcoor,blen);

}

void get_xyz(string line,double& x,double& y,double& z,int& no,char& seq)
{
	x = atof(line.substr(30,8).c_str());
	y = atof(line.substr(38,8).c_str());
	z = atof(line.substr(46,8).c_str());
	no = atoi(line.substr(22,4).c_str());
	seq = line[19];

}
void get_xyz(string line,double& x,double& y,double& z)
{
	x = atof(line.substr(30,8).c_str());
	y = atof(line.substr(38,8).c_str());
	z = atof(line.substr(46,8).c_str());
}

void read_PDB(ifstream& infile, double** &coor,string achain,int *&resno,string & seq,char* &chid,string item,bool c_opt)
{
	infile.seekg(ios::beg);
	string temp;
	int length = 0, no ;
	char aseq ; 
	double x,y,z;
	if (c_opt){
		while (infile.good()){
			getline(infile,temp);
			if (temp.compare(0,4,"ATOM") == 0){
				if (temp.compare(13,item.length(),item) == 0){
					char chainid = temp.at(21);
					string::size_type pos;
					pos = achain.find(chainid);
					if (pos != achain.npos){
						get_xyz(temp,x,y,z,no,aseq);
						coor[length][0] = x ;
						coor[length][1] = y ;
						coor[length][2] = z ;
						resno[length] = no;
						seq.push_back(aseq);
						chid[length] = temp.at(21);
						length ++;
					}
				}
			}
		}
	}
	else{
		while (infile.good()){
			getline(infile,temp);
			if (temp.compare(0,4,"ATOM") == 0){
				if (!temp.compare(13,item.size(),item)){
					get_xyz(temp,x,y,z,no,aseq);
					coor[length][0] = x ;
					coor[length][1] = y ;
					coor[length][2] = z ;
					resno[length] = no;
					seq.push_back(aseq);
					chid[length] = temp.at(21);
					length ++;
				}
			}
		}
	}
	if (length == 0)
		cout << "read 0 atom " << endl;

}


int get_APDB_length(ifstream& infile, bool c_opt,string& item, string & achain){
	infile.seekg(ios::beg);
	int length = 0;
	string temp;
	if (c_opt){
		while (infile.good()){
			getline(infile,temp);
			if (temp.compare(0,4,"ATOM") == 0){
				if (temp.compare(13,item.length(),item) == 0){
					char chainid = temp.at(21);
					string::size_type pos;
					pos = achain.find(chainid);
					if (pos != achain.npos){
						length ++;
					}
				}
			}
		}
	}
	else{
		while (infile.good()){
			getline(infile,temp);
			if (temp.compare(0,4,"ATOM") == 0){
				if (!temp.compare(13,item.length(),item)){
					length ++;
				}
			}
		}
	}

	if (length == 0){
		cout << "read 0 atom from PDB file." << endl;
		exit(EXIT_FAILURE);
	}
	return length;
}

int get_BPDB_length(ifstream& infile, bool c_opt,string &item,string &bchain){
	infile.seekg(ios::beg);
	int length = 0;
	string temp;
	if (c_opt){
		while (infile.good()){
			getline(infile,temp);
			if (temp.compare(0,4,"ATOM") == 0){
				if (temp.compare(13,item.size(),item) == 0){
					char chainid = temp.at(21);
					string::size_type pos;
					pos = bchain.find(chainid);
					if (pos != bchain.npos){
						length ++;
					}
				}
			}
		}
	}
	else{
		while (infile.good()){
			getline(infile,temp);
			if (temp.compare(0,4,"ATOM") == 0){
				if (temp.compare(13,item.size(),item) == 0){
					length ++;
				}
			}
		}
	}

	if (length == 0){
		cout << "read 0 atom from PDB file." << endl;
		exit(EXIT_FAILURE);
	}
	return length;
}


void print_2Array(double** & array2,int n1,int n2){
	for(int i = 0;i < n1 ;i++){
		for (int j = 0 ; j < n2; j++){
			cout << array2[i][j] << "   ";
		}
		cout << endl;
	}
}

int getmin(int alen,int blen){

	if (alen>blen)
		return blen;
	else
		return alen;
};


void NWDP_TM(int len1, int len2, double gap_open, int j2i[])
{
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
    //Input: score[1:len1, 1:len2], and gap_open
    //Output: j2i[1:len2] \in {1:len1} U {-1}
	//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical

	int i, j;
	double h, v, d;

	//initialization
	val[0][0]=0;
	for(i=0; i<=len1; i++)
	{
		val[i][0]=0;
		path[i][0]=false; //not from diagonal
	}

	for(j=0; j<=len2; j++)
	{
		val[0][j]=0;
		path[0][j]=false; //not from diagonal
		j2i[j]=-1;	//all are not aligned, only use j2i[1:len2]
	}      

	
	//decide matrix and path
	for(i=1; i<=len1; i++)	
	{	
		for(j=1; j<=len2; j++)
		{
			d=val[i-1][j-1]+score[i][j]; //diagonal   score is the score 

			//symbol insertion in horizontal (= a gap in vertical)
			h=val[i-1][j];
			if(path[i-1][j]) //aligned in last position
				h += gap_open;				

			//symbol insertion in vertical
			v=val[i][j-1];
			if(path[i][j-1]) //aligned in last position
				v += gap_open;


			if(d>=h && d>=v)
			{
				path[i][j]=true; //from diagonal
				val[i][j]=d;
			}
			else 
			{
				path[i][j]=false; //from horizontal
				if(v>=h)
					val[i][j]=v;
				else					
					val[i][j]=h;
			}
		} //for i
	} //for j

	//trace back to extract the alignment
	i = len1;
    j = len2;
	while(i>0 && j>0)
	{
		if(path[i][j]) //from diagonal
		{
			j2i[j-1]=i-1;		
			i--;
			j--;
		}
		else 			
		{
			h=val[i-1][j];
			if(path[i-1][j]) h +=gap_open;

			v=val[i][j-1];
			if(path[i][j-1]) v +=gap_open;

			if(v>=h)
				j--;
			else
				i--;
		}
	}	
}

void NWDP_TM(double **x, double **y, int len1, int len2, double t[3], double u[3][3], double d02, double gap_open, int j2i[])
{
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
    //Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
    //Output: j2i[1:len2] \in {1:len1} U {-1}
	//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical

	int i, j;
	double h, v, d;

	//initialization
	val[0][0]=0;
	for(i=0; i<=len1; i++)
	{
		val[i][0]=0;
		path[i][0]=false; //not from diagonal
	}

	for(j=0; j<=len2; j++)
	{
		val[0][j]=0;
		path[0][j]=false; //not from diagonal
		j2i[j]=-1;	//all are not aligned, only use j2i[1:len2] // rigth
	}      
	double xx[3], dij;

	// path np algorithm got it
	//decide matrix and path
	for(i=1; i<=len1; i++)	
	{	
		transform(t, u, &x[i-1][0], xx);
		for(j=1; j<=len2; j++)
		{
			//d=val[i-1][j-1]+score[i][j]; //diagonal
			dij=dist(xx, &y[j-1][0]);    					
			d=val[i-1][j-1] +  1.0/(1+dij/d02);

			//symbol insertion in horizontal (= a gap in vertical)
			h=val[i-1][j];
			if(path[i-1][j]) //aligned in last position
				h += gap_open;				

			//symbol insertion in vertical
			v=val[i][j-1];
			if(path[i][j-1]) //aligned in last position
				v += gap_open;


			if(d>=h && d>=v)
			{
				path[i][j]=true; //from diagonal
				val[i][j]=d;
			}
			else 
			{
				path[i][j]=false; //from horizontal
				if(v>=h)
					val[i][j]=v;
				else					
					val[i][j]=h;
			}
		} //for i
	} //for j

	//trace back to extract the alignment
	i=len1;
    j=len2;
	while(i>0 && j>0)
	{
		if(path[i][j]) //from diagonal
		{
			j2i[j-1]=i-1;		
			i--;
			j--;
		}
		else 			
		{
			h=val[i-1][j];
			if(path[i-1][j]) h +=gap_open;

			v=val[i][j-1];
			if(path[i][j-1]) v +=gap_open;

			if(v>=h)
				j--;
			else
				i--;
		}
	}	
}

//+ss
void NWDP_TM(int *secx, int *secy, int len1, int len2, double gap_open, int j2i[])
{
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
    //Input: secondary structure secx, secy, and gap_open
    //Output: j2i[1:len2] \in {1:len1} U {-1}
	//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical

	int i, j;
	double h, v, d;

	//initialization
	val[0][0]=0;
	for(i=0; i<=len1; i++)
	{
		val[i][0]=0;
		path[i][0]=false; //not from diagonal
	}

	for(j=0; j<=len2; j++)
	{
		val[0][j]=0;
		path[0][j]=false; //not from diagonal
		j2i[j]=-1;	//all are not aligned, only use j2i[1:len2]
	}      
	
	//decide matrix and path
	for(i=1; i<=len1; i++)	
	{	
		for(j=1; j<=len2; j++)
		{
			//d=val[i-1][j-1]+score[i][j]; //diagonal			
			if(secx[i-1]==secy[j-1])
			{
				d=val[i-1][j-1] + 1.0;
			}
			else
			{
				d=val[i-1][j-1];
			}

			//symbol insertion in horizontal (= a gap in vertical)
			h=val[i-1][j];
			if(path[i-1][j]) //aligned in last position
				h += gap_open;				

			//symbol insertion in vertical
			v=val[i][j-1];
			if(path[i][j-1]) //aligned in last position
				v += gap_open;


			if(d>=h && d>=v)
			{
				path[i][j]=true; //from diagonal
				val[i][j]=d;
			}
			else 
			{
				path[i][j]=false; //from horizontal
				if(v>=h)
					val[i][j]=v;
				else					
					val[i][j]=h;
			}
		} //for i
	} //for j

	//trace back to extract the alignment
	i=len1;
    j=len2;
	while(i>0 && j>0)
	{
		if(path[i][j]) //from diagonal
		{
			j2i[j-1]=i-1;		
			i--;
			j--;
		}
		else 			
		{
			h=val[i-1][j];
			if(path[i-1][j]) h +=gap_open;

			v=val[i][j-1];
			if(path[i][j-1]) v +=gap_open;

			if(v>=h)
				j--;
			else
				i--;
		}
	}	
}


void free_memory()
{
	DeleteArray(&path, xlen+1);
	DeleteArray(&val, xlen+1);
	DeleteArray(&score, xlen+1);
	DeleteArray(&xa, xlen);
	DeleteArray(&xt, xlen);
	DeleteArray(&ya, ylen);
	DeleteArray(&r1, minlen);
	DeleteArray(&r2, minlen);
	DeleteArray(&xtm, minlen);
	DeleteArray(&ytm, minlen);
    
 	system("rm -f dssr-*");   
    delete [] seqx;
    delete [] seqy;
    delete [] secx;
    delete [] secy;
    delete [] xresno;
    delete [] yresno;
    
}


//     1, collect those residues with dis<d;
//     2, calculate TMscore
int score_fun8( double **xa, 
                double **ya, 
                int n_ali,
                double d,
                int i_ali[], 
                double *score1,
                int score_sum_method
              )
{
    double score_sum=0, di;
    double d_tmp=d*d;
	double d02=d0*d0;
	double score_d8_cut = score_d8*score_d8;
    
    int i, n_cut, inc=0;

    while(1)
    {
		n_cut=0;
        score_sum=0;
        for(i=0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if(di<d_tmp)
            {
                i_ali[n_cut]=i;
                n_cut++;
            }
            if(score_sum_method==8)
            {				
                if(di<=score_d8_cut)
                {					
                    score_sum += 1/(1+di/d02);
                }                
            }
            else
            {
				score_sum += 1/(1+di/d02);
            }
        }
        //there are not enough feasible pairs, reliefe the threshold 		
        if(n_cut<3 && n_ali>3)
        {
			inc++;
			double dinc=(d+inc*0.5);
			d_tmp = dinc * dinc;
        }
        else
        {
            break;
        }

    }  

    *score1=score_sum/Lnorm;

    return n_cut;
}

// TMscore search engine
// input:   two aligned vector sets: x, y
//          scale parameter d0
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8
//                                  
//          
// output:  the best rotaion matrix t0, u0 that results in highest TMscore
double TMscore8_search( double **xtm, 
                        double **ytm,
                        int Lali, // aligned length of residue
                        double t0[3],
                        double u0[3][3],
                        int simplify_step,
                        int score_sum_method,
                        double *Rcomm
                       )
{   
    int i, m;
    double score_max, score, rmsd;    
    const int kmax=Lali;    
    int k_ali[kmax], ka, k;
    double t[3];
    double u[3][3];
	double d;
	

	//iterative parameters
	int n_it=20;            //maximum number of iterations
    const int n_init_max=6; //maximum number of different fragment length  fragemnt type
    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   Lali/2^i 
    int L_ini_min=4;  // but the fragmetn must >= 4
    if(Lali<4) L_ini_min=Lali;   // if the paired length > 4 , L_ini_min = 4, else = 4,or,3,2
    int n_init=0, i_init;
    for(i=0; i<n_init_max-1; i++) // for loop for 0,1,2,3... search fragment type
    {
        n_init++;
        L_ini[i]=(int) (Lali/pow(2.0, (double) i));// for fragemnt length
        if(L_ini[i]<=L_ini_min)
        {
            L_ini[i]=L_ini_min;
            break;
        }
    }
    if(i==n_init_max-1) // i mean the number of  fragment 
    {
        n_init++;
        L_ini[i]=L_ini_min;
    }
    // all the above is set the L_inin   the fragment lebgth of searching
    score_max=-1;
    //find the maximum score starting from local structures superposition
    int i_ali[kmax], n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment
    for(i_init=0; i_init<n_init; i_init++) //loop for the fragment
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag; // biduishang - pianduan  ye jiushi shengyu de canji geshu
      
        i=0;   
        while(1)
        {
            //extract the fragment starting from position i 
            ka=0;
            for(k=0; k<L_frag; k++)
            {
				int kk=k+i; // k , i  if L_frag = the pair residue
                r1[k][0]=xtm[kk][0];  
                r1[k][1]=xtm[kk][1]; 
                r1[k][2]=xtm[kk][2];   
                
                r2[k][0]=ytm[kk][0];  
                r2[k][1]=ytm[kk][1]; 
                r2[k][2]=ytm[kk][2];
                
                k_ali[ka]=kk; // kk for what?
                ka++; // for the coorditate
            }
            
            //extract rotation matrix based on the fragment
            Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
            if(i_init==0)
            {
                *Rcomm=sqrt(rmsd/Lali);// normalized the rmsd
            }
			//print_2Array(xtm,Lali,3) ; cout << endl;
			//print_2Array(xt,Lali,3); cout << endl;
			//cout << Lali << endl;
            do_rotation(xtm, xt, Lali, t, u);
            d=d0_search-1;
            n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, score_sum_method);
            if(score>score_max)
            {
                score_max=score;
                *Rcomm=sqrt(rmsd/Lali);// normalized the rmsd
                
                //save the rotation matrix
                for(k=0; k<3; k++)
                {
                    t0[k]=t[k];
                    u0[k][0]=u[k][0];
                    u0[k][1]=u[k][1];
                    u0[k][2]=u[k][2];
                }
            }
            
            //try to extend the alignment iteratively            
            d=d0_search+1; //d set the residue distance 
            for(int it=0; it<n_it; it++)            
            {
                ka=0;
                for(k=0; k<n_cut; k++)
                {
                    m=i_ali[k];
                    r1[k][0]=xtm[m][0];  
                    r1[k][1]=xtm[m][1]; 
                    r1[k][2]=xtm[m][2];
                    
                    r2[k][0]=ytm[m][0];  
                    r2[k][1]=ytm[m][1]; 
                    r2[k][2]=ytm[m][2];
                    
                    k_ali[ka]=m;
                    ka++;
                } 
                //extract rotation matrix based on the fragment                
                Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
                do_rotation(xtm, xt, Lali, t, u);
                n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, score_sum_method);
                if(score>score_max)
                {
                	*Rcomm=sqrt(rmsd/Lali);// normalized the rmsd
                    score_max=score;

                    //save the rotation matrix
                    for(k=0; k<3; k++)
                    {
                        t0[k]=t[k];
                        u0[k][0]=u[k][0];
                        u0[k][1]=u[k][1];
                        u0[k][2]=u[k][2];
                    }                     
                }
                
                //check if it converges                 
			
                if(n_cut==ka)
                {				
                    for(k=0; k<n_cut; k++)
                    {
                        if(i_ali[k]!=k_ali[k])
						{
							break;
						}
                    }
                    if(k==n_cut)
                    {						
                        break; //stop iteration
                    }
                }                                                               
            } //for iteration            

			if(i<iL_max)
			{
				i=i+simplify_step; //shift the fragment		
				if(i>iL_max) i=iL_max;  //do this to use the last missed fragment
			}
			else if(i>=iL_max)
			{
				break;
			}
        }//while(1)
        //end of one fragment
		//cout << "flagment: " << L_frag << " Align lengt: " << Lali << " RMscore: "<< score << endl;
    }//for(i_init
    return score_max;
}

//Comprehensive TMscore search engine
// input:   two vector sets: x, y
//          an alignment invmap0[] between x and y
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8          
// output:  the best rotaion matrix t, u that results in highest TMscore
double detailed_search( double **x,
                        double **y, 
                        int x_len, 
                        int y_len, 
                        int invmap0[],
                        double t[3],
                        double u[3][3],
                        int simplify_step,
                        int score_sum_method                        
                       )
{
    //x is model, y is template, try to superpose onto y
    int i, j, k;     
    double tmscore;
    double rmsd(0);



    k=0;
    for(i=0; i<y_len; i++) 
    {
        j=invmap0[i];
        if(j>=0) //aligned
        {
            xtm[k][0]=x[j][0];
            xtm[k][1]=x[j][1];
            xtm[k][2]=x[j][2];
                
            ytm[k][0]=y[i][0];
            ytm[k][1]=y[i][1];
            ytm[k][2]=y[i][2];
            k++;
        }
    }
    
    //detailed search 40-->1
    tmscore=TMscore8_search(xtm, ytm, k, t, u, simplify_step, score_sum_method, &rmsd);  
   // add for calculating the distance of corresponding resiudes
	k=0;
    for(i=0; i<y_len; i++) 
    {
        j=invmap0[i];
        if(j>=0) //aligned
        {
            xtm[k][0]=x[j][0];
            xtm[k][1]=x[j][1];
            xtm[k][2]=x[j][2];
                
            ytm[k][0]=y[i][0];
            ytm[k][1]=y[i][1];
            ytm[k][2]=y[i][2];
            k++;
        }
    }
	rmsd = 0;
	//Kabsch(xtm,ytm,k,1,&rmsd,t,u);
	do_rotation(xtm,xt,k,t,u);
	//printf("ditance: ");
    //for(i=0;i<k;i++)
	//	printf("%5.2f ",sqrt(dist(ytm[i],xt[i])));
	//printf("%i\n",k);
	// tmscore and rmsd
    return tmscore;
}



//compute the score quickly in three iterations
double get_score_fast(double **x, double **y, int x_len, int y_len, int invmap[])
{
    double rms, tmscore, tmscore1, tmscore2; // three TMscore
    int i, j, k;

    k=0;
    for(j=0; j<y_len; j++) // loop for y structure 
    {
        i=invmap[j];  // i means the order of the x structure
        if(i>=0) // aligned  copy coordinate
        {
            r1[k][0]=x[i][0]; 
            r1[k][1]=x[i][1];
            r1[k][2]=x[i][2]; // orignal x 

            r2[k][0]=y[j][0];
            r2[k][1]=y[j][1];
            r2[k][2]=y[j][2]; // orignal y
            
            xtm[k][0]=x[i][0];
            xtm[k][1]=x[i][1];
            xtm[k][2]=x[i][2]; // x
            
            ytm[k][0]=y[j][0];
            ytm[k][1]=y[j][1];
            ytm[k][2]=y[j][2];  // y             
            
            k++;
        }
        else if(i!=-1)
        {
            cout << "Wrong map!" << endl;
        }       
    }
    Kabsch(r1, r2, k, 1, &rms, t, u); // the import parameter is rms and t and u
    // evlate the rmsd
    //evaluate score   
    double di;
	const int len=k; // k = the length of aligned 
    //double *dis;   // dis 
	double dis[len];
	double d00=d0_search ; // set for calkulate the tmscore second iterator
	double d002=d00*d00;
	double d02=d0*d0; //  set
	
    int n_ali=k;
	double xrot[3];
	tmscore=0;
	
	for(k=0; k<n_ali; k++)
	{
        transform(t, u, &xtm[k][0], xrot);  // tranform the origial coordinate to xrot      
        di=dist(xrot, &ytm[k][0]); // 
        dis[k]=di; // the aligned reside distance
        tmscore += 1/(1+di/d02); // tmscore   so that the TMscore corresnate to the distance
    }


	//tmscore = TMSCORE(t,u,n_ali,dis);
	//print_1Array(dis,n_ali); cout << tmscore << endl;//exit(0);
   
   
   //second iteration  what this mean?
    double d002t=d002;
    while(1)
    {
		j=0;
        for(k=0; k<n_ali; k++)
        {            
            if(dis[k]<=d002t) // manzu diejia hou de juli de canji  jixu diejie
            {
                r1[j][0]=xtm[k][0];
                r1[j][1]=xtm[k][1];
                r1[j][2]=xtm[k][2];
                
                r2[j][0]=ytm[k][0];
                r2[j][1]=ytm[k][1];
                r2[j][2]=ytm[k][2];
                
                j++;
            }
        }
        //there are not enough feasible pairs, relieve the threshold 
        if(j<3 && n_ali>3)
        {
            d002t += 0.5; // if the aligend is  less tha 3 then add the distance 
        }
        else
        {
            break; // break   j >=3 and 
        }
    }
    
    if(n_ali!=j) // zengjia juli tiaojian zhihou  bidui shang de can ji bianxiaole 
    {
        Kabsch(r1, r2, j, 1, &rms, t, u);// calkulate the rmsd again
    	tmscore1=0; // another tmsocre
    	for(k=0; k<n_ali; k++)
    	{
            transform(t, u, &xtm[k][0], xrot);        
            di=dist(xrot, &ytm[k][0]);
            dis[k]=di;
            tmscore1 += 1/(1+di/d02);
        }
        
        //third iteration
        d002t=d002+1; // add the d002 agian
       
        while(1)
        {
			j=0;
            for(k=0; k<n_ali; k++)
            {            
                if(dis[k]<=d002t)
                {
                    r1[j][0]=xtm[k][0];
                    r1[j][1]=xtm[k][1];
                    r1[j][2]=xtm[k][2];
                    
                    r2[j][0]=ytm[k][0];
                    r2[j][1]=ytm[k][1];
                    r2[j][2]=ytm[k][2];
                                        
                    j++;
                }
            }
            //there are not enough feasible pairs, relieve the threshold 
            if(j<3 && n_ali>3)
            {
                d002t += 0.5;
            }
            else
            {
                break;
            }
        }

        //evaluate the score
        Kabsch(r1, r2, j, 1, &rms, t, u);
        tmscore2=0;
        for(k=0; k<n_ali; k++)
        {
            transform(t, u, &xtm[k][0], xrot);
            di=dist(xrot, &ytm[k][0]);
            tmscore2 += 1/(1+di/d02);
        }    
    }
    else
    {
        tmscore1=tmscore;
        tmscore2=tmscore;
    }
    
      
    if(tmscore1>=tmscore) tmscore=tmscore1;
    if(tmscore2>=tmscore) tmscore=tmscore2;

    // just calkulate a tmscore
    return tmscore; // no need to normalize this score because it will not be used for latter scoring
}


//perform gapless threading to find the best initial alignment
//input: x, y, x_len, y_len
//output: y2x0 stores the best alignment: e.g., 
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
double get_initial( double **x, 
                    double **y, 
                    int x_len,
                    int y_len, 
                    int *y2x
                   )
{
    int min_len=getmin(x_len, y_len);//  if x_len, y_len = 8 min_len =8
    if(min_len<=5) cout << "Sequence is too short <=5!" << endl;
    
    int min_ali= min_len/2;              //minimum size of considered fragment   min_ali = 4  
    if(min_ali<=5)  min_ali=5;    // min_ali = 5
    int n1, n2;
    n1 = -y_len+min_ali;  // n1 = -8 + 4 = -4
    n2 = x_len-min_ali; // n2 = 8- 4 = 4

    int i, j, k, k_best;
    double tmscore, tmscore_max=-1;

    k_best=n1;  // k_best = -4
    for(k=n1; k<=n2; k++) // for -4,-3,-2,-1,0,1,2,3,4
    {
        //get the map  get the map 
        for(j=0; j<y_len; j++) // for j =0,1,2,3,4,5,6,7,8  in y_len
        {
            i=j+k; //  -y_len+min_ali),x_len-min_ali+y_len]
            if(i>=0 && i<x_len)//  y2x   i mean the order of the x structure sequence 
            {
                y2x[j]=i; // j mean the order of the y structure sequence  y(j) align to x(i)
            }
            else
            {
                y2x[j]=-1; // -1 mean this not aligned  y(j)  not aligned to  
            }
        }
        
        //evaluate the map quickly in three iterations
		//this is not real tmscore, it is used to evaluate the goodness of the initial alignment
        tmscore=get_score_fast(x, y, x_len, y_len, y2x);  //y2x means the sequence has aligned 
		// three tmscore return the max score
        if(tmscore>=tmscore_max)
        {
            tmscore_max=tmscore;
            k_best=k; // what k mean?
        }
    }
    
    //extract the best map
    k=k_best;
    for(j=0; j<y_len; j++)
    {
        i=j+k;
        if(i>=0 && i<x_len)
        {
            y2x[j]=i;
        }
        else
        {
            y2x[j]=-1;
        }
    }    

    return tmscore_max;
}




//1->coil, 2->helix, 3->turn, 4->strand




//get initial alignment from secondary structure alignment
//input: x, y, x_len, y_len
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1


// get_initial5 in TMalign
//get initial alignment of local structure superposition
//input: x, y, x_len, y_len
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
bool get_initial_local(  double **x, 
						 double **y, 
						 int x_len,
						 int y_len, 
						 int *y2x
						 )
{
    double GL, rmsd;    
    double t[3];
    double u[3][3];

	double d01=d0+1.5;
	if(d01 < D0_MIN) d01=D0_MIN;
	double d02=d01*d01;
	
	double GLmax=0;
	int n_frag=20; //length of fragment for superposition
	int ns=20; //tail length to discard
	int *invmap=new int[y_len+1];
      
	int aL=getmin(x_len, y_len);
	if(aL>250)
	{
		n_frag=50;
	}
	else if(aL>200)
	{
		n_frag=40;
	}
	else if(aL>150)
	{
		n_frag=30;
	}
	else
	{
		n_frag=20;
	}
	
	int smallest=aL/3; // I change here from aL/2 to aL/3

	if(n_frag>smallest) n_frag=smallest; 
	if(ns>smallest) ns=smallest;
    
	int m1=x_len-n_frag-ns;
	int m2=y_len-n_frag-ns;

	bool flag=false;

	for(int ii=0; ii<y_len; ii++) 
	{                
		y2x[ii]=-1;
	}

	int count=0;
	for(int i=ns-1; i<m1; i=i+n_frag) //index starts from 0, different from FORTRAN
	{
		for(int j=ns-1; j<m2; j=j+n_frag)
		{
			for(int k=0; k<n_frag; k++) //fragment in y
			{ 
				r1[k][0]=x[k+i][0];  
				r1[k][1]=x[k+i][1]; 
				r1[k][2]=x[k+i][2]; 													
                        
                r2[k][0]=y[k+j][0];  
                r2[k][1]=y[k+j][1]; 
                r2[k][2]=y[k+j][2];
			}


			Kabsch(r1, r2, n_frag, 1, &rmsd, t, u);			
			count++;

			double gap_open=0.0;			
			NWDP_TM(x, y, x_len, y_len, t, u, d02, gap_open, invmap);
			GL=get_score_fast(x, y, x_len, y_len, invmap);
			if(GL>GLmax)
			{
				GLmax=GL;
				for(int ii=0; ii<y_len; ii++) 
				{                
					y2x[ii]=invmap[ii];
				}
				flag=true;
			}
		}
	}	


	delete [] invmap;
	return flag;

}


// score_matrix_rmsd  mean without second structure information
//with invmap(i) calculate score(i,j) using RMSD rotation
void score_matrix_rmsd(  double **x, 
						 double **y, 
						 int x_len,
						 int y_len,
						 int *y2x
						 )
{
	double t[3], u[3][3];
	double rmsd, dij;
	double d01=d0+1.5;
	if(d01 < D0_MIN) d01=D0_MIN;
	double d02=d01*d01;

	double xx[3];
	int i, k=0;
	for(int j=0; j<y_len; j++)
	{
		i=y2x[j];
		if(i>=0)
		{
			r1[k][0]=x[i][0];  
			r1[k][1]=x[i][1]; 
			r1[k][2]=x[i][2];   
            
			r2[k][0]=y[j][0];  
			r2[k][1]=y[j][1]; 
			r2[k][2]=y[j][2];
			
			k++;
		}
	}
	Kabsch(r1, r2, k, 1, &rmsd, t, u);
	//do_rotation(x, xt, x_len, t, u);
	
	
	for(int ii=0; ii<x_len; ii++)
	{		
		transform(t, u, &x[ii][0], xx);
		for(int jj=0; jj<y_len; jj++)
		{
			//dij=dist(&xt[ii][0], &y[jj][0]);   
			dij=dist(xx, &y[jj][0]); 
			score[ii+1][jj+1] = 1.0/(1+dij/d02);
			//	cout << ii+1 << " " << jj+1 << " " << score[ii+1][jj+1]<< endl;
		}
	}		
}


void score_matrix_rmsd_sec(  double **x, 
							 double **y, 
							 int x_len,
							 int y_len,
							 int *y2x
							 )
{
	double t[3], u[3][3];
	double rmsd, dij;
	double d01=d0+1.5;
	if(d01 < D0_MIN) d01=D0_MIN;
	double d02=d01*d01;

	double xx[3];
	int i, k=0;
	for(int j=0; j<y_len; j++)
	{
		i=y2x[j];
		if(i>=0)
		{
			r1[k][0]=x[i][0];  
			r1[k][1]=x[i][1]; 
			r1[k][2]=x[i][2];   
            
			r2[k][0]=y[j][0];  
			r2[k][1]=y[j][1]; 
			r2[k][2]=y[j][2];
			
			k++;
		}
	}
	Kabsch(r1, r2, k, 1, &rmsd, t, u);

	
	for(int ii=0; ii<x_len; ii++)
	{		
		transform(t, u, &x[ii][0], xx);
		for(int jj=0; jj<y_len; jj++)
		{
			dij=dist(xx, &y[jj][0]); 
			if(secx[ii]==secy[jj])
			{
				score[ii+1][jj+1] = 1.0/(1+dij/d02) + 0.5;
			}
			else
			{
				score[ii+1][jj+1] = 1.0/(1+dij/d02);
			}		
		}
	}		
}


//get initial alignment from secondary structure and previous alignments
//input: x, y, x_len, y_len
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ssplus( double **x, 
						 double **y, 
						 int x_len,
						 int y_len,
						 int *y2x0,
						 int *y2x						
						 )
{

	//create score matrix for DP
	score_matrix_rmsd_sec(x, y, x_len, y_len, y2x0);
	
	double gap_open=-1.0;
	NWDP_TM(x_len, y_len, gap_open, y2x);
}


void get_initial_splus(double **x, double **y, int x_len, int y_len, int* y2x0 ,int* y2x)
{

	score_matrix_rmsd(x, y, x_len, y_len, y2x0);
	double gap_open = -1.0;
	NWDP_TM(x_len,y_len,gap_open,y2x);
}

void find_max_frag(double **x, int *resno, int len, int *start_max, int *end_max)
{
	int r_min, fra_min=4;           //minimum fragment for search
	double d;
	int start;
	int Lfr_max=0, flag;

	r_min= (int) (len*1.0/3.0); //minimum fragment, in case too small protein
	if(r_min > fra_min) r_min=fra_min;
	
	int inc=0;
	double dcu0_cut=dcu0*dcu0;;
	double dcu_cut=dcu0_cut;

	while(Lfr_max < r_min)
	{		
		Lfr_max=0;			
		int j=1;    //number of residues at nf-fragment
		start=0;
		for(int i=1; i<len; i++)
		{			
			d = dist(x[i-1], x[i]);
			flag=0;
			if(dcu_cut>dcu0_cut)
			{
				if(d<dcu_cut)
				{
					flag=1;
				}
			}
			else if(resno[i] == (resno[i-1]+1)) //necessary??
			{
				if(d<dcu_cut)
				{
					flag=1;
				}
			}

			if(flag==1)
			{
				j++;

				if(i==(len-1))
				{
					if(j > Lfr_max) 
					{
						Lfr_max=j;
						*start_max=start;
						*end_max=i;						
					}
					j=1;
				}
			}
			else
			{
				if(j>Lfr_max) 
				{
					Lfr_max=j;
					*start_max=start;
					*end_max=i-1;										
				}

				j=1;
				start=i;
			}
		}// for i;
		
		if(Lfr_max < r_min)
		{
			inc++;
			double dinc=pow(1.1, (double) inc) * dcu0;
			dcu_cut= dinc*dinc;
		}
	}//while <;	
}

//perform fragment gapless threading to find the best initial alignment
//input: x, y, x_len, y_len
//output: y2x0 stores the best alignment: e.g., 
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
double get_initial_fgt( double **x, 
						double **y, 
						int x_len,
						int y_len, 
						int *xresno,
						int *yresno,
						int *y2x
						)
{
	int fra_min=4;           //minimum fragment for search
	int fra_min1=fra_min-1;  //cutoff for shift, save time

	int xstart=0, ystart=0, xend=0, yend=0;

	find_max_frag(x, xresno, x_len,  &xstart, &xend);
	find_max_frag(y, yresno, y_len, &ystart, &yend);


	int Lx = xend-xstart+1;
	int Ly = yend-ystart+1;
	int *ifr, *y2x_;
	int L_fr=getmin(Lx, Ly);
	ifr= new int[L_fr];
	y2x_= new int[y_len+1];

	//select what piece will be used (this may araise ansysmetry, but
	//only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1
	//if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1

	if(Lx<Ly || (Lx==Ly && x_len<=y_len))
	{		
		for(int i=0; i<L_fr; i++)
		{
			ifr[i]=xstart+i;
		}
	}
	else if(Lx>Ly || (Lx==Ly && x_len>y_len))
	{		
		for(int i=0; i<L_fr; i++)
		{
			ifr[i]=ystart+i;
		}	
	}

	
	int L0=getmin(x_len, y_len); //non-redundant to get_initial1
	if(L_fr==L0)
	{
		int n1= (int)(L0*0.1); //my index starts from 0
		int n2= (int)(L0*0.89);

		int j=0;
		for(int i=n1; i<= n2; i++)
		{
			ifr[j]=ifr[i];
			j++;
		}
		L_fr=j;
	}


	//gapless threading for the extracted fragment
	double tmscore, tmscore_max=-1;

	if(Lx<Ly || (Lx==Ly && x_len<=y_len))
	{
		int L1=L_fr;
	    int min_len=getmin(L1, y_len);    
		int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment 
		if(min_ali<=fra_min1)  min_ali=fra_min1;    
		int n1, n2;
		n1 = -y_len+min_ali; 
		n2 = L1-min_ali;

		int i, j, k;
		for(k=n1; k<=n2; k++)
		{
			//get the map
			for(j=0; j<y_len; j++)
			{
				i=j+k;
				if(i>=0 && i<L1)
				{				
					y2x_[j]=ifr[i];
				}
				else
				{
					y2x_[j]=-1;
				}
			}

			//evaluate the map quickly in three iterations
			tmscore=get_score_fast(x, y, x_len, y_len, y2x_);

			if(tmscore>=tmscore_max)
			{
				tmscore_max=tmscore;
				for(j=0; j<y_len; j++)
				{
					y2x[j]=y2x_[j];
				}
			}
		}
	}
	else
	{
		int L2=L_fr;
	    int min_len=getmin(x_len, L2);    
		int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment 
		if(min_ali<=fra_min1)  min_ali=fra_min1;    
		int n1, n2;
		n1 = -L2+min_ali; 
		n2 = x_len-min_ali;

		int i, j, k;	

		for(k=n1; k<=n2; k++)
		{
			//get the map
			for(j=0; j<y_len; j++)
			{
				y2x_[j]=-1;
			}

			for(j=0; j<L2; j++)
			{
				i=j+k;
				if(i>=0 && i<x_len)
				{
					y2x_[ifr[j]]=i;
				}
			}
        
			//evaluate the map quickly in three iterations
			tmscore=get_score_fast(x, y, x_len, y_len, y2x_);
			if(tmscore>=tmscore_max)
			{
				tmscore_max=tmscore;
				for(j=0; j<y_len; j++)
				{
					y2x[j]=y2x_[j];
				}
			}
		}
	}    


	delete [] ifr;
	delete [] y2x_;
    return tmscore_max;
}





//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double DP_iter( double **x,
                double **y, 
                int x_len, 
                int y_len, 
                double t[3],
                double u[3][3],
				int invmap0[],
				int g1,
				int g2,
				int iteration_max                                   
				)
{
	// g1 = 0 , g2 = 2 , iterator_max = 30
    double gap_open[2]={-0.6, 0}; // for gap open
    double rmsd; 
    int *invmap=new int[y_len+1]; // the change imap
    
    int iteration, i, j, k;
    double tmscore, tmscore_max, tmscore_old=0;    
    int score_sum_method=8, simplify_step=40;
    tmscore_max=-1;

	//double d01=d0+1.5;
    double d02=d0*d0;
    for(int g=g1; g<g2; g++) // loop 2 0,1
    {
        for(iteration=0; iteration<iteration_max; iteration++) // loop 30
        {           
			NWDP_TM(x, y, x_len, y_len, t, u, d02, gap_open[g], invmap);
            
            k=0;
            for(j=0; j<y_len; j++) 
            {
                i=invmap[j];

                if(i>=0) //aligned
                {
                    xtm[k][0]=x[i][0];
                    xtm[k][1]=x[i][1];
                    xtm[k][2]=x[i][2];
                    
                    ytm[k][0]=y[j][0];
                    ytm[k][1]=y[j][1];
                    ytm[k][2]=y[j][2];
                    k++;
                }
            }

            tmscore=TMscore8_search(xtm, ytm, k, t, u, simplify_step, score_sum_method, &rmsd);
			// calculate the TMscroe
           
            if(tmscore>tmscore_max)
            {
                tmscore_max=tmscore;
                for(i=0; i<y_len; i++) 
                {                
                    invmap0[i]=invmap[i]; // base on the  tmscore to calculate the rmsd                                     
                }				
            }
    
            if(iteration>0)
            {
                if(fabs(tmscore_old-tmscore)<0.000001) // yi jin  bubian le jiu tuichu , yaoshi bian le  na me jiu 30 
                {     
                    break;       
                }
            }
            tmscore_old=tmscore;
        }// for iteration           
        
    }//for gapopen   different gap
    
    
    delete []invmap;
    return tmscore_max;
}

void parameter_set4rowsearch(int xlen, int ylen)
{
	D0_MIN = 0.5;
	dcu0 = 4.25;
	Lnorm = getmin(xlen,ylen);
	d0 = 5.0 ;
	d0_search = d0;
	//score_d8 = 1.5*pow(Lnorm*1.0, 0.3)+3.5;
	score_d8 = 2.5*pow(Lnorm*1.0,0.45);
}

void parameter_set4search(int xlen, int ylen)
{
	//parameter initilization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
	D0_MIN=0.5; 
	dcu0=4.25;                       //update 3.85-->4.25
 
	Lnorm=getmin(xlen, ylen);        //normaliz TMscore by this in searching
	if(Lnorm<=8)                    //update 15-->19
    {
        d0=0.168;                   //update 0.5-->0.168
    }
    else
    {
        d0=(1.032*pow((Lnorm*1.0-8.111), 1.0/3)-0.440);
        //d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    }

	D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search    


	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;

    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void parameter_set4final(double len)
{
	D0_MIN=0.5; 
 
	Lnorm=len;            //normaliz TMscore by this in searching
    //*
	if(Lnorm<=8) // for the mininus d0         
    {
        d0=0.5;          
    }
    else
    {
        d0=(1.032*pow((Lnorm*1.0-8.111), 1.0/3)-0.440);
    }
	// do was update for RNA with experiment, 2016-12-08
    
	if(d0<D0_MIN) d0=D0_MIN;   

	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  

}

void parameter_set4final(double len,double do_min)
{
	D0_MIN=do_min; 
 
	Lnorm=len;            //normaliz TMscore by this in searching
    //*
	if(Lnorm<=8) // for the mininus d0         
    {
        d0=do_min;          
    }
    else
    {
        d0=(1.032*pow((Lnorm*1.0-8.111), 1.0/3)-0.440) -0.6;
    }
	// do was update for RNA with experiment, 2016-12-08
	//*/ // this code is for TMalign
	/*
	if ( Lnorm > 100)
		d0 = 19.6490;
	else
		d0 = 2.394*pow(Lnorm*1.0,0.437);

	*/
    if(d0<D0_MIN) d0=D0_MIN;   

	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  
}


void parameter_set4scale(int len, double d_s)
{
 
	d0=d_s;          
	Lnorm=len;            //normaliz TMscore by this in searching

	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  

}

double TMSCORE(double t[3], double u[3][3], int n_ali, double * &distance){
	double xrot[3];
	//cout << "AA" << endl;
	double d02 = d0*d0;
	double di,tmscore(0);
	distance = new double[n_ali];
	if (!distance)
		cout << "alocating memory fail in 1712 " << endl;
	for ( int k =0 ; k < n_ali ; k ++){
		transform(t , u ,& xtm[k][0] , xrot);
		di  = dist(xrot,&ytm[k][0]);
		distance [k] = di;
		tmscore += 1/(1 + di /d02);
	}
	delete [] distance;
	return tmscore;
}


void output_alignment(string xname,
					 string yname,
					 int x_len,
					 int y_len,
					 double t[3],
					 double u[3][3],
					 double TM1,
					 double TM2,
					 double rmsd,
					 double d0_out,
					 int m1[], 
					 int m2[],
					 int n_ali8,
					 int n_ali,
					 double TM_0,
					 double Lnorm_0,
					 double d0_0,
					 int * & ano,
					 int * & bno,
					 bool a_opt,
					 bool u_opt,
					 bool d_opt,
					 bool o_opt,
					 char superposed_file[100]
					 )
{
    double seq_id;          
    int i, j, k;
    double d;
    int ali_len=x_len+y_len; //maximum length of alignment
	char *seqM, *seqxA, *seqyA;
	seqM=new char[ali_len];
	seqxA=new char[ali_len];
	seqyA=new char[ali_len];
	do_rotation(xa, xt, x_len, t, u);
	seq_id=0;
	int kk=0, i_old=0, j_old=0;
	//print_1Array(seqy,y_len); cout << endl;
	for(k=0; k<n_ali8; k++)
	{
		for(i=i_old; i<m1[k]; i++)
		{
			//align x to gap
			seqxA[kk]=aseq[i];
			seqyA[kk]='-';
			seqM[kk]=' ';					
			kk++;
		}

		for(j=j_old; j<m2[k]; j++)
		{
			//align y to gap
			seqxA[kk]='-';
            seqyA[kk]=bseq[j];
            seqM[kk]=' ';
            kk++;
		}

		seqxA[kk]=aseq[m1[k]];
		seqyA[kk]=bseq[m2[k]];
		if(seqxA[kk]==seqyA[kk])
		{
			seq_id++;
		}
		d=sqrt(dist(&xt[m1[k]][0], &ya[m2[k]][0]));
		if(d<d0_out)
		{
			seqM[kk]=':';
		}
		else
		{
			seqM[kk]='.';
		} 
		kk++;  
		i_old=m1[k]+1;
		j_old=m2[k]+1;
	}

	//tail
	for(i=i_old; i<x_len; i++)
	{
		//align x to gap
		seqxA[kk]=aseq[i];
		seqyA[kk]='-';
		seqM[kk]=' ';					
		kk++;
	}    
	for(j=j_old; j<y_len; j++)
	{
		//align y to gap
		seqxA[kk]='-';
		seqyA[kk]=bseq[j];
		seqM[kk]=' ';
		kk++;
	}
 
    seqxA[kk]='\0';
    seqyA[kk]='\0';
    seqM[kk]='\0';
	

	seq_id=seq_id/( n_ali8+0.00000001); //what did by TMalign, but not reasonable, it should be n_ali8    
	
	/*cout <<endl;	
	cout << " *****************************************************************************" << endl
		 << " * TM-align (Version "<< version <<"): A protein structural alignment algorithm     *" << endl
		 << " * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)       *" << endl
		 << " * Please email your comments and suggestions to Yang Zhang (zhng@umich.edu) *" << endl
		 << " *****************************************************************************" << endl;	

	*/

	cout << "*******************************************************************************" << endl
		<<  "**************  RNAalign : a RNA structural alignment algortithm **************" << endl
		<<  "*******************************************************************************" << endl;
	
	//printf("\nName of structure A: %s\n", yname); 
	cout << "Name of structure A : " << yname << endl;
	cout << "Name of structure B : " << xname << " (to be superposed onto structure A)" << endl;
	//printf("Name of structure B: %s (to be superimposed onto structure A)\n", xname); 
	printf("Length of structure A: %d residues\n", y_len);
	printf("Length of structure B: %d residues\n\n", x_len);

	printf("Aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n\n", n_ali8, rmsd, seq_id); 
	
	printf("RMscore=%6.5f (if normalized by length of structure B, i.e., LN=%d, d0=%.2f)\n", TM2, x_len, d0B);
	printf("RMscore=%6.5f (if normalized by length of structure A, i.e., LN=%d, d0=%.2f)\n", TM1, y_len, d0A);
	if(a_opt)
	{
		double L_ave=(x_len+y_len)*0.5;
		printf("RMscore=%6.5f (if normalized by average length of two structures, i.e., LN=%.2f, d0=%.2f)\n", TM3, L_ave, d0a);
	}
	if(u_opt)
	{		
		printf("RMscore=%6.5f (if normalized by user-specified LN=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
	}
	if(d_opt)
	{		
		printf("RMscore=%6.5f (if scaled by user-specified d0=%.2f, and LN=%.2f)\n", TM5, d0_scale, Lnorm_0);
	}
	printf("(You should use RMscore normalized by length of the reference RNA)\n");

   

    printf("\n----- The rotation matrix to rotate Structure B to Structure A -----\n");
    
    printf("i\t%18s %15s %15s %15s\n", "t[i]", "u[i][0]", "u[i][1]", "u[i][2]");
    for(k=0; k<3; k++)
    {
        printf("%d\t%18.10f %15.10f %15.10f %15.10f\n",\
                k, t[k], u[k][0], u[k][1] ,u[k][2]);
    }
    printf("\nCode for rotating Structure B from (x,y,z) to (X,Y,Z):\n");
    printf("for(k=0; k<L; k++)\n");
    printf("{\n");
    printf("   X[k] = t[0] + u[0][0]*x[k] + u[0][1]*y[k] + u[0][2]*z[k]\n");
    printf("   Y[k] = t[1] + u[1][0]*x[k] + u[1][1]*y[k] + u[1][2]*z[k]\n");
    printf("   Z[k] = t[2] + u[2][0]*x[k] + u[2][1]*y[k] + u[2][2]*z[k]\n");    
    printf("}\n"); 
    
    
    //output structure alignment
    printf("\n(\":\" denotes residue pairs of d < %4.1f Angstrom, ", d0_out);
    printf("\".\" denotes other aligned residues)\n");
	printf("%s\n", seqxA);
    printf("%s\n", seqM);
    printf("%s\n", seqyA);

	// output residue number to ensure the actual alignment 
	//cout << endl <<" ----begin to output residue number of alignment---- " << endl;
	//print_1Array(ano,x_len);
	//print_1Array(bno,y_len);
	cout << endl;
    // the following code set output file parameter
	if (o_opt)
	{
		FILE *fp = fopen(superposed_file,"w");
		if (!fp)
		{
			cout << " open superposed file fail , programm exit." << endl;
			exit(1);
		}
		fprintf(fp , "REMARK RNAalign Version %s \n", version.c_str());
		fprintf(fp , "REMARK Structure A:%s  Length = %d \n" , yname.c_str() ,y_len);
		fprintf(fp , "REMARK Structure B:%s  Length = %d \n" , xname.c_str() ,x_len);
		fprintf(fp , "RRMARK RMRscore normalized by Length = %.2f, d0 = %.2f \n", Lnorm_0,d0_0);
		fprintf(fp , "REMARK Aligned Length = %d , RMSD = %.2f , RMscore = %.5f , ID = %.3f \n" , n_ali8, rmsd , TM_0 , seq_id);
		//fclose(fp);
		// next code need to been extended by superpose all the structure instead of only one type atom 
		// with the transport matrix and original coordinate
		ifstream in(xname.c_str());
		if (!in)
			cout << "Cann't open file to write the superposed structure!" << endl;
		string line;double x,y,z;
		while(getline(in,line))
		{
			if (line.size() >= 6) // for "ATOM"
			{
				if (line.substr(0,4) == "ATOM" || line.substr(0,6) == "HETATM")
				{
					double tx,ty,tz;
					get_xyz(line,x,y,z);
					tx = t[0] + u[0][0] * x + u[0][1] * y + u[0][2] * z;
					ty = t[1] + u[1][0] * x + u[1][1] * y + u[1][2] * z;
					tz = t[2] + u[2][0] * x + u[2][1] * y + u[2][2] * z;
					//fprintf(fp,"%s\n",line.c_str());
					fprintf(fp,"%s%8.3f%8.3f%8.3f%s\n",line.substr(0,30).c_str(),tx,ty,tz,
							line.substr(54,line.size()-54).c_str());
				}
			}
		}
		/*
		for (int i=0; i < x_len; i ++)
		{
			fprintf(fp, "ATOM  %5d  %3s %3c %5d    %8.3f%8.3f%8.3f\n",
					i+1,item.c_str(),aseq[i],xresno[i],xt[i][0],xt[i][1],xt[i][2]);
		}
		*/
		fclose(fp);
	}




/*	if(o_opt)
	{
		output_superpose(xname, yname, x_len, y_len, t, u, rmsd, d0_out, m1, m2, n_ali8, seq_id, TM_0, Lnorm_0, d0_0);
	}

*/
	delete [] seqM;
	delete [] seqxA;
	delete [] seqyA;
    
}

///
//  internal loop,stem,hairpin,other
int read_dssra(string name,string chain_id)
{

	int state(1);
	state = access( (name +".x3dna-dssr").c_str(),F_OK);
	if (state == -1)
	{
		cout << "current dir don't exits RNA second structure information,please check dssr-3dna" << endl;
		exit(0);
	}

	ifstream in((name + ".3dna-dssr").c_str());

	if (!in)
	{
		cout << "open file: " << name << ".3dna-dssr fail" << endl;
		exit(0);
	}

	string tag =  ">" + name + chain_id ; 

	string temp;
	while (getline(in,temp))
	{
		if (temp.substr(0,5) == tag){
			string seq,ss_state;
			getline(in,seq);
			getline(in,ss_state);
			if ( xlen != seq.size())
				cout << "Warning " << name << " may miss backbone " << item << endl;
			for (int i =0;i < xlen; i++)
			{
				if (ss_state[i] == '.')
					secx[i] = 0;
				if (ss_state[i] == '(' || ss_state[i] == ')')
					secx[i] = 1;
			}
			break;
		}
	}
	return 0;
}

int read_dssrb(string name,string chain_id)
{

	int state(1);
	state = access( (name +".x3dna-dssr").c_str(),F_OK);
	if (state == -1)
	{
		cout << "current dir don't exits RNA second structure information,please check dssr-3dna" << endl;
		exit(0);
	}

	ifstream in((name + ".3dna-dssr").c_str());

	if (!in)
	{
		cout << "open file: " << name << ".3dna-dssr fail" << endl;
		exit(0);
	}

	string tag =  ">" + name + "-" + chain_id ; 

	string temp;
	while (getline(in,temp))
	{
		if (temp.substr(0,5) == tag){
			string seq,ss_state;
			getline(in,seq);
			getline(in,ss_state);
			if ( ylen != seq.size())
				cout << "Warning " << name << " may miss backbone " << item << endl;
			for (int i =0;i < ylen; i++)
			{
				if (ss_state[i] == '.')
					secy[i] = 0;
				if (ss_state[i] == '(' || ss_state[i] == ')')
					secy[i] = 1;
			}
			break;
		}

	}
	return 0;
}

int read_a_second_structual(string name,string chainid)
{
	int state(1);
	state = access( (name +".x3dna-dssr").c_str(),F_OK);
	if (state == -1)
	{
		cout << "current dir don't exits RNA second structure information,please check dssr-3dna" << endl;
		exit(0);
	}

	ifstream in((name + ".x3dna-dssr").c_str());

	if (!in)
	{
		cout << "open file: " << name << ".x3dna-dssr fail" << endl;
		exit(0);
	}

	string list("List of");
	//string tag_array[3]={"internal","hairpin","stem"};
	string temp;
	for(int i=0;i<xlen;i++)
		secx[i] = 0; // for single-stranded
	while(getline(in,temp)) // deal x3dna out file 
	{
		if (temp.substr(0,7) == list) // extract string "List of "
		{
			stringstream os(temp.c_str());
			int number(0);
			string field;
			os >> field >> field >> number >> field; // second structural element number
			string special_file;
			if (number == 1)
				special_file = "stem";
			if (number >= 2)
				special_file = "stems";
			string bugles;
			if (number == 1)
				bugles = "bulge";
			if (number >= 2)
				bugles = "bulges";
			if (field == bugles)
			{
				for(int j=1;j<=number;j++)
				{
					getline(in,temp);
					getline(in,temp);
					stringstream ss;
					ss << temp;
					string frag1,frag2;
					ss >> frag1 >> frag2 ;
					int loop_time = frag2.size();
					ss >> frag2;
					for (int k=0;k<(loop_time-1);k++)
					{
						int index = frag2.find(',');
						string cell = frag2.substr(0,(index+1));
						char f2(frag2.at(0));
						ss.str("");
						ss.clear();
						ss << cell.substr(3,cell.size()-3);
						frag2.erase(0,index+1); // erase the first field divieded by comma
						ss >> index;
						for(int m=0;m<xlen;m++) // assign second structure information
						{
							if (xresno[m] == index && chidx[m] == f2)
							{
								secx[m] =4;
								break;
							}
						}
					}
					int index(0);
					string cell = frag2;
					char f2(frag2.at(0));
					ss.str("");
					ss.clear();
					ss << frag2.substr(3,cell.size()-3);
					ss >> index;
					for(int m=0;m<xlen;m++) // assign second structure information
					{
						if (xresno[m] == index && chidx[m] == f2)
						{
							secx[m] = 4; // 3 for hairpin loop
							break;
						}
					}
					getline(in,temp);
					getline(in,temp); // the rest lines
				}
			}
			if (field == special_file) // for tag stem
			{
				for(int j=1;j<=number;j++)// for stem number
				{ 
					stringstream ss;
					string number_tag;
					ss << j ;
					ss >> number_tag;
					string tag("stem#" + number_tag);
					int stem_number(0);
					while(getline(in,temp))
					{
						if (temp.find(tag) != temp.npos)
						{
							stringstream ss1;
							ss1 << temp;
							string kk;
							ss1 >> kk >> kk;
							kk = kk.substr(4,kk.length()-4);
							ss1.clear();
							ss1.str("");
							ss1  <<  kk;
							ss1 >> stem_number ; // for residue number in one stem
							break;
						}
					}
					getline(in,temp);
					getline(in,temp);
					getline(in,temp);
					getline(in,temp);
					for(int k=0;k<stem_number;k++) // for the resiude 
					{
						getline(in,temp);
						stringstream ss1;
						ss1 << temp;
						string frag1,frag2; // need to deal information
						ss1 >> frag1 >> frag1 >> frag2;
						int size1(frag1.length()),size2(frag2.length());
						int res_no1(0),res_no2(0);
						string kk(frag1.substr(3,(size1-3)));
						char f1,f2;
						f1 = frag1.at(0);
						f2 = frag2.at(0);
						ss1.str("");
						ss1 << kk; // if the residue is not a standard , bug
						ss1 >> res_no1;
						ss1.clear();
						ss1 << frag2.substr(3,(size2-3)); // if the residue is not a standard , bug
						ss1 >> res_no2;
						for(int n=0;n<xlen;n++) // for assigning the second  information
						{
							if (xresno[n] == res_no1 && chidx[n] == f1)
							{
								secx[n] = 1; // 1 for stem
							}
							if (xresno[n] == res_no2 && chidx[n] == f2)
							{
								secx[n] = 1; // 1 for stem
							}
						}
					}
				}
			}
			if (field == "internal")
			{
				for (int j=1;j<=number;j++) // loop the number of internal loop 
				{
					stringstream ss;
					getline(in,temp); // line1
					getline(in,temp); // line2 
					getline(in,temp); // line3
					ss.clear();
					ss <<  temp;
					string frat1,frat2;
					ss >> frat1 >> frat2 >> frat2;
					ss.clear();
					ss << frat1.substr(4,frat1.length()-4);
					int loop_time(0);
					ss >> loop_time;
					for(int n=0;n<(loop_time-1);n++)
					{
						int index = frat2.find(',');
						char f2(frat2.at(0));
						string cell = frat2.substr(0,(index+1));
						ss.clear();
						ss.str("");
						ss << cell.substr(3,cell.size()-3);
						frat2.erase(0,index+1); // erase the first field divieded by comma
						ss >> index;
						for(int m=0;m<xlen;m++) // assign second structure information
						{
							if (xresno[m] == index && chidx[m] == f2)
							{
								secx[m] =2;
								break;
							}
						}
					}
					int index(0);
					string cell = frat2;
					char f2(frat2.at(0));
					ss.str("");
					ss.clear();
					ss << frat2.substr(3,cell.size()-3);
					ss >> index;
					for(int m=0;m<xlen;m++) // assign second structure information
					{
						if (xresno[m] == index && chidx[m] == f2)
						{
							secx[m] = 2; // 3 for hairpin loop
							break;
						}
					}
					// for another size of internal loop 
					getline(in,temp); // line4
					ss.str("");
					ss.clear();
					ss << temp;
					ss >> frat1 >> frat2 >> frat2;
					ss.str("");
					ss.clear();
					ss << frat1.substr(4,frat1.length()-4);
					ss >> loop_time;
					for(int n=0;n<(loop_time-1);n++) // loop number - number ',' = 1
					{
						int index = frat2.find(',');
						string cell = frat2.substr(0,(index+1));
						char f2(frat2.at(0));
						ss.str("");
						ss.clear();
						ss << cell.substr(3,cell.size()-3);
						frat2.erase(0,index+1); // erase the first field divieded by comma
						ss >> index;
						for(int m=0;m<xlen;m++) // assign second structure information
						{
							if (xresno[m] == index && chidx[m] == f2)
							{
								secx[m] = 2; // 2 for internal loop
								break;
							}
						}
					}
					cell = frat2;
					f2 = frat2.at(0);
					ss.str("");
					ss.clear();
					ss << frat2.substr(3,cell.size()-3);
					ss >> index;
					for(int m=0;m<xlen;m++) // assign second structure information
					{
						if (xresno[m] == index && chidx[m] == f2)
						{
							secx[m] = 2; // 3 for hairpin loop
							break;
						}
					}

				}

			}
			if (field == "hairpin")
			{
				for(int j=1;j<=number;j++) // the total number  hairpin loop 
				{
					stringstream ss;
					string number_tag;
					getline(in,temp);
					getline(in,temp);
					getline(in,temp);
					ss.clear();
					ss << temp;
					string frat1,frat2;
					ss >> frat1 >> frat2 >> frat2;
					int loop_time(0);
					ss.str("");ss.clear();
					ss << frat1.substr(4,(frat1.size()-4));
					ss >> loop_time;
					for(int k=0;k<(loop_time-1);k++)
					{
						int index = frat2.find(',');
						char f2(frat2.at(0));
						string cell = frat2.substr(0,(index+1));
						ss.str("");
						ss.clear();
						ss << cell.substr(3,cell.size()-3);
						frat2.erase(0,index+1); // erase the first field divieded by comma
						ss >> index;
						for(int m=0;m<xlen;m++) // assign second structure information
						{
							if (xresno[m] == index && chidx[m] == f2)
							{
								secx[m] = 3; // 3 for hairpin loop
								break;
							}
						}
					}

					// for the last hairloop loop
					int index(0);
					string cell = frat2;
					char f2(frat2.at(0));
					ss.str("");
					ss.clear();
					ss << frat2.substr(3,cell.size()-3);
					ss >> index;
					for(int m=0;m<xlen;m++) // assign second structure information
					{
						if (xresno[m] == index && chidx[m] == f2)
						{
							secx[m] = 3; // 3 for hairpin loop
							break;
						}
					}
				}
			}
		}
	}
	return 0;
}

int read_b_second_structual(string name,string chainid)
{
	int state(1);
	state = access( (name +".x3dna-dssr").c_str(),F_OK);
	if (state == -1)
	{
		cout << "current dir don't exits RNA second structure information,please check dssr-3dna" << endl;
		exit(0);
	}

	ifstream in((name + ".x3dna-dssr").c_str());

	if (!in)
	{
		cout << "open file: " << name << ".x3dna-dssr fail" << endl;
		exit(0);
	}

	string list("List of");
	//string tag_array[3]={"internal","hairpin","stem"};
	string temp;
	for(int i=0;i<ylen;i++)
		secy[i] = 0; // for single-stranded
	while(getline(in,temp)) // deal x3dna out file 
	{
		if (temp.substr(0,7) == list) // extract string "List of "
		{
			stringstream os(temp.c_str());
			int number(0);
			string field;
			os >> field >> field >> number >> field; // second structural element number
			string special_file;
			if (number == 0)
				special_file = "stem";
			if (number >= 1)
				special_file = "stems";
			string bugles;
			if (number == 1)
				bugles = "bulge";
			if (number >= 2)
				bugles = "bulges";
			if (field == bugles)
			{
				for(int j=1;j<=number;j++)
				{
					getline(in,temp);
					getline(in,temp);
					stringstream ss;
					ss << temp;
					string frag1,frat2;
					ss >> frag1 >> frat2 ;
					int loop_time = frat2.size();
					ss >> frat2;
					for (int k=0;k<(loop_time-1);k++)
					{
						int index = frat2.find(',');
						string cell = frat2.substr(0,(index+1));
						char f2(frat2.at(0));
						ss.str("");
						ss.clear();
						ss << cell.substr(3,cell.size()-3);
						frat2.erase(0,index+1); // erase the first field divieded by comma
						ss >> index;
						for(int m=0;m<ylen;m++) // assign second structure information
						{
							if (yresno[m] == index && chidy[m] == f2)
							{
								secy[m] =4;  // 4 for bugles
								break;
							}
						}
					}
					// for the last bulge
					int index(0);
					string cell = frat2;
					char f2(frat2.at(0));
					ss.str("");
					ss.clear();
					ss << frat2.substr(3,cell.size()-3);
					ss >> index;
					for(int m=0;m<ylen;m++) // assign second structure information
					{
						if (yresno[m] == index && chidy[m] == f2)
						{
							secy[m] = 4; // 3 for hairpin loop
							break;
						}
					}
					getline(in,temp);
					getline(in,temp); // the rest lines
				}
			}
			if (field == special_file) // for tag stem
			{
				for(int j=1;j<=number;j++)// for stem number
				{ 
					stringstream ss;
					string number_tag;
					ss << j ;
					ss >> number_tag;
					string tag("stem#" + number_tag);
					int stem_number(0);
					while(getline(in,temp))
					{
						if (temp.find(tag) != temp.npos)
						{
							stringstream ss1;
							ss1 << temp;
							string kk;
							ss1 >> kk >> kk;
							kk = kk.substr(4,kk.length()-4);
							ss1.clear();
							ss1.str("");
							ss1  <<  kk;
							ss1 >> stem_number ; // for residue number in one stem
							break;
						}
					}
					getline(in,temp);
					getline(in,temp);
					getline(in,temp);
					getline(in,temp);
					for(int k=0;k<stem_number;k++) // for the resiude 
					{
						getline(in,temp);
						stringstream ss1;
						ss1 << temp;
						string frag1,frag2; // need to deal information
						ss1 >> frag1 >> frag1 >> frag2;
						int size1(frag1.length()),size2(frag2.length());
						int res_no1(0),res_no2(0);
						string kk(frag1.substr(3,(size1-3)));
						char f1,f2;
						f1 = frag1.at(0);
						f2 = frag2.at(0);
						ss1.str("");
						ss1 << kk; // if the residue is not a standard , bug
						ss1 >> res_no1;
						ss1.clear();
						ss1 << frag2.substr(3,(size2-3)); // if the residue is not a standard , bug
						ss1 >> res_no2;
						for(int n=0;n<ylen;n++) // for assigning the second  information
						{
							if (yresno[n] == res_no1 && chidy[n] == f1)
							{
								secy[n] = 1; // 1 for stem
							}
							if (yresno[n] == res_no2 && chidy[n] == f2)
							{
								secy[n] = 1; // 1 for stem
							}
						}
					}
				}
			}
			if (field == "internal")
			{
				for (int j=1;j<=number;j++) // loop the number of internal loop 
				{
					stringstream ss;
					getline(in,temp); // line1
					getline(in,temp); // line2 
					getline(in,temp); // line3
					ss.clear();
					ss <<  temp;
					string frat1,frat2;
					ss >> frat1 >> frat2 >> frat2;
					ss.clear();
					ss << frat1.substr(4,frat1.length()-4);
					int loop_time(0);
					ss >> loop_time;
					for(int n=0;n<(loop_time-1);n++)
					{
						int index = frat2.find(',');
						char f2(frat2.at(0));
						string cell = frat2.substr(0,(index+1));
						ss.clear();
						ss.str("");
						ss << cell.substr(3,cell.size()-3);
						frat2.erase(0,index+1); // erase the first field divieded by comma
						ss >> index;
						for(int m=0;m<ylen;m++) // assign second structure information
						{
							if (yresno[m] == index && chidy[m] == f2)
							{
								secy[m] =2;
								break;
							}
						}
					}
					int index(0);
					string cell = frat2;
					char f2(frat2.at(0));
					ss.str("");
					ss.clear();
					ss << frat2.substr(3,cell.size()-3);
					ss >> index;
					for(int m=0;m<ylen;m++) // assign second structure information
					{
						if (yresno[m] == index && chidy[m] == f2)
						{
							secy[m] = 2; // 3 for hairpin loop
							break;
						}
					}
					// for another size of internal loop 
					getline(in,temp); // line4
					ss.str("");
					ss.clear();
					ss << temp;
					ss >> frat1 >> frat2 >> frat2;
					ss.str("");
					ss.clear();
					ss << frat1.substr(4,frat1.length()-4);
					ss >> loop_time;
					for(int n=0;n<(loop_time-1);n++) // loop number - number ',' = 1
					{
						int index = frat2.find(',');
						string cell = frat2.substr(0,(index+1));
						char f2(frat2.at(0));
						ss.str("");
						ss.clear();
						ss << cell.substr(3,cell.size()-3);
						frat2.erase(0,index+1); // erase the first field divieded by comma
						ss >> index;
						for(int m=0;m<ylen;m++) // assign second structure information
						{
							if (yresno[m] == index && chidy[m] == f2)
							{
								secy[m] = 2; // 2 for internal loop
								break;
							}
						}
					}
					cell = frat2;
					f2 = frat2.at(0);
					ss.str("");
					ss.clear();
					ss << frat2.substr(3,cell.size()-3);
					ss >> index;
					for(int m=0;m<ylen;m++) // assign second structure information
					{
						if (yresno[m] == index && chidy[m] == f2)
						{
							secy[m] = 2; // 3 for hairpin loop
							break;
						}
					}

				}

			}
			if (field == "hairpin")
			{
				for(int j=1;j<=number;j++) // the total number  hairpin loop 
				{
					stringstream ss;
					string number_tag;
					getline(in,temp);
					getline(in,temp);
					getline(in,temp);
					ss.clear();
					ss << temp;
					string frat1,frat2;
					ss >> frat1 >> frat2 >> frat2;
					int loop_time(0);
					ss.str("");ss.clear();
					ss << frat1.substr(4,(frat1.size()-4));
					ss >> loop_time;
					for(int k=0;k<(loop_time-1);k++)
					{
						int index = frat2.find(',');
						char f2(frat2.at(0));
						string cell = frat2.substr(0,(index+1));
						ss.str("");
						ss.clear();
						ss << cell.substr(3,cell.size()-3);
						frat2.erase(0,index+1); // erase the first field divieded by comma
						ss >> index;
						for(int m=0;m<ylen;m++) // assign second structure information
						{
							if (yresno[m] == index && chidy[m] == f2)
							{
								secy[m] = 3; // 3 for hairpin loop
								break;
							}
						}
					}

					// for the last hairloop loop
					int index(0);
					string cell = frat2;
					char f2(frat2.at(0));
					ss.str("");
					ss.clear();
					ss << frat2.substr(3,cell.size()-3);
					ss >> index;
					for(int m=0;m<ylen;m++) // assign second structure information
					{
						if (yresno[m] == index && chidy[m] == f2)
						{
							secy[m] = 3; // 3 for hairpin loop
							break;
						}
					}
				}
			}
		}
	}
	return 0;
}


void free_momery()
{
	DeleteArray(&path, xlen+1);
	DeleteArray(&val, ylen+1);
	DeleteArray(&score, xlen+1);
	DeleteArray(&xa, xlen);
	DeleteArray(&xt, xlen);
	DeleteArray(&ya, ylen);
	DeleteArray(&r1, minlen);
	DeleteArray(&r2, minlen);
	DeleteArray(&xtm, minlen);
	DeleteArray(&ytm, minlen);
	system("rm -f dssr-*");
	delete [] secx;
	delete [] secy;
}
