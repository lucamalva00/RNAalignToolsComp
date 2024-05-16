#include <iostream>
#include <string>
#include "Kabsch.h"
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
using namespace std;

template <class A> void NewArray(A*** array, int Narray1,int Narray2);
template <class A> void DeleteArray(A ** array, int Narray);

int main(int argc,char* argv[])
{
	string temp;
	ifstream in(argv[1]);
	if (!in)
	{
		cout << "can't can open file " << argv[1] << endl;
		exit(1);
	}
	double x(0),y(0),z(0);
	int n = atoi(argv[3]);
	double **set1,**set2;
	NewArray(&set1, n, 3);
	NewArray(&set2, n, 3);
	for (int i=0;i<n;i++)
	{
		getline(in,temp);
		stringstream stringio(temp);
		stringio >> x >> y >> z;
		//cout  << x << " " << y << " "<< z << endl;
		set1[i][0] = x;
		set1[i][1] = y;
		set1[i][2] = z;
	}
	/*
	for (int i=0;i<35;i++)
	{
		cout << set1[i][0] << " " << set1[i][1] << " " << set1[i][2] << endl;
	}
	*/
	in.close();
	ifstream in2(argv[2]);
	if (!in2)
	{
		cout << "can't open file " << argv[2] << endl;
		exit(1);
	}
	for (int i=0;i<n;i++)
	{
		getline(in2,temp);
		stringstream stringio(temp);
		stringio >> x >> y >> z;
		//cout  << x << " " << y << " "<< z << endl;
		set2[i][0] = x;
		set2[i][1] = y;
		set2[i][2] = z;
	}
	/*
	for (int i=0;i<35;i++)
	{
		cout << set2[i][0] << " " << set2[i][1] << " " << set2[i][2] << endl;
	}
	*/
	in2.close();
	int mode = 0;
	double rmsd(0);
	double *rms = &rmsd;
	double t[3],r[3][3]; 
	Kabsch(set1,set2,n,mode,rms,t,r);
	rmsd = sqrt(rmsd/n);
	cout << "rmsd: " << rmsd <<" "<< *rms << endl;
	DeleteArray(&set1, n);
	DeleteArray(&set2, n);
	return 0;
}

template <class A> void NewArray(A*** array, int Narray1,int Narray2)
{                                                            
	*array = new A* [Narray1];
	for (int i=0;i<Narray1;i++)
		*(*array + i) = new A [Narray2];
};
template <class A> void DeleteArray(A ** array, int Narray){
	for(int i=0; i<Narray; i++)
		if(*(*array+i)) delete [] *(*array+i);
	if(Narray) delete [] (*array);
	(*array)=NULL;
};  
