#ifndef NEEDLE_H_HH
#define NEEDLE_H_HH
#include <string>
#include <iostream>
#include <cstdlib>
#include "function.h"
using namespace std;
int output_fasta(string& name,string& fasta,string chain);
string needle_align(string& name1,string& name2,string chain1,string chain2);
void read_needle_align(int* align,string& result,string& name1,string& name2);
float read_needle_align_return(int* align,string& result,string& name1,string& name2);
void parameter_set4rowsearch(int xlen, int ylen);
#endif
