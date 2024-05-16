#include"needle.h"
#include<fstream>
#include "sstream"

#include "unistd.h"
//#include "cio"
int output_fasta(string& name,string& fasta,string chain)
{
	if (!access((name+chain+".fasta").c_str(),0))
		return 0;
	ofstream out((name+chain+".fasta").c_str());
	out << ">" << name << endl;
	out << fasta[0];
	for(int i=1;i<fasta.size();i++)
		if (!(i%60))
			out << endl << fasta[i];
		else
			out << fasta[i];
	return 0;
}

string needle_align(string& name1,string& name2,string chain1,string chain2)
{
	ostringstream os;
	string outfile(name1+chain1+name2+chain2+".align");
	os << "needle -asequence " << name1+chain1 << ".fasta -bsequence " << name2+chain2 << ".fasta -gapopen 10.0 -gapextend 2.0 "
		<< " -outfile " << outfile << " 2>/dev/null";
	if (!access(outfile.c_str(),0))
		return outfile;
	system(os.str().c_str());
	return outfile;
}

void read_needle_align(int* align,string& result,string& name1,string& name2)
{
	ifstream in(result.c_str());
	string temp;
	int length(0);
	while(getline(in,temp))
		if (temp.compare(0,9,"# Length:") == 0)
			break;
	istringstream istring;
	istring.str(temp);
	istring >> temp >> temp >> length;
	while(getline(in,temp))
		if (temp.compare(0,11,"# Identity:") == 0)
			break;
	string skip;
	istring.clear();
	istring.str(temp);
	istring >> skip >> skip >> skip;
	int posi;
	posi = skip.find('/');
	double identity1,identity2,identity;
	string temp1(skip.substr(0,posi));
	string temp2(skip.substr(posi+1,skip.size()));
	identity1 = atof(temp1.c_str());identity2 = atof(temp2.c_str());
	identity  = identity1/identity2;
	string tag1(name1),tag2(name2);
	int middle_signal =0;
	char up_char;
	string up_line,middle_line,down_line;
	if (tag1 != tag2)
	{
		while(getline(in,temp))
		{
			if(temp.compare(0,tag1.size(),tag1) == 0 )
			{
				istring.clear();
				istring.str(temp);
				istring >> skip >> skip >> skip;
				for(int i=0;i<skip.size();i++)
				{
					up_char = skip.at(i);
					up_line.push_back(up_char);
				}
				middle_signal = 1;
				continue;
			}
			if (middle_signal ==1)
			{
				skip = temp;
				for(int i=0;i<skip.size();i++)
				{
					up_char = skip.at(i);
					if (up_char == '.' || up_char == '|' || up_char == ':')
						middle_line.push_back(up_char);
				}
				middle_signal = 0;
			}
			if(temp.compare(0,tag2.size(),tag2) == 0 ){
				istring.clear();
				istring.str(temp);;
				istring >> skip >> skip >> skip ;;
				for (int i=0;i<skip.size();i++){
					up_char = skip.at(i);
					down_line.push_back(up_char);
				}
				continue;
			}
		}
	}
	if (tag1 == tag2)
	{
		int mark_again(0);
		int push_ready(0);
		while(getline(in,temp))
		{
			if (temp.compare(0,tag1.size(),tag1) == 0 );
			{
				mark_again ++ ;
				if ( mark_again % 2 )
					push_ready ++ ;
				istring.clear();
				istring.str(temp);
				istring >> skip >> skip >> skip ;
				if ( push_ready == 1){
					push_ready -- ;
					for (int i=0;i<skip.size();i++){
						up_char = skip.at(i);
						up_line.push_back(up_char);
					}
					middle_signal = 1;
					continue;
				}
				if ( middle_signal == 1 ){
					istring.clear();
					istring.str(temp);
					istring >> skip;
					for (int i=0;i<skip.size();i++){
						up_char = skip.at(i);
						if ( up_char == '.' || up_char == '|' || up_char == ':')
							middle_line.push_back(up_char);
					}
					middle_signal =0;
				}
			}
		}
		down_line = up_line ;
	}
	if ( down_line.size() != up_line.size() )
		cout << "error happen when read sequence from the alignment file " << endl;
	if ( !middle_line.size() )
		cout << "there is not aligned residuce at the sequecne" << endl;
	
	in.close();
	//system(("rm "+result).c_str());
	//cout << up_line << endl << down_line << endl;
	int up_gap(0),down_gap(0);
	for(int i=0;i<down_line.size();i++)
	{
		if(down_line.at(i) != '-' && up_line.at(i) !=  '-')
		{
			//cout << down_line.at(i) << " " << up_line.at(i) << " " << i-down_gap << " " << i - up_gap << endl;
			align[i - down_gap] = i - up_gap;
		}
		else
		{
			align[i - down_gap] = -1;
			if (down_line.at(i) == '-')
				down_gap ++;
			if (up_line.at(i) == '-')
				up_gap ++ ;
		}
	}
	//cout << "upline: " <<up_line << up_gap << endl << "downline" << down_line << down_gap << endl;
}

float read_needle_align_return(int* align,string& result,string& name1,string& name2)
{
	ifstream in(result.c_str());
	string temp;
	int length(0);
	while(getline(in,temp))
		if (temp.compare(0,9,"# Length:") == 0)
			break;
	istringstream istring;
	istring.str(temp);
	istring >> temp >> temp >> length;
	while(getline(in,temp))
		if (temp.compare(0,11,"# Identity:") == 0)
			break;
	string skip;
	istring.clear();
	istring.str(temp);
	istring >> skip >> skip >> skip;
	int posi;
	posi = skip.find('/');
	double identity1,identity2,identity;
	string temp1(skip.substr(0,posi));
	string temp2(skip.substr(posi+1,skip.size()));
	identity1 = atof(temp1.c_str());identity2 = atof(temp2.c_str());
	identity  = identity1/identity2;
	string tag1(name1),tag2(name2);
	int middle_signal =0;
	char up_char;
	string up_line,middle_line,down_line;
	if (tag1 != tag2)
	{
		while(getline(in,temp))
		{
			if(temp.compare(0,tag1.size(),tag1) == 0 )
			{
				istring.clear();
				istring.str(temp);
				istring >> skip >> skip >> skip;
				for(int i=0;i<skip.size();i++)
				{
					up_char = skip.at(i);
					up_line.push_back(up_char);
				}
				middle_signal = 1;
				continue;
			}
			if (middle_signal ==1)
			{
				skip = temp;
				for(int i=0;i<skip.size();i++)
				{
					up_char = skip.at(i);
					if (up_char == '.' || up_char == '|' || up_char == ':')
						middle_line.push_back(up_char);
				}
				middle_signal = 0;
			}
			if(temp.compare(0,tag2.size(),tag2) == 0 ){
				istring.clear();
				istring.str(temp);;
				istring >> skip >> skip >> skip ;;
				for (int i=0;i<skip.size();i++){
					up_char = skip.at(i);
					down_line.push_back(up_char);
				}
				continue;
			}
		}
	}
	if (tag1 == tag2)
	{
		int mark_again(0);
		int push_ready(0);
		while(getline(in,temp))
		{
			if (temp.compare(0,tag1.size(),tag1) == 0 );
			{
				mark_again ++ ;
				if ( mark_again % 2 )
					push_ready ++ ;
				istring.clear();
				istring.str(temp);
				istring >> skip >> skip >> skip ;
				if ( push_ready == 1){
					push_ready -- ;
					for (int i=0;i<skip.size();i++){
						up_char = skip.at(i);
						up_line.push_back(up_char);
					}
					middle_signal = 1;
					continue;
				}
				if ( middle_signal == 1 ){
					istring.clear();
					istring.str(temp);
					istring >> skip;
					for (int i=0;i<skip.size();i++){
						up_char = skip.at(i);
						if ( up_char == '.' || up_char == '|' || up_char == ':')
							middle_line.push_back(up_char);
					}
					middle_signal =0;
				}
			}
		}
		down_line = up_line ;
	}
	if ( down_line.size() != up_line.size() )
		cout << "error happen when read sequence from the alignment file " << endl;
	if ( !middle_line.size() )
		cout << "there is not aligned residuce at the sequecne" << endl;
	
	in.close();
	//system(("rm "+result).c_str());
	//cout << up_line << endl << down_line << endl;
	int up_gap(0),down_gap(0);
	for(int i=0;i<down_line.size();i++)
	{
		if(down_line.at(i) != '-' && up_line.at(i) !=  '-')
		{
			//cout << down_line.at(i) << " " << up_line.at(i) << " " << i-down_gap << " " << i - up_gap << endl;
			align[i - down_gap] = i - up_gap;
		}
		else
		{
			align[i - down_gap] = -1;
			if (down_line.at(i) == '-')
				down_gap ++;
			if (up_line.at(i) == '-')
				up_gap ++ ;
		}
	}
	return identity;
}
