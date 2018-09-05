// ConsoleApplication8.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#define nl cout<<endl;
using namespace std;

ifstream open_file(string filepath)//open the file with given filepath and returns the file handle
{
	ifstream inFile;
	inFile.open(filepath);
	if (inFile.fail())
	{
		cerr << "Error opening file" << endl;
		exit(1);
	}
	return inFile;
}
vector<double> get_Ris(vector<double>signals, long long int p)//gets the Ris from the signals
{
	vector<double> Ri;
	double value = 0;
	for (long long int i = 0; i <= p; i++)
	{
		value = 0;
		for (long long int j = 0; j < signals.size() - i; j++)
		{
			value +=signals[j] *signals[j + i];//check for the overflows once-----------------------------------------------------------------------------
			
		}
		Ri.push_back(value);
		//cout << value << endl;
	}
	return Ri;
	
}
vector<double> get_ais(vector<double>signals, long long int p)//phi(k)=R(k)=summationo_over_all_samples(signals(n)*s(n+k)) where n is the sample number
{
	vector<double>ai(p + 1, 0), Ri, E(p + 1, 0), k, bi;
	double value;
	Ri=get_Ris(signals,p);//gets the Ris from the signals and store in the vector Ri;
	E[0]=Ri[0];
	k.push_back(-1);//k[0] is invalid
	bi = ai;
	for (long long int i = 1; i <= p; i++)
	{
		value = 0;
		for (long long int j = 1; j < i; j++)
		{
			value += ai[i-1]*Ri[i - j];
		}
		value=(Ri[i]-value)/ E[i - 1];
		
		k.push_back(value);
		bi[i] = value;
		for (int j = 1; j < i; j++)
		{
			bi[j] = ai[j] - value*ai[i - j];
		}
		ai = bi;
		E[i] = (1 - value*value)*E[i - 1];

	}
	

	for (int i = 0; i < ai.size(); i++)cout << ai[i] << " ";
	cout << endl;
	
	return ai;
}

int _tmain(int argc, _TCHAR* argv[])
{
	ifstream infile;
	ofstream ofs;
	long long int item,p=12;
	string folder_path = "vowel_data",filepath;
	vector <string> vowels = { "a", "e", "i", "o", "u" };//all the vowels stored in the vector vowels
	vector <double> data;
	for (long long int i = 0; i < 1; i++)
	{
		for (long long int j = 0; j < vowels.size(); j++)
		{
			filepath = folder_path + "/150101010_" + vowels[j] + "_" + to_string(i + 1) + ".txt";
			infile=open_file(filepath);
			cout << filepath << endl;
			while (!infile.eof())//storing input file in store vector
			{
				
				infile >> item;
				data.push_back(item);
				
			}
			get_ais(data, p);
			data.clear();
		}
				
	}
	return 0;
}

