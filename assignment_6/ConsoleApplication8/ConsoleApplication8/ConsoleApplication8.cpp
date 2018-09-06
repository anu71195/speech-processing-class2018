// ConsoleApplication8.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#define nl cout<<endl;
using namespace std;

double abs_max(double a, double b)//return the value of absolute max of given parameters;
{
	if (abs(a) > abs(b))return abs(a);
	return abs(b);
}
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
void store_values(vector<vector <double> >input,string filename)
{
	ofstream ofs;
	ofs.open(filename);
	for (int i = 0; i < input.size(); i++)
	{
		for (int j = 0; j < input[i].size(); j++)
		{
			ofs << to_string(input[i][j]) + " ";
		}
		ofs << endl;
	}

}
vector<double>hamming_window(vector<double>signals)
{
	double pi = 3.14159265358979323846;
	for (int i = 0; i < signals.size(); i++)
	{
		signals[i] = signals[i] * (0.54-0.46*cos(2*pi*(i/(signals.size()-1))));
	}
	return signals;
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
vector<double> get_ais(vector<double>signals, vector<double> Ri,long long int p)//phi(k)=R(k)=summationo_over_all_samples(signals(n)*s(n+k)) where n is the sample number
{
	vector<double>ai(p + 1, 0), E(p + 1, 0), k, bi;
	double value;
	
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
double get_gain_square(vector<double> Ri,vector<double>  ai,int p)
{
	double G2=Ri[0];
	for (int i = 1; i <= p; i++)
	{
//		G2 -= ai[1] * Ri[1];  //???????????????????????
	}
	return G2;
	
}
vector <double> get_cis(vector<double> ai, vector<double> Ri, double G2,int p)
{
	vector<double>ci(p+1,0);
	double value;
	ci[0] = log(G2);
	cout << " c0 is " << ci[0] << endl;
	for (int i = 1; i <= p; i++)
	{
		value = 0;
		for (int j = 1; j < i; j++)
		{
			value += ((double)j / (double)i)*ci[j] * ai[i - j];
		}
		ci[i] = ai[i] + value;
	}
	cout << endl;
	return ci;

}
vector<double> dc_shift_normalize_vac(vector<double> data)//dc shift then voice activity detection and then normalize
{
	vector<double> store_temp;
	double a = 0, b = 0, counta = 0, countb = 0, dc_shift_bound = 1000, threshold, first_high, last_high, max_v = INT_MIN, min_v = INT_MAX, amplitude=10;
	for (int i = 0; i < data.size(); i++)//finding the ambient sound for positive side and negative of the signals only at the beginning and ending of the signal (dc_shift_bound number of samples on both sides)
	{
		if (data[i] < 0 && (i <= dc_shift_bound || i >= (data.size() - dc_shift_bound)))
		{
			a += data[i];
			counta++;
		}
		else if (data[i]>0 && (i <= dc_shift_bound || i >= (data.size() - dc_shift_bound)))
		{
			b += data[i];
			countb++;
		}
	}
	a = a / counta;//average ambient sound for positive and negative sides
	b = b / countb;
	threshold = 2 * abs_max(a, b);//defining threshold with the help of average ambient sound for positive and negative sides
	for (int i = 0; i < data.size(); i++)//getting the first index of the input signal where the value is greater than threshold
	{
		if (data[i]>threshold)
		{
			first_high = i;
			break;
		}
	}
	for (int i = data.size() - 1; i >= 0; i--)//getting the last index of the input signal where the value is greater than threshold
	{
		if (data[i] > threshold)
		{
			last_high = i;
			break;
		}
	}
	//the first and last index where the values are greater than threshold are stored in first_high and last_high variable;
	for (int i = first_high; i <= last_high; i++)//removing ambient sound i.e. taking values between first high and last high  and dc shift i.e. subtracting average ambient sound on positive and negative side of the signals and finally storing in store_temp
	{
		if (data[i]>0)		store_temp.push_back(data[i] - a);
		else store_temp.push_back(data[i] - b);
		if (store_temp[store_temp.size() - 1] > max_v)max_v = store_temp[store_temp.size() - 1];
		if (store_temp[store_temp.size() - 1] < min_v)min_v = store_temp[store_temp.size() - 1];
	}
	for (int i = 0; i < store_temp.size(); i++)//normalize between -A and +A amplitude
	{
		if (store_temp[i]>0)store_temp[i] = (amplitude*((double)store_temp[i] / (double)max_v));
		else store_temp[i] =  (-amplitude*((double)store_temp[i] / (double)min_v));
	}
	data = store_temp;
	return data;
}
int _tmain(int argc, _TCHAR* argv[])
{
	ifstream infile;
	ofstream ofs;
	long long int item,p=12,num_samples=320;
	double G2;
	string folder_path = "vowel_data",filepath,s_filename;//store filename=s_filename
	vector <string> vowels = { "a", "e", "i", "o", "u" };//all the vowels stored in the vector vowels
	vector <double> data,ai,ci,Ri,all_data;
	vector<vector <double> >all_ais,all_cis,all_Ris;
	for (long long int i = 0; i < 1; i++)
	{
		for (long long int j = 0; j < vowels.size(); j++)
		{
			filepath = folder_path + "/150101010_" + vowels[j] + "_" + to_string(i + 1) + ".txt";
			infile = open_file(filepath);
			cout << filepath << endl;
			while (!infile.eof())//storing input file in store vector
			{

				infile >> item;
				data.push_back(item);

			}
			all_data = dc_shift_normalize_vac(data);
			for (int k = 0; k < 5; k++)
			{
				data.assign(all_data.begin() + k * 25 * num_samples / 100, all_data.begin() + k * 25 * num_samples / 100+num_samples);
				data = hamming_window(data);
				Ri = get_Ris(data, p);//gets the Ris from the signals and store in the vector Ri;
				for (int i = 0; i < Ri.size(); i++)cout << "Ri:" << Ri[i] << "\t";
				cout << endl;
				ai = get_ais(data, Ri, p);
				G2 = get_gain_square(Ri, ai, p);
				ci = get_cis(ai, Ri, G2, p);
				all_ais.push_back(ai);
				all_cis.push_back(ci);
				all_Ris.push_back(Ri);

				data.clear();
			}
			s_filename = "Ai_" + vowels[j] + "_.txt";
			store_values(all_ais, s_filename);
			s_filename = "Ci_" + vowels[j] + "_.txt";
			store_values(all_cis, s_filename);
			s_filename = "Ri_" + vowels[j] + "_.txt";
			store_values(all_Ris, s_filename);
			all_data.clear();
			all_ais.clear();
			all_cis.clear();
		}
		
	}


	return 0;
}

