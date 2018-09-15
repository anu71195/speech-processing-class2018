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

void print_vector(vector<long double> input)
{
	for (int i = 0; i<input.size(); i++)
	{
		cout << input[i] << " ";
	}
	cout << endl;
}
long double abs_max(long double a, long double b)//return the value of absolute max of given parameters;
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
void store_values(vector<vector <long double> >input, string filename)
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
vector<long double>hamming_window(vector<long double>signals)
{
	long double pi = 3.14159;//265358979323846;
	for (int i = 0; i < signals.size(); i++)		signals[i] = signals[i] * (0.54 - 0.46*cos(2 * pi*((long double)i / (signals.size() - 1))));
	return signals;
}
vector<long double> get_Ris(vector<long double>signals, long long int p)//gets the Ris from the signals
{
	vector<long double> Ri;
	long double value = 0;
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
vector<long double> get_ais(vector<long double>signals, vector<long double> Ri, long long int p)//phi(k)=R(k)=summationo_over_all_samples(signals(n)*s(n+k)) where n is the sample number
{
	vector<long double>ai(p + 1, 0), E(p + 1, 0), k, bi, output_ai;
	long double value;
	
	E[0]=Ri[0];
	k.push_back(-1);//k[0] is invalid
	bi = ai;
	output_ai.push_back(0);
	for (long long int i = 1; i <= p; i++)
	{
		value = 0;
		for (long long int j = 1; j < i; j++)
		{
			value += ai[j]*Ri[i - j];
		}
		value=(Ri[i]-value)/ E[i - 1];
		output_ai.push_back(value);
		k.push_back(value);
		bi[i] = value;
		for (int j = 1; j < i; j++)
		{
			bi[j] = ai[j] - value*ai[i - j];
		}
		ai = bi;
		E[i] = (1 - value*value)*E[i - 1];

	}
	//print_vector(Ri);
	//print_vector(ai);
	//for (int i = 0; i < ai.size(); i++)cout << ai[i] << " ";
//	cout << endl;
	//cout << output_ai.size() << endl;
	//return output_ai;
	return ai;
}
long double get_gain_square(vector<long double> Ri, vector<long double>  ai, long long int p)
{
	double G2=Ri[0];
//	for (int i = 1; i <= p; i++)
//	{
//		G2 -= ai[1] * Ri[1];  //???????????????????????
	//}
	return G2;
	
}
vector <long double> get_cis(vector<long double> ai, vector<long double> Ri, long double G2, long long int p)
{
	vector<long double>ci(p + 1, 0);
	long double value;
	ci[0] = log(G2);
	//cout << " c0 is " << ci[0] << endl;
	//print_vector(ai);
	//print_vector(ci);
	for (int i = 1; i <= p; i++)
	{
		value = 0;
		for (int j = 1; j < i; j++)
		{
			value += ((long double)j / (long double)i)*ci[j] * ai[i - j];
			//cout << value << endl;
		}
		ci[i] = ai[i] + value;
		//print_vector(ci);
	}
	//cout << endl;
	return ci;

}
vector<long double> dc_shift_normalize_vac(vector<long double> data)//dc shift then voice activity detection and then normalize
{
	vector<long double> store_temp;
	long double a = 0, b = 0, counta = 0, countb = 0, dc_shift_bound = 1000, threshold, first_high, last_high, max_v = INT_MIN, min_v = INT_MAX, amplitude = 5000;
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
		if (store_temp[i]>0)store_temp[i] = (amplitude*((long double)store_temp[i] / (long double)max_v));
		else store_temp[i] = (-amplitude*((long double)store_temp[i] / (long double)min_v));
	}
	data = store_temp;
	return data;
}
int _tmain(int argc, _TCHAR* argv[])
{
	ifstream infile;
	ofstream ofs;
	long long int p=12,num_samples=320;
	long double G2,item;
	string folder_path = "vowel_data",filepath,s_filename;//store filename=s_filename
	vector <string> vowels = { "a", "e", "i", "o", "u" };//all the vowels stored in the vector vowels
	vector <long double> data, ai, ci, Ri, all_data;
	vector<vector <long double> >all_ais, all_cis, all_Ris;
	for (long long int i = 0; i < 1; i++)
	{
		for ( int j = 0; j < vowels.size(); j++)
		{
			filepath = folder_path + "/150101010_" + vowels[j] + "_" + to_string(i + 1) + ".txt";
			//filepath="trimmed.txt";
			infile = open_file(filepath);
			cout << filepath << endl;
			while (!infile.eof())//storing input file in store vector
			{

				infile >> item;
				data.push_back(item);

			}
			all_data = data;
			all_data = dc_shift_normalize_vac(data);
			
			for (int k = 0; k < 5; k++)
			{
				
				data.assign(all_data.begin() + k * 25 * num_samples / 100, all_data.begin() + k * 25 * num_samples / 100 + num_samples); 
				
				data = hamming_window(data);
				Ri = get_Ris(data, p);//gets the Ris from the signals and store in the vector Ri;
				////for (int i = 0; i < Ri.size(); i++)cout << "Ri:" << Ri[i] << "\t";
				//cout << endl;
				ai = get_ais(data, Ri, p);
				G2 = get_gain_square(Ri, ai, p);
				//cout << "G2 is " << G2 << endl;
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

