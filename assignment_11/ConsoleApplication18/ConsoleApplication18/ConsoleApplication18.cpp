// ConsoleApplication8.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#define nl cout<<endl;
#define ll long long int
#define lld long double
using namespace std;

void print_vector(vector<lld> input)
{
	for (int i = 0; i<input.size(); i++)
	{
		cout << input[i] << " ";
	}
	cout << endl;
}
void print_matrix(vector < vector <lld> > input)
{
	for (int i = 0; i<input.size(); i++)
	{
		print_vector(input[i]);
	}
	nl
}
lld abs_max(lld a, lld b)//return the value of absolute max of given parameters;
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
void store_values(vector<vector <lld> >input, string filename)
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
vector<lld>hamming_window(vector<lld>signals)
{
	lld pi = 3.14159;//265358979323846;
	for (int i = 0; i < signals.size(); i++)		signals[i] = signals[i] * (0.54 - 0.46*cos(2 * pi*((lld)i / (signals.size() - 1))));
	return signals;
}
vector<lld> get_Ris(vector<lld>signals, ll p)//gets the Ris from the signals
{
	vector<lld> Ri;
	lld value = 0;
	for (ll i = 0; i <= p; i++)
	{
		value = 0;
		for (ll j = 0; j < signals.size() - i; j++)
		{
			value += signals[j] * signals[j + i];//check for the overflows once-----------------------------------------------------------------------------

		}
		Ri.push_back(value);
		//cout << value << endl;
	}
	return Ri;

}
vector<lld> get_ais(vector<lld>signals, vector<lld> Ri, ll p)//phi(k)=R(k)=summationo_over_all_samples(signals(n)*s(n+k)) where n is the sample number
{
	vector<lld>ai(p + 1, 0), E(p + 1, 0), k, bi, output_ai;
	lld value;

	E[0] = Ri[0];
	k.push_back(-1);//k[0] is invalid
	bi = ai;
	output_ai.push_back(0);
	for (ll i = 1; i <= p; i++)
	{
		value = 0;
		for (ll j = 1; j < i; j++)
		{
			value += ai[j] * Ri[i - j];
		}
		value = (Ri[i] - value) / E[i - 1];
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
lld get_gain_square(vector<lld> Ri, vector<lld>  ai, ll p)
{
	double G2 = Ri[0];
	//	for (int i = 1; i <= p; i++)
	//	{
	//		G2 -= ai[1] * Ri[1];  //???????????????????????
	//}
	return G2;

}
vector <lld> get_cis(vector<lld> ai, vector<lld> Ri, lld G2, ll p)
{
	vector<lld>ci(p + 1, 0);
	lld value;
	ci[0] = log(G2);
	//cout << " c0 is " << ci[0] << endl;
	//print_vector(ai);
	//print_vector(ci);
	for (int i = 1; i <= p; i++)
	{
		value = 0;
		for (int j = 1; j < i; j++)
		{
			value += ((lld)j / (lld)i)*ci[j] * ai[i - j];
			//cout << value << endl;
		}
		ci[i] = ai[i] + value;
		//print_vector(ci);
	}
	//cout << endl;
	return ci;

}
vector<lld> dc_shift_normalize_vac(vector<lld> data)//dc shift then voice activity detection and then normalize
{
	
	vector<lld> store_temp;
	lld a = 0, b = 0, counta = 0, countb = 0, dc_shift_bound = 1000, threshold, first_high, last_high, max_v = INT_MIN, min_v = INT_MAX, amplitude = 10000, dc_shift_value = 0, count = 0;
	
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
		if (i <= dc_shift_bound || i >= (data.size() - dc_shift_bound))
		{
			dc_shift_value += data[i];
			count++;
		}
	}
	
	
	a = a / counta;//average ambient sound for positive and negative sides
	b = b / countb;
	//a = dc_shift_value / count;
	//b = a;
	threshold = 2*max(a,b);//defining threshold with the help of average ambient sound for positive and negative sides
	
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
		if (store_temp[i]>0)store_temp[i] = (amplitude*((lld)store_temp[i] / (lld)max_v));
		else store_temp[i] = (-amplitude*((lld)store_temp[i] / (lld)min_v));
	}
	
	data = store_temp;
	return data;
}
vector<vector <lld> > get_input_cis(string filepath, lld num_samples, lld p)//get the input from the file and calculate the ci's which is the input to the program
{
	lld item, G2;
	vector<lld> data, all_data, ai, Ri, ci;
	vector<vector <lld> >all_ais, all_cis, all_Ris;
	ifstream infile;
	infile = open_file(filepath);
	cout << "INPUT FILE: " << filepath << endl;
	while (!infile.eof())//storing input file in store vector
	{

		infile >> item;
		data.push_back(item);
	}
	all_data = data;
	all_data = dc_shift_normalize_vac(data);
	for (int k = 0; k < 5; k++)
	{

		data.assign(all_data.begin() + k * 25 * num_samples / 100, all_data.begin() + k * 25 * num_samples / 100 + num_samples);//taking 320 samples ==num_samples by doing 25% shift each time

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
	store_values(all_ais, "input_ais.txt");//storing values into the files
	store_values(all_cis, "input_cis.txt");
	store_values(all_Ris, "input_ris.txt");
	return all_cis;
}
lld tokura_distance(vector<vector <lld> > a, vector<vector <lld> > b, vector <lld> weights)//finding tokhura's distance
{
	lld distance = 0;
	vector <lld> vec_distance;
	for (int i = 0; i < a.size(); i++)
	{
		distance = 0;
		for (int j = 1; j < a[i].size(); j++)
		{
			distance += weights[j - 1] * (a[i][j] - b[i][j])*(a[i][j] - b[i][j]);//finding distance;
		}
		vec_distance.push_back(distance);
	}
	distance = 0;
	for (int i = 0; i < vec_distance.size(); i++)
	{
		distance += vec_distance[i];//finding the average distance for the 5 frames
	}
	distance = distance / ((lld)vec_distance.size());
	return distance;
}
int find_row_with_min_dist(vector < vector <lld> > input)//findingin the index with the least distance where index tells the vowels number
{
	int row = 0;
	lld distance = DBL_MAX;
	for (int i = 0; i < input.size(); i++)
	{
		for (int j = 0; j < input[i].size(); j++)
		{
			if (input[i][j] < distance)
			{
				distance = input[i][j];//conidition if distance found is less then update hte  min_distance and the corresponding index
				row = i;
			}
		}
	}
	return row;
}

int _tmain(int argc, _TCHAR* argv[])
{
	ifstream infile;
	ofstream ofs;
	int iterations = 20, recog_index,nframes=5;
	lld p = 12, num_samples = 320, distance;
	lld G2, item;
	string folder_path = "data", filepath, s_filename;//store filename=s_filename
	vector <lld> data, ai, ci, Ri, all_data, tokura_weights = { 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61 }, temp, digits = {0,1,2,3,4,5,6,7,8,9};
	vector<vector <lld> >all_ais, all_cis, all_Ris, input_cis,universe;

	for (ll i = 0; i < iterations; i++)
	{
		for (ll j = 0; j < digits.size(); j++)
		{
			
			filepath =folder_path+"/150101010_" + to_string((int)digits[j]) + "_" + to_string(i + 1) + ".txt";
			cout << filepath << endl;
			infile = open_file(filepath);
			
			while (!infile.eof())//storing input file in store vector
			{

				infile >> item;
				data.push_back(item);

			}
			data = dc_shift_normalize_vac(data);//does dc shift, normalize the data and do voice activity detecting
			all_data = data;
			for (ll k = 0; k < nframes; k++)
			{
				
				data.assign(all_data.begin() + k * 25 * num_samples / 100, all_data.begin() + k * 25 * num_samples / 100 + num_samples);//takes num_samples at a time and then shift 25% subsequently while taking data
				data = hamming_window(data);//apply hamming window on the data
				Ri = get_Ris(data, p);//gets the Ris from the signals and store in the vector Ri;
				ai = get_ais(data, Ri, p);
				G2 = get_gain_square(Ri, ai, p);
				ci = get_cis(ai, Ri, G2, p);
				universe.push_back(ci);
				//all_ais.push_back(ai);
				//all_cis.push_back(ci);
				//all_Ris.push_back(Ri);
				data.clear();
			}
			
			//all_data.clear();
			//all_ais.clear();
			//all_cis.clear();
		}
	}
	cout << universe.size() << endl;

	return 0;
}
