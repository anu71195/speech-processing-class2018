// ConsoleApplication5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#define nl cout<<endl;
using namespace std;
ifstream open_file(string filepath)
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
int abs_max(int a, int b)
{
	if (abs(a) > abs(b))return abs(a);
	return abs(b);
}
vector <long long int> find_zcr(string filename)
{
	ifstream inFile;
	int item, a = 0, b = 0, counta = 0, countb = 0, threshold, first_high, last_high, max = 0, min = 0, amplitude = 5000, zcr = 0, samples = 1000, current_zcr = 0,dc_shift_bound=1000;
	vector<int> store, store_temp;
	inFile = open_file(filename);
	while (!inFile.eof())
	{
		inFile >> item;
		store.push_back(item);
	}
	for (int i = 0; i < store.size(); i++)
	{
		if (store[i] < 0 &&(i<=dc_shift_bound || i>=(store.size()-dc_shift_bound)))
		{
			a += store[i];
			counta++;
		}
		else if (store[i]>0 &&  (i <= dc_shift_bound || i >= (store.size() - dc_shift_bound)))
		{
			b += store[i];
			countb++;
		}
	}
	a = a / counta;
	b = b / countb;
	
	threshold = 2 * abs_max(a, b);
	for (int i = 0; i < store.size(); i++)
	{
		if (store[i]>threshold)
		{
			first_high = i;
			break;
		}
	}
	for (int i = store.size() - 1; i >= 0; i--)
	{
		if (store[i] > threshold)
		{
			last_high = i;
			break;
		}
	}
	for (int i = first_high; i <= last_high; i++)//removing ambient sound and dc shift
	{
		if (store[i]>0)		store_temp.push_back(store[i] - a);
		else store_temp.push_back(store[i] - b);
		if (store_temp[store_temp.size() - 1] > max)max = store_temp[store_temp.size() - 1];
		if (store_temp[store_temp.size() - 1] < min)min = store_temp[store_temp.size() - 1];
	}
	for (int i = 0; i < store_temp.size(); i++)//normalize
	{
		if (store_temp[i]>0)store_temp[i] = (int)(amplitude*((float)store_temp[i] / (float)max));
		else store_temp[i] = (int)(-amplitude*((float)store_temp[i] / (float)min));
	}

	int frequency = 0;
	long long int energy = 0, total_energy = 0;
	for (int j = 1; j < store_temp.size(); j++)
	{
		energy += store_temp[j] * store_temp[j];
		if (store_temp[j] * store_temp[j - 1] <= 0) current_zcr++;
		if (j%samples == 0)
		{
			frequency++;
			zcr += current_zcr;
			current_zcr = 0;
			total_energy += energy;
			energy = 0;
		}
	}
	vector<long long int> store_return;
	store_return.push_back(zcr / frequency);
	store_return.push_back(total_energy / frequency);
	return store_return;
	
}
int _tmain(int argc, _TCHAR* argv[])
{
	long long int zcr1 = 108, zcr6 = 202, zcr,energy;
	vector <long long int>en_zcr;
	string filename = "data/150101010_";
	for (int i = 0; i<10; i++)
	{
		en_zcr = find_zcr(filename + "1_" + to_string(i + 1) + ".txt");
		cout << en_zcr[0] << " " << en_zcr[1] << endl;
	}
	nl;
	for (int i = 0; i<10; i++)
	{
		en_zcr = find_zcr(filename + "6_" + to_string(i + 1) + ".txt");
		cout << en_zcr[0] << " " << en_zcr[1] << endl;
	}
	zcr = en_zcr[0], energy = en_zcr[1];
	if (abs(zcr - zcr1)>abs(zcr - zcr6))cout << "the sound is of 6"<<endl;
	else cout << "the sound is of 1"<<endl;
}