// ConsoleApplication3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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
int _tmain(int argc, _TCHAR* argv[])
{
	ifstream inFile;
	string end_line_text, filename = "data/150101010_";
	int item;
	vector <vector <int> > store;
	vector <int> store_temp;
	//file format somethingi_j.txt where i and j digit and its iteration number something part must be given in filename;see number_of_iterations variable
	int  countp = 0, countn = 0, bitspersample = 16, channels = 1, samplerate = 16000, normalized = 0, number_of_iterations = 10, ambient_sound_average = 0, a = 0, b = 0,max=0,min=INT_MAX;
	int amplitude = 5000,threshold,first_high,last_high,zcr=0,zcr1,zcr6,samples=1000;
	//number of iterations which represents j
	inFile=open_file("data/ambient_sound.txt");
	while (!inFile.eof())
	{
		inFile >> item;
		if (item>0)
		{
			a += item;
			countp++;
		}
		else
		{
			b += item;
			countn++;
		}
	}
	a = a / countp, b = b / countn;
	cout << a << " " << b << endl;//a is positive side average and b is negative side average;

	//reading files of 1 at a time and storing dc shift value in the store variable;
	for (int i = 0; i<10; i++)
	{
	//	cout << filename + "1_" + to_string(i + 1)+".txt" << endl;
		inFile = open_file(filename + "1_" + to_string(i + 1) + ".txt");
		max = 0,min=INT_MAX;
		while (!inFile.eof())
		{
			inFile >> item;
			if (item>0)item=item - a;
			else item = item - b;
			if (item > max)max = item;
			if (item < min)min = item;
			store_temp.push_back(item);
		}
		threshold = (int)(2 * ((float)abs_max(a, b)*(float)5000 / (float)abs_max(max, min)));
		cout << "threshold is" << threshold << endl;
		for (int j = 0; j < store_temp.size(); j++)//normalizing each file between -amplitude and +amplitude;
		{
			if (store_temp[j]>0)store_temp[j] = (int)(amplitude*((float)store_temp[j] /(float) max));
			else store_temp[j] = (int)(-amplitude*((float)store_temp[j]/(float)min));
		}
		for (int j = 0; j <store_temp.size(); j++)
		{
			if (abs(store_temp[j])>threshold)
			{
				first_high = j; break;
			}
		}
		for (int j = store_temp.size()-1; j >= 0; j--)
		{
			if (abs(store_temp[j])>threshold)
			{
				last_high = j; break;
			}
		}
		vector<int> update_store_temp;
		int current_zcr = 0;
		cout << "first high and last high" << endl;
		cout << first_high << endl;
		cout << last_high << endl;
		cout << store_temp.size() << endl;
		cout << threshold << endl;
		for (int j = first_high ; j < last_high ; j++)
		{
			update_store_temp.push_back(store_temp[j]);
			//if (store_temp[j] * store_temp[j - 1] <= 0)current_zcr++;
		}
		int cur_zcr = 0,frequency=0;
		for (int j = first_high; j < last_high ; j++)
		{
			if (store_temp[j] * store_temp[j - 1] <= 0)current_zcr++;
			if (j%samples == 0)
			{
				frequency++;
				cur_zcr += current_zcr;
				current_zcr = 0;
			}
		}
		zcr += cur_zcr/frequency;
		store.push_back(update_store_temp);
		cout << "zcr is"<<zcr /(i+1)<< endl;
		store_temp.clear();
	}
	zcr = zcr / 10;
	zcr1 = zcr;


	store_temp.clear();
	store.clear();
	zcr = 0;
	cout << endl;
	//reading files of 6 at a time and storing dc shift value in the store variable;
	for (int i = 0; i<10; i++)
	{
		//	cout << filename + "6_" + to_string(i + 1)+".txt" << endl;
		inFile = open_file(filename + "6_" + to_string(i + 1) + ".txt");
		max = 0, min = INT_MAX;
		while (!inFile.eof())
		{
			inFile >> item;
			if (item>0)item = item - a;
			else item = item - b;
			if (item > max)max = item;
			if (item < min)min = item;
			store_temp.push_back(item);
		}
		threshold = (int)(2 * ((float)abs_max(a, b)*(float)5000 / (float)abs_max(max, min)));
		cout << "threshold is" << threshold << endl;
		for (int j = 0; j < store_temp.size(); j++)//normalizing each file between -amplitude and +amplitude;
		{
			if (store_temp[j]>0)store_temp[j] = (int)(amplitude*((float)store_temp[j] / (float)max));
			else store_temp[j] = (int)(-amplitude*((float)store_temp[j] / (float)min));
		}
		for (int j = 0; j <store_temp.size(); j++)
		{
			if (abs(store_temp[j])>threshold)
			{
				first_high = j; break;
			}
		}
		for (int j = store_temp.size() - 1; j >= 0; j--)
		{
			if (abs(store_temp[j])>threshold)
			{
				last_high = j; break;
			}
		}
		vector<int> update_store_temp;
		int current_zcr = 0;
		cout << "first high and last high" << endl;
		cout << first_high << endl;
		cout << last_high << endl;
		cout << store_temp.size() << endl;
		cout << threshold << endl;
		for (int j = first_high; j < last_high; j++)
		{
			update_store_temp.push_back(store_temp[j]);
			
		}
		int cur_zcr = 0, frequency = 0;
		for (int j = first_high; j < last_high; j++)
		{
			if (store_temp[j] * store_temp[j - 1] <= 0)current_zcr++;
			if (j%samples == 0)
			{
				frequency++;
				cur_zcr += current_zcr;
				current_zcr = 0;
			}
		}
		zcr += cur_zcr / frequency;
		store.push_back(update_store_temp);
		cout << "zcr is" << zcr / (i + 1) << endl;
		store_temp.clear();
	}
	zcr = zcr / 10;
	zcr6 = zcr;
	cout << "zcr 1 and 6 are " << endl;
	cout << zcr1 << " " << zcr6 << endl;
	cout << a << " " << b << endl;

	return 0;
}