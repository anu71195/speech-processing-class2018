// ConsoleApplication5.cpp : Defines the entry point for the console application.
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
int abs_max(int a, int b)//return the value of absolute max of given parameters;
{
	if (abs(a) > abs(b))return abs(a);
	return abs(b);
}
vector <long long int> find_zcr(string filename)
{
	ifstream inFile;
	int item, a = 0, b = 0, counta = 0, countb = 0, threshold, first_high, last_high, max = 0, min = 0, amplitude = 5000, zcr = 0, samples = 1000, current_zcr = 0, dc_shift_bound = 1000;
	vector<int> store, store_temp;
	inFile = open_file(filename);
	while (!inFile.eof())//storing input file in store vector
	{
		inFile >> item;
		store.push_back(item);
	}
	for (int i = 0; i < store.size(); i++)//finding the ambient sound for positive side and negative of the signals only at the beginning and ending of the signal (dc_shift_bound number of samples on both sides)
	{
		if (store[i] < 0 && (i <= dc_shift_bound || i >= (store.size() - dc_shift_bound)))
		{
			a += store[i];
			counta++;
		}
		else if (store[i]>0 && (i <= dc_shift_bound || i >= (store.size() - dc_shift_bound)))
		{
			b += store[i];
			countb++;
		}
	}
	a = a / counta;//average ambient sound for positive and negative sides
	b = b / countb;
	threshold = 2 * abs_max(a, b);//defining threshold with the help of average ambient sound for positive and negative sides
	for (int i = 0; i < store.size(); i++)//getting the first index of the input signal where the value is greater than threshold
	{
		if (store[i]>threshold)
		{
			first_high = i;
			break;
		}
	}
	for (int i = store.size() - 1; i >= 0; i--)//getting the last index of the input signal where the value is greater than threshold
	{
		if (store[i] > threshold)
		{
			last_high = i;
			break;
		}
	}
	//the first and last index where the values are greater than threshold are stored in first_high and last_high variable;
	for (int i = first_high; i <= last_high; i++)//removing ambient sound i.e. taking values between first high and last high  and dc shift i.e. subtracting average ambient sound on positive and negative side of the signals and finally storing in store_temp
	{
		if (store[i]>0)		store_temp.push_back(store[i] - a);
		else store_temp.push_back(store[i] - b);
		if (store_temp[store_temp.size() - 1] > max)max = store_temp[store_temp.size() - 1];
		if (store_temp[store_temp.size() - 1] < min)min = store_temp[store_temp.size() - 1];
	}
	for (int i = 0; i < store_temp.size(); i++)//normalize between -A and +A amplitude
	{
		if (store_temp[i]>0)store_temp[i] = (int)(amplitude*((float)store_temp[i] / (float)max));
		else store_temp[i] = (int)(-amplitude*((float)store_temp[i] / (float)min));
	}
	int frequency = 0;
	long long int energy = 0, total_energy = 0;
	for (int j = 1; j < store_temp.size(); j++)//finding the total energy and total zcr
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
	store_return.push_back(zcr / frequency);//storing the average zcr and average energy of each frame in store//frame is taken of length samples
	store_return.push_back(total_energy / frequency);
	return store_return;

}
vector<long long int> train(string filename)//train the system to find the average zcr1 and zcr 6
{
	
	vector<long long int> zcr_output;
	long long int zcr1 = 0, zcr6 = 0, zcr, energy;
	vector <long long int>en_zcr;
	for (int i = 0; i<10; i++)//finding  zcr for 1
	{
		en_zcr = find_zcr(filename + "1_" + to_string(i + 1) + ".txt");
		zcr1 += en_zcr[0]; 
	}
	nl;
	for (int i = 0; i<10; i++)//finding  zcr of 6
	{
		en_zcr = find_zcr(filename + "6_" + to_string(i + 1) + ".txt");
		zcr6 += en_zcr[0];
	}
	zcr1 = zcr1 / 10;//finding average zcr for 1 and 6 finding dividing it with number of zcrs calculated i.e. number of files
	zcr6 = zcr6 / 10;
	zcr = en_zcr[0], energy = en_zcr[1];
	zcr_output.push_back(zcr1);
	zcr_output.push_back(zcr6);
	return zcr_output;
}
void record_sound()
{
	cout << "NOTE:- It is expected that input at the beginning and ending has ambient sound of about 1 second" << endl << endl;
	system("Recording_Module.exe 3 out.wav out.txt");
}
int _tmain(int argc, _TCHAR* argv[])
{
	record_sound();//function to record the sound of the person
	cout << "training to find zcr for 1 and 6..."<<endl;
	cout << "files for training are in the location \"./data/\""<<endl;
	string filename = "data/150101010_",input_file="out.txt";
	ofstream ofs;
	ifstream inFile;
	vector<long long int> en_zcr,vec_zcr;
	long long int zcr,zcr1,zcr6;
	vec_zcr=train(filename);//traing the system to find the average value of zcr1 and zcr 6
	ofs.open("zcr.txt");
	ofs << to_string(vec_zcr[0])+"\n"+to_string(vec_zcr[1]);//writing the average values of zcr1 and zcr6 to the file zcr.txt in the current directory
	zcr1 = vec_zcr[0];
	zcr6 = vec_zcr[1];
	cout << "average zcr1 = " << zcr1 << endl;
	cout << "average zcr6 = " << zcr6 << endl<<endl;
	cout << "the same above value for average zcr of 1 and 6 are stored in the file zcr.txt in the same directory" << endl;
	en_zcr = find_zcr(input_file);//finding zcr of the given input file
	zcr = en_zcr[0];
	cout << "zcr of the file given is :- " << zcr << endl;
	if (abs(zcr - vec_zcr[0])>abs(zcr - vec_zcr[1]))cout << "the sound classified is of 6" << endl<<endl;//decision writing whether it is 1 or 6
	else cout << "the sound is of 1" << endl<<endl;
}