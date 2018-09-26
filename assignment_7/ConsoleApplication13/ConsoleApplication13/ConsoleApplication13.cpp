// ConsoleApplication13.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include  <iostream>
#include  <fstream>
#include  <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <algorithm>
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
void print_vector(vector<double>input)
{
	for (int i = 0; i < input.size(); i++)
	{
		cout << input[i] << " ";
	}
	nl
}
void print_vector(vector<int>input)
{
	for (int i = 0; i < input.size(); i++)
	{
		cout << input[i] << " ";
	}
	nl
}
void print_matrix(vector <vector <double> > input)
{
	for (int i = 0; i < input.size(); i++)
	{
		print_vector(input[i]);
	}
	cout << endl;
}
double doubleRand() {
	return double(rand()) / (double(RAND_MAX) + 1.0);
}
vector<double> parse_line(string line, char delimit)
{
	int first = 0;
	vector < double>output;
	for (int i = 0; i < line.size(); i++)
	{
		if (line[i] == delimit)
		{
			output.push_back(stod(line.substr(first, i - first)));
			first = i + 1;
		}
	}
	output.push_back(stod(line.substr(first, line.size() - first)));
	return output;
}
vector <vector <double > > get_input_csv(string filename)
{
	vector<vector<double> > X;
	ifstream ifs(filename);
	string line;
	while (ifs.good())
	{
		vector<double>temp_data;
		getline(ifs, line, '\n');
		if (line.size() == 0)break;
		temp_data = parse_line(line, ',');
		X.push_back(temp_data);
	}
	return X;
}
vector <vector <double > > get_input_txt(string filename)
{
	vector <vector <double> > X;
	vector <double> temp;
	ifstream infile;
	infile = open_file(filename);
	double item;
	int count = 0;
	while (!infile.eof())//storing input file in store vector
	{
		infile >> item;
		temp.push_back(item);
		count++;
		if (count == 13)
		{
			temp.erase(temp.begin());
			X.push_back(temp);
			temp.clear();
			count = 0;
		}
	}
	return X;
}
double tokhura(vector<double> a, vector<double> b, vector <double> weights)//finding tokhura's distance
{
	double distance = 0;
	vector <long double> vec_distance;
	for (int i = 0; i < a.size(); i++)
	{
		distance += weights[i] * (a[i] - b[i])*(a[i] - b[i]);//finding distance;	
	}
	return distance;
}
double find_distortion(vector<vector<double> > X, vector<int>  Y, vector< vector< double> >  V, vector<double>tokura_weights)
{
	double distortion = 0;
	for (int i = 0; i < X.size(); i++)
	{
		distortion += tokhura(X[i], V[Y[i]], tokura_weights);
	}
	distortion = distortion / X.size();
	return distortion;
}
vector <double> get_centroid(vector <vector <double > > X)
{
	vector <double> output(X[0].size(),0);
	for (int i = 0; i < X.size(); i++)
	{
		for (int j = 0; j < X[0].size(); j++)
		{
			output[j] += X[i][j];
		}
	}
	for (int j = 0; j < output.size(); j++)
	{
		output[j] = output[j] / (double)X.size();
	}
	return output;
}
vector<int>  k_means(vector<vector<double> > X, vector< vector<double> > &V, vector<double>tokura_weights)
{
	//	double max_val = get_max(X);
	//double min_val = get_min(X);
	vector<int> Y(X.size(), 0);
	double min_dist = numeric_limits<double>::max(), dist, distortion = numeric_limits<double>::max(), old_distortion = numeric_limits<double>::max();
	int assigned_class = -1, max_iterations = 100;

	//classification 
	for (int iter = 0; iter < max_iterations; iter++)
	{
		cout << "iterations number=" << iter << endl;
		old_distortion = distortion;
		for (int i = 0; i < X.size(); i++)
		{
			for (int j = 0; j < V.size(); j++)
			{
				dist = tokhura(X[i], V[j], tokura_weights);
				if (min_dist>dist)
				{
					min_dist = dist;
					assigned_class = j;
				}
			}
			min_dist = numeric_limits<double>::max();
			Y[i] = assigned_class;
		}

		//finding centroids;
		for (int i = 0; i < V.size(); i++)
		{
			vector<double>temp(V[0].size(), 0);
			int count = 0;
			for (int j = 0; j < X.size(); j++)
			{
				if (Y[j] == i)
				{
					count++;
					for (int k = 0; k < X[0].size(); k++)
					{
						temp[k] += X[j][k];
					}
				}
			}
			for (int k = 0; k < X[0].size(); k++)
			{
				if (count == 0)break;
				temp[k] = temp[k] / (double)count;
			}
			if (count)V[i] = temp;
		}
		distortion = find_distortion(X, Y, V, tokura_weights);
		cout << "old distortion=" << old_distortion << endl;
		cout << "new distortion=" << distortion << endl;
		cout << "change in distortion=" << old_distortion - distortion << endl;
		if (abs(old_distortion - distortion) < 0.00000000000001)
		{
			break;
		}
		nl

	}
	print_vector(Y);
	return Y;

}
vector < vector  <double> >split_codebook(vector <vector <double> > V, double epsilon)
{
	vector <vector<double> >output;
	vector <double> temp;
	for (int i = 0; i < V.size(); i++)
	{
		for (int j = 0; j < V[i].size(); j++)
		{
			temp.push_back(V[i][j] + epsilon);
		}
		output.push_back(temp);
		temp.clear();
		for (int j = 0; j < V[i].size(); j++)
		{
			temp.push_back(V[i][j] - epsilon);
		}
		output.push_back(temp);
		temp.clear();
	}
	return output;
}
void LBG(vector <vector <double >  > X,int k,vector <double> tokura_weights,double epsilon)
{
	vector< vector<double> > V;
	vector <int> Y;
	int num_clusters = 0;
	//initialization
	V.push_back(get_centroid(X));
	print_matrix(V);

	while (V.size() < k)
	{

		V = split_codebook(V, epsilon);
		print_matrix(V);
		Y = k_means(X, V, tokura_weights);
		print_matrix(V);
	}
	print_matrix(V);
	print_vector(Y);

}
int _tmain(int argc, _TCHAR* argv[])
{
	vector<vector<double> > X;
	vector<double>tokura_weights = { 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61 };
	//string filename = "My_Universe.txt";
	string filename = "Universe.csv";
	int k = 8;
	double epsilon = 0.03;
	//X = get_input_csv(filename);
	X = get_input_txt(filename);
	LBG(X,k,tokura_weights,epsilon);
	return 0;
}

