// ConsoleApplication12.cpp : Defines the entry point for the console application.
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
void print_vector(vector<double>input)//printing a single dimension vector with datatype double
{
	for (int i = 0; i < input.size(); i++)
	{
		cout << input[i] << " ";
	}
	nl
}
void print_vector(vector<int>input)//printing a single dimension vector with datatype int
{
	for (int i = 0; i < input.size(); i++)
	{
		cout << input[i] << " ";
	}
	nl
}
void print_matrix(vector <vector <double> > input)//printing a matrix with datatype double
{
	for (int i = 0; i < input.size(); i++)
	{
		print_vector(input[i]);
	}
	cout << endl;
}
double doubleRand() //generate a random number between 0 and 1
{
	return double(rand()) / (double(RAND_MAX) + 1.0);
}
vector<double> parse_line(string line, char delimit)//parse line from csv file with delimiter delimit
{
	int first = 0;
	vector < double>output;//stores the splitted numbers from line delimited by delimit
	for (int i = 0; i < line.size(); i++)
	{
		if (line[i] == delimit)
		{
			output.push_back(stod(line.substr(first, i - first)));//stores the number by converting string to double
			first = i + 1;
		}
	}
	output.push_back(stod(line.substr(first, line.size() - first)));
	return output;//return the output to the parent function
}
vector <vector <double > > get_input_csv(string filename)//get input from csv file
{
	vector<vector<double> > X;//stores all the data from filename
	ifstream ifs(filename);
	string line;
	while (ifs.good())
	{
		vector<double>temp_data;
		getline(ifs, line, '\n');//stores the line separated by \n and stores in line
		if (line.size() == 0)break;
		temp_data = parse_line(line, ',');//parse line delimited by ,
		X.push_back(temp_data);//stores it in X
	}
	return X;//returns the X to the parent function
}
vector <vector <double > > get_input_txt(string filename)//get input from txt file
{
	vector <vector <double> > X;//stores all the data
	vector <double> temp;
	ifstream infile;
	infile = open_file(filename);
	double item;
	int count = 0;
	while (!infile.eof())//storing input file in store vector
	{
		infile >> item;//take input from the file
		temp.push_back(item);//store it in temp
		count++;
		if (count == 13)
		{
			temp.erase(temp.begin());//ignore first value in the vector
			X.push_back(temp);//and store it in X so temp is of 12 dimensions in our case
			temp.clear();
			count = 0;
		}
	}
	return X;//return X to the parent function
}
double tokhura(vector<double> a, vector<double> b, vector <double> weights)//finding tokhura's distance
{
	double distance = 0;
	vector <long double> vec_distance;
	for (int i = 0; i < a.size(); i++)
	{
		distance += weights[i] * (a[i] - b[i])*(a[i] - b[i]);//finding distance using tokhura's weights;	
	}
	return distance;
}
double get_max(vector<vector<double> > X)//get max of all the values in X
{
	double max_val=0;
	for (int i = 0; i < X.size(); i++)
	{
		for (int j = 0; j < X[i].size(); j++)
		{
			max_val = max(max_val,X[i][j]);
		}
	}
	return max_val;
}
double get_min(vector<vector<double> > X)//get min of all the values in X
{
	double min_val = numeric_limits<double>::max();
	for (int i = 0; i < X.size(); i++)
	{
		for (int j = 0; j < X[i].size(); j++)
		{
			min_val = min(min_val, X[i][j]);
		}
	}
	return min_val;
}
double find_distortion(vector<vector<double> > X, vector<int>  Y, vector< vector< double> >  V, vector<double>tokura_weights)//finding the distorting in our centroids with respect to the X using tokhura's distance
{
	double distortion = 0;
	for (int i = 0; i < X.size(); i++)
	{
		distortion+=tokhura(X[i], V[Y[i]], tokura_weights);//finding distorting using tokura distance and then sum it in distortion X[i] is the vector Y[i]is its label and V[Y[i]] is the centroid it belongs to
	}
	distortion = distortion / X.size();//average distortion
	return distortion;//return distortion to the parent function
}
void store_values(vector <int> input, string filename)//store thevector in the file given by the filename
{
	ofstream ofs;
	ofs.open(filename);
	for (int i = 0; i < input.size(); i++)
	{
		ofs << to_string(input[i]) ;
		ofs << endl;
	}

}
void store_values(vector<vector < double> >input, string filename)
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
void k_means(vector<vector<double> > X, int k, vector<double>tokura_weights)
{
//	double max_val = get_max(X);
	//double min_val = get_min(X);

	vector < vector<double> > V;
	vector<int> Y(X.size(),0);
	double min_dist = numeric_limits<double>::max(), dist, distortion = numeric_limits<double>::max(), old_distortion = numeric_limits<double>::max();
	int assigned_class = -1,max_iterations=100;//max iterations 100
	for(int i=0;i<k;i++)//initalization of the centroids
	{
		V.push_back(X[(i/(k-1))*(X.size()-1)]);//initializing centroids taking some point from X and then use it as initialcentroid
	}

	//classification 
	for (int iter = 0; iter < max_iterations; iter++)//run it for max max_iterations
	{
		cout << "iterations number="<<iter << endl;
		old_distortion = distortion;//storing the old distortion in old distortion and then subsequently in the further code finding distortion new one
		for (int i = 0; i < X.size(); i++)
		{
			for (int j = 0; j < k; j++)
			{
				dist = tokhura(X[i], V[j], tokura_weights);//finding tokhura's distance
				if (min_dist>dist)
				{
					min_dist = dist;//finding the least tokhura's distance
					assigned_class = j;
				}
			}
			min_dist = numeric_limits<double>::max();
			Y[i] = assigned_class;//using the distance cluster it to one of the class which has min distance
		}

		//finding centroids;
		for (int i = 0; i < V.size(); i++)//finding centroids
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
						temp[k] += X[j][k];//sum over all points belong to every cluster one by one
					}
				}
			}
			for (int k = 0; k < X[0].size(); k++)
			{
				if (count == 0)break;
				temp[k] = temp[k] / (double)count;//taking average
			}
			if (count)V[i] = temp;//make sure the centroid changing has atleast one point in its cluster
		}
		distortion=	find_distortion(X,Y,V,tokura_weights);//finding the distorting in our centroids with respect to the X using tokhura's distance
		cout << "old distortion="<<old_distortion << endl;
		cout << "new distortion=" << distortion << endl;
		cout << "change in distortion="<<old_distortion - distortion << endl;//the difference in distortion
		if (abs(old_distortion - distortion) < 0.00000000000001)//checking if change in distortion is very low if yes  then break
		{
			break;
		}	
		nl
		
	}nl
	cout << endl;
	cout << "CENTROIDS" << endl;
	print_matrix(V);
	cout << "LABELS" << endl;
	print_vector(Y);//get all the Y values
	cout << endl;
	cout << "LABELS AND CENTROIDS ARE STORED IN THE FILE labels.txt AND centroids.txt RESPECTIVELY IN THE SAME DIRECTORY IN WHICH CPP FILE IS THERE";nl
	store_values(Y, "labels.txt");
	store_values(V,"centroids.txt");
}
int _tmain(int argc, _TCHAR* argv[])
{
	vector<vector<double> > X;
	vector<double>tokura_weights = { 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61 };//tokura's weights initialized
	//string filename = "Universe.csv";//uncomment this and comment below line to run universe.csv
	string filename = "My_Universe.txt";//datafile
	int k = 8;//number of clusters
	//X = get_input_csv(filename);//uncomment this and comment below line to run universe.csv
	X = get_input_txt(filename);//get input from the file Myuniverse.txt
	k_means(X,k,tokura_weights);//run k m eans on the dataset
	return 0;
}

