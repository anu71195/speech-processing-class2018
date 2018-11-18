// ConsoleApplication8.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#define nl cout<<endl;
#define ll long long int
#define lld long double

using namespace std;
vector<string> filename_vector;
int iteration_num,avg_num;

lld vector_sum(vector<lld>input)
{
	lld output = 0;
	for (ll i = 0; i < input.size(); i++)output += input[i];
	return output;
}
void print_vector(vector<lld>input)//printing a single dimension vector with datatype double
{
	for (ll i = 0; i < input.size(); i++)		cout << input[i] << " ";	nl
}
void print_vector(vector<ll>input)//printing a single dimension vector with datatype int
{
	for (ll i = 0; i < input.size(); i++)	cout << input[i] << " ";	nl
}
void print_matrix(vector<ll>input)//printing a single dimension vector with datatype int
{
	for (ll i = 0; i < input.size(); i++)	cout << input[i] << " ";	nl
}
void print_matrix(vector<int>input)//printing a single dimension vector with datatype int
{
	for (ll i = 0; i < input.size(); i++)	cout << input[i] << " ";	nl
}
void print_matrix(vector<lld>input)//printing a single dimension vector with datatype int
{
	for (ll i = 0; i < input.size(); i++)	cout << input[i] << " ";	nl
}
void print_matrix(vector <vector <lld> > input)//printing a matrix with datatype double
{
	for (ll i = 0; i < input.size(); i++)		print_vector(input[i]); nl
}
void print_matrix(vector <vector <int> > input)//printing a matrix with datatype double
{
	for (ll i = 0; i < input.size(); i++)		print_matrix(input[i]); nl
}
void print_matrix(vector < vector <vector <lld> > > input)
{
	for (ll i = 0; i < input.size(); i++)cout << i << " " << endl, print_matrix(input[i]); nl
}
vector<lld >addone(vector<lld >input)
{
	for (int i = 0; i < input.size(); i++)input[i]++;
	return input;
}

bool isfileexists(const string & filename)
{
	ifstream ifile(filename.c_str());
	return (bool)ifile;
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
vector < vector <ll> > get_matrix_file_input(string filename, ll N,ll option)//getting matrix from file filename where each row has N columns
{
	vector <vector <ll> > X;
	vector <ll> temp;
	ll item;
	ll count = 0;
	ifstream infile = open_file(filename);
	while (!infile.eof())
	{
		infile >> item;
		temp.push_back(item);
		count++;
		if (count == N)		count = 0, X.push_back(temp), temp.clear();
	}
	return X;
}
vector < vector <lld> > get_matrix_file_input(string filename, ll N)//getting matrix from file filename where each row has N columns
{
	vector <vector <lld> > X;
	vector <lld> temp;
	lld item;
	ll count = 0;
	ifstream infile = open_file(filename);
	while (!infile.eof())
	{
		infile >> item;
		temp.push_back(item);
		count++;
		if (count == N)		count = 0, X.push_back(temp), temp.clear();
	}
	return X;
}
vector <lld> get_vector_file_input(string filename)//getting vector from file
{
	vector <lld> output;
	lld item;
	ifstream infile = open_file(filename);
	struct stat buffer;

	while (!infile.eof())
	{
		infile >> item;
		output.push_back(item);
	}
	return output;
}

vector <lld> get_vector_file_input(string filename,int & flag)//getting vector from file
{
	vector <lld> output;
	lld item;
	
	
	flag = isfileexists(filename);
	if (flag == 0)return output;
	ifstream infile = open_file(filename);
	while (!infile.eof())
	{
		infile >> item;
		output.push_back(item);
	}
	return output;
}vector < vector <lld> >  get_A(string filename, ll N)//getting transition matrix from input file initial model
{
	return get_matrix_file_input(filename, N);
}
vector< vector <lld> > get_B(string filename, ll M)//getting observation probability from input file initial model
{
	return get_matrix_file_input(filename, M);
}
vector <lld> get_pi_matrix(string filename)//getting state matrix from input file initial model
{
	return get_vector_file_input(filename);
}
vector<lld> get_observation_sequence(string filename)//getting observation sequence vector from input file initial model
{
	return get_vector_file_input(filename);
}
lld abs_max(lld a, lld b)//return the value of absolute max of given parameters;
{
	if (abs(a) > abs(b))return abs(a);
	return abs(b);
}

void store_values(lld input, string filename)//store thevector in the file given by the filename
{
	ofstream ofs;
	ofs.open(filename);
	ofs << to_string(input);
	ofs << endl;


}
void store_values(vector <lld> input, string filename)//store thevector in the file given by the filename
{
	ofstream ofs;
	ofs.open(filename);
	for (ll i = 0; i < input.size(); i++)
	{
		ofs << to_string(input[i]);
		ofs << endl;
	}

}
void store_values(vector<vector < lld> >input, string filename)//STORES MATRIX IN THE FILE GIVEN BY THENAME filename
{
	ofstream ofs;
	ofs.open(filename);
	for (ll i = 0; i < input.size(); i++)
	{
		for (ll j = 0; j < input[i].size(); j++)
		{
			ofs << to_string(input[i][j]) + " ";
		}ofs << endl;
	}

}
void store_values(vector < vector < vector <lld > > > input, string filename)
{
	ofstream ofs;
	ofs.open(filename);
	for (ll i = 0; i < input.size(); i++)
	{
		for (ll j = 0; j < input[i].size(); j++)
		{
			for (ll k = 0; k < input[i][j].size(); k++)
			{
				ofs << to_string(input[i][j][k]) + " ";
			}ofs << endl;
		}ofs << endl;
	}
}
void store_values(vector<vector < ll> >input, string filename)//STORES MATRIX IN THE FILE GIVEN BY THENAME filename
{
	ofstream ofs;
	ofs.open(filename);
	for (ll i = 0; i < input.size(); i++)
	{
		for (ll j = 0; j < input[i].size(); j++)
		{
			ofs << to_string(input[i][j]) + " ";
		}ofs << endl;
	}

}
vector<lld>hamming_window(vector<lld>signals)
{
	lld pi = 3.14159;//265358979323846;
	for (ll i = 0; i < signals.size(); i++)		signals[i] = signals[i] * (0.54 - 0.46*cos(2 * pi*((lld)i / (signals.size() - 1))));
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
	lld G2 = Ri[0];
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
	for (ll i = 1; i <= p; i++)
	{
		value = 0;
		for (ll j = 1; j < i; j++)
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
	
	for (ll i = 0; i < data.size(); i++)//finding the ambient sound for positive side and negative of the signals only at the beginning and ending of the signal (dc_shift_bound number of samples on both sides)
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
	//a = 1000, b = 1000;
	//a = dc_shift_value / count;
	//b = a;
	threshold = max(a,b);//defining threshold with the help of average ambient sound for positive and negative sides
	threshold = 1000;
	
	for (ll i = 0; i < data.size(); i++)//getting the first index of the input signal where the value is greater than threshold
	{
		if (data[i]>threshold)
		{
			first_high = i;
			break;
		}
	}
	for (ll i = data.size() - 1; i >= 0; i--)//getting the last index of the input signal where the value is greater than threshold
	{
		if (data[i] > threshold)
		{
			last_high = i;
			break;
		}
	}
	
	//the first and last index where the values are greater than threshold are stored in first_high and last_high variable;
	for (ll i = first_high; i <= last_high; i++)//removing ambient sound i.e. taking values between first high and last high  and dc shift i.e. subtracting average ambient sound on positive and negative side of the signals and finally storing in store_temp
	{
		if (data[i]>0)		store_temp.push_back(data[i] - a);
		else store_temp.push_back(data[i] - b);
		if (store_temp[store_temp.size() - 1] > max_v)max_v = store_temp[store_temp.size() - 1];
		if (store_temp[store_temp.size() - 1] < min_v)min_v = store_temp[store_temp.size() - 1];
	}
	
	for (ll i = 0; i < store_temp.size(); i++)//normalize between -A and +A amplitude
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
	for (ll k = 0; k < 5; k++)
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

	return all_cis;
}
lld tokura_distance(vector<vector <lld> > a, vector<vector <lld> > b, vector <lld> weights)//finding tokhura's distance
{
	lld distance = 0;
	vector <lld> vec_distance;
	for (ll i = 0; i < a.size(); i++)
	{
		distance = 0;
		for (ll j = 1; j < a[i].size(); j++)
		{
			distance += weights[j - 1] * (a[i][j] - b[i][j])*(a[i][j] - b[i][j]);//finding distance;
		}
		vec_distance.push_back(distance);
	}
	distance = 0;
	for (ll i = 0; i < vec_distance.size(); i++)
	{
		distance += vec_distance[i];//finding the average distance for the 5 frames
	}
	distance = distance / ((lld)vec_distance.size());
	return distance;
}
ll find_row_with_min_dist(vector < vector <lld> > input)//findingin the index with the least distance where index tells the vowels number
{
	ll  row = 0;
	lld distance = DBL_MAX;
	for (ll i = 0; i < input.size(); i++)
	{
		for (ll j = 0; j < input[i].size(); j++)
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
vector < vector <lld> > create_universe(ll  &iterations, ll  &nframes, lld &p, lld &num_samples,ll offset)
{
	ifstream infile;
	ofstream ofs;
	lld G2, item;
	string folder_path = "data", filepath;//store filename=s_filename
	vector <lld> data, ai, ci, Ri, all_data, temp, digits = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	vector<vector <lld> >all_ais, all_cis, all_Ris, input_cis, universe;
	filename_vector.clear();
	for (ll i = 0+offset; i < iterations+offset; i++)
	{
		for (ll j = 0; j < digits.size(); j++)
		{

			filepath = folder_path + "/150101010_" + to_string((ll)digits[j]) + "_" + to_string(i + 1) + ".txt";
			filename_vector.push_back(filepath);
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
				//cout <<k<<" "<< all_data.size() << " " << k * 25 * num_samples / 100 + num_samples << endl;
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
	return universe;
}
lld tokhura(vector<lld> a, vector<lld> b, vector <lld> weights)//finding tokhura's distance
{
	lld distance = 0;
	vector <lld> vec_distance;
	for (ll i = 1; i < a.size(); i++)
	{
		distance += weights[i-1] * (a[i] - b[i])*(a[i] - b[i]);//finding distance using tokhura's weights;	
	}
	return distance;
}
lld find_distortion(vector<vector<lld> > X, vector<ll >  Y, vector< vector< lld> >  V, vector<lld>tokura_weights)//finding the distorting in our centroids with respect to the X using tokhura's distance
{
	lld distortion = 0;
	for (ll i = 0; i < X.size(); i++)
	{
		distortion += tokhura(X[i], V[Y[i]], tokura_weights);//finding distorting using tokura distance and then sum it in distortion X[i] is the vector Y[i]is its label and V[Y[i]] is the centroid it belongs to
	}
	distortion = distortion / X.size();//average distortion
	return distortion;//return distortion to the parent function
}
vector <lld> get_centroid(vector <vector <lld > > X)
{
	vector <lld> output(X[0].size(), 0);
	for (ll i = 0; i < X.size(); i++)
	{
		for (ll j = 0; j < X[0].size(); j++)
		{
			output[j] += X[i][j];
		}
	}
	for (ll j = 0; j < output.size(); j++)
	{
		output[j] = output[j] / (lld)X.size();
	}
	return output;
}
vector<ll >  k_means(vector<vector<lld> > X, vector< vector<lld> > &V, vector<lld>tokura_weights,lld threshold)
{
	//	lld max_val = get_max(X);
	//lld min_val = get_min(X);
	vector<ll > Y(X.size(), 0);
	lld min_dist = numeric_limits<lld>::max(), dist, distortion = numeric_limits<lld>::max(), old_distortion = numeric_limits<lld>::max();
	ll  assigned_class = -1, max_iterations = 100;//max iterations 100

	//classification 
	for (ll iter = 0; iter < max_iterations; iter++)//run it for max max_iterations
	{
		cout << "iterations number=" << iter << endl;
		old_distortion = distortion;//storing the old distortion in old distortion and then subsequently in the further code finding distortion new one
		for (ll i = 0; i < X.size(); i++)
		{
			for (ll j = 0; j < V.size(); j++)
			{
				dist = tokhura(X[i], V[j], tokura_weights);//finding tokhura's distance
				if (min_dist>dist)
				{
					min_dist = dist;//finding the least tokhura's distance
					assigned_class = j;
				}
			}
			min_dist = numeric_limits<lld>::max();
			Y[i] = assigned_class;//using the distance cluster it to one of the class which has min distance
		}

		//finding centroids;
		for (ll i = 0; i < V.size(); i++)//finding centroids
		{
			vector<lld>temp(V[0].size(), 0);
			ll  count = 0;
			for (ll j = 0; j < X.size(); j++)
			{
				if (Y[j] == i)
				{
					count++;
					for (ll k = 0; k < X[0].size(); k++)
					{
						temp[k] += X[j][k];//sum over all points belong to every cluster one by one
					}
				}
			}
			for (ll k = 0; k < X[0].size(); k++)
			{
				if (count == 0)break;
				temp[k] = temp[k] / (lld)count;//taking average
			}
			if (count)V[i] = temp;//make sure the centroid changing has atleast one point in its cluster
		}
		distortion = find_distortion(X, Y, V, tokura_weights);//finding the distorting in our centroids with respect to the X using tokhura's distance
		cout << "old distortion=" << old_distortion << endl;
		cout << "new distortion=" << distortion << endl;
		cout << "change in distortion=" << old_distortion - distortion << endl;//the difference in distortion
		if (abs(old_distortion - distortion) < threshold)//checking if change in distortion is very low if yes  then break
		{
			break;
		}
		nl

	}
	//get all the Y values
	return Y;

}
vector < vector  <lld> >split_codebook(vector <vector <lld> > V, lld epsilon)//split codebook i.e. divide the centroidinto by adding and subtracting epsilon so codebook size llds
{
	vector <vector<lld> >output;
	vector <lld> temp;
	for (ll i = 0; i < V.size(); i++)
	{
		for (ll j = 0; j < V[i].size(); j++)
		{
			temp.push_back(V[i][j] + epsilon);//get centroid  by adding epsilon
		}
		output.push_back(temp);
		temp.clear();
		for (ll j = 0; j < V[i].size(); j++)
		{
			temp.push_back(V[i][j] - epsilon);//get centroid by subtracting epsilon
		}
		output.push_back(temp);//store this in the output which is the new codebook
		temp.clear();
	}
	return output;//return new codebook
}
void LBG(vector<vector<lld>>&V, vector<vector<lld> >&Y_all, vector<ll> & Y, vector <vector <lld >  >& X, ll  &k, vector <lld>& tokura_weights, lld& epsilon, ll  nframes,lld threshold)
{
	lld min_dist = numeric_limits<lld>::max(), dist;
	vector<lld > temp;
	ll  num_clusters = 0, assigned_class = -1;
	//initialization
	V.push_back(get_centroid(X));//initialize the first centroid which is the mean of all the points
	//print_matrix(V);

	while (V.size() < k)
	{
		cout << "code book size is " << V.size() << endl;
		V = split_codebook(V, epsilon);//split codebook i.e. divide the centroidinto by adding and subtracting epsilon so codebook size llds
		Y = k_means(X, V, tokura_weights,threshold);//apply kmeans on this using this new centroids and dataset as X
	}
	cout << "Y size is " << Y.size() << endl;
	for (ll i = 0; i < X.size(); i++)
	{
		for (ll j = 0; j < V.size(); j++)
		{
			dist = tokhura(X[i], V[j], tokura_weights);//finding tokhura's distance
			if (min_dist>dist)
			{
				min_dist = dist;//finding the least tokhura's distance
				assigned_class = j;
			}
		}
		min_dist = numeric_limits<lld>::max();
		Y[i] = assigned_class;//using the distance cluster it to one of the class which has min distance
		temp.push_back(assigned_class);
		if (temp.size() == nframes)			Y_all.push_back(temp), temp.clear();
		
	}

}
void check_stochastic(vector<lld> input)
{
	lld temp = 0;
	for (ll i = 0; i < input.size(); i++)
	{
		temp += input[i];

	}
	cout << temp << endl;
}

int check_stochastic(vector<vector<lld> > input)
{
	lld temp = 0, total = 0;
	for (ll i = 0; i < input.size(); i++)
	{
		temp = 0;
		for (ll j = 0; j < input[i].size(); j++)
		{
			temp += input[i][j];
		}
		total += temp;
		if (temp == 1)cout << "yes it is 1" <<" "<<1-temp<< endl;
		else cout << "no it is not 1" << " " << 1-temp << endl;
	
	}

	if (input.size() == total)cout << "yes" << endl;
	else cout << "NO "<<total << endl;
	cout << endl;
	return 0;
}
void forward_initialization(vector<vector<lld> >&alpha, vector < lld> pi_matrix, vector <vector <lld> > B, vector<lld> O)//initialization of forward procedure
{
	vector < lld> temp;
	for (ll i = 0; i < pi_matrix.size(); i++)		temp.push_back(pi_matrix[i] * B[i][O[0] - 1]);
	alpha.push_back(temp);
}
void forward_induction(vector<vector <lld> > &alpha, vector < vector<lld> >& A, vector<vector<lld> >& B, vector <lld>& O, ll &T)//induction step in forward procedure
{
	vector<lld> alpha_temp;
	lld temp;
	for (ll t = 0; t < T - 1; t++)
	{
		alpha_temp.clear();
		for (ll j = 0; j < A.size(); j++)
		{
			temp = 0;
			for (ll i = 0; i < A.size(); i++)
			{
				temp += alpha[alpha.size() - 1][i] * A[i][j];//summation term
			}
			alpha_temp.push_back(temp*B[j][O[t + 1] - 1]);//creating vector for every time t is updated
		}
		alpha.push_back(alpha_temp);//storing vector for time t in alpha 
	}
}
void termination(lld &prob, vector<vector<lld> > & alpha)//finding probability of observation given model
{
	prob = 0;
	for (ll i = 0; i < alpha[0].size(); i++)prob += alpha[alpha.size() - 1][i];
}
void forward_procedure(vector < vector <lld> > A, vector < vector < lld> > B, vector <lld> pi_matrix, vector <lld> O, ll T, vector <vector <lld> > &alpha, lld &P_O_lambda)//forward procedure
{
	
	//initialization step
	
	alpha.clear();
	forward_initialization(alpha, pi_matrix, B, O);
	
	//induction step
	forward_induction(alpha, A, B, O, T);

	//terminationa step
	termination(P_O_lambda, alpha);

	
}
void backward_initialization(vector<vector<lld> > &beta, ll N)//initialization of backward procedure
{
	beta.push_back(vector<lld>(N, 1));//at time t =T all the states have value 1
}
int check_all_zeros(vector<vector<lld> > beta)
{
	lld total = 0;
	for (ll i = 0; i < beta.size(); i++)
	{
		for (ll j = 0; j < beta[i].size(); j++)
		{
			total += beta[i][j];
		}
	}
	if (total == 5)return 1;
	return 0;
}
void backward_induction(vector<vector<lld> > & beta, vector<vector<lld> >& A, vector < vector < lld> > &B, vector <lld> &O, ll& T)//induction step in backward procedure
{
	vector<lld> b_temp;
	vector<vector<lld> > beta_temp = beta;
	lld temp = 0;
	
	for (ll t = T - 1; t > 0; t--)
	{
		b_temp.clear();
		for (ll i = 0; i < A.size(); i++)
		{
			temp = 0;
			for (ll j = 0; j < A.size(); j++)
			{
				temp += A[i][j] * B[j][O[t] - 1] * beta[0][j];//summation term
				//cout <<i<<" "<<j<<" "<< A[i][j] << " " << " " << B[j][O[t] - 1] << " " << beta[0][j] << endl;
			}
			b_temp.push_back(temp);//storing values from summation term in vector for every time t for every state
		}
		beta.insert(beta.begin(), b_temp);//making a matrix beta by storing values from b_temp vector for each time and each state 1 to T and 1 to N
	}
	//print_matrix(beta);
	
}
void backward_termination(vector<vector<lld> > & beta, lld &P)
{
	P = 0;
	for (ll i = 0; i<beta[0].size(); i++) //this is for self assessment and not required as there is no further use
	{
		P += beta[0][i];
		
	}
}
void backward_procedure(vector < vector <lld> > A, vector < vector < lld> > B, vector <lld> O, ll T, vector<vector<lld> > &beta, lld& P_beta)//backward procedure
{

	//initialization step
	beta.clear();
	backward_initialization(beta, A.size());

	//induction step
	backward_induction(beta, A, B, O, T);
	

	
	backward_termination(beta, P_beta);
	//cout << "beta is " << endl;
	//print_matrix(beta);


}
void viterbi_initialization(vector<vector<lld> > &delta, vector< vector < lld> > & psi, vector<vector<lld> >& b, vector<lld>& pi_matrix, vector<lld>& O)//viterbi algorithm initialization
{
	vector<lld> temp;
	for (ll i = 0; i < pi_matrix.size(); i++)		temp.push_back(pi_matrix[i] * b[i][O[0] - 1]);
	delta.push_back(temp), psi.push_back(vector<lld>(pi_matrix.size(), 0));
}
void viterbi_induction(vector<vector<lld> > &delta, vector < vector < lld> >& psi, vector < vector <lld> > A, vector < vector < lld> > b, vector<lld>O, ll T)//induction step in viterbi algorithm
{
	vector<lld> delta_temp, psi_temp;
	lld temp = LONG_MIN, temp2;
	ll index;
	
	for (ll t = 1; t < T; t++)
	{
		delta_temp.clear(), psi_temp.clear();
		for (ll j = 0; j < A.size(); j++)
		{
			temp = LONG_MIN;

			for (ll i = 0; i < A.size(); i++)
			{
				temp2 = delta[delta.size() - 1][i] * A[i][j];
				//cout << temp << " " << temp2 << endl;
			
				if (temp < temp2)
				{
					temp = temp2, index = i;			//max value step and its corresponding index
				}
			}
			delta_temp.push_back(temp*b[j][O[t] - 1]);	//storring max value for each time t and state in delta_temp
			psi_temp.push_back(index);//and its corresponding index psi_temp
		}
		delta.push_back(delta_temp);//storing delta_temp vector in delta matrix
		psi.push_back(psi_temp);//storing psi_temp vector in psi matrix
	}
	//cout << psi.size() << endl;
}
void viterbi_termination(vector< vector< lld> > delta, lld &p_star, lld &q_star)//viterbi terimination step
{
	p_star = delta[delta.size() - 1][0];
	for (ll i = 0; i < delta[delta.size() - 1].size(); i++)
	{
		if (p_star < delta[delta.size() - 1][i])//getting the biggest probability from delta and storing it in p_star and its corresponding index in q_star
		{
			p_star = delta[delta.size() - 1][i];
			q_star = i;
		}
	}
}
vector<lld> viterbi_backtracking(vector<vector<lld> > psi, lld q_star)//back tracking best state sequence
{
	vector< lld > output;
	for (ll t = psi.size() - 2; t >= 0; t--)
	{
		output.insert(output.begin(), q_star);//storing q_star or state sequence values in output
		q_star = psi[t + 1][q_star];//from the q_star obtained from the last step getting q_star for previous step
	}
	output.insert(output.begin(), q_star);
	return output;
}
void viterbi(vector<vector<lld> > A, vector < vector <lld> > b, vector<lld> pi_matrix, vector<lld> O, ll T, vector<vector<lld> > &delta, vector<vector<lld> > &psi, vector<lld> &state_sequence, lld &  p_star, lld & q_star)
{

	//initialization
	delta.clear(), psi.clear();
	viterbi_initialization(delta, psi, b, pi_matrix, O);

	//induction step
	viterbi_induction(delta, psi, A, b, O, T);

	//termination step
	viterbi_termination(delta, p_star, q_star);


	//backtracking step
	state_sequence = viterbi_backtracking(psi, q_star);
	
}
lld get_polambda_at_t(vector<vector<lld> > & A, vector < vector < lld> > & B, vector < vector <lld> > & alpha, vector < vector <lld> > & beta, vector<lld> & O, ll t)
{
	lld polambda = 0;
	ll N = A[0].size();
	
	for (ll i = 0; i < N; i++)
	{

		for (ll j = 0; j < N; j++)
		{

			polambda += alpha[t][i] * A[i][j] * B[j][O[t + 1] - 1] * beta[t + 1][j];
			
		}
	}
	
	return polambda;
	
	
}
void scale(vector<vector<lld>> & A)
{
	for (int i = 0; i < A.size(); i++)
	{
		lld sum_value = 0;
		for (int j = 0; j < A[i].size(); j++)
		{
			sum_value += A[i][j];
		}
		if (sum_value == 0)continue;
		for (int j = 0; j < A[i].size(); j++)
		{
			A[i][j] = A[i][j] / sum_value;
		}
	}
}
lld matrix_sum(vector<vector<lld> > & zhi_t)
{
	lld output = 0;
	for (ll i = 0; i < zhi_t.size(); i++)
	{
		for (ll j = 0; j < zhi_t[i].size(); j++)
		{
			output += zhi_t[i][j];
		}
	}
	return output;

}
void scale2d(vector<vector<lld> > & zhi_t)
{
	lld divisor = matrix_sum(zhi_t);
	for (ll i = 0; i < zhi_t.size(); i++)
	{
		for (ll j = 0; j < zhi_t[i].size(); j++)
		{
			zhi_t[i][j] = zhi_t[i][j] / divisor;
		}
	}
}
vector <vector <lld> > get_zhi_at_t(vector<vector<lld> > & A, vector < vector < lld> > & B, vector < vector <lld> > & alpha, vector < vector <lld> > & beta, vector<lld> & O, ll t)
{
	vector<vector<lld> > zhi_t(A[0].size(), vector<lld>(A[0].size(), 0));
	ll N = A[0].size();
	lld polambda = get_polambda_at_t(A, B, alpha, beta, O, t);
	lld MIN_VALUE = 1e-150;
	for (ll i = 0; i < N; i++)
	{

		for (ll j = 0; j < N; j++)
		{
			if (polambda)
			{
				zhi_t[i][j] = alpha[t][i] * A[i][j] * B[j][O[t + 1] - 1] * beta[t + 1][j] / polambda;
				
			}
			else
			{
				zhi_t[i][j] = 0;


			}
			if (zhi_t[i][j] >0&& zhi_t[i][j] < MIN_VALUE) zhi_t[i][j] = MIN_VALUE;
			
		}
	}
	
	//print_matrix(A);// , print_matrix(B), print_matrix(alpha), print_matrix(beta), print_matrix(zhi_t);
	scale2d(zhi_t);
	return zhi_t;

}
void findgamma(vector <vector < vector <lld> > > &zhi, vector < vector <lld> > & gamma)
{
	vector < vector < lld > > output(zhi.size(), vector<lld>(zhi[0].size(), 0));
	//cout << output.size() << endl;
	//cout << zhi.size() << endl;
	//cout << zhi[0].size() << endl;
	
	for (ll t = 0; t < zhi.size(); t++)
	{

		for (ll i = 0; i < zhi[0].size(); i++)
		{
			output[t][i] = vector_sum(zhi[t][i]);
			
		}
	}
	
	
	gamma = output;
}
void findzhi(vector <vector < vector <lld> > > &zhi, vector<vector<lld> > & A, vector < vector < lld> > & B, vector < vector <lld> > & alpha, vector < vector <lld> > & beta, vector<lld> & O, ll T)
{

	for (ll t = 0; t < T - 1; t++)
		zhi.push_back(get_zhi_at_t(A, B, alpha, beta, O, t));
	
	
	
}
void update_pi_matrix(vector<vector <lld> > & gamma, vector <lld> &pi_matrix)
{
	
	for (ll i = 0; i < pi_matrix.size(); i++)		pi_matrix[i] = gamma[0][i];

}
lld get_sum_gamma_i(vector< vector<lld> > & gamma, ll i)
{
	lld output = 0;
	for (ll t = 0; t < gamma.size(); t++)		output += gamma[t][i];
	return output;

}
lld get_sum_zhi_ij(vector<vector<vector<lld> > > & zhi, ll i, ll j)
{
	lld output = 0;
	for (ll t = 0; t < zhi.size(); t++)output += zhi[t][i][j];
	return output;

}

void update_A_matrix(vector <vector < vector <lld> > > &zhi, vector<vector<lld> > & gamma, vector < vector < lld> > & A)
{
	
	
	
	
	
	bool flag = 0;
	for (ll i = 0; i < A.size(); i++)
	{
		lld sum_gamma_ti = get_sum_gamma_i(gamma, i);
		

	
		for (ll j = 0; j < A[i].size(); j++)
		{
			
			A[i][j] = get_sum_zhi_ij(zhi, i, j);// / sum_gamma_ti;
			//if (sum_gamma_ti)A[i][j] = get_sum_zhi_ij(zhi, i, j) / sum_gamma_ti;
			//A[i][j] = 0;
	
			
			
		}
	}
	
//	print_matrix(A);

	
	scale(A);
	
}
lld get_bjk_num(vector<vector< lld> > & gamma, ll j, vector< lld> & O, ll K)
{
	lld output = 0;
	int flag = 0;
	for (ll t = 0; t < gamma.size(); t++)
	{
		if ((O[t] - 1) == K)
		{
			flag = 1;
			output += gamma[t][j];

		}
	}
	
//	if (output == 0)cout << "0 ";
//	else cout << output<<" ";
	return output;
}
lld get_bjk_den(vector<vector<lld> > & gamma, ll j)
{
	lld output = 0;
	for (ll t = 0; t < gamma.size(); t++)output += gamma[t][j];
	return output;
}
void adjust_stochastic(vector<vector<lld> > & B)
{
	lld temp = 0,max_value,max_index;
	for (ll i = 0; i < B.size(); i++)
	{
		max_value = B[i][0], max_index = 0;;
		for (ll j = 0; j < B[i].size(); j++)
		{
			temp += B[i][j];
			if (B[i][j]>max_value)
			{
				max_value = B[i][j];
				max_index = j;
			}
		}
		B[i][max_index] += 1 - temp;
		
		temp = 0;
	}
}
void update_B_matrix(vector< vector <lld> > & gamma, vector<lld> & O, vector < vector<lld> > & B)
{
	//cout << "---------------------------------------------------------------------------------------------------------------------------------" << endl;

	lld temp = 0;
	
	for (ll j = 0; j < B.size(); j++)
	{
		lld bjk_den = get_bjk_den(gamma, j);
		for (ll k = 0; k < B[j].size(); k++)
		{
			B[j][k] = get_bjk_num(gamma, j, O, k);// / bjk_den;
			
			if (B[j][k] == 0)B[j][k] = 1e-16;
			//if (bjk_den)B[j][k] = get_bjk_num(gamma, j, O, k);// / bjk_den;
			//else B[j][k] = 0;
		
				
		}
		
		

	}
	//adjust_stochastic(B
	
	
	scale(B);
	
	
	//cout << "---------------------------------------------------------------------------------------------------------------------------------" << endl;
	

}
void EM(vector <vector < vector <lld> > > &zhi, vector<vector<lld> > & A, vector < vector < lld> > & B, vector < vector <lld> > & alpha, vector < vector <lld> > & beta, vector<vector<lld> > &gamma, vector<lld> & O, ll T, vector <lld> &pi_matrix)
{
	zhi.clear(), gamma.clear();

	findzhi(zhi, A, B, alpha, beta, O, T);
	
	
	findgamma(zhi, gamma);
	


	//update_pi_matrix(gamma, pi_matrix);
	


	update_A_matrix(zhi, gamma, A);
	

	update_B_matrix(gamma, O, B);
	





}


void inner_wrapper(vector < vector <lld> >& A, vector < vector < lld> >& B, vector <lld>& pi_matrix, vector <lld> O, ll T)
{


	for (ll i = 0; i < 20; i++)
	{
		vector <vector < vector <lld> > > zhi;
		vector <vector <lld> > alpha, beta, delta, psi, gamma;
		vector<lld> state_sequence;
		lld P_O_lambda, p_star, q_star, P_beta;

		
		
		forward_procedure(A, B, pi_matrix, O, T, alpha, P_O_lambda);//running forward procedure
		
		backward_procedure(A, B, O, T, beta, P_beta);//running backward procedure

		
		
		viterbi(A, B, pi_matrix, O, T, delta, psi, state_sequence, p_star, q_star);//running viterbi algorithm		
		
		
		EM(zhi, A, B, alpha, beta, gamma, O, T, pi_matrix);
		//cout << "state sequence "<< state_sequence.size() << endl;
		//print_matrix(state_sequence);
		//cout << P_O_lambda << endl;
		
		
		

	}
}
void outer_wrapper(vector<vector<lld> > & universe, vector<vector<lld> > & V, vector<vector<lld> > & Y_all, vector<lld>&tokura_weights, vector<ll>&Y, lld &p, lld &num_samples, lld&epsilon, ll &iterations, ll &nframes, ll &k,lld &threshold,ll offset)
{
	universe = create_universe(iterations, nframes, p, num_samples,offset);
	LBG(V,Y_all, Y, universe, k, tokura_weights, epsilon,nframes,threshold);//run LBG algorithm on the dataset
	
}
void make_zero(vector < vector < vector <lld> > > & input)
{
	for (ll i = 0; i < input.size(); i++)
	{
		for (ll j = 0; j < input[i].size(); j++)
		{
			for (ll k = 0; k < input[i][j].size();k++)
			{
				input[i][j][k] = 0;
			}
		}
	}
}
void add_values_matrix(vector<vector<lld> > &old_input, vector<vector<lld> >&new_input)
{
	for (ll i = 0; i < new_input.size(); i++)
	{
		for (ll j = 0; j < new_input[i].size(); j++)
		{
			old_input[i][j] += new_input[i][j];
		}
	}
}
void divide_by(vector<vector <vector <lld> > > & input, lld divisor)
{
	for (ll i = 0; i < input.size(); i++)
	{
		for (ll j = 0; j < input[i].size(); j++)
		{
			for (ll k = 0; k < input[i][j].size(); k++)
			{
				input[i][j][k] =input[i][j][k]/divisor;
			}
		}
	}
}
void create_take_average(vector <vector<vector<lld> > >&new_A, vector <vector<vector<lld> > >&all_A)
{
	
	make_zero(all_A);
	for (ll i = 0; i < new_A.size(); i++)
	{
		add_values_matrix(all_A[i%all_A.size()], new_A[i]);
	}
	
	divide_by(all_A,new_A.size()/all_A.size());
}
void get_universe_and_observation_sequence(vector<vector<lld> > & universe, vector<vector<lld> > & V, vector<vector<lld> > & Y_all, vector<lld>&tokura_weights, vector<ll>&Y, lld &p, lld &num_samples, lld&epsilon, ll &iterations, ll &nframes, ll &k, lld &threshold,ll old_ones,ll offset)
{
	if (old_ones==0)
	{
		outer_wrapper(universe, V, Y_all, tokura_weights, Y, p, num_samples, epsilon, iterations, nframes, k, threshold, offset);
		store_values(Y_all, "Y_all.txt");
	}
	else
	{
		Y_all = get_matrix_file_input("Y_all.txt", nframes);
	}
}
void create_initial_AB(vector <vector<vector<lld> > > &all_A, vector <vector<vector<lld> > >& all_B,vector<lld>& digits,ll &N_states,ll &N_codebook, string& file_A_matrix,string& file_B_matrix)
{
	vector<vector<lld> >  A, B;
	A = get_A(file_A_matrix, N_states);//getting transition matrix from input file initial model
	B = get_B(file_B_matrix, N_codebook);//getting observation probability from input file initi
	for (int i = 0; i < digits.size(); i++)
	{
		all_A.push_back(A);
		all_B.push_back(B);
	}
}
void training(vector <vector<vector<lld> > > &all_A, vector <vector<vector<lld> > > &all_B, vector<vector<lld> >& Y_all, vector<lld> & pi_matrix, vector<lld> &digits, ll T, string &file_pi_matrix)
{
	vector <vector<vector<lld> > >  new_A, new_B;
	vector<vector<lld> > A, B;
	for (int j = 0; j < 3; j++)
	{
		avg_num = j;
		cout << "average number " << j << " --------------------------------------------------------------------------------------------------------------------" << endl;
		//for (int i = 0; i < Y_all.size(); i++)
		for (int i = 0; i <Y_all.size(); i++)
		{
			cout << i << " is the iteration number" << endl;
			iteration_num = i;
			A = all_A[i%digits.size()];
			B = all_B[i%digits.size()];
			pi_matrix = get_pi_matrix(file_pi_matrix);//getting state matrix from input file initial mod
			inner_wrapper(A, B, pi_matrix, addone(Y_all[i]), T);//runs forward, backward and viterbi algorithm and output the result
			new_A.push_back(A);
			new_B.push_back(B);
			//check_stochastic(A);
			//check_stochastic(B);

		}

		create_take_average(new_A, all_A);
		create_take_average(new_B, all_B);
		cout << new_A.size() << " " << new_B.size() << endl;
		cout << all_A.size() << " " << all_B.size() << endl;
		new_A.clear(), new_B.clear();
		
	}
	store_values(all_A,"all_A.txt");
	store_values(all_B, "all_B.txt");
}
lld find_observation_probability(vector<vector<lld> > & A,vector<vector<lld> > & B, vector<lld> &pi_matrix, vector<lld> & O, vector<lld>& digits, ll & T)
{
	vector <vector < vector <lld> > > zhi;
	vector <vector <lld> > alpha, beta, delta, psi, gamma;
	vector<lld> state_sequence;
	lld P_O_lambda, p_star, q_star, P_beta;



	forward_procedure(A, B, pi_matrix, O, T, alpha, P_O_lambda);//running forward procedure

	backward_procedure(A, B, O, T, beta, P_beta);//running backward procedure
	//cout << " -------------------------------------------------------------------------------------------------------------------------";
	//print_matrix(A); print_matrix(B), print_matrix(pi_matrix), print_matrix(O), print_matrix(alpha), print_matrix(beta);
	//cout << T << " " << P_O_lambda << " " << P_beta;
	//cout << " -------------------------------------------------------------------------------------------------------------------------";
	viterbi(A, B, pi_matrix, O, T, delta, psi, state_sequence, p_star, q_star);//running viterbi algorithm		
	return p_star;
}
ll classification(vector <vector<vector<lld> > > &all_A, vector <vector<vector<lld> > > &all_B,vector<lld> & pi_matrix,vector<lld>& O, vector<lld> & digits,ll& T)
{
	vector<lld> probabilities;
	lld max_value;
	ll max_index;
	for (ll i = 0; i < digits.size(); i++)
	{
		probabilities.push_back(find_observation_probability(all_A[i],all_B[i],pi_matrix,addone(O),digits,T));
	}
	max_value = probabilities[0], max_index = 0;
	for (ll i = 0; i < probabilities.size(); i++)
	{
		if (max_value < probabilities[i])
		{
			max_value = probabilities[i];
			max_index = i;
		}
	}
	print_matrix(probabilities);
	return digits[max_index];

}
void testing(vector <vector<vector<lld> > > &all_A, vector <vector<vector<lld> > > &all_B, vector<vector<lld> >& input_observations, vector<lld> & pi_matrix, vector<lld> &digits, ll T, string &file_pi_matrix)
{
	print_matrix(all_A);
	print_matrix(all_B);
	ll digit;
	for (ll i = 0; i < input_observations.size(); i++)
	{

		pi_matrix = get_pi_matrix(file_pi_matrix);//getting state matrix from input file initial mod
		digit=classification(all_A,all_B,pi_matrix,input_observations[i],digits,T);
		cout << filename_vector[i] << endl;
		cout << "digit classified is " << digit << endl;
	}
}
void outer_wrapper2(vector<vector<lld> > & X, vector<vector<lld> > & V, vector<vector<lld> > & Y_all, vector<lld>&tokura_weights, vector<ll>&Y, lld &p, lld &num_samples, lld&epsilon, ll &iterations, ll &nframes, ll &k, lld &threshold, ll offset)
{

	
	
	vector<vector<lld> > universe;
	Y_all.clear();
	for (int i = 0; i < X.size(); i++)Y.push_back(-1);
	universe = create_universe(iterations, nframes, p, num_samples, offset);
	X = universe;
	lld min_dist = numeric_limits<lld>::max(), dist;
	ll  num_clusters = 0, assigned_class = -1;
	vector<lld > temp;
	cout << "isze of X is "<< X.size() << endl;
	cout << "isze of V is"<<V.size() << endl;
	for (ll i = 0; i < X.size(); i++)
	{
		//cout << "is is" << i << endl;
		for (ll j = 0; j < V.size(); j++)
		{
			//cout << j << " ";
			dist = tokhura(X[i], V[j], tokura_weights);//finding tokhura's distance
			//cout << " aaaaaaaaa";
			if (min_dist>dist)
			{
				min_dist = dist;//finding the least tokhura's distance
				assigned_class = j;
			}
		}
		min_dist = numeric_limits<lld>::max();
		//cout << " aaaaaaaaa";
		//Y[i] = assigned_class;//using the distance cluster it to one of the class which has min distance
		//cout << " aaaaaaaaa";
		temp.push_back(assigned_class);
		if (temp.size() == nframes)			Y_all.push_back(temp), temp.clear();

	}
	cout << "Here";
	//LBG(V, Y_all, Y, universe, k, tokura_weights, epsilon, nframes, threshold);//run LBG algorithm on the dataset

}
void outer_wrapper3(vector<vector<lld> > & X, vector<vector<lld> > & V, vector<vector<lld> > & Y_all, vector<lld>&tokura_weights, vector<ll>&Y, lld &p, lld &num_samples, lld&epsilon, ll &iterations, ll &nframes, ll &k, lld &threshold, ll offset)
{

	lld min_dist = numeric_limits<lld>::max(), dist;
	ll  num_clusters = 0, assigned_class = -1;
	Y_all.clear();
	vector<lld > temp;
	cout << "isze of X is " << X.size() << endl;
	cout << "isze of V is" << V.size() << endl;
	for (ll i = 0; i < X.size(); i++)
	{
		//cout << "is is" << i << endl;
		for (ll j = 0; j < V.size(); j++)
		{
			//cout << j << " ";
			dist = tokhura(X[i], V[j], tokura_weights);//finding tokhura's distance
			//cout << " aaaaaaaaa";
			if (min_dist>dist)
			{
				min_dist = dist;//finding the least tokhura's distance
				assigned_class = j;
			}
		}
		min_dist = numeric_limits<lld>::max();
		//cout << " aaaaaaaaa";
		//Y[i] = assigned_class;//using the distance cluster it to one of the class which has min distance
		//cout << " aaaaaaaaa";
		temp.push_back(assigned_class);
		if (temp.size() == nframes)			Y_all.push_back(temp), temp.clear();

	}
	cout << "Here";
	//LBG(V, Y_all, Y, universe, k, tokura

	

}
vector<vector <lld> > recording_module_process(ll num_samples,ll nframes,ll p,string filename,int & flag)
{
	
		vector<lld>data, all_data;
		
		data = get_vector_file_input(filename,flag);
		vector<vector<lld> > dummy;
		if (flag == 0)return dummy;
		data = dc_shift_normalize_vac(data);//does dc shift, normalize the data and do voice activity detecting
		all_data = data;
		lld G2, item;
		vector <lld>  ai, ci, Ri, temp, digits = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		vector<vector <lld> >all_ais, all_cis, all_Ris, input_cis, universe;
		for (ll k = 0; k < nframes; k++)
		{
			//cout <<k<<" "<< all_data.size() << " " << k * 25 * num_samples / 100 + num_samples << endl;
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
		return universe;
	
}
void record_sound()
{
	cout << "NOTE:- It is expected that input at the beginning and ending has ambient sound of about 1 second" << endl << endl;
	cout << "The output of the recorded sound will be stored in out.wav and out.txt for sound and its numerical samples respectively" << endl << endl;
	system("Recording_Module.exe 4 out.wav out.txt");
}
int _tmain(ll  argc, _TCHAR* argv[])
{
	
	vector <vector<vector<lld> > > all_A, all_B;
	vector<vector<lld> > training_universe, training_V,testing_universe,testing_V;
	vector<vector<lld > > Y_all,input_observations;
	vector<lld>pi_matrix, tokura_weights = { 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61 }, digits = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };;
	vector<ll >training_Y,testing_Y;
	string file_A_matrix = "data/A_matrix_initialmodel.txt", file_B_matrix = "data/B_matrix_initialmodel.txt", file_pi_matrix = "data/pi_matrix_initialmodel.txt";
	lld vector_length = 12, num_samples_perframe = 320, epsilon = 0.03,threshold=0.1;//vector length is p
	ll  iterations_perdigit = 20, nframes = 40, k = 32, T=nframes, N_states = 5, N_codebook = 32,old_universe=0,train=1,test=1,option;
	
	get_universe_and_observation_sequence(training_universe, training_V, Y_all, tokura_weights, training_Y, vector_length, num_samples_perframe, epsilon, iterations_perdigit, nframes, k, threshold, old_universe, 0);
	cout << " training universe is " << training_universe.size() << endl;
	if (train)
	{
		create_initial_AB(all_A, all_B, digits, N_states, N_codebook, file_A_matrix, file_B_matrix);
		training(all_A, all_B, Y_all, pi_matrix, digits, T, file_pi_matrix);
	}
	while (1)
	{
		cout << "Press 0 to stop the algorithm | ";
		cout << "Press 1 to test against testing data already in the file | ";
		cout << "Press 2 to give recording module | ";
		cout << "Press 3 to give filename";
		cin >> option;
		if (option == 0)break;
		if (test&&option==1)
		{
			ll iterations_perdigit = 10;
			cout << training_V.size() << endl;
			outer_wrapper2(testing_universe, training_V, input_observations, tokura_weights, testing_Y, vector_length, num_samples_perframe, epsilon, iterations_perdigit, nframes, k, threshold, 20);
			print_matrix(input_observations);
			testing(all_A, all_B, input_observations, pi_matrix, digits, T, file_pi_matrix);
		}
		else if (option == 2||option==3)
		{
			string filename;
			if (option == 2)
			{
				record_sound();
				 filename = "out.txt";
			}
			else if (option == 3)
			{

				cin >> filename;
			}
			else
			{
				cout << "WRONG OPTION CHOOSE AGAIN" << endl;
				continue;
			}
			filename_vector.clear();
			filename_vector.push_back(filename);
			
			vector<vector <lld> > universe;
			int flag;
			universe=recording_module_process(num_samples_perframe,nframes,vector_length,filename,flag);
			if (flag == 0)
			{
				cout << "WRONG FILENAME TRY AGAIN" << endl;
				continue;
			}
			outer_wrapper3(universe, training_V, input_observations, tokura_weights, testing_Y, vector_length, num_samples_perframe, epsilon, iterations_perdigit, nframes, k, threshold, 20);
		
			print_matrix(input_observations);
			testing(all_A, all_B, input_observations, pi_matrix, digits, T, file_pi_matrix);
		}
		else
		{
			cout << "WRONG OPTION CHOOSE AGAIN" << endl;
		}
		
	}
	
	for (int i = 0; i < all_A.size(); i++)
	{
		check_stochastic(all_A[i]);
	}
	for (int i = 0; i < all_B.size(); i++)
	{
		check_stochastic(all_B[i]);
	}
	//print_matrix(all_A);
	//print_matrix(all_B);
	//inner_wrapper(A, B, pi_matrix, Y_all[0], T);//runs forward, backward and viterbi algorithm and output the result
	

	return 0;
}