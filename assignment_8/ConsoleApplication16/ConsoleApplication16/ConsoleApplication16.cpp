// ConsoleApplication16.cpp : Defines the entry point for the console application.
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
#define ll long long int
#define lld long double
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
void print_vector(vector<lld>input)//printing a single dimension vector with datatype double
{
	for (ll i = 0; i < input.size(); i++)		cout << input[i] << " ";	nl
}
void print_vector(vector<ll>input)//printing a single dimension vector with datatype int
{
	for (ll i = 0; i < input.size(); i++)	cout << input[i] << " ";	nl
}
void print_matrix(vector <vector <lld> > input)//printing a matrix with datatype double
{
	for (ll i = 0; i < input.size(); i++)		print_vector(input[i]); nl
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
	while (!infile.eof())
	{
		infile >> item;
		output.push_back(item);
	}
	return output;
}
vector < vector <lld> >  get_A(string filename,ll N)//getting transition matrix from input file initial model
{
	return get_matrix_file_input(filename,N);
}
vector< vector <lld> > get_B(string filename, ll M)//getting observation probability from input file initial model
{
	return get_matrix_file_input(filename,M);
}
vector <lld> get_pi_matrix(string filename)//getting state matrix from input file initial model
{
	return get_vector_file_input(filename);
}
vector<lld> get_observation_sequence(string filename)//getting observation sequence vector from input file initial model
{
	return get_vector_file_input(filename);
}
void forward_initialization(vector<vector<lld> >&alpha,vector < lld> pi_matrix, vector <vector <lld> > B,vector<lld> O)
{
	vector < lld> temp;	
	for (ll i = 0; i < pi_matrix.size(); i++)		temp.push_back(pi_matrix[i] * B[i][O[0]-1]);
	alpha.push_back(temp);
}
void forward_induction(vector<vector <lld> > &alpha, vector < vector<lld> >& A, vector<vector<lld> >& B, vector <lld>& O,ll &T)
{
	vector<lld> alpha_temp;
	lld temp;
	for (ll t = 0; t < T - 1; t++)
	{
		alpha_temp.clear();
		for (ll j = 0; j < A.size(); j++)
		{
			temp = 0;
			for (ll i = 0; i < A.size();i++)
			{
				temp+=alpha[alpha.size() - 1][i] * A[i][j];
			}
			alpha_temp.push_back(temp*B[j][O[t+1]-1]);
		}
		alpha.push_back(alpha_temp);
	}
}
void termination(lld &prob,vector<vector<lld> > & alpha)
{
	prob = 0;
	for (ll i = 0; i < alpha[0].size(); i++)prob += alpha[alpha.size()-1][i];
}
void forward_procedure(vector < vector <lld> > A, vector < vector < lld> > B, vector <lld> pi_matrix, vector <lld> O,ll T)
{
	vector <vector <lld> > alpha;
	lld P_O_lambda=0;
	//initialization step
	forward_initialization(alpha,pi_matrix,B,O);
	
	//induction step
	forward_induction(alpha,A,B,O,T);
	
	//terminationa step
	termination(P_O_lambda,alpha);

	cout << "output of forward procedure is = ";
	cout << P_O_lambda << endl;
}
void backward_initialization(vector<vector<lld> > &beta,ll N)
{
	beta.push_back(vector<lld> (N,1));
}
void backward_procedure(vector < vector <lld> > A, vector < vector < lld> > B,  vector <lld> O, ll T)
{
	vector<vector<lld> > beta;
	//initialization step
	backward_initialization(beta,A.size());
	
}
void wrapper(vector < vector <lld> > A, vector < vector < lld> > B, vector <lld> pi_matrix, vector <lld> O, ll T)
{
	forward_procedure(A, B, pi_matrix, O, T);
	backward_procedure(A, B, O, T);
}
int _tmain(int argc, _TCHAR* argv[])
{
	string file_A_matrix = "data/A_matrix_initialmodel.txt", file_B_matrix = "data/B_matrix_initialmodel.txt", file_pi_matrix = "data/pi_matrix_initialmodel.txt";
	string file_observation1 = "data/observation_sequence_1_initialmodel.txt", file_observation2 = "data/observation_sequence_2_initialmodel.txt";
	vector <vector <lld> > A, B;//A is transition matrix , B is observation probability
	vector <lld> pi_matrix,Observation_sequence1,Observation_sequence2;//pi_matrix is the state matrix
	ll N_states = 5,N_codebook=32,T=85;//N_states is number of states , N_codebook is the size of codebook
	A=get_A(file_A_matrix,N_states);//getting transition matrix from input file initial model
	B = get_B(file_B_matrix,N_codebook);//getting observation probability from input file initial model
	pi_matrix=get_pi_matrix(file_pi_matrix);//getting state matrix from input file initial model
	Observation_sequence1 = get_observation_sequence(file_observation1);
	Observation_sequence2 = get_observation_sequence(file_observation2);
	wrapper(A,B,pi_matrix,Observation_sequence1,T);
	return 0;
}
