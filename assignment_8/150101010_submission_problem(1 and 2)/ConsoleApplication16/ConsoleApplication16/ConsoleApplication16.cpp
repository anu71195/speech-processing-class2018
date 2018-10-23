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
void forward_initialization(vector<vector<lld> >&alpha,vector < lld> pi_matrix, vector <vector <lld> > B,vector<lld> O)//initialization of forward procedure
{
	vector < lld> temp;	
	for (ll i = 0; i < pi_matrix.size(); i++)		temp.push_back(pi_matrix[i] * B[i][O[0]-1]);
	alpha.push_back(temp);
}
void forward_induction(vector<vector <lld> > &alpha, vector < vector<lld> >& A, vector<vector<lld> >& B, vector <lld>& O,ll &T)//induction step in forward procedure
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
				temp+=alpha[alpha.size() - 1][i] * A[i][j];//summation term
			}
			alpha_temp.push_back(temp*B[j][O[t+1]-1]);//creating vector for every time t is updated
		}
		alpha.push_back(alpha_temp);//storing vector for time t in alpha 
	}
}
void termination(lld &prob,vector<vector<lld> > & alpha)//finding probability of observation given model
{
	prob = 0;
	for (ll i = 0; i < alpha[0].size(); i++)prob += alpha[alpha.size()-1][i];
}
void forward_procedure(vector < vector <lld> > A, vector < vector < lld> > B, vector <lld> pi_matrix, vector <lld> O,ll T)//forward procedure
{
	vector <vector <lld> > alpha;
	lld P_O_lambda=0;
	//initialization step
	forward_initialization(alpha,pi_matrix,B,O);
	
	//induction step
	forward_induction(alpha,A,B,O,T);
	
	//terminationa step
	termination(P_O_lambda,alpha);

	cout << "P(O/lambda) = " << P_O_lambda; nl nl
	cout << "alpha matrix is"; nl;
	print_matrix(alpha);
}
void backward_initialization(vector<vector<lld> > &beta,ll N)//initialization of backward procedure
{
	beta.push_back(vector<lld> (N,1));//at time t =T all the states have value 1
}
void backward_induction(vector<vector<lld> > & beta, vector<vector<lld> > A, vector < vector < lld> > B, vector <lld> O, ll T)//induction step in backward procedure
{
	vector<lld> b_temp;
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
			}
			b_temp.push_back(temp);//storing values from summation term in vector for every time t for every state
		}
		beta.insert(beta.begin(), b_temp);//making a matrix beta by storing values from b_temp vector for each time and each state 1 to T and 1 to N
	}
}
void backward_procedure(vector < vector <lld> > A, vector < vector < lld> > B,  vector <lld> O, ll T)//backward procedure
{
	vector<vector<lld> > beta;
	//initialization step
	backward_initialization(beta,A.size());
	
	//induction step
	backward_induction(beta,A,B,O,T);
	cout << "alpha matrix is"; nl;
	print_matrix(beta);
	

}
void viterbi_initialization(vector<vector<lld> > &delta, vector< vector < lld> > & psi, vector<vector<lld> > b,vector<lld> pi_matrix,vector<lld> O)//viterbi algorithm initialization
{
	vector<lld> temp;
	for (ll i = 0; i < pi_matrix.size(); i++)		temp.push_back(pi_matrix[i]*b[i][O[0]-1]);
	delta.push_back(temp), psi.push_back(vector<lld> (pi_matrix.size(),0));
}
void viterbi_induction(vector<vector<lld> > &delta, vector < vector < lld> >& psi, vector < vector <lld> > A, vector < vector < lld> > b, vector<lld>O, ll T)//induction step in viterbi algorithm
{
	vector<lld> delta_temp,psi_temp;
	lld temp=INT_MIN,temp2;
	ll index;
	for (ll t = 1; t < T; t++)
	{
		delta_temp.clear(),psi_temp.clear();
		for (ll j = 0; j < A.size(); j++)
		{
			temp = INT_MIN;
			
			for (ll i = 0; i < A.size(); i++)
			{
				temp2 = delta[delta.size() - 1][i] * A[i][j];
				if (temp < temp2)
				{
					temp = temp2,	index = i;			//max value step and its corresponding index
				}
			}
			delta_temp.push_back(temp*b[j][O[t] - 1]);	//storring max value for each time t and state in delta_temp
			psi_temp.push_back(index);//and its corresponding index psi_temp
		}
		delta.push_back(delta_temp);//storing delta_temp vector in delta matrix
		psi.push_back(psi_temp);//storing psi_temp vector in psi matrix
	}
}
void viterbi_termination(vector< vector< lld> > delta,lld &p_star,lld &q_star)//viterbi terimination step
{
	p_star = INT_MIN;
	for (ll i = 0; i < delta[delta.size() - 1].size(); i++)
	{
		if (p_star < delta[delta.size() - 1][i])//getting the biggest probability from delta and storing it in p_star and its corresponding index in q_star
		{
			p_star = delta[delta.size()-1][i];
			q_star = i;
		}
	}
}
vector<lld> viterbi_backtracking(vector<vector<lld> > psi, lld q_star)//back tracking best state sequence
{
	vector< lld > output;
	for (int t = psi.size() - 2; t >= 0; t--)
	{
		output.insert(output.begin(), q_star);//storing q_star or state sequence values in output
		q_star = psi[t+1][q_star];//from the q_star obtained from the last step getting q_star for previous step
	}
	output.insert(output.begin(), q_star);
	return output;
}
void viterbi(vector<vector<lld> > A,vector < vector <lld> > b, vector<lld> pi_matrix, vector<lld> O, ll T)
{
	vector<vector<lld> > delta, psi;
	vector<lld> state_sequence;
	lld p_star, q_star;
	//initialization
	viterbi_initialization(delta,psi,b,pi_matrix,O);
	
	//induction step
	viterbi_induction(delta,psi,A,b,O,T);
	
	//termination step
	viterbi_termination(delta,p_star,q_star);

	
	//backtracking step
	state_sequence=viterbi_backtracking(psi,q_star);
	cout << "delta is"; nl
	print_matrix(delta);
	cout << "psi is"; nl
	print_matrix(psi);
	cout << "state sequence is"; nl
	print_vector(state_sequence);
	cout << "p* = " << p_star << endl;
}
void wrapper(vector < vector <lld> > A, vector < vector < lld> > B, vector <lld> pi_matrix, vector <lld> O, ll T)
{
	forward_procedure(A, B, pi_matrix, O, T);//running forward procedure
	backward_procedure(A, B, O, T);//running backward procedure
	viterbi(A,B,pi_matrix,O,T);//running viterbi algorithm
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
	Observation_sequence1 = get_observation_sequence(file_observation1);//reading observation sequence 1 from file input
	Observation_sequence2 = get_observation_sequence(file_observation2);//reading observation sequence 2 from file input
	wrapper(A,B,pi_matrix,Observation_sequence1,T);//runs forward, backward and viterbi algorithm and output the result
	return 0;
}

