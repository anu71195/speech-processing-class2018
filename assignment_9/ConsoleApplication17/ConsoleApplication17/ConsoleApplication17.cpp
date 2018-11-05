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
void print_matrix(vector<lld>input)//printing a single dimension vector with datatype int
{
	for (ll i = 0; i < input.size(); i++)	cout << input[i] << " ";	nl
}
void print_matrix(vector <vector <lld> > input)//printing a matrix with datatype double
{
	for (ll i = 0; i < input.size(); i++)		print_vector(input[i]); nl
}
void print_matrix(vector < vector <vector <lld> > > input)
{
	for (ll i = 0; i < input.size(); i++)cout << i << " " << endl , print_matrix(input[i]); nl
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
vector < vector <lld> >  get_A(string filename, ll N)//getting transition matrix from input file initial model
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
void check_stochastic(vector<lld> input)
{
	lld temp = 0;
	for (ll i = 0; i < input.size(); i++)
	{
		temp += input[i];

	}
	cout << temp << endl;
}

void check_stochastic(vector<vector<lld> > input)
{
	lld temp = 0;
	for (ll i = 0; i < input.size(); i++)
	{
		temp = 0;
		for (ll j = 0; j < input[i].size(); j++)
		{
			temp += input[i][j];
		}
		cout << temp << endl;
	}
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
void forward_procedure(vector < vector <lld> > A, vector < vector < lld> > B, vector <lld> pi_matrix, vector <lld> O, ll T, vector <vector <lld> > &alpha,lld &P_O_lambda)//forward procedure
{
	//print_matrix(B);
	//check_stochastic(B);
	//initialization step
	//print_matrix(pi_matrix);
	alpha.clear();
	forward_initialization(alpha, pi_matrix, B, O);
	//print_matrix(alpha);
	//induction step
	forward_induction(alpha, A, B, O, T);

	//terminationa step
	termination(P_O_lambda, alpha);

	//cout << "P(O/lambda) = " << P_O_lambda; nl nl
		//		cout << "alpha matrix is"; nl;
		//print_matrix(alpha);
	//	cout << alpha.size() << endl;
	//print_matrix(alpha);
	//cout << alpha.size() << endl;
	store_values(P_O_lambda, "POlambda.txt");
	store_values(alpha, "alpha.txt");
}
void backward_initialization(vector<vector<lld> > &beta, ll N)//initialization of backward procedure
{
	beta.push_back(vector<lld>(N, 1));//at time t =T all the states have value 1
}
void backward_induction(vector<vector<lld> > & beta, vector<vector<lld> >& A, vector < vector < lld> > &B, vector <lld> &O, ll& T)//induction step in backward procedure
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
void backward_termination(vector<vector<lld> > & beta,lld &P)
{
	P = 0;
	for (int i = 0; i<beta[0].size(); i++) //this is for self assessment and not required as there is no further use
	{
		P += beta[0][i];
		//cout << beta[0][i] << endl;
	}
}
void backward_procedure(vector < vector <lld> > A, vector < vector < lld> > B, vector <lld> O, ll T, vector<vector<lld> > &beta,lld& P_beta)//backward procedure
{
	
	//initialization step
	beta.clear();
	backward_initialization(beta, A.size());

	//induction step
	backward_induction(beta, A, B, O, T);
	//cout << "beta matrix is"; nl;
//	print_matrix(beta);
	backward_termination(beta,P_beta);
	store_values(beta, "beta.txt");


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
	for (int t = psi.size() - 2; t >= 0; t--)
	{
		output.insert(output.begin(), q_star);//storing q_star or state sequence values in output
		q_star = psi[t + 1][q_star];//from the q_star obtained from the last step getting q_star for previous step
	}
	output.insert(output.begin(), q_star);
	return output;
}
void viterbi(vector<vector<lld> > A, vector < vector <lld> > b, vector<lld> pi_matrix, vector<lld> O, ll T, vector<vector<lld> > &delta, vector<vector<lld> > &psi,vector<lld> &state_sequence,lld &  p_star, lld & q_star)
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
//	cout << "delta is"; nl
//		print_matrix(delta);
//	cout << "psi is"; nl
//		print_matrix(psi);
	//cout << "state sequence is"; nl
	//	print_vector(state_sequence);
//	cout << "p* = " << p_star << endl;nl
	store_values(delta, "delta.txt");
	store_values(psi, "psi.txt");
	store_values(state_sequence, "state_sequence.txt");
	store_values(p_star, "p_star.txt");
}
lld get_polambda_at_t(vector<vector<lld> > & A, vector < vector < lld> > & B, vector < vector <lld> > & alpha, vector < vector <lld> > & beta, vector<lld> & O, ll t)
{
	lld polambda = 0;
	ll N = A[0].size();
	for (ll i = 0; i < N;i++)
		for (ll j = 0; j < N;j++)
			polambda += alpha[t][i] * A[i][j] * B[j][O[t + 1] - 1] * beta[t + 1][j];
	return polambda;
}
vector <vector <lld> > get_zhi_at_t(vector<vector<lld> > & A, vector < vector < lld> > & B, vector < vector <lld> > & alpha, vector < vector <lld> > & beta, vector<lld> & O, ll t)
{
	vector<vector<lld> > zhi_t(A[0].size(),vector<lld>(A[0].size(),0));
	ll N = A[0].size();
	lld polambda = get_polambda_at_t(A, B, alpha, beta, O, t);
	for (ll i = 0; i < N; i++)
		for (ll j = 0; j < N; j++)
		{
			if(polambda)zhi_t[i][j] = alpha[t][i] * A[i][j] * B[j][O[t + 1] - 1] * beta[t + 1][j] / polambda;
			else zhi_t[i][j] = 0;
		}
	return zhi_t;

}
void findgamma(vector <vector < vector <lld> > > &zhi, vector < vector <lld> > & gamma)
{
	vector < vector < lld > > output(zhi.size(), vector<lld>(zhi[0].size(), 0));
	for (ll t = 0; t < zhi.size(); t++)
		for (ll i = 0; i < zhi[0].size(); i++)
			output[t][i] = vector_sum(zhi[t][i]);
	gamma = output;
}
void findzhi(vector <vector < vector <lld> > > &zhi, vector<vector<lld> > & A, vector < vector < lld> > & B, vector < vector <lld> > & alpha, vector < vector <lld> > & beta, vector<lld> & O, ll T)
{
	//print_matrix(beta);
	for (ll t = 0; t < T-1; t++)
		zhi.push_back(get_zhi_at_t(A,B,alpha,beta,O,t));	
}
void update_pi_matrix(vector<vector <lld> > & gamma,vector <lld> &pi_matrix)
{
	//print_matrix(gamma);
	for (ll i = 0; i < pi_matrix.size(); i++)		pi_matrix[i]=gamma[0][i];

}
lld get_sum_gamma_i(vector< vector<lld> > & gamma, int i)
{
	lld output = 0;
	for (ll t = 0; t < gamma.size(); t++)		output += gamma[t][i];
	return output;

}
lld get_sum_zhi_ij(vector<vector<vector<lld> > > & zhi, int i, int j)
{
	lld output = 0;
	for (ll t = 0; t < zhi.size(); t++)output += zhi[t][i][j];
	return output;

}
void update_A_matrix(vector <vector < vector <lld> > > &zhi, vector<vector<lld> > & gamma,vector < vector < lld> > & A)
{
	//print_matrix(gamma);
	for (ll i = 0; i < A.size(); i++)
	{
		lld sum_gamma_ti = get_sum_gamma_i(gamma,i);
		//cout << "sumgamma is"<<sum_gamma_ti << endl;
		for (ll j = 0; j < A[i].size(); j++)
		{
			if(sum_gamma_ti)A[i][j] = get_sum_zhi_ij(zhi, i, j) / sum_gamma_ti;
			else A[i][j] = 0;
			//cout << A[i][j] << endl;
		}
	}
}
lld get_bjk_num(vector<vector< lld> > & gamma, ll j, vector< lld> & O,ll K)
{
	lld output = 0;
	for (ll t = 0; t < gamma.size(); t++)if((O[t]-1)==K)output += gamma[t][j];
	return output;
}
lld get_bjk_den(vector<vector<lld> > & gamma, ll j)
{
	lld output = 0;
	for (ll t = 0; t < gamma.size(); t++)output += gamma[t][j];
	return output;
}
void update_B_matrix(vector< vector <lld> > & gamma, vector<lld> & O,vector < vector<lld> > & B)
{
	for (ll j = 0; j < B.size(); j++)
	{
		lld bjk_den = get_bjk_den(gamma, j);
		for (ll k = 0; k < B[j].size(); k++)
		{
			if(bjk_den)B[j][k] = get_bjk_num(gamma, j, O, k) / bjk_den;
			else B[j][k] = 0;
		}
	}	

}
void EM(vector <vector < vector <lld> > > &zhi, vector<vector<lld> > & A, vector < vector < lld> > & B, vector < vector <lld> > & alpha, vector < vector <lld> > & beta, vector<vector<lld> > &gamma, vector<lld> & O, ll T, vector <lld> &pi_matrix)
{
	zhi.clear(), gamma.clear();
	
	findzhi(zhi, A, B, alpha, beta, O, T);
	//cout << "zhi is" << endl;
	//print_matrix(zhi); 
	
	findgamma(zhi,gamma);
	//cout << "gamma is " << endl;
	//print_matrix(zhi);

	
	update_pi_matrix(gamma,pi_matrix);
//	cout << "update pi_matr is" << endl;
//	print_vector(pi_matrix);

	update_A_matrix(zhi, gamma, A);
	//cout << "update A matrix is " << endl;
	//print_matrix(A);

	update_B_matrix(gamma,O,B);
	//cout << "updated B Matrix is " << endl;
	//print_matrix(B);





}


void wrapper(vector < vector <lld> >& A, vector < vector < lld> >& B, vector <lld>& pi_matrix, vector <lld> O, ll T)
{
	
	
	for (int i = 0; i < 20; i++)
	{
		vector <vector < vector <lld> > > zhi;
		vector <vector <lld> > alpha, beta, delta, psi, gamma;
		vector<lld> state_sequence;
		lld P_O_lambda, p_star, q_star,P_beta;

		//cout << "Old" << endl;
		//print_matrix(gamma);
		forward_procedure(A, B, pi_matrix, O, T, alpha, P_O_lambda);//running forward procedure
		//print_matrix(pi_matrix);
		//print_matrix(alpha);
		backward_procedure(A, B, O, T, beta,P_beta);//running backward procedure

//		cout << beta.size() << endl;
	//	print_matrix(beta);
		viterbi(A, B, pi_matrix, O, T, delta, psi, state_sequence, p_star, q_star);//running viterbi algorithm		
		
		
		
		
		EM(zhi, A, B, alpha, beta, gamma, O, T, pi_matrix);
		//print_matrix(zhi);
		//print_matrix(state_sequence);
		cout << P_O_lambda << endl;
		
	}
}
int _tmain(int argc, _TCHAR* argv[])
{
	string file_A_matrix = "data/A_matrix_initialmodel.txt", file_B_matrix = "data/B_matrix_initialmodel.txt", file_pi_matrix = "data/pi_matrix_initialmodel.txt";
	string file_observation1 = "data/observation_sequence_1_initialmodel.txt", file_observation2 = "data/observation_sequence_2_initialmodel.txt";
	vector <vector <lld> > A, B;//A is transition matrix , B is observation probability
	vector <lld> pi_matrix, Observation_sequence1, Observation_sequence2;//pi_matrix is the state matrix
	ll N_states = 5, N_codebook = 32, T = 85;//N_states is number of states , N_codebook is the size of codebook
	A = get_A(file_A_matrix, N_states);//getting transition matrix from input file initial model
	B = get_B(file_B_matrix, N_codebook);//getting observation probability from input file initial model(emission matrix)
	pi_matrix = get_pi_matrix(file_pi_matrix);//getting state matrix from input file initial model
	Observation_sequence1 = get_observation_sequence(file_observation1);//reading observation sequence 1 from file input
	Observation_sequence2 = get_observation_sequence(file_observation2);//reading observation sequence 2 from file input
	wrapper(A, B, pi_matrix, Observation_sequence1, T);//runs forward, backward and viterbi algorithm and output the result
	return 0;
}
