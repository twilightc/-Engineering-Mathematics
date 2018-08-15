#pragma once
#include<string>
#include<vector>
#include<stack>
#include "DotNetUtilities.h"
#include<iostream>
#include<fstream>	
#include"DataManager.h"

class Mmatrix {
private:
public:
	Mmatrix();

	bool to_find(System::String^);
	bool search_func(std::string);
	void FindFunc(std::vector<Matrix>&, std::string,Ans*);

	bool Mmatrix::Deal_outstack(std::vector<std::string> &out, std::vector<Matrix> &matrix, std::vector<int> &rec);
	bool isDigit(std::string);
	bool search_menu(std::string);
	int search_Ms(std::string, std::vector<Matrix>&);

	std::vector<std::vector<double>> muti(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
	std::vector<std::vector<double>> Trsp(Matrix&);
	std::vector<std::vector<double>> Trsp(std::vector<std::vector<double> > &);
	std::vector<std::vector<double>> Utri(std::vector<std::vector<double> > &v);
	std::vector<std::vector<double>> Ltri(std::vector<std::vector<double> > &v);
	std::vector<std::vector<double>> inverse(std::vector<std::vector<double>>&, bool&);
	void Ach(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, int);
	void ch(std::vector<std::vector<double>>&);
	void chL(std::vector<std::vector<double>>&);
	std::vector<double> solve_linear_s(std::vector<std::vector<double>>&, std::vector<std::vector<double>> &);
	int rank(std::vector<std::vector<double>>&);
	
	int n_0_len(std::vector<double> &);
	int R_n0_len(std::vector<double> &);
	bool check_size_m(std::deque<Matrix> &);
	bool check_size_m_dot(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);

	bool C_size(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, int);

	Matrix vectorToM(std::vector<std::vector<double>>);		//連續使用函數時可用
	std::vector<std::vector<double>>  MTovector(Matrix *);
	std::vector<std::vector<double>>  MTovector(Matrix);		//連續使用函數時可用
	double  scalar(std::vector<std::vector<double> > &);		//matrix to scalar
	Matrix lu_decomposition(std::vector<std::vector<double>> &);
	std::vector<double> solvelineareq(Matrix &, std::vector<std::vector<double> > &);
	double determinant(std::vector <std::vector<double>>&);
	std::vector<std::vector<double> > InverseM(Matrix &);
	std::vector<std::vector<double> > adjoint(Matrix &);
	std::vector<std::vector<double> > least_suqare(Matrix &, Matrix &);

	void power_method(Matrix &,std::vector<double>*);
	void find_MatrixName(std::vector<Matrix> &, std::deque<Matrix> &, std::string &);	//finding certain matrix in vector<Matrix>
	std::string mathematics_operation(std::vector<int> &, std::vector<Matrix> &, std::vector<std::string> &);	//函數標題就是四則運算的英文
	std::string Mk_str(Ans&);
};