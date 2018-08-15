#pragma once
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<stack>
#include"DataManager.h"
class Mvector {
private:
	Mvector *mv;	//convenient to operate function
public:	
	Mvector();

	bool to_find(System::String^);
	bool search_func(std::string);
	void FindFunc(std::vector<Vector>&, std::string,Ans*);

	int n_0_len(std::vector<double>&);
	bool isDigit(std::string);
	bool check_size(std::vector<Vector>&);

	double D_Product(Vector&, Vector&);
	double D_Product(std::vector<double>&, std::vector<double>&);
	std::vector<double> C_product(Vector&, Vector&);

	double Norm(std::vector<double>&);
	double Compo_a_on_b(std::vector<double>&, std::vector<double>&);
	double Tri_area(Vector&, Vector&);
	double Angle_V(Vector&, Vector&);
	std::vector<double> plus(std::vector<double>&, std::vector<double>&,int);
	std::vector<double> Proj(Vector&,Vector&);
	std::vector<double> Proj(std::vector<double>&, std::vector<double>&);
	std::vector<double> Normal(std::vector<double>&);
	bool Para_Jud(Vector&,Vector&);
	std::vector<double> Pl_normal(Vector&, Vector&);	
	std::vector<std::vector<double>> Orth_basis(std::vector<Vector>&);
	bool Linear_Jud(Vector&,Vector&);
	bool Linear_Jud(std::vector<Vector>&);
	bool Orth_Jud(Vector&, Vector&);

	void ch(std::vector<Vector>&);
	bool check_size_v(std::deque<Vector> &);
	bool find_operator(std::string);		//finding cretain operator
	void find_VectorName(std::vector<Vector> &, std::deque<Vector> &, std::string &);	//finding certain vector in vector<Vector>
	std::string mathematics_operation(std::vector<Vector> &, std::vector<std::string> &,std::vector<int>&);	//函數標題就是四則運算的英文
	std::string Mk_str(Ans&);
};