#pragma once
#include<iostream>
#include<stdlib.h>
#include<vector>
#include<queue>
#include<algorithm>
#include"MathOp.h"
#include"DataManager.h"
#define eps 0.0001

typedef double(*func)(double, double);

typedef double(*func2)(double, double);

class Newton {
	MathOp *mathop;
	std::vector<std::string> fx;
	std::string out;
	std::vector<std::vector<double>>xs;
public:
	Newton();
	Newton(std::vector<std::string>);
	void recycle(std::vector<std::vector<double>>&);
	void recycle(std::vector<double>&);
	void clearOut();
	void NewtonAlgorism(double, double, double, double, std::vector<std::string>);
	void NewtonAlgorism(double, double, double, std::vector<std::string>);
	void storePs(std::vector<std::vector<double>>, std::vector<std::string>);

	std::string Outcome();
	std::string nextline();

	std::vector<double> solvelineareq(Matrix&, std::vector<std::vector<double>>&);
	std::vector<std::vector<double>> inv_hessian(std::vector<std::vector<double>>&);
	std::vector<std::vector<double>> InverseM(Matrix&);
	std::vector<std::vector<double>> hessian_2d(double, double, double, double, std::vector<std::string>);
	std::vector<std::vector<double>> hessian_1d(double, double, double, std::vector<std::string>);
	std::vector<std::vector<double>> getXs();

	Matrix lu_decomposition(std::vector<std::vector<double> > &);
};