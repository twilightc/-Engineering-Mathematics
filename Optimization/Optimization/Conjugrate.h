#pragma once
#pragma once
#include<iostream>
#include<stdlib.h>
#include<vector>
#include<queue>
#include<algorithm>
#include"MathOp.h"
#include<algorithm>
#include"DataManager.h"
#include"DotNetUtilities.h"
#define eps 0.0001
#define eps2 0.0001
#define eps3 0.0001

class Conjugrate {
public:
	Conjugrate();
	std::vector<std::vector<double>> getPs();
	std::vector<std::vector<double>> trans(std::vector<std::vector<double>>);
	std::vector<std::vector<double>> hessian_2d(double, double, double, double, std::vector<std::string>);
	std::string nextline();
	std::string Outcome();

	void conjugrateAgorism(double, double, double, std::vector<std::string>);
	void conjugrateAgorism(double, double, double,double, std::vector<std::string>);
	void storePs(std::vector<std::vector<double>>, std::vector<std::string>);
	void recycle(std::vector<std::vector<double> > &);
	void recycle(std::vector<double > &);
	double find_a(double, double, double, double, std::vector<std::string>);

private:
	MathOp *mathop;
	std::string out;
	std::vector<std::vector<double>> steps;
};