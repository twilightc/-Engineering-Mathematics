#pragma once
#pragma once
#include<iostream>
#include<stdlib.h>
#include<vector>
#include<queue>
#include<algorithm>
#include"MathOp.h"
#include"DataManager.h"
#include"DotNetUtilities.h"

class Qnewton {
	std::vector<std::string> fx;
	std::vector<std::vector<double>> xs;
	std::string out;
	MathOp *mathop;
public:
	Qnewton();
	Qnewton(std::vector<std::string>);
	void QnewtonAlgorism(double, double, std::vector<std::string>);
	void QnewtonAlgorism(double, std::vector<std::string>);
	void storeXs(std::vector<std::vector<double>>, std::vector<std::string>);
	void clean();

	double getAlpha(std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::string>);
	double getAlpha(std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>);

	std::string nextline();
	std::string getOut();
	std::vector<std::vector<double>> trans(std::vector<std::vector<double>>);
	std::vector<std::vector<double>> getPs();
};
