#pragma once
#include<vector>
#include<iostream>
#include<string>

class MathOp {
	std::vector<std::vector<std::string>> fs;
	std::vector<std::vector<double>> steps;
public:
	MathOp();
	std::vector<std::string> split_to_format(std::string);
	std::vector<std::string> format(std::vector<std::string>);
	std::vector<std::vector<double>> getPs();

	bool ch_get_fs(int, std::vector<std::string>&);

	void storeFx(std::vector<std::string>);
	void storePs(std::vector < std::vector < double >>, std::vector<std::string>);
	void recycle(std::vector<std::vector<double>>&);
	void recycle(std::vector<double>&);
	void print();

	double scalar(std::vector<std::vector<double>>&);
	double ComputeFx(double, double, std::vector<std::string>);
	double ComputeFx(double, std::vector<std::string>);
	double f_prime(double, double, std::vector<std::string>);
	double s_prime(double, double, std::vector<std::string>);
	double fx(double, double, double, std::vector<std::string>);
	double fy(double x, double y, double k, std::vector<std::string>);
	double fxx(double, double, double, std::vector<std::string>);
	double fxy(double, double, double, double, std::vector<std::string>);
	double fyy(double, double, double, std::vector<std::string>);
	bool isSymbol(std::string);
	bool isDigit(std::string);

private:

	std::vector<std::string> pow_backet(std::vector<std::string>);
	std::vector<std::string> mut_backet(std::vector<std::string>);
	std::vector<std::string> sub_backet(std::vector<std::string>);
	std::vector<std::string> div_backet(std::vector<std::string>);
};
