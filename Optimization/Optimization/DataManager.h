#pragma once
#include<vector>
#include<iostream>
#include<fstream>
#include<string>

struct Matrix {
	int row_len, col_len;
	std::string index;
	std::vector<std::vector<double>> Datas;
	std::string Name;
	std::vector<std::vector<double>> U;			//uppertriangle
	std::vector<std::vector<double>> L;			//lowertriangle
	std::vector<std::vector<double>> P;
	Matrix Matrix::operator=(Matrix &);
	Matrix Matrix::operator=(Matrix *);
	Matrix* Matrix::operator+(Matrix &);
	Matrix* Matrix::operator+(Matrix *);

	friend Matrix* Matrix::operator*(double, Matrix &);								//四則運算用

	friend std::vector<std::vector<double>> Matrix::operator+(std::vector<std::vector<double>>, std::vector<std::vector<double>>); //new
	friend std::vector<std::vector<double>> Matrix::operator*(std::vector<std::vector<double>>, double);	//new
	friend std::vector<std::vector<double>> Matrix::operator*(std::vector<std::vector<double>>, std::vector<std::vector<double>>);

	friend Matrix Matrix::operator*(double, std::vector<std::vector<double>> &);

	friend Matrix* Matrix::operator*(std::vector<std::vector<double>> &, Matrix &);	//四則運算用
																					//std::vector<std::vector<double>>  Matrix::operator=(std::vector<std::vector<double>> &);
	friend Matrix* Matrix::operator*(Matrix &, Matrix &);
};

//*****************************************************************

//定義控管資料class
class DataManager
{
private:
	//儲存方程式資料
	std::vector<std::string> Equations;
	//紀錄向量ID，用於控管
	int EquationIndex;
	//紀錄檔案路徑名稱
	std::string FileName;
public:
	DataManager();
	//讀取向量資料
	bool LoadEquationData();
	//取得向量資料
	std::vector<std::string> GetEquations();
	//設置檔案路徑名稱
	void SetFileName(std::string fileName);
};