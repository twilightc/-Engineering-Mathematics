#pragma once
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
//定義向量資料結構
struct Vector
{
	int len;
	std::string index;
	std::string Name;
	std::vector<double> Data;

	//the following are operator overloading
	bool Vector::operator==(Vector &);
	Vector Vector::operator=(Vector &);
	Vector Vector::operator+(Vector &);
	Vector Vector::operator*(Vector);
	friend Vector Vector::operator*(double, Vector &);
	friend Vector Vector::operator*(std::vector<double>, Vector &);
};
struct Eigen {
	std::vector<std::vector<double> > e_vector;
	std::vector<double> e_value;
};
struct Matrix {
	int row_len, col_len;
	std::string index;
	std::vector<std::vector<double>> Datas;
	std::string Name;
	std::vector<std::vector<double>> U;			//uppertriangle
	std::vector<std::vector<double>> L;			//lowertriangle
	std::vector<std::vector<double>> P;
	Eigen eigen;								//eigen set
	Matrix Matrix::operator=(Matrix &);
	Matrix Matrix::operator=(Matrix *);
	Matrix* Matrix::operator+(Matrix &);
	Matrix* Matrix::operator+(Matrix *);
	friend Matrix* Matrix::operator*(double, Matrix &);
	friend Matrix Matrix::operator*(double, std::vector<std::vector<double>> &);
	friend Matrix* Matrix::operator*(std::vector<std::vector<double>> &, Matrix &);
	//std::vector<std::vector<double>>  Matrix::operator=(std::vector<std::vector<double>> &);
	friend Matrix* Matrix::operator*(Matrix &, Matrix &);
	friend std::vector<std::vector<double>> Matrix::operator*(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
	//friend std::vector<std::vector<double>> operator*(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
};
struct Ans {
	std::vector<std::vector<double>> vvAns;
	std::vector<double> vAns;
	std::vector<double> AnsData;
	std::string strAns;
	Matrix Mans;
	Vector Vans;
	bool bAns;
	int jud;
	double dAns;
};


//定義控管資料class
class DataManager
{
private:
	//儲存向量資料
	std::vector<Vector> Vectors;
	//紀錄向量ID，用於控管
	int  VectorVariableIndex;
	//紀錄檔案路徑名稱
	std::string FileName;
	//儲存矩陣資料
	std::vector<Matrix> matrix;
	//紀錄向量ID，用於控管
	int  MatrixVariableIndex;
	//紀錄檔案路徑名稱
	std::string MfileName;
public:
	DataManager();
	//讀取向量資料
	bool LoadVectorData();
	//取得向量資料
	std::vector<Vector> GetVectors();
	//設置檔案路徑名稱
	void SetFileName(std::string fileName);
	//**************************************************
	bool LoadData();
	std::vector<Matrix> getMatrix();
	void setName(std::string);
};
