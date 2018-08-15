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

	friend Matrix* Matrix::operator*(double, Matrix &);								//�|�h�B���

	friend std::vector<std::vector<double>> Matrix::operator+(std::vector<std::vector<double>>, std::vector<std::vector<double>>); //new
	friend std::vector<std::vector<double>> Matrix::operator*(std::vector<std::vector<double>>, double);	//new
	friend std::vector<std::vector<double>> Matrix::operator*(std::vector<std::vector<double>>, std::vector<std::vector<double>>);

	friend Matrix Matrix::operator*(double, std::vector<std::vector<double>> &);

	friend Matrix* Matrix::operator*(std::vector<std::vector<double>> &, Matrix &);	//�|�h�B���
																					//std::vector<std::vector<double>>  Matrix::operator=(std::vector<std::vector<double>> &);
	friend Matrix* Matrix::operator*(Matrix &, Matrix &);
};

//*****************************************************************

//�w�q���޸��class
class DataManager
{
private:
	//�x�s��{�����
	std::vector<std::string> Equations;
	//�����V�qID�A�Ω󱱺�
	int EquationIndex;
	//�����ɮ׸��|�W��
	std::string FileName;
public:
	DataManager();
	//Ū���V�q���
	bool LoadEquationData();
	//���o�V�q���
	std::vector<std::string> GetEquations();
	//�]�m�ɮ׸��|�W��
	void SetFileName(std::string fileName);
};