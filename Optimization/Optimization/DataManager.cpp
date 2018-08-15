#include "DataManager.h"

DataManager::DataManager()
{
	EquationIndex = 0;
}

bool DataManager::LoadEquationData()
{
	std::fstream fin;
	//開啟檔案，傳入open函數的參數有兩個，欲開起的檔案名稱，開啟檔案的模式參數(這邊std::ios::in為讀取(輸入)狀態)
	fin.open(FileName, std::ios::in);
	//讀取失敗回傳false
	if (!fin)
	{
		return false;
	}
	else
	{
		//定義讀取檔案字串暫存變數
		std::string tempSring;

		//執行讀檔迴圈，並在讀到檔案結尾時結束
		while (!fin.eof())
		{
			//從檔案讀取字串
			fin >> tempSring;
			//解析到向量標記"V"
			Equations.push_back(tempSring);
			//遞增EquationIndex，標記到當前讀取向量ID
			EquationIndex++;
		}
		return true;
	}
}

std::vector<std::string> DataManager::GetEquations()
{
	return Equations;
}

void DataManager::SetFileName(std::string fileName)
{
	FileName = fileName;
}


Matrix Matrix::operator=(Matrix &mat) {
	std::vector<double> d;

	//info copy
	row_len = mat.row_len;
	col_len = mat.col_len;
	Name = mat.Name;
	index = mat.index;
	if (!mat.Datas.empty()) {
		for (int i = 0; i < mat.Datas.size(); i++) {
			Datas.push_back(d);
			for (int j = 0; j < mat.Datas[0].size(); j++) {
				Datas[i].push_back(mat.Datas[i][j]);
			}
		}
	}
	else if (!mat.L.empty()) {
		d.clear();
		for (int i = 0; i < mat.L.size(); i++) {
			L.push_back(d);
			for (int j = 0; j < mat.L[0].size(); j++) {
				L[i].push_back(mat.L[i][j]);
			}
		}
	}
	else if (!mat.U.empty()) {
		d.clear();
		for (int i = 0; i < mat.U.size(); i++) {
			U.push_back(d);
			for (int j = 0; j < mat.U[0].size(); j++) {
				U[i].push_back(mat.U[i][j]);
			}
		}
	}
	else if (!mat.P.empty()) {
		d.clear();
		for (int i = 0; i < mat.P.size(); i++) {
			P.push_back(d);
			for (int j = 0; j < mat.P[0].size(); j++) {
				P[i].push_back(mat.P[i][j]);
			}
		}
	}

	return *this;
}

Matrix Matrix::operator=(Matrix *mat) {
	std::vector<double> d;
	row_len = mat->row_len;
	col_len = mat->col_len;
	Name = mat->Name;
	index = mat->index;

	if (!mat->Datas.empty()) {
		for (int i = 0; i < mat->Datas.size(); i++) {
			Datas.push_back(d);
			for (int j = 0; j < mat->Datas[0].size(); j++) {
				Datas[i].push_back(mat->Datas[i][j]);
			}
		}
	}
	else if (!mat->L.empty()) {
		d.clear();
		for (int i = 0; i < mat->L.size(); i++) {
			L.push_back(d);
			for (int j = 0; j < mat->L[0].size(); j++) {
				L[i].push_back(mat->L[i][j]);
			}
		}
	}
	else if (!mat->U.empty()) {
		d.clear();
		for (int i = 0; i < mat->U.size(); i++) {
			U.push_back(d);
			for (int j = 0; j < mat->U[0].size(); j++) {
				U[i].push_back(mat->U[i][j]);
			}
		}
	}
	else if (!mat->P.empty()) {
		d.clear();
		for (int i = 0; i < mat->P.size(); i++) {
			P.push_back(d);
			for (int j = 0; j < mat->P[0].size(); j++) {
				P[i].push_back(mat->P[i][j]);
			}
		}
	}
	return *this;
}


Matrix* Matrix::operator+(Matrix &mat) {
	Matrix *ma = new Matrix();
	std::vector<double> d;		//輔助進行加法，因為兩個二維vector容器無法直接進行加減
								//初始化，避免賦值時產生錯誤
	for (int i = 0; i < mat.row_len; i++) {
		ma->Datas.push_back(d);
		for (int j = 0; j < mat.col_len; j++) {
			ma->Datas[i].push_back(0);
		}
	}
	//資訊複製
	ma->row_len = mat.row_len;
	ma->col_len = mat.col_len;
	ma->Name = mat.Name;
	ma->index = mat.index;

	for (int i = 0; i < mat.row_len; i++) {
		for (int j = 0; j < col_len; j++)
			ma->Datas[i][j] = Datas[i][j] + mat.Datas[i][j];
	}
	return ma;
}
Matrix* Matrix::operator+(Matrix *mat) {
	Matrix *ma = new Matrix();
	std::vector<double> d;		//輔助進行加法，因為兩個二維vector容器無法直接進行加減

								//初始化，避免賦值時產生錯誤
	for (int i = 0; i < mat->row_len; i++) {
		ma->Datas.push_back(d);
		for (int j = 0; j < mat->col_len; j++) {
			ma->Datas[i].push_back(0);
		}
	}
	//資訊複製
	ma->row_len = mat->row_len;
	ma->col_len = mat->col_len;
	ma->Name = mat->Name;
	ma->index = mat->index;

	for (int i = 0; i < mat->row_len; i++) {
		for (int j = 0; j < col_len; j++)
			ma->Datas[i][j] = Datas[i][j] + mat->Datas[i][j];
	}

	return ma;
}


//constant multi, needn't to check dimension
Matrix* operator*(double input_constant, Matrix &mat) {
	Matrix *ma = new Matrix();
	std::vector<double> d;

	ma->row_len = mat.row_len;
	ma->col_len = mat.col_len;
	ma->Name = mat.Name;
	ma->index = mat.index;

	for (int i = 0; i < mat.row_len; i++) {
		ma->Datas.push_back(d);
		for (int j = 0; j < mat.col_len; j++) {
			ma->Datas[i].push_back(input_constant * mat.Datas[i][j]);
		}
	}
	return ma;
}

//const multi(in vector<vector<double>> form)
Matrix* operator*(std::vector<std::vector<double>> &input_constant, Matrix &mat) {
	std::vector<double> d;
	Matrix *ma = new Matrix();

	ma->row_len = mat.row_len;
	ma->col_len = mat.col_len;
	ma->Name = mat.Name;
	ma->index = mat.index;


	std::cout << input_constant[0][0] << " ";
	for (int i = 0; i < mat.row_len; i++) {
		ma->Datas.push_back(d);
		for (int j = 0; j < mat.col_len; j++) {
			ma->Datas[i].push_back(input_constant[0][0] * mat.Datas[i][j]);
		}
	}
	return ma;
}

Matrix operator*(double input_constant, std::vector<std::vector<double>> &mat) {
	Matrix ma;
	std::vector<double> d;

	for (int i = 0; i < mat.size(); i++) {
		ma.Datas.push_back(d);
		for (int j = 0; j < mat[0].size(); j++) {
			ma.Datas[i].push_back(input_constant * mat[i][j]);
		}
	}
	return ma;
}

//new
std::vector<std::vector<double>> operator+(std::vector<std::vector<double>> m1, std::vector<std::vector<double>> m2) {
	for (int i = 0; i < m1.size(); i++) {
		for (int j = 0; j < m1[0].size(); j++) {
			m1[i][j] = m1[i][j] + m2[i][j];
		}
	}
	/*
	std::cout << "m1" << std::endl;
	for (int i = 0; i < m1.size(); i++) {
	for (int j = 0; j < m1[0].size(); j++) {
	std::cout << m1[i][j] << " ";
	}
	std::cout << "\n";
	}*/
	return m1;
}

std::vector<std::vector<double>> operator*(std::vector<std::vector<double>> m, double input_constant) {
	std::vector<std::vector<double>> result;
	std::vector<double> assist;
	for (int i = 0; i < m.size(); i++) {
		result.push_back(assist);
		for (int j = 0; j < m[0].size(); j++) {
			result[i].push_back(input_constant * m[i][j]);
		}
	}
	return result;
}
//caution:no check dimension here
std::vector<std::vector<double>> operator*(std::vector<std::vector<double>> m2, std::vector<std::vector<double>> m1) {
	std::vector<std::vector<double> >result;
	std::vector<double> assist_result;
	//initialize
	for (int i = 0; i < m2.size(); i++) {
		result.push_back(assist_result);
		for (int j = 0; j < m1[0].size(); j++) {
			result[i].push_back(0);
		}
	}
	//processing
	for (int i = 0; i < m2.size(); i++)
		for (int j = 0; j < m1[0].size(); j++)
			for (int k = 0; k < m1.size(); k++) {
				//std::cout << "m1:"<< m1[k][j] << " ";
				//std::cout << "m2:" << m2[i][k] << "\n";
				result[i][j] += m2[i][k] * m1[k][j];
			}
	return result;
}

Matrix* operator*(Matrix &m2, Matrix &m1) {
	std::vector<double> d;
	double save = 0;		//暫存結果
	Matrix *ma = new Matrix();

	ma->row_len = m2.row_len;
	ma->col_len = m2.col_len;
	ma->Name = m2.Name;
	ma->index = m2.index;

	//initialize
	for (int i = 0; i < m2.Datas.size(); i++) {
		ma->Datas.push_back(d);
		for (int j = 0; j < m1.Datas[0].size(); j++) {
			ma->Datas[i].push_back(0);
		}
	}

	//processing

	for (int i = 0; i < m2.Datas.size(); i++)
		for (int j = 0; j < m1.Datas[0].size(); j++)
			for (int k = 0; k < m1.Datas.size(); k++)
				ma->Datas[i][j] += m2.Datas[i][k] * m1.Datas[k][j];


	//test data
	/*
	std::cout << "multi(vector<vector>,Matrix)" << std::endl;
	for (int i = 0; i < ma->Datas.size(); i++) {
	for (int j = 0; j < ma->Datas[0].size(); j++) {
	std::cout << ma->Datas[i][j] << " ";
	}
	std::cout << "\n";
	}
	*/
	return ma;
}