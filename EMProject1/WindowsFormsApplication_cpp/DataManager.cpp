#include"DataManager.h"
#include"Mvector.h"			//overloading operator for vector

//以下都是operator overloading(有些傳入型態不同)
//for vector(useless now)
inline bool Vector::operator==(Vector& v1) {
	return (Data.size() == v1.Data.size()) ? true : false;
}

Vector Vector::operator=(Vector& vec) {
	for (int i = 0; i<vec.Data.size(); i++) {
		Data.push_back(vec.Data.at(i));
	}
	return *this;
}

Vector Vector::operator+(Vector& vec) {
	Vector ve;

	if (Data.size() == vec.Data.size()) {
		for (int i = 0; i < Data.size(); i++) {
			ve.Data.push_back(Data.at(i) + vec.Data.at(i));
		}
	}
	else {
		std::cout << "size doesn't match" << std::endl;

		getchar();
		exit(1);
	}
	return ve;
}
Vector operator*(double input_constant, Vector &vec) {
	Vector ve;
	for (int i = 0; i < vec.Data.size(); i++) {
		ve.Data.push_back(input_constant * vec.Data.at(i));
	}
	return ve;
}
Vector operator*(std::vector<double> dou, Vector &vec) {
	Vector ve;
	double b = dou.at(dou.size() - 1);
	for (int i = 0; i < vec.Data.size(); i++) {
		ve.Data.push_back(b * vec.Data.at(i));
	}
	return ve;
}
std::vector<double> operator*(double input_constant, std::vector<double> &vec) {
	std::vector<double> ve;
	for (int i = 0; i < vec.size(); i++)
		ve.push_back(input_constant * vec[i]);
	return ve;
}

/*
------------------------------------------------------------
*/

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

	std::cout << "equal" << std::endl;
	for (int i = 0; i < mat->row_len; i++) {
		for (int j = 0; j < mat->col_len; j++) {
			std::cout << Datas[i][j] << " ";
		}
		std::cout << "\n";
	}

	return *this;
}


Matrix* Matrix::operator+(Matrix &mat) {
	Matrix *ma = new Matrix();
	std::vector<double> d;		//輔助進行加法，因為兩個二維vector容器無法直接進行加減

	std::cout << row_len << std::endl;


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

	//test data
	std::cout << "multi(double,2Darray)" << std::endl;
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			std::cout << ma.Datas[i][j] << " ";
		}
		std::cout << "\n";
	}
	//test data

	return ma;

}

std::vector<std::vector<double>> operator*(std::vector<std::vector<double>> &m2, std::vector<std::vector<double>> &m1) {
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
			for (int k = 0; k < m1.size(); k++)
				result[i][j] += m2[i][k] * m1[k][j];

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
	std::cout << "multi(vector<vector>,Matrix)" << std::endl;
	for (int i = 0; i < ma->Datas.size(); i++) {
		for (int j = 0; j < ma->Datas[0].size(); j++) {
			std::cout << ma->Datas[i][j] << " ";
		}
		std::cout << "\n";
	}

	return ma;
}

DataManager::DataManager()
{
	VectorVariableIndex = 0;
	MatrixVariableIndex = 0;
}

bool DataManager::LoadVectorData()
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
		//標記當前讀取向量ID
		int currentLoadVectorID = 0;
		//定義向量資料暫存變數
		std::vector<double> tempVectorData;
		//定義讀取檔案字串暫存變數
		std::string tempSring;
		//從檔案讀取字串，解析掉向量總數
		fin >> tempSring;		
		//執行讀檔迴圈，並在讀到檔案結尾時結束
		while (!fin.eof())
		{
			//從檔案讀取字串
			fin >> tempSring;
			//解析到向量標記"V"
			if (tempSring == "V")
			{
				if (currentLoadVectorID != 0)
				{
					//定義暫存向量資料結構
					Vector tempVector;
					//存入向量資料
					tempVector.Data = tempVectorData;
					tempVector.len = tempVectorData.size();
					//定義向量變數名稱，依VectorVariableIndex變數作名稱的控管
					tempVector.index = std::to_string(VectorVariableIndex);
					std::string vectorVariableTemp = "v" + std::to_string(VectorVariableIndex);
					//存入向量變數名稱
					tempVector.Name = vectorVariableTemp;
					//存入向量
					Vectors.push_back(tempVector);
					//遞增VectorVariableIndex，以確保變數名稱不重複
					VectorVariableIndex++;
					//清除向量資料暫存
					tempVectorData.clear();
				}
				//遞增currentLoadVectorID，標記到當前讀取向量ID
				currentLoadVectorID++;
				//從檔案讀取字串，解析掉向量維度
				fin >> tempSring;
			}
			else
			{
				//讀取向量資料，並將string轉為double
				double value;
				value = (double)strtod(tempSring.c_str(), NULL);
				//將向量資料存入暫存
				tempVectorData.push_back(value);
			}
		}
		
		//讀入輸入檔案中最後一個向量資訊
		Vector tempVector;
		//tempVectorData.pop_back();
		tempVector.Data = tempVectorData;
		tempVector.len = tempVectorData.size();
		std::string vectorVariableTemp = "v" + std::to_string(VectorVariableIndex);
		tempVector.Name = vectorVariableTemp;
		Vectors.push_back(tempVector);
		VectorVariableIndex++;
		//讀取成功回傳false
		return true;
	}
}

std::vector<Vector> DataManager::GetVectors(){
	return Vectors;
}
void DataManager::SetFileName(std::string fileName){
	FileName = fileName;
}
//****************************************************
bool DataManager::LoadData() {
	std::fstream fin;
	//開啟檔案，傳入open函數的參數有兩個，欲開起的檔案名稱，開啟檔案的模式參數(這邊std::ios::in為讀取(輸入)狀態)
	fin.open(MfileName, std::ios::in);
	//讀取失敗回傳false
	if (!fin)
		return false;
	else {
		int row = 0, col = 0, cnt = 0;
		//標記當前讀取矩陣ID
		int currentLoadmatrixID = 0;
		//定義矩陣資料暫存變數
		std::vector<double> temp_row;
		std::vector<std::vector<double>> tempMatrixData;
		//定義讀取檔案字串暫存變數
		std::string tempSring;
		//從檔案讀取字串，解析掉矩陣總數
		fin >> tempSring;

		//執行讀檔迴圈，並在讀到檔案結尾時結束
		while (!fin.eof()) {
			//從檔案讀取字串
			fin >> tempSring;
			//解析到矩陣標記"M"
			if (tempSring == "M")
			{
				if (currentLoadmatrixID != 0)
				{
					//定義暫存向量資料結構
					Matrix tempMatrix;
					//存入向量資料
					tempMatrix.Datas = tempMatrixData;
					//定義向量變數名稱，依VectorVariableIndex變數作名稱的控管
					tempMatrix.index = std::to_string(MatrixVariableIndex);
					std::string matrixVariableTemp = "m" + std::to_string(MatrixVariableIndex);
					//存入向量變數名稱
					tempMatrix.Name = matrixVariableTemp;
					//存入row & column
					tempMatrix.row_len = row;
					tempMatrix.col_len = col;
					//存入向量
					matrix.push_back(tempMatrix);
					//遞增VectorVariableIndex，以確保變數名稱不重複
					MatrixVariableIndex++;
					//清除向量資料暫存;
					tempMatrixData.clear();
					row = 0; col = 0; cnt = 0;
				}
				//遞增currentLoadVectorID，標記到當前讀取矩陣ID
				currentLoadmatrixID++;
				//從檔案讀取字串，解析掉向量維度
				fin >> tempSring;
				row = std::stoi(tempSring);
				fin >> tempSring;
				col = std::stoi(tempSring);
			}
			else {
				//讀取矩陣資料，並將string轉為double
				double value;
				value = (double)strtod(tempSring.c_str(), NULL);
				//將矩陣資料存入暫存
				temp_row.push_back(value);
				cnt++;
				if (cnt == col) {
					tempMatrixData.push_back(temp_row);
					temp_row.clear();
					cnt = 0;
				}
			}
		}
		//讀入輸入檔案中最後一個向量資訊
		Matrix tempMatrix;
		tempMatrix.Datas = tempMatrixData;
		std::string vectorVariableTemp = "m" + std::to_string(MatrixVariableIndex);
		tempMatrix.Name = vectorVariableTemp;
		tempMatrix.row_len = row;
		tempMatrix.col_len = col;
		matrix.push_back(tempMatrix);
		MatrixVariableIndex++;
		//讀取成功回傳false
		return true;
	}
}
std::vector<Matrix> DataManager::getMatrix() { return matrix; }
void DataManager::setName(std::string name) { MfileName = name; }
