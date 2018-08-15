#include"DataManager.h"
#include"Mvector.h"			//overloading operator for vector

//�H�U���Ooperator overloading(���ǶǤJ���A���P)
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
	std::vector<double> d;		//���U�i��[�k�A�]����ӤG��vector�e���L�k�����i��[��

	std::cout << row_len << std::endl;


	//��l�ơA�קK��Ȯɲ��Ϳ��~
	for (int i = 0; i < mat.row_len; i++) {
		ma->Datas.push_back(d);
		for (int j = 0; j < mat.col_len; j++) {
			ma->Datas[i].push_back(0);
		}
	}
	//��T�ƻs
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
	std::vector<double> d;		//���U�i��[�k�A�]����ӤG��vector�e���L�k�����i��[��



	//��l�ơA�קK��Ȯɲ��Ϳ��~
	for (int i = 0; i < mat->row_len; i++) {
		ma->Datas.push_back(d);
		for (int j = 0; j < mat->col_len; j++) {
			ma->Datas[i].push_back(0);
		}
	}
	//��T�ƻs
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
	double save = 0;		//�Ȧs���G
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
	//�}���ɮסA�ǤJopen��ƪ��ѼƦ���ӡA���}�_���ɮצW�١A�}���ɮת��Ҧ��Ѽ�(�o��std::ios::in��Ū��(��J)���A)
	fin.open(FileName, std::ios::in);
	//Ū�����Ѧ^��false
	if (!fin)
	{
		return false;
	}
	else
	{
		//�аO��eŪ���V�qID
		int currentLoadVectorID = 0;
		//�w�q�V�q��ƼȦs�ܼ�
		std::vector<double> tempVectorData;
		//�w�qŪ���ɮצr��Ȧs�ܼ�
		std::string tempSring;
		//�q�ɮ�Ū���r��A�ѪR���V�q�`��
		fin >> tempSring;		
		//����Ū�ɰj��A�æbŪ���ɮ׵����ɵ���
		while (!fin.eof())
		{
			//�q�ɮ�Ū���r��
			fin >> tempSring;
			//�ѪR��V�q�аO"V"
			if (tempSring == "V")
			{
				if (currentLoadVectorID != 0)
				{
					//�w�q�Ȧs�V�q��Ƶ��c
					Vector tempVector;
					//�s�J�V�q���
					tempVector.Data = tempVectorData;
					tempVector.len = tempVectorData.size();
					//�w�q�V�q�ܼƦW�١A��VectorVariableIndex�ܼƧ@�W�٪�����
					tempVector.index = std::to_string(VectorVariableIndex);
					std::string vectorVariableTemp = "v" + std::to_string(VectorVariableIndex);
					//�s�J�V�q�ܼƦW��
					tempVector.Name = vectorVariableTemp;
					//�s�J�V�q
					Vectors.push_back(tempVector);
					//���WVectorVariableIndex�A�H�T�O�ܼƦW�٤�����
					VectorVariableIndex++;
					//�M���V�q��ƼȦs
					tempVectorData.clear();
				}
				//���WcurrentLoadVectorID�A�аO���eŪ���V�qID
				currentLoadVectorID++;
				//�q�ɮ�Ū���r��A�ѪR���V�q����
				fin >> tempSring;
			}
			else
			{
				//Ū���V�q��ơA�ñNstring�ରdouble
				double value;
				value = (double)strtod(tempSring.c_str(), NULL);
				//�N�V�q��Ʀs�J�Ȧs
				tempVectorData.push_back(value);
			}
		}
		
		//Ū�J��J�ɮפ��̫�@�ӦV�q��T
		Vector tempVector;
		//tempVectorData.pop_back();
		tempVector.Data = tempVectorData;
		tempVector.len = tempVectorData.size();
		std::string vectorVariableTemp = "v" + std::to_string(VectorVariableIndex);
		tempVector.Name = vectorVariableTemp;
		Vectors.push_back(tempVector);
		VectorVariableIndex++;
		//Ū�����\�^��false
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
	//�}���ɮסA�ǤJopen��ƪ��ѼƦ���ӡA���}�_���ɮצW�١A�}���ɮת��Ҧ��Ѽ�(�o��std::ios::in��Ū��(��J)���A)
	fin.open(MfileName, std::ios::in);
	//Ū�����Ѧ^��false
	if (!fin)
		return false;
	else {
		int row = 0, col = 0, cnt = 0;
		//�аO��eŪ���x�}ID
		int currentLoadmatrixID = 0;
		//�w�q�x�}��ƼȦs�ܼ�
		std::vector<double> temp_row;
		std::vector<std::vector<double>> tempMatrixData;
		//�w�qŪ���ɮצr��Ȧs�ܼ�
		std::string tempSring;
		//�q�ɮ�Ū���r��A�ѪR���x�}�`��
		fin >> tempSring;

		//����Ū�ɰj��A�æbŪ���ɮ׵����ɵ���
		while (!fin.eof()) {
			//�q�ɮ�Ū���r��
			fin >> tempSring;
			//�ѪR��x�}�аO"M"
			if (tempSring == "M")
			{
				if (currentLoadmatrixID != 0)
				{
					//�w�q�Ȧs�V�q��Ƶ��c
					Matrix tempMatrix;
					//�s�J�V�q���
					tempMatrix.Datas = tempMatrixData;
					//�w�q�V�q�ܼƦW�١A��VectorVariableIndex�ܼƧ@�W�٪�����
					tempMatrix.index = std::to_string(MatrixVariableIndex);
					std::string matrixVariableTemp = "m" + std::to_string(MatrixVariableIndex);
					//�s�J�V�q�ܼƦW��
					tempMatrix.Name = matrixVariableTemp;
					//�s�Jrow & column
					tempMatrix.row_len = row;
					tempMatrix.col_len = col;
					//�s�J�V�q
					matrix.push_back(tempMatrix);
					//���WVectorVariableIndex�A�H�T�O�ܼƦW�٤�����
					MatrixVariableIndex++;
					//�M���V�q��ƼȦs;
					tempMatrixData.clear();
					row = 0; col = 0; cnt = 0;
				}
				//���WcurrentLoadVectorID�A�аO���eŪ���x�}ID
				currentLoadmatrixID++;
				//�q�ɮ�Ū���r��A�ѪR���V�q����
				fin >> tempSring;
				row = std::stoi(tempSring);
				fin >> tempSring;
				col = std::stoi(tempSring);
			}
			else {
				//Ū���x�}��ơA�ñNstring�ରdouble
				double value;
				value = (double)strtod(tempSring.c_str(), NULL);
				//�N�x�}��Ʀs�J�Ȧs
				temp_row.push_back(value);
				cnt++;
				if (cnt == col) {
					tempMatrixData.push_back(temp_row);
					temp_row.clear();
					cnt = 0;
				}
			}
		}
		//Ū�J��J�ɮפ��̫�@�ӦV�q��T
		Matrix tempMatrix;
		tempMatrix.Datas = tempMatrixData;
		std::string vectorVariableTemp = "m" + std::to_string(MatrixVariableIndex);
		tempMatrix.Name = vectorVariableTemp;
		tempMatrix.row_len = row;
		tempMatrix.col_len = col;
		matrix.push_back(tempMatrix);
		MatrixVariableIndex++;
		//Ū�����\�^��false
		return true;
	}
}
std::vector<Matrix> DataManager::getMatrix() { return matrix; }
void DataManager::setName(std::string name) { MfileName = name; }
