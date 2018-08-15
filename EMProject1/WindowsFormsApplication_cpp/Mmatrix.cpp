#include "Mmatrix.h"

Mmatrix::Mmatrix() {
}

bool Mmatrix::Deal_outstack(std::vector<std::string> &out, std::vector<Matrix> &matrix,
	std::vector<int> &rec) {
	std::vector<std::string> tempOut;
	std::stack<std::string> Opstack;
	std::vector<int> index;
	bool check = false;
	for (int i = 0; i < out.size(); i++) {
		if (this->search_menu(out[i])) {
			Opstack.push(out[i]);
			check = true;
		}
		else if (this->search_func(out[i])) {
			Opstack.push(out[i]);
			check = true;
		}
		else if (out[i] == ")")
			for (int j = 0; !Opstack.empty(); j++)
				if (Opstack.top() == "(") {
					Opstack.pop();
					break;
				}
				else {
					tempOut.push_back(Opstack.top());
					index.push_back(-1);
					Opstack.pop();
				}
		else {
			for (unsigned int j = 0; j < matrix.size(); j++)
				if (out[i] == matrix[j].Name.c_str()) {
					tempOut.push_back(matrix[j].Name);
					index.push_back(1);
					check = true;
					break;
				}
			if (this->isDigit(out[i])) {
				check = true;
				tempOut.push_back(out[i]);
				index.push_back(0);
			}
		}
	}
	if (!Opstack.empty())
		while (!Opstack.empty()) {
			tempOut.push_back(Opstack.top());
			index.push_back(-1);
			Opstack.pop();
		}
	if (!check)
		return false;
	else {
		out = tempOut;
		rec = index;
		return true;
	}
}
int Mmatrix::rank(std::vector < std::vector < double >> &matrix) {
	double check = 0;
	std::cout << "Utri " << "\n";
	std::vector<std::vector<double>> output = Utri(matrix);
	std::cout << "rank:" << "\n";
	int count = 0;
	for (int i = 0; i< matrix.size(); i++) {
		int len = this->n_0_len(output[i]);
		if (len > 0) count++;
	}
	return count;
}
std::vector<std::vector<double>> Mmatrix::Trsp(Matrix &m) {
	std::vector<double> r;
	std::vector<std::vector<double>> temp;
	for (int i = 0; i < m.col_len; i++) {
		for (int j = 0; j < m.row_len; j++) {
			r.push_back(m.Datas[j][i]);
		}
		temp.push_back(r);
		r.clear();
	}
	return temp;
}

std::vector<double> Mmatrix::solve_linear_s(std::vector<std::vector<double>> &variable, std::vector<std::vector<double>> &agument) {
	std::vector<std::vector<double>> output;
	std::vector<double> s;
	output = Utri(variable);
	std::cout << output.size() << std::endl;
	s.push_back(determinant(output));


	int ii = 0, jj = 0;		//count
	for (int i = 0; i < output[0].size(); i++) {
		for (int j = 0; j < output.size(); j++) {
			variable[j][i] = agument[j][ii];
		}
		output = Utri(variable);
		s.push_back(determinant(output));
	}

	for (int i = 0; i < s.size(); i++) {
		std::cout << s.at(i) << std::endl;
	}

	return s;
}

std::vector<std::vector<double>> Mmatrix::muti(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b) {
	std::vector<std::vector<double>> temp;
	std::vector<double> t;
	if (!this->C_size(a, b, 1)) return temp;
	double sum = 0.0;

	for (int i = 0; i < a.size(); i++) {
		for (int l = 0; l < b[0].size(); l++) {
			for (int j = 0; j < a[0].size(); j++) {
				sum += a[i][j] * b[j][l];
			}
			t.push_back(sum);
			sum = 0.0;
		}
		temp.push_back(t);
		t.clear();
	}
	return temp;
}

std::vector<std::vector<double>> Mmatrix::inverse(std::vector<std::vector<double>> &v, bool &jud) {
	std::vector<std::vector<double>> t = v;
	std::vector<std::vector<double>> temp;
	if (t.size() != t[0].size()) {
		jud = false;
		return temp;
	}
	//單位矩陣
	for (int i = 0; i < t.size(); i++) {
		std::vector<double> a(t[0].size(),0);
		temp.push_back(a);
	}
	for (int i = 0; i < t.size(); i++) {
		temp[i][i] = 1;
	}
	//*******************上三角****************************
	int start = 0;
	//排列順序改變
	Ach(t,temp,1);
	for (int i = 0; i < t.size()&&R_n0_len(t[i])!=0; i++) {
		if (t[i][start] == 0.0) {
			start = 1;
			while (t[i][start] == 0.0 && start<t[0].size())
				start++;
		}
		std::cout << i << "\n";
		for (int j = i + 1; j < t.size()&&start <t[0].size(); j++) {
			double k = t[j][start] / t[i][start];
			int len = n_0_len(t[j]);
			for (int l = start, r = 0; r < len && fabs(k)>pow(10,-8); l++, r++) {
				//原本矩陣操作
				double te = t[j][l] - k * t[i][l];
				if (fabs(te) < pow(10, -8)) t[j][l] = 0.0;
				else t[j][l] = te;
				//形成反矩陣
				double te2 = temp[j][l] - k*temp[i][l];
				if(fabs(te2) < pow(10, -8)) temp[j][l] = 0.0;
				else temp[j][l] = te2;
			}
		}
		Ach(t,temp,1);
		for (int x = 0; x < t.size(); x++) {
			for (int y = 0; y < t[0].size(); y++) {
				std::cout << t[x][y]<<" ";
			}
			std::cout <<"\n";
		}
		start = 0;
	}
	//************************向上消***********************************
	//若沒反矩陣
	bool next = true;
	for (int i = 0; i < v.size(); i++)
		if (t[i][i] == 0) next = false;

	if (next) {
		int Len = t[0].size();
		start = Len - 1;
		for (int i = t.size() - 1; i > 0; i--) {
				if (t[i][start] == 0.0) {
					while (t[i][start] == 0.0 && start>0)
						start--;
				}
			for (int j = i - 1; j >= 0&&start>0; j--) {
				double k = t[j][start] / t[i][start];
				int len = R_n0_len(t[j]);
				for (int l = start, r = 0; r < len && k != 0; l--, r++) {
					double te = t[j][l] - k * t[i][l];
					if (fabs(te) < pow(10, -8)) t[j][l] = 0.0;
					else t[j][l] = te;
					//形成反矩陣
					double te2 = temp[j][l] - k* temp[i][l];
					if (fabs(te2) < pow(10, -8)) temp[j][l] = 0.0;
					else temp[j][l] = te2;
				}
			}
			start = Len - 1;
		}
	}
	for (int i = 0; i < t.size(); i++) {
		double k = 1 / t[i][i];
		t[i][i] = 1;
		for (int j = 0; j < t[0].size(); j++)
			temp[i][j] *= k;
	}
	std::cout << "end*******************\n\n";
	for (int i = 0; i < t.size() && next; i++) {
		for (int j = 0; j < t[0].size(); j++)
			std::cout << temp[i][j] << " ";
		std::cout << "\n";
	}
	std::cout << "********************************\n\n";
	for (int i = 0; i < t.size() && next; i++) {
		for (int j = 0; j < t[0].size(); j++)
			std::cout << t[i][j] << " ";
		std::cout << "\n";
	}
	if (!next)
		std::cout << "\n\nerror" << "\n";
	else std::cout << "\n\n good \ >o< /" << "\n\n";
	jud = next;

	return temp;
}
//矩陣的排列
void Mmatrix::Ach(std::vector<std::vector<double>> &v, std::vector<std::vector<double>> &s,int cho) {
	if (cho == 1) {
		for (int i = 0; i < v.size(); i++) {
			//revise at here
			for (int j = i; j < v.size(); j++) {
				if (n_0_len(v[j]) > n_0_len(v[i])) {
					std::vector<double> t = v[j];
					v[j] = v[i];
					v[i] = t;
					t.clear();
					t = s[j];
					s[j] = s[i];
					s[i] = t;
				}
			}
		}
	}
}
int Mmatrix::search_Ms(std::string s,std::vector<Matrix> &datas) {
	int i = 0;
	while (i < datas.size()) {
		if (s == datas[i].Name);
			return i;
		i++;
	}
	return -1;
}
bool Mmatrix::search_menu(std::string in) {
	std::string menu[6] = { "(", ",", "+", "-", ",","*" };
	for (int i = 0; i < 6; i++)
		if (menu[i] == in)
			return true;
	return false;
}
bool Mmatrix::isDigit(std::string s) {
	for (int i = 0; i < s.size(); i++)
		if (s.at(i)<'0' || s.at(i)>'9')
			return false;
	return true;
}
bool Mmatrix::C_size(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b,int jud) {
	switch (jud) {
	case 0:
		if (a.size() == b.size() && a[0].size() == b[0].size())
			return true;
		break;
	case 1:
		if (a[0].size() == b.size())
			return true;
		break;
	default:
		return false;
		break;
	}
	return false;
}
//排序
void Mmatrix::chL(std::vector<std::vector<double> > &v) {
	for (int i = v.size()-1;i >=0; i--) {
		//revise at here
		for (int j = i; j >=0; j--) {
			if (R_n0_len(v[j]) > R_n0_len(v[i])) {
				std::vector<double> t = v[i];
				v[i] = v[j];
				v[j] = t;
			}
		}
	}
}
void Mmatrix::ch(std::vector<std::vector<double> > &v) {
	for (int i = 0; i < v.size(); i++) {
		//revise at here
		for (int j = i; j < v.size() ; j++) {
			if (n_0_len(v[j]) > n_0_len(v[i])) {
				std::vector<double> t = v[j];
				v[j] = v[i];
				v[i] = t;
			}
		}
	}
}
//找長度
int Mmatrix::n_0_len(std::vector<double> &v) {
	for (int i = 0; i < v.size(); i++)
		if (v[i] != 0.0)
			return v.size() - i;
	return 0;
}
int Mmatrix::R_n0_len(std::vector<double> &v) {
	for (int i = v.size()-1; i >= 0; i--)
		if (fabs(v[i]) > 0.0000001)
			return i+1;
	return 0;
}
//上下 三角
std::vector<std::vector<double>> Mmatrix::Utri(std::vector<std::vector<double> > &v) {
	std::vector<std::vector<double> > t = v;
	int start = 0;
	ch(t);
	for (int i = 0; i < t.size(); i++) {
		if (t[i][start] == 0.0) {
			start = 1;
			while (t[i][start] == 0.0 && start < t[0].size()) {
				start++;
				if (start >= t[0].size()) break;
			}
		}
		for (int j = i + 1; j < t.size()&&start<t[0].size(); j++) {
			double k = t[j][start] / t[i][start];
			int len = n_0_len(t[j]);
			for (int l = start, r = 0; r < len && k != 0; l++, r++) {
				double te = t[j][l] - k * t[i][l];
				if (fabs(te) < pow(10,-8)) t[j][l] = 0.0;
				else t[j][l] = te;
			}
		}
		ch(t);
		start = 0;
	}
	return t;
}	
std::vector<std::vector<double>> Mmatrix::Ltri(std::vector<std::vector<double> > &v) {
	std::vector<std::vector<double> > t = v;
	int Len = t[0].size();
	int start = Len- 1;
	chL(t);
	for (int i = t.size()-1; i >= 0; i--) {
		if (t[i][start] == 0.0) {
			while (t[i][start] == 0.0 && start >0)
				start--;
		}
		for (int j = i - 1; j >=0&&start>0; j--) {
			double k = t[j][start] / t[i][start];
			int len = R_n0_len(t[j]);
			for (int l = start, r = 0; r < len && k!=0; l--, r++) {
				double te = t[j][l] - k * t[i][l];
				if (fabs(te) < pow(10, -8)) t[j][l] = 0.0;
				else t[j][l] = te;
			}
		}
		chL(t);	
		start = Len-1;
	}
	return t;
}
Matrix Mmatrix::vectorToM(std::vector<std::vector<double> > vector) {
	Matrix m;
	m.Datas = vector;
	return m;
}

//Lu decomposition
//don't use it to print out the result directly
//PA=LU
Matrix Mmatrix::lu_decomposition(std::vector<std::vector<double> > &m) {
	Matrix result;
	std::vector<double> assist_L;
	std::vector<double> assist_U;
	int dim = m.size();
	double rmax;			//maximum value of row
	int row_pos;			//max value position of rmax
	int* pivot = new int[dim];		//if pivot[i] = k -> m[i][k] = 1;

	for (int i = 0; i < dim; i++) pivot[i] = i;	//exchange order of each vector(i = row)

	for (int r = 0; r < dim; r++) {

		//find the ith row(maximum)
		rmax = 0;
		for (int i = r; i < dim; i++) {
			//else if (fabs(m[i][r] - rmax) < 0.00000001) throw "singular matrix";  //singular matrix
			if (m[i][r] != 0 && fabs(m[i][r] - rmax) > 0.00000001) {
				//if (fabs(m[i][r]) > fabs(rmax)) { //consider the max/min value always first, then use this.
				rmax = m[i][r];
				row_pos = i;
				break;
				//if (fabs(m[i][r] - rmax) > 0.00000001)    break;
			}
		}


		//exchange info between pivot (order of matrix elimination)
		//that is, recording the order during we processing the matrix.
		int temp_order = pivot[r];
		pivot[r] = pivot[row_pos];
		pivot[row_pos] = temp_order;

		//exchange line
		m.at(row_pos).swap(m.at(r));

		//lu decomposition by gauss elimination
		for (int i = r + 1; i < dim; i++) {
			m[i][r] /= m[r][r];
			for (int j = r + 1; j < dim; j++) {
				m[i][j] -= m[i][r] * m[r][j];
			}
		}
	}

	//seperate L and U
	for (int row = 0; row < dim; row++) {
		result.L.push_back(assist_L);
		result.U.push_back(assist_U);
		for (int col = 0; col < dim; col++) {

			result.L[row].push_back(0);
			result.U[row].push_back(0);
			//determine whether it's U or L, or on diagonal
			if (row > col) {
				result.U[row][col] = 0;
				result.L[row][col] = m[row][col];
			}
			else if (row < col) {
				result.L[row][col] = 0;
				result.U[row][col] = m[row][col];
			}
			else {
				result.U[row][col] = m[row][col];
				result.L[row][col] = 1;
			}
		}
	}

	assist_L.clear();
	//copy p
	for (int i = 0; i < dim; i++) {
		result.P.push_back(assist_L);
		for (int j = 0; j < dim; j++) {
			result.P[i].push_back(0);
		}
	}
	for (int i = 0; i < dim; i++)	result.P[i][pivot[i]] = 1;


	//print out result for test
	for (int i = 0; i < result.U.size(); i++) {
		for (int j = 0; j <result.U[0].size(); j++) {
			std::cout << result.U[i][j] << " ";
		}
		std::cout << "\n";
	}

	delete pivot;
	return result;
}

//solve linear equation by back substitution
std::vector<double> Mmatrix::solvelineareq(Matrix &m, std::vector<std::vector<double> > &augment) {
	Matrix temp = m;
	//get L, U , and P(maybe useless in the example)
	Matrix lup = lu_decomposition(temp.Datas);
	int dim = lup.U.size();				//range = uppertri/lowertri.size
	double* y = new double[dim];
	double* x = new double[dim];
	std::vector<double> result;

	//here we don't consider the order
	//Ly=b (forward substution)
	for (int i = 0; i < dim; i++) {
		y[i] = augment[i][0];
		for (int j = 0; j < i; j++)
			y[i] -= lup.L[i][j] * y[j];
	}

	//Ux=y (back substution)
	for (int i = dim - 1; i >= 0; i--) {
		x[i] = 0;
		//iterative substract y 
		for (int j = dim - 1; j > i; j--)
			y[i] -= lup.U[i][j] * x[j];
		//devide coefficient of x(get answer of x[i])
		x[i] = y[i] / lup.U[i][i];
	}

	for (int i = 0; i < dim; i++)
		result.push_back(x[i]);

	std::cout << "result:" << std::endl;
	std::cout << result.size() << std::endl;

	for (int i = 0; i < dim ; i++) {
	std::cout << result.at(i) << "\n";
	}
	std::cout << "\n";
	
	delete y;
	delete x;
	return result;
}

//after get linear system's solve, do that : AX=I 
//and then seperate X by column [x1 x2 x3... xn]
//calculate each row vector's inverse
//finally, row vector combination = inverse matrix
//caution: be careful for the "difference" when using this to calculate the other function.
std::vector<std::vector<double> > Mmatrix::InverseM(Matrix &m) {
	std::vector<std::vector<double> > inverse;
	std::vector<std::vector<double> > iden;
	std::vector<std::vector<double> > iden_row;
	std::vector<double> x(m.Datas.size(), 0);
	std::vector<double> assist_iden;

	if (m.Datas.size() != m.Datas[0].size()) throw "not square matrix";

	//initialize
	for (int i = 0; i < m.Datas.size(); i++) {
		iden.push_back(assist_iden);
		inverse.push_back(assist_iden);
		for (int j = 0; j < m.Datas.size(); j++) {
			iden[i].push_back(0);
			inverse[i].push_back(0);
		}
	}

	//turn 30x30 -> 30x1
	assist_iden.clear();
	for (int i = 0; i < m.Datas.size(); i++) {
		iden_row.push_back(assist_iden);
		for (int j = 0; j < 1; j++)
			iden_row[i].push_back(0);
	}

	//processing
	for (int i = 0; i < m.Datas.size(); i++) {
		iden[i][i] = 1;

		for (int j = 0; j < m.Datas.size(); j++) {
			iden_row[j][0] = iden[j][i];
		}

		//since we solve it column(row vector) by column, before we move to the next column
		//we need to clear it
		//if (i != 0)iden[i - 1][i] = 0;
		if (i != 0) iden_row[i - 1][0] = 0;

		//load solution
		x = solvelineareq(m, iden_row);

		//solution correspond to the each row vector
		for (int j = 0; j < inverse.size(); j++)
			inverse[j][i] = x[j];
	}

	/*
	std::cout << "inverse:" << inverse.size() << std::endl;
	for (int i = 0; i < inverse.size(); i++) {
	for (int j = 0; j < inverse[0].size(); j++) {
	std::cout << inverse[i][j] << " ";
	}
	std::cout << "\n";
	}*/

	return inverse;
}

//approximating eigenvalue
//at most 3D
void Mmatrix::power_method(Matrix &A,std::vector<double> *fetch) {
	//set initial guess
	std::vector<std::vector<double> > X;
	std::vector<std::vector<double> > temp;
	std::vector<double> assist;

	//initialize(row vector), set X = 1
	for (int i = 0; i < A.Datas.size(); i++) {
		X.push_back(assist);
		for (int j = 0; j < 1; j++) {
			if (i == 0)X[i].push_back(1);
			else X[i].push_back(0);
		}

	}
	//initialize result to save
	for (int i = 0; i < A.Datas.size(); i++) {
		A.eigen.e_vector.push_back(assist);
		for (int j = 0; j < 1; j++)
			A.eigen.e_vector[i].push_back(0);

	}
	//initialize temp
	for (int i = 0; i < A.Datas.size(); i++) {
		temp.push_back(assist);
		for (int j = 0; j < 1; j++)
			temp[i].push_back(0);
	}

	//for (int i = 0; i < 20; i++)
	//A.eigen.e_value.push_back(0);


	//iterate result, here we approximating by 6 times
	for (int t = 0; t < 200; t++) {
		temp = A.Datas * X;

		//Rayleigh quotient
		A.eigen.e_value.clear();
		A.eigen.e_value.push_back((scalar(Trsp(temp) * X) / scalar(Trsp(X) * X)));

		//find maximum
		double scale = temp[0][0];
		for (int i = 1; i < temp.size(); i++)
			scale = (fabs(temp[i][0]) > fabs(scale)) ? temp[i][0] : scale;

		//scaling by maximum
		for (int i = 0; i < temp.size(); i++)
			temp[i][0] /= scale;

		//save result(not finished yet)
		std::cout << "eigenvalue:" << std::endl;
		for (int i = 0; i < A.eigen.e_value.size(); i++) {
			std::cout << A.eigen.e_value[i] << " ";
		}
		std::cout << "\n";


		for (int i = 0; i < temp.size(); i++) {
			X[i][0] = temp[i][0];
		}
	}
	for (int i = 0; i < A.eigen.e_value.size(); i++) {
		fetch->push_back(A.eigen.e_value[i]);
		std::cout << A.eigen.e_value[i] << " ";
	}
	std::cout << "\n";
}

std::vector<std::vector<double> > Mmatrix::adjoint(Matrix &m) {
	Matrix result = determinant(m.Datas) * InverseM(m);
	
	std::cout << "\n" << "adjoint" << std::endl;
	for (int i = 0; i < result.Datas.size(); i++) {
		for (int j = 0; j < result.Datas[0].size(); j++) {
			std::cout << result.Datas[i][j] << " ";
		}
		std::cout << "\n";
	}

	return result.Datas;
}

//transpose for 2d vector
std::vector<std::vector<double>> Mmatrix::Trsp(std::vector<std::vector<double> > &m) {
	std::vector<double> r;
	std::vector<std::vector<double>> temp;

	//std::cout << "m:" << m.size() << "x" << m[0].size() << std::endl;
	for (int i = 0; i < m[0].size(); i++) {
		for (int j = 0; j < m.size(); j++) {
			r.push_back(m[j][i]);
		}
		temp.push_back(r);
		r.clear();
	}
	return temp;
}

inline double  Mmatrix::scalar(std::vector<std::vector<double> > &m) {
	return ((m.size() == 1) && (m[0].size() == 1)) ? m[0][0] : throw "cannot convert to scalar";
}

inline std::vector<std::vector<double>>  Mmatrix::MTovector(Matrix m) {
	return m.Datas;
}

inline std::vector<std::vector<double>>  Mmatrix::MTovector(Matrix *m) {
	return m->Datas;
}

std::vector<std::vector<double> > Mmatrix::least_suqare(Matrix &A, Matrix &Y) {
	Matrix X;
	std::vector<std::vector<double> > a;

	X.Datas = Trsp(A) * A.Datas;
	//a = Trsp(A) * Y.Datas;
	a = InverseM(X) * Trsp(A) * Y.Datas;
	/*
	std::cout << "a:" << std::endl;
	std::cout << a.size() << std::endl;
	std::cout << a[0].size() << std::endl;
	for (int i = 0; i < a.size(); i++) {
	for (int j = 0; j < a[0].size(); j++) {
	std::cout << a[i][j] << " ";
	}
	std::cout << "\n";
	}*/

	return a;
}
//尋找運算符號
bool Mmatrix::to_find(System::String ^in) {
	//Tri後(含Tri)沒用到
	std::string menu[17] = { "(", ",", "+","-",",", "*",
		"Trsp", "Utri" ,"Ltri", "inverse", "S_sol",
		"det", "rank", "adj","L_s", "power","LU" };

	std::string temp = "";
	MarshalString(in, temp);
	for (int i = 0; i < 17; i++)
		if (temp == menu[i].c_str()) {
			return true;
		}
	return false;
}
//search function
bool Mmatrix::search_func(std::string func) {
	//same as above
	std::string menu[12] = { "Trsp", "Utri" ,"Ltri", "inverse", "S_sol",
		"det", "Tri", "rank","adj", "L_s","power", "LU" };
	for (int i = 0; i < 12; i++)
		if (func == menu[i])
			return true;
	return false;
}
//find_function    1是向量資料 -1是運算符號 0是純量
void Mmatrix::FindFunc(std::vector<Matrix> &m, std::string func,Ans *ans) {
	if (func == "muti") {
		ans->vvAns = this->muti(m[0].Datas,m[1].Datas);
		ans->jud = 1;
	}
	else if (func == "det") {
		std::vector<std::vector<double>> tri = this->Utri(m[0].Datas);
		ans->dAns = this->determinant(tri);
		ans->jud = 0;
	}
	else if (func == "rank") {
		std::cout << func<<"\n";
		ans->dAns = this->rank(m[0].Datas);
		ans->jud = 0;
		std::cout << ans->dAns << "\n";
	}
	else if (func == "Trsp") {
		ans->vvAns = this->Trsp(m[0]);
		ans->jud = 1;
	}
	else if (func == "Utri") {
		ans->vvAns = this->Utri(m[0].Datas);
		ans->jud = 1;
	}
	else if (func == "Ltri") {
		ans->vvAns = this->Ltri(m[0].Datas);
		ans->jud = 1;
	}
	else if (func == "inverse") {
		if (m[0].col_len == m[0].row_len) {
			ans->vvAns = this->InverseM(m[0]);
			if (!ans->vvAns.empty())
				ans->jud = 1;
			else ans->jud = -1000;
		}
		else ans->jud = -1000;
	}
	else if (func == "S_sol") {
		for (int j = 0; j < m.size(); j++) {
			std::cout << "m name: " << m[j].Name << "\n";
		}
		ans->vAns = this->solvelineareq(m[1], m[0].Datas);
		ans->jud = 2;
	}
	else if (func=="adj") {
		ans->vvAns = this->adjoint(m[0]);
		if (!ans->vvAns.empty())
			ans->jud = 1;
		else ans->jud = -1000;
	}
	else if (func=="L_s") {
		if (m[1].col_len == m[0].row_len) {
			ans->vvAns = this->least_suqare(m[1], m[0]);
			if (!ans->vvAns.empty()) ans->jud = 1;
			else ans->jud = -1000;
		}
		else ans->jud = -1000;
	}
	else if (func=="power") {
		this->power_method(m.back(),&ans->vAns);
		if (!ans->vAns.empty()) ans->jud = 2;
		else ans->jud = -1000;
	}
	else if (func == "LU") {
		std::vector<std::vector<double>> utr = this->Utri(m.back().Datas);
		ans->Mans = this->lu_decomposition(m.back().Datas);
		ans->jud = 3;
		ans->Mans.U = utr;
	}
}

//經過上三角處理後，只要將對角數值作乘積即可得到答案
double Mmatrix::determinant(std::vector <std::vector<double>> &m) {
	double result = 1;	//initialize
	Matrix temp = lu_decomposition(m);
		
	for (int i = 0; i < m.size(); i++) {
		result *= temp.U.at(i).at(i);
		if (i >= 1)
			if (fabs(result - temp.U.at(i - 1).at(i - 1)) < 0.000000001) throw "determinant doesn't exist";
	}
	return result;
}

//size check matrix
bool Mmatrix::check_size_m(std::deque<Matrix> &datas) {
	double row = datas[0].Datas.size();
	double col = datas[0].Datas[0].size();
	if (datas.size() > 1) {
		for (int i = 1; i < datas.size(); i++) {
			if (datas[i].Name == "operator") continue;
			if (row != datas[i].Datas.size() || col != datas[i].Datas[i].size())
				return false;
		}
	}
	return true;
}

//only for dot(remember to check transpose)
inline bool Mmatrix::check_size_m_dot(std::vector<std::vector<double>> &data1, std::vector<std::vector<double>> &data2) {
	return (data1.size() == data2[0].size()) ? true : false;
}

//finding matrix name
void Mmatrix::find_MatrixName(std::vector<Matrix> &matrixs, std::deque<Matrix> &Matrix_stack, std::string &str) {
	for (int i = 0; i < matrixs.size(); i++) {
		if (str == matrixs[i].Name) {
			//data check
			for (int r = 0; r < matrixs[i].row_len; r++) {
				for (int k = 0; k < matrixs[i].col_len; k++) {
					std::cout << matrixs[i].Datas[r][k] << " ";
				}
				std::cout << "\n";
			}
			//data check
			Matrix_stack.push_back(matrixs[i]);
			break;
		}
	}
}

std::string Mmatrix::mathematics_operation(std::vector<int> &index, std::vector<Matrix> &Matrixs, std::vector<std::string> &formula) {
	Ans ans;
	Matrix mat;
	Matrix *ptrmat = &mat;			//預備用
	Matrix const_for_oprate;
	const_for_oprate.Name = "operator";
	double d;		//save constant
	bool sizecheck;
	bool sizecheck_dot=false;
	bool done = false;
	std::deque<Matrix> Matrix_stack;	//save Vector
	std::vector<Matrix> m;

	for (int i = 0; i < formula.size(); i++) {
		//std::cout << Matrixs[0].Name << " " << std::endl;

		std::cout <<formula[i] << " ";
	}
	std::cout <<std::endl;
	for (int i = 0; i < formula.size(); i++) {
		std::cout <<"\n" << formula.at(i) << std::endl;
		if (index.at(i) == -1) {
			if (formula.at(i) == "+") {

				if (!done) {	//檢查維度只需一次
					sizecheck = check_size_m(Matrix_stack);
					if (!sizecheck) return "size doesn't match";
					done = true;
				}

				mat = Matrix_stack.at(Matrix_stack.size() - 2) + Matrix_stack.at(Matrix_stack.size() - 1);

				//運算完成，處理剩餘資源
				for (int i = 0; i < 2; i++) Matrix_stack.pop_back();
				Matrix_stack.push_back(mat);
			}
			else if (formula.at(i) == "-") {
				if (!done) {	//檢查維度只需一次
					sizecheck = check_size_m(Matrix_stack);
					if (!sizecheck) {
						return "size doesn't match";
					}
					done = true;
				}
				mat = Matrix_stack.at(Matrix_stack.size() - 2) + (-1) * Matrix_stack.at(Matrix_stack.size() - 1);
				
				for (int i = 0; i < 2; i++) Matrix_stack.pop_back();
				Matrix_stack.push_back(mat);
			}
			else if (formula.at(i) == "*") {
				if ((Matrix_stack.at(Matrix_stack.size() - 2).Name) == "operator" || (Matrix_stack.at(Matrix_stack.size() - 1).Name) == "operator")
					mat = Matrix_stack.at(Matrix_stack.size() - 2).Datas * Matrix_stack.at(Matrix_stack.size() - 1);
				else {
					if (!sizecheck_dot) {	//乘法維度檢查需不只一次
						sizecheck = check_size_m_dot(Matrix_stack.at(Matrix_stack.size() - 2).Datas, Matrix_stack.at(Matrix_stack.size() - 1).Datas);
						if (!sizecheck) {
							return "size doesn't match";
						}
					}
					mat = Matrix_stack.at(Matrix_stack.size() - 2) * Matrix_stack.at(Matrix_stack.size() - 1);
				}
					
				//mat.Datas = Matrix_stack.at(Matrix_stack.size() - 2).Datas * Matrix_stack.at(Matrix_stack.size() - 1).Datas;

				for (int i = 0; i < 2; i++) Matrix_stack.pop_back();
				Matrix_stack.push_back(mat);
			}
			else if (formula[i] == ",") {
				m.push_back(Matrix_stack.back());
				Matrix_stack.pop_back();
			}
			else if(this->search_func(formula[i])){
				m.push_back(Matrix_stack.back());
				for (int j = 0; j < m.size();j++) {
					std::cout <<"m name: "<< m[j].Name << "\n";
				}
				Matrix_stack.clear();
				this->FindFunc(m, formula[i], &ans);
				std::cout<<"2:" << ans.jud << "\n";
			}
			else {
				Matrix_stack.clear();
				std::cout << "error" << "\n";
			}
		}
		else if (index[i] == 0) {
			std::vector<double> save_double;
			save_double.push_back(std::stod(formula.at(i), NULL));
			const_for_oprate.Datas.push_back(save_double);
			Matrix_stack.push_back(const_for_oprate);
		}
		else if (index.at(i) == 1) {
			find_MatrixName(Matrixs, Matrix_stack, formula.at(i));
		}
		const_for_oprate.Datas.clear();
		mat.Datas.clear();
	}
	//Matrix Mans
	if (!Matrix_stack.empty()) {
		std::cout << Matrix_stack[0].Datas.size() << std::endl;
		ans.vvAns = Matrix_stack[0].Datas;
		ans.jud = 1;
	}
	Matrix_stack.clear();
	return this->Mk_str(ans);
}

std::string Mmatrix::Mk_str(Ans &ans) {
	std::string result = "[ ";
	std::string n=""; 
	MarshalString(System::Environment::NewLine,n);
	if (ans.jud == 0) {
		result += std::to_string(ans.dAns);
	}
	else if (ans.jud == 1) {
		for (int i = 0; i < ans.vvAns.size(); i++) {
			for (int j = 0; j < ans.vvAns[0].size(); j++) {
				result += std::to_string(ans.vvAns[i][j]) + ", ";
			}
			result += n;
		}
		ans.vvAns.clear();
	}
	else if (ans.jud == 2) {
		if (!ans.vAns.empty()) {
			for (int i = 0; i < ans.vAns.size(); i++) {
				result += std::to_string(ans.vAns[i])+n;
			}
		}
		else {
			result += "fault"+n;
		}
		ans.vAns.clear();
	}
	else if (ans.jud== 3) {
		result = " ";
		result += "Upper tri: [" + n;
		for (int i = 0; i < ans.Mans.U.size(); i++) {
			for (int j = 0; j < ans.Mans.U[0].size(); j++) {
				result += std::to_string(ans.Mans.U[i][j])+", ";
			}
			result += n;
		}
		result += "]"+n+" Lower trir: [" +n; 
		for (int i = 0; i < ans.Mans.L.size(); i++) {
			for (int j = 0; j < ans.Mans.L[0].size(); j++) {
				result += std::to_string(ans.Mans.L[i][j]) + ", ";
			}
			result += n;
		}
	}
	else if (ans.jud == -1) {
		if (ans.bAns)
			ans.strAns = "yes";
		else ans.strAns = "no";
		result += ans.strAns;
	}
	else result += "error!!!";
	result += " ]";
	ans.jud = NULL;
	return result;
}
