#include"Newton.h"
#include"DotNetUtilities.h"

Newton::Newton() {
	mathop = new MathOp();
}
Newton::Newton(std::vector<std::string> f) {
	mathop = new MathOp();
	this->fx = mathop->format(f);
}

void Newton::NewtonAlgorism(double start, double end, double x, double y, std::vector<std::string> f) {
	double h = 0.00001, k = 1;
	Matrix m;
	std::vector<std::vector<double>> X;
	std::vector<std::vector<double>> X_1;
	std::vector<std::vector<double>> hessian;
	std::vector<std::vector<double>> g;
	std::vector<double> assist;
	//2d = 2 
	for (int i = 0; i < 2; i++) {
		X.push_back(assist);
		//X_1.push_back(assist);
		g.push_back(assist);
		for (int j = 0; j < 1; j++) {
			//X_1[i].push_back(0);
			X[i].push_back(0);
		}
	}

	X[0][0] = x;
	X[1][0] = y;
	this->storePs(X, f);
	for (int i = 0; i < 50; i++) {
		//gradient
		g[0].push_back(mathop->fx(X[0][0], X[1][0], h, f));
		g[1].push_back(mathop->fy(X[0][0], X[1][0], k, f));

		//initial F
		hessian = this->hessian_2d(X[0][0], X[1][0], h, k, f);
		this->out += "hessian:" + this->nextline() + "[ ";

		for (int j = 0; j < hessian.size(); j++) {
			for (int l = 0; l<hessian.back().size(); l++) {
				this->out += std::to_string(hessian[j][l]) + " , ";
			}
			this->out += this->nextline();
		}
		this->out += "]" + this->nextline();

		//inv
		hessian = this->inv_hessian(hessian);
		this->out += "inv_hessian:" + this->nextline() + "[ ";
		for (int j = 0; j < hessian.size(); j++) {
			for (int l = 0; l< hessian.back().size(); l++)
				this->out += std::to_string(hessian[j][l]) + " , ";

			this->out += this->nextline();
		}
		this->out += "]" + this->nextline();

		//set iterative condition(be careful of order that mathematic operation, it will effect the answer)
		X_1 = hessian * g * (-1) + X;

		this->out += "x: [" + this->nextline();
		this->out += std::to_string(X_1[0][0]) + " , " + this->nextline() +
			std::to_string(X_1[1][0]) + this->nextline() + "]" + this->nextline();

		//judge if the condition  will be satisfied
		if (fabs(X_1[0][0] + X[0][0] * (-1)) < eps && fabs(X_1[1][0] + X[1][0] * (-1)) < eps) break;
		else X = X_1;
		this->storePs(X, f);

		g[0].clear();
		g[1].clear();
	}

	this->out += "x,x1: [ " + std::to_string(X[0][0]) + " , " + std::to_string(X[1][0]) +
		" ]" + this->nextline();
	this->out += "f(x,y): " + std::to_string(mathop->ComputeFx(X[0][0], X[1][0], f)) +
		this->nextline();

	std::cout << "answer" << "\n";
	std::cout << "x:" << X[0][0] << "\n";
	std::cout << "x1:" << X[1][0] << "\n";
	std::cout << "f(x,y):" << mathop->ComputeFx(X[0][0], X[1][0], f) << "\n";
}

void Newton::NewtonAlgorism(double start, double end, double x, std::vector<std::string> f) {
	double h = 0.00001, k = 1;
	Matrix m;
	std::vector<std::vector<double>> X;
	std::vector<std::vector<double>> X_1;
	std::vector<std::vector<double>> hessian;
	std::vector<std::vector<double>> g;
	std::vector<double> assist;

	//1d = 1 
	for (int i = 0; i < 1; i++) {
		X.push_back(assist);
		//X_1.push_back(assist);
		g.push_back(assist);
		for (int j = 0; j < 1; j++) {
			//X_1[i].push_back(0);
			X[i].push_back(0);
		}
	}

	X[0][0] = x;
	this->storePs(X, f);
	for (int i = 0; i < 50; i++) {
		g[0].push_back(mathop->fx(X[0][0], 0, h, f));

		hessian = this->hessian_1d(X[0][0], 0, h, f);

		this->out += "hessain:" + this->nextline() + "[ ";
		for (int j = 0; j < hessian.size(); j++) {
			for (int l = 0; l < hessian[0].size(); l++) {
				this->out += std::to_string(hessian[j][l]) + " , ";
			}
			this->out += this->nextline();
		}
		this->out += " ]" + this->nextline();

		hessian = this->inv_hessian(hessian);

		this->out += "inv_hessain:" + this->nextline() + "[ ";
		for (int j = 0; j < hessian.size(); j++) {
			for (int l = 0; l < hessian[0].size(); l++) {
				this->out += std::to_string(hessian[j][l]) + " , ";
			}
			this->out += this->nextline();
		}
		this->out += " ]" + this->nextline();

		X_1 = hessian * g * (-1) + X;

		if (fabs(X_1[0][0] + X[0][0] * (-1)) < eps) break;
		else X = X_1;
		this->storePs(X, f);
		g[0].clear();
	}
	this->out += this->nextline();
	this->out += "x: [ " + std::to_string(X[0][0]) + " ]" + this->nextline();
	this->out += "f(x): " + std::to_string(mathop->ComputeFx(X[0][0], f))
		+ this->nextline();

	std::cout << "answer" << "\n";
	std::cout << "x:" << X[0][0] << "\n";
	std::cout << "f(x):" << mathop->ComputeFx(X[0][0], f) << "\n";
}
//1x1
std::vector<std::vector<double>> Newton::hessian_1d(double x, double y, double h, std::vector<std::string> f) {
	std::vector<std::vector<double>> hessian;
	std::vector<double> assist;
	double result = mathop->fxx(x, y, h, f);
	hessian.push_back(assist);
	hessian[0].push_back(result);

	for (int i = 0; i < hessian.size(); i++) {
		for (int j = 0; j < hessian[0].size(); j++) {
			std::cout << hessian[i][j];
		}
		std::cout << "\n";
	}
	return hessian;
}

//2x2
std::vector<std::vector<double>> Newton::hessian_2d(double x, double y, double h, double k, std::vector<std::string> f) {
	std::vector<std::vector<double>> hessian;
	std::vector<double> assist;
	std::queue<double> s_d;

	s_d.push(mathop->fxx(x, y, h, f));
	s_d.push(mathop->fxy(x, y, h, k, f));
	s_d.push(mathop->fxy(y, x, k, h, f));
	s_d.push(mathop->fyy(x, y, k, f));

	for (int i = 0; i < 2; i++) {
		hessian.push_back(assist);
		for (int j = 0; j < 2; j++) {
			hessian[i].push_back(s_d.front());
			s_d.pop();
		}
	}
	std::cout << "\n" << "hessian" << std::endl;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			std::cout << hessian[i][j] << " ";
		}
		std::cout << "\n";
	}
	return hessian;
}

std::vector<std::vector<double>> Newton::inv_hessian(std::vector<std::vector<double>> &m) {
	Matrix temp;
	std::vector<std::vector<double>> result;
	temp.Datas = m;
	result = InverseM(temp);
	std::cout << "inv_hessian" << "\n";

	for (int i = 0; i < m.size(); i++) {
		for (int j = 0; j < m[0].size(); j++) {
			std::cout << result[i][j] << " ";
		}
		std::cout << "\n";
	}
	return result;
}

//Lu decomposition
//don't use it to print out the result directly
//PA=LU
Matrix  Newton::lu_decomposition(std::vector<std::vector<double> > &m) {
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
		rmax = m[r][r];
		row_pos = r;
		for (int i = r + 1; i < dim; i++) {
			if (fabs(m[i][r] - rmax) > 0.00000001) break;
			else if (fabs(m[i][r] - rmax) < 0.00000001) throw "singular matrix";	//singular matrix
			else {
				//if (fabs(m[i][r]) > fabs(rmax)) {	//consider the max/min value always first, then use this.
				rmax = m[i][r];
				row_pos = i;
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

	delete pivot;
	return result;
}
std::vector<double> Newton::solvelineareq(Matrix &m, std::vector<std::vector<double>> &augment) {
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

	delete y;
	delete x;
	return result;
}

std::vector<std::vector<double>> Newton::InverseM(Matrix &m) {
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
	return inverse;
}

void Newton::recycle(std::vector<std::vector<double> > &m) {
	std::vector<std::vector<double>> remove;
	m.swap(remove);
}
void Newton::recycle(std::vector<double> &v) {
	std::vector<double> remove;
	v.swap(remove);
}

std::string Newton::Outcome() {
	return this->out;
}

std::string Newton::nextline() {
	std::string t;
	MarshalString(System::Environment::NewLine, t);
	return t;
}
void Newton::clearOut() {
	this->out.clear();
	//this->xs.clear();
}
void Newton::storePs(std::vector<std::vector<double>> x, std::vector<std::string> f) {
	std::vector<double> t;
	if (x.size()>1) {
		t.push_back(x[0][0]);
		t.push_back(x[0][1]);
		t.push_back(mathop->ComputeFx(x[0][0], x[1][0], f));
	}
	else {
		t.push_back(x[0][0]);
		t.push_back(0);
		t.push_back(mathop->ComputeFx(x[0][0], f));
	}
	this->xs.push_back(t);
	mathop->recycle(t);
}
std::vector<std::vector<double>> Newton::getXs() {
	return this->xs;
}