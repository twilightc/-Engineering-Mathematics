#include"Conjugrate.h"

Conjugrate::Conjugrate() {
	this->mathop = new MathOp();
}
void Conjugrate::conjugrateAgorism(double x, double y, double start, double end, std::vector<std::string> f) {
	double h = 0.0001, k = 0.0001;
	std::vector<std::vector<double>> X;
	std::vector<std::vector<double>> d;
	std::vector<std::vector<double>> Q;
	std::vector<std::vector<double>> alpha;
	std::vector<std::vector<double>> beta;
	std::vector<std::vector<double>> g;
	std::vector<std::vector<double>> g_1;
	std::vector<double> assist;


	//initialze
	for (int i = 0; i < 2; i++) {
		X.push_back(assist);
		d.push_back(assist);
		g.push_back(assist);
		g_1.push_back(assist);
		Q.push_back(assist);
		for (int j = 0; j < 1; j++) g_1[i].push_back(0);
	};
	for (int i = 0; i < 1; i++) {
		alpha.push_back(assist);
		beta.push_back(assist);
	}
	this->out += "X0: [ "+std::to_string(x)+" , "+std::to_string(y)+
				 " ]"+this->nextline();
	X[0].push_back(x);
	X[1].push_back(y);
	this->storePs(X,f);
	g[0].push_back(mathop->fx(X[0][0], X[1][0], h, f));
	g[1].push_back(mathop->fy(X[0][0], X[1][0], k, f));
	this->out += "g0: [ " + std::to_string(g[0][0]) + " , " + std::to_string(g[1][0])
		+ " ]" + this->nextline();
	d = g*(-1);
	//compute hessian
	Q = hessian_2d(X[0][0], X[1][0], h, k, f);
	std::cout << "\n";
	
	//conjugate gradient assuring that the ans will converge in n step(upperbound)
	for (int i = 0; i < X.size();i++){
		if (fabs(g[0][0]) < eps && fabs(g[1][0]) < eps) break;

		//calculate Xn+1
		alpha = (this->trans(g)*d*(-1))*(1 / mathop->scalar(this->trans(d)*Q*d));

		X = X + d*alpha;
		this->storePs(X, f);
		this->out += "X: [ " + std::to_string(X[0][0]) + " , " + std::to_string(X[1][0])
					  + " ]" + this->nextline();
		
		//renew gradient
		g_1[0][0] = mathop->fx(X[0][0], X[1][0], h, f);
		g_1[1][0] = mathop->fy(X[0][0], X[1][0], k, f);


		if (fabs(g[0][0]) < eps && fabs(g[1][0]) < eps) break;
		
		//calculate beta
		beta = this->trans(g_1)*g_1 * (1 / mathop->scalar(this->trans(g)*g));

		//next direction
		d = g_1*(-1) + d*beta;
		g = g_1;
		this->out += "g: [ " + std::to_string(g[0][0]) + " , " + std::to_string(g[1][0])
					  + " ]" + this->nextline();
	}
		this->out += "(x,y): [ " + std::to_string(X[0][0]) + " , " + std::to_string(X[1][0])
						+ " ]" + this->nextline();
		this->out += "f(x,y): [ " + std::to_string(mathop->ComputeFx(X[0][0],X[1][0],f))
					 +" ]"+this->nextline();
	std::cout << "-----\n";
	std::cout << "X[0][0]:" << X[0][0] << ", X[1][0]:" << X[1][0] << std::endl;
	std::cout << "f:" << mathop->ComputeFx(X[0][0], X[1][0],f) << std::endl;

	recycle(d);
	recycle(alpha);
	recycle(beta);
	recycle(g);
	recycle(g_1);
	recycle(Q);
}
std::vector<std::vector<double>> Conjugrate::hessian_2d(double x, double y, double h, double k, std::vector<std::string> f) {
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
	/*
	std::cout << "\n" << "hessian" << std::endl;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			std::cout << hessian[i][j] << " ";
		}
		std::cout << "\n";
	}*/
	return hessian;
}

void Conjugrate::conjugrateAgorism(double x, double start, double end, std::vector<std::string> f) {	
	double interval1 = start, interval2 = end;		//set interval
	double h = 0.000001;
	std::vector<std::vector<double>> X;
	std::vector<std::vector<double>> X1;
	std::vector<std::vector<double>> d;
	std::vector<std::vector<double>> alpha;
	std::vector<std::vector<double>> beta;
	std::vector<std::vector<double>> Q;
	std::vector<double> assist;


	for (int i = 0; i < 1; i++) {
	X.push_back(assist);
	X1.push_back(assist);
	d.push_back(assist);
	alpha.push_back(assist);
	beta.push_back(assist);
	Q.push_back(assist);
	}
	X[0].push_back(x);
	this->out += "X0: [ " + std::to_string(x) +" ]" + this->nextline();

	d[0].push_back((-1)*mathop->fx(X[0][0],0, h, f));
	this->out += "g0: [ " + std::to_string(d[0][0])+ " ]" + this->nextline();


	//conjugate gradient assuring that the ans will converge in n step(upperbound)
	for (int k = 0; k < X.size(); k++) {

	if (fabs(X[0][0]) < eps)break;
	//備胎
	//alpha[0].push_back(mathop->scalar(trans(d)*d * (1 / mathop->scalar(d*d*mathop->fxx(X[0][0], 0, h, f)))));
	//預先使用
	alpha[0].push_back(find_a(0.0, 2.0, X[0][0], d[0][0], f));
	std::cout << alpha[0].back()<<std::endl;

	//calculate Xn+1
	X1[0].push_back(mathop->scalar(d * mathop->scalar(alpha)+ X));
	this->out += "X: [ " + std::to_string(X1[0][0]) + " ]" + this->nextline();
	this->storePs(X1,f);

	//calculate beta
	beta[0].push_back((mathop->fx(X1[0][0],0 ,h, f)*mathop->fx(X1[0][0],0 ,h, f)) /
	(mathop->fx(X[0][0],0, h, f)*mathop->fx(X[0][0],0, h, f)));
	
	//next direction
	d[0][0] = (-1) * mathop->scalar((d * mathop->scalar(beta)))+ (mathop->f_prime(X1[0][0], h, f));
	
	//setting stop condition(3 types)
	this->out += "g0: [ " + std::to_string(d[0][0]) + " ]" + this->nextline();
	
	//assign new X(basiclly useless,since it's 1D)
	X[0][0] = X1[0][0];

	alpha[0].clear();
	beta[0].clear();
	X1[0].clear();		//clear the former recording -- it may cause dimension doesn't match
	}
	this->out += "x: [ " + std::to_string(X[0][0]) +" ]" + this->nextline();
	this->out += "f(x): [ " + std::to_string(mathop->ComputeFx(X[0][0],0, f))+ " ]" + this->nextline();
	std::cout << "X: [" << X[0][0]<<" ]" << std::endl;
	std::cout << mathop->ComputeFx(X[0][0], 0,f) << std::endl;

	//recycle(d);
	//recycle(alpha);
	//recycle(beta);
}
//decenting by spec direction, 1 dimension
double Conjugrate::find_a(double interval1, double interval2, double X, double d, std::vector<std::string> f) {
	double a = 0;
	double n_min = mathop->ComputeFx(X + d*a, 0, f);
	//I don't know why but i should be increased by 0.001, the ans will accurate.
	for (double i = 0.0; i < 2.0; i += 0.001) {
		double current = mathop->ComputeFx(X + d*i,0, f);
		if (current < n_min) {
			a = i;
			n_min = current;
		}
	}
	std::cout << "  alpha: " << a << std::endl;			//print out current alpha
	return a;
}

std::vector<std::vector<double>> Conjugrate::trans(std::vector<std::vector<double>> x) {
	std::vector<double> r;
	std::vector<std::vector<double>> temp;
	for (int i = 0; i < x.back().size(); i++) {
		for (int j = 0; j < x.size(); j++) {
			r.push_back(x[j][i]);
		}
		temp.push_back(r);
	}
	return temp;
}

std::string Conjugrate::Outcome() {
	return this->out;
}
void Conjugrate::recycle(std::vector<std::vector<double> > &m) {
	std::vector<std::vector<double>> remove;
	m.swap(remove);
}
void Conjugrate::recycle(std::vector<double> &v) {
	std::vector<double> remove;
	v.swap(remove);
}
std::vector<std::vector<double>> Conjugrate::getPs() {
	return this->steps;
}
void Conjugrate::storePs(std::vector<std::vector<double>> xn, std::vector<std::string> f) {
	std::vector<double> t;
	if (xn.size() > 1) {
		t.push_back(xn[0][0]);
		t.push_back(xn[0][1]);
		t.push_back(mathop->ComputeFx(xn[0][0], xn[0][1], f));
	}
	else {
		t.push_back(xn[0][0]);
		t.push_back(0);
		t.push_back(mathop->ComputeFx(xn[0][0], 0, f));
	}
	this->steps.push_back(t);
}

std::string Conjugrate::nextline() {
	std::string t;
	MarshalString(System::Environment::NewLine, t);
	return t;
}