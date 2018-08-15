#include"Qnewton.h"
Qnewton::Qnewton() {
	mathop = new MathOp();
}
Qnewton::Qnewton(std::vector<std::string> s) {
	mathop = new MathOp();
	this->fx = s;
}
void Qnewton::QnewtonAlgorism(double x, double y, std::vector<std::string> f) {
	double k = 0.0000001, h = 0.0000001, eps = 0.000001;
	std::vector<std::vector<double>> g;
	std::vector<std::vector<double>> F;
	std::vector<std::vector<double>> Xn;
	std::vector<double> t(2, 0);
	F.push_back(t);
	F.push_back(t);
	F[0][0] = 1;
	F[1][1] = 1;
	t.clear();
	g.push_back(t);
	g.push_back(t);

	t.push_back(x);
	Xn.push_back(t);
	t[0] = y;
	Xn.push_back(t);
	this->storeXs(Xn, f);
	this->out += "X:[ ";
	this->out += std::to_string(Xn[0][0]) + " , " +
		std::to_string(Xn[1][0]) + " ]" + this->nextline();

	g[0].push_back(mathop->fx(Xn[0][0], Xn[1][0], h, f));
	g[1].push_back(mathop->fy(Xn[0][0], Xn[1][0], k, f));
	this->out += "g:[ " + std::to_string(g[0][0]) + " , " + std::to_string(g[1][0])
		+ " ]" + this->nextline();
	std::vector<std::vector<double>> d;

	for (int i = 0; i < 1500; i++) {
		if (fabs(g[0][0]) < eps && fabs(g[1][0]) < eps)
			break;
		else
			d = F*g*-1;

		double alpha = this->getAlpha(Xn, d, f);
		//double alpha = this->getAlpha(g,d,F)*-1;
		Xn = Xn + d*alpha;

		this->out += this->nextline() + "X:[ ";
		this->out += std::to_string(Xn[0][0]) + " , " +
			std::to_string(Xn[1][0]) + " ]" + this->nextline();
		this->storeXs(Xn, f);

		std::vector<std::vector<double>> g_next;
		std::vector<double> val(1, 0);
		g_next.push_back(val);
		g_next.push_back(val);

		g_next[0].back() = mathop->fx(Xn[0][0], Xn[1][0], h, f);
		g_next[1].back() = mathop->fy(Xn[0][0], Xn[1][0], k, f);

		std::vector<std::vector<double>> x_step = d*alpha;
		std::vector<std::vector<double>> g_step = g_next + g*-1;

		g[0][0] = g_next[0].back();
		g[1][0] = g_next[1].back();

		this->out += this->nextline() + "g:[ " + std::to_string(g[0][0]) + " , " + std::to_string(g[1][0])
			+ " ]" + this->nextline();

		std::vector<std::vector<double>> temp;
		std::vector<std::vector<double>> ft1 = x_step*this->trans(x_step);
		double ftd1 = mathop->scalar(this->trans(x_step)*g_step);
		temp = ft1*(1 / ftd1);

		std::vector<std::vector<double>>ft2 = (F*g_step)*this->trans(F*g_step);
		double ftd2 = -1 * mathop->scalar(this->trans(g_step)*F*g_step);
		temp = temp + ft2*(1 / ftd2);

		F = F + temp;

		this->out += this->nextline() + "F: [" + this->nextline();
		this->out += std::to_string(F[0][0]) + " , " + std::to_string(F[0][1]) +
			this->nextline() + std::to_string(F[1][0]) + "  , " +
			std::to_string(F[1][1]) + " ]" + this->nextline();

	}

	std::cout << "(x,y): " << Xn[0][0] << " , " << Xn[1][0] << "\n";
	std::cout << "f(x,y )=" << mathop->ComputeFx(Xn[0][0], Xn[1][0], f);

	this->out += this->nextline() + "x,y=[ " + std::to_string(Xn[0][0]) + " , " + std::to_string(Xn[1][0]) + " ] " + this->nextline();
	this->out += "f(x,y)=[ " + std::to_string(mathop->ComputeFx(Xn[0][0], Xn[1][0], f)) +
		" ]" + this->nextline();

}
//use this function to search alpha
double Qnewton::getAlpha(std::vector<std::vector<double>> xn, std::vector<std::vector<double>> d, std::vector<std::string> f) {
	double minalpha = -1, nowval = 0;
	std::vector<std::vector<double>> x = xn;
	std::vector<std::vector<double>> t = xn + d*minalpha;
	double minval = mathop->ComputeFx(t[0][0], t[1][0], f);
	for (double i = -1; i < 1; i += 0.001) {
		t = x + d*i;
		nowval = mathop->ComputeFx(t[0][0], t[1][0], f);
		//std::cout << nowval<<"\n";
		if (nowval<minval) {
			minalpha = i;
			minval = nowval;
		}
		t.clear();
	}return minalpha;
}
//quit it
double Qnewton::getAlpha(std::vector<std::vector<double>> g, std::vector<std::vector<double>> d, std::vector<std::vector<double>> F) {
	double son = mathop->scalar(this->trans(g)*d);
	double mother = mathop->scalar(this->trans(d)*F*d);
	return son / mother;
}

void Qnewton::QnewtonAlgorism(double x, std::vector<std::string> f) {
	double h = 0.000001, g, xn = x;
	g = mathop->fx(xn, 0, h, f);

}
std::vector<std::vector<double>> Qnewton::trans(std::vector<std::vector<double>> v) {
	std::vector<double> r;
	std::vector<std::vector<double>> temp;

	for (int i = 0; i < v[0].size(); i++) {
		for (int j = 0; j < v.size(); j++) {
			r.push_back(v[j][i]);
		}
		temp.push_back(r);
		r.clear();
	}
	return temp;
}

void Qnewton::storeXs(std::vector<std::vector<double>> x, std::vector<std::string>f) {
	std::vector<double> t;
	if (x.size() > 1) {
		t.push_back(x[0][0]);
		t.push_back(x[0][1]);
		t.push_back(mathop->ComputeFx(x[0][0], x[0][1], f));
	}
	else {
		t.push_back(x[0][0]);
		t.push_back(0);
		t.push_back(mathop->ComputeFx(x[0][0], f));
	}
	this->xs.push_back(t);
}

std::vector<std::vector<double>> Qnewton::getPs() {
	return this->xs;
}
std::string Qnewton::getOut() {
	return this->out;
}
std::string Qnewton::nextline() {
	std::string t;
	MarshalString(System::Environment::NewLine, t);
	return t;
}
void Qnewton::clean() {
	this->fx.clear();
	this->out.clear();
}