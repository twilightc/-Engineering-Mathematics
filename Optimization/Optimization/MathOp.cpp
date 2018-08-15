#include "MathOp.h"
#include"DotNetUtilities.h"
#include<sstream>
#include <stack>
#include<deque>

MathOp::MathOp() {
}
//input the computed f(x)
void MathOp::storeFx(std::vector<std::string> f) {
	this->fs.push_back(f);
}
//f(x) to format
std::vector<std::string> MathOp::split_to_format(std::string input) {
	std::vector<std::string> temp;
	std::stringstream ss;
	std::string s;
	int cnt = 0;
	bool check = false;
	for (int i = 0; i < input.size(); i++) {
		//let a char to a string
		std::string t;
		ss << input.at(i);
		ss >> t;
		ss.clear();
		if (this->isSymbol(t)) {
			// if t is an operation , let s into temp
			if (!s.empty())
				temp.push_back(s);
			temp.push_back(t);
			s.clear();
		}
		// if this char is not a dot and an operation , judge that it is a number
		else  s += t;
		t.clear();
	}
	if (!s.empty()) temp.push_back(s);
	//add backets to fx
	std::vector<std::string> t = this->pow_backet(temp);
	t = this->div_backet(t);
	t = this->mut_backet(t);
	t = this->sub_backet(t);
	return t;
}

//add backets to f(x) (pow) 
std::vector<std::string> MathOp::pow_backet(std::vector<std::string> t) {
	std::deque<std::string> out;
	for (int i = 0; i < t.size(); i++) {
		if (t[i] == "^") {
			// first, add "("
			// if "?^" this condition
			if (out[out.size() - 1] != ")") {
				int j = out.size() - 1;
				while (out[j] != "+" && out[j] != "-" && out[j] != "*" &&out[j]!="/"&& j - 1 >= 0)  j--;
				if (j - 1 < 0) out.push_front("(");
				else out.insert(out.begin() + j + 1, "(");
			}
			else {
				// to count the numbers of the "()"
				int cnt = 0;
				int j = out.size() - 1;
				while (j >= 0) {
					if (out[j] == ")") cnt++;
					else if (out[j] == "(") cnt--;

					if (cnt == 0 && j != t.size() - 1) break;
					j--;
				}
				if (j - 1 < 0) out.push_front("(");
				else out.insert(out.begin() + j, "(");
			}
			out.push_back(t[i]);
			//and then add )
			if (t[i + 1] == "(") {
				//if ^( this condition
				int j = i + 1;
				while (t[j] != ")"&&j < t.size()) {
					out.push_back(t[j]);
					j++;
				}
				out.push_back(t[j]);
				out.push_back(")");
				i = j;
			}
			else {
				//if ^? this condition. ? is a number
				int  j = i + 1;
				//meet the + - * , so add )
				while ((t[j] != "+" && t[j] != "-"&&t[j] != "*"&&t[j]!="/")) {
					out.push_back(t[j]);
					j++;
					if (j >= t.size()) break;
				}
				out.push_back(")");
				i = j - 1;
			}
		}
		else out.push_back(t[i]);
	}

	std::vector<std::string> ans;
	for (int i = 0; i < out.size(); i++) {
		ans.push_back(out[i]);
	}
	out.clear();
	return ans;
}

//add parentheses while encounter substract 
std::vector<std::string> MathOp::sub_backet(std::vector<std::string> f) {
	std::deque<std::string> t;
	for (int i = 0; i < f.size(); i++) {
		if (f[i] == "-") {
			// if meet * , add "(" first.
			if (t[t.size() - 1] != ")") {
				int j = t.size() - 1;
				while (t[j] != "+"&&t[j] != "-"&&j - 1 >= 0)  j--;

				if (j - 1 < 0) t.push_front("(");
				else t.insert(t.begin() + j + 1, "(");
			}
			else {
				// to count the numbers of the "()"
				int cnt = 0;
				int j = t.size() - 1;
				while (j >= 0) {
					if (t[j] == ")") cnt++;
					else if (t[j] == "(") cnt--;

					if (cnt == 0 && j != t.size() - 1) break;
					j--;
				}
				if (j - 1 < 0) t.push_front("(");
				else t.insert(t.begin() + j, "(");
			}
			t.push_back(f[i]);
			if (f[i + 1] != "(") {
				int j = i + 1;
				while (f[j] != "+"&&f[j] != "-"&&f[j] != "*"&&f[j] != ")"&&f[j]!="/"&&f[j] != "(") {
					t.push_back(f[j]);
					j++;
					if (j >= f.size()) break;
				}
				t.push_back(")");
				i = j - 1;
			}
			else {
				int j = i + 1;
				// to count the numbers of the "()"
				int cnt = 0;
				while (j < f.size()) {
					if (f[j] == "(") cnt++;
					else if (f[j] == ")") cnt--;

					if (cnt == 0 && j != i + 1) break;
					t.push_back(f[j]);
					j++;
				}
				t.push_back(f[j]);
				t.push_back(")");
				i = j;
			}
		}
		else t.push_back(f[i]);
	}
	std::vector<std::string > out;
	for (int i = 0; i < t.size(); i++) {
		out.push_back(t[i]);
	}
	return out;

}


//add backets to f(x) (*)
std::vector<std::string> MathOp::mut_backet(std::vector<std::string> f) {
	std::deque<std::string> t;
	for (int i = 0; i < f.size(); i++) {
		if (f[i] == "*") {
			// if meet * , add "(" first.
			t.push_back(f[i]);
			if (f[i - 1] != ")") {
				int j = t.size() - 1;
				while (t[j] != "/"&&t[j] != "+"&&t[j] != "-"&&j - 1 >= 0)  j--;

				if (j - 1 < 0) t.push_front("(");
				else t.insert(t.begin() + j + 1, "(");
			}
			else {
				// to count the numbers of the "()"
				int cnt = 0;
				int j = t.size() - 1;
				while (j >= 0) {
					if (t[j] == ")") cnt++;
					else if (t[j] == "(") cnt--;

					if (cnt == 0 && j != t.size() - 1) break;
					j--;
				}
				if (j - 1 < 0) t.push_front("(");
				else t.insert(t.begin() + j, "(");
			}
			if (f[i + 1] != "(") {
				int j = i + 1;
				while (f[j] != "+"&&f[j] != "-"&&f[j] != "*"&&f[j] != ")"&&f[j] != "("&&f[j]!="/") {
					t.push_back(f[j]);
					j++;
					if (j >= f.size()) break;
				}
				t.push_back(")");
				i = j - 1;
			}
			else {
				int j = i + 1;
				// to count the numbers of the "()"
				int cnt = 0;
				while (j < f.size()) {
					if (f[j] == "(") cnt++;
					else if (f[j] == ")") cnt--;

					if (cnt == 0 && j != i + 1) break;
					t.push_back(f[j]);
					j++;
				}
				t.push_back(f[j]);
				t.push_back(")");
				i = j;
			}
		}
		else t.push_back(f[i]);
	}
	std::vector<std::string > out;
	for (int i = 0; i < t.size(); i++) {
		out.push_back(t[i]);
	}
	return out;
}
std::vector<std::string> MathOp::div_backet(std::vector<std::string> f) {
	std::deque<std::string> t;
	for (int i = 0; i < f.size(); i++) {
		if (f[i] == "/") {
			// if meet / , add "(" first.
			if (f[i - 1] != ")") {
				int j = t.size() - 1;
				while (t[j] != "/"&&t[j] != "+"&&t[j] != "-"&&t[j] != "*" && j - 1 >= 0)  j--;

				if (j - 1 < 0) t.push_front("(");
				else t.insert(t.begin() + j + 1, "(");
			}
			else {
				// to count the numbers of the "()"
				int cnt = 0;
				int j = t.size() - 1;
				while (j >= 0) {
					if (t[j] == ")") cnt++;
					else if (t[j] == "(") cnt--;

					if (cnt == 0 && j != t.size() - 1) break;
					j--;
				}
				if (j - 1 < 0) t.push_front("(");
				else t.insert(t.begin() + j, "(");
			}
			t.push_back(f[i]);
			if (f[i + 1] != "(") {
				int j = i + 1;
				while (f[j] != "+"&&f[j] != "-"&&f[j] != "*"&&f[j] != "/"&&f[j] != ")"&&f[j] != "(") {
					t.push_back(f[j]);
					j++;
					if (j >= f.size()) break;
				}
				t.push_back(")");
				i = j - 1;
			}
			else {
				int j = i + 1;
				// to count the numbers of the "()"
				int cnt = 0;
				while (j < f.size()) {
					if (f[j] == "(") cnt++;
					else if (f[j] == ")") cnt--;

					if (cnt == 0 && j != i + 1) break;
					t.push_back(f[j]);
					j++;
				}
				t.push_back(f[j]);
				t.push_back(")");
				i = j;
			}
		}
		else t.push_back(f[i]);
	}
	std::vector<std::string > out;
	for (int i = 0; i < t.size(); i++) {
		out.push_back(t[i]);
	}
	return out;
}
//to postfix
std::vector<std::string> MathOp::format(std::vector<std::string> f) {
	std::vector<std::string> tempOut;
	std::stack<std::string> Opstack;
	for (int i = 0; i < f.size(); i++) {
		if (f[i] == "x" || f[i] == "y") {
			tempOut.push_back(f[i]);
		}
		else if (this->isSymbol(f[i]) && f[i] != ")") {
			Opstack.push(f[i]);
		}
		else if (f[i] == ")")
			for (int j = 0; !Opstack.empty(); j++)
				if (Opstack.top() == "(") {
					Opstack.pop();
					break;
				}
				else {
					tempOut.push_back(Opstack.top());
					Opstack.pop();
				}
		else {
			tempOut.push_back(f[i]);
		}
	}

	if (!Opstack.empty())
		while (!Opstack.empty()) {
			tempOut.push_back(Opstack.top());
			Opstack.pop();
		}
	return tempOut;
}

//(to two variables fx)
//use postfix to compute the output of the fx inputted 
double MathOp::ComputeFx(double a, double b, std::vector<std::string> fx) {
	std::vector<std::string> f = this->format(fx);
	std::vector<double> t;
	double ans = 0.0;

	for (int i = 0; i < f.size(); i++) {
		if (f[i] == "x") {
			t.push_back(a);
		}
		else if (f[i] == "y") {
			t.push_back(b);
		}
		else if (isSymbol(f[i])) {
			if (f[i] == "+") {
				double temp = t[t.size() - 2] + t[t.size() - 1];
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
			else if (f[i] == "-") {
				double temp = t[t.size() - 2] + (-1) * t[t.size() - 1];
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
			else if (f[i] == "*") {
				double temp = t[t.size() - 2] * t[t.size() - 1];
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
			else if (f[i] == "^") {
				double temp = pow(t[t.size() - 2], t[t.size() - 1]);
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
			else if (f[i] == "/") {
				double temp = t[t.size() - 2] / t[t.size() - 1];
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
		}
		else {
			t.push_back(std::stod(f[i]));
		}

	}
	ans = t.back();
	return ans;
}
// to only one variable f(x)
double MathOp::ComputeFx(double a, std::vector<std::string> fx) {
	std::vector<std::string> f = this->format(fx);
	std::vector<double> t;
	double ans = 0.0;

	for (int i = 0; i < f.size(); i++) {

		if (f[i] == "x") {
			t.push_back(a);
		}
		else if (isSymbol(f[i])) {
			if (f[i] == "+") {
				double temp = t[t.size() - 2] + t[t.size() - 1];
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
			else if (f[i] == "-") {
				double temp = t[t.size() - 2] + (-1) * t[t.size() - 1];
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
			else if (f[i] == "*") {
				double temp = t[t.size() - 2] * t[t.size() - 1];
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
			else if (f[i] == "^") {
				double temp = pow(t[t.size() - 2], t[t.size() - 1]);
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
			else if (f[i] == "/") {
				double temp = t[t.size() - 2] / t[t.size() - 1];
				t.pop_back();
				t.pop_back();
				t.push_back(temp);
			}
		}
		else {
			t.push_back(std::stod(f[i]));
		}
		ans = t.back();
	}
	return ans;
}
// to choose which fx we use into the optimization algorism
bool MathOp::ch_get_fs(int ch, std::vector<std::string> &get) {
	if (ch < this->fs.size()) {
		get = this->fs.at(ch);
		return true;
	}
	return false;
}
// assure that string s is a number
bool MathOp::isDigit(std::string s) {
	for (int i = 0; i < s.size(); i++)
		if (s.at(i) < '0' || s.at(i) > '9') {
			if (s.at(i) != '.')
				return false;
		}

	return true;
}
// assure that string s is an operation
bool MathOp::isSymbol(std::string s) {
	std::string op[7] = { "+","-","^","(",")","*","/" };
	for (int i = 0; i < 7; i++) {
		if (s == op[i])
			return true;
	}
	return false;
}
// to test if the fs is true
void MathOp::print() {
	for (int i = 0; i < this->fs.size(); i++) {
		for (int j = 0; j < fs[i].size(); j++) {
			std::cout << fs[i][j] << " ";
		}
		std::cout << "\n";
	}
}

//____________________________deriave______________________________________
double MathOp::f_prime(double x, double h, std::vector<std::string> f) {
	return 0.5*(fx(x + h, 0, h, f) -fx(x - h, 0, h, f)) / h;
}

double MathOp::s_prime(double x, double h, std::vector<std::string> f) {
	return 0.5*(this->fx(x + h, 0, h, f) - this->fx(x - h, 0, h, f));
}

double MathOp::fx(double x, double y, double h, std::vector<std::string> f) {
	std::cout<<"" << this->ComputeFx(x + h, y, f) << "   " << this->ComputeFx(x - h, y, f) << std::endl;
	return 0.5*(this->ComputeFx(x + h, y, f) - this->ComputeFx(x - h, y, f)) / h;
}

double MathOp::fy(double x, double y, double k, std::vector<std::string> f) {
	return 0.5*(this->ComputeFx(x, y + k, f) - this->ComputeFx(x, y - k, f)) / k;
}

double MathOp::fxx(double x, double y, double h, std::vector<std::string> f) {
	double f1 = this->ComputeFx(x + h, y, f);
	double fm = 2 * this->ComputeFx(x, y, f);
	double f2 = this->ComputeFx(x - h, y, f);
	return (f1 - fm + f2) / pow(h, 2);
}

double MathOp::fyy(double x, double y, double k, std::vector<std::string> f) {
	double f1 = this->ComputeFx(x, y + k, f);
	double fm = 2 * this->ComputeFx(x, y, f);
	double f2 = this->ComputeFx(x, y - k, f);
	return (f1 - fm + f2) / pow(k, 2);
}
double MathOp::fxy(double x, double y, double h, double k, std::vector<std::string> f) {
	double f1 = this->ComputeFx(x + h, y + k, f);
	double f2 = this->ComputeFx(x + h, y - k, f);
	double f3 = this->ComputeFx(x - h, y + k, f);
	double f4 = this->ComputeFx(x - h, y - k, f);
	return (f1 - f2 - f3 + f4) / (4 * h*k);
}
double MathOp::scalar(std::vector<std::vector<double>> &m) {
	return ((m.size() == 1) && (m[0].size() == 1)) ? m[0][0] : throw "cannot convert to scalar";
}

void MathOp::recycle(std::vector<std::vector<double> > &m) {
	std::vector<std::vector<double>> remove;
	m.swap(remove);
}

void MathOp::recycle(std::vector<double> &v) {
	std::vector<double> remove;
	v.swap(remove);
}
void MathOp::storePs(std::vector<std::vector<double>> x, std::vector<std::string> f) {
	std::vector<double> t;
	t.push_back(x[0][0]);
	t.push_back(x[1][0]);
	t.push_back(this->ComputeFx(x[0][0], x[0][1], f));
	this->steps.push_back(t);
	this->recycle(t);
}

std::vector<std::vector<double>> MathOp::getPs() {
	return this->steps;
}