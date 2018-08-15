#include "Mvector.h"
#include <math.h>
#include "DotNetUtilities.h"

Mvector::Mvector() {}
//dot運算
double Mvector::D_Product(Vector& v1, Vector& v2) {
	double sum = 0.0,product=0.0;
	std::vector<Vector> cal;
	cal.push_back(v1); cal.push_back(v2);
	//卻認維度是否依樣
	if (!this->check_size(cal)) 
		return 0;
	for (int i = 0; i < v1.len; i++) {
		product=v1.Data[i]*v2.Data[i];
		sum += product;
	}
	return sum;
}
double Mvector::D_Product(std::vector<double>& v1, std::vector<double> &v2) {
	double sum = 0.0, product = 0.0;
	//卻認維度是否依樣
	if (v1.size() != v2.size())
		return NULL;

	for (int j = 0; j < v1.size(); j++) {
		product = 1.0;
		product = v1[j] * v2[j];
		sum += product;
	}
	return sum;
}
//cross product
std::vector<double> Mvector::C_product(Vector &v1,Vector &v2) {
	//確認維度是否一樣
	if (v1.len != v2.len||v1.len != 3) {
		std::vector<double> a; a.push_back(-1);
		return a;
	}
	std::vector<double> temp;
	double i = v1.Data[1] * v2.Data[2] - v1.Data[2] * v2.Data[1];
	double j = -1*(v1.Data[0] * v2.Data[2] - v1.Data[2] * v2.Data[0]);
	double k = v1.Data[0] * v2.Data[1] - v1.Data[1] * v2.Data[0];
	temp.push_back(i); temp.push_back(j); temp.push_back(k);
	return temp;
}
// norm 長度
double Mvector::Norm(std::vector<double> &datas) {
	double sum = 0.0;
	for (int i = 0; i < datas.size(); i++)
		sum += pow(datas[i],2);
	return pow(sum,0.5);
}
// 正規劃normalization
std::vector<double> Mvector::Normal(std::vector<double> &datas) {
	double leng = this->Norm(datas);
	std::vector<double> temp=datas;
	for (int i = 0; i < temp.size(); i++)
		temp[i]=temp[i] / leng;
	return temp;
}

//component a on b
double Mvector::Compo_a_on_b(std::vector<double> &a, std::vector<double> &b) {
	double dot=this->D_Product(a,b);
	double divide = this->Norm(b);
	return dot/divide;
}

//projection a on b
std::vector<double> Mvector::Proj(Vector &a, Vector &b) {
	//確認維度
	if (a.Data.size() != b.Data.size()) {
		std::vector<double> a; a.push_back(-1);
		return a;
	}
	double muti = 0;
	if(this->D_Product(b, b)!=0)
		muti = this->D_Product(a, b) / this->D_Product(b, b);
	std::vector<double> temp;
	for (int i = 0; i < b.len; i++)
		temp.push_back(muti*b.Data[i]);
	return temp;
}
std::vector<double> Mvector::Proj(std::vector<double> &a, std::vector<double> &b) {
	//確認維度
	if (a.size() != b.size()) {
		std::vector<double> a; a.push_back(-1);
		return a;
	}
	double muti = 0;
	if (this->D_Product(b, b) != 0)
		muti = this->D_Product(a, b) / this->D_Product(b, b);
	std::vector<double> temp(a.size(),0);
	for (int i = 0; i < b.size(); i++)
		temp[i]=muti*b[i];
	return temp;
}
//兩向量形成三角形的面積
double Mvector::Tri_area(Vector &v1,Vector &v2) {
	if (v1.Data.size()!=v2.Data.size()&&v1.Data.size()!=3)
		return -1;

	double area = 0.0;
	double vec = this->Norm(this->C_product(v1,v2));
	area = 0.5*vec;
	return area;
}
//平行判斷 0是指維度不同，-1是指不是平行 1表示平行
bool Mvector::Para_Jud(Vector &v1,Vector &v2) {
	//確認維度相同
	std::vector<Vector> check;
	check.push_back(v1); check.push_back(v2);
	if (!this->check_size(check))
		return false;
	//算出相差幾倍
	double comp = v1.Data[0] / v2.Data[0];	
	for (int i = 1; i < v1.len; i++)
		//若某值倍數為否，則回傳-1
		if (comp != v1.Data[i] / v2.Data[i])
			return false;
	return true;
}

//0是維度不同 1是 0否
bool Mvector::Orth_Jud(Vector &v1, Vector &v2) {
	//確認維度
	std::vector<Vector> check;
	check.push_back(v1); check.push_back(v2);
	if (!this->check_size(check))
		return false;
	//內積是否為0
	if (this->D_Product(v1, v2) == 0)
		return true;
	else
		return false;
}
// 兩向量角度
double Mvector::Angle_V(Vector &v1, Vector &v2) {
	double cos = 0.0;
	cos = this->D_Product(v1, v2) /
		 (this->Norm(v1.Data)*this->Norm(v2.Data));
	double angle = std::acos(cos);
	return angle;
}
// 法向量平面
std::vector<double> Mvector::Pl_normal(Vector &v1, Vector &v2) {
	if (v1.len != v2.len||v1.len != 3) {
		std::vector<double> a; a.push_back(-1);
		return a;
	}
	std::vector<double> temp = this->C_product(v2,v1);
	return temp;
}

bool Mvector::Linear_Jud(std::vector<Vector> &v) {
	bool size = false,check=false;
	std::vector<Vector> t=v;
	if (v.size() == 1) return true;
	if (!this->check_size(v)) { 
		check = true;
		size = false;
		return true;
	}
	int start = 0;
	ch(t);
	for (int i = 0; i < t.size()-1; i++) {
		if (t[i].Data[start] == 0.0) {
			start = 1;
			while (t[i].Data[start] == 0.0) {
				start++;
				if (start >= t[i].Data.size())
					break;
			}
		}
		for (int j = i + 1; j < t.size()&&start<t[0].Data.size(); j++) {
			double k = t[j].Data[start] / t[i].Data[start];
			int len = n_0_len(t[j].Data);
			for (int l = start, r = 0; r < len && k != 0; l++, r++) {
				double te = t[j].Data[l] - k * t[i].Data[l];
				if (fabs(te) < pow(10, -8)) t[j].Data[l] = 0.0;
				else t[j].Data[l] = te;
			}
		}
		ch(t);
		start = 0;
	}
	int cnt = 0;
	for (int i = 0; i < t.size(); i++) {
		if (this->n_0_len(t[i].Data)!=0)
			cnt++;
	}
	if(cnt==t.size())
		check= true;
	else check= false;
	return check;
}
void Mvector::ch(std::vector<Vector> &v) {
	for (int i = 0; i < v.size(); i++) {
		//revise at here
		for (int j = i; j < v.size(); j++) {
			if (n_0_len(v[j].Data) > n_0_len(v[i].Data)) {
				std::vector<double> t = v[j].Data;
				v[j] = v[i];
				v[i].Data = t;
			}
		}
	}
}
bool Mvector::Linear_Jud(Vector &v1, Vector &v2) {
	if (v1.len != v2.len) return true;
	int cnt = 1;
	double k = 0;
	if (v2.Data[0] != 0) 
		k = v1.Data[0] / v2.Data[0];
	else k = 0;	
	for (int i = 1; i < v1.Data.size(); i ++)
		if (k == v1.Data[i] / v2.Data[i])
			cnt++;
	cnt++;
	if (cnt != v1.len)
		return true;
	else
		return false;
}
std::vector<std::vector<double>> Mvector::Orth_basis(std::vector<Vector> &u) {
	int len = u[0].len;
	std::vector<std::vector<double>> v;
	if (!this->check_size(u)) return v;

	v.push_back(u[0].Data);
	for (int i = 1; i < u.size(); i++) {
		std::vector<double> proj_t= this->Proj(u[i].Data, v[0]);

		for (int j = 1; j < i; j++) {
			 proj_t = this->plus(proj_t, this->Proj(u[i].Data, v[j]),1);
		}	
		v.push_back(this->plus(u[i].Data,proj_t,-1));
	}
	for (int i = 0; i < u.size(); i++)  v[i] = this->Normal(v[i]);

	std::cout << " below:  \n";
	for (int i = 0; i < v.size(); i++) {
		for (int j = 0; j < len; j++) {
			std::cout <<v[i][j]<< " ";
		}
		std::cout << "\n";
	}
	return v;
}
std::vector<double> Mvector::plus(std::vector<double>&a,std::vector<double>&b,int mut) {
	std::vector<double> t;
	for (int i = 0; i < a.size(); i++) t.push_back(a[i]+(b[i]*mut));
	return t;
}
//確認是否為純量
bool Mvector::isDigit(std::string s) {
	for (int i = 0; i < s.size(); i++)
		if (s.at(i)<'0' || s.at(i)>'9') 
			if(s.at(i)!='.')
				return false;
	return true;
}
//確認維度
int Mvector::n_0_len(std::vector<double> &v) {
	for (int i = 0; i < v.size(); i++)
		if (fabs(v[i]) > 0.0000000001)
			return v.size() - i;
	return 0;
}
//sizecheck
bool Mvector::check_size(std::vector<Vector> &datas) {
	double len = datas[0].Data.size();
	if (datas.size() > 1) {
		for (int i = 1; i < datas.size(); i++) {
			if (len != datas[i].Data.size()) {
				return false;
			}
		}
	}
	return true;
}
//尋找運算符號
bool Mvector::to_find(System::String ^in) {
	std::string menu[19] = { "(", ",", "+","-",",", "*",
							"dot", "Norm" ,"Normal", "Crossp", "Compo",
							"Proj", "Tri", "Paral","Oj", "Angle",
							"P_n", "Lij", "Ob"};
	std::string temp = "";
	MarshalString(in, temp);
	for (int i = 0; i < 19; i++)
		if (temp == menu[i].c_str()) {
			return true;
		}
	return false;
}
//search function
bool Mvector::search_func(std::string func) {
	std::string menu[13] = { "dot","Norm" ,"Normal", "Crossp", "Compo", 
							 "Proj", "Tri", "Paral","Oj", "Angle", 
							 "P_n", "Lij", "Ob" };
	for (int i = 0; i < 13; i++)
		if (func == menu[i])
			return true;
	return false;
}
//find function
void Mvector::FindFunc(std::vector<Vector> &v,std::string func,Ans *ans) {
	/*for (int i = 0; i < v.size(); i++) {
		for (int j = 0; j < v[0].Data.size(); j++)
			std::cout << v[i].Data[j] << " ";
		std::cout << "\n";
	}*/
	std::cout << "\n"<<func<<"\n";
	if (func == "Norm") {
		ans->dAns=this->Norm(v[0].Data);
		ans->jud = 0;
	}
	else if (func == "dot") {
		ans->dAns=this->D_Product(v[0],v[1]);
		ans->jud = 0;
		std::cout << "dot: " << ans->dAns << "\n";
	}
	else if (func == "Normal") {
		ans->AnsData=this->Normal(v[0].Data);
		ans->jud= 1;
	}
	else if (func == "Crossp"&&v.size()==2) {
		if (v[0].len == 3 && v[1].len == 3)
			ans->AnsData=this->C_product(v[0],v[1]);
		ans->jud= 1;
	}
	else if (func=="Compo"&&v.size()==2) {
		ans->dAns=this->Compo_a_on_b(v[1].Data,v[0].Data);
		ans->jud= 0;
	}
	else if (func == "Proj") {
		ans->AnsData=this->Proj(v[1], v[0]);
		ans->jud= 1;
	}
	else if (func == "Tri") {
		ans->dAns=this->Tri_area(v[1],v[0]);
		ans->jud= 0;
	}
	else if (func == "Paral") {
		ans->bAns=this->Para_Jud(v[0],v[1]);
		ans->jud= -1;
	}
	else if (func == "Oj") {
		ans->bAns=this->Orth_Jud(v[0],v[1]);
		ans->jud= -1;
	}
	else if (func=="Angle") {
		ans->dAns=this->Angle_V(v[0],v[1]);
		ans->jud= 0;
	} 
	else if (func=="P_n") {
		ans->AnsData=this->Pl_normal(v[0],v[1]);
		ans->jud= 1;
	}
	else if (func=="Lij") {
		if(v.size()==2) ans->bAns=this->Linear_Jud(v[0],v[1]);
		else if (v.size() == 1) ans->bAns = true;
		else ans->bAns = this->Linear_Jud(v);
		ans->jud= -1;
	}
	else if (func == "Ob") {
		std::vector<Vector> v2;
		while(!v.empty()){
			v2.push_back(v.back());
			v.pop_back();
		}
		ans->vvAns= this->Orth_basis(v2);
		ans->jud= 2;
		if (ans->vvAns.empty()) ans->jud = -1000;
	}
}
//finding certain operator
bool Mvector::find_operator(std::string op) {
	std::string menu[4] = { "+", "-", "*", "/", };
	for (int i = 0; i < 4; i++)
		if (op == menu[i].c_str())
			return true;
	return false;
}

//finding certain vector
void Mvector::find_VectorName(std::vector<Vector> &vectors, std::deque<Vector> &Vetcor_stack, std::string &str) {
	for (int i = 0; i < vectors.size(); i++) {
		if (str == vectors[i].Name) {
			Vetcor_stack.push_back(vectors[i]);
		}
	}
}

//sizecheck for vector
bool Mvector::check_size_v(std::deque<Vector> &datas) {
	double len = datas[0].len;
	if (datas.size() > 1) {
		for (int i = 1; i < datas.size(); i++){
			if (datas[i].Name == "operator") continue;
			if (len != datas[i].len)
				return false;
		}
	}
	return true;
}

//calculate result of postfix
std::string Mvector::mathematics_operation(std::vector<Vector> &vecs, std::vector<std::string> &formula,std::vector<int> &index) {
	Mvector mv;
	Ans ans;
	std::vector<std::string> op_stack;		//save operation
	std::deque<Vector> Vector_stack;	//save Vector
	std::vector<Vector> v;
	Vector const_for_oprate;
	Vector result;
	double d;		//save constant
	bool sizecheck;
	bool done = false;
	result.Name = "result";		//convenient to check
	const_for_oprate.Name = "operator";

	//bool Mvector::check_size(std::vector<Vector> &datas)

	for (int i = 0; i < formula.size(); i++) {
		//std::cout << formula.at(i) << std::endl;
		if (index[i]==-1) {
			//push into operator_stack
			if (formula.at(i) == "+") {
				if (!done) {	//檢查維度只需一次
					sizecheck = check_size_v(Vector_stack);
					if (!sizecheck) {
						return "size doesn't match";
					}
					done = true;
				}
				result = Vector_stack.at(Vector_stack.size() - 2) + Vector_stack.at(Vector_stack.size() - 1);
				std::cout << result.Data.at(0) << std::endl;
				//運算完成，處理剩餘資源
				for (int i = 0; i < 2; i++) Vector_stack.pop_back();
				Vector_stack.push_back(result);
			}
			else if (formula.at(i) == "-") {
				if (!done) {	//檢查維度只需一次
					sizecheck = check_size_v(Vector_stack);
					if (!sizecheck) {
						return "size doesn't match";
					}
					done = true;
				}
				result = Vector_stack.at(Vector_stack.size() - 2) + (-1) * Vector_stack.at(Vector_stack.size() - 1);

				std::cout << result.Data.at(0) << std::endl;
				for (int i = 0; i < 2; i++) Vector_stack.pop_back();
				Vector_stack.push_back(result);
			}
			else if (formula.at(i) == "*") {
				if ((Vector_stack.at(Vector_stack.size() - 2).Name) == "operator"||
					 (Vector_stack.at(Vector_stack.size()-2).Data.size()==1))
					result = Vector_stack.at(Vector_stack.size() - 2).Data * Vector_stack.at(Vector_stack.size() - 1);
				
				std::cout << result.Data.at(0) << std::endl;

				for (int i = 0; i < 2; i++) Vector_stack.pop_back();
					Vector_stack.push_back(result);
			}
			else if (formula[i] == ",") {
				Vector_stack.back().len = Vector_stack.back().Data.size();
				v.push_back(Vector_stack.back());
				Vector_stack.pop_back();
			}
			else if (this->search_func(formula[i])) {
				v.push_back(Vector_stack[0]);
				this->FindFunc(v,formula[i],&ans);
				Vector_stack.clear();
			}
		}
		else if (index[i]==0) {
			const_for_oprate.Data.push_back(std::stod(formula.at(i), NULL));
			Vector_stack.push_back(const_for_oprate);
		}
		else {
			find_VectorName(vecs, Vector_stack, formula.at(i));
		}
		const_for_oprate.Data.clear();
		result.Data.clear();
	}
	//std::cout << Vector_stack.at(0).Data.size() << std::endl;
	if (!Vector_stack.empty()) {
		std::cout << Vector_stack[0].Data.size() << std::endl;
		ans.Vans = Vector_stack[0];
		ans.jud = 3;
	}
	if (ans.jud == 0) {
		std::cout<<ans.dAns<<"\n";
	}
	else if (ans.jud == 1) {
		for (int i = 0; i < ans.AnsData.size(); i++)
			std::cout<<ans.AnsData[i]<<" ";
	}
	else if (ans.jud == -1) {
		if (ans.bAns)
			ans.strAns = "yes";
		else ans.strAns = "no";
		std::cout << ans.strAns << "\n";
	}
	else if (ans.jud == 3) {
		for (int i = 0; i < ans.Vans.Data.size(); i++)
			std::cout << ans.Vans.Data[i] << "  ";
	}

	return this->Mk_str(ans);
}
std::string Mvector::Mk_str(Ans &ans) {
	std::string result="[ ";
	if (ans.jud == 0) {
		result += std::to_string(ans.dAns);
	}
	else if (ans.jud == 1) {
		for (int i = 0; i < ans.AnsData.size(); i++)
			result += std::to_string(ans.AnsData[i])+", ";
	}
	else if (ans.jud == -1) {
		if (ans.bAns)
			ans.strAns = "yes";
		else ans.strAns = "no";
		result += ans.strAns;
	}
	else if (ans.jud == 3) {
		for (int i = 0; i < ans.Vans.Data.size(); i++)
			result += std::to_string(ans.Vans.Data[i])+", ";
	}
	else if (ans.jud == 2) {
		std::string n = "";
		MarshalString(System::Environment::NewLine, n);

		for (int i = 0; i < ans.vvAns.size(); i++) {
			for (int j = 0; j < ans.vvAns[0].size(); j++) {
				result += std::to_string(ans.vvAns[i][j])+", ";
			}
			result += n;
		}
	}
	else result+= "error!!";
	result += " ]";
	return result;
}