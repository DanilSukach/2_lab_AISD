#include <iostream>
#include <locale>
#include <cmath>
#include <iomanip>
#include <complex>
using namespace std;

template<typename T>
class Polynomial {
	struct list {
		size_t kof;
		T data;
		list* link;
		list(const T& data, const size_t kof, list* link) {
			this->data = data;
			this->kof = kof;
			this->link = link;
		}
	}*_head;
	double eps;
public:
	Polynomial() {
		_head = nullptr;
		eps = 0.0001;
	}
	Polynomial(const size_t n) {
		eps = 0.0001;
		list* p = new list(1, 0, nullptr);
		_head = p;
		list* tmp = _head;
		for (size_t i = 1; i <= n; i++) {
			list* p1 = new list(1, i, nullptr);
			tmp->link = p1;
			tmp = tmp->link;
		}
	}
	bool operator ==(const Polynomial& v) {
		Polynomial temp;
		list* p1 = _head;
		list* p2 = v._head;
		while (p1 || p2) {
			if (p1 && p2) {
				if (abs(p1->data - p2->data) > eps) {
					return false;
				}
				else {
					p1 = p1->link;
					p2 = p2->link;
				}
			}
			else {
				if (p1) {
					if (p1->data == 0) {
						p1 = p1->link;
					}
					else {
						return false;
					}
				}
				if (p2) {
					if (p2->data == 0) {
						p2 = p2->link;
					}
					else {
						return false;
					}
				}
			}
		}
		return true;
	}
	void Set(const T& value, const size_t ind) {
		list* p = _head;
		while (p) {
			if (p->kof == ind) {
				p->data = value;
				break;
			}
			p = p->link;
			if (!p) throw "invalid index";
		}
	}
	T& operator[](const size_t ind) const {
		if (ind < 0) throw "invalid index";
		list* p = _head;
		for (size_t i = 0; i < ind; i++) {
			p = p->link;
			if (!p) throw "invalid index";
		}
		return p->data;
	}
	friend ostream& operator<<(ostream& os, const Polynomial& v) {
		list* p = v._head;
		while (p) {
			if (p->data == 0) {
				p = p->link;
			}
			else {
				if (p->data > 0) {
					if (p->kof == 0) {
						os << p->data;
					}
					else {
						if (p->data == 1) {
							os << "+x^" << p->kof;
						}
						else {
							os << "+" << p->data << "x^" << p->kof;
						}
					}
				}
				else {
					if (p->kof == 0) {
						os << p->data;
					}
					else {
						if (p->data == -1) {
							os << "-x^" << p->kof;
						}
						else {
							os << p->data << "x^" << p->kof;
						}
					}
				}
				p = p->link;
			}
		}
		return os;
	}

	Polynomial operator +(const Polynomial& v) const {
		Polynomial temp;
		list* p1 = _head;
		list* p2 = v._head;
		T temp_value;
		while (p1 || p2) {
			if (p1 && p2) {
				temp_value = p1->data + p2->data;
				push_back(temp, temp_value, p1->kof);
				p1 = p1->link;
				p2 = p2->link;
			}
			else {
				if (p1) {
					temp_value = p1->data;
					push_back(temp, temp_value, p1->kof);
					p1 = p1->link;
				}
				if (p2) {
					temp_value = p2->data;
					push_back(temp, temp_value, p2->kof);
					p2 = p2->link;
				}
			}
		}
		return temp;
	}
	Polynomial operator -(const Polynomial& v) const {
		Polynomial temp;
		list* p1 = _head;
		list* p2 = v._head;
		T temp_value;
		while (p1 || p2) {
			if (p1 && p2) {
				temp_value = p1->data - p2->data;
				push_back(temp, temp_value, p1->kof);
				p1 = p1->link;
				p2 = p2->link;
			}
			else {
				if (p1) {
					temp_value = p1->data;
					push_back(temp, temp_value, p1->kof);
					p1 = p1->link;
				}
				if (p2) {
					temp_value = -p2->data;
					push_back(temp, temp_value, p2->kof);
					p2 = p2->link;
				}
			}
		}
		return temp;
	}
	Polynomial& operator *=(const T& value) {
		list* p = _head;
		while (p) {
			p->data *= value;
			p = p->link;
		}
		return *this;
	}
	Polynomial operator *(const T& value) const {
		Polynomial temp;
		list* p = _head;
		T value_temp;
		while (p) {
			value_temp = p->data * value;
			push_back(temp, value_temp, p->kof);
			p = p->link;
		}
		return temp;
	}
	friend Polynomial operator *(const T& value, const Polynomial& v) {
		Polynomial temp;
		list* p = v._head;
		T value_temp;
		while (p) {
			value_temp = p->data * value;
			temp.push_back(temp, value_temp, p->kof);
			p = p->link;
		}
		return temp;
	}
	double calculation(const T& x) const {
		double result = 0;
		list* p = _head;
		while (p) {
			result += p->data * pow(x, p->kof);
			p = p->link;
		}
		return result;
	}
	void push_back(Polynomial& temp, const T& value, const size_t kof)const {
		if (temp._head) {
			list* p1 = temp._head;
			while (p1->link) {
				p1 = p1->link;
			}
			list* p2 = new list(value, kof, nullptr);
			p1->link = p2;
		}
		else
		{
			list* p1 = new list(value, kof, nullptr);
			temp._head = p1;
		}
	}
	void formula_Kardano(double& x1, double& x2, double& x3, double& Q) const {
		list* tmp = _head;
		T d = tmp->data;
		tmp = tmp->link;
		if (!tmp)throw"invalid Polinomial";
		T c = tmp->data;
		tmp = tmp->link;
		if (!tmp)throw"invalid Polinomial";
		T b = tmp->data;
		tmp = tmp->link;
		if (!tmp)throw"invalid Polinomial";
		T a = tmp->data;
		tmp = tmp->link;
		if (tmp)throw"invalid Polinomial";
		double p = (3 * static_cast<double>(a) * static_cast<double>(c) - pow(b, 2)) / (3 * pow(a, 2));
		double q = (2 * pow(b, 3) - 9 * static_cast<double>(a) * static_cast<double>(b) * static_cast<double>(c) + 27 * pow(a, 2) * d) / (27 * pow(a, 3));
		Q = pow((p / 3), 3) + pow((q / 2), 2);
		cout << "Q = " << Q << endl;
		double u = cbrt((-q / 2) + sqrt((pow(q, 2) / 4) + (pow(p, 3) / 27)));
		double v = cbrt((-q / 2) - sqrt((pow(q, 2) / 4) + (pow(p, 3) / 27)));
		if (Q > 0) {
			double y = u + v;
			x1 = y - (b / (3 * a));
		}
		if (abs(Q) == 0) {
			double y1 = 2 * cbrt(-q / 2);
			double y2 = -cbrt(-q / 2);
			x1 = y1 - (b / (3 * a));
			x2 = x3 = y2 - (b / (3 * a));
		}
		if (Q < 0) {
			double fi = 0;
			if (abs(q) < 0.1) {
				fi = 3.1415926535 / 2;
			}
			if (q > 0) {
				fi = atan(sqrt(-Q) / (-q / 2)) + 3.1415926535;
			}
			if (q < 0) {
				fi = atan(sqrt(-Q) / (-q / 2));
			}
			x1 = (2 * sqrt(-p / 3) * cos(fi / 3)) - (b / (3 * a));
			x2 = (2 * sqrt(-p / 3) * cos((fi + 2 * 3.1415926535) / 3)) - (b / (3 * a));
			x3 = (2 * sqrt(-p / 3) * cos((fi + 4 * 3.1415926535) / 3)) - (b / (3 * a));
		}
	}
};
template <typename T>

class Polynomial<complex<T>> {
	struct list {
		size_t kof;
		complex<T> data;
		list* link;
		list(const complex<T>& data, const size_t kof, list* link) {
			this->data = data;
			this->kof = kof;
			this->link = link;
		}
	}*_head;
	double eps;
public:
	Polynomial() {
		_head = nullptr;
		eps = 0.0001;
	}
	Polynomial(const size_t n) {
		eps = 0.0001;
		complex<float> x(1.f, 1.f);
		list* p = new list(x, 0, nullptr);
		_head = p;
		list* tmp = _head;
		for (size_t i = 1; i <= n; i++) {
			list* p1 = new list(x, i, nullptr);
			tmp->link = p1;
			tmp = tmp->link;
		}
	}
	complex<T>& operator[](const size_t ind) const {
		if (ind < 0) throw "invalid index";
		list* p = _head;
		for (size_t i = 0; i < ind; i++) {
			p = p->link;
			if (!p) throw "invalid index";
		}
		return p->data;
	}
	void Set(const complex<T>& value, const size_t ind) {
		list* p = _head;
		while (p) {
			if (p->kof == ind) {
				p->data = value;
				break;
			}
			p = p->link;
			if (!p) throw "invalid index";
		}
	}
	bool operator ==(const Polynomial<complex<T>>& v) {
		Polynomial<complex<T>> temp;
		list* p1 = _head;
		list* p2 = v._head;
		while (p1 || p2) {
			if (p1 && p2) {
				if (abs(p1->data - p2->data) > eps) {
					return false;
				}
				else {
					p1 = p1->link;
					p2 = p2->link;
				}
			}
			else {
				if (p1) {
					complex<T> x(0, 0);
					if (p1->data == x) {
						p1 = p1->link;
					}
					else {
						return false;
					}
				}
				if (p2) {
					complex<T> x(0, 0);
					if (p2->data == x) {
						p2 = p2->link;
					}
					else {
						return false;
					}
				}
			}
		}
		return true;
	}
	friend ostream& operator<<(ostream& os, const Polynomial<complex<T>>& v) {
		list* p = v._head;
		while (p) {
			if (p->kof == 0) {
				os << p->data << "+";
			}
			else if (!p->link) {
				os << p->data << "x^" << p->kof;
			}
			else {
				os << p->data << "x^" << p->kof << "+";
			}
			p = p->link;
		}
		return os;
	}
	Polynomial operator +(const Polynomial<complex<T>>& v) const {
		Polynomial<complex<T>> temp;
		list* p1 = _head;
		list* p2 = v._head;
		complex<T> temp_value;
		while (p1 || p2) {
			if (p1 && p2) {
				temp_value = p1->data + p2->data;
				push_back(temp, temp_value, p1->kof);
				p1 = p1->link;
				p2 = p2->link;
			}
			else {
				if (p1) {
					temp_value = p1->data;
					push_back(temp, temp_value, p1->kof);
					p1 = p1->link;
				}
				if (p2) {
					temp_value = p2->data;
					push_back(temp, temp_value, p2->kof);
					p2 = p2->link;
				}
			}
		}
		return temp;
	}
	Polynomial operator -(const Polynomial<complex<T>>& v) const {
		Polynomial<complex<T>> temp;
		list* p1 = _head;
		list* p2 = v._head;
		complex<T> temp_value;
		while (p1 || p2) {
			if (p1 && p2) {
				temp_value = p1->data - p2->data;
				push_back(temp, temp_value, p1->kof);
				p1 = p1->link;
				p2 = p2->link;
			}
			else {
				if (p1) {
					temp_value = p1->data;
					push_back(temp, temp_value, p1->kof);
					p1 = p1->link;
				}
				if (p2) {
					temp_value = -p2->data;
					push_back(temp, temp_value, p2->kof);
					p2 = p2->link;
				}
			}
		}
		return temp;
	}
	Polynomial<complex<T>>& operator *=(const complex<T>& value) {
		list* p = _head;
		while (p) {
			p->data *= value;
			p = p->link;
		}
		return *this;
	}
	Polynomial<complex<T>> operator *(const complex<T>& value) const {
		Polynomial<complex<T>> temp;
		list* p = _head;
		complex<T> value_temp;
		while (p) {
			value_temp = p->data * value;
			push_back(temp, value_temp, p->kof);
			p = p->link;
		}
		return temp;
	}
	friend Polynomial<complex<T>> operator *(const complex<T>& value, const Polynomial<complex<T>>& v) {
		Polynomial<complex<T>> temp;
		list* p = v._head;
		complex<T> value_temp;
		while (p) {
			value_temp = p->data * value;
			temp.push_back(temp, value_temp, p->kof);
			p = p->link;
		}
		return temp;
	}
	void push_back(Polynomial<complex<T>>& temp, const complex<T>& value, const size_t kof)const {
		if (temp._head) {
			list* p1 = temp._head;
			while (p1->link) {
				p1 = p1->link;
			}
			list* p2 = new list(value, kof, nullptr);
			p1->link = p2;
		}
		else
		{
			list* p1 = new list(value, kof, nullptr);
			temp._head = p1;
		}
	}
	complex<T> calculation(const complex<T>& x) const {
		complex<T> result(0, 0);
		list* p = _head;
		while (p) {
			complex<T> temp(pow(x, p->kof));
			result += p->data * temp;
			p = p->link;
		}
		return result;
	}
	void formula_Kardano() const {
		cout << "Нет действительных корней" << endl;
	}
};


void menu_1() {
	cout << "1.int" << endl;
	cout << "2.float" << endl;
	cout << "3.double" << endl;
	cout << "4.complex<float>" << endl;
	cout << "5.complex<double>" << endl;
	cout << "0.Выход" << endl;
}

void menu_2() {
	cout << "1.Установить коэффициент" << endl;
	cout << "2.Вывести многочлен" << endl;
	cout << "3.Вывести коэффициент при заданной степени" << endl;
	cout << "4.Сложение" << endl;
	cout << "5.Вычитание" << endl;
	cout << "6.Умножение на скаляр" << endl;
	cout << "7.Вычисление значения в Х" << endl;
	cout << "8.Нахождение действительных корней" << endl;
	cout << "9.Проверка на равенство" << endl;
	cout << "0.Выход" << endl;
}
int main() {
	setlocale(LC_ALL, "Russian");
	while (true) {
		system("cls");
		menu_1();
		size_t flag_1;
		cin >> flag_1;
		switch (flag_1)
		{
		case 1:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial<int> p1(number_1);
			Polynomial <int>p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2_1;
				cin >> flag_2_1;
				switch (flag_2_1)
				{
				case 1:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						int number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p1.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						int number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p2.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial <int> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_2;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						Polynomial <int> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial <int> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					size_t flag_2;
					int value;
					cout << "На что умножить?" << endl;
					cin >> value;
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						Polynomial <int> p3 = p1 * value;
						Polynomial <int> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					int value;
					cout << "Введите Х" << endl;
					cin >> value;
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					double x1, x2, x3, Q;
					try {
						p1.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						p2.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;
				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}
				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 2:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial<float> p1(number_1);
			Polynomial <float>p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2_2;
				cin >> flag_2_2;
				switch (flag_2_2)
				{
				case 1:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						float number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p1.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						float number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p2.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial <float> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_2;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						Polynomial <float> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial <float> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					size_t flag_2;
					float value;
					cout << "На что умножить?" << endl;
					cin >> value;
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						Polynomial <float> p3 = p1 * value;
						Polynomial <float> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					float value;
					cout << "Введите Х" << endl;
					cin >> value;
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					double x1, x2, x3, Q;
					try {
						p1.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						p2.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;
				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}
				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 3:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial<double> p1(number_1);
			Polynomial <double>p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2_3;
				cin >> flag_2_3;
				switch (flag_2_3)
				{
				case 1:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						double number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p1.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						double number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p2.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial <double> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_2;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						Polynomial <double> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial <double> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					size_t flag_2;
					double value;
					cout << "На что умножить?" << endl;
					cin >> value;
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						Polynomial <double> p3 = p1 * value;
						Polynomial <double> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					double value;
					cout << "Введите Х" << endl;
					cin >> value;
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					double x1, x2, x3, Q;
					try {
						p1.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						p2.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;
				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}
				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 4:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial <complex<float>> p1(number_1);
			Polynomial <complex<float>> p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2;
				cin >> flag_2;
				switch (flag_2)
				{
				case 1:
				{
					system("cls");
					size_t flag_5;
					cout << "У какого многочлена?" << endl;
					cin >> flag_5;
					switch (flag_5)
					{
					case 1:
					{
						size_t ind;
						size_t flag = 1;
						while (flag) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							float real;
							cout << "Введите действительную часть" << endl;
							cin >> real;
							float im;
							cout << "Введите мнимую часть" << endl;
							cin >> im;
							complex<float>x(real, im);
							try {
								p1.Set(x, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						size_t flag = 1;
						while (flag) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							float real;
							cout << "Введите действительную часть" << endl;
							cin >> real;
							float im;
							cout << "Введите мнимую часть" << endl;
							cin >> im;
							complex<float>x(real, im);
							try {
								p2.Set(x, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_0;
					cout << "У какого многочлена?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial<complex<float>> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_0;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						Polynomial<complex <float>> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial<complex <float>> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					size_t flag_0;
					float value;
					cout << "На что умножить?" << endl;
					cin >> value;
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						Polynomial<complex <float>> p3 = p1 * value;
						Polynomial<complex <float>> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					float real;
					cout << "Введите действительную часть" << endl;
					cin >> real;
					float im;
					cout << "Введите мнимую часть" << endl;
					cin >> im;
					complex<float> value(real, im);
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					try {
						p1.formula_Kardano();
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						p2.formula_Kardano();
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;

				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}

				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 5:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial <complex<double>> p1(number_1);
			Polynomial <complex<double>> p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2;
				cin >> flag_2;
				switch (flag_2)
				{
				case 1:
				{
					system("cls");
					size_t flag_0;
					cout << "У какого многочлена?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						size_t ind;
						size_t flag = 1;
						while (flag) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							double real;
							cout << "Введите действительную часть" << endl;
							cin >> real;
							double im;
							cout << "Введите мнимую часть" << endl;
							cin >> im;
							complex<double>x(real, im);
							try {
								p1.Set(x, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						size_t flag = 1;
						while (flag) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							float real;
							cout << "Введите действительную часть" << endl;
							cin >> real;
							float im;
							cout << "Введите мнимую часть" << endl;
							cin >> im;
							complex<double>x(real, im);
							try {
								p2.Set(x, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_0;
					cout << "У какого многочлена?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial<complex<double>> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_0;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						Polynomial<complex <double>> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial<complex <double>> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					size_t flag_0;
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						double real;
						cout << "Введите действительную часть" << endl;
						cin >> real;
						double im;
						cout << "Введите мнимую часть" << endl;
						cin >> im;
						complex<double> value(real, im);
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						double real;
						cout << "Введите действительную часть" << endl;
						cin >> real;
						double im;
						cout << "Введите мнимую часть" << endl;
						cin >> im;
						complex<double> value(real, im);
						Polynomial<complex <double>> p3 = p1 * value;
						Polynomial<complex <double>> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					double real;
					cout << "Введите действительную часть" << endl;
					cin >> real;
					double im;
					cout << "Введите мнимую часть" << endl;
					cin >> im;
					complex<double> value(real, im);
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					try {
						p1.formula_Kardano();
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						p2.formula_Kardano();
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;

				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}

				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 0:
		{
			return 0;
		}
		}
	}
}