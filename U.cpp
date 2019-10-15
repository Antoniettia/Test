#include<iostream>
#include<iomanip>
#include<fstream>
#include <vector>
#include <cmath>
using namespace std;

double abs(double a, double b) {									//取两数之差的绝对值
	if (a > b)	return a - b;
	else return b - a;
}
double getD(vector<double> list, int j) {							//累加计算D
	double result = 0;

	for (int i = 0; i < j; i++) {
		result += list[i];
	}

	result /= j * list[j - 1];
	return result;
}

double getft(vector <double> t, double a, double b, int i, double M) {			//依照公式，计算Fi
	double f = 0;
	f = 1 - exp(a * (exp(-(b * (t[i - 1] + M))) - exp(-(b * M))));
	return f;
}

int main() {
	fstream reader;											//初始化参数
	reader.open("SYS.txt", ios::in);
	double d1, d2, d3, M, n1, xr, xl, xm, a, b, f, xi;
	int n;
	vector<double> t;										//用double类型的vector存储数据ti
	vector <double> ft;										//用于存储横坐标
	vector <double> D;										//用于存储各组数据对应的D
	vector <double> len;									//用于存储纵坐标

	while (!reader.eof()) {									//读入ti到vector t中
		reader >> d1;
		reader >> d2;
		reader >> d3;
		t.push_back(d3);
		if (reader.fail())	break;
	}

	reader.close();
	t.pop_back();
	n = t.size();
	M = t[n - 1];

	for (int i = 1; i <= n; i++) {
		D.push_back(getD(t, i));								//存入每组数据得出的D
		cout << "D	" << D[i - 1] << endl;

		if (D[i - 1] >= 0.5) {									//G-O模型步骤1
			cout << "参数估计无解，停止计算。" << endl;
			continue;
		}

		xl = (1 - 2 * D[i - 1]) / 2;							//G-O模型步骤2
		xr = 1 / D[i - 1];
		xm = (xr + xl) / 2;
		xi = 0.0001;

		while (abs(xr, xl) > xi) {
			f = (1 - D[i - 1] * xm) * exp(xm) + (D[i - 1] - 1) * xm - 1;
			if (f > xi) {										//G-O模型步骤3
				xl = xm;
				xm = (xr + xl) / 2;
			}
			else {
				if (f < (-xi)) {
					xr = xm;
					xm = (xr + xl) / 2;
				}
				else {
					break;
				}
			}
		}

		b = xm / t[i - 1];										//G-O模型步骤4
		a = i / (1 - exp(-(b * t[i - 1])));

		cout << endl << i << "    a	" << setprecision(16) << a;
		cout << "    b	" << setprecision(16) << b << endl << endl;

		ft.push_back(getft(t, a, b, i, M));						//得到每一组数据的ft并存入vector ft中

	}

	n1 = ft.size();												//升序排列ft得到uj
	for (int i = 0; i < n1; i++) {
		for (int j = i; j < n1; j++) {
			if (ft[i] > ft[j]) {
				double temp = ft[i];
				ft[i] = ft[j];
				ft[j] = temp;
			}
		}
	}

	f = 1;
	for (int i = 0; i < n1; i++) {								//依照公式得到纵坐标并输出坐标信息
		double l = (i + f) / n1;
		len.push_back(l);
		cout << i + 1 << "	横坐标为" << setprecision(16) << ft[i] << "		纵坐标为";
		cout << setprecision(16) << len[i] << endl;
	}

	return 0;

}
