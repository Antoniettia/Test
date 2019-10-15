#include<iostream>
#include<iomanip>
#include<fstream>
#include <vector>
#include <cmath>
using namespace std;

double abs(double a, double b) {							//取两数之差绝对值
	if (a > b)	return a - b;
	else return b - a;
}
double getD(vector<double> list, int j) {					//累加计算得到D
	double result = 0;

	for (int i = 0; i < j; i++) {
		result += list[i];
	}

	result /= j * list[j - 1];

	return result;
}
double getx(vector <double> t, double a, double b, int i, double N) {			//化简公式后依照公式计算xi
	double x = 0;
	x = -a * (exp(-(b * (t[i - 1] + N))) - exp(-(b * N)));
	return x;
}

int main() {
	fstream reader;											//初始化参数
	reader.open("SYS.txt", ios::in);
	double d1, d2, d3, N, xl, xr, xm, xi, a, b, f, sum = 0, n1;
	vector<double> t;										//用double类型的vector存储数据ti
	vector <double> D;										//用于存储多个D
	vector <double> x, x_sum;								//用于存储计算中间值x的数值，x累加的数值
	vector <double> y;										//用于存储横坐标
	vector <double> len;									//用于存储纵坐标
	int n;

	while (!reader.eof()) {									//读入ti
		reader >> d1;
		reader >> d2;
		reader >> d3;
		t.push_back(d3);
		if (reader.fail())	break;
	}
	reader.close();
	t.pop_back();

	n = t.size();
	N = t[n - 1];

	for (int i = 1; i <= n; i++) {						//G-O模型
		D.push_back(getD(t, i));
		cout << "D	" << D[i - 1] << endl;
		if (D[i - 1] >= 0.5) {							//步骤1
			cout << "参数估计无解，停止计算。" << endl;
			continue;
		}

		xl = (1 - 2 * D[i - 1]) / 2;					//步骤2
		xr = 1 / D[i - 1];
		xm = (xr + xl) / 2;
		xi = 0.0001;

		while (abs(xr, xl) > xi) {
			f = (1 - D[i - 1] * xm) * exp(xm) + (D[i - 1] - 1) * xm - 1;		//步骤3
			if (f > xi) {
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
		b = xm / t[i - 1];							//步骤4
		a = i / (1 - exp(-(b * t[i - 1])));
		cout << endl << i << "    a	" << setprecision(16) << a;
		cout << "    b	" << setprecision(16) << b << endl << endl;

		x.push_back(getx(t, a, b, i, N));				//调用函数计算得到xi并存入vector x中

	}

	n1 = x.size();
	for (int i = 0; i < n1; i++) {						//升序排列x
		for (int j = i; j < n1; j++) {
			if (x[i] > x[j]) {
				double temp = x[i];
				x[i] = x[j];
				x[j] = temp;
			}
		}
		sum += x[i];									//计算并存储xj的累加值
		x_sum.push_back(sum);
	}

	f = 1;
	for (int i = 0; i < n1; i++) {
		double l = x_sum[i] / x_sum[n1 - 1];			//累加后相除得到yj

		y.push_back(l);
		l = (i + f) / n1;								//计算得到纵坐标
		len.push_back(l);

		cout << i + 1 << "	横坐标为" << setprecision(16) << y[i] << "		纵坐标为";			//输出坐标信息
		cout << setprecision(16) << len[i] << endl;
	}
	return 0;

}
