#include<iostream>
#include<iomanip>
#include<fstream>
#include <vector>
#include <cmath>
using namespace std;

double abs(double a, double b) {							//取两数之差的绝对值
	if (a > b)	return a - b;
	else return b - a;
}
double getD(vector<double> list) {							//用于计算D 
	double result = 0;
	int n = list.size();

	for (int i = 0; i < n; i++) {
		result += list[i];
	}

	result /= n * list[n - 1];
	return result;
}


int main() {
	fstream reader;
	reader.open("SYS.txt", ios::in);
	double d1, d2, d3, xl, xr, xm, xi, a, b, f,D = 0;			//初始化参数
	int n = 0;
	vector<double> t;										//用double类型的vector存储数据ti
	while (!reader.eof()) {
		reader >> d1;
		reader >> d2;
		reader >> d3;
		t.push_back(d3);									//仅读入第三列数据，即累加故障间隔时间
		if (reader.fail())	break;
	}
	reader.close();
	t.pop_back();
	n = t.size();
	D = getD(t);											//步骤1，计算得到D并判断它与0.5的大小关系
	if (D >= 0.5) {
		cout << "参数估计无解，停止计算。" << endl;
		exit(0);
	}
	xl = (1 - 2 * D) / 2;
	xr = 1 / D;
	xm = (xr + xl) / 2;										//步骤2，计算得到xl和xr
	xi = 0.1;											//设定xi

	while (abs(xr, xl) > xi) {							//步骤3，判断绝对值与xi的大小关系后循环
		f = (1 - D * xm) * exp(xm) + (D - 1) * xm - 1;
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
				break;										//以上两种情况均不是，跳转步骤4
			}
		}
	}
	b = xm / t[n - 1];										//步骤4，计算得到a，b
	a = n / (1 - exp(-(b * t[n - 1])));

	cout << "a	" << setprecision(20) << a << endl;			//输出计算结果
	cout << "b	" << setprecision(20) << b << endl;
	return 0;

}
