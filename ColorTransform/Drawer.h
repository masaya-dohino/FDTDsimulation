#ifndef _DRAWER_H

#define _DRAWER_H

#include<string>

#include<vector>

#include<math.h>

#include<fstream>

#include<complex>

#include<algorithm>

#include<iostream>


using namespace std;

const int WINDOW_W = 700;

const int WINDOW_H = 700;

class MyColor;



class Drawer {

public:

	virtual void draw() = 0;

};



void drawBitmapString(void* font, string s, double x, double y);



template<class T> string to_s(const T& n) {

	string rt;

	stringstream ss;

	ss << n;

	ss >> rt;

	return rt;

}



class MyColor {

public:

	double r, g, b;

	MyColor(double _r, double _g, double _b) :r(_r), g(_g), b(_b)

	{

	}

	MyColor() :r(0), g(0), b(0)

	{

	}

};



MyColor operator+ (MyColor& a, MyColor& b);



MyColor operator- (MyColor& a, MyColor& b);



MyColor operator* (const double& a, const MyColor& b);



MyColor operator/ (const MyColor& b, const double& a);

bool operator== (MyColor& a, MyColor& b);





template<class T> class Range {

	T _min, _max, _interval;

	bool _minF, _maxF;

public:

	Range(T a, T b, T inter) :_interval(inter)

	{

		_min = min(a, b);

		_max = max(a, b);

		_minF = _maxF = true;

		cout << _minF << endl;



	}



	Range()

	{

		_minF = _maxF = false;

	}



	T getMax()

	{

		if (_maxF)

			return _max;

		else

			throw "not max";

	}



	T getMin()

	{

		if (_maxF)

			return _min;

		else

			throw "not min";

	}



	T getInterVal() {

		return _interval;

	}



	void setMax(T a) {

		if (_maxF) {

			_max = max(_max, a);

		}

		else {

			_max = a;

			_maxF = true;

		}

		if (!_maxF) cout << "max---" << endl;

	}



	void setMin(T a) {

		if (_minF) {

			_min = min(_min, a);

		}

		else {

			_min = a;

			_minF = true;

		}

		if (!_minF) cout << "min---" << endl;

	}



	void set(T a)

	{

		setMin(a);

		setMax(a);

	}

};


namespace ColorList {

	const MyColor Red = MyColor(1.0, 0.0, 0.0);

	const MyColor Green = MyColor(0.0, 1.0, 0.0);

	const MyColor Blue = MyColor(0.0, 0.0, 1.0);

	const MyColor White = MyColor(1.0, 1.0, 1.0);

	const MyColor Black = MyColor(0.0, 0.0, 0.0);

	MyColor getColor(int);

};

template<class T> class DataSet {

	static unsigned int ID;

public:

	vector< complex<T> > Data;	//todo Xで昇順になるようにすべき

	Range<double> X, Y;

	string dataName;

	MyColor  color;

	unsigned int dataNum;



	void Initialize() {

		X = Range<T>();

		Y = Range<T>();

		color = ColorList::getColor(ID++);

	}



	DataSet() :dataNum(0)

	{

		Initialize();

	}



	DataSet(string name) :dataNum(0)

	{

		Initialize();

		ifstream fp(name);	if (!fp) throw "no file";

		T tmpX, tmpY;

		while (fp) {

			fp >> tmpX >> tmpY;

			setData(tmpX, tmpY);

		}

	}



	//X座標のデータあり

	DataSet(string name, T* X) :dataNum(0)

	{

		Initialize();

		ifstream fp(name);	if (!fp) throw "no file";

		T tmp;

		for (int i = 0; fp; i++) {

			fp >> tmp;

			setData(X[i], tmp);

		}

	}



	//X座標の範囲あり

	DataSet(string name, Range<T> _X) {

		Initialize();

		X = _X;

		ifstream fp(name);	if (!fp) throw "no file";

		T tmp;



		for (int i = 0; fp; i++) {

			fp >> tmp;

			setData(_X.getMin() + i * _X.getInterVal(), tmp);

		}

	}



	void setName(string n)

	{

		dataName = n;

	}



	string getName()

	{

		return dataName;

	}



	void setData(T x, T y)

	{

		Data.push_back(complex<T>(x, y));

		X.set(x);

		Y.set(y);

		dataNum++;

	}



	Range<T> getRangeX() {

		return X;

	}



	Range<T> getRangeY() {

		return Y;

	}

};




#endif
