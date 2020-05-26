#ifndef _GRAPH_H_

#define _GRAPH_H_

#include"Drawer.h"

#include<vector>

#include<complex>

#include<string>

#include<map>

using namespace std;



typedef complex<double> Complex;

class Graph :public Drawer {

	vector< DataSet<double> > Solid, Dash;

	Range<double> X, Y;

	double offset;

public:

	Graph();

	Graph(int dataNum);

	~Graph();

	void draw();

	void draw_data(DataSet<double>, bool f);

	void draw_axis();

	void GraphVertex(Complex);

	void GraphVertex(double, double);

	void addSolidData(string);

	void addLineData(string);

};



#endif //_GRAPH_H_
