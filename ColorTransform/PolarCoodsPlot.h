#ifndef _POLAR_COODS_PLOT_H

#define _POLAR_COODS_PLOT_H



#define _USE_MATH_DEFINES

#include <stdlib.h>

#include <math.h>

#include <string>

#include <GL/glut.h>

#include <iostream>

#include <sstream>

#include <fstream>

#include <vector>

#include "RGBTransform.h"

#include "Drawer.h"

using namespace std;



typedef vector<double> DataArray;



class PolarCoodsPlot : public Drawer {

	int dataNum;

	DataArray data_r, data_g, data_b;

	double range_r, range_g, range_b;

	double v_max, v_min;

public:

	PolarCoodsPlot(string name) {

		//b440

		ifstream fpb(name + "(WL=440nm,AOI=135deg).txt");

		v_max = -100, v_min = 100;

		double tmpb;

		dataNum = 0;

		while (fpb) {

			fpb >> tmpb;

			if (dataNum >= 90 && dataNum < 270) {

				data_b.push_back(tmpb);

				v_min = min(v_min, tmpb);

			}

			dataNum++;

		}

		dataNum = data_b.size();

		for (int i = 0; i < dataNum; i++) {

			//			data_b[i] -= v_min;

			v_max = max(v_max, data_b[i]);

		}

		range_b = 10 / v_max;



		//g550

		ifstream fpg(name + "(WL=550nm,AOI=135deg).txt");

		v_max = -100, v_min = 100;

		double tmpg;

		dataNum = 0;

		while (fpg) {

			fpg >> tmpg;

			if (dataNum >= 90 && dataNum < 270) {

				data_g.push_back(tmpg);

				v_min = min(v_min, tmpg);

			}

			dataNum++;

		}

		dataNum = data_g.size();

		for (int i = 0; i < dataNum; i++) {

			//			data_g[i] -= v_min;

			v_max = max(v_max, data_g[i]);

		}

		range_g = 10 / v_max;



		//r700

		ifstream fpr(name + "(WL=700nm,AOI=135deg).txt");

		v_max = -100, v_min = 100;

		double tmpr;

		dataNum = 0;

		while (fpr) {

			fpr >> tmpr;

			if (dataNum >= 90 && dataNum < 270) {

				data_r.push_back(tmpr);

				v_min = min(v_min, tmpr);

			}

			dataNum++;

		}

		dataNum = data_r.size();

		for (int i = 0; i < dataNum; i++) {

			//			data_r[i] -= v_min;

			v_max = max(v_max, data_r[i]);

		}

		range_r = 10 / v_max;



		cout << dataNum << endl;

		cout << v_max << endl;

		cout << v_min << endl;

	};



	void draw() {

		glClear(GL_COLOR_BUFFER_BIT);

		glColor3d(255, 255, 255);		//”wŒi‚ð”’‚Å•`‰æ

		glRectd(-1, -1, 1, 1);



		glLineWidth(1);

		glColor3d(0, 0, 0);



		//‰~‚Ì•`‰æ

		glBegin(GL_LINE_LOOP);

		for (int i = 0; i < 360; i++) {

			double rad = i / 180.0 * M_PI;

			glVertex2d(0.9 * log10(11) * cos(rad), 0.9 * log10(11) * sin(rad));

		}

		glEnd();



		//í—p‘Î”•\Ž¦

		double a5 = log10(5 + 1);

		double a1 = log10(1 + 1);

		double a05 = log10(0.5 + 1);

		double a01 = log10(0.1 + 1);



		glColor3d(0.8, 0.8, 0.8);

		glBegin(GL_LINE_LOOP);

		for (int i = 0; i < 360; i++) {

			double rad = i / 180.0 * M_PI;

			glVertex2d(0.9 * a5 * cos(rad), 0.9 * a5 * sin(rad));

		}

		glEnd();

		glBegin(GL_LINE_LOOP);

		for (int i = 0; i < 360; i++) {

			double rad = i / 180.0 * M_PI;

			glVertex2d(0.9 * a1 * cos(rad), 0.9 * a1 * sin(rad));

		}

		glEnd();

		glBegin(GL_LINE_LOOP);

		for (int i = 0; i < 360; i++) {

			double rad = i / 180.0 * M_PI;

			glVertex2d(0.9 * a05 * cos(rad), 0.9 * a05 * sin(rad));

		}

		glEnd();

		glBegin(GL_LINE_LOOP);

		for (int i = 0; i < 360; i++) {

			double rad = i / 180.0 * M_PI;

			glVertex2d(0.9 * a01 * cos(rad), 0.9 * a01 * sin(rad));

		}

		glEnd();



		//ŠDF‚Å–Ú·‚è‚ð•`‰æ	DB•\‹L

		glColor3d(0.8, 0.8, 0.8);

		for (int i = 0; i < 12; i++) {

			double rad = i / 6.0 * M_PI;

			glBegin(GL_LINE_LOOP);

			glVertex2d(0, 0);

			glVertex2d(0.9 * cos(rad), 0.9 * sin(rad));

			glEnd();

		}

		glColor3d(0, 0, 0);

		for (int i = 0; i < 5; i++) {

			double rad45 = i / 4.0 * M_PI;

			glBegin(GL_LINE_LOOP);

			glVertex2d(0, 0);

			glVertex2d(cos(rad45), sin(rad45));

			glEnd();

		}



		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "10", 0, 0.9 * log10(11));

		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "5", 0, 0.9 * a5);

		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "1", 0, 0.9 * a1);

		drawBitmapString(GLUT_BITMAP_HELVETICA_12, ".5", 0, 0.9 * a05);

		drawBitmapString(GLUT_BITMAP_HELVETICA_12, ".1", 0, 0.9 * a01);



		/*

		//Ô‚Åƒf[ƒ^‚ð•`‰æ

		glColor3d(255,0,0);

		for(int i=0; i < 5; i++){

			double v = i*v_max/5;

			double s = v*range;



			drawBitmapString(GLUT_BITMAP_HELVETICA_18, to_s(s+v_min), 0, s);

		}

		*/

		glLineWidth(4.0);

		glColor3d(255, 0, 0);

		glBegin(GL_LINE_LOOP);

		for (int i = 0; i < dataNum; i++) {

			double rad = i / 180.0 * M_PI;

			double x = 0.9 * log10(range_r * data_r[i] + 1) * cos(rad);

			double y = 0.9 * log10(range_r * data_r[i] + 1) * sin(rad);

			glVertex2d(x, y);

		}

		glEnd();



		glLineWidth(3.0);

		glColor3d(0, 255, 0);

		glBegin(GL_LINE_LOOP);

		for (int i = 0; i < dataNum; i++) {

			double rad = i / 180.0 * M_PI;

			double x = 0.9 * log10(range_g * data_g[i] + 1) * cos(rad);

			double y = 0.9 * log10(range_g * data_g[i] + 1) * sin(rad);

			glVertex2d(x, y);

		}

		glEnd();



		glLineWidth(3.0);

		glColor3d(0, 0, 255);

		glBegin(GL_LINE_LOOP);

		for (int i = 0; i < dataNum; i++) {

			double rad = i / 180.0 * M_PI;

			double x = 0.9 * log10(range_b * data_b[i] + 1) * cos(rad);

			double y = 0.9 * log10(range_b * data_b[i] + 1) * sin(rad);

			glVertex2d(x, y);

		}

		glEnd();

	};



};



#endif //_POLAR_COODS_PLOT_H
