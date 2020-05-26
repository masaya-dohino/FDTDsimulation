
#define _CRTDBG_MAP_ALLOC

#define _USE_MATH_DEFINES

#include <stdlib.h>

#include <crtdbg.h>

#include <math.h>

#include <stdio.h>

#include <string>

#include <GL/glut.h>

#include <iostream>

#include <sstream>

#include <fstream>

#include <windows.h>

#include"Drawer.h"

#include "PolarCoodsPlot.h"

#include "RGBTransform.h"

#include "Graph.h"

#include"RGBPlot.h"

using namespace std;



//----------------------------------------------------

// 変数の宣言

//----------------------------------------------------

const int WindowPositionX = 200;  //生成するウィンドウ位置のX座標

const int WindowPositionY = 200;  //生成するウィンドウ位置のY座標

//const int WindowWidth = 512 + 200;      //生成するウィンドウの幅

//const int WindowHeight = 512;     //生成するウィンドウの高さ

const char WindowTitle[] = "シミュレーション";  //ウィンドウのタイトル



//-----------------------------------------------

//シミュレーション関係の変数

//----------------------------------------------



Drawer* drawer;

//----------------------------------------------------

// 関数プロトタイプ（後に呼び出す関数名と引数の宣言）

//----------------------------------------------------

void Initialize(void);   //初期設定時に呼び出す関数

void Idle(void);         //アイドル時に呼び出す関数

void Display(void);      //画面描画時に呼び出す関数



//----------------------------------------------------

// メイン関数

//----------------------------------------------------



int main(int argc, char* argv[]) {

    cout.setf(ios::fixed, ios::floatfield);	//固定小数点, 18桁制度で指定

    cout.precision(18);



    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);	//メモリリーク検出用

    glutInit(&argc, argv);                                     //環境の初期化

    glutInitWindowPosition(WindowPositionX, WindowPositionY);  //ウィンドウの位置の指定

    glutInitWindowSize(WINDOW_W, WINDOW_H);					//ウィンドウサイズの指定

    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE); //ディスプレイモードの指定

    glutCreateWindow(WindowTitle);                             //ウィンドウの作成

    glutIdleFunc(Idle);                                        //プログラムアイドル状態時に呼び出される関数

    glutDisplayFunc(Display);                                  //描画時に呼び出される関数を指定する（関数名：Display）

    Initialize();                                              //初期設定の関数を呼び出す

    glutMainLoop();



    return 0;

}



//----------------------------------------------------

// 初期設定の関数

//-----//---------------------------------------------------------------------------------------------------

void Initialize() {

    //drawer = new PolarCoodsPlot("(500nm,0deg)NTF(1deg).txt"); //極座標プロット

   // drawer = new PolarCoodsPlot("./rgbData/HairSimulation/incidence/"); //極座標プロット

    //drawer = new RGBTransform("Red_Refraction.txt"); 

    //drawer = new Graph();

    drawer = new RGBPlot();  //角反射率から, 色をだす

      //drawer = new PolarCoodsPlot("(WL=490nm,AOI=-80deg).txt");



}

//----------------------------------------------------

// アイドル時に呼び出される関数



void Idle() {

    Sleep(10);

    glutPostRedisplay(); //glutDisplayFunc()を１回実行する

}



//----------------------------------------------------

// 描画の関数

//----------------------------------------------------

void Display(void) {

    glClear(GL_COLOR_BUFFER_BIT);

    glColor3d(255, 255, 255);		//背景を白で描画

    glRectd(-1, -1, 1, 1);

    drawer->draw();



    glutSwapBuffers();

}

