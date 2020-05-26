
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

// �ϐ��̐錾

//----------------------------------------------------

const int WindowPositionX = 200;  //��������E�B���h�E�ʒu��X���W

const int WindowPositionY = 200;  //��������E�B���h�E�ʒu��Y���W

//const int WindowWidth = 512 + 200;      //��������E�B���h�E�̕�

//const int WindowHeight = 512;     //��������E�B���h�E�̍���

const char WindowTitle[] = "�V�~�����[�V����";  //�E�B���h�E�̃^�C�g��



//-----------------------------------------------

//�V�~�����[�V�����֌W�̕ϐ�

//----------------------------------------------



Drawer* drawer;

//----------------------------------------------------

// �֐��v���g�^�C�v�i��ɌĂяo���֐����ƈ����̐錾�j

//----------------------------------------------------

void Initialize(void);   //�����ݒ莞�ɌĂяo���֐�

void Idle(void);         //�A�C�h�����ɌĂяo���֐�

void Display(void);      //��ʕ`�掞�ɌĂяo���֐�



//----------------------------------------------------

// ���C���֐�

//----------------------------------------------------



int main(int argc, char* argv[]) {

    cout.setf(ios::fixed, ios::floatfield);	//�Œ菬���_, 18�����x�Ŏw��

    cout.precision(18);



    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);	//���������[�N���o�p

    glutInit(&argc, argv);                                     //���̏�����

    glutInitWindowPosition(WindowPositionX, WindowPositionY);  //�E�B���h�E�̈ʒu�̎w��

    glutInitWindowSize(WINDOW_W, WINDOW_H);					//�E�B���h�E�T�C�Y�̎w��

    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE); //�f�B�X�v���C���[�h�̎w��

    glutCreateWindow(WindowTitle);                             //�E�B���h�E�̍쐬

    glutIdleFunc(Idle);                                        //�v���O�����A�C�h����Ԏ��ɌĂяo�����֐�

    glutDisplayFunc(Display);                                  //�`�掞�ɌĂяo�����֐����w�肷��i�֐����FDisplay�j

    Initialize();                                              //�����ݒ�̊֐����Ăяo��

    glutMainLoop();



    return 0;

}



//----------------------------------------------------

// �����ݒ�̊֐�

//-----//---------------------------------------------------------------------------------------------------

void Initialize() {

    //drawer = new PolarCoodsPlot("(500nm,0deg)NTF(1deg).txt"); //�ɍ��W�v���b�g

   // drawer = new PolarCoodsPlot("./rgbData/HairSimulation/incidence/"); //�ɍ��W�v���b�g

    //drawer = new RGBTransform("Red_Refraction.txt"); 

    //drawer = new Graph();

    drawer = new RGBPlot();  //�p���˗�����, �F������

      //drawer = new PolarCoodsPlot("(WL=490nm,AOI=-80deg).txt");



}

//----------------------------------------------------

// �A�C�h�����ɌĂяo�����֐�



void Idle() {

    Sleep(10);

    glutPostRedisplay(); //glutDisplayFunc()���P����s����

}



//----------------------------------------------------

// �`��̊֐�

//----------------------------------------------------

void Display(void) {

    glClear(GL_COLOR_BUFFER_BIT);

    glColor3d(255, 255, 255);		//�w�i�𔒂ŕ`��

    glRectd(-1, -1, 1, 1);

    drawer->draw();



    glutSwapBuffers();

}

