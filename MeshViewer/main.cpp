#pragma once

#include <process.h>
#include <windows.h>
#include <iostream>

#include <vectormatrix.h>
#include <transfermatrix.h>
#include <glhandler.h>
#include <glutilities.h>
#include <stlhandler.h>
#include <modelhandler.h>

EITS::GLHandler *gl;
EITS::StlMesh *stl;
EITS::VolumeMesh *fem;

int argc;
char **argv;
int cur_label = 0;

//Key interface
void thread1(LPVOID pParam)
{
	std::cout<<"HELLO:)"<<std::endl;
	std::cout<<"********************* MeshViewer v0.0.0 *********************"<<std::endl;
	std::cout<<"Command list"<<std::endl;
	std::cout << "[LOAD DATA] load:filename(*.stl, *.fem)" << std::endl;
	std::cout << "[SAVE DATA] save:filename(*.stl, *.fem)" << std::endl;
	std::cout << "[CONVERT DATA] cvt:(stl2fem))" << std::endl;
	std::cout << "[EXIT] exit" << std::endl;
	std::cout<<"*******************************************************"<<std::endl;

	char command[256];
	char t_filename[256];
	while(true){
		std::cin>>command;
		if(strstr(command,"load:")){
			gl->setIsRun(false);
			Sleep(10);
			sscanf_s(command,"load:%s", t_filename, 256);
			if(strstr(command,".stl")){
				if(stl->load(t_filename)){
					stl->centerize();
					stl->getInfo(std::cout);
				}
			}
			if (strstr(command, ".fem")) {
				if (fem->load(t_filename)) {
					fem->getInfo(std::cout);
				}
			}
			gl->setIsRun(true);
		}
		else if (strstr(command, "save:")) {
			gl->setIsRun(false);
			Sleep(10);
			sscanf_s(command, "save:%s", t_filename, 256);
			if (strstr(command, ".stl")) {
				if (stl->save(t_filename)) {
				}
			}
			if (strstr(command, ".fem")) {
				if (fem->save(t_filename)) {
				}
			}
			gl->setIsRun(true);
		}
		else if (strstr(command, "cvt:")) {
			gl->setIsRun(false);
			Sleep(10);
			if (strstr(command, "stl2fem")) {
				EITS::convertSTL2FEM(stl, fem);
			}
			gl->setIsRun(true);
		}
		else if(strstr(command,"exit")){
			gl->setIsRun(false);
			std::cout<<"GOOD BYE:("<<std::endl;
			Sleep(2000);
			exit(0);
		}
		else {
			std::cout<<"UNDEFINED COMMAND"<<std::endl;
		}
	}
}


//Window
void thread2(LPVOID pParam)
{
	gl->init(&argc,argv);
	glutInitWindowSize(400 , 300);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition(100 , 50);
	glutCreateWindow(argv[0]);
	glutDisplayFunc(gl->displayGL);
	glutMotionFunc(gl->mouseMoveEventGL);
	glutMouseFunc(gl->mousePressEventGL);
	glutKeyboardFunc(gl->keyboardEventGL);
	glutMouseWheelFunc( gl->mouseWheelGL );
	glutIdleFunc(gl->idle);
	glutReshapeFunc(gl->resize);

	gl->initGL();
	glutMainLoop();
}


void object()
{
	if(gl->getIsRun()){
		stl->render();
		fem->render();
	}
}

void selectObject(EITS::Vector3d *_selected_coord, int _mode)
{
	//stl->select(_selected_coord, _mode);
}

int main(int _argc, char *_argv[])
{
	stl = new EITS::StlMesh;
	stl->setIsAutoScale(false);
	stl->setIsViewNode(true);
	stl->setIsIdentity(true);

	fem = new EITS::VolumeMesh;
	fem->setIsAutoScale(false);
	fem->setIsViewLine(true);
	fem->setIsViewFacet(true);

	argc=_argc;
	argv=_argv;
	gl = new EITS::GLHandler; 
	gl->getObject()->setObject(object);
	gl->getObject()->setSelectObject(selectObject);
	gl->setis_viewGrid(true);

	HANDLE hMutex;
	HANDLE hThread[2];
	hMutex = CreateMutex(NULL,FALSE,NULL);
	hThread[0] = (HANDLE)_beginthread(thread1,0,NULL);
	hThread[1] = (HANDLE)_beginthread(thread2,0,NULL);
	WaitForMultipleObjects(2,hThread,TRUE,INFINITE);
	CloseHandle(hMutex);

	delete gl;
	delete stl;
	delete fem;
}

