#pragma once

#include <process.h>
#include <windows.h>
#include <iostream>

#include <vectormatrix.h>
#include <transfermatrix.h>
#include <glhandler.h>
#include <glutilities.h>
#include <stlhandler.h>
#include <objhandler.h>
#include <voxelhandler.h>
#include <modelhandler.h>
#include <interface.h>

EITS::GLHandler *gl;
EITS::StlMesh *stl;
EITS::VolumeMesh *fem;
EITS::ObjMesh *obj;
EITS::VoxModel *vox;

int argc;
char **argv;

//Key interface
void thread1(LPVOID pParam)
{
	std::cout<<"HELLO:)"<<std::endl;
	std::cout<<"********************* MeshViewer v1.0.0 *********************"<<std::endl;
	std::cout<<"Command list: "<<std::endl;
	std::cout << "[LOAD DATA] load:filename(*.obj, *.stl, *.fem, *.vox)" << std::endl;
	std::cout << "[SAVE DATA] save:filename(*.stl, *.vox, (0:x,1:y,2:z):*.slice, *.fem)" << std::endl;
	std::cout << "[CONVERT DATA] cvt:(obj2stl, stl2obj, stl2fem, obj2vox:#res, stl2vox:#res)" << std::endl;
	std::cout << "[EXIT] exit" << std::endl;
	std::cout<<"*******************************************************"<<std::endl;

	char command[256];
	char t_filename[256];
	while(true){
		std::cin>>command;
		if(strstr(command,"load:")){
			gl->setIsRun(false);
			Sleep(100);
			sscanf_s(command,"load:%s", t_filename, 256);
			if(strstr(command,".stl")){
				if(stl->load(t_filename)){
					stl->centerize();
					stl->getInfo(std::cout);
				}
			}
			else if (strstr(command, ".fem")) {
				if (fem->load(t_filename)) {
					fem->getInfo(std::cout);
				}
			}
			else if (strstr(command, ".obj")) {
				if (obj->load(t_filename)) {
				}
			}
			else if (strstr(command, ".vox")) {
				if (vox->load(t_filename)) {
					vox->setSlicePos(EITS::Vector3i(vox->getWidth() / 2.0, vox->getHeight() / 2.0, vox->getDepth() / 2.0));
					vox->setSlice();
					vox->setIsAutoScale(true);
				}
			}
			else {
				std::cout << "-l:INVALID OPTION" << std::endl;
			}
			gl->setIsRun(true);
		}
		else if (strstr(command, "save:")) {
			gl->setIsRun(false);
			Sleep(100);
			sscanf_s(command, "save:%s", t_filename, 256);

			if (strstr(command, ".obj")) {
				//TODO
				//obj->save(t_filename);
			}
			else if (strstr(command, ".stl")) {
				stl->save(t_filename);
			}
			else if (strstr(command, ".vox")) {
				vox->save(t_filename);
			}
			else if (strstr(command, ".fem")) {
				if (fem->save(t_filename)) {
				}
			}
			else if (strstr(command, ".slice")) {
				int t_axis;
				sscanf_s(command, "save:%d:%s", &t_axis, t_filename, (unsigned int)sizeof(t_filename));
				if (t_axis == EITS::X_AXIS || t_axis == EITS::Y_AXIS || t_axis == EITS::Z_AXIS)
					vox->saveSlice(t_filename, t_axis);
				else {
					std::cout << "save:slice:INVALID OPTION" << std::endl;
				}
			}
			else {
				std::cout << "save:INVALID OPTION" << std::endl;
			}

			gl->setIsRun(true);
		}
		else if (strstr(command, "cvt:")) {
			gl->setIsRun(false);
			Sleep(100);
			if (strstr(command, "stl2fem")) {
				EITS::convertSTL2FEM(stl, fem);
			}
			else if (strstr(command, "fem2stl")) {
				//EITS::cnver:convertFME2STL(stl, fem);
			}
			else if (strstr(command, "obj2stl")) {
				EITS::convertOBJ2STL(obj, stl);
			}
			else if (strstr(command, "stl2obj")) {
				EITS::convertSTL2OBJ(stl, obj);
			}
			else if (strstr(command, "obj2vox:")) {
				int num = 0;
				sscanf_s(command, "cvt:obj2vox:%d", &num);
				if (num>0 && num <= 1024) {
					if (EITS::convertOBJ2VOX(obj, vox, num)) {
						vox->setSlicePos(EITS::Vector3i(num / 2, num / 2, num / 2));
						vox->setSlice();
					}
				}
			}
			else if (strstr(command, "stl2vox:")) {
				int num = 0;
				sscanf_s(command, "cvt:stl2vox:%d", &num);
				if (num>0 && num <= 1024) {
					if (EITS::convertSTL2VOX(stl, vox, num)) {
						vox->setSlicePos(EITS::Vector3i(num / 2, num / 2, num / 2));
						vox->setSlice();
					}
				}
				else {
					std::cout << "cvt:stl2vox:INVALID RESOLUTION" << std::endl;
				}
			}
			else {
				std::cout << "cvt:INVALID OPTION" << std::endl;
			}
			gl->setIsRun(true);
		}
		else if (strstr(command, "prc:")) {
			gl->setIsRun(false);
			Sleep(100);
			if (strstr(command, "fill")) {
				std::cout << "Filling color inside object...";
				if (vox->getIsLoaded()) {
					vox->setLabel();
					vox->setGradMapWithinSurf();
					vox->setSlice();
					std::cout << "[OK]" << std::endl;
				}
				else std::cout << "[FAIL]" << std::endl;
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
		obj->render();
		stl->render();
		fem->render();
		vox->render();
		vox->renderGuid();
		vox->renderSlicer();
	}
}

void selectObject(EITS::Vector3d *_selected_coord, int _mode)
{
	EITS::select(_selected_coord, _mode, EITS::SELECT_NODE, stl);
	EITS::select(_selected_coord, _mode, EITS::SELECT_FACET, stl);

//	stl->select(_selected_coord, _mode);
}

int main(int _argc, char *_argv[])
{
	obj = new EITS::ObjMesh;
	stl = new EITS::StlMesh;
	fem = new EITS::VolumeMesh;
	vox = new EITS::VoxModel;

	argc=_argc;
	argv=_argv;

	gl = new EITS::GLHandler; 
	gl->getObject()->setObject(object);
	gl->getObject()->setSelectObject(selectObject);
	gl->setIsViewGrid(true);

//	fem->is_view.label = true;
//	fem->is_view.facet = false;

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
	delete obj;
	delete vox;
}

