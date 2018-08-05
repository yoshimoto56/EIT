#pragma once
#define NUM_LIGHT_MAX 9
#define NUM_CAMERA_MAX 9
#define NUM_CAMERA_DEFAULT 5
#define NUM_LIGHT_DEFAULT 3
#define GL_REAL_VIRTUAL_SCALE 0.1
#define GL_VIRTUAL_REAL_SCALE 10.0
#define GL_CAMERA_INIPOS Vector3d(35,15,0)
#define GL_CAMERA_INIPOS_MOUSE Vector3d(0,0,0)

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <GL/freeglut.h>
#include <utilities.h>
#include <glutilities.h>
#include <transfermatrix.h>
#include <glutilities.h>

namespace EITS{

	class Object{
	public:
		Object(){}
		~Object(){}
		void (*render)();
		void setObject(void (*_render)()){render=_render;}

		void (*selectObject)(Vector3d *_selected_coord, int _mode);
		void setSelectObject(void(*_selectObject)(Vector3d *_selected_coord, int _mode)){selectObject=_selectObject;}
	};

	class GLHandler
	{
	private:
		static Window window;

		static MouseSelection select;
		static Vector3d pos_clicked;
		static int key;

		static int camera_mode;
		static Camera *camera;
		static int num_camera;

		static Light *light;
		static int num_light;

		static utilities dt; 
		static bool is_run;
		static bool is_view_grid;

		static double returnDepth(int,int);
		static Vector3d returnCameraCo(Vector3d);
		static Vector3d returnWorldCo(Vector3d);

		static void startFrame();
		static void endFrame();

		static Object *object;
	public:
		GLHandler();
		~GLHandler(void);

		void init(int *_argc, char **_argv);

		void setIsRun(bool _is_run){this->is_run=_is_run;}
		bool getIsRun(){return is_run;}
		void setIsViewGrid(bool _is_view_grid){this->is_view_grid=_is_view_grid;}
		bool getIsViewGrid(){return this->is_view_grid;}

		void setKey(int _key){this->key=_key;}
		int getKey(){return key;}
		Vector3d getPosClicked(){return pos_clicked;}

		void setNumLight(int _num_light){this->num_light=_num_light;}
		int getNumLight(){return num_light;}
		Light* getLightPointer(){return light;}
		Light* getLightPointerAt(int _index){return &light[_index];}

		void setCameraMode(int _camera_mode){camera_mode=_camera_mode;}
		void setNumCamera(int _num_camera){this->num_camera=_num_camera;}
		int getNumCamera(){return num_camera;}
		Camera* getCameraPointer(){return camera;}
		Camera* getCameraPointerAt(int _index){return &camera[_index];}

		void setWindowSize(Vector2d _size){window.size=_size;}
		Vector2d getWindowSize(){return window.size;}

		void initGL();
		static void displayGL();
		static void mousePressEventGL(int _button, int _state, int _x, int _y);
		static void mouseMoveEventGL(int _x, int _y);
		static void keyboardEventGL(unsigned char _key, int _x, int _y);
		static void idle();
		static void resize(int _width, int _height);
		static void mouseWheelGL(int _wheel, int _direction, int _x, int _y);

		Object* getObject(){return object;}

	};
};
