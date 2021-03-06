#pragma once
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
//#include <GL/glut.h>
#include <GL/freeglut.h>
#include <math.h>
#include <vectormatrix.h>
#include <transfermatrix.h>
#include <colormap.h>

#define VC_INIT_ZOOM 20.0
#define VC_DEFAULT_ZOOM 20.0


namespace EITS
{
	enum MState{
		MOUSE_NONE,
		MOUSE_PUSH,
		MOUSE_RELEASE
	};

	struct Window{
		char title[256];
		double fps;
		Vector2d size;
	};

	class Material{
	private:
		char name[256];
		Vector4f color;
		Vector4f ambient;
		Vector4f diffuse;
		Vector4f specular;
		Vector4f emission;
		float shininess;
		int illum;
	public:
		Material(char _name[256] = "default"){
			strcpy_s(name, _name);
			color = Vector4f(0.0, 0.0, 0.0, 1.0);
			ambient = Vector4f(0.0, 0.0, 0.0, 1.0);
			diffuse = Vector4f(1.0, 1.0, 1.0, 1.0);
			specular = Vector4f(0.5, 0.5, 0.5, 1.0);
			emission = Vector4f(0.0, 0.0, 0.0, 1.0);
			shininess = 5.0f;
			illum=4;
		}
		~Material(){}

		Material &operator = (const Material &_mat){ 
			strcpy_s(name, _mat.name); 
			color=_mat.color; 
			ambient=_mat.ambient;
			diffuse=_mat.diffuse;
			specular=_mat.specular;
			emission=_mat.emission;
			shininess=_mat.shininess;
			illum=_mat.illum;
			return (*this); }
		void setName(const char *_name){strcpy_s(name,_name);}
		void setColor(Vector4f _color){color = _color;}
		void setAmbient(Vector4f _ambient){ambient = _ambient;}
		void setDiffuse(Vector4f _diffuse){diffuse = _diffuse;}
		void setSpecular(Vector4f _specular){specular = _specular;}
		void setEmission(Vector4f _emission){emission = _emission;}
		void setShininess(float _shininess){shininess = _shininess;}
		void setColor(Vector3f _color, float _alpha=1.0){for(int i=0;i<3;i++){color.X[i] = _color.X[i];}color.w=_alpha;}
		void setAmbient(Vector3f _ambient, float _alpha=1.0){for(int i=0;i<3;i++){ambient.X[i] = _ambient.X[i];}ambient.w=_alpha;}
		void setDiffuse(Vector3f _diffuse, float _alpha=1.0){for(int i=0;i<3;i++){diffuse.X[i] = _diffuse.X[i];}diffuse.w=_alpha;}
		void setSpecular(Vector3f _specular, float _alpha=1.0){for(int i=0;i<3;i++){specular.X[i] = _specular.X[i];}specular.w=_alpha;}
		void setEmission(Vector3f _emission, float _alpha=1.0){for(int i=0;i<3;i++)emission.X[i] = _emission.X[i];emission.w=_alpha;}
		void setColorA(float _a){color.w=_a;}
		void setAmbientA(float _a){ambient.w=_a;}
		void setDiffuseA(float _a){diffuse.w=_a;}
		void setEmissionA(float _a){emission.w=_a;}
		void setSpecularA(float _a){specular.w=_a;}
		void setIllum(int _illum){this->illum=_illum;}

		char* getName(){return name;}
		Vector4f getColor(){return color;}
		Vector4f getAmbient(){return ambient;}
		Vector4f getDiffuse(){return diffuse;}
		Vector4f getEmission(){return emission;}
		Vector4f getSpecular(){return specular;}
		float getShininess(){return shininess;}
		int getIllum(){return illum;}

		void clear();
		void enable();
		void set();
		void disable();
	};

	class Light{
	private:
		Vector4f position;
		Vector4f ambient;
		Vector4f diffuse;
		Vector4f specular;
		int id;
		bool is_on;
	public:
		Light();
		~Light();

		void init();
		void set();
		void enable();
		void disable();

		void setIsOnLight(bool _is_on){this->is_on=_is_on;}
		void setID(int _id){this->id=_id;}
		void setPosition(Vector4f _position){position = _position;}
		void setAmbient(Vector4f _ambient){ambient = _ambient;}
		void setDiffuse(Vector4f _diffuse){diffuse = _diffuse;}
		void setSpecular(Vector4f _specular){specular = _specular;}

		bool getIsOn(){return this->is_on;}
		int getID(){return this->id;}
		Vector4f getAmbient(){return this->ambient;}
		Vector4f getDiffuse(){return this->diffuse;}
		Vector4f getSpecular(){return this->specular;}
		Vector4f getPosition(){return this->position;}
	};

	class MouseButton
	{
	public:
		Vector2f before;
		Vector2f current;
		Vector2f after;
		MState state;
		MouseButton();
		~MouseButton();
		void reset();
		void update();
	};
	class MouseSelection
	{
	protected:
		bool is_draw_select_region;
		Vector2d left_top;
		Vector2d right_bottom;
	public:
		MouseSelection();
		~MouseSelection();
		void setIsDrawSelectRegion(bool _is_draw_select_region){this->is_draw_select_region=_is_draw_select_region;}
		void setLeftTop(Vector2d _left_top){this->left_top=_left_top;}
		void setRightBottom(Vector2d _right_bottom){this->right_bottom=_right_bottom;}
		void render();
		bool getIsDrawSelectRegion(){return this->is_draw_select_region;}
	};

	enum CAM_MODE{
		CAM_PERSP,
		CAM_FRONT,
		CAM_SIDE,
		CAM_TOP
	};

	class Camera
	{
	private:
		double zoom;
		Vector3d angle;
		Vector3d position;
		Vector3d trans;
		Vector3d target;
		Vector3d upvector;
		transferMatrixd Tproj;
		transferMatrixd Tmodel;
		transferMatrixd Trot;
		transferMatrixd Ttran;

		int model;

		void setModelViewMatrix(bool _is_scale=true);

	public:

		Camera();
		~Camera();
		void init(int _cameraModel=CAM_PERSP);

		void setModel(int _model=CAM_PERSP){this->model=_model;}
		void resizePersp(Vector2d _size);
		void resizeOrtho(Vector2d _size, bool _is_scale=true);
		void attachTransMat(bool _is_fix=false);
		void attachRotMat(bool _is_scale=true);
		void attach(Vector2d _size, bool _is_scale=true);

		MouseButton rml[3];
		void rmlReset();
		void rmlUpdate();
		void update();
		void MouseInput(int _button, Vector2f _position, bool ml=true, bool mm=true, bool mr=true);
		void MouseMotion(Vector2f _position, bool ml=true, bool mm=true, bool mr=true);

		void setZoom(double _zoom){zoom=_zoom;}
		void setAngle(Vector3d _angle){angle=_angle;}
		void setPosition(Vector3d _position){position=_position;}
		void setTrans(Vector3d _trans){trans=_trans;}
		void setTarget(Vector3d _target){target=_target;}
		void setUpvector(Vector3d _upvector){upvector=_upvector;}
		void setTrot(transferMatrixd _Trot){this->Trot=_Trot;}
		void setTtran(transferMatrixd _Ttran){this->Ttran=_Ttran;}
		double getZoom(){return zoom;}
		Vector3d getPosition(){return position;}
		Vector3d getTrans(){return trans;}
		Vector3d getTarget(){return target;}
		Vector3d getUpvector(){return upvector;}
		Vector3d getAngle(){return angle;}

		void setTproj(transferMatrixd _Tproj){Tproj=_Tproj;}
		void setTmodel(transferMatrixd _Tmodel){Tmodel=_Tmodel;}
	};

	void glutSolidRectangle(Vector3d _size, bool _is_vert_shade=false);
	void glutLineArrows(double _scale=10);
	void glutAxsis();
	void glutAxsis(Vector3d _nx, Vector3d _ny, Vector3d _nz);
	void glutGridGround(double _size, int _num_grid=10);
	void glutDrawTextAt(std::string _text, Vector3d _position);
};