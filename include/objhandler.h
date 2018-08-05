#pragma once

#include <fstream>
#include <iostream>
#include <GL/freeglut.h>
#include <vectormatrix.h>
#include <mesh.h>

namespace EITS
{
	class ObjMesh: public SurfMesh
	{
	protected:
		bool loadOBJFile(const char *filename);
		bool loadMTLFile(const char *filename);
		char filename_obj[256];
		char filename_mtl[256];
		char directoryname[256];
	public:
		ObjMesh(SurfMesh *parent = 0){}
		~ObjMesh(){}
		char* getName(){return this->filename_obj;}
		bool load(const char *_filename);
		bool save(const char *_filename);
		void setup();
	};
};