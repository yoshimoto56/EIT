#pragma once

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <GL/freeglut.h>
#include <vectormatrix.h>
#include <glutilities.h>
#include <mesh.h>

namespace EITS
{
	enum{
		STL_ASCII,
		STL_BINARY
	};
	class StlMesh: public SurfMesh
	{
	public:
		StlMesh(SurfMesh *parent=0){
			this->num_material=1;this->newMaterial();
		}
		~StlMesh(){this->deleteMaterial();}
		bool load(const char *_filename);
		bool loadASCII(const char *_filename);
		bool loadBinary(const char *_filename);
		bool save(const char *_filename, int _type=STL_BINARY);
		bool saveASCII(const char *_filename);
		bool saveBinary(const char *_filename);
		void setup();
	};
};