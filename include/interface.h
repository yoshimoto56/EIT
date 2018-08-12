#pragma once

#include <GL/freeglut.h>
#include <pointlineface.h>
#include <voxelhandler.h>
#include <stlhandler.h>
#include <objhandler.h>
#include <mesh.h>

namespace EITS {
	void select(Vector3d *_coord, int _mode, int _type, StlMesh *_stl);
	void select(Vector3d *_coord, int _mode, int _type, VolumeMesh *_fem);
};

