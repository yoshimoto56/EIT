#pragma once
#include <pointlineface.h>
#include <triangleboxtest.h>
#include <voxelhandler.h>
#include <stlhandler.h>
#include <objhandler.h>
#include <mesh.h>
#include <tetgen.h>

namespace EITS
{
	bool convertSTL2VOX(StlMesh *_stl, VoxModel *_vox, int _num=64);
	bool convertOBJ2VOX(ObjMesh *_obj, VoxModel *_vox, int _num=64);
	bool convertSTL2OBJ(StlMesh *_stl, ObjMesh *_obj);
	bool convertOBJ2STL(ObjMesh *_obj, StlMesh *_stl);
	bool convertSTL2FEM(StlMesh *_stl, VolumeMesh *_fem, double _val = 1.2);
	bool convertSTL2TetGen(StlMesh *_stl, tetgenio *_in);
	bool convertTetGen2FEM(tetgenio *_out, VolumeMesh *_fem);

	void testIntersection(Facet *_mesh, VoxModel *_vox, int _map=255, int _thickness=3);
	void testIntersection(Facet *_mesh, VoxModel *_vox, Vector3d _color=Vector3d(0,0,1), int _thickness=3);
	bool calOBJandVOXcollision(ObjMesh *_obj, VoxModel *_vox, int *_index); 
//	double calBLADEandVOXcollision(transferMatrixd _Tw2o, BLADE *_blade, VoxModel *_vox, int *_index); 
//	bool calBLADEandVOXcollision(transferMatrixd _Tw2o, BLADE *_blade, VoxModel *_vox, int *_index, Vector3d *_points, int *_num_points); 
//	double calOBJandVOXdistance(ObjMesh *_obj, VoxModel *_vox, int *_index); 
	void millingVOXbyOBJ(ObjMesh *_obj, VoxModel *_vox);
	void millingVOXbyPOINTS(VoxModel *_vox, Vector3d *_points, double *depth, int _num_points=1);
//	void clippingVOXbyPLANE(VoxModel *_vox, PLANE *_plane, int _numPlane=1);

	void spanLabelVOX2OBJ(VoxModel *_vox, ObjMesh *_obj);
	void setLabelOBJNearest(ObjMesh *_obj);

	namespace ICP{
		enum MatchingType{
			ICP2P,ICP2S
		};
		transferMatrixd getOptimalRegistration(ObjMesh *_objA, ObjMesh *_objB, int _mode=ICP2P, double minError=0.00001);
		void matching(int _num_point, Vector3d *_pointA, Vector3d *_pointB, ObjMesh *_objB, int _mode);
		double errorMetric(double *_error, Vector3d *_pointB, ObjMesh *_objA);
		Vector3d centering(int _num_point, Vector3d *_point);
		Vector3d translating(int _num_point, Vector3d *_pointA, Vector3d *_pointB);
		transferMatrixd rotating(int _num_point, Vector3d *_pointA, Vector3d *_pointB);
		transferMatrixd minimizing(int _num_point, Vector3d *_pointA, Vector3d *_pointB);
	};
	transferMatrixd calTransMatBetween3(Vector3d *_p1, Vector3d *_p2);
	double calLongAxisFromPointCloud(int _num, Vector3d *_point, int *_index);
};
