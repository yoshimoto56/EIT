#pragma once
#include <vectormatrix.h>
#include <mesh.h>

namespace EITS{
	class Facet;

	//点と面の距離を計算するMethod
	double getFacePointDistance(Vector3d _point, Vector3d _normal, Vector3d _coord);
	//点を面に垂直に下ろした交点を計算するMethod
	Vector3d getFacePointProjection(Vector3d _point, Vector3d _lNormal, Vector3d _sNormal, Vector3d _coord);
	Vector3d getFacePointProjection(Vector3d _point, Vector3d _normal, Vector3d _coord);
	//線と面の交点を計算するMethod
	Vector3d getIntersectPointFace(Vector3d *_point, Vector3d _normal, Vector3d _coord);
	//2点で与えられる平行面の間に点があるかを計算するMethod
	bool isFacesSandPoint(Vector3d _coord1, Vector3d _coord2, Vector3d _point);
	//2点が面をはさんでいるかどうか計算するMethod
	bool isViewNodesSandFace(Vector3d *_point, Vector3d _normal, Vector3d _coord);
	//ある点が3角形の内部にあることを計算するMethod
	bool isViewNodeOnTriangle(Vector3d _point, Vector3d _normal, Vector3d *_coord);
	//TomasMollerの交差判定
	bool isIntersectLineTriangle(Vector3d *_point, Vector3d *_coord, Vector3d *_intersection);
	//面の内側か外側かを判定するMethod
	bool isViewNodeWithinTriangle(Vector3d _point, Vector3d _normal, Vector3d _coord);
	//点が三角形内に投影されるかを判定するMethod
	bool isProjectedPointOnTriangle(Vector3d _point, Vector3d _normal, Vector3d *_coord);
	Vector3d getProjectedPointOnTriangle(Vector3d _point, Vector3d _normal, Vector3d _coord);
	//点が有限面上にあるかを判定するMethod
	bool isViewNodeOnFiniteFace(Vector3d _point, Facet *_facet);
	bool isProjectedPointOnFiniteFace(Vector3d _point, Facet *_facet);
	//点と線の最短距離を求めるMethod
	double getPointLineDistance(Vector3d _point, Vector3d _dir, Vector3d _coord);
};
