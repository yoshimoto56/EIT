#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <GL/freeglut.h>
#include <vectormatrix.h>

namespace EITS{
	class ColorMap{
	private:
		Vector3d colormap[256];
		Vector3d currentColor;
		double min,max;
		int offset,rate;
		int rgb[3];
	public:
		ColorMap();
		~ColorMap();
		void setParam(double _min, double _max, int _offset, int _rate=200);
		bool load(const char*);
		void setColorMap();
		Vector3d glColorMap(double);
		Vector3d getColorMap(double);
		int* getColorMapip(double);
		double getNormalizedValue(double _value);
	};

	Vector3d HSV2RGB(Vector3d HSV);
	void MinMax(Vector3d RGB, double *_min, double *_max);
	Vector3d RGB2HSV(Vector3d RGB);

	double interpolate( double val, double y0, double x0, double y1, double x1 );
	double base( double val );
	double red( double gray );
	double green( double gray );
	double blue( double gray );
	Vector3d getJetColor(double _value);
};