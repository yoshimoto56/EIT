#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <GL/freeglut.h>
#include <vectormatrix.h>

namespace EITS{

	//TYPE OF COLOR MAP
	enum {
		MAP_NORMAL,
		MAP_JET,
		MAP_GRAY
	};
	class ColorMap{
	private:
		Vector3d color;
		double min,max;
	public:
		ColorMap();
		~ColorMap();
		//set minimum and maximum values
		void setParam(double _min, double _max);
		//opengl color setting
		void glColorMap(double, int _type = MAP_NORMAL);
		//get color vector from value
		Vector3d getColorMap(double _vale, int _type = MAP_NORMAL);
		//get color vector from value for opencv
		Vector3d getColorMapForCV(double _vale, int _type = MAP_NORMAL);
		//get normalized value from -1 to 1.
		double getNormalizedValue(double _value);
	};

	//functions for HSV and RGB color space convertion
	Vector3d hsv2rgb(Vector3d _hsv);
	void minMax(Vector3d _rgb, double *_min, double *_max);
	Vector3d rgb2hsv(Vector3d _rgb);

	//function for jet color
	double interpolate( double val, double y0, double x0, double y1, double x1 );
	double base( double val );
	double red( double gray );
	double green( double gray );
	double blue( double gray );
	//value = [-1:1]
	Vector3d getJetColor(double _value);
};