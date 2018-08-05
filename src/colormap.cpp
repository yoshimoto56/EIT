#include "colormap.h"

namespace EITS{

	ColorMap::ColorMap()
	{
		min=0;
		max=1;
	}

	ColorMap::~ColorMap()
	{
	}

	void ColorMap::setParam(double _min, double _max)
	{
		min=_min;
		max=_max;
	}

	void ColorMap::glColorMap(double _value, int _type)
	{
		glColor3dv(getColorMap(_value, _type).X);
	}

	Vector3d ColorMap::getColorMap(double _value, int _type)
	{
		if (_type == MAP_NORMAL)
			color = hsv2rgb(Vector3d(360 * (0.7 - 0.7 * (_value - min) / (max - min)), 1, 1));
		else if (_type == MAP_JET)
			color = getJetColor( getNormalizedValue(_value) );
		else {
			color.x = color.y = color.z = 0.5+0.5*getNormalizedValue(_value);
		}
		return color;
	}
	Vector3d ColorMap::getColorMapForCV(double _value, int _type)
	{
		if (_type == MAP_NORMAL) {
			color.x = hsv2rgb(Vector3d(360 * (0.7 - 0.7 * (_value - min) / (max - min)), 1, 1)).z;
			color.y = hsv2rgb(Vector3d(360 * (0.7 - 0.7 * (_value - min) / (max - min)), 1, 1)).y;
			color.z = hsv2rgb(Vector3d(360 * (0.7 - 0.7 * (_value - min) / (max - min)), 1, 1)).x;
		}
		else if (_type == MAP_JET) {
			color.x = getJetColor(getNormalizedValue(_value)).z;
			color.y = getJetColor(getNormalizedValue(_value)).y;
			color.z = getJetColor(getNormalizedValue(_value)).x;
		}
		else {
			color.x = color.y = color.z = 0.5 + 0.5*getNormalizedValue(_value);
		}
		return 255.0 * color;
	}
	double ColorMap::getNormalizedValue(double _value)
	{
		double t_value = (_value - min) / (max - min);
		if(t_value<0)t_value = 0;
		if(t_value>1)t_value = 1;
		return -1+2*t_value;
	}

	Vector3d hsv2rgb(Vector3d _hsv) {
		Vector3d rgb;
		int H_i;
		double f, p, q, t;
		if (fabs(_hsv.x)>360)
			_hsv.x = _hsv.x - ((int)(_hsv.x / 360)) * 360;
		if (_hsv.x<0)
			_hsv.x = 360 - _hsv.x;

		H_i = (int)floor(_hsv.X[0] / 60.0) % 6;
		f = _hsv.X[0] / 60.0 - (double)H_i;
		p = _hsv.X[2] * (1 - _hsv.X[1]);
		q = _hsv.X[2] * (1 - f*_hsv.X[1]);
		t = _hsv.X[2] * (1 - (1 - f)*_hsv.X[1]);
		if (H_i == 0) { rgb.X[0] = _hsv.X[2]; rgb.X[1] = t; rgb.X[2] = p; }
		else if (H_i == 1) { rgb.X[0] = q; rgb.X[1] = _hsv.X[2]; rgb.X[2] = p; }
		else if (H_i == 2) { rgb.X[0] = p; rgb.X[1] = _hsv.X[2]; rgb.X[2] = t; }
		else if (H_i == 3) { rgb.X[0] = p; rgb.X[1] = q; rgb.X[2] = _hsv.X[2]; }
		else if (H_i == 4) { rgb.X[0] = t; rgb.X[1] = p; rgb.X[2] = _hsv.X[2]; }
		else if (H_i == 5) { rgb.X[0] = _hsv.X[2]; rgb.X[1] = p; rgb.X[2] = q; }
		return rgb;
	}

	void minMax(Vector3d _rgb, double *_min, double *_max){
		*_min=*_max=0.0;
		for(int i=0;i<3;i++){
			if(*_min>_rgb.X[i])*_min=_rgb.X[i];
			if(*_max<_rgb.X[i])*_max=_rgb.X[i];
		}
	}

	Vector3d rgb2hsv(Vector3d _rgb){
		Vector3d hsv;
		double min,max;
		minMax(_rgb, &min, &max);
		if(max!=min){
			for(int i=0;i<3;i++){
				if(max==_rgb.X[i]){
					hsv.X[0]=60.0*(_rgb.X[(i+1)%3]-_rgb.X[(i+2)%3])/(max-min)+120.0*(double)i;
				}
			}
			if(hsv.X[0]<0)hsv.X[0]+=360.0;
		}
		else hsv.X[0]=-1;
		hsv.X[1]=max-min;
		hsv.X[2]=max;
		return hsv;
	}

	double interpolate( double val, double y0, double x0, double y1, double x1 ) {
		return (val-x0)*(y1-y0)/(x1-x0) + y0;
	}

	double base( double val ) {
		if ( val <= -0.75 ) return 0;
		else if ( val <= -0.25 ) return interpolate( val, 0.0, -0.75, 1.0, -0.25 );
		else if ( val <= 0.25 ) return 1.0;
		else if ( val <= 0.75 ) return interpolate( val, 1.0, 0.25, 0.0, 0.75 );
		else return 0.0;
	}

	double red( double gray ) {
		return base( gray - 0.5 );
	}

	double green( double gray ) {
		return base( gray );
	}

	double blue( double gray ) {
		return base( gray + 0.5 );
	}

	Vector3d getJetColor(double _value)
	{
		return Vector3d(red(_value),green(_value),blue(_value));
	}
}