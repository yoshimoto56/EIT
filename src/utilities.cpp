#include "utilities.h"

namespace EITS
{

	utilities::utilities(void)
	{
	}

	utilities::~utilities(void)
	{
	}

	void utilities::setStartTime()
	{
		QueryPerformanceFrequency(&freq);
		QueryPerformanceCounter(&start);
	}

	void utilities::setEndTime()
	{
		QueryPerformanceCounter(&end);
	}

	double utilities::getDeltaTime()
	{
		return (double)(end.QuadPart-start.QuadPart)/(double)freq.QuadPart;
	}


	void movingAverage(int _width, double *_input, double *_output, int _size)
	{
		if(_width < _size && _width > 1){
			for(int i = 0; i < _size; i++){
				_output[i]=0;
				for(int w =0;w < _width;w++){
					if(i-w>=0)
						_output[i] += _input[i-w];
				}
				if(i >= _width-1)
					_output[i] /= _width;
				else
					_output[i] /=  i + 1;
			}
		}
	}

};