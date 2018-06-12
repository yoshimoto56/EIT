#pragma once
#include <windows.h>


namespace EITS
{
	class utilities
	{
	private:
		LARGE_INTEGER freq;
		LARGE_INTEGER start;
		LARGE_INTEGER end;
	public:
		utilities(void);
		~utilities(void);
		void setStartTime();
		void setEndTime();
		double getDeltaTime();
	};
	void movingAverage(int _width, double *_input, double *_output, int _size);
};