#pragma once
#include <vector> 
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <random>

#include <vectormatrix.h>
#include <pointlineface.h>
#include <mesh.h>
#include <utilities.h>
#include <colormap.h>
#include <glutilities.h>

namespace EITS {
	enum NMODE
	{
		NONE,
		FIXED,
		UNFIXED,
		CONTACT,
		NONCONTACT,
		GROUNDED,
		UNGROUNDED,
		ELECTRODE,
		SENSING,
		DETECTOR
	};
}