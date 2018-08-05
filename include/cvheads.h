#pragma once
#include <opencv/cv.h>
#include <opencv/cxcore.h>
#include <opencv/highgui.h>
#include <opencv/ml.h>
#include <opencv2/core/core.hpp>
#include <opencv2/video/tracking.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
//#include <opencv2/nonfree/features2d.hpp>
#include <opencv/cxcore.h>
//#include <opencv2/nonfree/nonfree.hpp>
//#include <opencv2/legacy/legacy.hpp>

#define CV_VERSION_STR CVAUX_STR(CV_MAJOR_VERSION) CVAUX_STR(CV_MINOR_VERSION) CVAUX_STR(CV_SUBMINOR_VERSION)
#ifdef _DEBUG
#define CV_EXT_STR "d.lib"
#else
#define CV_EXT_STR ".lib"
#endif
#pragma comment(lib, "opencv_core"            CV_VERSION_STR CV_EXT_STR)
#pragma comment(lib, "opencv_highgui"        CV_VERSION_STR CV_EXT_STR)
#pragma comment(lib, "opencv_imgproc"  CV_VERSION_STR CV_EXT_STR)
#pragma comment(lib, "opencv_imgcodecs"  CV_VERSION_STR CV_EXT_STR)
#pragma comment(lib, "opencv_calib3d"  CV_VERSION_STR CV_EXT_STR)
//#pragma comment(lib, "opencv_gpu"   CV_VERSION_STR CV_EXT_STR)
#pragma comment(lib, "opencv_videoio"   CV_VERSION_STR CV_EXT_STR)
#pragma comment(lib, "opencv_objdetect"  CV_VERSION_STR CV_EXT_STR)
//#pragma comment(lib, "opencv_features2d" CV_VERSION_STR CV_EXT_STR)
//#pragma comment(lib, "opencv_nonfree" CV_VERSION_STR CV_EXT_STR)
//#pragma comment(lib, "opencv_flann"   CV_VERSION_STR CV_EXT_STR)
//#pragma comment(lib, "opencv_ffmpeg"  CV_VERSION_STR CV_EXT_STR)
//#pragma comment(lib, "opencv_ts"   CV_VERSION_STR CV_EXT_STR)
//#pragma comment(lib, "opencv_contrib"  CV_VERSION_STR CV_EXT_STR)
//#pragma comment(lib, "opencv_ml"   CV_VERSION_STR CV_EXT_STR)
//#pragma comment(lib, "opencv_legacy"  CV_VERSION_STR CV_EXT_STR)

namespace cv{
	static Size SIZE_QUADVGA = Size(1280, 960);
	static Size SIZE_SUPERVGA = Size(800, 600);
	static Size SIZE_VGA = Size(640, 480);
	static Size SIZE_QUARTERVGA = Size(320, 240);

	static Scalar SCALAR_BLAK = Scalar(0, 0, 0);
	static Scalar SCALAR_WHITE = Scalar(255, 255, 255);
	static Scalar SCALAR_RED = Scalar(0, 0, 255);
	static Scalar SCALAR_GREEN = Scalar(0, 255, 0);
	static Scalar SCALAR_BLUE = Scalar(255, 0, 0);
};


enum{
	IMAGE,
	CAMERA
};