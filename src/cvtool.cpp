#include "cvtool.h"

namespace CVTOOL{
	/*=== 関数 ===*/
	cv::Point warpPoint(cv::Point _p, cv::Mat _G)
	{
#ifdef _DEBUG
		if(_G.empty()){
			puts("Error: warpPoint");
			puts("- empty G");
			return cv::Point();
		}
#endif

		double h = _G.at<double>(2,0)*_p.x+_G.at<double>(2,1)*_p.y+_G.at<double>(2,2);
		return cv::Point(
					(_G.at<double>(0,0)*_p.x+_G.at<double>(0,1)*_p.y+_G.at<double>(0,2)) / h,
					(_G.at<double>(1,0)*_p.x+_G.at<double>(1,1)*_p.y+_G.at<double>(1,2)) / h
					);
	}
	cv::Point derivation(cv::Mat _img, cv::Point _p, cv::Mat _G, DERIVATION_TYPE _TYPE)
	{
#ifdef _DEBUG
		if(_img.empty()){
			puts("Warning: derivation");
			puts("- empty image");
			return cv::Point(0,0);
		}
#endif
		cv::Rect roi(0, 0, _img.size().width, _img.size().height);
		cv::Point p = warpPoint(_p, _G);
		if(!p.inside(roi)){
			return cv::Point(0, 0);
		}
		cv::Point drv = cv::Point(0,0);
		switch(_TYPE){
			case DERIVATION_TYPE::SOBEL:
				cv::Point p1 = warpPoint(_p + cv::Point(-1, -1), _G);
				cv::Point p2 = warpPoint(_p + cv::Point(0, -1), _G);
				cv::Point p3 = warpPoint(_p + cv::Point(1, -1), _G);
				cv::Point p4 = warpPoint(_p + cv::Point(-1, 0), _G);
				cv::Point p6 = warpPoint(_p + cv::Point(1, 0), _G);
				cv::Point p7 = warpPoint(_p + cv::Point(-1, 1), _G);
				cv::Point p8 = warpPoint(_p + cv::Point(0, 1), _G);
				cv::Point p9 = warpPoint(_p + cv::Point(1, 1), _G);
				if(p1.inside(roi) && p2.inside(roi) && p3.inside(roi) && p4.inside(roi) && p6.inside(roi) && p7.inside(roi) && p8.inside(roi) && p9.inside(roi)){
					drv.x = _img.at<uchar>(p3) + 2*_img.at<uchar>(p6) + _img.at<uchar>(p9)
							- _img.at<uchar>(p1) - 2*_img.at<uchar>(p4) - _img.at<uchar>(p7);
					drv.y = _img.at<uchar>(p7) + 2*_img.at<uchar>(p8) + _img.at<uchar>(p9)
							- _img.at<uchar>(p1) - 2*_img.at<uchar>(p2) - _img.at<uchar>(p3);
				}
				break;
		}
		return drv;
	}
	cv::Point derivation(cv::Mat _img, cv::Point _p, DERIVATION_TYPE _TYPE)
	{
#ifdef _DEBUG
		if(_img.empty()){
			puts("Warning: derivation");
			puts("- empty image");
			return cv::Point(0,0);
		}
#endif
		cv::Rect roi(0, 0, _img.size().width, _img.size().height);
		if(!_p.inside(roi)){
			return cv::Point(0, 0);
		}
		cv::Point drv = cv::Point(0,0);
		switch(_TYPE){
			case DERIVATION_TYPE::SOBEL:
				cv::Point p1 = _p + cv::Point(-1, -1);
				cv::Point p2 = _p + cv::Point(0, -1);
				cv::Point p3 = _p + cv::Point(1, -1);
				cv::Point p4 = _p + cv::Point(-1, 0);
				cv::Point p6 = _p + cv::Point(1, 0);
				cv::Point p7 = _p + cv::Point(-1, 1);
				cv::Point p8 = _p + cv::Point(0, 1);
				cv::Point p9 = _p + cv::Point(1, 1);
				if(p1.inside(roi) && p2.inside(roi) && p3.inside(roi) && p4.inside(roi) && p6.inside(roi) && p7.inside(roi) && p8.inside(roi) && p9.inside(roi)){
					drv.x = _img.at<uchar>(p3) + 2*_img.at<uchar>(p6) + _img.at<uchar>(p9)
							- _img.at<uchar>(p1) - 2*_img.at<uchar>(p4) - _img.at<uchar>(p7);
					drv.y = _img.at<uchar>(p7) + 2*_img.at<uchar>(p8) + _img.at<uchar>(p9)
							- _img.at<uchar>(p1) - 2*_img.at<uchar>(p2) - _img.at<uchar>(p3);
				}
				break;
		}
		return drv;
	}
	cv::Mat getHomography(cv::Mat _x, double _n)
	{
#ifdef _DEBUG
		if((_x.size().width * _x.size().height)!=8){
			puts("Error: getHomography");
			puts("- invalid x size");
			return cv::Mat();
		}
#endif

		//パラメータ_xが表すリー代数のsl(3)の要素計算:A(x)=SUM xi*Gi
		cv::Mat_<double> A;
		A = _x.at<double>(0) * G1+
			_x.at<double>(1) * G2+
			_x.at<double>(2) * G3+
			_x.at<double>(3) * G4+
			_x.at<double>(4) * G5+
			_x.at<double>(5) * G6+
			_x.at<double>(6) * G7+
			_x.at<double>(7) * G8;

		//H(x)=exp(A(x))の計算：指数関数行列なのでテーラー展開で計算
		cv::Mat_<double> H=(cv::Mat_<double>(3,3)<<
							1, 0, 0,
							0, 1, 0,
							0, 0, 1);
		cv::Mat_<double> exponentA = (cv::Mat_<double>(3,3)<<
							1, 0, 0,
							0, 1, 0,
							0, 0, 1);
		double factorialI = 1.0;
		for(int i = 1; i <= _n; i++){
			exponentA *= A;
			factorialI *= i;
			H += exponentA / factorialI;
		}
		return H;
	}

	cv::Mat trimImageWithHomography(cv::Size _size, cv::Mat _src, cv::Mat _H)
	{
		if(_size.width < 2 || _size.height < 2 || _src.empty() || _H.cols != 3 || _H.rows != 3) return cv::Mat();
		cv::Mat dst;
		if(_src.channels() == 3){
			dst = cv::Mat_<cv::Vec3b>(_size);
		}else{
			dst = cv::Mat_<uchar>(_size);
		}
		cv::Point2i p_src;
		int i, j;
		double h ;
		int t_count = 0;
		int area = _size.area();
		for(j = 0; j < _size.height; j ++){
			for(i = 0; i < _size.width; i ++){
				//ホモグラフィ変換：関数より直接書いたほうが速い
				//p_src = cvTransformFunction::warpPoint(p_tmp, _H);
				h = _H.at<double>(2, 0) * i + _H.at<double>(2, 1) * j + _H.at<double>(2, 2);
				p_src.x = (_H.at<double>(0,0) * i + _H.at<double>(0,1) * j + _H.at<double>(0, 2)) / h;
				p_src.y = (_H.at<double>(1,0) * i + _H.at<double>(1,1) * j + _H.at<double>(1, 2)) / h;
				if(0 <= p_src.x && p_src.x < _src.cols && 0 <= p_src.y && p_src.y < _src.rows){
					*(uchar*)(dst.data + dst.step * j + dst.channels() * i)	= *(uchar*)(_src.data + _src.step * p_src.y + _src.channels() * p_src.x);
					if(_src.channels() == 3){
						*(uchar*)(dst.data + dst.step * j + dst.channels() * i + 1) = *(uchar*)(_src.data + _src.step * p_src.y + _src.channels() * p_src.x +1);
						*(uchar*)(dst.data + dst.step * j + dst.channels() * i + 2) = *(uchar*)(_src.data + _src.step * p_src.y + _src.channels() * p_src.x +2);
					}
					t_count++;
				}
			}
		}
		if(t_count<1)return cv::Mat();
		return dst;
	}

	void drawHistgram(cv::Mat *_img_src, std::vector<double> _hist)
	{
		cv::Scalar color(255,255,255);
		for(int i = 1; i < _hist.size(); i ++){
			cv::line(*_img_src, cv::Point((i-1) * _img_src->cols / _hist.size(), _img_src->rows * (1 - _hist[i - 1]))
								,cv::Point(i * _img_src->cols / _hist.size(), _img_src->rows * (1 - _hist[i])), color);
		}
	}
	cv::Point2f getMoment(cv::Mat _img_src, cv::Point2i _dim)
	{
		double total = 0;
		cv::Point2f m = cv::Point2f(0,0);
		for(int j = 0; j < _img_src.rows; j ++){
			for(int i = 0; i < _img_src.cols ; i ++){
				if(_img_src.channels() == 3){
					m.x += pow(i, (double)_dim.x) * (double)_img_src.at<cv::Vec3b>(j,i)[0];
					m.y += pow(i, (double)_dim.y) * (double)_img_src.at<cv::Vec3b>(j,i)[0];
					total +=(double)_img_src.at<cv::Vec3b>(j,i)[0];
				}
				else{
					m.x += pow(i, (double)_dim.x) * (double)_img_src.at<uchar>(j,i);
					m.y += pow(i, (double)_dim.y) * (double)_img_src.at<uchar>(j,i);
					total +=(double)_img_src.at<uchar>(j,i);
				}
			}
		}
		if(total != 0)m *= 1.0/total;
		return m;
	}

	std::vector<cv::Mat> getHSV(cv::Mat _img_src)
	{
		std::vector<cv::Mat> img_hsv(3);
		img_hsv[0] = cv::Mat_<uchar>(_img_src.rows, _img_src.cols);
		img_hsv[1] = cv::Mat_<uchar>(_img_src.rows, _img_src.cols);
		img_hsv[2] = cv::Mat_<uchar>(_img_src.rows, _img_src.cols);
		if(_img_src.channels() == 3){
			for(int j = 0; j < _img_src.rows; j ++){
				for(int i = 0; i < _img_src.cols ; i ++){
					img_hsv[0].at<uchar>(j,i) = _img_src.at<cv::Vec3b>(j,i)[0];
					img_hsv[1].at<uchar>(j,i) = _img_src.at<cv::Vec3b>(j,i)[1];
					img_hsv[2].at<uchar>(j,i) = _img_src.at<cv::Vec3b>(j,i)[2];
				}
			}
		}
		return img_hsv;
	}

	cv::Mat getNormalizedImage(cv::Mat _img_src, int _max)
	{
		cv::Mat dst;
		_img_src.copyTo(dst);
		cv::Vec3b maxval(0,0,0);
		cv::Vec3b minval(255,255,255);
		
		for(int j = 0; j < _img_src.rows; j ++){
			for(int i = 0; i < _img_src.cols ; i ++){
				if(dst.channels() == 3){
					if(maxval[0] < _img_src.at<cv::Vec3b>(j, i)[0]){
						maxval[0] = _img_src.at<cv::Vec3b>(j, i)[0];
					}
					if(maxval[1] < _img_src.at<cv::Vec3b>(j, i)[1]){
						maxval[1] = _img_src.at<cv::Vec3b>(j, i)[1];
					}
					if(maxval[2] < _img_src.at<cv::Vec3b>(j, i)[2]){
						maxval[2] = _img_src.at<cv::Vec3b>(j, i)[2];
					}
					if(_img_src.at<cv::Vec3b>(j, i)[0] != 0 && minval[0] > _img_src.at<cv::Vec3b>(j, i)[0]){
						minval[0] = _img_src.at<cv::Vec3b>(j, i)[0];
					}
					if(_img_src.at<cv::Vec3b>(j, i)[1] != 0 && minval[1] > _img_src.at<cv::Vec3b>(j, i)[1]){
						minval[1] = _img_src.at<cv::Vec3b>(j, i)[1];
					}
					if(_img_src.at<cv::Vec3b>(j, i)[2] != 0 && minval[2] > _img_src.at<cv::Vec3b>(j, i)[2]){
						minval[2] = _img_src.at<cv::Vec3b>(j, i)[2];
					}
				}else{
					if(maxval[0] < _img_src.at<uchar>(j, i)){
						maxval[0] = _img_src.at<uchar>(j, i);
					}
					if(_img_src.at<uchar>(j, i) != 0 && minval[0] > _img_src.at<uchar>(j, i)){
						minval[0] = _img_src.at<uchar>(j, i);
					}
				}
			}
		}
		for(int j = 0; j < _img_src.rows; j ++){
			for(int i = 0; i < _img_src.cols ; i ++){
				if(dst.channels() == 3){
					if(_img_src.at<cv::Vec3b>(j, i)[0] != 0)
						dst.at<cv::Vec3b>(j, i)[0] = std::max(0.0, std::min((double)_max, (double)_max * (dst.at<cv::Vec3b>(j, i)[0] - minval[0]) / (double)(maxval[0] - minval[0])));
					if(_img_src.at<cv::Vec3b>(j, i)[1] != 0)
						dst.at<cv::Vec3b>(j, i)[1] = std::max(0.0, std::min((double)_max, (double)_max * (dst.at<cv::Vec3b>(j, i)[1] - minval[1]) / (double)(maxval[1] - minval[1])));
					if(_img_src.at<cv::Vec3b>(j, i)[2] != 0)
						dst.at<cv::Vec3b>(j, i)[2] = std::max(0.0, std::min((double)_max, (double)_max * (dst.at<cv::Vec3b>(j, i)[2] - minval[2]) / (double)(maxval[2] - minval[2])));
				}else{
					if(_img_src.at<uchar>(j, i) != 0)
						dst.at<uchar>(j, i) = std::max(0.0, std::min((double)_max, (double)_max * (dst.at<uchar>(j, i) - minval[0]) / (double)(maxval[0] - minval[0])));
				}
			}
		}
		return dst;
	}

	std::vector<double> getHistgram(cv::Mat _img_src)
	{
		std::vector<double> dst(256, 0);
		if(_img_src.empty())return dst;

		//ヒストグラムの計算
		for(int j = 0; j < _img_src.rows; j ++){
			for(int i = 0; i < _img_src.cols; i ++){
				if(_img_src.channels() ==3){
					if(_img_src.at<cv::Vec3b>(j,i)[2] != 0)
						dst[ _img_src.at<cv::Vec3b>(j,i)[2] ] ++; 
				}else{
					if(_img_src.at<uchar>(j,i) != 0)
						dst[ _img_src.at<uchar>(j,i) ] ++; 

				}
			}
		}

		//正規化
		double max = *std::max_element(dst.begin(),dst.end());
		for(int i = 0; i < 256; i ++)
			dst[i] /= max;

		return dst;
	}

	void fastResize(cv::Mat *_src, cv::Mat *_dst, cv::Size _size){
		*_dst=cv::Mat(_size, _src->type());
		int i,x1,y1,x2,y2;
#ifdef _OPENMP
#pragma omp parallel private(i, x1, y1, x2, y2)
#pragma omp for
#endif
		for(i=0;i<_size.height*_size.width;i++){
			x1 = i % _size.width; x2 = x1 * ((double)_src->cols/_size.width);
			y1 = i / _size.width; y2 = y1 * ((double)_src->rows/_size.height);
			_dst->at<cv::Vec3b>(i)=_src->at<cv::Vec3b>(y2,x2);
		}
	}
}
