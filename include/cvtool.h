#pragma once
#include "cvheads.h"

namespace CVTOOL
{
	static cv::Mat_<double> G1=(cv::Mat_<double>(3,3)<<0,0,1,0,0,0,0,0,0);//���[�㐔�̊��s��
	static cv::Mat_<double> G2=(cv::Mat_<double>(3,3)<<0,0,0,0,0,1,0,0,0);//���[�㐔�̊��s��
	static cv::Mat_<double> G3=(cv::Mat_<double>(3,3)<<0,1,0,0,0,0,0,0,0);//���[�㐔�̊��s��
	static cv::Mat_<double> G4=(cv::Mat_<double>(3,3)<<0,0,0,1,0,0,0,0,0);//���[�㐔�̊��s��
	static cv::Mat_<double> G5=(cv::Mat_<double>(3,3)<<1,0,0,0,-1,0,0,0,0);//���[�㐔�̊��s��
	static cv::Mat_<double> G6=(cv::Mat_<double>(3,3)<<0,0,0,0,-1,0,0,0,1);//���[�㐔�̊��s��
	static cv::Mat_<double> G7=(cv::Mat_<double>(3,3)<<0,0,0,0,0,0,1,0,0);//���[�㐔�̊��s��
	static cv::Mat_<double> G8=(cv::Mat_<double>(3,3)<<0,0,0,0,0,0,0,1,0);//���[�㐔�̊��s��
	static cv::Mat_<double> J_G=(cv::Mat_<double>(9,8) <<
						G1(0),G2(0),G3(0),G4(0),G5(0),G6(0),G7(0),G8(0),
						G1(1),G2(1),G3(1),G4(1),G5(1),G6(1),G7(1),G8(1),
						G1(2),G2(2),G3(2),G4(2),G5(2),G6(2),G7(2),G8(2),
						G1(3),G2(3),G3(3),G4(3),G5(3),G6(3),G7(3),G8(3),
						G1(4),G2(4),G3(4),G4(4),G5(4),G6(4),G7(4),G8(4),
						G1(5),G2(5),G3(5),G4(5),G5(5),G6(5),G7(5),G8(5),
						G1(6),G2(6),G3(6),G4(6),G5(6),G6(6),G7(6),G8(6),
						G1(7),G2(7),G3(7),G4(7),G5(7),G6(7),G7(7),G8(7),
						G1(8),G2(8),G3(8),G4(8),G5(8),G6(8),G7(8),G8(8));//���R�r�A��

	enum G_TYPE{HOMOGRAPHY,};//�ϊ��s��̃^�C�v
	enum DERIVATION_TYPE{SOBEL,};//�摜�̌��z�̌v�Z���@
	/*=== �֐� ===*/
	//�ϊ�_G�ɂ��__p�̕ϊ�
	cv::Point warpPoint(cv::Point _p, cv::Mat _G);
	//�摜_img�̓__p�ɂ�������z���v�Z
	cv::Point derivation(cv::Mat _img, cv::Point _p, DERIVATION_TYPE _TYPE=SOBEL);
	//�摜_img�̓_warpPoint(_p,_G)�ɂ�������z���v�Z
	cv::Point derivation(cv::Mat _img, cv::Point _p, cv::Mat _H, DERIVATION_TYPE _TYPE=SOBEL);
	//�W�����̃p�����[�^�x�N�g��_x����C�z���O���t�B�ϊ��s����v�Z
	cv::Mat getHomography(cv::Mat _x, double _n = 170);//�w���v�Z��DBL_MAX�̃I�[�_�Ɠ����ɂȂ�悤�f�t�H���g�J��Ԃ��񐔂��w�肵�Ă���
	//�T�C�Y�ƃz���O���t�B�s���^���C�摜���g���~���O����֐�
	cv::Mat trimImageWithHomography(cv::Size _size, cv::Mat _src, cv::Mat _H);

	void drawHistgram(cv::Mat *_img_src, std::vector<double> _hist);
	std::vector<cv::Mat> getHSV(cv::Mat _img_src);
	cv::Mat getNormalizedImage(cv::Mat _img_src, int _max = 255);
	std::vector<double> getHistgram(cv::Mat _img_src);
	cv::Point2f getMoment(cv::Mat _img_src, cv::Point2i _dim);


	//�������T�C�Y�֐�
	void fastResize(cv::Mat *_src, cv::Mat *_dst, cv::Size _size);
};

