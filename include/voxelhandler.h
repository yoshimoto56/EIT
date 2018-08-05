#pragma once
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <GL/freeglut.h>
#include <pointlineface.h>
#include <vectormatrix.h>
#include <transfermatrix.h>
#include <colormap.h>
#include <glutilities.h>
#include <imageprocessing.h>

#define GAIN_MARGIN 1.2

namespace EITS
{
	#define FREE_COL Vector3d(0,0,0)
	#define SURF_COL Vector3d(0,0,1)
	#define INTR_COL Vector3d(0,1,1)
	#define INVIS_COL Vector3d(1,1,1)

	enum COORDINATE{
		COORD_LOCAL,
		COORD_WORLD
	};
	enum{
		X_AXIS,
		Y_AXIS,
		Z_AXIS
	};

	class VoxModel
	{
	private:
		class Voxel
		{
		public:
			int neibor[27];
			bool is_surface;
			bool is_free_space;
			int label;
			int label_index;
			Vector3d color;
			Vector3d center;
			Voxel &operator=(const Voxel &parent){
				label=parent.label;
				label_index=parent.label_index;
				color=parent.color;
				center=parent.center;
				is_surface=parent.is_surface;
				is_free_space=parent.is_free_space;
				memcpy(neibor, parent.neibor, sizeof(int)*27);
				return(*this);
			}
			Voxel(){
				label=-1;
				center=Vector3d(0,0,0); 
				color=Vector3d(0,0,0); 
				is_surface=false; 
				is_free_space=false;
				label_index=-1;
				memset(neibor, -1, sizeof(int) * 27);
			}
			~Voxel(){}
		}***data;
		ColorMap colormap;

		Material mat;
		bool is_loaded;
		bool is_view_data;
		bool is_draw_wire_frame;
		bool is_view_slice;
		bool is_auto_scale;
		bool is_milling;
		bool is_affine;
		int num_voxel;
		int map_depth;
		int width;
		int height;
		int depth;
		int free_label;
		double scale;
		Vector3d size;
		Vector3d center;
		double length;
		unsigned char *pixels;
		unsigned char *buffer;
		IMG *slice;
		Vector3i slice_pos;
		GLuint name[3];

		//モデル-ボクセル変換行列の計算
		void calModel2Voxel();
		//スケール位置自動調整行列の計算
		void calAutoAffine();
		//ボクセル中心の計算
		void calVoxelCenter();
		//近傍ボクセルの設定
		void calNeiborVoxel();

	public:
		VoxModel(void);
		~VoxModel(void);

		//座標系の行列
		localObject Tr;

		void initTexture();
		//スケールは最大サイズ
		void setScale(double _scale){scale=_scale * GAIN_MARGIN;}
		double getScale(){return scale;}
		//サイズは3次元
		void setSize(Vector3d _size){size=_size;}
		void setSize(double _size){size=Vector3d(_size,_size,_size);}
		void setVoxelNum(int _width, int _height, int _depth);
		void setVoxelSameNum(int _sizeVoxel);
		void setVoxelValue(int x, int y, int z, double _val){
			data[x][y][z].color=Vector3d(0,0,1);
			data[x][y][z].is_surface=true;}
		void setVoxelColor(int x, int y, int z, Vector3d _color){
			data[x][y][z].color=_color;
			data[x][y][z].is_surface=true;}
		void setVoxelValue(int _index, double _val){
			int x,y,z;
			z=_index%depth;
			y=(int)(_index/depth)%height;
			x=(int)(_index/depth)/height;
			data[x][y][z].color=Vector3d(0,0,1);
			data[x][y][z].is_surface=true;}

		Vector3d getSize(){return size;}
		int getWidth(){return width;}
		int getHeight(){return height;}
		int getDepth(){return depth;}
		int getNumVoxel(){return num_voxel;}
		Vector3d getCenter(int x, int y, int z){return data[x][y][z].center;}
		Vector3d getCenter(int _index){
			int x,y,z;
			z=_index%depth;
			y=(int)(_index/depth)%height;
			x=(int)(_index/depth)/height;
			return data[x][y][z].center;
		}

		//座標設定の計算
		void calCoordinateSetting();

		//変換行列の設定及び取得

		int getIndexOfVoxelAt(Vector3d _Puni, int _coordinate=COORD_WORLD);
		void deleteVoxel();
		void newVoxel();
		void render(bool _is_scaled=true);
		void renderGuid();
		void clear();
		bool getIsLoaded(){return this->is_loaded;}
		bool load(const char* _filename);
		bool save(const char* _filename);
		void setCenter(Vector3d _center){this->center=_center;}
		Vector3d getCenter(){return this->center;}
		void setLength(double _length){this->length=_length;}
		double getLength(){return length;}
		void setIsDrawWireFrame(bool _is_draw_wire_frame){this->is_draw_wire_frame=_is_draw_wire_frame;}
		void setIsLoaded(bool _is_loaded){this->is_loaded = _is_loaded;}
		//描画モードの設定
		void setIsViewData(bool _is_view_data){this->is_view_data=_is_view_data;}
		void setIsViewSlice(bool _is_view_slice){this->is_view_slice=_is_view_slice;}
		void setIsAutoScale(bool _is_auto_scale){this->is_auto_scale=_is_auto_scale;}
		bool getIsAutoScale(){return this->is_auto_scale;}

		Vector3d getColorAt(int _width, int _height ,int _depth);
		Vector3d getColorAt(int _index);
		int getLabelIndexAt(int _index);
		void setFreeSpaceAt(int _index);
		void setSlicePos(Vector3i _slice_pos){this->slice_pos=_slice_pos;}
		void setSliceSize(int _size);
		void setSlice(bool _invert = false);
		void renderGrid();
		void renderSlicer();
		void renderImages(int _width, int _height);

		bool saveSlice(const char *_filename, int _axis = X_AXIS);

		//クロージングをスライス画像に対して適用する
		void closingSlice(int _times);

		//特定の色のボクセル近傍を，特定の色の条件であれば指定した色に変える(戻り値は変えたボクセルの個数)
		int setNeiborColorWith(Vector3d _tColor, Vector3d _nColor, Vector3d _cColor=Vector3d(0,0,0));
		//点で挟む領域（注意：_point1は物体の内部に設定しなければならない）
		void setLabelBetween(Vector3d _point1, Vector3d _point2, Vector3d _color, int _label);
		//ラベルをセットする
		int setLabel();
		//連結空間を特定の色で塗りつぶす
		void setColorMap(int _num);
		//特定のボクセルの色と同じボクセルをすべて指定した色で塗り替える
		void resetColorAt(int _index=0, Vector3d _color=Vector3d(0,0,0), bool _is_free_space=true);
		//特定のボクセルの色以外のボクセルをSurfaceと内部にわけて塗りかえる
		void resetColorSurfAndInExceptAt(int _index, Vector3d _sColor, Vector3d _iColor);
		//物体の内部にグラデーションマップを生成する
		void setGradMapWithinSurf();
		//物体の内部を表面と同じ色で塗りつぶす
		void setColorWithinSurf();

		//特定のindexにおけるSurfaceボクセルを削り、Surfaceを設定し直す
		bool millingSurfaceWith(int _index);
		bool coloringSurfaceWith(int _index, Vector3d _col);
		//これは単に特定の位置におけるボクセルを自由空間にするMethod
//		void millingVoxels(Vector3d _Pworld);

	};
};
