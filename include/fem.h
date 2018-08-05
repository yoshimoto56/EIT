#pragma once
#define MODEL_SCALE 20.0
#define TARGET_CONDUCTIVITY 100

#include <cpplapack.h>//To solve Ax=y
#include <GL/freeglut.h>
#include <eitsheads.h>

namespace EITS{
	enum REGULARISATION
	{
		GAUSSIAN,
		LAPLACIAN,
		NEWTON_TIKHONOV
	};
	enum VIEWMODE
	{
		POTENTIAL,
		CONDUCTIVITY
	};
	enum ESTMETHOD
	{
		EIT,
		EMBIT
	};
	enum EMODE
	{
		MEASURED,
		TARGET
	};

	class eNode: public Node{
	public:
		eNode(){this->clear();}
		~eNode(){}
		//電位
		double potential;
		std::vector<int> index_list_facet;
		void clear(){
			this->potential = 0;
			//this->index = -1;
			this->state = 0;
			//this->vertex = Vector3d(0,0,0);
			//this->normal = Vector3d(0,0,0);
		}
	};

	class eFacet:public Facet{
	public:
		eFacet();
		~eFacet();

		Facet &operator = (const eFacet &_facet);

		//導電率
		double conductivity;
		//電位
		VectorNd potential;
		//局所電位の取得
		double getLocalPotentialAt(Vector3d);
		//電場
		VectorNd electric_field;
		//電場の計算
		VectorNd calElectricField();
		//電流密度
		VectorNd current_density;
		//電流密度の計算
		VectorNd calCurrentDensity();
		//電流密度の大きさ
		double abs_current_density;
		//電流密度の大きさの計算
		double calAbsCurrentDensity();
		//構成行列に関する変数
		Matrixd K;
		Matrixd dN;

		void clear();
	};

	class eTetrahedra:public Tetrahedra{
	public:
		eTetrahedra();
		~eTetrahedra();

		//導電率
		double conductivity;
		//電位
		VectorNd potential;
		//局所電位の取得
		double getLocalPotentialAt(Vector3d);
		//電場
		VectorNd electric_field;
		//電場の計算
		VectorNd calElectricField();
		//電流密度
		VectorNd current_density;
		//電流密度の計算
		VectorNd calCurrentDensity();
		//電流密度の大きさ
		double abs_current_density;
		//電流密度の大きさの計算
		double calAbsCurrentDensity();

		//構成行列に関する変数
		Matrixd K;
		Matrixd dN;

		void clear();

	};

	class FEME2D:public SurfMesh{
	private:

		//剛性行列が利用可能か
		bool is_matrix_available;
		//カラーマップ
		ColorMap color_map;
		//解析モード(EIT or EMBIT)
		int mode;

		//構成行列
		Matrixd K;
		Matrixd K_A;
		Matrixd K_AA;
		Matrixd K_AB;
		Matrixd K_BA;
		Matrixd K_BB;
		Matrixd invK_A;
		Matrixd invK_AB;
		Matrixd invK_BA;
		Matrixd invK_AA;
		Matrixd invK_BB;

		std::vector<int> renum_A;
		std::vector<int> renum_AA;

		//境界の点数
		int num_grounded;
		//接触点の数
		int num_electrode;
		//電極の数
		int num_sensing;
		//全電極を構成するノードの数
		int num_sensing_node;
		//候補の数
		int num_candidate;

		eNode *node;
		eFacet *facet;

		//各ノード電位に関する変数
		VectorNd potential;
		VectorNd potential_A;
		VectorNd potential_B;
		//各ノード電流に関する変数
		VectorNd current;
		VectorNd current_A;
		VectorNd current_B;

		//解析用の電極の情報
		int curElectIndex;
		int curInputIndex;

		std::vector<int> index_list_input;
		std::vector<int> index_list_sensing;
		std::vector<std::vector<int>> index_list_sensing_each;
		std::vector<int> index_list_free_element;

		std::vector<double> mvoltage;
		std::vector<double> meanf;
		std::vector<double> error;
		std::vector<int> mvIndexList;
		std::vector<int> enodeIndexList;

		//size = num_sensing*num_candidate*num_sensing
		std::vector<double>simulatedPotential;
		//size = num_sensing*num_candidate
		std::vector<double>min_potential_sim;
		std::vector<double>max_potential_sim;
		//size = num_sensing*num_sensing
		std::vector<double>measuredPotential;
		double min_potential_meas;
		double max_potential_meas;
		//接触点
		Vector3d position_touch;
		//計測準備完了か
		bool is_sensing_ready;
		//表示モード
		int mode_view;

		//EIT用の変数
		//EITの全条件の電極ノードにおける電位
		std::vector<double> potential_sim;
		std::vector<double> potentialm_sim;
		std::vector<double> potential_meas;
		std::vector<double> potentialm_meas;
		std::vector<double> V1_c;
		std::vector<double> V2_c;
		std::vector<double> J_cm_free;
		std::vector<double> J_cm;
		std::vector<double> J_c;
		bool is_jacobian_loaded;
		double lambda;
/*
		CPPL::dgematrix A_J;
		CPPL::dgematrix J;
		CPPL::dgematrix tJ;
		CPPL::dgematrix Q;
		CPPL::dcovector dV;
		CPPL::dcovector dSigma;
*/
		Matrixd J;
		Matrixd tJ;
		Matrixd A_J;
		Matrixd invA_J;
		Matrixd invA_JtJ;
		Matrixd Q;
		VectorNd dV;
		VectorNd dSigma;
		Vector2d mSigma;

	public:
		FEME2D();
		~FEME2D();
		void newMesh();
		void deleteMesh();
		void clear();
		void render();
		bool load(const char *_filename);
		void init();
		void calCenter();
		void calLine();
		void calVertex();
		void calFacetList();

		int getNumElectrode(){return this->num_electrode;}


		//剛性行列の計算
		void calStiffnessMatrix();
		//逆行列の事前計算
		void calInverseMatrix();
		//Set mode of the node at
		void setNodeMode(Vector3d _coord);
		//特定のマテリアル番号のノードを固定点にする
		void setNodeFixedAtMaterial(int _mat);
		//Set node selection mode as 
		void setNodeSelectionMode(int _mode);
		void setNodeFlagFunc(int);

		/*各種設定ファイルの読込出力*/
		//物理特性の設定に関する関数
		void loadConductivity(const char* _filename);
		void saveConductivity(const char* _filename);
		bool loadPreCalMatrix(const char* _filename);
		bool savePreCalMatrix(const char* _filename);
		bool loadElectrodeSetting(const char* _filename);
		bool saveElectrodeSetting(const char* _filename);
		bool loadElectrodeIndexList(const char* _filename);
		bool savePotentialList(const char* _filename);
		bool saveSimulatedPotentialForAll(const char* _filename);
		bool loadSimulatedPotentialForAll(const char* _filename);
		bool savePotentialSimulatedFor(int _condition, const char* _filename);
		bool saveMeanPotentialSimulatedFor(int _condition, const char* _filename);
		bool loadMeanPotentialSimulatedFor(int _condition, const char* _filename);
		bool loadPotentialSimulatedFor(int _condition, const char* _filename);
		bool loadPotentialMeasuredFor(int _condition, const char* _filename);
		bool loadMeanPotentialMeasuredFor(int _condition, const char* _filename);
		bool saveJacobianFor(const char* _filename);
		bool saveMeanJacobianFor(const char* _filename);
		bool loadJacobian(const char* _filename);
		bool saveInvJacobian(const char* _filename);
		bool loadInvJacobian(const char* _filename);
		bool loadMeanJacobian(const char* _filename);
		bool saveDSigma(const char* _filename);
		bool loadDSigma(const char* _filename);

		void clearSimulatedPotential();

		bool getIsMatrixAvailable(){return this->is_matrix_available;}	
		void setIsMatrixAvailable(bool _is_matrix_available){this->is_matrix_available = _is_matrix_available;}

		ColorMap setColorMap(ColorMap _color_map){color_map = _color_map;}
		ColorMap getColorMap(){return color_map;}
		void setModeView(int _mode){this->mode_view = _mode;
			if(_mode == POTENTIAL){
				color_map.setParam(0, 3);
				std::cout<<"Vis-Mode: POTENTIAL"<<std::endl;
			}
			else if(_mode == CONDUCTIVITY){
				color_map.setParam(0, 150);
				std::cout<<"Vis-Mode: CONDUCTIVITY"<<std::endl;
			}
		};

		/*境界条件に関する関数*/
		//特定のラベルの電極をグランドにする
		void selectStateAtLabelOf(int _state, int _label);
		//labelにGNDを，label+1にElectrodeをセットする
		void setGroundAndElectrode(int _label);
		//labelにGNDを，elemの周囲にElectrodeをセットする
		void setGroundAndElectrodeAt(int _label, int _elem);
		//labelにGNDを，選択された要素のノードにElectrodeをセットする
		void setGroundAndElectrodeAtSelected(int _label);
		//listで指定したノードに電位を与える
		void setPotential(double *_potential, std::vector<int> _list);
		//Electrodeラベルに電位を与える
		void setPotential(double _potential);

		/*ソルバに関する関数*/
		//選択された面の導電率を変更する
		void setConductivityForSelected(double _conductivity);
		//読み込んだ導電率変化を物性値に反映させる
		void updateConductivityWithDSigma();
		//ディリクレ条件から電位を求める
		void Potential2Potential();
		//電流密度を計算する
		void calCurrentDensity();
		//計算結果を更新する
		void update();
		//Ku=jを未知パラメータ既知パラメータに分解して解く
		void solve();
		//特定の条件に対する解を求める
		void solveFor(int _condition, double _voltage = 1);
		//境界条件推定法のサンプルデータ生成用
		void solveBoundFor(int _condition, double _voltage = 1);

		double getPotentialAt(Vector3d _position);
		double getSigmaAt(Vector3d _position);
		void selectObject(Vector3d *_selected_coord, int _mode);
		void selectObject2(Vector3d *_selected_coord, int _mode);
		void setLabelAndStateToSelected(int _label, int _state);
		void clearLabel(int _state = NONE);

		std::vector<double> getLineProfile(int _resolution, Vector3d _v1, Vector3d _v2);

		/*EIT用の関数*/
		void setMode(int _mode);
		int getMode();
		/*逆問題推定法*/
		void calJacobianFor(int _condition);
		//電極条件conditionについて，要素elemのヤコビアンを求める
		void calJacobianForAt(int _condition, int _elem);
		//電極条件conditionについて，要素elemの境界ヤコビアンを求める
		void calBJacobianForAt(int _condition, int _elem);
		//電極ごとに平均化する
		void calMeanJacobian();
		//電極ごとに電位の平均値を求める
		void calMeanPotential();
		//摂動を与えて剛性行列を計算する
		void setConductivityPurtabationAt(double _d_sigma, int _elem);
		//逆問題を解くための行列を求めておく
		void calAJacobian(double _lambda = 0, int _method = NEWTON_TIKHONOV);
		//逆問題を解く
		void solveInverseProblem(bool _is_set_cmap = true);
		//||Ax-y||
		double calResidualVector();
		//||Rx||
		double calRegularizedSolution();
		//接触中心
		Vector3d calCenterOfContact(int _mode = MEASURED, double _ratio = 0.5);
		//計測データと接続するための関数
		void setMeasuredPotential(int _condition, double *_voltage);
		//最大値を更新
		void updateMinMax(bool _is_set_cmap);

		/*境界条件推定法*/
		//与えた電位パターンからもっともらしい入力条件を見つける
		int getOptimizedInputFor(int _condition, double *_voltage);
		void calErrorBtwSimMeas();
		void calMinMaxPotentialSim();
		void calMinMaxPotentialMeas();

		//接触点座標を描画座標系に変更
		Vector3d getPositionTouch(){
			if(this->is_auto_scale){
				double scale = 10.0 / this->size.abs();
				return scale * (this->position_touch + center);
			}
			return this->position_touch;
		}
		bool getIsSensingReady(){return this->is_sensing_ready;}

		//評価に関する関数(See Tawil et al., 2011)
		//面積比
		double calSpatialResolution(double _rate = 0.5);
		//Silvera2012_PhD_eq(4.2)
		double calShapeDeformation();
		//面積誤差
		double calSizeError();
		//位置誤差
		double calPositionError();
		//TODO:電位を再計算
		double calPotentialError();
	};

	class FEME3D:public VolumeMesh{
	private:

		//剛性行列が利用可能か
		bool is_matrix_available;
		//カラーマップ
		ColorMap color_map;

		//構成行列
		Matrixd K;
		Matrixd K_A;
		Matrixd K_AA;
		Matrixd K_AB;
		Matrixd K_BA;
		Matrixd invK_A;
		Matrixd invK_AA;
		Matrixd invK_AB;
		Matrixd invK_BA;

		std::vector<int> renum_A;
		std::vector<int> renum_AA;

		//境界の点数
		int num_grounded;
		int num_electrode;

		eNode *node;
		eTetrahedra *elem;

		//各ノード電位に関する変数
		VectorNd potential;
		VectorNd potential_A;
		VectorNd potential_B;
		//各ノード電流に関する変数
		VectorNd current;
		VectorNd current_A;
		VectorNd current_B;

		//解析用の電極の情報
		int curElectIndex;
		int curInputIndex;

		std::vector<int> index_list_input;
		std::vector<int> index_list_electrode;

		std::vector<double> mvoltage;
		std::vector<double> meanf;
		std::vector<double> error;
		std::vector<int> mvIndexList;
		std::vector<int> enodeIndexList;

		//size=numClip*tNumCandidate*numClip
		std::vector<double>simulatedPotential;
	public:
		FEME3D();
		~FEME3D();
		void newMesh();
		void deleteMesh();
		void clear();
		void render();
		bool load(const char *_filename);
		void init();
		void calCenter();
		void calLine();
		void calFacet();

		//物理特性の設定に関する関数
		void loadParameter(const char* _filename);
		//剛性行列の計算
		void calStiffnessMatrix();
		//逆行列の事前計算
		void calInverseMatrix();
		//Set mode of the node at
		void setNodeMode(Vector3d _coord);
		//特定のマテリアル番号のノードを固定点にする
		void setNodeFixedAtMaterial(int _mat);
		//Set node selection mode as 
		void setNodeSelectionMode(int _mode);
		void setNodeFlagFunc(int);

		//各種設定ファイルの読込出力
		bool loadPreCalMatrix(const char* _filename);
		bool savePreCalMatrix(const char* _filename);
		bool loadElectrodeSetting(const char* _filename);
		bool saveElectrodeSetting(const char* _filename);
		bool loadElectrodeIndexList(const char* _filename);
		bool savePotentialList(const char* _filename);
		bool saveSimulatedPotentialForAll(const char* _filename);
		bool loadSimulatedPotentialForAll(const char* _filename);

		bool getIsMatrixAvailable(){return this->is_matrix_available;}	
		void setIsMatrixAvailable(bool _is_matrix_available){this->is_matrix_available = _is_matrix_available;}

		ColorMap setColorMap(ColorMap _color_map){color_map = _color_map;}
		ColorMap getColorMap(){return color_map;}


		//電圧セットの読込
		bool loadMeasuredVoltageSet(const char* _filename);

		//電位解析に関するSolver
		void update();
		void setPotential(double *_potential, std::vector<int> _list);
		void setPotential(double _potential);
		//ディリクレ条件から電位を求める
		void Potential2Potential();
		double getPotentialAt(Vector3d _position);

		void calCurrentDensity();

		void selectObject(Vector3d *_selected_coord, int _mode);

		/*
		void setGroundElectrodeList(bool _isGroundOnly=false);
		void selectGroundElectrodesAt(int _index);
		void selectActiveElectrodesAt(int _index);
		bool nextInputAs(double _voltage);
		void evaluateInput();
		double getErrorAt(int _condition, int _index, int _numVoltage, double *_voltage); 
		int getOptimizedInputFor(int _condition, int _numVoltage, double *_voltage);
		void getOptimizedInput(double _initialValue);
		double getMaximumSimulatedPotentialFor(int _condition);
		double getMinimumSimulatedPotentialFor(int _condition);
		*/
	};
};
