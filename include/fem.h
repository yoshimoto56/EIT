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
		//�d��
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

		//���d��
		double conductivity;
		//�d��
		VectorNd potential;
		//�Ǐ��d�ʂ̎擾
		double getLocalPotentialAt(Vector3d);
		//�d��
		VectorNd electric_field;
		//�d��̌v�Z
		VectorNd calElectricField();
		//�d�����x
		VectorNd current_density;
		//�d�����x�̌v�Z
		VectorNd calCurrentDensity();
		//�d�����x�̑傫��
		double abs_current_density;
		//�d�����x�̑傫���̌v�Z
		double calAbsCurrentDensity();
		//�\���s��Ɋւ���ϐ�
		Matrixd K;
		Matrixd dN;

		void clear();
	};

	class eTetrahedra:public Tetrahedra{
	public:
		eTetrahedra();
		~eTetrahedra();

		//���d��
		double conductivity;
		//�d��
		VectorNd potential;
		//�Ǐ��d�ʂ̎擾
		double getLocalPotentialAt(Vector3d);
		//�d��
		VectorNd electric_field;
		//�d��̌v�Z
		VectorNd calElectricField();
		//�d�����x
		VectorNd current_density;
		//�d�����x�̌v�Z
		VectorNd calCurrentDensity();
		//�d�����x�̑傫��
		double abs_current_density;
		//�d�����x�̑傫���̌v�Z
		double calAbsCurrentDensity();

		//�\���s��Ɋւ���ϐ�
		Matrixd K;
		Matrixd dN;

		void clear();

	};

	class FEME2D:public SurfMesh{
	private:

		//�����s�񂪗��p�\��
		bool is_matrix_available;
		//�J���[�}�b�v
		ColorMap color_map;
		//��̓��[�h(EIT or EMBIT)
		int mode;

		//�\���s��
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

		//���E�̓_��
		int num_grounded;
		//�ڐG�_�̐�
		int num_electrode;
		//�d�ɂ̐�
		int num_sensing;
		//�S�d�ɂ��\������m�[�h�̐�
		int num_sensing_node;
		//���̐�
		int num_candidate;

		eNode *node;
		eFacet *facet;

		//�e�m�[�h�d�ʂɊւ���ϐ�
		VectorNd potential;
		VectorNd potential_A;
		VectorNd potential_B;
		//�e�m�[�h�d���Ɋւ���ϐ�
		VectorNd current;
		VectorNd current_A;
		VectorNd current_B;

		//��͗p�̓d�ɂ̏��
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
		//�ڐG�_
		Vector3d position_touch;
		//�v������������
		bool is_sensing_ready;
		//�\�����[�h
		int mode_view;

		//EIT�p�̕ϐ�
		//EIT�̑S�����̓d�Ƀm�[�h�ɂ�����d��
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


		//�����s��̌v�Z
		void calStiffnessMatrix();
		//�t�s��̎��O�v�Z
		void calInverseMatrix();
		//Set mode of the node at
		void setNodeMode(Vector3d _coord);
		//����̃}�e���A���ԍ��̃m�[�h���Œ�_�ɂ���
		void setNodeFixedAtMaterial(int _mat);
		//Set node selection mode as 
		void setNodeSelectionMode(int _mode);
		void setNodeFlagFunc(int);

		/*�e��ݒ�t�@�C���̓Ǎ��o��*/
		//���������̐ݒ�Ɋւ���֐�
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

		/*���E�����Ɋւ���֐�*/
		//����̃��x���̓d�ɂ��O�����h�ɂ���
		void selectStateAtLabelOf(int _state, int _label);
		//label��GND���Clabel+1��Electrode���Z�b�g����
		void setGroundAndElectrode(int _label);
		//label��GND���Celem�̎��͂�Electrode���Z�b�g����
		void setGroundAndElectrodeAt(int _label, int _elem);
		//label��GND���C�I�����ꂽ�v�f�̃m�[�h��Electrode���Z�b�g����
		void setGroundAndElectrodeAtSelected(int _label);
		//list�Ŏw�肵���m�[�h�ɓd�ʂ�^����
		void setPotential(double *_potential, std::vector<int> _list);
		//Electrode���x���ɓd�ʂ�^����
		void setPotential(double _potential);

		/*�\���o�Ɋւ���֐�*/
		//�I�����ꂽ�ʂ̓��d����ύX����
		void setConductivityForSelected(double _conductivity);
		//�ǂݍ��񂾓��d���ω��𕨐��l�ɔ��f������
		void updateConductivityWithDSigma();
		//�f�B���N����������d�ʂ����߂�
		void Potential2Potential();
		//�d�����x���v�Z����
		void calCurrentDensity();
		//�v�Z���ʂ��X�V����
		void update();
		//Ku=j�𖢒m�p�����[�^���m�p�����[�^�ɕ������ĉ���
		void solve();
		//����̏����ɑ΂���������߂�
		void solveFor(int _condition, double _voltage = 1);
		//���E��������@�̃T���v���f�[�^�����p
		void solveBoundFor(int _condition, double _voltage = 1);

		double getPotentialAt(Vector3d _position);
		double getSigmaAt(Vector3d _position);
		void selectObject(Vector3d *_selected_coord, int _mode);
		void selectObject2(Vector3d *_selected_coord, int _mode);
		void setLabelAndStateToSelected(int _label, int _state);
		void clearLabel(int _state = NONE);

		std::vector<double> getLineProfile(int _resolution, Vector3d _v1, Vector3d _v2);

		/*EIT�p�̊֐�*/
		void setMode(int _mode);
		int getMode();
		/*�t��萄��@*/
		void calJacobianFor(int _condition);
		//�d�ɏ���condition�ɂ��āC�v�felem�̃��R�r�A�������߂�
		void calJacobianForAt(int _condition, int _elem);
		//�d�ɏ���condition�ɂ��āC�v�felem�̋��E���R�r�A�������߂�
		void calBJacobianForAt(int _condition, int _elem);
		//�d�ɂ��Ƃɕ��ω�����
		void calMeanJacobian();
		//�d�ɂ��Ƃɓd�ʂ̕��ϒl�����߂�
		void calMeanPotential();
		//�ۓ���^���č����s����v�Z����
		void setConductivityPurtabationAt(double _d_sigma, int _elem);
		//�t�����������߂̍s������߂Ă���
		void calAJacobian(double _lambda = 0, int _method = NEWTON_TIKHONOV);
		//�t��������
		void solveInverseProblem(bool _is_set_cmap = true);
		//||Ax-y||
		double calResidualVector();
		//||Rx||
		double calRegularizedSolution();
		//�ڐG���S
		Vector3d calCenterOfContact(int _mode = MEASURED, double _ratio = 0.5);
		//�v���f�[�^�Ɛڑ����邽�߂̊֐�
		void setMeasuredPotential(int _condition, double *_voltage);
		//�ő�l���X�V
		void updateMinMax(bool _is_set_cmap);

		/*���E��������@*/
		//�^�����d�ʃp�^�[����������Ƃ��炵�����͏�����������
		int getOptimizedInputFor(int _condition, double *_voltage);
		void calErrorBtwSimMeas();
		void calMinMaxPotentialSim();
		void calMinMaxPotentialMeas();

		//�ڐG�_���W��`����W�n�ɕύX
		Vector3d getPositionTouch(){
			if(this->is_auto_scale){
				double scale = 10.0 / this->size.abs();
				return scale * (this->position_touch + center);
			}
			return this->position_touch;
		}
		bool getIsSensingReady(){return this->is_sensing_ready;}

		//�]���Ɋւ���֐�(See Tawil et al., 2011)
		//�ʐϔ�
		double calSpatialResolution(double _rate = 0.5);
		//Silvera2012_PhD_eq(4.2)
		double calShapeDeformation();
		//�ʐό덷
		double calSizeError();
		//�ʒu�덷
		double calPositionError();
		//TODO:�d�ʂ��Čv�Z
		double calPotentialError();
	};

	class FEME3D:public VolumeMesh{
	private:

		//�����s�񂪗��p�\��
		bool is_matrix_available;
		//�J���[�}�b�v
		ColorMap color_map;

		//�\���s��
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

		//���E�̓_��
		int num_grounded;
		int num_electrode;

		eNode *node;
		eTetrahedra *elem;

		//�e�m�[�h�d�ʂɊւ���ϐ�
		VectorNd potential;
		VectorNd potential_A;
		VectorNd potential_B;
		//�e�m�[�h�d���Ɋւ���ϐ�
		VectorNd current;
		VectorNd current_A;
		VectorNd current_B;

		//��͗p�̓d�ɂ̏��
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

		//���������̐ݒ�Ɋւ���֐�
		void loadParameter(const char* _filename);
		//�����s��̌v�Z
		void calStiffnessMatrix();
		//�t�s��̎��O�v�Z
		void calInverseMatrix();
		//Set mode of the node at
		void setNodeMode(Vector3d _coord);
		//����̃}�e���A���ԍ��̃m�[�h���Œ�_�ɂ���
		void setNodeFixedAtMaterial(int _mat);
		//Set node selection mode as 
		void setNodeSelectionMode(int _mode);
		void setNodeFlagFunc(int);

		//�e��ݒ�t�@�C���̓Ǎ��o��
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


		//�d���Z�b�g�̓Ǎ�
		bool loadMeasuredVoltageSet(const char* _filename);

		//�d�ʉ�͂Ɋւ���Solver
		void update();
		void setPotential(double *_potential, std::vector<int> _list);
		void setPotential(double _potential);
		//�f�B���N����������d�ʂ����߂�
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
