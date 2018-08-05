#include "fem.h"

namespace EITS{
	/****************eFacet******************/
	eFacet::eFacet()
	{
		this->setFacetTypeAsTriangle();
		this->K = Matrixd(3,3);
		this->dN = Matrixd(2,3);
		this->conductivity = 1;
		this->potential = VectorNd(3);
		this->current_density = VectorNd(2);
		this->electric_field = VectorNd(2);
		this->clear();
	}

	eFacet::~eFacet()
	{
		this->clear();
		this->K.free();
		this->dN.free();
		this->potential.free();
		this->current_density.free();
		this->electric_field.free();
	}

	void eFacet::clear()
	{
		this->K.identity();
		memset(this->current_density.X, 0, sizeof(double) * this->current_density.n);
		memset(this->electric_field.X, 0, sizeof(double) * this->electric_field.n);
		memset(this->potential.X, 0, sizeof(double) * this->potential.n);
		this->abs_current_density = 0;

		this->index_material = 0;
		this->area = 0;

		for(int i=0;i<6;i++){
			this->normal[i]=Vector3d(0,0,0);
		}
	}

	VectorNd eFacet::calElectricField()
	{
		Matrixd temp_dN(3, 2);
		temp_dN = dN.trn();
		electric_field = temp_dN * potential;
		return electric_field;
	}

	VectorNd eFacet::calCurrentDensity()
	{
		calElectricField();
		current_density = conductivity * electric_field;
		abs_current_density = current_density.abs();
		return current_density; 
	}

	Facet& eFacet::operator= (const eFacet &_facet){
		type=_facet.type;
		this->is_selected=_facet.is_selected;
		this->is_enabled=_facet.is_enabled;
		index_facet=_facet.index_facet;
		index_material=_facet.index_material;
		num_node=_facet.num_node;
		num_line=_facet.num_line;
		num_normal=_facet.num_normal;
		normal_type = _facet.normal_type;
		index_elem = _facet.index_elem;
		area = _facet.area;
		memcpy(index_node,_facet.index_node,sizeof(int)*num_node);
		memcpy(vertex,_facet.vertex,sizeof(Vector3d)*num_node);
		memcpy(index_normal,_facet.index_normal,sizeof(int)*num_normal);
		memcpy(normal,_facet.normal,sizeof(Vector3d)*num_normal);
		memcpy(line,_facet.line,sizeof(Line)*num_line);
		memcpy(K.X, _facet.K.X, sizeof(double) * _facet.K.n * _facet.K.m);
		memcpy(dN.X, _facet.dN.X, sizeof(double) * _facet.dN.n * _facet.dN.m);
		memcpy(potential.X, _facet.potential.X, sizeof(double) * _facet.potential.n);
		memcpy(electric_field.X, _facet.electric_field.X, sizeof(double) * _facet.electric_field.n);
		memcpy(current_density.X, _facet.current_density.X, sizeof(double) * _facet.current_density.n);

		conductivity = _facet.conductivity;
		return (*this);
	}
	double eFacet::getLocalPotentialAt(Vector3d _position)
	{
		double value;
		Matrixd temp_inv = Matrixd(4,4);
		VectorNd temp_result = VectorNd(4);
		VectorNd temp_position = VectorNd(4);
		temp_position.X[0] = 1;
		temp_position.X[1] = _position.x;
		temp_position.X[2] = _position.y;
		temp_position.X[3] = _position.z;

		//TODO
//		temp_inv = calVolume();
		temp_inv = calArea();
		temp_result = temp_inv * potential;
		temp_result = temp_position * temp_result;
		value = temp_result.X[0]; 
		temp_inv.free();
		temp_result.free();
		temp_position.free();
		return value;
	}
	/****************eFacet******************/

	/****************Tetrahedra******************/
	eTetrahedra::eTetrahedra()
	{
		this->K = Matrixd(4,4);
		this->position.malloc(12);
		this->potential = VectorNd(4);
		this->dN = Matrixd(3,4);
		this->conductivity = 1;
		this->current_density = VectorNd(3);
		this->electric_field = VectorNd(3);
		this->clear();
	}

	eTetrahedra::~eTetrahedra()
	{
		this->clear();
		this->K.free();
		this->position.free();
		this->potential.free();
		this->dN.free();
		this->current_density.free();
		this->electric_field.free();
	}

	void eTetrahedra::clear()
	{
		this->K.identity();
		memset(this->position.X, 0, sizeof(double) * this->position.n);
		memset(this->current_density.X, 0, sizeof(double) * this->current_density.n);
		memset(this->electric_field.X, 0, sizeof(double) * this->electric_field.n);
		memset(this->potential.X, 0, sizeof(double) * this->potential.n);
		this->abs_current_density = 0;

		this->index_material = 0;
		this->center = Vector3d(0,0,0);
		memset(this->index_node, -1, sizeof(int) * 4);
		this->volume = 0;
		order_node[0]=2;order_node[1]=1;order_node[2]=0;order_node[3]=3;
		order_node[4]=2;order_node[5]=3;order_node[6]=1;order_node[7]=0;
		order_node[8]=0;order_node[9]=3;order_node[10]=2;order_node[11]=1;
		order_node[12]=0;order_node[13]=1;order_node[14]=3;order_node[15]=2;
		for(int i=0;i<4;i++){
			this->normal[i]=Vector3d(0,0,0);
		}
	}

	VectorNd eTetrahedra::calElectricField()
	{
		Matrixd temp_dN(4, 3);
		temp_dN = dN.trn();
		electric_field = temp_dN * potential;
		return electric_field;
	}

	VectorNd eTetrahedra::calCurrentDensity()
	{
		calElectricField();
		current_density = conductivity * electric_field;
		abs_current_density = current_density.abs();
		return current_density; 
	}

	double eTetrahedra::getLocalPotentialAt(Vector3d _position)
	{
		double value;
		Matrixd temp_inv = Matrixd(4,4);
		VectorNd temp_result = VectorNd(4);
		VectorNd temp_position = VectorNd(4);
		temp_position.X[0] = 1;
		temp_position.X[1] = _position.x;
		temp_position.X[2] = _position.y;
		temp_position.X[3] = _position.z;

		temp_inv = calVolume();
		temp_result = temp_inv * potential;
		temp_result = temp_position * temp_result;
		value = temp_result.X[0]; 
		temp_inv.free();
		temp_result.free();
		temp_position.free();
		return value;
	}
	/****************Tetrahedra******************/

	/****************FEM2D******************/
	FEME2D::FEME2D()
	{
		this->is_sensing_ready = false;
		this->num_node = 0;
		this->num_facet = 0;
		this->num_grounded = 0;
		this->num_electrode = 0;
		this->num_sensing = 0;
		this->num_sensing_node = 0;
		this->num_candidate = 0;
		this->num_line = 0;
		this->newMesh();
		renum_A.push_back(0);
		renum_AA.push_back(0);

		this->mode_view = POTENTIAL;
		this->mode = EMBIT;
		this->is_jacobian_loaded = false;
		this->lambda = 1;

		this->clear();
		this->clearView();

		this->position_touch = Vector3d(0, 0, 0);

		color_map.setParam(0, 3);//TODO:アクセッサを用意してマップ範囲を変えられるように
	}

	FEME2D::~FEME2D()
	{
		this->deleteMesh();
	}

	void FEME2D::newMesh()
	{
		/*Memory domain for the FEM analisys*/
		this->node = new eNode[this->num_node + 2];
		this->line = new Line[this->num_line + 2];
		this->facet = new eFacet[this->num_facet + 2];
		this->normal = new Vector3d[this->num_normal+2];
		this->vertex = new Vector3d[this->num_node + 2];
/*
		this->color = new Vector3f[this->num_node + 2];
		this->labelIndex = new int[this->num_node + 2];
		this->is_selecteded = new bool[this->num_node+2];

		for(int i=0;i<num_node+2;i++){
			is_selecteded[i]=false;
			labelIndex[i]=-1;
		}
*/

//		this->init();
	}

	void FEME2D::deleteMesh()
	{
		this->is_loaded = false;
		delete []node;
		delete []line;
		delete []facet;
		delete []normal;
		delete []vertex;
/*
		delete []labelIndex;
		delete []color;
		delete []is_selecteded;
*/
		this->clear();
		if(num_node > 0){
			K.free();
			K_A.free();
			K_AA.free();
			K_AB.free();
			K_BA.free();
			K_BB.free();
			invK_AB.free();
			invK_BA.free();
			invK_A.free();
			invK_AA.free();

			potential.free();
			potential_A.free();
		}
	}

	void FEME2D::clear()
	{
		this->is_jacobian_loaded = false;
		this->is_sensing_ready = false;
		this->is_matrix_available = false;
		this->is_loaded = false;
		this->center = Vector3d(0, 0, 0);
		this->size = Vector3d(1, 1, 1);
	}

	void FEME2D::init()
	{
		this->curInputIndex = 0;
		if(num_node > 0){
			memset(this->potential.X,0,sizeof(double)*this->potential.n);
			memset(this->potential_A.X,0,sizeof(double)*this->potential_A.n);
			memset(this->potential_B.X,0,sizeof(double)*this->potential_B.n);
		}
/*
		for(int i=0;i<num_node;i++){
			this->node[i].clear();
			renum_A[i]=0;
			renum_AA[i]=0;
		}
*/
	}

	void FEME2D::calCenter()
	{
		std::cout<<"Calculating center...";
		Vector3d min(INT_MAX, INT_MAX, INT_MAX);
		Vector3d max(-INT_MAX, -INT_MAX, -INT_MAX);
		for(int i=0;i<this->num_node;i++){
			for(int k=0;k<3;k++){
				if(min.X[k]>this->node[i].vertex.X[k])
					min.X[k]=this->node[i].vertex.X[k];
				if(max.X[k]<this->node[i].vertex.X[k])
					max.X[k]=this->node[i].vertex.X[k];
			}
		}
		this->size.x = fabs(max.x-min.x);
		this->size.y = fabs(max.y-min.y);
		this->size.z = fabs(max.z-min.z);
		this->center = (max+min)/2.0;
		std::cout<<"[OK]"<<std::endl;
	}

	void FEME2D::calLine()
	{
		std::cout<<"Extracting lines...";
		if(!this->is_loaded){
			this->is_view_line_extracted = false;
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}

		Line t_line;
		delete []this->line;
		this->num_line = 0;
		this->line = new Line[this->num_facet * 3 + 2];//三角形なので3本

		for(int i=0;i<this->num_facet;i++){
			for(int j=0;j<3;j++){
				t_line.index_node[0] = facet[i].index_node[j];
				t_line.index_node[1] = facet[i].index_node[(j+1)%3];
				t_line.vertex[0] = facet[i].vertex[j];
				t_line.vertex[1] = facet[i].vertex[(j+1)%3];
				if(!isSharedLine(t_line, this->line, this->num_line)){
					this->num_line++;
					this->line[this->num_line-1] = t_line;
				}
			}
		}
		this->is_view_line_extracted = true;
		std::cout<<"[OK]"<<std::endl;
	}	

	void FEME2D::calVertex()
	{
		std::cout<<"Calculating vertices...";
		int count=1;
		this->node[0].vertex=this->facet[0].vertex[0];
		for(int i=0;i<this->num_facet;i++){
			for(int j=0;j<3;j++){
				for(int k=0;k<count;k++){
					if(this->facet[i].vertex[j]==this->node[k].vertex){
						this->facet[i].index_node[j]=k;
						k=count;
					}
					else if(k==count-1){
						this->node[count].vertex=this->facet[i].vertex[j];
						this->facet[i].index_node[j]=count;
						count++;
					}
				}
			}
		}
		this->num_node = count;
		std::cout<<"[OK]"<<std::endl;
/*
		eNode *t_node = new eNode[this->num_node];
		memcpy(t_node, node, sizeof(eNode)*this->num_node);
		delete []node;
		this->node = t_node;
*/
	}

	void FEME2D::calFacetList()
	{
		std::cout<<"Generating index facet list...";
		for(int i=0;i<this->num_node;i++){
			this->node[i].index_list_facet.clear();
		}
		for(int i=0;i<this->num_facet;i++){
			for(int j=0;j<3;j++){
				this->node[this->facet[i].index_node[j]].index_list_facet.push_back(i);
			}
		}
		std::cout<<"[OK]"<<std::endl;
	}

	bool FEME2D::load(const char *_filename)
	{
		this->is_loaded = false;
		std::cout <<"Loading .fem2d data...";
		int count=0;
		float *pf;
		int *pd;
		std::ifstream file;
		eFacet sample;
		char buf[256];
		file.open(_filename, std::ios::in|std::ios::binary);
		if (!file.is_open()){
			std::cout << "[FAIL]" << std::endl;
			return false;
		}
		file.read(buf,84);
		pd = (int *)&buf[80];
		this->deleteMesh();//TODO:２回目で落ちるので要チェック
		this->num_facet=*pd;
		this->num_normal=this->num_facet;
		this->num_node=3*this->num_facet;
		this->num_line=3*this->num_facet;
		this->newMesh();
		sample.setFacetTypeAsTriangle();
		while (file.read(buf,50)){
			pf = (float *)buf;
			sample.normal[0].x=*(pf+ 0);
			sample.normal[0].y=*(pf+ 1);
			sample.normal[0].z=*(pf+ 2);
			sample.vertex[0].x=*(pf+ 3);
			sample.vertex[0].y=*(pf+ 4);
			sample.vertex[0].z=*(pf+ 5);
			sample.vertex[1].x=*(pf+ 6);
			sample.vertex[1].y=*(pf+ 7);
			sample.vertex[1].z=*(pf+ 8);
			sample.vertex[2].x=*(pf+ 9);
			sample.vertex[2].y=*(pf+10);
			sample.vertex[2].z=*(pf+11);
			sample.index_facet=count;
			sample.index_normal[0]=count;
			sample.index_normal[0]=count;
			sample.index_normal[0]=count;
			normal[count]=sample.normal[0];
			facet[count]=sample;
			count++;
		}
		std::cout<<"[OK]"<<std::endl;
		file.close();

		this->is_loaded = true;
		this->calVertex();
		this->calCenter();
		this->calTransform();
//		this->calFacetNormal();
		this->calScale();
		this->calLine();
		this->calFacetList();
		std::cout<<"Number of node: "<<this->num_node<<std::endl;
		std::cout<<"Number of element: "<<this->num_facet<<std::endl;
		return true;
	}

	bool FEME2D::loadPreCalMatrix(const char* _filename)
	{
		this->is_matrix_available = false;
		Sleep(10);

		std::ifstream file;
		int n;
		std::cout<<"Loading electrical matrix data...";
		file.open(_filename, std::ios::in|std::ios::binary);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		file.read((char*)&n,4);
		if(this->num_node != n){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		this->num_node = n;
		file.read((char*)&this->num_grounded,4);
		K_A = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		invK_A = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		K_AA = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		K_AB = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		invK_AA = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		invK_AB = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);

		for(int i=0;i<this->num_node;i++)
			file.read((char*)&this->node[i].state,sizeof(int));
		for(int i=0;i<this->num_node-this->num_grounded;i++)
			file.read((char*)&this->renum_A[i],sizeof(int));
		file.read((char*)invK_A.X,invK_A.n*invK_A.m*sizeof(double));

		if(K.n!=0){
			K.free();
			K.n=0;
		}
		this->is_matrix_available = true;
		std::cout<<"[OK]"<<std::endl;
		file.clear();
		file.close();
		return true;
	}

	bool FEME2D::savePreCalMatrix(const char* _filename)
	{
		std::ofstream file;
		std::cout<<"Saving electrical matrix data...";
		file.open(_filename, std::ios::out|std::ios::binary);
		if (!this->is_loaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		file.write((char*)&this->num_node, 4);
		file.write((char*)&this->num_grounded, 4);
		for(int i=0;i<this->num_node;i++)
			file.write((char*)&this->node[i].state, sizeof(int));
		for(int i=0;i<this->num_node-this->num_grounded;i++){
			file.write((char*)&this->renum_A[i], sizeof(int));
		}
		file.write((char*)invK_A.X, invK_A.n * invK_A.m * sizeof(double));
		std::cout<<"[OK]"<<std::endl;
		file.clear();
		file.close();
		return true;
	}

	void FEME2D::loadConductivity(const char* _filename)
	{
		std::ifstream file;
		file.open(_filename, std::ios::in);

		std::cout<<"Loading parameters...";
		if(!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}
		for(int i=0; i<num_facet; i++){
			file>>this->facet[i].conductivity;
		}
		file.close();
		std::cout<<"[OK]"<<std::endl;
	}

	void FEME2D::saveConductivity(const char* _filename)
	{
		std::ofstream file;
		file.open(_filename, std::ios::out);

		std::cout<<"Saving parameters...";
		if(!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}
		for(int i=0; i<num_facet; i++){
			file<<this->facet[i].conductivity;
		}
		file.close();
		std::cout<<"[OK]"<<std::endl;
	}

	bool FEME2D::loadElectrodeSetting(const char* _filename)
	{
		this->is_matrix_available = false;
		Sleep(10);
		std::ifstream file;
		std::cout<<"Loading electrode setting...";
		file.open(_filename, std::ios::in);
		if (!this->is_loaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}

		//状態を一度クリアする
		this->num_grounded = 0;
		this->num_electrode = 0;
		this->num_sensing = 0;
		this->num_sensing_node = 0;
		this->num_candidate = 0;
		for(int i=0;i<num_node;i++){
			this->node[i].state = NONE;
			this->node[i].label = -1;
		}

		//ファイルの状態をセットする
		char buf[256];
		int tindex, tstate, tlabel;
		//最初の行には電極の数，状態番号，それらを構成するノードの数が入っている
		file>>buf;
		sscanf_s(buf, "%d,%d,%d", &tindex, &tstate,&tlabel);
		index_list_sensing.clear();
		index_list_sensing_each.clear();
		this->num_sensing = tindex;
		this->num_sensing_node = tlabel;
		index_list_sensing_each.resize(tindex);
		V1_c.resize(this->num_facet * num_sensing_node, 0);
		V2_c.resize(this->num_facet * num_sensing_node, 0);
		J_c.resize(this->num_facet * num_sensing_node, 0);
		J_cm.resize(this->num_facet * num_sensing, 0);
		potential_sim.resize(this->num_sensing * num_sensing_node, 0);
		potentialm_sim.resize(this->num_sensing * num_sensing, 0);
		potential_meas.resize(this->num_sensing * num_sensing_node, 0);
		potentialm_meas.resize(this->num_sensing * num_sensing, 0);

		while(!file.eof()){
			file>>buf;
			sscanf_s(buf, "%d,%d,%d", &tindex, &tstate,&tlabel);
			if(tindex >= 0 && tindex < this->num_node){
				this->node[tindex].state = tstate;
				this->node[tindex].label = tlabel;
				if(tstate == EITS::GROUNDED){
					this->num_grounded++;
				}
				else if(tstate== EITS::SENSING){
					index_list_sensing.push_back(tindex);
					index_list_sensing_each[tlabel].push_back(tindex);
				}
				else if(tstate== EITS::ELECTRODE){
					this->num_electrode++;
				}
			}
		}//確認：自由要素かつ電極となっていないかどうか
		index_list_free_element.clear();
		for(int i =0;i<this->num_facet;i++){
			//ノードがすべて検出設定であれば自由要素
			if(node[facet[i].index_node[0]].state == DETECTOR &&
				node[facet[i].index_node[1]].state == DETECTOR &&
				node[facet[i].index_node[2]].state == DETECTOR)
					index_list_free_element.push_back(i);
		}
#ifdef _DEBUG
		for(int i =0;i<this->index_list_sensing_each.size();i++){
			std::cout<<i<<": "<<this->index_list_sensing_each[i].size()<<std::endl;
		}
#endif
		J_cm_free.resize(index_list_free_element.size() * num_sensing, 0);

		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		std::cout<<"Number of elecrodes: "<<num_sensing<<"=="<<this->index_list_sensing_each.size()<<std::endl;
		std::cout<<"Number of sensing nodes: "<<num_sensing_node<<"=="<<this->index_list_sensing.size()<<std::endl;
		std::cout<<"Number of detector elements: "<<index_list_free_element.size()<<std::endl;

		return true;
	}

	bool FEME2D::saveElectrodeSetting(const char* _filename)
	{
		std::ofstream file;
		std::cout<<"Saving electrode setting...";
		file.open(_filename, std::ios::out);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		int num_es = 0;
		int max_lab = 0;
		for(int i=0;i<this->num_node;i++){
			if(this->node[i].state==SENSING){
				if(max_lab<=this->node[i].label){
					max_lab = this->node[i].label;
				}
				num_es++;
			}
		}

		file<<max_lab + 1<<","<<SENSING<<","<<num_es;
		for(int i=0;i<this->num_node;i++){
			if(this->node[i].state!=NONE){
				file<<std::endl<<i<<","<<this->node[i].state<<","<<this->node[i].label;
			}
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	//全ての候補についてシミュレーションしてデータを保存する
	bool FEME2D::saveSimulatedPotentialForAll(const char *_filename)
	{
		std::ofstream file;
		double t_potential;
		char buf[256];
		std::cout<<"Saving simulated potential...";
		file.open(_filename, std::ios::out);
		if (!this->is_loaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		this->num_candidate = this->num_node - this->num_sensing_node;
		this->num_electrode = 1;

		file<<this->num_sensing<<std::endl;
		file<<this->num_candidate<<std::endl;

		for(int i=0;i<this->num_sensing;i++){
			this->selectStateAtLabelOf(GROUNDED, i);
			this->calInverseMatrix();
			for(int j=0;j<this->num_node;j++){//接触点の候補
				if(this->node[j].state != GROUNDED && this->node[j].state != SENSING){
					for(int k=0;k<this->num_node;k++)
						this->node[k].state = NONE;
					this->node[j].state = ELECTRODE;
					//シミュレーションを実行
					this->setPotential(3.0);//TODO:今回は3.0V，計測システムの電源電圧にしたがって決める
					this->Potential2Potential();

					for(int k=0;k<this->num_sensing;k++){

						//電極の電圧は該当ノード電圧の平均値で与える
						t_potential = 0;
						for(int l=0;l<index_list_sensing_each[k].size();l++){
							t_potential+=node[index_list_sensing_each[k][l]].potential;
						}
						//ファイルに出力する
						file<<t_potential/index_list_sensing_each[k].size();
						if(k < this->num_sensing - 1)file<<",";
					}
					file<<std::endl;
				}
			}
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool FEME2D::loadSimulatedPotentialForAll(const char *_filename)
	{
		std::ifstream file;
		std::cout<<"Loading voltage set...";
		this->is_sensing_ready = false;
		file.open(_filename, std::ios::in);
		if (!this->is_loaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		std::string buf;
		std::string tmp1;

		simulatedPotential.clear();
		index_list_input.clear();
		file>>this->num_sensing;
		file>>this->num_candidate;
		for(int i=0;i<this->num_sensing;i++){
			for(int j=0;j<this->num_node;j++){
				if(this->node[j].state != GROUNDED && this->node[j].state != SENSING){
					file>>buf;
					std::stringstream tmp2(buf);
					while(std::getline(tmp2,tmp1,','))//カンマ区切りを分割
					{
						simulatedPotential.push_back(std::stod(tmp1));
					}
					if(i==0)index_list_input.push_back(j);
				}
			}
		}
		this->calMinMaxPotentialSim();

		measuredPotential.resize(this->num_sensing * this->num_sensing);
		for(int i =0;i<measuredPotential.size();i++)
			measuredPotential[i] = 0;
		this->error.resize(this->num_candidate);
		this->is_sensing_ready = true;

		std::cout<<"[OK]"<<std::endl;
		std::cout<<"Num Sensing: "<<this->num_sensing<<std::endl;
		std::cout<<"Num Candidate: "<<this->num_candidate<<std::endl;
		std::cout<<"Size: "<<this->num_sensing*this->num_sensing*this->num_candidate<<"=="<<simulatedPotential.size()<<std::endl;
		return true;
	}

	bool FEME2D::savePotentialSimulatedFor(int _condition, const char* _filename)
	{
		std::ofstream file;
		std::cout<<"Saving simulated potential...";
		file.open(_filename, std::ios::out);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i=0;i<this->index_list_sensing.size();i++){
			file<<this->potential_sim[this->index_list_sensing.size() * _condition + i]<<std::endl;
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool FEME2D::saveMeanPotentialSimulatedFor(int _condition, const char* _filename)
	{
		std::ofstream file;
		std::cout<<"Saving simulated mean potential...";
		file.open(_filename, std::ios::out);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		this->calMeanPotential();
		for(int i=0;i<this->num_sensing;i++){
			file<<this->potentialm_sim[this->num_sensing * _condition + i]<<std::endl;
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool FEME2D::loadMeanPotentialSimulatedFor(int _condition, const char* _filename)
	{
		std::ifstream file;
		std::cout<<"Loading simulated mean potential...";
		file.open(_filename, std::ios::in);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i=0;i<this->num_sensing;i++){
			file>>this->potentialm_sim[this->num_sensing * _condition + i];
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool FEME2D::loadPotentialSimulatedFor(int _condition, const char* _filename)
	{
		std::ifstream file;
		std::cout<<"Loading simulated potential...";
		file.open(_filename, std::ios::in);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i=0;i<this->index_list_sensing.size();i++){
			file>>this->potential_sim[this->index_list_sensing.size() * _condition + i];
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	bool FEME2D::loadPotentialMeasuredFor(int _condition, const char* _filename)
	{
		std::ifstream file;
		std::cout<<"Loading measured potential...";
		file.open(_filename, std::ios::in);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i=0;i<this->index_list_sensing.size();i++){
			file>>this->potential_meas[this->index_list_sensing.size() * _condition + i];
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	bool FEME2D::loadMeanPotentialMeasuredFor(int _condition, const char* _filename)
	{
		std::ifstream file;
		std::cout<<"Loading measured mean potential...";
		file.open(_filename, std::ios::in);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i=0;i<this->num_sensing;i++){
			file>>this->potentialm_meas[this->num_sensing * _condition + i];
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	bool FEME2D::saveJacobianFor(const char* _filename)
	{
		std::ofstream file;
		std::cout<<"Saving Jacobian...";
		file.open(_filename, std::ios::out);
		if (!this->is_loaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int j=0;j<this->num_facet;j++){
			for(int i=0;i<this->num_sensing_node;i++){
				file<<J_c[this->num_sensing_node * j + i]<<std::endl;
			}
		}
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool FEME2D::saveMeanJacobianFor(const char* _filename)
	{
		std::ofstream file;
		std::cout<<"Saving Mean Jacobian...";
		file.open(_filename, std::ios::out);
		if (!this->is_loaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		this->calMeanJacobian();
		for(int j=0;j<this->index_list_free_element.size();j++){
			for(int i=0;i<this->num_sensing;i++){
				file<<J_cm_free[this->num_sensing * j + i]<<std::endl;
			}
		}
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool FEME2D::loadJacobian(const char *_filename)
	{
		std::ifstream file;
		char filename[256];
		double t_J;
		this->is_jacobian_loaded = false;
		std::cout<<"Loading jacobian...";
		if(this->num_sensing_node < 1 || this->num_sensing < 1 || this->num_facet < 1 ){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		J.malloc(this->index_list_free_element.size(), this->num_sensing_node * this->num_sensing);
		tJ.malloc(this->num_sensing_node * this->num_sensing, this->index_list_free_element.size());
		Q.malloc(this->index_list_free_element.size(), this->index_list_free_element.size());
		Q.identity();
		dV.malloc(this->num_sensing_node * this->num_sensing);
		dSigma.malloc(this->index_list_free_element.size());

		for(int k=0;k<this->num_sensing;k++){
			std::cout<<".";
			sprintf_s(filename, "%s%d.csv", _filename, k);
			file.open(filename, std::ios::in);
			if (!this->is_loaded||!file.is_open()){
				std::cout<<"[FAIL]"<<std::endl;
				return false;
			}
			int count = 0;
			for(int j=0;j<this->num_facet;j++){
				for(int i=0;i<this->num_sensing_node;i++){
					file>>t_J;
					J(this->num_sensing_node * k + i, j) = t_J;
					tJ(j, this->num_sensing_node * k + i) = t_J;
					//std::cout<<t_J<<",";
					if(j == this->index_list_free_element[count]){
						tJ.X[this->num_sensing_node * k + i + this->num_sensing_node* this->num_sensing*count] = t_J;
						J.X[count+ this->index_list_free_element.size()*(this->num_sensing_node * k + i)] = t_J;
						if(i == this->num_sensing_node-1)count++;
					}
				}
			}
			file.close();
			file.clear();
		}
		std::cout<<"[OK]"<<std::endl;
		this->is_jacobian_loaded = true;
		return true;
	}
	
	bool FEME2D::loadInvJacobian(const char *_filename)
	{
		std::ifstream file;
		double t_J;

		this->is_sensing_ready = false;
		std::cout<<"Loading inverse jacobian...";
		if(this->num_sensing_node < 1 || this->num_sensing < 1 || this->num_facet < 1 ){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		int dim_t = this->index_list_free_element.size();
		int dim_s;
		if(mode==EIT)
			dim_s = this->num_sensing * (this->num_sensing-2);
		else if(mode==EMBIT)
			dim_s = this->num_sensing * (this->num_sensing-1);
		this->A_J.malloc(dim_t, dim_t);
		this->invA_J.malloc(dim_t, dim_t);
		this->invA_JtJ.malloc(dim_s, dim_t);

		file.open(_filename, std::ios::in|std::ios::binary);
		if(!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		file.read((char*)this->A_J.X,sizeof(double) * dim_t * dim_t);
		file.read((char*)this->invA_J.X,sizeof(double) * dim_t * dim_t);
		file.read((char*)this->invA_JtJ.X,sizeof(double) * dim_t * dim_s);
		file.close();
		this->is_sensing_ready = true;
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool FEME2D::saveInvJacobian(const char *_filename)
	{
		std::ofstream file;
		double t_J;
		std::cout<<"Saving inverse jacobian...";
		if(!this->is_sensing_ready){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		int dim_t = this->index_list_free_element.size();
		int dim_s;
		if(mode==EIT)
			dim_s = this->num_sensing * (this->num_sensing-2);
		else if(mode==EMBIT)
			dim_s = this->num_sensing * (this->num_sensing-1);
		file.open(_filename, std::ios::out|std::ios::binary);
		if(!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i = 0;i<dim_t * dim_t;i++)
			file.write((char*)&this->A_J.X[i],sizeof(double));
		for(int i = 0;i<dim_t * dim_t;i++)
			file.write((char*)&this->invA_J.X[i],sizeof(double));
		for(int i = 0;i<dim_t * dim_s;i++)
			file.write((char*)&this->invA_JtJ.X[i],sizeof(double));
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool FEME2D::loadMeanJacobian(const char *_filename)
	{
		std::ifstream file;
		char filename[256];
		double t_J;
		this->is_jacobian_loaded = false;
		std::cout<<"Loading mean jacobian...";
		if(this->num_sensing_node < 1 || this->num_sensing < 1 || this->num_facet < 1 ){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		//モードによって-2か-1かを変えないといけない
		int t_num_sensing;
		if(mode == EIT)t_num_sensing = this->num_sensing-2;
		else if(mode == EMBIT)t_num_sensing = this->num_sensing-1;
		J.malloc(this->index_list_free_element.size(), this->num_sensing * t_num_sensing);
		tJ.malloc(this->num_sensing * t_num_sensing, this->index_list_free_element.size());
		Q.malloc(this->index_list_free_element.size(), this->index_list_free_element.size());
		Q.identity();
		dV.malloc(this->num_sensing * t_num_sensing);
		dSigma.malloc(this->index_list_free_element.size());

		for(int k=0;k<this->num_sensing;k++){
			std::cout<<".";
			sprintf_s(filename, "%s%d.csv", _filename, k);
			file.open(filename, std::ios::in);
			if (!this->is_loaded||!file.is_open()){
				std::cout<<"[FAIL]"<<std::endl;
				return false;
			}
			for(int j=0;j<this->index_list_free_element.size();j++){
				int count = 0;
				for(int i=0;i<this->num_sensing;i++){
					file>>t_J;
					if(this->mode == EIT){
						if(i!=k&&i!=((k+1)%this->num_sensing)){//入力電極をスキップ
							tJ.X[t_num_sensing * k + count + this->num_sensing * t_num_sensing * j] = t_J;
							J.X[j+ this->index_list_free_element.size() * (t_num_sensing * k + count)] = t_J;
							count++;
						}
					}else{
						if(i!=k){//入力電極をスキップ
							tJ.X[t_num_sensing * k + count + this->num_sensing * t_num_sensing * j] = t_J;
							J.X[j+ this->index_list_free_element.size() * (t_num_sensing * k + count)] = t_J;
							count++;
						}
					}
				}
			}
			file.close();
			file.clear();
		}
		std::cout<<"[OK]"<<std::endl;
		this->is_jacobian_loaded = true;
		return true;
	}
	bool FEME2D::saveDSigma(const char* _filename)
	{
		std::ofstream file;
		char filename[256];
		std::cout<<"Savig defference of conductivity...";
		file.open(_filename, std::ios::out);
		if (!this->is_loaded||!file.is_open()||this->dSigma.n==0){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i=0;i<this->dSigma.n;i++){
			file<<this->dSigma.X[i]<<std::endl;
		}
		file.close();
		file.clear();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	bool FEME2D::loadDSigma(const char* _filename)
	{
		std::ifstream file;
		char filename[256];
		double t_dsigma;
		double max_sigma = -1000;
		double min_sigma = 1000;
		std::cout<<"Loading defference of conductivity...";
		file.open(_filename, std::ios::in);
		if (!this->is_loaded||!file.is_open()||this->dSigma.n==0){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i=0;i<this->dSigma.n;i++){
			file>>t_dsigma;
			if(max_sigma<t_dsigma)max_sigma = t_dsigma;
			if(min_sigma>t_dsigma)min_sigma = t_dsigma;
			this->dSigma.X[i] = t_dsigma;
		}
//		this->color_map.setParam(min_sigma, max_sigma, 0, 1);
		mSigma.x = min_sigma;
		mSigma.y = max_sigma;
		file.close();
		file.clear();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	void FEME2D::clearSimulatedPotential()
	{
		potentialm_sim.resize(this->num_sensing * num_sensing, 0);
	}

	void FEME2D::calStiffnessMatrix()
	{
		std::cout<<"Calculating stiffness matrix...";
		if(!this->is_loaded){
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}
		if(num_node > 1 ){
			K = Matrixd(num_node, num_node);
			K_A = Matrixd(num_node, num_node);
			K_AA = Matrixd(num_node, num_node);
			K_AB = Matrixd(num_node,num_node);
			K_BA = Matrixd(num_node,num_node);
			K_BB = Matrixd(num_node,num_node);
			invK_A = Matrixd(num_node, num_node);
			invK_AB = Matrixd(num_node, num_node);
			invK_BA = Matrixd(num_node,num_node);
			invK_AA = Matrixd(num_node, num_node);
			potential = VectorNd(num_node);
			potential_A = VectorNd(num_node);

			renum_A.resize(num_node);
			renum_AA.resize(num_node);

			this->init();
		}

		Matrixd temp_inv(3,3);
		for(int i=0;i<num_facet;i++){
			temp_inv=facet[i].calArea();
			for(int j=0;j<3;j++)
				for(int k=0;k<2;k++)
					facet[i].dN.X[2*j+k] = temp_inv.X[3*(k+1)+j];//dNを求める
			for(int j=0;j<3;j++){
				for(int k=0;k<3;k++){
					facet[i].K.X[3*j+k] = 0;
					for(int l=0;l<2;l++)
						facet[i].K.X[3*j+k] += facet[i].dN.X[2*j+l]*facet[i].dN.X[2*k+l];//要素剛性マトリクスに値を代入
					facet[i].K.X[3*j+k] *= facet[i].area * facet[i].conductivity;
				}
			}
			for(int j=0;j<3;j++)
				for(int k=0;k<3;k++)
					K.X[num_node * facet[i].index_node[j] + facet[i].index_node[k]] += facet[i].K.X[3*j+k];//剛性マトリクスに要素マトリクスを挿入
		}
		std::cout<<"[OK]"<<std::endl;
	}

	void FEME2D::setConductivityPurtabationAt(double _d_sigma, int _elem)
	{
		if(!this->is_loaded || this->num_node < 1 || num_facet < 1){
			return;
		}
		this->init();

		Matrixd temp_inv(3,3);
		memset(K.X, 0, sizeof(double) * this->num_node *this->num_node );
		for(int i=0;i<num_facet;i++){
			temp_inv = facet[i].calArea();
			for(int j = 0; j < 3; j++)
				for(int k=0; k < 2; k++)
					facet[i].dN.X[2 * j + k] = temp_inv.X[3 * (k + 1) + j];
			for(int j = 0; j < 3; j++){
				for(int k = 0; k < 3; k++){
					facet[i].K.X[3 * j + k] = 0;
					for(int l = 0; l < 2; l++)
						facet[i].K.X[3*j+k] += facet[i].dN.X[2*j+l]*facet[i].dN.X[2*k+l];
					if(i == _elem)//該当要素の導電率に摂動を与える
						facet[i].K.X[3*j+k] *= facet[i].area * (facet[i].conductivity + _d_sigma);
					else
						facet[i].K.X[3*j+k] *= facet[i].area * facet[i].conductivity;
				}
			}
			for(int j=0;j<3;j++)
				for(int k=0;k<3;k++)
					K.X[num_node * facet[i].index_node[j] + facet[i].index_node[k]] += facet[i].K.X[3*j+k];
		}
	}
	void FEME2D::calAJacobian(double _lambda, int _method)
	{
		std::cout<<"Calculating A-Jacobian...";
		this->is_sensing_ready = false;
		Matrixd JtJ = this->tJ*this->J;
		if(this->is_jacobian_loaded){
			this->lambda = _lambda;
			if(_method == GAUSSIAN){
				Q.identity();
			}
			else if(_method == LAPLACIAN){
				Q.identity();
			}
			else if(_method == NEWTON_TIKHONOV){
				for(int i=0;i<JtJ.m;i++)
					Q(i,i) = JtJ(i,i);
			}
			this->A_J = JtJ + (this->lambda * this->Q);
#ifdef _DEBUG
			for(int i=0;i<A_J.m;i++){
				for(int j=0;j<A_J.n;j++){
					if(A_J(i,j)==0){
						std::cout<<i<<","<<j<<std::endl;
					}
				}
			}
#endif
			//MEMO:メッシュサイズが大きすぎるとメモリリークが起きる
			this->invA_J = this->A_J.inv();
			this->invA_JtJ =  this->invA_J * this->tJ;
			this->is_sensing_ready = true;
			std::cout<<"[OK]"<<std::endl;
		}else{
			std::cout<<"[FAIL]"<<std::endl;
		}
	}

	void FEME2D::updateMinMax(bool _is_set_cmap){
		double min_sigma = 100000;
		double max_sigma = -100000;
		for(int i = 0;i<this->dSigma.n;i++){
			if(this->facet[this->index_list_free_element[i]].is_selected){
				if(max_sigma < dSigma(i))
					max_sigma = dSigma(i);
				if(min_sigma > dSigma(i))
					min_sigma = dSigma(i);
			}
		}
		mSigma.x = min_sigma;
		mSigma.y = max_sigma;
#ifdef _DEBUG
		std::cout<<"MAX: "<<max_sigma<<std::endl;
		std::cout<<"MIN: "<<min_sigma<<std::endl;
		std::cout<<"[OK]"<<std::endl;
#endif
		if(_is_set_cmap){
			color_map.setParam(min_sigma, max_sigma);
		}
	}

	void FEME2D::solveInverseProblem(bool _is_set_cmap)
	{
#ifdef _DEBUG
		std::cout<<"Solving inverse problem...";
#endif
		if(invA_JtJ.m > 0 && invA_JtJ.n > 0 ){
			for(int i = 0;i<this->num_sensing;i++){
				int count =0;
				for(int j = 0;j<this->num_sensing;j++){
					if(mode == EIT){
						if(j!=i&&j!=((i+1)%this->num_sensing)){
							dV.X[(this->num_sensing-2)*i+count]
								= this->potentialm_meas[this->num_sensing*i+j] - this->potentialm_sim[this->num_sensing*i+j];
							count++;
						}
					}
					else if(mode == EMBIT){
						if(j!=i){
							dV.X[(this->num_sensing-1)*i+count]
								= this->potentialm_meas[this->num_sensing*i+j];
							count++;
						}
					}
				}
			}
			dSigma = invA_JtJ * dV;
			this->updateMinMax(_is_set_cmap);
		}
	}

	double FEME2D::calResidualVector()
	{
		VectorNd Ax_y = this->J * this->dSigma - this->dV;
		return Ax_y.abs();
	}

	double FEME2D::calRegularizedSolution()
	{
		VectorNd Rx = this->Q * this->dSigma;
		return Rx.abs();
	}

	Vector3d FEME2D::calCenterOfContact(int _mode, double _ratio)
	{
		double S_total = 0;
		double x = 0;
		double y = 0;
		double z = 0;
		double sigma;

		/*//重心による方法
		for(int e = 0; e < this->index_list_free_element.size(); e++){
			if(_mode == MEASURED)sigma = (this->dSigma(e) - mSigma.x)/(mSigma.y - mSigma.x);
			else if(_mode == TARGET)sigma = this->facet[this->index_list_free_element[e]].conductivity;//Testデータを導電率に反映させた場合は使える
			if(sigma<0)sigma=0;
			for(int j = 0; j < 3; j++){
				x += this->facet[this->index_list_free_element[e]].vertex[j].x * sigma/3.0;
				y += this->facet[this->index_list_free_element[e]].vertex[j].y * sigma/3.0;
				z += this->facet[this->index_list_free_element[e]].vertex[j].z * sigma/3.0;
			}
			S_total += sigma;
		}
		if(_mode == MEASURED){
			this->position_touch = Vector3d(x/S_total, y/S_total, z/S_total);
#ifdef _DEBUG
			std::cout<<this->position_touch<<std::endl;
#endif
		}
		return Vector3d(x/S_total, y/S_total, z/S_total);
		*/
		/*//最大値による方法
		sigma = mSigma.x;
		int t_index_center = 0;
		for(int e = 0; e < this->index_list_free_element.size(); e++){
			if(sigma < this->dSigma(e)){
				t_index_center = this->index_list_free_element[e];
				sigma = this->dSigma(e);
			}
		}
		for(int j = 0; j < 3; j++){
			x += this->facet[t_index_center].vertex[j].x/3.0;
			y += this->facet[t_index_center].vertex[j].y/3.0;
			z += this->facet[t_index_center].vertex[j].z/3.0;
		}*/
		//最大値の_ratio*100%に入る領域の重心
		int t_count = 0;
		for(int e = 0; e < this->index_list_free_element.size(); e++){
			if(_mode == MEASURED)sigma = this->dSigma(e);
			else if(_mode == TARGET)sigma = this->facet[this->index_list_free_element[e]].conductivity;//Testデータを導電率に反映させた場合は使える
			if(mSigma.y * _ratio < sigma){
				for(int j = 0; j < 3; j++){
					x += this->facet[this->index_list_free_element[e]].vertex[j].x * sigma / 3.0;
					y += this->facet[this->index_list_free_element[e]].vertex[j].y * sigma / 3.0;
					z += this->facet[this->index_list_free_element[e]].vertex[j].z * sigma / 3.0;
				}
				S_total += sigma;
			}
		}
		this->position_touch = Vector3d(x/S_total, y/S_total, z/S_total);
		return this->position_touch;
	}

	void FEME2D::calInverseMatrix()
	{
		this->is_matrix_available = false;
		if(K.n==0)
			this->calStiffnessMatrix();
		int num_node_grounded = num_node-num_grounded;
		K_A = Matrixd(num_node_grounded, num_node_grounded);
		invK_A = Matrixd(num_node_grounded, num_node_grounded);
//		std::cout<<Kp<<std::endl;
		int count1 = 0;
		int count2 = 0;
		for(int i=0;i<num_node;i++){
			if(node[i].state != GROUNDED ){
				renum_A[count1]=i;
				for(int j=0;j<num_node;j++){
					if(node[j].state != GROUNDED ){
						K_A.X[num_node_grounded*count1+count2] = K.X[num_node*i+j];
						count2++;
					}				
				}
				count1++;
				count2=0;
			}
		}
		std::cout<<"Calculating inverse electrical matrix...";
		if(!this->is_loaded){
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}

		invK_A = K_A.inv();
		//std::cout<<invKp_A<<std::endl;
		if(1)
	//	if(MathCalMethod::Inverse_Gauss(&invKp_A,&Kp_A)!=0)
		{
			K_AA = Matrixd(num_node,num_node);
			invK_AA = Matrixd(num_node,num_node);
			K_AB = Matrixd(num_node,num_node);
			invK_AB = Matrixd(num_node,num_node);
			this->is_matrix_available = true;
			std::cout<<"[OK]"<<std::endl;
		}
		else
		{
			invK_A.free();
			K_A.free();
			this->is_matrix_available = false;
			std::cout<<"[FAIL]"<<std::endl;
		}
	}

	double FEME2D::getPotentialAt(Vector3d _position){
		double result=0;
		if(this->is_matrix_available){
			Matrixd t_position = Matrixd(3,4);
			for(int i=0;i<this->num_facet;i++){
				int j=0;
				if(EITS::isProjectedPointOnTriangle(_position, this->facet[i].normal[0], this->facet[i].vertex)){//面の内外判定
					result=facet[i].getLocalPotentialAt(_position);
					return result;
				}
			}
		}
		return 0;
	}

	double FEME2D::getSigmaAt(Vector3d _position){
		int count = 0;
		double result=0;
		Matrixd t_position = Matrixd(3,4);
		for(int i=0;i<this->num_facet;i++){
			if(i==this->index_list_free_element[count]){
				if(EITS::isProjectedPointOnTriangle(_position, this->facet[i].normal[0], this->facet[i].vertex)){//面の内外判定
					return this->dSigma[count];
				}
				count++;
			}
		}
		return 0;
	}

	void FEME2D::selectObject(Vector3d *_selected_coord, int _mode)
	{
		if(!this->is_loaded)
			return;
		Vector3d select_normal;
		bool is_po_wi_tr;
		static int index_current = 0;

		/*
		if(_mode==GLUT_ACTIVE_SHIFT){
			this->is_matrix_available = false;
			this->num_grounded = 0;
		}
		if(_mode==GLUT_ACTIVE_CTRL){
			this->num_electrode = 0;
		}
		for(int i=0;i<this->num_node;i++){
			if(_mode==GLUT_ACTIVE_SHIFT&&this->node[i].state!=ELECTRODE&&this->node[i].state!=SENSING)
				this->node[i].state = NONE;
			if(_mode==GLUT_ACTIVE_CTRL&&this->node[i].state!=GROUNDED&&this->node[i].state!=SENSING)
				this->node[i].state = NONE;
		}
		*/

		//TODO:ALTキーを押しているときはリセットしない
		if(_mode!=GLUT_ACTIVE_ALT){
			for(int i=0;i<this->num_facet;i++){
				this->facet[i].is_selected = false;
			}
		}

		for(int i=0;i<this->num_node;i++){
			for(int j=0;j<4;j++){
				select_normal=(_selected_coord[(j+1)%4]-_selected_coord[j])
					%(_selected_coord[j+4]-_selected_coord[j]);
				select_normal/=select_normal.abs();

				if(this->is_auto_scale)
					is_po_wi_tr = EITS::isViewNodeWithinTriangle(10.0*(this->node[i].vertex - this->center)/this->size.abs(), select_normal, _selected_coord[j]);
				else 
					is_po_wi_tr = EITS::isViewNodeWithinTriangle((this->node[i].vertex), select_normal, _selected_coord[j]);
//					is_po_wi_tr = EITS::is_view_nodeWithinTriangle((this->node[i].vertex + this->center)/10.0, select_normal, _selected_coord[j]);

				if(!is_po_wi_tr){
					if(j==3){
					/*
						if(_mode==GLUT_ACTIVE_SHIFT){
							this->node[i].state = GROUNDED;
							num_grounded++;
						}
						if(_mode==GLUT_ACTIVE_CTRL){
							this->node[i].state = ELECTRODE;
							this->num_electrode++;
						}
						*/
						if(_mode==GLUT_ACTIVE_ALT){
							for(int k=0;k<this->node[i].index_list_facet.size();k++){
								this->facet[this->node[i].index_list_facet[k]].is_selected = true;
							}
						}
					}
				}
				else{
					break;
				}
			}
		}
		if(_mode==GLUT_LEFT_BUTTON){
				index_current=0;
		}else{
				index_current++;
		}
	}

	void FEME2D::selectObject2(Vector3d *_selected_coord, int _mode)
	{
		if(!this->is_loaded)
			return;
		Vector3d select_normal;
		bool is_po_wi_tr;
		static int index_current = 0;

		//SHIFTキーを押しているときはリセットしない
		if(_mode != GLUT_ACTIVE_SHIFT){
			for(int i=0;i<this->num_facet;i++){
				this->facet[i].is_selected = false;
			}
		}

		for(int i=0;i<this->num_node;i++){
			for(int j=0;j<4;j++){
				select_normal=(_selected_coord[(j+1)%4]-_selected_coord[j])
					%(_selected_coord[j+4]-_selected_coord[j]);
				select_normal/=select_normal.abs();

				if(this->is_auto_scale)
					is_po_wi_tr = EITS::isViewNodeWithinTriangle(10.0*(this->node[i].vertex - this->center)/this->size.abs(), select_normal, _selected_coord[j]);
				else
					is_po_wi_tr = EITS::isViewNodeWithinTriangle((this->node[i].vertex), select_normal, _selected_coord[j]);
//					is_po_wi_tr = EITS::is_view_nodeWithinTriangle((this->node[i].vertex + this->center)/10.0, select_normal, _selected_coord[j]);

				if(!is_po_wi_tr){
					if(j==3){
						for(int k=0;k<this->node[i].index_list_facet.size();k++){
							this->facet[this->node[i].index_list_facet[k]].is_selected = true;
						}
					}
				}
				else{
					break;
				}
			}
		}
	}

	void FEME2D::setLabelAndStateToSelected(int _label, int _state)
	{
		for(int i=0;i<this->num_facet;i++){
			if(this->facet[i].is_selected){
				for(int j=0;j<3;j++){
					this->node[this->facet[i].index_node[j]].label = _label;
					this->node[this->facet[i].index_node[j]].state = _state;
				}
			}
		}
	}

	void FEME2D::clearLabel(int _state)
	{
		for(int j=0;j<this->num_node;j++){
			if(_state == NONE){
				this->node[j].label = -1;
				this->node[j].state = EITS::NONE;
			}
			else{
				if(_state == this->node[j].state){
					this->node[j].label = -1;
					this->node[j].state = EITS::NONE;
				}
			}
		}
	}

	std::vector<double> FEME2D::getLineProfile(int _resolution, Vector3d _v1, Vector3d _v2)
	{
		std::vector<double> dst;
		Vector3d t_pos;
		for(int i=0;i<_resolution;i++){
			t_pos = (_v2-_v1)*(i/(double)(_resolution-1))+_v1;
			dst.push_back(this->getSigmaAt(t_pos));
		}
		return dst;
	}
	void FEME2D::setMode(int _mode)
	{
		this->mode = _mode;
	}

	int FEME2D::getMode( )
	{
		return this->mode;
	}
	void FEME2D::calJacobianFor(int _condition)
	{
		EITS::utilities t;

		//全ての自由要素について：確認：自由要素はどうやって決めたか
		for(int i=0;i<this->index_list_free_element.size();i++){
			if(i == 0)
				t.setStartTime();
			//(0) 
			//境界条件_conditionにセット
			if(mode == EIT){
				this->setGroundAndElectrode(_condition);
				this->calJacobianForAt(_condition, this->index_list_free_element[i]);
			}else if(mode == EMBIT){
				this->setGroundAndElectrodeAt(_condition, this->index_list_free_element[i]);
				this->calBJacobianForAt(_condition, this->index_list_free_element[i]);
			}
			if(i == 0){
				t.setEndTime();
				std::cout<<"Estimated calculation time: "<<t.getDeltaTime() * this->index_list_free_element.size()<<" s"<<std::endl;
			}
		}
	}

	void FEME2D::calJacobianForAt(int _condition, int _elem)
	{
		int i;

		//(1) 摂動あり:TODO摂動の値
		double d_sigma = 1e-6;
		this->setConductivityPurtabationAt(d_sigma, _elem);
		//駆動境界の電圧値を3.3Vにセット
		this->setPotential(3.3);
		//電位分布を求める
		this->solve();
		//計測点の電位のみ取り出す
#ifdef _OPENMP
#pragma omp parallel private( i )
#pragma omp for
#endif
		for(i = 0; i<this->index_list_sensing.size();i++){
			V2_c[this->num_sensing_node * _elem + i] = this->node[this->index_list_sensing[i]].potential;
		}

		//(2) 摂動なし
		this->setConductivityPurtabationAt(0, _elem);
		//駆動境界の電圧値を3.3Vにセット
		this->setPotential(3.3);
		//電位分布を求める
		this->solve();
		//計測点の電位のみ取り出す
#ifdef _OPENMP
#pragma omp parallel private( i )
#pragma omp for
#endif
		for(i = 0; i < this->index_list_sensing.size(); i++){
			V1_c[this->num_sensing_node * _elem + i] = this->node[this->index_list_sensing[i]].potential;
			this->potential_sim[this->num_sensing_node * _condition + i] =  this->node[this->index_list_sensing[i]].potential;
		}
	

		//(3) ヤコビアンの計算
#ifdef _OPENMP
#pragma omp parallel private( i )
#pragma omp for
#endif
		for(i = 0; i < this->index_list_sensing.size(); i++){
			J_c[this->num_sensing_node * _elem + i] = 
				(this->V2_c[this->num_sensing_node * _elem + i] - this->V1_c[this->num_sensing_node * _elem + i] ) / d_sigma;
		}

		if(_elem == 0){
			double tJcm[16];
			//電極番号ごとに加算
			memset(tJcm,0,16*sizeof(double));
			for(int i=0;i<this->index_list_sensing.size();i++){
				tJcm[this->node[this->index_list_sensing[i]].label] += J_c[this->num_sensing_node * this->index_list_free_element[_elem] + i];
			}
			//各電極ごとのインデックスの数で割って平均化
			for(int i=0;i<this->index_list_sensing_each.size();i++){
				tJcm[i] /= this->index_list_sensing_each[i].size();
				std::cout<<tJcm[i]<<std::endl;
			}
		}

	}

	void FEME2D::calBJacobianForAt(int _condition, int _elem)
	{
		int i;
/*
		//(1) 摂動あり
		this->setConductivityPurtabationAt(0, _elem);
		//駆動境界の電圧値を2Vにセット
		this->setPotential(2.0);
		//電位分布を求める
		this->solve();
		//計測点の電位のみ取り出す
#ifdef _OPENMP
#pragma omp parallel private( i )
#pragma omp for
#endif
		for(i = 0; i<this->index_list_sensing.size();i++){
			V2_c[this->num_sensing_node * _elem + i] = this->node[this->index_list_sensing[i]].potential;
		}
*/
		//(2)摂動なし
		this->setConductivityPurtabationAt(0, _elem);
		//駆動境界の電圧値を1Vにセット
		this->setPotential(1.0);
		//電位分布を求める
		this->solve();
		//計測点の電位のみ取り出す
#ifdef _OPENMP
#pragma omp parallel private( i )
#pragma omp for
#endif
		for(i = 0; i < this->index_list_sensing.size(); i++){
			V1_c[this->num_sensing_node * _elem + i] = this->node[this->index_list_sensing[i]].potential;
			this->potential_sim[this->num_sensing_node * _condition + i] =  this->node[this->index_list_sensing[i]].potential;
		}

		//(3) ヤコビアンの計算
#ifdef _OPENMP
#pragma omp parallel private( i )
#pragma omp for
#endif
		for(i = 0; i < this->index_list_sensing.size(); i++){
			J_c[this->num_sensing_node * _elem + i] = this->V1_c[this->num_sensing_node * _elem + i] ;
				//this->V2_c[this->num_sensing_node * _elem + i]  - this->V1_c[this->num_sensing_node * _elem + i] ;
		}
	}

	void FEME2D::calMeanJacobian(){
		for(int j=0;j<this->index_list_free_element.size();j++){
			memset(&J_cm_free[this->num_sensing * j ], 0, sizeof(double) * this->num_sensing);
			//電極番号ごとに加算
			for(int i=0;i<this->index_list_sensing.size();i++){
				J_cm_free[this->num_sensing * j + this->node[this->index_list_sensing[i]].label] += J_c[this->num_sensing_node * this->index_list_free_element[j] + i];
			}
			//各電極ごとのインデックスの数で割って平均化
			for(int i=0;i<this->index_list_sensing_each.size();i++){
				J_cm_free[this->num_sensing * j + i] /= this->index_list_sensing_each[i].size();
			}
		}
	}

	void FEME2D::calMeanPotential(){
		int count = 0;
		for(int j=0;j<this->num_sensing;j++){
			memset(&this->potentialm_sim[this->num_sensing * j ], 0, sizeof(double) * this->num_sensing);
			for(int i=0;i<this->index_list_sensing.size();i++){
				this->potentialm_sim[this->num_sensing * j + this->node[this->index_list_sensing[i]].label] += this->potential_sim[this->num_sensing_node * j + i];
			}
			for(int i=0;i<this->index_list_sensing_each.size();i++){
				this->potentialm_sim[this->num_sensing * j + i] /= this->index_list_sensing_each[i].size();
			}
		}
	}

	void FEME2D::setPotential(double *_potential, std::vector<int> _list)
	{
		memset(this->potential.X, 0, sizeof(double) * this->potential.n);
		for(int i=0;i<_list.size();i++){
			this->node[_list[i]].state = ELECTRODE;
			this->potential.X[_list[i]] = _potential[i];
		}
	}

	void FEME2D::setPotential(double _potential)
	{
		int i;
#ifdef _OPENMP
#pragma omp parallel private( i )
#pragma omp for
#endif
		for(i = 0; i < this->num_node; i++){
			if(this->node[i].state == ELECTRODE)
				this->potential.X[i] = _potential;
			else this->potential.X[i] = 0;
		}
	}

	void FEME2D::setConductivityForSelected(double _conductivity)
	{
		for(int i = 0; i < this->num_facet; i++){
			if(this->facet[i].is_selected){
				this->facet[i].conductivity = _conductivity;
			}else{
				this->facet[i].conductivity = 1;
			}
		}
	}

	void FEME2D::updateConductivityWithDSigma()
	{
		int count = 0;
		for(int i=0;i<this->num_facet;i++){
			if(i==this->index_list_free_element[count]){
				this->facet[i].conductivity += this->dSigma(count);
				count++;
			}
		}
	}


	void FEME2D::Potential2Potential(){
		if(!this->is_matrix_available || !this->is_loaded)
			return;
		int count1=0;
		int count2=0;
		int count3=0;
		int count4=0;
		for(int i=0;i<num_node-num_grounded;i++){
			if(node[renum_A[i]].state != ELECTRODE ){	
				for(int j=0;j<num_node-num_grounded;j++){
					if(node[renum_A[j]].state == ELECTRODE){
						invK_BA.X[num_electrode*count1+count2]=invK_A.X[(num_node-num_grounded)*i+j];
						count2++;
					}
				}
				renum_AA[count1] = renum_A[i];			
				count1++;
				count2=0;
			}
			else{
				for(int j=0;j<num_node-num_grounded;j++){
					if(node[renum_A[j]].state == ELECTRODE){
						invK_AA.X[num_electrode * count3 + count4]=invK_A.X[(num_node - num_grounded) * i + j];
						count4++;				
					}
				}
				potential_A.X[count3]=potential.X[renum_A[i]];
				count3++;
				count4=0;
			}
		}
		if(num_electrode > 0){
			potential_A.n = num_electrode;
			invK_AA.n = invK_AA.m = num_electrode;
			K_AA = invK_AA.inv();
			invK_BA.m = num_node - num_grounded - num_electrode;
			invK_BA.n = num_electrode;
			potential_B = invK_BA*K_AA*potential_A;
			this->update();
/*
			if(this->is_viewVectorField)
				this->calCurCurrentDensity();
*/
		}
	}
	void FEME2D::solveFor(int _condition, double _voltage)
	{
		int i;
		//境界条件_conditionにセット
		//this->setGroundAndElectrodeAtSelected(_condition);
		this->setGroundAndElectrode( _condition );
		this->setConductivityPurtabationAt(0, 0);
		//駆動境界の電圧値を_voltageにセット
		this->setPotential(_voltage);
		//電位分布を求める
		this->solve();
		//計測点の電位のみ取り出す
#ifdef _OPENMP
#pragma omp parallel private( i )
#pragma omp for
#endif
		for(i = 0; i < this->index_list_sensing.size(); i++){
			this->potential_sim[this->num_sensing_node * _condition + i] =  this->node[this->index_list_sensing[i]].potential;
		}

	}

	void FEME2D::solveBoundFor(int _condition, double _voltage)
	{
		int i;
		//境界条件_conditionにセット
		this->setGroundAndElectrodeAtSelected(_condition);
		this->setConductivityPurtabationAt(0, 0);
		//駆動境界の電圧値を1Vにセット
		this->setPotential(_voltage);
		//電位分布を求める
		this->solve();
		//計測点の電位のみ取り出す
#ifdef _OPENMP
#pragma omp parallel private( i )
#pragma omp for
#endif
		for(i = 0; i < this->index_list_sensing.size(); i++){
			this->potential_sim[this->num_sensing_node * _condition + i] =  this->node[this->index_list_sensing[i]].potential;
		}
	}

	void FEME2D::solve()
	{
		if(!this->is_loaded)
			return;
		int count1=0;
		int count2=0;
		int count3=0;
		int count4=0;
		int t_size_a = num_electrode + num_grounded;
		int t_size_b = num_node - num_electrode - num_grounded;
		CPPL::dgematrix A(t_size_b, t_size_b);
		CPPL::dcovector y(t_size_b);
		CPPL::dgematrix t_B(t_size_b, t_size_a);
		CPPL::dcovector t_y(t_size_a);
		for(int i=0;i<num_node;i++){
			if(node[i].state != ELECTRODE && node[i].state != GROUNDED){	
				for(int j=0;j<num_node;j++){
					if(node[j].state != ELECTRODE && node[j].state != GROUNDED){
						A(count1, count2) = K.X[num_node * i + j];
						count2++;
					}
					else{
						t_B(count1 , count3) = K.X[num_node * i + j];
						count3++;				
					}
				}
				renum_AA[count1] = i;
				count1++;
				count2 = 0;
				count3 = 0;
			}
			else{
				t_y(count4) = potential.X[i];
				count4++;
			}
		}
		if(num_electrode > 0){
			y = - t_B * t_y;
			A.dgesv(y);
			int i;
			for(i = 0; i < t_size_b; i++)
				potential.X[renum_AA[i]] = y(i);
			for(i = 0; i < num_node; i++)
				node[i].potential=potential.X[i];
			for(i = 0; i < num_facet; i++){
				facet[i].potential.X[0] = potential.X[facet[i].index_node[0]];
				facet[i].potential.X[1] = potential.X[facet[i].index_node[1]];
				facet[i].potential.X[2] = potential.X[facet[i].index_node[2]];
			}
		}
	}

	void FEME2D::calCurrentDensity()
	{
		if(!this->is_matrix_available || !this->is_loaded)
			return;
		for(int i=0;i<this->num_facet;i++){
			this->facet[i].calCurrentDensity();
		}
	}

	//使っていない
	void FEME2D::selectStateAtLabelOf(int _state, int _label)
	{
		if(_state != GROUNDED && _state !=ELECTRODE)//UNSUPPORTED
			return;

		if(_state == GROUNDED)
			num_grounded = 0;
		else if(_state == ELECTRODE)
			num_electrode = 0;

		for(int i=0;i<this->num_node;i++){
			if(this->node[i].label != -1){
				if(_state == GROUNDED)
					this->node[i].state = SENSING;
			}
			if(this->node[i].label  == _label){
				if(_state == GROUNDED){
					this->node[i].state = _state;
					num_grounded ++;
				}
				else if(_state == ELECTRODE){
					if(this->node[i].state != GROUNDED){
						this->node[i].state = _state;
						num_electrode ++;
					}
				}
			}
		}
	}

	//EITで使う
	void FEME2D::setGroundAndElectrode(int _label)
	{
		num_grounded = 0;
		num_electrode = 0;

//		index_list_input.clear();
		for(int i=0;i<this->num_node;i++){
			this->node[i].state = NONE;
			if(this->node[i].label != -1){
				this->node[i].state = SENSING;
			}
			if(this->node[i].label == _label){
				this->node[i].state = GROUNDED;
				num_grounded ++;
			}
			else if(this->node[i].label == (_label + 1)%this->num_sensing){
				this->node[i].state = ELECTRODE;
				num_electrode ++;
			}
		}
	}

	//labelと一致するところをGNDに，要素のノードをELEに
	void FEME2D::setGroundAndElectrodeAt(int _label, int _elem)
	{
		num_grounded = 0;
		num_electrode = 0;

//		index_list_input.clear();

		for(int i=0;i<this->num_node;i++){
			this->node[i].state = NONE;
			if(this->node[i].label != -1){
				this->node[i].state = SENSING;
			}
			if(this->node[i].label == _label){
				this->node[i].state = GROUNDED;
				num_grounded ++;
			}
		}
		for(int i=0;i<3;i++){
			this->node[this->facet[_elem].index_node[i]].state = ELECTRODE;
			this->num_electrode++;
		}
	}

	void FEME2D::setGroundAndElectrodeAtSelected(int _label)
	{
		num_grounded = 0;
		num_electrode = 0;

//		index_list_input.clear();

		for(int i=0;i<this->num_node;i++){
			this->node[i].state = NONE;
			if(this->node[i].label != -1){
				this->node[i].state = SENSING;
			}
			if(this->node[i].label == _label){
				this->node[i].state = GROUNDED;
				num_grounded ++;
			}
		}

		for(int j=0;j<this->num_facet;j++){
			if(this->facet[j].is_selected){
				for(int i=0;i<3;i++){
					if(this->node[this->facet[j].index_node[i]].state != ELECTRODE)
						this->num_electrode++;
					this->node[this->facet[j].index_node[i]].state = ELECTRODE;
				}
			}
		}
	}
	void FEME2D::update()
	{
		if(!this->is_loaded)
			return;
		for(int i=0;i<num_node-num_grounded-num_electrode;i++)
			potential.X[renum_AA[i]] = potential_B.X[i];
		for(int i=0;i<num_node;i++)
			node[i].potential=potential.X[i];
		for(int i=0;i<num_facet;i++)
			for(int j=0;j<3;j++)
				facet[i].potential.X[j] = potential.X[facet[i].index_node[j]];
	}

	void FEME2D::render(void)
	{
		if(!this->is_loaded)
			return;

		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);
		glPushMatrix();

		if(this->is_auto_scale){
			glScaled(10.0/this->size.abs(), 10.0/this->size.abs(), 10.0/this->size.abs());
			glTranslated(-this->center.x, -this->center.y, -this->center.z);
		}
		/*
		else//モデルスケールによる調整なしの描画ではcm単位で描画する
			glScaled(1/10.0, 1/10.0, 1/10.0);
		*/

		//Lineの描画
		if(this->is_view_line && this->is_view_line_extracted){
			glDisable(GL_LIGHTING);
			glColor3d(0.1, 0.1, 0.1);
			glLineWidth(0.5);
			for(int i=0;i<num_line;i++){
				//接触点を含む線は描画しない
				//if(!node[this->line[i].nodeIndex[0]].isNodeContact&&!node[this->line[i].nodeIndex[1]].isNodeContact){
					glBegin(GL_LINE_STRIP);
					glVertex3dv(this->line[i].vertex[0].X);	
					glVertex3dv(this->line[i].vertex[1].X);	
					glEnd();
				//}
			}
			glEnable(GL_LIGHTING);
		}

		if(this->is_view_facet){
			//if(this->is_viewMisesStress)
			glDisable(GL_LIGHTING);
			int count = 0;
			for(int i=0;i<this->num_facet;i++)
			{

	//			if(this->surf[i].isCollidSurface){//接触面の表示
	//				this->faceMaterial.setAmbient(Vector3D(0.8,0,0.8));
	//				glColor3dv(Vector3D(1, 0, 1).X);
	//			}
				glBegin(GL_POLYGON);
	//			if(!this->is_view_smooth)
	//				glNormal3dv(facet[i].normal[0].X);
				for(int j=0;j<3;j++){
	//TODO::変位等の可視化
	//				if(this->is_viewMisesStress){
	//					this->map.glColorMap(this->surf[i].deform[j].abs());
	//					this->map.glColorMap(this->elem[this->surf[i].elemIndex].strain.X[3*j+0]*this->elem[this->surf[i].elemIndex].strain.X[3*j+0]+this->elem[this->surf[i].elemIndex].strain.X[3*j+1]*this->elem[this->surf[i].elemIndex].strain.X[3*j+1]+this->elem[this->surf[i].elemIndex].strain.X[3*j+2]*this->elem[this->surf[i].elemIndex].strain.X[3*j+2]);
	//					this->map.glColorMap(this->node[this->surf[i].nodeIndex[j]].deformation.abs());
	//				}
	//頂点法線
	//				if(this->is_view_smooth)
	//					glNormal3dv(this->facet[i].normal[j].X);

					glColor3dv(getJetColor(this->color_map.getNormalizedValue(facet[i].potential.X[j])).X);
					if(this->facet[i].is_selected)glColor3d(1,1,1);
					if(this->mode_view == CONDUCTIVITY ){
						if(i==this->index_list_free_element[count]){
							if(this->dSigma.n>0)
								glColor3dv(getJetColor(this->color_map.getNormalizedValue(this->dSigma(count))).X);
							if(j==2&&this->index_list_free_element.size()-1>count)count++;
						}
						else{
							glColor3d(0.5,0.5,0.5);
							//if(node[this->facet[i].index_node[0]].label != -1 )
							if(node[this->facet[i].index_node[0]].state == SENSING &&
								node[this->facet[i].index_node[1]].state == SENSING &&
								node[this->facet[i].index_node[2]].state == SENSING)
								glColor3d(1,1,1);
						}
					}
					glVertex3dv(this->facet[i].vertex[j].X);
				}
				glEnd();
				/*
				char text[256];
				sprintf(text,"%04d", i);
				glutDrawTextAt(text, (this->facet[i].vertex[0]+this->facet[i].vertex[1]+this->facet[i].vertex[2])/3.0);
				*/
			}
			glEnable(GL_LIGHTING);

		}
		
		if(this->is_view_node){
			glPointSize(5);
			glDisable(GL_LIGHTING );
			glBegin(GL_POINTS);
			for(int i=0;i<num_node;i++){
				if(node[i].state == GROUNDED ){glColor3d(1,0,0);}
				else if(node[i].state == ELECTRODE ){glColor3d(1,1,0);}
				else if(node[i].state == SENSING ){glColor3d(0,1,0);}
				else {glColor3d(0,1,1);}
				glVertex3dv(this->node[i].vertex.X);
			}
			glEnd();
			glEnable(GL_LIGHTING );
		}
		glPopMatrix();
		glDisable(GL_DEPTH_TEST );
		glDisable(GL_ALPHA_TEST);
	}

	void FEME2D::setMeasuredPotential(int _condition, double *_voltage)
	{
		//計測値を更新
		memcpy(&this->potentialm_meas[this->num_sensing * _condition], _voltage, sizeof(double) * this->num_sensing);
	}

	//この関数は使っていない（古い）
	int FEME2D::getOptimizedInputFor(int _condition, double *_voltage)
	{
		if(this->num_candidate>0 && this->index_list_input.size() == this->num_candidate && measuredPotential.size()!=0){
			//計測値を更新
			memcpy(&measuredPotential[this->num_sensing * _condition], _voltage, sizeof(double) * this->num_sensing);
			//最大最小値の算出
			this->calMinMaxPotentialMeas();
			//全ての候補に対して誤差を計算
			this->calErrorBtwSimMeas();
			//誤差が小さい候補トップNをスタックする
			if(error.size()>0){
				int t_num_top = 10;//TODO: N=10
				//最小要素を取得する
				std::vector<std::pair<double, int> > t_list;
				for (int i = 0; i < error.size(); i++) {
					t_list.push_back(std::pair<double, int>(error[i], i));
				}
				std::sort(t_list.begin(), t_list.end());

				Vector3d t_position=Vector3d(0,0,0);
				for(int i=0;i<t_num_top;i++){
					t_position+=2.0*(t_num_top - i)*node[index_list_input[t_list[i].second]].vertex/((double)(t_num_top+1)*t_num_top);
				}
				position_touch = t_position;
				if(this->is_auto_scale)
					position_touch*=10.0/this->scale;
				else
					position_touch*=1/10.0;
				return index_list_input[t_list[0].second];
			}
		}
	/*
				cout << showpos<<setprecision(4)<< setiosflags(ios::fixed);
				std::cout << "Error at min input node is " << meanf[index] << ","<< error[index] << std::endl;
				for(int i=0;i<this->numClip;i++){
					cout<<simulatedPotential[numClip*(numCandidate*_condition+index)+i]<<",";
				}
				cout<<endl;
				for(int i=0;i<this->numClip;i++){
					cout<<_voltage[i]<<",";
				}
				cout<<endl;
	*/
		return -1;
	}

	void FEME2D::calErrorBtwSimMeas()
	{
		int t_index1;
		double t_error;
		for(int i = 0; i < this->num_candidate; i++){
			t_error = 0;
			if(max_potential_meas!=0){
				for(int j = 0; j < this->num_sensing; j++){
					for(int k = 0; k < this->num_sensing; k++){
						t_index1 = this->num_sensing * (num_candidate * j + i) + k;
						t_error +=pow((simulatedPotential[t_index1] - min_potential_sim[i])/(max_potential_sim[i] - min_potential_sim[i])
							-(measuredPotential[this->num_sensing * j + k] - min_potential_meas)/(max_potential_meas - min_potential_meas), 2.0);
					}
				}
			}
			t_error /=  (double)this->num_sensing*this->num_sensing;
			error[i] = t_error;
		}
	}

	void FEME2D::calMinMaxPotentialSim()
	{
		int t_index1;
		min_potential_sim.resize(this->num_candidate);
		max_potential_sim.resize(this->num_candidate);
		for(int i = 0; i < this->num_candidate; i++){
			min_potential_sim[i] = 10000;
			max_potential_sim[i] = 0;
			for(int j = 0; j < this->num_sensing; j++){
				for(int k = 0; k < this->num_sensing; k++){
					t_index1 = this->num_sensing * (this->num_candidate * j + i ) + k;
					if(min_potential_sim[i] > simulatedPotential[t_index1] )
						min_potential_sim[i] = simulatedPotential[t_index1];
					if(max_potential_sim[i] < simulatedPotential[t_index1] )
						max_potential_sim[i] = simulatedPotential[t_index1];
				}
			}
		}
	}

	void FEME2D::calMinMaxPotentialMeas()
	{
		int t_index1;
		min_potential_meas = 10000;
		max_potential_meas = 0;
		for(int i = 0; i < this->num_sensing; i++){
			for(int j = 0; j < this->num_sensing; j++){
				t_index1 = this->num_sensing * j + i;
				if(min_potential_meas > measuredPotential[t_index1] )
					min_potential_meas = measuredPotential[t_index1];
				if(max_potential_meas < measuredPotential[t_index1] )
					max_potential_meas = measuredPotential[t_index1];
			}
		}
	}

	double FEME2D::calSpatialResolution(double _rate)
	{
		double RES = 0;
		double TOTAL = 0;
		//_rate * 100% of the maximum amplitude
		double HL = _rate * (mSigma.y - mSigma.x);
		for(int i=0; i<this->index_list_free_element.size(); i++){
//			if(this->facet[this->index_list_free_element[i]].conductivity == TARGET_CONDUCTIVITY){
				TOTAL+=this->facet[this->index_list_free_element[i]].area;
				if(this->dSigma[i] - mSigma.x > HL){
					RES+=this->facet[this->index_list_free_element[i]].area;
				}
//			}
		}
		if(TOTAL == 0)return 0;
		return RES/TOTAL;
	}

	double FEME2D::calShapeDeformation()
	{
		double RES50 = this->calSpatialResolution(0.5);
		double RES75 = this->calSpatialResolution(0.75);
		double HL = 0.25 * (mSigma.y - mSigma.x);
		double beta = atan(2 * HL / (RES75 - RES50));
		return beta;		
	}

	double FEME2D::calSizeError()
	{
		double RES = 0;
		double TOTAL = 0;
		//75% of the maximum amplitude
		double HL = 0.75 * (mSigma.y - mSigma.x);
		for(int i=0; i<this->index_list_free_element.size(); i++){
			if(this->facet[this->index_list_free_element[i]].conductivity == TARGET_CONDUCTIVITY){
				TOTAL+=this->facet[this->index_list_free_element[i]].area;
			}
			if(this->dSigma[i] - mSigma.x > HL){
				RES+=this->facet[this->index_list_free_element[i]].area;
			}
		}
		if(TOTAL == 0)return 0;
		return RES/TOTAL;
	}

	double FEME2D::calPositionError()
	{
		Vector3d P_m = this->calCenterOfContact();
		Vector3d P_t = this->calCenterOfContact(TARGET);//この部分は有効ではないかもしれない
		double error = (P_m - P_t).abs();
		return error;
	}
	double FEME2D::calPotentialError()
	{
		//未実装
		return 0;
	}


	/****************FEM3D******************/
	FEME3D::FEME3D()
	{
		this->num_elem = 0;
		this->num_node = 0;
		this->num_facet = 0;
		this->num_grounded = 0;
		this->num_electrode = 0;
		this->num_line = 0;
		this->newMesh();
		this->clear();
		this->clearView();

		color_map.setParam(0, 5);//TODO:アクセッサを用意してマップ範囲を変えられるように
	}

	FEME3D::~FEME3D()
	{
		this->deleteMesh();
	}

	void FEME3D::newMesh()
	{
		/*Memory domain for the FEM analisys*/
		this->node = new eNode[num_node + 2];
		this->elem = new eTetrahedra[num_elem + 2];
		this->line = new Line[this->num_line + 2];
		this->facet = new Facet[this->num_facet + 2];

		K = Matrixd(num_node,num_node);
		K_A = Matrixd(num_node,num_node);
		K_AA = Matrixd(num_node,num_node);
		invK_AB = Matrixd(num_node/2,num_node);
		invK_A = Matrixd(num_node,num_node);
		invK_AA = Matrixd(num_node/2,num_node/2);
		potential = VectorNd(num_node);
		potential_A = VectorNd(num_node/2);

		if(num_node > 1){
			renum_A.resize(num_node);
			renum_AA.resize(num_node);
		}
		this->init();
	}

	void FEME3D::deleteMesh()
	{
		this->is_loaded = false;
		delete []elem;
		delete []node;
		delete []line;
		delete []facet;

		this->clear();

		K.free();
		K_A.free();
		K_AA.free();
		K_BA.free();
		invK_AB.free();
		invK_BA.free();
		invK_A.free();
		invK_AA.free();
		potential.free();
		potential_A.free();
	}

	void FEME3D::clear()
	{
		this->is_matrix_available = false;
		this->is_loaded = false;
		this->is_auto_scale = true;
		this->is_view_line_extracted = false;
		this->is_view_facet_extracted = false;
		this->center = Vector3d(0, 0, 0);
		this->size = Vector3d(1, 1, 1);
/*
		this->curElectIndex=0;
		this->isMemLock=false;
		this->isAdjustedScale=true;
		this->isObjectLoaded=false;
		this->isBCSLoaded=false;
		this->isElecMatrixAvailable=false;
		this->isMechMatrixAvailable=false;
		this->is_view_lineExtracted=false;
		this->isElementNormal=false;
*/
	}

	void FEME3D::init()
	{
		this->curInputIndex=0;
		if(this->potential.n>0)
			memset(this->potential.X,0,sizeof(double)*this->potential.n);
		if(this->potential_A.n>0)
			memset(this->potential_A.X,0,sizeof(double)*this->potential_A.n);
		if(this->potential_B.n>0)
			memset(this->potential_B.X,0,sizeof(double)*this->potential_B.n);
		for(int i=0;i<num_elem;i++){
			this->elem[i].clear();
		}
		for(int i=0;i<num_node;i++){
			this->node[i].clear();
			renum_A[i]=0;
			renum_AA[i]=0;
		}
	}

	void FEME3D::calCenter()
	{
		Vector3d min(INT_MAX,INT_MAX,INT_MAX);
		Vector3d max(-INT_MAX,-INT_MAX,-INT_MAX);
		for(int i=0;i<this->num_node;i++){
			for(int k=0;k<3;k++){
				if(min.X[k]>this->node[i].vertex.X[k])
					min.X[k]=this->node[i].vertex.X[k];
				if(max.X[k]<this->node[i].vertex.X[k])
					max.X[k]=this->node[i].vertex.X[k];
			}
		}
		this->size.x=fabs(max.x-min.x);
		this->size.y=fabs(max.y-min.y);
		this->size.z=fabs(max.z-min.z);
		this->center=(max+min)/2.0;
	}

	void FEME3D::calLine()
	{
		std::cout<<"Extracting lines...";
		if(!this->is_loaded){
			this->is_view_facet_extracted = false;
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}

		Line t_line;
		delete []this->line;
		this->num_line = 0;
		this->line = new Line[this->num_elem * 6 +2];//4面体なので6本

		for(int i=0;i<this->num_elem;i++){
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					if(j!=k){
						t_line.index_node[0] = elem[i].index_node[j];
						t_line.index_node[1] = elem[i].index_node[k];
						t_line.vertex[0] = node[elem[i].index_node[j]].vertex;
						t_line.vertex[1] = node[elem[i].index_node[k]].vertex;
						if(!isSharedLine(t_line, this->line, this->num_line)){
							this->num_line++;
							this->line[this->num_line-1] = t_line;
						}
					}
				}
			}
		}
		this->is_view_line_extracted = true;
		std::cout<<"[OK]"<<std::endl;
	}

	void FEME3D::calFacet()
	{
		std::cout<<"Extracting facets...";
		if(!this->is_loaded){
			this->is_view_facet_extracted = false;
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}
		Facet t_facet;
		t_facet.setFacetTypeAsTriangle();
		delete []this->facet;
		this->num_facet = 0;
		this->facet= new Facet[this->num_elem * 4];
		for(int i=0; i<this->num_elem; i++){
			for(int j=0; j<4; j++){
				t_facet.normal[j] = elem[i].normal[j];
				t_facet.index_node[0] = elem[i].index_node[(j+0)%4];
				t_facet.index_node[1] = elem[i].index_node[(j+1)%4];
				t_facet.index_node[2] = elem[i].index_node[(j+2)%4];
				t_facet.vertex[0] = node[elem[i].index_node[(j+0)%4]].vertex;
				t_facet.vertex[1] = node[elem[i].index_node[(j+1)%4]].vertex;
				t_facet.vertex[2] = node[elem[i].index_node[(j+2)%4]].vertex;

				//std::cout<<t_facet.vertex[0]<<std::endl;

				if(!isSharedFacet(t_facet, this->facet, this->num_facet)){
					this->num_facet++;
					t_facet.index_facet = this->num_facet-1;
					t_facet.index_elem = i;
					this->facet[this->num_facet-1] = t_facet;
				}
			}
		}
		this->is_view_facet_extracted = true;
		std::cout<<"[OK]"<<std::endl;
	}

	bool FEME3D::load(const char *_filename)
	{
		std::ifstream file;
		char buf[512];
		char dummy[256];
		int num;
		std::cout<<"Loading mesh...";
		file.open(_filename, std::ios::in);
		if(!file.is_open())
		{
			std::cout<<"[FAIL]"<<_filename<<std::endl;
			return false;
		}
		this->deleteMesh();
		while(file.getline(buf, sizeof(buf))){
			if(strstr(buf, "nNodes ") || strstr(buf, "nVertex ")){
				sscanf_s(buf, "%s %d", dummy, (unsigned)sizeof(dummy), &num);
				num_node = num;
				continue;
			}
			if(strstr(buf, "nTetrahedra ") || strstr(buf, "nTetrahedron ")){
				sscanf_s(buf, "%s %d", dummy, (unsigned)sizeof(dummy), &num);
				num_elem = num;
				break;
			}
		}
		this->newMesh();

		while(file.getline(buf, sizeof(buf)))
		{
			if(strstr(buf, "# Data ")) break;
		}
		while(!strstr(buf, "@1"))
		{
			file >> buf;
			if(file.eof()){
				std::cout<<"[FAIL]"<<std::endl;
				return false;
			}
		}
		for(int i=0; i< num_node; i++)
		{
			file 
				>> node[i].vertex.x 
				>> node[i].vertex.y 
				>> node[i].vertex.z;
		}
		while(!strstr(buf, "@2"))
		{
			file >> buf;
			if(file.eof()){
				std::cout<<"[FAIL]"<<std::endl;
				return false;
			}
		}
		for(int i=0; i< num_elem; i++){
			for(int j=0; j<4; j++){
				file >> elem[i].index_node[j];
				elem[i].index_node[j]--;
			}
		}
		for(int i=0;i<num_elem;i++)
			for(int j=0;j<4;j++)
				for(int k=0;k<3;k++)
					elem[i].position.X[3*j+k] = node[elem[i].index_node[j]].vertex.X[k];

		while(!strstr(buf, "@3")){
			file >> buf;
			if(file.eof()){
				std::cout<<"[FAIL]"<<std::endl;
				return false;
			}
		}
		for(int i=0; i< num_elem; i++)
			file >> elem[i].index_material;
		file.clear();
		file.close();

		this->is_loaded = true;
		this->calCenter();

		this->calLine();
		this->calFacet();

		//外部関数にする
		/*
		this->calObjectLine();
		this->calObjectSurfacePolygon();
		this->calObjectSurfaceNormal();
		this->calObjectSurfaceCenter();
		this->setNeiborSurface();
		this->setNeiborSurfaceNode();
		this->setNeiborElement();
		this->calObjectNewSurfaceNormal();
		this->calObjectNewNodeNormal();
*/
		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool FEME3D::loadPreCalMatrix(const char* _filename)
	{
		this->is_matrix_available = false;
		Sleep(10);

		std::ifstream file;
		int n;
		std::cout<<"Loading electrical matrix data...";
		file.open(_filename, std::ios::in|std::ios::binary);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		file.read((char*)&n,4);
		if(this->num_node != n){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		this->num_node = n;
		file.read((char*)&this->num_grounded,4);

		K_A = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		invK_A = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		K_AA = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		K_AB = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		K_BA = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		invK_AA = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		invK_AB = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);
		invK_BA = Matrixd(this->num_node-this->num_grounded,this->num_node-this->num_grounded);

		for(int i=0;i<this->num_node;i++){
			file.read((char*)&this->node[i].state,sizeof(int));
		}
		for(int i=0;i<this->num_node-this->num_grounded;i++){
			file.read((char*)&this->renum_A[i],sizeof(int));
		}
		file.read((char*)invK_A.X,invK_A.n*invK_A.m*sizeof(double));

		if(K.n!=0){
			K.free();
			K.n=0;
		}
		this->is_matrix_available = true;
		std::cout<<"[OK]"<<std::endl;
		file.clear();
		file.close();
		return true;
	}

	bool FEME3D::savePreCalMatrix(const char* _filename)
	{
		std::ofstream file;
		std::cout<<"Saving electrical matrix data...";
		file.open(_filename, std::ios::out|std::ios::binary);
		if (!this->is_loaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		file.write((char*)&this->num_node,4);
		file.write((char*)&this->num_grounded,4);
		for(int i=0;i<this->num_node;i++)
			file.write((char*)&this->node[i].state,sizeof(int));
		for(int i=0;i<this->num_node - this->num_grounded;i++)
			file.write((char*)&this->renum_A[i], sizeof(int));
		file.write((char*)invK_A.X,invK_A.n*invK_A.m*sizeof(double));
		std::cout<<"[OK]"<<std::endl;
		file.clear();
		file.close();
		return true;
	}

	bool FEME3D::loadElectrodeSetting(const char* _filename)
	{
		this->is_matrix_available = false;
		Sleep(10);
		std::ifstream file;
		std::cout<<"Loading electrode setting...";
		file.open(_filename, std::ios::in|std::ios::binary);
		if (!this->is_loaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}

		//状態を一度クリアする
		this->num_grounded = 0;
		this->num_electrode = 0;
		for(int i=0;i<num_node;i++){
			this->node[i].state = NONE;
		}

		//ファイルの状態をセットする
		char buf[256];
		int tindex,tstate;
		while(!file.eof()){
			file>>buf;
			sscanf_s(buf, "%d,%d", &tindex, &tstate);
			this->node[tindex].state = tstate;
			if(tstate == EITS::GROUNDED){
				this->num_grounded++;
			}
			else if(tstate== EITS::ELECTRODE){
				this->num_electrode++;
			}
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;

		return true;
	}

	bool FEME3D::saveElectrodeSetting(const char* _filename)
	{
		std::ofstream file;
		std::cout<<"Saving electrode setting...";
		file.open(_filename, std::ios::out|std::ios::binary);
		if (!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i=0;i<this->num_node;i++){
			if(this->node[i].state!=NONE){
				file<<i<<","<<this->node[i].state<<std::endl;
			}
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}


/*

	bool FEME3D::loadElectrodeIndexList(const char* _filename)
	{
		this->isElecMatrixAvailable=false;
		Sleep(10);
		ifstream file;
		std::cout<<"Loading electrode index list...";
		file.open(_filename, ios::in|ios::binary);
		if (!this->isObjectLoaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		char buf[256];
		int tindex1,tindex2;
		//ファイルの状態をセットする
		this->electIndexList.clear();
		this->groundIndexList.clear();
		while(!file.eof()){
			file>>buf;
			sscanf(buf, "%d,%d", &tindex1, &tindex2);
			this->electIndexList.push_back(tindex1);
			this->groundIndexList.push_back(tindex2);
		}
		file.clear();
		file.close();
		this->curElectIndex=0;
		std::cout<<"[OK]"<<std::endl;
	//	this->setGroundElectrodeList();
	//	this->calInverseElecMatrix();
		return true;
	}
	bool FEME3D::savePotentialList(const char* _filename)
	{
		ofstream file;
		char buf[256];
		std::cout<<"Saving results...";
		file.open(_filename, ios::out|ios::binary);
		if (!this->isObjectLoaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		for(int i=0;i<this->num_node;i++){
			if(this->node[i].state==SensingNodeON){
				file<<i<<","<<this->node[i].electrodeIndex<<","<<this->getPotentialAt(this->node[i].coord)<<endl;
			}
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	bool FEME3D::loadMeasuredVoltageSet(const char *_filename)
	{
		this->isElecMatrixAvailable=false;
		Sleep(10);
		ifstream file;
		std::cout<<"Loading voltage set...";
		file.open(_filename, ios::in|ios::binary);
		if (!this->isObjectLoaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		char buf[256];
		int tindex1,tindex2;
		double tvoltage;
		//ファイルの状態をセットする
		mvoltage.clear();
		mvIndexList.clear();
		while(!file.eof()){
			file>>buf;
			sscanf(buf, "%d,%d,%lf", &tindex1, &tindex2, &tvoltage);
			mvoltage.push_back(tvoltage);
			mvIndexList.push_back(tindex1);
		}
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		this->isElecMatrixAvailable=true;
		return true;
	}
	bool FEME3D::saveSimulatedPotentialForAll(const char *_filename)
	{
		ofstream file;
		double *tPotential;
		int *tNumPotential;
		char buf[256];
		std::cout<<"Saving simulated potential..."<<endl;
		file.open(_filename, ios::out|ios::binary);
		if (!this->isObjectLoaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		tPotential=new double[this->numClip]();
		tNumPotential=new int[this->numClip]();
		file<<this->numClip<<endl;
		file<<this->num_node-this->electIndexList.size()<<endl;
		this->numElectrode=1;
		for(int i=0;i<this->numClip;i++){
			this->selectGroundElectrodesAt(i);
			this->calInverseElecMatrix();
			for(int j=0;j<this->num_node;j++){//TODO:いずれはノードの削減を行う必要がある（表面を抽出した結果や裏面などを考慮する）
				if(this->node[j].state!=GroundNodeON&&this->node[j].state!=ElectrodeNodeON&&this->node[j].state!=SensingNodeON){
					//シミュレーションを実行
					for(int k=0;k<this->num_node;k++)
						this->node[k].isNodeElectrode=false;
					this->node[j].isNodeElectrode=true;
					this->setPotential(5);//TODO:今回は5V，計測システムの電源電圧にしたがって決める
					this->Potential2Potential();
					//各電極における電圧の平均値を求める
					for(int k=0;k<this->numClip;k++){
						tPotential[k]=0;
						tNumPotential[k]=0;
					}
					for(int k=0;k<this->electIndexList.size();k++){
						tPotential[node[this->electIndexList[k]].electrodeIndex]+=node[this->electIndexList[k]].potential;
						tNumPotential[node[this->electIndexList[k]].electrodeIndex]++;
					}
					//ファイルに出力する
					for(int k=0;k<this->numClip;k++){
						file<<tPotential[k]/tNumPotential[k];
						if(k<this->numClip-1)file<<",";
					}
					file<<endl;
				}
			}
		}
		delete []tPotential;
		delete []tNumPotential;
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	bool FEME3D::loadSimulatedPotentialForAll(const char *_filename)
	{
		ifstream file;
		std::cout<<"Loading voltage set...";
		file.open(_filename, ios::in|ios::binary);
		if (!this->isObjectLoaded||!file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		std::string buf;
		std::string tmp1;

		simulatedPotential.clear();
		inputIndexList.clear();
		file>>this->numClip;
		file>>this->numCandidate;
		for(int i=0;i<this->numClip;i++){
			for(int j=0;j<this->num_node;j++){//TODO:いずれはノードの削減を行う必要がある（表面を抽出した結果や裏面などを考慮する）
				if(this->node[j].state!=GroundNodeON&&this->node[j].state!=ElectrodeNodeON&&this->node[j].state!=SensingNodeON){
					file>>buf;
					std::is_tringstream tmp2(buf);
					while(std::getline(tmp2,tmp1,','))//カンマ区切りを分割
					{
						simulatedPotential.push_back(std::stod(tmp1));
					}
					if(i==0)inputIndexList.push_back(j);
				}
			}
		}
		std::cout<<"[OK]"<<std::endl;
		return true;
	}
*/
	void FEME3D::loadParameter(const char* _filename)
	{
		std::ifstream file;
		file.open(_filename, std::ios::in);

		std::cout<<"Loading parameters...";
		if(!this->is_loaded || !file.is_open()){
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}
//TODO:ラベルを利用してデータ量を減らす
		for(int i=0;i<num_elem;i++)
		{
			file>>this->elem[i].conductivity;
		}
		file.close();
		std::cout<<"[OK]"<<std::endl;
	}

	void FEME3D::calStiffnessMatrix()
	{
		std::cout<<"Calculating stiffness matrix...";
		if(!this->is_loaded){
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}

		this->K = Matrixd(num_node,num_node);
		Matrixd temp_inv(4,4);
		for(int i=0;i<num_elem;i++){
			temp_inv=elem[i].calVolume();
			for(int j=0;j<4;j++)
				for(int k=0;k<3;k++)elem[i].dN.X[3*j+k]=temp_inv.X[4*(k+1)+j];//dNを求める
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					elem[i].K.X[4*j+k]=0;
					for(int l=0;l<3;l++)elem[i].K.X[4*j+k]+=elem[i].dN.X[3*j+l]*elem[i].dN.X[3*k+l];//要素剛性マトリクスに値を代入
					elem[i].K.X[4*j+k]*=elem[i].volume*elem[i].conductivity;
				}
			}
			for(int j=0;j<4;j++)
				for(int k=0;k<4;k++)
					K.X[num_node*elem[i].index_node[j]+elem[i].index_node[k]]+=elem[i].K.X[4*j+k];//剛性マトリクスに要素マトリクスを挿入
		}
		std::cout<<"[OK]"<<std::endl;
	}

	void FEME3D::calInverseMatrix()
	{
		this->is_matrix_available = false;
		if(K.n==0)
			this->calStiffnessMatrix();

		K_A = Matrixd(num_node-num_grounded, num_node-num_grounded);
		invK_A = Matrixd(num_node-num_grounded, num_node-num_grounded);
//		std::cout<<Kp<<std::endl;
		int count1 = 0;
		int count2 = 0;
		for(int i=0;i<num_node;i++){
			if(node[i].state != GROUNDED ){
				renum_A[count1]=i;
				for(int j=0;j<num_node;j++){
					if(node[j].state != GROUNDED ){
						K_A.X[(num_node-num_grounded)*count1+count2]=K.X[num_node*i+j];
						count2++;
					}				
				}
				count1++;
				count2=0;
			}
		}
		std::cout<<"Calculating inverse electrical matrix...";
		if(!this->is_loaded){
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}

		invK_A = K_A.inv();
		//std::cout<<invKp_A<<std::endl;
		if(1)
	//	if(MathCalMethod::Inverse_Gauss(&invKp_A,&Kp_A)!=0)
		{
			K_AA = Matrixd(num_node,num_node);
			invK_AA = Matrixd(num_node,num_node);
			K_AB = Matrixd(num_node,num_node);
			invK_AB = Matrixd(num_node,num_node);
			K_BA = Matrixd(num_node,num_node);
			invK_BA = Matrixd(num_node,num_node);
			this->is_matrix_available = true;
			std::cout<<"[OK]"<<std::endl;
		}
		else
		{
			invK_A.free();
			K_A.free();
			this->is_matrix_available = false;
			std::cout<<"[FAIL]"<<std::endl;
		}
	}

	//電圧入力の場合の計算
	void FEME3D::setPotential(double *_potential, std::vector<int> _list)
	{
		memset(this->potential.X, 0, sizeof(double) * this->potential.n);
		for(int i=0;i<_list.size();i++){
			this->node[_list[i]].state = ELECTRODE;
			this->potential.X[_list[i]] = _potential[i];
		}
	}
	void FEME3D::setPotential(double _potential)
	{
		for(int i=0;i<this->num_node;i++){
			if(this->node[i].state == ELECTRODE)
				this->potential.X[i] = _potential;
			else this->potential.X[i] = 0;
		}
	}

	void FEME3D::Potential2Potential(){
		if(!this->is_matrix_available || !this->is_loaded)
			return;
		int count1=0;
		int count2=0;
		int count3=0;
		int count4=0;
		for(int i=0;i<num_node-num_grounded;i++){
			if(node[renum_A[i]].state != ELECTRODE ){	
				for(int j=0;j<num_node-num_grounded;j++){
					if(node[renum_A[j]].state == ELECTRODE){
						invK_BA.X[num_electrode*count1+count2]=invK_A.X[(num_node-num_grounded)*i+j];
						count2++;
					}
				}
				renum_AA[count1] = renum_A[i];			
				count1++;
				count2=0;
			}
			else{
				for(int j=0;j<num_node-num_grounded;j++){
					if(node[renum_A[j]].state == ELECTRODE){
						invK_AA.X[num_electrode * count3 + count4]=invK_A.X[(num_node - num_grounded) * i + j];
						count4++;				
					}
				}
				potential_A.X[count3]=potential.X[renum_A[i]];
				count3++;
				count4=0;
			}
		}
		if(num_electrode>0){
			potential_A.n=num_electrode;
			invK_AA.n=invK_AA.m=num_electrode;
			K_AA=invK_AA.inv();
			invK_BA.m=num_node-num_grounded-num_electrode;
			invK_BA.n=num_electrode;
			potential_B=invK_BA*K_AA*potential_A;
			this->update();
/*
			if(this->is_viewVectorField)
				this->calCurCurrentDensity();
*/
		}
	}
	void FEME3D::calCurrentDensity()
	{
		if(!this->is_matrix_available || !this->is_loaded)
			return;
		if(this->is_view_facet){
			for(int i=0;i<this->num_facet;i++){
				this->elem[this->facet[i].index_elem].calCurrentDensity();
			}
		}else{
			for(int i=0;i<this->num_elem;i++){
				this->elem[i].calCurrentDensity();
			}
		}
	}

	void FEME3D::update()
	{
		if(!this->is_matrix_available||!this->is_loaded)
			return;
		for(int i=0;i<num_node-num_grounded-num_electrode;i++)
			potential.X[renum_AA[i]]=potential_B.X[i];
		for(int i=0;i<num_node;i++)
			node[i].potential=potential.X[i];
		for(int i=0;i<num_elem;i++)
			for(int j=0;j<4;j++)
				elem[i].potential.X[j]=potential.X[elem[i].index_node[j]];
/*
		for(int i=0;i<num_facet;i++)
			for(int j=0;j<3;j++)
				facet[i].potential[j] = potential.X[surf[i].nodeIndex[j]];
*/
	}

	void FEME3D::selectObject(Vector3d *_selected_coord, int _mode)
	{
		if(!this->is_loaded)
			return;
		Vector3d select_normal;
		bool is_po_wi_tr;
		static int index_current = 0;

		if(_mode==GLUT_ACTIVE_SHIFT){
			this->num_grounded = 0;
		}
		if(_mode==GLUT_ACTIVE_CTRL){
			this->num_electrode = 0;
		}

		for(int i=0;i<this->num_node;i++){
			if(_mode==GLUT_ACTIVE_SHIFT&&this->node[i].state!=ELECTRODE)
				this->node[i].state = NONE;
			if(_mode==GLUT_ACTIVE_CTRL&&this->node[i].state!=GROUNDED)
				this->node[i].state = NONE;
		}

		for(int i=0;i<this->num_node;i++){
			for(int j=0;j<4;j++){
				select_normal=(_selected_coord[(j+1)%4]-_selected_coord[j])
					%(_selected_coord[j+4]-_selected_coord[j]);
				select_normal/=select_normal.abs();

				if(this->is_auto_scale)
					is_po_wi_tr = EITS::isViewNodeWithinTriangle(10.0*(this->node[i].vertex + this->center)/this->size.abs(), select_normal, _selected_coord[j]);
				else 
					is_po_wi_tr = EITS::isViewNodeWithinTriangle((this->node[i].vertex + this->center)/10.0, select_normal, _selected_coord[j]);

				if(!is_po_wi_tr){
					if(j==3){
						if(_mode==GLUT_ACTIVE_SHIFT){
							this->node[i].state = GROUNDED;
							num_grounded++;
						}
						if(_mode==GLUT_ACTIVE_CTRL){
							this->node[i].state = ELECTRODE;
							this->num_electrode++;
						}
					}
				}
				else{
					break;
				}
			}
		}
		if(_mode==GLUT_LEFT_BUTTON){
				index_current=0;
		}else{
				index_current++;
		}
	}

/*
	void FEME3D::setNodeSelectionMode(int _mode)
	{
		this->nodeSelectionMode=_mode;
	}
	void FEME3D::setNodeMode(Vector3D _coord)
	{
		if(!this->isObjectLoaded)
			return;

		Vector3D temp_vector;
		int i=0;
		while(i<num_node){
			if(this->isAdjustedScale)
				temp_vector=_coord-10*(node[i].coord+node[i].deformation+this->worldCenter)/this->worldSize;
			else
				temp_vector=_coord-(node[i].coord+node[i].deformation+this->worldCenter)/10.0;

			if(temp_vector.abs()<0.05){
				if(this->nodeSelectionMode==FixedNodeON&&!node[i].isNodeFixed){//0->1::前カウントしたのを除く
					this->clickedCoord=Vector3D(0,0,0);
					numFixed++;
					node[i].isNodeFixed=true;
					cout << "Fixed node #" << i<< endl;
				}
				if(this->nodeSelectionMode==FixedNodeOFF&&node[i].isNodeFixed){//1->0::前カウントしたのを除く
					this->clickedCoord=Vector3D(0,0,0);
					numFixed--;
					node[i].isNodeFixed=false;
					cout << "Free node #" << i<< endl;
				}
				if(this->nodeSelectionMode==ContactNodeON&&!node[i].isNodeContact){//0->1::前カウントしたのを除く
					this->clickedCoord=10*(node[i].coord+this->worldCenter)/this->worldSize;
					numContact++;
					node[i].isNodeContact=true;
					cout << "Contact node #" << i<< endl;
				}
				if(this->nodeSelectionMode==ContactNodeOFF&&node[i].isNodeContact){//1->0::前カウントしたのを除く
					this->clickedCoord=Vector3D(0,0,0);
					numContact--;
					node[i].isNodeContact=false;
					cout << "Free node #" << i<< endl;
				}
				break;
			}
			else i++;
		}
	}
	void FEME3D::selectNodes(int _mode, int _button)
	{
		Vector3D selectNormal;
		bool isPoWiTr;
		this->setNodeSelectionMode(_mode);
		static int currentIndex=0;

		//特定の条件だけリセットしたほうがよい
		if(_button==GLUT_LEFT_BUTTON){
			if(this->nodeSelectionMode==FixedNodeON)
				this->numFixed=0;
			else if(this->nodeSelectionMode==ContactNodeON)
				this->numContact=0;
			else if(this->nodeSelectionMode==GroundNodeON)
				this->num_grounded=0;
			else if(this->nodeSelectionMode==ElectrodeNodeON)
				this->numElectrode=0;
		}

		for(int i=0;i<this->num_node;i++){
			for(int j=0;j<4;j++){
				selectNormal=(selectCoord[(j+1)%4]-selectCoord[j])%(selectCoord[j+4]-selectCoord[j]);
				selectNormal/=selectNormal.abs();

				if(this->isAdjustedScale)
					isPoWiTr=PointLineFace::is_view_nodeWithinTriangle(10*(this->node[i].coord+this->worldCenter)/this->worldSize,selectNormal,selectCoord[j]);
				else 
					isPoWiTr=PointLineFace::is_view_nodeWithinTriangle((this->node[i].coord+this->worldCenter)/10.0,selectNormal,selectCoord[j]);

				if(!isPoWiTr){
					if(j==3&&this->nodeSelectionMode==FixedNodeON){
						this->node[i].isNodeFixed=true;
						this->node[i].state=FixedNodeON;
						numFixed++;
					}
					else if(j==3&&this->nodeSelectionMode==ContactNodeON){
						this->node[i].isNodeContact=true;
						this->node[i].state=ContactNodeON;
						numContact++;
					}
					else if(j==3&&this->nodeSelectionMode==GroundNodeON){
						this->node[i].isNodeGround=true;
						this->node[i].state=GroundNodeON;
						this->node[i].electrodeIndex=currentIndex;
						num_grounded++;
					}
					else if(j==3&&this->nodeSelectionMode==ElectrodeNodeON){
						this->node[i].isNodeElectrode=true;
						this->node[i].state=ElectrodeNodeON;
						this->node[i].electrodeIndex=currentIndex;
						numElectrode++;
					}
				}
				else{
					if(_button==GLUT_LEFT_BUTTON){
						if(!this->node[i].isNodeFixed&&!this->node[i].isNodeContact&&!this->node[i].isNodeGround&&!this->node[i].isNodeElectrode)
							this->node[i].state=NONE;
						if(this->nodeSelectionMode==FixedNodeON)
							this->node[i].isNodeFixed=false;
						else if(this->nodeSelectionMode==ContactNodeON)
							this->node[i].isNodeContact=false;
						else if(this->nodeSelectionMode==GroundNodeON)
							this->node[i].isNodeGround=false;
						else if(this->nodeSelectionMode==ElectrodeNodeON)
							this->node[i].isNodeElectrode=false;
					}
					break;
				}
			}
		}
		if(_button==GLUT_LEFT_BUTTON){
				currentIndex=0;
		}else{
				currentIndex++;
		}
	}
*/
	double FEME3D::getPotentialAt(Vector3d _position){
		double result=0;
		if(this->is_matrix_available){
			Matrixd t_position = Matrixd(3,4);
			for(int i=0;i<this->num_elem;i++){
				int j=0;
				while(j<4){
					//どの四面体の内側かを判定する
					for(int k=0;k<3;k++)
						for(int l=0;l<4;l++)
							t_position.X[3*l+k]=elem[i].position.X[3*elem[i].order_node[4*j+l]+k];
					if(elem[i].normal[j].x*(_position.x-t_position.X[3*0+0])
						+elem[i].normal[j].y*(_position.y-t_position.X[3*0+1])
						+elem[i].normal[j].z*(_position.z-t_position.X[3*0+2])<0){//面の内外判定
						break;
					}
					else if(j==3){
						result=elem[i].getLocalPotentialAt(_position);
						return result;
					}
					j++;
				}
			}
		}
		return 0;
	}

	/*
	void FEME3D::setGroundElectrodeList(bool _isGroundOnly)
	{
		this->isElecMatrixAvailable=false;
		this->num_grounded=0;
		if(!_isGroundOnly)
			this->numElectrode=0;
		for(int i=0;i<this->num_node;i++){
			if(!_isGroundOnly)
				this->node[i].isNodeElectrode=false;
			this->node[i].isNodeGround=false;
		}
		for(int i=0;i<this->num_node;i++){
			if(!_isGroundOnly&&this->node[i].electrodeIndex==this->electIndexList[this->curElectIndex]){
				this->node[i].state=SensingNodeON;
				this->node[i].isNodeElectrode=true;
				this->numElectrode++;
			}
			else if(this->node[i].electrodeIndex==this->groundIndexList[this->curElectIndex]){
				this->node[i].state=SensingNodeON;
				this->node[i].isNodeGround=true;
				this->num_grounded++;
			}
			else if(this->node[i].state!=NONE){
				this->node[i].state=SensingNodeON;
			}
		}
		if(!_isGroundOnly){
			this->curElectIndex++;
			if(this->curElectIndex==this->electIndexList.size())
				this->curElectIndex=0;
		}
		this->isElecMatrixAvailable=true;
	}

	void FEME3D::selectGroundElectrodesAt(int _index)
	{
		if(this->isBCSLoaded){
			this->isElecMatrixAvailable=false;
			this->num_grounded=0;
			for(int i=0;i<this->num_node;i++)
				this->node[i].isNodeGround=false;
			for(int i=0;i<this->num_node;i++){
				if(this->node[i].electrodeIndex==_index){
					this->node[i].state=SensingNodeON;
					this->node[i].isNodeGround=true;
					this->num_grounded++;
				}
			}
			this->isElecMatrixAvailable=true;
		}
	}
	void FEME3D::selectActiveElectrodesAt(int _index)
	{
		if(this->isBCSLoaded){
			this->isElecMatrixAvailable=false;
			this->numElectrode=0;
			for(int i=0;i<this->num_node;i++)
				this->node[i].isNodeElectrode=false;
			for(int i=0;i<this->num_node;i++){
				if(this->node[i].electrodeIndex==_index){
					this->node[i].state=SensingNodeON;
					this->node[i].isNodeElectrode=true;
					this->numElectrode++;
				}
			}
			this->isElecMatrixAvailable=true;
		}
	}

	bool FEME3D::nextInputAs(double _voltage)
	{
		this->numElectrode=0;
		this->node[curElectIndex].isNodeElectrode=false;
		do{
			this->curElectIndex++;
			if(this->curElectIndex==this->num_node){
				this->curElectIndex=0;
				return false;
			}
		}while(this->node[curElectIndex].state==SensingNodeON&&this->electIndexList[0]!=this->node[curElectIndex].electrodeIndex);
		enodeIndexList.push_back(curElectIndex);
		this->numElectrode=1;
		this->node[curElectIndex].isNodeElectrode=true;
		this->potential.X[curElectIndex]=_voltage;
		return true;
	}
	void FEME3D::evaluateInput()
	{
		double mf=0;
		double vf=0;
		for(int i=0;i<this->mvIndexList.size();i++){
			mf+=this->potential.X[this->mvIndexList[i]]/this->mvoltage[i];
		}
		mf/=this->mvIndexList.size();
		for(int i=0;i<this->mvIndexList.size();i++){
			vf+=pow((this->potential.X[this->mvIndexList[i]]/this->mvoltage[i])-mf,2);
		}
		vf/=this->mvIndexList.size();
		meanf.push_back(mf);
		error.push_back(vf);
		cout<<this->curElectIndex<<":"<<mf<<","<<vf<<endl;
	}
	double FEME3D::getErrorAt(int _condition, int _index, int _numVoltage, double *_voltage)
	{
		double mf=0;
		double vf=0;
		double maxSimVoltage=0;
		double maxSenVoltage=0;
		double minSimVoltage=100;
		double minSenVoltage=100;
		int indexSim;
		int indexSen;
		for(int i=0;i<this->numClip;i++){
			if(maxSimVoltage<=simulatedPotential[numClip*(numCandidate*_condition+_index)+i]){
				maxSimVoltage=simulatedPotential[numClip*(numCandidate*_condition+_index)+i];
				indexSim=i;
			}
			if(maxSenVoltage<=_voltage[i]){
				maxSenVoltage=_voltage[i];
				indexSen=i;
			}
			if(i!=_condition){
				if(minSimVoltage>=simulatedPotential[numClip*(numCandidate*_condition+_index)+i]){
					minSimVoltage=simulatedPotential[numClip*(numCandidate*_condition+_index)+i];
				}
				if(minSenVoltage>=_voltage[i]){
					minSenVoltage=_voltage[i];
				}
			}
		}
		for(int i=0;i<this->numClip;i++){
			if(i!=_condition)
				mf+=pow((simulatedPotential[numClip*(numCandidate*_condition+_index)+i]-minSimVoltage)/(maxSimVoltage-minSimVoltage)
				-(_voltage[i]-minSenVoltage)/(maxSenVoltage-minSenVoltage),2.0);
		}
		mf/=numClip-1;
		if(_condition==0){
			meanf.push_back(mf);
			error.push_back(mf);
		}else{
			meanf[_index]+=mf;
			error[_index]+=mf;
		}
		return mf;
	}
	int FEME3D::getOptimizedInputFor(int _condition, int _numVoltage, double *_voltage)
	{
		if(this->numCandidate>0 && this->inputIndexList.size()==this->numCandidate){
			if(_condition==0){
				meanf.clear();
				error.clear();
			}
			if(_condition==0||meanf.size()>0){
				for(int i=0;i<this->numCandidate;i++){
						this->getErrorAt(_condition,i,_numVoltage,_voltage);
				}
			}
			if(error.size()>0&&_condition==this->numClip-1){
				//cout<<endl;
				//最小要素を取得する
				std::vector<double>::iterator iter = std::min_element(error.begin(), error.end());
				size_t index = std::distance(error.begin(), iter);
		//		std::cout << "min input node is " << index << std::endl;

				std::vector<std::pair<double, int> > list;
				for (int i = 0; i < error.size(); i++) {
					//平面裏半分は候補に入れない
					if(this->node[inputIndexList[i]].coord.y<0)
						list.push_back(std::pair<double, int>(error[i], i));
				}
				std::sort(list.begin(), list.end());
				Vector3D tposition=Vector3D(0,0,0);
				int tNumTop=10;
				for(int i=0;i<tNumTop;i++){//TODO:ひとまずデータ点を10確保
					tposition+=2*(tNumTop-i)*node[inputIndexList[list[i].second]].coord/((tNumTop+1)*tNumTop);
				}
				touchPos=tposition;
				if(this->isAdjustedScale)
					touchPos*=10.0/this->worldSize;
				else
					touchPos*=1/10.0;


				//ひとまず簡易的にindexがどのノードかを調べる
				for(int j=0;j<this->num_node;j++)
					this->node[j].isNodeElectrode=false;
				for(int j=0;j<inputIndexList.size();j++){//TODO:いずれはノードの削減を行う必要がある（表面を抽出した結果や裏面などを考慮する）
					if(j==list[0].second){//index?
						this->numElectrode=1;
						this->node[inputIndexList[j]].isNodeElectrode=true;
						//cout<<"Error:"<<list[0].first<<endl;
						double val_in=0;
						if(list[0].first<1.8){
							for(int k=0;k<16;k++)val_in=max(val_in,_voltage[k]);
							this->setPotential(15*val_in);
						}
						else
							this->setPotential(0);
						this->Potential2Potential();
						this->curInputIndex=inputIndexList[j];
						return this->curInputIndex;
					}
				}
			}
		}
		return -1;
	}
	double FEME3D::getMaximumSimulatedPotentialFor(int _condition)
	{
		double maxVs=0;
		for(int i=0;i<this->numClip;i++){
			if(i!=_condition){
				if(maxVs<simulatedPotential[numClip*(numCandidate*_condition+curInputIndex)+i]){
					maxVs=simulatedPotential[numClip*(numCandidate*_condition+curInputIndex)+i];
				}
			}
		}
		return maxVs;
	}
	double FEME3D::getMinimumSimulatedPotentialFor(int _condition)
	{
		double minVs=100;
		for(int i=0;i<this->numClip;i++){
			if(i!=_condition){
				if(minVs>simulatedPotential[numClip*(numCandidate*_condition+curInputIndex)+i]){
					minVs=simulatedPotential[numClip*(numCandidate*_condition+curInputIndex)+i];
				}
			}
		}
		return minVs;
	}

	void FEME3D::getOptimizedInput(double _initialValue)
	{
		std::vector<double>::iterator iter = std::min_element(error.begin(), error.end());
		size_t index = std::distance(error.begin(), iter);
		std::cout << "min input node is " << enodeIndexList[index] << std::endl;

		this->numElectrode=0;
		for(int i=0;i<this->num_node;i++)
			this->node[i].isNodeElectrode=false;
		this->node[ enodeIndexList[index] ].isNodeElectrode=true;
		this->numElectrode=1;
		this->potential.X[ enodeIndexList[index] ]=_initialValue/meanf[index];
		std::cout << "Potential condition is " << _initialValue/meanf[index] << std::endl;

		meanf.clear();
		error.clear();
		enodeIndexList.clear();
	}
*/

	void FEME3D::render(void)
	{
		if(!this->is_loaded)
			return;

		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);
		glPushMatrix();

		if(this->is_auto_scale)
			glScaled(10.0/this->size.abs(), 10.0/this->size.abs(), 10.0/this->size.abs());
		else//モデルスケールによる調整なしの描画ではcm単位で描画する
			glScaled(1/10.0, 1/10.0, 1/10.0);
		glTranslated(this->center.x, this->center.y, this->center.z);

		//Lineの描画
		if(this->is_view_line && this->is_view_line_extracted){
			glDisable(GL_LIGHTING);
			glColor3d(0.1, 0.1, 0.1);
			glLineWidth(0.5);
			for(int i=0;i<num_line;i++){
				//接触点を含む線は描画しない
				//if(!node[this->line[i].nodeIndex[0]].isNodeContact&&!node[this->line[i].nodeIndex[1]].isNodeContact){
					glBegin(GL_LINE_STRIP);
					glVertex3dv(this->line[i].vertex[0].X);	
					glVertex3dv(this->line[i].vertex[1].X);	
					glEnd();
				//}
			}
			glEnable(GL_LIGHTING);
		}

		if(this->is_view_facet && this->is_view_facet_extracted){
			//if(this->is_viewMisesStress)
			glDisable(GL_LIGHTING);
			for(int i=0;i<this->num_facet;i++)
			{

	//			if(this->surf[i].isCollidSurface){//接触面の表示
	//				this->faceMaterial.setAmbient(Vector3D(0.8,0,0.8));
	//				glColor3dv(Vector3D(1, 0, 1).X);
	//			}
				glBegin(GL_POLYGON);
				if(!this->is_view_smooth)
					glNormal3dv(facet[i].normal[0].X);
				for(int j=0;j<3;j++){
	//TODO::変位等の可視化
	//				if(this->is_viewMisesStress){
	//					this->map.glColorMap(this->surf[i].deform[j].abs());
	//					this->map.glColorMap(this->elem[this->surf[i].elemIndex].strain.X[3*j+0]*this->elem[this->surf[i].elemIndex].strain.X[3*j+0]+this->elem[this->surf[i].elemIndex].strain.X[3*j+1]*this->elem[this->surf[i].elemIndex].strain.X[3*j+1]+this->elem[this->surf[i].elemIndex].strain.X[3*j+2]*this->elem[this->surf[i].elemIndex].strain.X[3*j+2]);
	//					this->map.glColorMap(this->node[this->surf[i].nodeIndex[j]].deformation.abs());
	//				}
	//頂点法線
					if(this->is_view_smooth)
						glNormal3dv(this->facet[i].normal[j].X);

					color_map.glColorMap(node[facet[i].index_node[j]].potential, MAP_JET);
					glVertex3dv(this->facet[i].vertex[j].X);
				}
				glEnd();
			}
			glEnable(GL_LIGHTING);

		}
/*
		if(this->is_viewVectorField){
			glDisable(GL_LIGHTING);
			double ang;
			for(int i=0;i<this->numTetrahedra;i++){
				glPushMatrix();
				glTranslated(elem[i].center.x,elem[i].center.y,elem[i].center.z);
				if(elem[i].absCurrentDensity>0){
					if(elem[i].currentDensity.X[0]>=0)
						ang=acos(elem[i].currentDensity.X[2]/elem[i].absCurrentDensity)/M_PI*180.0;
					else	ang=M_PI*180.0-acos(elem[i].currentDensity.X[2]/elem[i].absCurrentDensity)/M_PI*180.0;
				}
				this->map.glColorMap(elem[i].absCurrentDensity);
				glRotated( ang, -elem[i].currentDensity.X[1], elem[i].currentDensity.X[0], 0.0);
				glutSolidArrow(0.5,0,0.05,0.5,8);
				glPopMatrix();
			}
			glEnable(GL_LIGHTING);
		}
*/
		if(this->is_view_node){
			glPointSize(5);
			glDisable(GL_LIGHTING );
			glBegin(GL_POINTS);
			for(int i=0;i<num_node;i++){
				if(node[i].state == GROUNDED ){glColor3d(1,0,0);}
				else if(node[i].state == ELECTRODE ){glColor3d(1,1,0);}
				else {glColor3d(0,1,1);}
				glVertex3dv(this->node[i].vertex.X);
			}
			glEnd();
			glEnable(GL_LIGHTING );
		}
		glPopMatrix();
		glDisable(GL_DEPTH_TEST );
		glDisable(GL_ALPHA_TEST);
	}
}