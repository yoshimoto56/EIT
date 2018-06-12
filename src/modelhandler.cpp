#include "modelhandler.h"

namespace EITS
{

	bool convertSTL2VOX(StlMesh *_stl, VoxModel *_vox, int _num)
	{
		std::cout<<"Converting STL data to VOX data ...";

		_vox->setIsLoaded(false);
		_vox->deleteVoxel();
		_vox->setVoxelSameNum(_num);
		_vox->newVoxel();
		_vox->setScale(_stl->getScale());
		_vox->setCenter(_stl->getCenter());
		_vox->calCoordinateSetting();
		_vox->setIsLoaded(true);

		for(int i=0;i<_stl->getNumFacet();i++){
			testIntersection(_stl->getFacetPointer(i), _vox, 255);
		}

		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	bool convertOBJ2VOX(ObjMesh *_obj, VoxModel *_vox, int _num)
	{
		std::cout<<"Converting OBJ data to VOX data ...";

		_vox->setIsLoaded(false);
		_vox->deleteVoxel();
		_vox->setVoxelSameNum(_num);
		_vox->newVoxel();
		_vox->setScale(_obj->getScale());
		_vox->setCenter(_obj->getCenter());
		_vox->calCoordinateSetting();
		_vox->setIsLoaded(true);

		for(int i=0;i<_obj->getNumFacet();i++){
			Vector3d t_color;
			t_color.x=_obj->getMaterial(_obj->getFacet(i).index_material).getDiffuse().x;
			t_color.y=_obj->getMaterial(_obj->getFacet(i).index_material).getDiffuse().y;
			t_color.z=_obj->getMaterial(_obj->getFacet(i).index_material).getDiffuse().z;
			testIntersection(_obj->getFacetPointer(i),_vox,t_color,1);
		}

		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	bool convertOBJ2STL(ObjMesh *_obj, StlMesh *_stl)
	{
		std::cout<<"Converting OBJ data to STL data ...";

		_stl->deleteMesh();
		_stl->deleteMaterial();
		_stl->setNumNode(_obj->getNumNode());
		_stl->setNumFacet(_obj->getNumFacet());
		_stl->setNumNormal(_obj->getNumNormal());
		_stl->setNumMaterial(_obj->getNumMaterial());
		_stl->setNumLine(_obj->getNumLine());
		_stl->newMesh();
		_stl->newMaterial();
		for(int i = 0; i < _obj->getNumFacet(); i++){
			_stl->setFacet(i, _obj->getFacet(i));
		}
		for(int i = 0; i < _obj->getNumNormal(); i++){
			_stl->setNormal(i, _obj->getNormal(i));
		}
		for(int i = 0; i < _obj->getNumNode(); i++){
			_stl->setVertex(i, _obj->getVertex(i));
		}
		_stl->setup();

		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	bool convertSTL2OBJ(StlMesh *_stl, ObjMesh *_obj)
	{
		std::cout<<"Converting STL data to OBJ data ...";

		_obj->deleteMesh();
		_obj->deleteMaterial();
		_obj->setNumNode(_stl->getNumNode());
		_obj->setNumFacet(_stl->getNumFacet());
		_obj->setNumNormal(_stl->getNumNormal());
		_obj->setNumMaterial(_stl->getNumMaterial());
		_obj->setNumLine(_stl->getNumLine());
		_obj->newMesh();
		_obj->newMaterial();
		for(int i = 0; i < _stl->getNumFacet(); i++){
			_obj->setFacet(i, _stl->getFacet(i));
			_obj->setNormal(i, _stl->getNormal(i));
		}
		for(int i = 0; i < _stl->getNumNode(); i++){
			_obj->setVertex(i, _stl->getVertex(i));
		}
		_obj->setup();

		std::cout<<"[OK]"<<std::endl;
		return true;
	}

	bool convertSTL2FEM(StlMesh *_stl, VolumeMesh *_fem, double val)
	{
		tetgenio in, out, addin, bgmin;
		tetgenbehavior test;
		char **argv;

		argv = (char**)malloc(sizeof(char*) * 20);
		argv[0] = (char *)malloc(256 * sizeof(char));
		argv[1] = (char *)malloc(256 * sizeof(char));
		argv[2] = (char *)malloc(256 * sizeof(char));

		int argc;
		_stl->save("temp.stl", STL_ASCII);
		sprintf_s(argv[0], 256, "tetGen");
		sprintf_s(argv[1], 256, "-pq%f", val);
		sprintf_s(argv[2], 256, "temp.stl");
		argc = 3;
		if (!test.parse_commandline(argc, argv)) {
			return false;
		}
		if (test.refine) {
			if (!in.load_tetmesh(test.infilename, tetgenbehavior::STL)) {
				terminatetetgen(NULL, 3);
			}
		}
		else {
			if (!in.load_plc(test.infilename, (int)test.object)) {
				terminatetetgen(NULL, 3);
			}
		}
		if (test.insertaddpoints) {
			if (!addin.load_node(test.addinfilename)) {
				addin.numberofpoints = 0l;
			}
		}
		if (test.metric) {
			if (!bgmin.load_tetmesh(test.bgmeshfilename, tetgenbehavior::STL)) {
				bgmin.numberoftetrahedra = 0l;
			}
		}
		//To Do オプションを設定可能にする
		//メッシュサイズ
		//領域分類
		if (bgmin.numberoftetrahedra > 0l) {
			tetrahedralize(&test, &in, &out, &addin, &bgmin);
		}
		else {
			tetrahedralize(&test, &in, &out, &addin, NULL);
		}
		//Delete temp file
		remove("temp.stl");
		if (!convertTetGen2FEM(&out, _fem)){
			return false;
		}
		free(argv);
		_fem->setup();
		return true;
	}
	bool convertSTL2TetGen(StlMesh *_stl, tetgenio *_in)
	{
		tetgenio::facet *f;
		tetgenio::polygon *p;

		if (_stl->getNumNode()<1)
			return false;
		_in->firstnumber = 1;
		_in->numberofpoints = _stl->getNumNode();
		_in->pointlist = new REAL[_in->numberofpoints * 3];
		for (int i = 0; i<_stl->getNumNode(); i++) {
			for (int j = 0; j<3; j++) {
				_in->pointlist[3 * i + j] = _stl->getVertex(i).X[j];
			}
		}
		_in->numberoffacets = _stl->getNumFacet();
		_in->facetlist = new tetgenio::facet[_in->numberoffacets];
		_in->facetmarkerlist = new int[_in->numberoffacets];
		for (int i = 0; i<_stl->getNumFacet(); i++) {
			f = &_in->facetlist[i];
			f->numberofpolygons = 1;
			f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
			f->numberofholes = 0;
			f->holelist = NULL;
			p = &f->polygonlist[0];
			p->numberofvertices = 3;
			p->vertexlist = new int[p->numberofvertices];
			for (int j = 0; j>3; j++) {
				p->vertexlist[j] = _stl->getFacet(i).index_node[j];
			}
			_in->facetmarkerlist[i] = 0;//ここが不明
		}
		return true;
	}
	bool convertTetGen2FEM(tetgenio *_out, VolumeMesh *_fem)
	{
		_fem->clear();
		_fem->deleteMesh();
		_fem->setNumNode(_out->numberofpoints);
		_fem->setNumElem(_out->numberoftetrahedra);
		_fem->setNumLine(4 * _out->numberoftetrahedra);
		_fem->setNumFacet(4 * _out->numberoftetrahedra);
		_fem->newMesh();
		for (int i = 0; i< _out->numberofpoints; i++) {
			for (int j = 0; j<3; j++) {
				_fem->getNodeAt(i)->vertex.X[j] = _out->pointlist[3 * i + j];
			}
		}
		for (int i = 0; i< _out->numberoftetrahedra; i++) {
			if (_out->numberofcorners != 4)
				return false;
			for (int j = 0; j<4; j++) {
				_fem->getElemAt(i)->index_node[j] = _out->tetrahedronlist[4 * i + j] - 1;
			}
			_fem->getElemAt(i)->index_material = 0;//ToDo
		}

		_fem->setIsLoaded(true);

		return true;
	}

	void testIntersection(Facet *_facet, VoxModel *_vox, int _map, int _thickness)
	{
		AABBtest isIntersectFaceBox;
		double vp[9];
		int j,k,l;
		Vector3d Wmin,Wmax;
		Vector3d Vmin,Vmax;
		for(int ivp=0;ivp<3;ivp++)
			for(int dim=0;dim<3;dim++)
				vp[3*ivp+dim]=_facet->vertex[ivp].X[dim];
		Wmin=_facet->getMinVertex();
		Wmax=_facet->getMaxVertex();
		Vmin=_vox->Tr.getModel2Voxel()*Wmin;
		Vmax=_vox->Tr.getModel2Voxel()*Wmax;
		for(j=Vmin.x-_thickness;j<Vmax.x+_thickness;j++){
			if(0<=j&&j<_vox->getWidth()){
				for(k=Vmin.y-_thickness;k<Vmax.y+_thickness;k++){
					if(0<=k&&k<_vox->getHeight()){
						for(l=Vmin.z-_thickness;l<Vmax.z+_thickness;l++){
							if(0<=l&&l<_vox->getDepth()){
								isIntersectFaceBox.setData(_vox->getSize().x/2.0,_vox->getCenter(j,k,l).X,vp);
								if(isIntersectFaceBox.triBoxOverlap()){
									_vox->setVoxelValue(j,k,l,_map);
								}
							}
						}
					}
				}
			}
		}
	}
	void testIntersection(Facet *_facet, VoxModel *_vox, Vector3d _color, int _thickness)
	{
		AABBtest isIntersectFaceBox;
		double vp[9];
		int j,k,l;
		Vector3d Wmin,Wmax;
		Vector3d Vmin,Vmax;
		for(int ivp=0;ivp<3;ivp++)
			for(int dim=0;dim<3;dim++)
				vp[3*ivp+dim]=_facet->vertex[ivp].X[dim];
		Wmin=_facet->getMinVertex();
		Wmax=_facet->getMaxVertex();
		Vmin=_vox->Tr.getModel2Voxel()*Wmin;
		Vmax=_vox->Tr.getModel2Voxel()*Wmax;
		for(j=Vmin.x-_thickness;j<Vmax.x+_thickness;j++){
			if(0<=j&&j<_vox->getWidth()){
				for(k=Vmin.y-_thickness;k<Vmax.y+_thickness;k++){
					if(0<=k&&k<_vox->getHeight()){
						for(l=Vmin.z-_thickness;l<Vmax.z+_thickness;l++){
							if(0<=l&&l<_vox->getDepth()){
								isIntersectFaceBox.setData(_vox->getSize().x/2.0,_vox->getCenter(j,k,l).X,vp);
								if(isIntersectFaceBox.triBoxOverlap()){
									_vox->setVoxelColor(j,k,l,_color);
								}
							}
						}
					}
				}
			}
		}
	}
	bool calOBJandVOXcollision(ObjMesh *_obj, VoxModel *_vox, int *_index)
	{
		int num = _obj->getNumNode();
		//ボクセルのフラグを開放する
		int index;
		double max=0;
		bool isCollide=false;
		*_index=-1;
		if(_vox->getIsLoaded()){
			for(int i=0;i<num;i++){
				//まず座標を変換、次にボクセルを参照、ボクセルにフラグを立てる
				index=_vox->getIndexOfVoxelAt(_obj->Tr.getWorld2Model()*_obj->getVertex(i));
				if(index!=-1){
					if(_vox->getColorAt(index).abs()!=0){
						isCollide=true;
					}
				}
			}
		}
		return isCollide;
	}
	/*
	double calBLADEandVOXcollision(transferMatrixd _Tw2o, BLADE *_blade, VoxModel *_vox, int *_index)
	{
		int num = _blade->getNumPc();
		int index;
		double max=0;
		bool isCollide=false;
		*_index=-1;
		if(_vox->getIsLoaded()){
			for(int i=0;i<num;i++){
				index=_vox->getIndexOfVoxelAt(_Tw2o*_blade->getPcAt(i));
				if(index!=-1){
					if(_vox->getColorAt(index).abs()!=0){
						isCollide=true;
						if(max<_vox->getMapAt(index,true,_Tw2o*_blade->getPcAt(i))){
							max=_vox->getMapAt(index,true,_Tw2o*_blade->getPcAt(i));
							*_index=i;
						}
					}
				}
			}
		}
		return max;
	}
	bool calBLADEandVOXcollision(transferMatrixd _Tw2o, BLADE *_blade, VoxModel *_vox, int *_index, Vector3d *_points, int *_numPoints)
	{
		*_numPoints=0;
		int numPoints=0;
		int num = _blade->getNumPc();
		int index=-1;
		bool isCollide=false;
		double max=0;
		*_index=-1;
		if(_vox->getIsLoaded()){
			for(int i=0;i<num;i++){
				index=_vox->getIndexOfVoxelAt(_Tw2o*_blade->getPcAt(i));
				if(index!=-1){
					if(_vox->getColorAt(index).abs()!=0){
						isCollide=true;
						if(_vox->getColorAt(index)!=INVIS_COL){
							_points[numPoints]=_vox->Tr.getModel2World()*_Tw2o*_blade->getPcAt(i);
							numPoints++;
						}
						if(max<_vox->getMapAt(index,true,_Tw2o*_blade->getPcAt(i))){
							max=_vox->getMapAt(index,true,_Tw2o*_blade->getPcAt(i));
							*_index=i;
						}
					}
				}
			}
		}
		*_numPoints=numPoints;
		return isCollide;
	}*/
	/*
	double calOBJandVOXdistance(ObjMesh *_obj, VoxModel *_vox, int *_index)
	{
		double min=INT_MAX;
		double dist=0;
		if(_vox->getIsLoaded()){
			for(int i=0;i<_vox->getNumVoxel();i++){
				if(_vox->getColorAt(i)==SURF_COL){
					dist=(_obj->Tr.getWorld2Model()*_obj->getTri()->getCenterCoord()-_vox->Tr.getWorld2Model()*_vox->getCenter(i)).abs();
					if(min>dist)
						min=dist;
				}
			}
		}
		return min;
	}
	*/
	void millingVOXbyOBJ(ObjMesh *_obj, VoxModel *_vox)
	{
		int num = _obj->getNumNode();
		//ボクセルのフラグを開放する
		int index;
		int max=0;
		bool isCollide=false;
		if(_vox->getIsLoaded()){
			for(int i=0;i<num;i++){
				//まず座標を変換、次にボクセルを参照、ボクセルにフラグを立てる
				index=_vox->getIndexOfVoxelAt(_obj->Tr.getWorld2Model()*_obj->getVertex(i));
				if(index!=-1)
					_vox->millingSurfaceWith(index);
			}
		}
	}
	void millingVOXbyPOINTS(VoxModel *_vox, Vector3d *_points, double *_depth, int _numPoints)
	{
		int index;
		if(_vox->getIsLoaded()){
			for(int i=0;i<_numPoints;i++){
				index=_vox->getIndexOfVoxelAt(_points[i],COORD_LOCAL);
				if(index!=-1){
					if(_vox->getColorAt(index)!=INVIS_COL){
						_vox->coloringSurfaceWith(index, INVIS_COL);
					}
					else _depth[i]=0;
				}
				else _depth[i]=0;
			}
		}
	}

	void spanLabelVOX2OBJ(VoxModel *_vox, ObjMesh *_obj)
	{
		Vector3f tcolor;
		for(int i=0;i<_obj->getNumNode();i++){
			tcolor.x=_vox->getColorAt(_vox->getIndexOfVoxelAt(_obj->getVertex(i))).x;
			tcolor.y=_vox->getColorAt(_vox->getIndexOfVoxelAt(_obj->getVertex(i))).y;
			tcolor.z=_vox->getColorAt(_vox->getIndexOfVoxelAt(_obj->getVertex(i))).z;
			_obj->setColorAt(i,tcolor);
			_obj->setLabelIndexAt(i,_vox->getLabelIndexAt(_vox->getIndexOfVoxelAt(_obj->getVertex(i))));
		}
	}
	void setLabelOBJNearest(ObjMesh *_obj)
	{
		double mdist;
		for(int i=0;i<_obj->getNumNode();i++){
			mdist=1000;
			if(_obj->getLabelIndex(i)==-1){
				for(int j=0;j<_obj->getNumNode();j++){
					if(_obj->getLabelIndex(j)!=-1){
						if(mdist>(_obj->getVertex(i)-_obj->getVertex(j)).abs()){
							mdist=(_obj->getVertex(i)-_obj->getVertex(j)).abs();
							_obj->setLabelIndexAt(i,_obj->getLabelIndex(j));
							_obj->setColorAt(i,_obj->getColor(j));
						}
					}
				}
			}
		}
	}
	/*
	void clippingVOXbyPLANE(VoxModel *_vox, PLANE *_plane, int _numPlane)
	{
		bool isInside=true;
		for(int i=0;i<_vox->getNumVoxel();i++){
			if(_vox->getColorAt(i)==SURF_COL){
				isInside=true;
				for(int j=0;j<_numPlane;j++)
					if(_plane[j].getNormal()*(_vox->getCenter(i)-_plane->getCenter())>=0)
						isInside=false;
				if(isInside)
					_vox->coloringSurfaceWith(i, INVIS_COL);
				//				_vox->setFreeSpaceAt(i);
			}
		}
	}
	*/
	transferMatrixd ICP::getOptimalRegistration(ObjMesh *_objA, ObjMesh *_objB, int _mode, double minError)
	{
		int numPoint=_objA->getNumNode();
		transferMatrixd Ta2b;
		Vector3d *pointA=new Vector3d[numPoint];
		Vector3d *pointB=new Vector3d[numPoint];
		double *error=new double[numPoint];
		double dE=1000;
		double pError=0;
		double cError=0;


		//TODO:LOOP
		while(dE>minError){

			//位置合わせを行う点群を設定
			for(int i=0;i<numPoint;i++)
				pointA[i]=_objA->Tr.getWorld2Model()*_objA->getVertex(i);

			//最近傍探索点を取得
			ICP::matching(numPoint, pointA, pointB, _objB, _mode);
			//誤差の取得
			pError=cError;
			cError=ICP::errorMetric(error, pointB, _objA);
			dE=fabs(cError-pError);
			//最適化を行うための変換行列を取得
			Ta2b=ICP::minimizing(numPoint, pointA, pointB);
			//変換行列の更新
			_objA->Tr.setWorld2Model(Ta2b*_objA->Tr.getWorld2Model());

			std::cout<<"dE:"<<dE<<std::endl;
		}
		std::cout<<"Mean Square Error:"<<cError<<std::endl;
		std::cout<<"Transfer Matrix:"<<std::endl;
		std::cout<<_objA->Tr.getWorld2Model()<<std::endl;

		delete []pointA;
		delete []pointB;
		delete []error;
		return _objA->Tr.getWorld2Model();
	}
	//最近傍点の取得
	void  ICP::matching(int _numPoint, Vector3d *_pointA, Vector3d *_pointB, ObjMesh *_objB, int _mode)
	{
		double dist=0;
		double mDist=INT_MAX;
		Vector3d normal;
		Vector3d vertex;
		for(int i=0;i<_numPoint;i++){
			if(_mode==ICP2P){//点と点の場合
				for(int j=0;j<_objB->getNumNode();j++){
				}
			}
			if(_mode==ICP2S){//点と面の場合
				mDist=INT_MAX;
				for(int j=0;j<_objB->getNumFacet();j++){
					normal=_objB->getFacetPointer(j)->normal[0];
					vertex=_objB->getFacetPointer(j)->vertex[0];
					dist=getFacePointDistance(_pointA[i],normal,vertex);
					if(mDist>dist){
						mDist=dist;
						_pointB[i]=getFacePointProjection(_pointA[i], normal, vertex);
					}
				}
			}
		}
	}
	double  ICP::errorMetric(double *_error, Vector3d *_pointB, ObjMesh *_objA)
	{
		double error=0;
		double max=0;
		double min=INT_MAX;
		for(int i=0;i<_objA->getNumNode();i++){
			_error[i]=(_pointB[i]-_objA->Tr.getWorld2Model()*_objA->getVertex(i)).abs();
			if(max<_error[i])
				max=_error[i];
			if(min>_error[i])
				min=_error[i];
			error+=_error[i]*_error[i];
		}
		for(int i=0;i<_objA->getNumNode();i++){
			if(max!=0){
				Vector3d color=HSV2RGB(Vector3d(240.0-240.0*(_error[i]-min)/(max-min),1,1));
				_objA->setColorAt(i,Vector3f(color.x, color.y, color.z));
			}
		}
		return sqrt(error)/_objA->getNumNode();
	}
	Vector3d  ICP::centering(int _numPoint, Vector3d *_point)
	{
		Vector3d center;
		for(int i=0;i<_numPoint;i++)
			center+=_point[i]/(double)_numPoint;
		return center;
	}
	Vector3d  ICP::translating(int _numPoint, Vector3d *_pointA, Vector3d *_pointB){
		Vector3d centerA=ICP::centering(_numPoint,_pointA);
		Vector3d centerB=ICP::centering(_numPoint,_pointB);
		return centerB-centerA;
	}
	transferMatrixd ICP::rotating(int _numPoint, Vector3d *_pointA, Vector3d *_pointB)
	{
		transferMatrixd Trot;
		Matrixd Xa = Matrixd(_numPoint, 3);
		Matrixd Xb = Matrixd(_numPoint, 3);
		Matrixd R = Matrixd(3, 3);
		for(int i=0;i<_numPoint;i++){
			Xa.X[_numPoint*0+i]=_pointA[i].x;
			Xa.X[_numPoint*1+i]=_pointA[i].y;
			Xa.X[_numPoint*2+i]=_pointA[i].z;

			Xb.X[_numPoint*0+i]=_pointB[i].x;
			Xb.X[_numPoint*1+i]=_pointB[i].y;
			Xb.X[_numPoint*2+i]=_pointB[i].z;
		}
		if((Xa*Xa.trn()).det()!=0){
			R=Xb*Xa.trn()*(Xa*Xa.trn()).inv();
			Vector3d nx=Vector3d(R.X[0],R.X[1],R.X[2]);
			Vector3d ny=Vector3d(R.X[3],R.X[4],R.X[5]);
			Vector3d nz=Vector3d(R.X[6],R.X[7],R.X[8]);
			Trot.setRotMatrix(nx/nx.abs(), ny/ny.abs(), nz/nz.abs());
		}
		Xa.free();
		Xb.free();
		R.free();
		return Trot;
	}
	transferMatrixd ICP::minimizing(int _numPoint, Vector3d *_pointA, Vector3d *_pointB)
	{
		transferMatrixd Ttran;
		transferMatrixd Trot;
		transferMatrixd TcenterA;
		transferMatrixd TcenterB;
		Vector3d *pointA=new Vector3d[_numPoint];
		Vector3d *pointB=new Vector3d[_numPoint];
		Ttran.setTranslate(ICP::translating(_numPoint, _pointA, _pointB));
		TcenterA.setTranslate(ICP::centering(_numPoint,_pointA));
		TcenterB.setTranslate(ICP::centering(_numPoint,_pointB));

		for(int i=0;i<_numPoint;i++){
			pointA[i]=TcenterA.inv()*_pointA[i];
			pointB[i]=TcenterB.inv()*_pointB[i];
		}
		Trot=ICP::rotating(_numPoint, pointA, pointB);

		delete []pointA;
		delete []pointB;

		return Ttran*TcenterA*Trot*TcenterA.inv();
	}

	transferMatrixd calTransMatBetween3(Vector3d *_p1, Vector3d *_p2)
	{
		Matrixd Pp(4,4);
		Matrixd Pn(4,4);
		Matrixd tTm2w;
		transferMatrixd Tdst;
		Vector3d temp1=(_p1[1]-_p1[0])%(_p1[2]-_p1[0]);
		Vector3d temp2=(_p2[1]-_p2[0])%(_p2[2]-_p2[0]);
		if(temp1.abs()!=0)temp1/=temp1.abs();
		if(temp2.abs()!=0)temp2/=temp2.abs();
		for(int j=0;j<4;j++){
			Pp.X[4*3+j]=1;
			Pn.X[4*3+j]=1;
		}
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				Pp.X[4*k+j]=_p1[j].X[k];
				Pn.X[4*k+j]=_p2[j].X[k];
			}
		}
		for(int k=0;k<3;k++){
			Pp.X[4*k+3]=temp1.X[k];
			Pn.X[4*k+3]=temp2.X[k];
		}
		tTm2w=Pn*Pp.trn()*(Pp*Pp.trn()).inv();
		memcpy(Tdst.X,tTm2w.X,sizeof(double)*4*4);
		return Tdst;
	}

	double calLongAxisFromPointCloud(int _num, Vector3d *_point, int *_index)
	{
		double dst=0;
		for(int i=0;i<_num;i++){
			for(int j=_num-1;j>=0;j--){
				if(i!=j){
					if(dst<(_point[i]-_point[j]).abs()){
						dst=(_point[i]-_point[j]).abs();
						_index[0]=i;
						_index[1]=j;
					}
				}
			}
		}
		return dst;
	}
};

/*

bool INTERACTION::intractINSTandVOX(Instrument *_inst, voxelHandler *_vox)
{
	bool isInteracted=false;
	if(!_inst->getTool()->getIsLoaded()||!_vox->getIsLoaded())return isInteracted;
	int num=_inst->getTool()->num_vertex;
	Vector3d Pworld;
	int max=0;
	for(int i=0; i<num; i++){
		Pworld=_inst->getTw2v()*_inst->getTi2w()*_inst->getTool()->vertex[i]*2;//VRスペースでは座標を2倍にする
		if(_vox->millingSurfaceWith(Pworld))
			isInteracted=true;
	}
	return isInteracted;
}


int INTERACTION::navigateINSTwithVOX(Instrument *_inst, voxelHandler *_vox)
{
	if(!_inst->getTool()->getIsLoaded()||!_vox->getIsLoaded())return 0;
	int num=_inst->getTool()->num_vertex;
	Vector3d Pworld;
	int max=0;
	int val;
	for(int i=0; i<num; i++){
//		Pworld=_inst->getTw2v()*_inst->getTi2w()*_inst->getTool()->vertex[i];//オリジナル
		Pworld=_inst->getTi2w()*_inst->getTool()->vertex[i];//仮に２倍にしている
		val=_vox->navigateWith(Pworld);
		if(max<val)
			max=val;
	}
	return max;
}



*/
