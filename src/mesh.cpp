#include "mesh.h"


namespace EITS{
	/****************Surface Mesh******************/
	SurfMesh::SurfMesh(void)
		:num_node(0),num_normal(0),num_facet(0),num_line(0),num_material(0),id(0)
		,is_view(true),is_view_line(false),is_view_node(false),is_view_facet(true)
		,is_auto_scale(true),is_loaded(false),is_tri(false),is_identity(false)
		,is_selected(false),is_alpha_blend(false),is_view_line_extracted(false)
		,is_listed(true),is_display_list(false),is_vertex_color_enabled(false),
		is_view_label(false)
	{
		this->newMesh();
		this->newMaterial();
	}


	SurfMesh::~SurfMesh(void)
	{
		this->deleteMesh();
		this->deleteMaterial();
	}

	void SurfMesh::deleteMesh()
	{
		this->is_loaded=false;
		delete []vertex;
		delete []facet;
		delete []normal;
		delete []label_index;
		delete []color;
		delete []line;
		delete []is_selected;
	}

	void SurfMesh::newMesh()
	{
		this->vertex=new Vector3d[this->num_node+2];
		this->color=new Vector3f[this->num_node+2];
		this->label_index=new int[this->num_node+2];
		this->facet=new Facet[this->num_facet+2];
		this->normal=new Vector3d[this->num_normal+2];
		this->is_selected=new bool[this->num_node+2];
		this->line=new Line[this->num_line+2];

		for(int i=0;i<num_node+2;i++){
			is_selected[i]=false;
			label_index[i]=-1;
		}
	}

	void SurfMesh::deleteMaterial()
	{
		delete []material;
	}

	void SurfMesh::newMaterial()
	{
		this->material=new Material[this->num_material+2];
		this->material[0].setDiffuse(Vector4f(0.5,0.5,0.5,0.5));
		this->material[0].setAmbient(Vector4f(0.0,0.0,0.0,0.5));
	}

	void SurfMesh::clearMesh()
	{
		this->is_listed=false;
		this->is_loaded=false;
		this->is_view_line_extracted=false;
		this->is_vertex_color_enabled=false;
		this->is_display_list=false;
		num_node = 0;
		num_normal = 0;
		num_facet = 0;
		num_material = 0;
		num_line = 0;
		this->scale=0;
		this->size=Vector3d(1,1,1);
		this->center=Vector3d(0,0,0);
	}
	void SurfMesh::clearView()
	{
		this->is_view=true;
		this->is_view_facet=true;
		this->is_view_line=true;
		this->is_view_node=true;
		this->is_view_label=true;
	}
	void SurfMesh::addVertex(Vector3d &_vertex)
	{
		num_node++;
		vertex = (Vector3d*)realloc(vertex, num_node*sizeof(Vector3d));
		vertex[num_node-1] = _vertex;
		is_selected = (bool*)realloc(is_selected, num_node*sizeof(bool));
		is_selected[num_node-1]=false;
	}
	void SurfMesh::addNormal(Vector3d &_normal)
	{
		num_normal++;
		normal = (Vector3d*)realloc(normal, num_normal*sizeof(Vector3d));
		normal[num_normal-1] = _normal;
	}
	void SurfMesh::addFacet(Facet &_facet)
	{
		num_facet++;
		facet = (Facet*)realloc(facet, num_facet*sizeof(Facet));
		facet[num_facet-1] = _facet;
	}
	void SurfMesh::addMaterial(Material &_material)
	{
		num_material++;
		material = (Material*)realloc(material, num_material*sizeof(Material));
		material[num_material-1] = _material;
	}
	std::ostream& SurfMesh::getInfo(std::ostream &stream)
	{
		stream<<"Number of vertices：" << num_node << std::endl;
		stream<<"Number of normals：" << num_normal << std::endl;
		stream<<"Number of facets：" << num_facet << std::endl;
		stream<<"Number of materials：" << num_material << std::endl;
		return stream;
	}

	void SurfMesh::calFacetNormal()
	{
		Vector3d temp[2];
		delete []normal;
		this->normal=new Vector3d[this->num_facet];
		for(int i=0;i<num_facet;i++){
			this->facet[i].normal_type=NORMAL_FACET;
			this->facet[i].num_normal=1;
			for(int j=1;j<3;j++)
				temp[j-1]=vertex[facet[i].index_node[j]]-vertex[facet[i].index_node[0]];
			for(int j=0;j<3;j++){
				this->normal[i].X[j] =temp[1].X[(j+1)%3]*temp[0].X[(j+2)%3]
						-temp[1].X[(j+2)%3]*temp[0].X[(j+1)%3];
			}
			if(this->normal[i].abs()!=0)this->normal[i]/=this->normal[i].abs();
			this->facet[i].index_normal[0]=i;
			this->facet[i].normal[0]=this->normal[i];
		}
	}
	void SurfMesh::calFacetVertex()
	{
		for(int i=0;i<this->num_facet;i++){
			for(int j=0;j<this->facet[i].num_node;j++){
				this->facet[i].vertex[j]=vertex[this->facet[i].index_node[j]];
			}
		}
	}
	void SurfMesh::calLine()
	{
		if(this->num_line>1)
			delete []this->line;
		this->num_line=0;
		this->line= new Line[this->num_facet*4];//四角形までしか想定していない

		for(int i=0;i<this->num_facet;i++){
			for(int j=0;j<facet[i].num_node;j++){
				facet[i].line[j].index_node[0]=facet[i].index_node[j];
				facet[i].line[j].index_node[1]=facet[i].index_node[(j+1)%facet[i].num_node];
				facet[i].line[j].vertex[0]=facet[i].vertex[j];
				facet[i].line[j].vertex[1]=facet[i].vertex[(j+1)%facet[i].num_node];
				if(!isSharedLine(facet[i].line[j],line,num_line)){
					this->line[this->num_line]=facet[i].line[j];
					this->num_line++;
				}
			}
		}
		this->is_view_line_extracted=true;
	}

	void SurfMesh::calVertex()
	{
		int count=1;
		this->vertex[0]=this->facet[0].vertex[0];
		for(int i=0;i<this->num_facet;i++){
			for(int j=0;j<3;j++){
				for(int k=0;k<count;k++){
					if(this->facet[i].vertex[j]==this->vertex[k]){
						this->facet[i].index_node[j]=k;
						k=count;
					}
					else if(k==count-1){
						this->vertex[count]=this->facet[i].vertex[j];
						this->facet[i].index_node[j]=count;
						count++;
					}
				}
			}
		}
		this->num_node=count;
		Vector3d *tempVertex=new Vector3d[this->num_node];
		memcpy(tempVertex,vertex,sizeof(Vector3d)*this->num_node);
		delete []vertex;
		this->vertex=tempVertex;
	}
	void SurfMesh::calCenter()
	{
		Vector3d min(INT_MAX,INT_MAX,INT_MAX);
		Vector3d max(-INT_MAX,-INT_MAX,-INT_MAX);
		for(int i=0;i<this->num_node;i++){
			for(int k=0;k<3;k++){
				if(min.X[k]>this->vertex[i].X[k])
					min.X[k]=this->vertex[i].X[k];
				if(max.X[k]<this->vertex[i].X[k])
					max.X[k]=this->vertex[i].X[k];
			}
		}
		this->size.x=fabs(max.x-min.x);
		this->size.y=fabs(max.y-min.y);
		this->size.z=fabs(max.z-min.z);
		this->center=(max+min)/2.0;
	}
	void SurfMesh::calScale()
	{
		double max=-INT_MAX;
		for(int i=0;i<3;i++)
			if(max<size.X[i])max=size.X[i];
		scale=max;
	}
	void SurfMesh::calTransform()
	{
		transferMatrixd Tscale;
		transferMatrixd Ttrans;
		Tscale.setScale(10/this->scale);
		Ttrans.setTranslate(Vector3d(-this->center.x, -this->center.y, -this->center.z));
		Tr.setModel2Uni(Tscale*Ttrans);
	}

	//ファイルフォーマットごとにオーバーロードする
	bool SurfMesh::load(const char *_filename)
	{
		return true;
	}

	bool SurfMesh::save(const char *_filename)
	{
		return true;
	}


	void SurfMesh::render()
	{
		if(!this->is_loaded||!this->is_view)
			return;
		glDisable(GL_CULL_FACE); 
		if(this->is_alpha_blend){
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
			glEnable(GL_ALPHA_TEST);
		}
		glEnable(GL_DEPTH_TEST);

		glPushMatrix();

		if(is_auto_scale){
			glMultMatrixd(Tr.getModel2Uni().getTr4GL());
		}
		if(this->is_display_list){
			if(this->is_listed)
				glCallList(id);
			else this->makeDisplayList(0);
		}else renderList();

		glPopMatrix();

		glDisable(GL_BLEND);
		glDisable(GL_ALPHA_TEST);
	}

	void SurfMesh::renderList()
	{
		Vector3d p1,p2;
		int pre_mat = -1, cur_mat = 0;
		if(this->is_view_facet){
			glEnable(GL_LIGHTING);
			if(this->is_view_label)
				glDisable(GL_LIGHTING);
			for ( int i=0; i<num_facet; i++ ){
				if(facet[i].is_enabled){
					cur_mat = facet[i].index_material;
					if ( pre_mat != cur_mat ){
						material[cur_mat].set();
						pre_mat = cur_mat;
					}
					if(facet[i].is_selected){
						material_s.set();
						pre_mat=-1;
					}
					if(this->facet[i].normal_type==NORMAL_FACET)
						glNormal3dv(this->normal[this->facet[i].index_normal[0]].X);
					glBegin(facet[i].type);
					for ( int j=0; j<facet[i].num_node; j++ ){
						if(this->facet[i].normal_type==NORMAL_POINT)
							glNormal3dv(this->normal[this->facet[i].index_normal[j]].X);
						glColor3fv(color[facet[i].index_node[j]].X);
						p1=Tr.getModel2World()*vertex[facet[i].index_node[j]];
						glVertex3dv(p1.X);
					}
					glEnd();
				}
			}
		}
		if(this->is_view_line){
			glDisable(GL_LIGHTING);
			glColor3d(0.1, 0.1, 0.1);
			glLineWidth(0.5);
			if(this->is_view_line_extracted){
				for(int i=0;i<num_line;i++)
					this->line[i].render();
			}else{
				for(int i=0;i<num_facet;i++){
					glBegin(GL_LINE_LOOP);
					for ( int j=0; j<facet[i].num_node; j++ )
						glVertex3dv(vertex[facet[i].index_node[j]].X);
					glEnd();
				}
			}
			glEnable(GL_LIGHTING);
		}
		if(this->is_view_node){
			glDisable(GL_LIGHTING);
			glPointSize(3);
			for(int i=0;i<num_node;i++){
				glPushMatrix();
				glMultMatrixd(Tr.getModel2World().getTr4GL());
				if(this->is_selected[i])
					glColor3d(1,0,1);
				if(this->is_view_label){
					glColor3fv(color[i].X);
				}
				else glColor3d(0,1,1);
				glBegin(GL_POINTS);
				glVertex3dv(vertex[i].X);
				glEnd();
				glPopMatrix();
			}
			glEnable(GL_LIGHTING);
		}
	}

	GLuint SurfMesh::makeDisplayList(int _id)
	{
		id = glGenLists(1);
		glNewList(id, GL_COMPILE);
		renderList();
		glEndList();
		this->is_listed=true;
		return id;
	}

	int SurfMesh::getVertexIndexNear(Vector3d _pos)
	{
		Vector3d error;
		Vector3d minError=Vector3d(1000,1000,1000);
		int min_error_index=-1;
		for(int i=0;i<this->num_node;i++){
			if(!this->is_identity){
				if(is_auto_scale)
					error=(Tr.getModel2Uni()*Tr.getWorld2Model()*this->vertex[i]-_pos);
				else
					error=(Tr.getWorld2Model()*this->vertex[i]-_pos);
			}
			else error=this->vertex[i]-_pos;
			if(error.abs()<minError.abs()&&error.abs()<10){
				minError=error;
				min_error_index=i;
			}
			this->is_selected[i]=false;
		}
		if(min_error_index!=-1){
			this->is_selected[min_error_index]=true;
			this->selected_vertex=this->vertex[min_error_index];
			this->selected_index=min_error_index;
		}
		return min_error_index;
	}
	int SurfMesh::getFacetIndexNear(Vector3d _pos)
	{
		double dist;
		double mDist=10;//クリックした点から面までの距離が最大でも10
		int mDistIndex=-1;
		Vector3d Pc_o;
		Vector3d Pf_o;
		if(!this->is_identity){
			if(is_auto_scale)
				Pc_o=Tr.getModel2World()*Tr.getUni2Model()*_pos;
			else
				Pc_o=Tr.getModel2World()*_pos;
		}
		else Pc_o=_pos;
		for(int i=0;i<this->num_facet;i++){
			if(facet[i].is_enabled){
				facet[i].setIsSelected(false);
				if(isProjectedPointOnFiniteFace(Pc_o, &facet[i])){
					dist=getFacePointDistance(Pc_o, facet[i].normal[0], facet[i].vertex[0]);
					if(mDist>fabs(dist)){
						mDist=fabs(dist);
						mDistIndex=i;
					}
				}
			}
		}
		if(mDistIndex!=-1)
			facet[mDistIndex].setIsSelected(true);
		return mDistIndex;
	}
	void SurfMesh::clearSelection()
	{
		for(int i=0;i<num_facet;i++){
			facet[i].is_selected=false;
			this->is_selected[i]=false;
		}
	}

	void SurfMesh::quadrize()
	{
		std::cout<<"Quadrizing surface...";
		for(int i=0;i<this->num_facet;i++)
			if(this->facet[i].is_enabled)
				for(int j=0;j<this->num_facet;j++)
					if(i!=j)
						convertTriangle2Quad(&this->facet[i],&this->facet[j]);
		this->calLine();
		std::cout<<"[OK]"<<std::endl;
	}
	void SurfMesh::centerize()
	{
		std::cout << "Centerizing mesh...";
		this->calCenter();
		for (int i = 0;i < this->num_node;i++) {
			this->vertex[i] -= this->center;
		}
		for (int i = 0;i < this->num_facet;i++) 
			for(int j=0;j<this->facet[i].num_node;j++)
				this->facet[i].vertex[j] -= this->center;
		std::cout << "[OK]" << std::endl;
	}
	/****************Surface Mesh******************/


	/****************Volume Mesh******************/
	VolumeMesh::VolumeMesh():num_elem(0),num_node(0),num_facet(0),num_line(0)
	{
		this->newMesh();
		this->clear();
		this->clearView();
	}

	VolumeMesh::~VolumeMesh()
	{
	}

	void VolumeMesh::newMesh()
	{
		this->elem = new Tetrahedra[this->num_elem + 2];
		this->node = new Node[this->num_node + 2];
		this->line = new Line[this->num_line + 2];
		this->facet = new Facet[this->num_facet + 2];
	}

	void VolumeMesh::deleteMesh()
	{
		this->is_loaded = false;
		delete []elem;
		delete []node;
		delete []line;
		delete []facet;
	}

	void VolumeMesh::clear()
	{
		this->is_loaded = false;
		this->is_auto_scale = false;
		this->is_view_line_extracted = false;
		this->is_view_facet_extracted = false;
		this->center = Vector3d(0, 0, 0);
		this->size = Vector3d(1, 1, 1);
		this->scale = 1;
	}

	void VolumeMesh::clearView()
	{
		this->is_view_line = true;
		this->is_view_node = true;
		this->is_view_facet = true;
		this->is_view_smooth = false;
	}

	bool VolumeMesh::load(const char *_filename)
	{
		std::ifstream file;
		char buf[512];
		char dummy[256];
		int num;

		file.open(_filename, std::ios::in);
		if(!file.is_open())
		{
			std::cout<<"Fail to load volume mesh "<<_filename<<std::endl;
			return false;
		}
		this->deleteMesh();
		while(file.getline(buf, sizeof(buf))){
			if(strstr(buf, "nNodes ") || strstr(buf, "nVertex ")){
				sscanf_s(buf, "%s %d", &dummy, 256, &num); 
				num_node = num;
				continue;
			}
			if(strstr(buf, "nTetrahedra ") || strstr(buf, "nTetrahedron ")){
				sscanf_s(buf, "%s %d", &dummy, 256, &num);
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

		this->setup();

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
*/		return true;
	}

	bool VolumeMesh::save(const char *_filename)
	{
		std::ofstream file;
		std::cout<<"Saving mesh data...";
		file.open(_filename, std::ios::out);
		if(!this->is_loaded || !file.is_open())
		{
			std::cout<<"[FAIL]"<<std::endl;
			return false;
		}
		file<<"nNodes "<<this->num_node<<std::endl;
		file<<"nTetrahedra "<<this->num_elem<<std::endl;

		file<<"Nodes { float[3] Coordinates } @1"<<std::endl;
		file<<"Tetrahedra { int[4] Nodes } @2"<<std::endl;
		file<<"TetrahedronData { byte Materials } @3"<<std::endl;
		file<<"# Data section follows"<<std::endl;
		file<<"@1"<<std::endl;
		for(int i=0; i< num_node; i++)
			file<<node[i].vertex.X[0]<<" "<<node[i].vertex.X[1]<<" "<<node[i].vertex.X[2]<<std::endl;
		file<<"@2"<<std::endl;
		for(int i=0; i< num_elem; i++){
			for(int j=0; j<4; j++){
				file << elem[i].index_node[j]+1<<" ";
			}
			file<<std::endl;
		}
		file<<"@3"<<std::endl;
		for(int i=0; i< num_elem; i++)
			file << elem[i].index_material<<std::endl;
		file.clear();
		file.close();
		std::cout<<"[OK]"<<std::endl;
		return true;
	}
	void VolumeMesh::calCenter()
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
	void VolumeMesh::calScale()
	{
		double max = -INT_MAX;
		for (int i = 0;i<3;i++)
			if (max<size.X[i])max = size.X[i];
		scale = max;
	}
	void VolumeMesh::calTransform()
	{
		transferMatrixd Tscale;
		transferMatrixd Ttrans;
		Tscale.setScale(10 / this->scale);
		Ttrans.setTranslate(Vector3d(-this->center.x, -this->center.y, -this->center.z));
		Tr.setModel2Uni(Tscale*Ttrans);
	}
	std::ostream& VolumeMesh::getInfo(std::ostream &stream)
	{
		stream << "Number of vertices：" << this->num_node << std::endl;
		stream << "Number of tetrahedras：" << this->num_elem << std::endl;
		stream << "Number of surfaces：" << this->num_facet << std::endl;
		stream << "Number of lines：" << this->num_line << std::endl;
		return stream;
	}
	void VolumeMesh::render(void)
	{
		if(!this->is_loaded)
			return;

		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);
		glPushMatrix();

		if (is_auto_scale) {
			glMultMatrixd(Tr.getModel2Uni().getTr4GL());
		}

		//Lineの描画
		if(this->is_view_line && this->is_view_line_extracted){
			glDisable(GL_LIGHTING);
			glColor3d(0.1, 0.1, 0.1);
			glLineWidth(0.5);
			for(int i = 0; i < num_line; i++){
				glBegin(GL_LINE_STRIP);
				glVertex3dv(this->line[i].vertex[0].X);	
				glVertex3dv(this->line[i].vertex[1].X);	
				glEnd();
			}
			glEnable(GL_LIGHTING);
		}

		if(this->is_view_facet && this->is_view_facet_extracted){
			//if(this->is_viewMisesStress)
			glDisable(GL_LIGHTING);
			glColor3d(1.0, 0.1, 0.1);
//				this->faceMaterial.enable();
//				this->faceMaterial.setAmbient(Vector3D(0,0,0));
			glBegin(GL_TRIANGLES);
			for(int i=0;i<this->num_facet;i++){
//				if(!this->is_view_smooth)
//					glNormal3dv(facet[i].normal[0].X);
					for(int j=0;j<3;j++){
						glVertex3dv(this->facet[i].vertex[j].X);
					}
			}
			glEnd();
			glEnable(GL_LIGHTING);
		}

		if(this->is_view_node){
			glPointSize(5);
			glDisable(GL_LIGHTING );
			glColor3d(0,1,1);
			glBegin(GL_POINTS);
			for(int i=0;i<num_node;i++){
				glVertex3dv(this->node[i].vertex.X);
			}
			glEnd();
			glEnable(GL_LIGHTING );
		}
		glPopMatrix();
		glDisable(GL_DEPTH_TEST );
	}

	int VolumeMesh::makeDisplayList(int _id){
		int tid;
		tid = glGenLists(_id+1);
		glNewList(tid, GL_COMPILE);
		this->render();
		glEndList();
		return tid;
	}

	void VolumeMesh::calLine()
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
				for(int k=j;k<4;k++){
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

	void VolumeMesh::calFacet()
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

	/****************Volume Mesh******************/

	bool isSharedLine(Line _new, Line *_old, int _num)
	{
		for(int i=0;i<_num;i++){
			if(_old[i].index_node[0]==_new.index_node[0]&&_old[i].index_node[1]==_new.index_node[1])
				return true;
			else if(_old[i].index_node[0]==_new.index_node[1]&&_old[i].index_node[1]==_new.index_node[0])
				return true;
		}
		return false;
	}

	bool isSharedFacet(Facet _new, Facet *_old, int _num)
	{
		int count=0;
		for(int i=0; i<_num; i++){
			if((_new.index_node[0]==_old[i].index_node[0] && _new.index_node[1] == _old[i].index_node[1] && _new.index_node[2] == _old[i].index_node[2]) ||
				(_new.index_node[0]==_old[i].index_node[1] && _new.index_node[1] == _old[i].index_node[2] && _new.index_node[2] == _old[i].index_node[0]) ||
				(_new.index_node[0]==_old[i].index_node[2] && _new.index_node[1] == _old[i].index_node[0] && _new.index_node[2] == _old[i].index_node[1]) ||
				(_new.index_node[0]==_old[i].index_node[0] && _new.index_node[1] == _old[i].index_node[2] && _new.index_node[2] == _old[i].index_node[1]) ||
				(_new.index_node[0]==_old[i].index_node[1] && _new.index_node[1] == _old[i].index_node[0] && _new.index_node[2] == _old[i].index_node[2]) ||
				(_new.index_node[0]==_old[i].index_node[2] && _new.index_node[1] == _old[i].index_node[1] && _new.index_node[2] == _old[i].index_node[0])){
					return true;
			}	
		}
		return false;
	}

	bool convertTriangle2Quad(Facet *_tri1, Facet *_tri2)
	{
		if(_tri1->type!=GL_TRIANGLES||_tri2->type!=GL_TRIANGLES)
			return false;
		if(fabs(_tri1->normal[0]*_tri2->normal[0])<0.999)
			return false;

		Vector3d temp_vertex[4];
		int temp_index_node[4];
		int temp_shared_index[2];
		Vector3d temp_shared_vertex[2];
		int count=0;
		Facet temp_facet;
		temp_facet=*_tri1;

		//共有されている2点を抽出する
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				if(_tri1->index_node[i]==_tri2->index_node[j]){
					temp_shared_index[count]=_tri1->index_node[i];
					temp_shared_vertex[count]=_tri1->vertex[i];
					count++;
				}
			}
		}
		if(count!=2)return false;
		//共有されていない2点を抽出する
		count=0;
		for(int i=0;i<3;i++){
			count=0;
			for(int j=0;j<2;j++){
				if(_tri1->index_node[i]!=temp_shared_index[j])
					count++;
				if(count==2){
					temp_index_node[0]=_tri1->index_node[i];
					temp_vertex[0]=_tri1->vertex[i];
				}
			}
		}
		count=0;
		for(int i=0;i<3;i++){
			count=0;
			for(int j=0;j<2;j++){
				if(_tri2->index_node[i]!=temp_shared_index[j])
					count++;
				if(count==2){
					temp_index_node[2]=_tri2->index_node[i];
					temp_vertex[2]=_tri2->vertex[i];
				}
			}
		}
		//順番を考慮して一時変数に格納する
		temp_index_node[1]=temp_shared_index[0];
		temp_index_node[3]=temp_shared_index[1];
		temp_vertex[1]=temp_shared_vertex[0];
		temp_vertex[3]=temp_shared_vertex[1];

		//ポリゴン種類の変更
		_tri1->setFacetTypeAsPolygon();
		_tri2->setFacetTypeAsPolygon();	

		//値の継承
		_tri1->index_facet=temp_facet.index_facet;
		_tri1->normal_type=temp_facet.normal_type;
		_tri1->index_material=temp_facet.index_material;
		memcpy(_tri1->normal,temp_facet.normal,sizeof(Vector3d)*4);
		memcpy(_tri1->index_normal,temp_facet.index_normal,sizeof(int)*4);
		memcpy(_tri1->vertex,temp_vertex,sizeof(Vector3d)*4);
		memcpy(_tri1->index_node,temp_index_node,sizeof(int)*4);

		_tri2->index_facet=temp_facet.index_facet;
		_tri2->normal_type=temp_facet.normal_type;
		_tri2->index_material=temp_facet.index_material;
		memcpy(_tri2->normal,temp_facet.normal,sizeof(Vector3d)*4);
		memcpy(_tri2->index_normal,temp_facet.index_normal,sizeof(int)*4);
		_tri2->is_enabled=false;
		memcpy(_tri2->vertex,temp_vertex,sizeof(Vector3d)*4);
		memcpy(_tri2->index_node,temp_index_node,sizeof(int)*4);

		return true;
	}

	void VolumeMesh::setup()
	{
		this->is_loaded = true;
		this->calCenter();
		this->calScale();
		this->calTransform();
		this->calLine();
		this->calFacet();
	}
}