#include "mesh.h"


namespace EITS{
	/****************Surface Mesh******************/
	SurfMesh::SurfMesh(void)
		:num_node(0),num_facet(0),num_line(0),num_material(0)
	{
		this->newMesh();
		this->newMaterial();
		this->clearEnable();
		this->clearView();
	}

	SurfMesh::~SurfMesh(void)
	{
		this->deleteMesh();
		this->deleteMaterial();
	}

	void SurfMesh::deleteMesh()
	{
		this->is_loaded = false;
		delete[]node;
		delete[]facet;
		delete[]color;
		delete[]line;
	}

	void SurfMesh::newMesh()
	{
		this->node = new Node[this->num_node + 2];
		this->color = new Vector3f[this->num_node + 2];
		this->facet = new Facet[this->num_facet + 2];
		this->line = new Line[this->num_line + 2];
	}

	void SurfMesh::deleteMaterial()
	{
		delete[]material;
	}

	void SurfMesh::newMaterial()
	{
		this->material = new Material[this->num_material + 2];
		this->material[0].setDiffuse(Vector4f(0.5, 0.5, 0.5, 0.0));
		this->material[0].setAmbient(Vector4f(0.0, 0.0, 0.0, 0.0));
		this->material_selected.setDiffuse(COLOR_SELECTED);
		this->material_selected.setAmbient(COLOR_SELECTED);
		this->material_selected.setSpecular(COLOR_SELECTED);
	}

	void SurfMesh::clearMesh()
	{
		this->is_listed=false;
		this->is_loaded=false;
		this->is_line_extracted=false;
		num_node = 0;
		num_facet = 0;
		num_line = 0;
		num_material = 0;
		this->scale=0;
		this->size=Vector3d(1,1,1);
		this->center=Vector3d(0,0,0);
	}

	void SurfMesh::clearEnable() {
		this->is_enable.auto_trans = true;
		this->is_enable.alpha_blend = false;
		this->is_enable.display_list = false;
		this->is_enable.vertex_color = false;
	}

	void SurfMesh::clearView()
	{
		this->is_view.object = true;
		this->is_view.node = true;
		this->is_view.line = true;
		this->is_view.facet = true;
		this->is_view.label = false;
	}

	void SurfMesh::addNode(Node &_node)
	{
		num_node++;
		node = (Node*)realloc(node, num_node*sizeof(Node));
		node[num_node-1] = _node;
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
		stream << "Number of nodes：" << num_node << std::endl;
		stream << "Number of facets：" << num_facet << std::endl;
		stream << "Number of lines：" << num_line << std::endl;
		stream << "Number of materials：" << num_material << std::endl;
		return stream;
	}

	void SurfMesh::calFacetNormal()
	{
		Vector3d t_normal;
		Vector3d temp[2];
		for(int i=0;i<num_facet;i++){
			this->facet[i].normal_type = NORMAL_FACET;
			this->facet[i].num_normal = 1;
			for(int j=1;j<3;j++)
				temp[j-1] = node[facet[i].index_node[j]].position - node[facet[i].index_node[0]].position;
			for(int j=0;j<3;j++){
				t_normal.X[j] = temp[1].X[(j+1)%3] * temp[0].X[(j+2)%3]
						-temp[1].X[(j+2)%3] * temp[0].X[(j+1)%3];
			}
			if(t_normal.abs()!=0)t_normal/=t_normal.abs();
			this->facet[i].index_normal[0] = i;
			this->facet[i].normal[0] = t_normal;
		}
	}

	void SurfMesh::calFacetVertex()
	{
		for(int i=0;i<this->num_facet;i++){
			for(int j=0;j<this->facet[i].num_node;j++){
				this->facet[i].position[j] = node[this->facet[i].index_node[j]].position;
			}
		}
	}

	void SurfMesh::calLine()
	{
		if(this->num_line > 0)
			delete []this->line;
		this->num_line = 0;
		this->line= new Line[this->num_facet*4];//CAUTION: Supported polygons are triangle and rectangle

		for(int i=0;i<this->num_facet;i++){
			for(int j=0;j<facet[i].num_node;j++){
				facet[i].line[j].index_node[0] = facet[i].index_node[j];
				facet[i].line[j].index_node[1] = facet[i].index_node[(j+1)%facet[i].num_node];
				facet[i].line[j].position[0] = facet[i].position[j];
				facet[i].line[j].position[1] = facet[i].position[(j+1)%facet[i].num_node];
				if(!isSharedLine(facet[i].line[j], line,num_line)){
					this->line[this->num_line] = facet[i].line[j];
					this->num_line++;
				}
			}
		}
		this->is_line_extracted=true;
	}

	void SurfMesh::calNode()
	{
		int count = 1;
		this->node[0].position = this->facet[0].position[0];
		for (int i = 0;i<this->num_facet;i++) {
			for (int j = 0;j<3;j++) {
				for(int k=0;k<count;k++){
					if(this->facet[i].position[j] == this->node[k].position){
						this->facet[i].index_node[j] = k;
						k = count;
					}
					else if(k == count-1){
						this->node[count].position = this->facet[i].position[j];
						this->facet[i].index_node[j] = count;
						count++;
					}
				}
			}
		}
		this->num_node = count;
		Node *t_node = new Node[this->num_node + 2];
		for (int i = 0;i < this->num_node;i++)
			t_node[i] = node[i];
		delete []node;
		this->node = t_node;
	}

	void SurfMesh::calCenter()
	{
		Vector3d min(INT_MAX,INT_MAX,INT_MAX);
		Vector3d max(-INT_MAX,-INT_MAX,-INT_MAX);
		for(int i=0;i<this->num_node;i++){
			for(int k=0;k<3;k++){
				if(min.X[k]>this->node[i].position.X[k])
					min.X[k]=this->node[i].position.X[k];
				if(max.X[k]<this->node[i].position.X[k])
					max.X[k]=this->node[i].position.X[k];
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
		if(!this->is_loaded||!this->is_view.object)
			return;
		glDisable(GL_CULL_FACE); 
		if(this->is_enable.alpha_blend){
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
			glEnable(GL_ALPHA_TEST);
		}
		glEnable(GL_DEPTH_TEST);

		glPushMatrix();

		if(this->is_enable.auto_trans){
			glMultMatrixd(Tr.getModel2Uni().getTr4GL());
		}
		else {
			glMultMatrixd(Tr.getModel2World().getTr4GL());
		}
		if(this->is_enable.display_list){
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
		int pre_mat = -1, cur_mat = 0;

		if (this->is_view.node) {
			glDisable(GL_LIGHTING);
			glPointSize(5);
			for (int i = 0;i<num_node;i++) {
				this->node[i].render();
			}
			glEnable(GL_LIGHTING);
		}

		if(this->is_view.line && this->is_line_extracted){
			glDisable(GL_LIGHTING);
			glLineWidth(0.5);
			for (int i = 0;i < num_line;i++) {
				this->line[i].render();
			}
			glEnable(GL_LIGHTING);
		}

		if (this->is_view.facet) {
			if (this->is_view.label)
				glDisable(GL_LIGHTING);
			for (int i = 0; i<num_facet; i++) {
				if (facet[i].is_enabled) {
					cur_mat = facet[i].index_material;
					if (pre_mat != cur_mat) {
						material[cur_mat].set();
						pre_mat = cur_mat;
					}
					if (facet[i].is_selected) {
						material_selected.set();
						pre_mat = -1;
					}
					this->facet[i].render();
				}
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

	void SurfMesh::selectAt(Vector3d _coord, int _type) 
	{
		double dist;
		double min_dist = 1000;
		int t_index = -1;
		Vector3d Pc_o;
		if (this->is_enable.auto_trans)
			Pc_o = Tr.getUni2Model()*_coord;
		else
			Pc_o = Tr.getModel2World()*_coord;
		if (_type == SELECT_NODE){
			for (int i = 0;i<this->num_node;i++) {
				dist = (this->node[i].position - Pc_o).abs();
				if (dist < min_dist && dist < 10) {
					min_dist = dist;
					t_index = i;
				}
			}
			if (0 <= t_index&&t_index<num_node)
				this->node[t_index].is_selected = !this->node[t_index].is_selected;
		}
		if (_type == SELECT_FACET) {
			for (int i = 0;i<this->num_facet;i++) {
				if (facet[i].is_enabled) {
					if (isProjectedPointOnFiniteFace(Pc_o, &facet[i])) {
						dist = getFacePointDistance(Pc_o, facet[i].normal[0], facet[i].position[0]);
						if (min_dist > fabs(dist)) {
							min_dist = fabs(dist);
							t_index = i;
						}
					}
				}
			}
			if (0 <= t_index&&t_index < num_facet)
				this->facet[t_index].is_selected = !this->facet[t_index].is_selected;
		}
	}

	void SurfMesh::selectWithin(Vector3d *_coord, int _type)
	{
		if (!this->is_loaded)
			return;
		Vector3d select_normal;
		bool is_po_wi_tr;
		for (int i = 0;i<this->num_node;i++) {
			for (int j = 0;j<4;j++) {
				select_normal = (_coord[(j + 1) % 4] - _coord[j])
					% (_coord[j + 4] - _coord[j]);
				select_normal /= select_normal.abs();

				if (this->is_enable.auto_trans)
					is_po_wi_tr = EITS::isViewNodeWithinTriangle(10.0*(this->node[i].position - this->center) / this->scale, select_normal, _coord[j]);
				else
					is_po_wi_tr = EITS::isViewNodeWithinTriangle((this->node[i].position), select_normal, _coord[j]);
				if (!is_po_wi_tr) {
					if (j == 3) {
						if (_type == SELECT_NODE)
							this->node[i].is_selected = true;
						if (_type == SELECT_FACET) {
							for (int k = 0;k<num_facet;k++)
								for(int l=0;l<3;l++)
									if(facet[k].index_node[l] == i)
										this->facet[k].is_selected = true;
						}
					}
				}
				else {
					break;
				}
			}
		}
	}

	void SurfMesh::clearSelection(int _type)
	{
		if (_type == SELECT_FACET || _type == SELECT_NONE) {
			for (int i = 0;i < num_facet;i++) {
				facet[i].is_selected = false;
			}
		}
		if (_type == SELECT_NODE || _type == SELECT_NONE) {
			for (int i = 0;i < num_node;i++) {
				this->node[i].is_selected = false;
			}
		}
		if (_type == SELECT_LINE || _type == SELECT_NONE) {
			for (int i = 0;i < num_line;i++) {
				this->line[i].is_selected = false;
			}
		}
	}

	void SurfMesh::quadrize()
	{
		std::cout<<"Quadrizing surface ...";
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
		std::cout << "Centerizing mesh ...";
		this->calCenter();
		for (int i = 0;i < this->num_node;i++) {
			this->node[i].position -= this->center;
		}
		for (int i = 0;i < this->num_line;i++)
			for (int j = 0;j<2;j++)
				this->line[i].position[j] -= this->center;
		for (int i = 0;i < this->num_facet;i++)
			for (int j = 0;j<this->facet[i].num_node;j++)
				this->facet[i].position[j] -= this->center;
		std::cout << "[OK]" << std::endl;
		this->calCenter();
		this->calTransform();
	}
	/****************Surface Mesh******************/


	/****************Volume Mesh******************/
	VolumeMesh::VolumeMesh():num_elem(0),num_node(0),num_facet(0),num_line(0)
	{
		this->newMesh();
		this->clear();
		this->clearEnable();
		this->clearView();
	}

	VolumeMesh::~VolumeMesh()
	{
	}

	void VolumeMesh::newMesh()
	{
		this->is_loaded = false;
		this->node = new Node[this->num_node + 2];
		this->line = new Line[this->num_line + 2];
		this->facet = new Facet[this->num_facet + 2];
		this->elem = new Tetrahedra[this->num_elem + 2];
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
		this->is_line_extracted = false;
		this->is_facet_extracted = false;
		this->center = Vector3d(0, 0, 0);
		this->size = Vector3d(1, 1, 1);
		this->scale = 1;
	}

	void VolumeMesh::clearView()
	{
		this->is_view.object = true;
		this->is_view.node = true;
		this->is_view.line = true;
		this->is_view.facet = true;
		this->is_view.label = false;
		this->is_view.element = false;
	}

	void VolumeMesh::clearEnable()
	{
		this->is_enable.auto_trans = true;
		this->is_enable.alpha_blend = false;
		this->is_enable.display_list = false;
		this->is_enable.vertex_color = false;
	}

	bool VolumeMesh::load(const char *_filename)
	{
		std::ifstream file;
		char buf[512];
		char dummy[256];
		int num;

		std::cout << "Loading volume mesh data...";
		file.open(_filename, std::ios::in);
		if(!file.is_open())
		{
			std::cout << "[FAIL]" << std::endl;
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
				>> node[i].position.x 
				>> node[i].position.y 
				>> node[i].position.z;
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
			}
		}
		for (int i = 0;i < num_elem;i++) {
			for (int j = 0;j<4;j++)
				for (int k = 0;k<3;k++)
					elem[i].position.X[3 * j + k] = node[elem[i].index_node[j]].position.X[k];
			elem[i].calNormal();
		}

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

		std::cout << "[OK]" << std::endl;
		this->setup();

		return true;
	}

	bool VolumeMesh::save(const char *_filename)
	{
		std::ofstream file;
		std::cout<<"Saving volume mesh data...";
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
			file<<node[i].position.X[0]<<" "<<node[i].position.X[1]<<" "<<node[i].position.X[2]<<std::endl;
		file<<"@2"<<std::endl;
		for(int i=0; i< num_elem; i++){
			file << elem[i].index_node[0] << " " << elem[i].index_node[1] << " " << elem[i].index_node[2] << " " << elem[i].index_node[3] << std::endl;
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
				if(min.X[k]>this->node[i].position.X[k])
					min.X[k]=this->node[i].position.X[k];
				if(max.X[k]<this->node[i].position.X[k])
					max.X[k]=this->node[i].position.X[k];
			}
		}
		this->size.x=fabs(max.x-min.x);
		this->size.y=fabs(max.y-min.y);
		this->size.z=fabs(max.z-min.z);
		this->center=(max+min)/2.0;
#ifdef _DEBUG
		std::cout << "Center: " << this->center << std::endl;
		std::cout << "Size: " << this->size << std::endl;
#endif
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
		stream << "Number of nodes：" << this->num_node << std::endl;
		stream << "Number of lines：" << this->num_line << std::endl;
		stream << "Number of facets：" << this->num_facet << std::endl;
		stream << "Number of tetrahedras：" << this->num_elem << std::endl;
		return stream;
	}

	void VolumeMesh::render(void)
	{
		if(!this->is_loaded)
			return;

		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);
		glPushMatrix();

		if (this->is_enable.auto_trans) {
			glMultMatrixd(Tr.getModel2Uni().getTr4GL());
		}
		else {
			glMultMatrixd(Tr.getModel2World().getTr4GL());
		}

		//Draw node
		if (this->is_view.node) {
			glDisable(GL_LIGHTING);
			glPointSize(5);
			for (int i = 0;i<num_node;i++) {
				this->node[i].render();
			}
			glEnable(GL_LIGHTING);
		}

		//Draw line
		if (this->is_view.line && this->is_line_extracted) {
			glDisable(GL_LIGHTING);
			glLineWidth(0.5);
			for (int i = 0;i < num_line;i++) {
				this->line[i].render();
			}
			glEnable(GL_LIGHTING);
		}

		//Draw facet
		if(this->is_view.facet && this->is_facet_extracted){
			glDisable(GL_LIGHTING);
			for(int i=0;i<this->num_facet;i++){
				this->facet[i].render();
			}
			glEnable(GL_LIGHTING);
		}

		glPopMatrix();
		glEnable(GL_CULL_FACE);

	}

	int VolumeMesh::makeDisplayList(int _id){
		int tid;
		tid = glGenLists(_id+1);
		glNewList(tid, GL_COMPILE);
		this->render();
		glEndList();
		return tid;
	}

	void VolumeMesh::calNode()
	{
		std::cout << "Deleting duplicated node ...";
		if (!this->is_loaded) {
			std::cout << "[FAIL]" << std::endl;
			return;
		}
		
		Node *t_node = new Node[this->num_node + 2];
		int *t_index_map = new int[this->num_node + 2];
		int count = 1;
		t_node[0] = this->node[0];
		t_index_map[0] = 0;
		for (int i = 1;i < this->num_node;i++) {
			for (int j = 0;j < count;j++) {
				if (t_node[j].position == this->node[i].position) {
					t_index_map[i] = j;
					break;
				}
				else {
					if (j == count - 1) {
						t_node[count] = node[i];
						t_index_map[i] = count;
						count++;
					}
				}
			}
		}
		this->num_node = count;
		delete[]node;
		node = t_node;
		for (int i = 0;i < this->num_elem;i++) {
			for(int j=0;j<4;j++)
				this->elem[i].index_node[j] = t_index_map[this->elem[i].index_node[j]];
		}
		delete[]t_index_map;
		std::cout << "[OK]" << std::endl;
	}


	void VolumeMesh::calLine()
	{
		std::cout<<"Extracting lines ...";
		if(!this->is_loaded){
			this->is_facet_extracted = false;
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}

		Line t_line;
		delete []this->line;
		this->num_line = 0;
		this->line = new Line[this->num_elem * 6 + 2];//4面体なので6本

		for(int i=0;i<this->num_elem;i++){
			for(int j=0;j<4;j++){
				for(int k=j;k<4;k++){
					if(j!=k){
						t_line.index_node[0] = elem[i].index_node[j];
						t_line.index_node[1] = elem[i].index_node[k];
						t_line.position[0] = node[elem[i].index_node[j]].position;
						t_line.position[1] = node[elem[i].index_node[k]].position;
						if(!isSharedLine(t_line, this->line, this->num_line)){
							this->num_line++;
							this->line[this->num_line-1] = t_line;
						}
					}
				}
			}
		}
		this->is_line_extracted = true;
		std::cout<<"[OK]"<<std::endl;
	}

	void VolumeMesh::calFacet()
	{
		std::cout<<"Extracting facets ...";
		if(!this->is_loaded){
			this->is_facet_extracted = false;
			std::cout<<"[FAIL]"<<std::endl;
			return;
		}
		Facet t_facet;
		t_facet.setFacetTypeAsTriangle();
		delete []this->facet;
		this->num_facet = 0;
		this->facet= new Facet[this->num_elem * 4 + 2];
		int t_index_duplicated;
		bool *t_is_duplicated = new bool[this->num_elem * 4 + 2];
		memset(t_is_duplicated, false, sizeof(bool)*(this->num_elem * 4 + 2));
		for(int i=0; i<this->num_elem; i++){
			for(int j=0; j<4; j++){
				t_facet.normal[j] = elem[i].normal[j];
				t_facet.index_node[0] = elem[i].index_node[(j+0)%4];
				t_facet.index_node[1] = elem[i].index_node[(j+1)%4];
				t_facet.index_node[2] = elem[i].index_node[(j+2)%4];
				t_facet.position[0] = node[elem[i].index_node[(j+0)%4]].position;
				t_facet.position[1] = node[elem[i].index_node[(j+1)%4]].position;
				t_facet.position[2] = node[elem[i].index_node[(j+2)%4]].position;

				t_index_duplicated = isSharedFacet(t_facet, this->facet, this->num_facet);
				if(t_index_duplicated == -1){
					t_facet.index_facet = this->num_facet;
					t_facet.index_elem = i;
					this->facet[this->num_facet] = t_facet;
					this->num_facet++;
				}
				else {
					t_is_duplicated[t_index_duplicated] = true;
				}
			}
		}
		int count = 0;
		for (int i = 0; i < this->num_facet; i++) {
			this->facet[count] = this->facet[i];
			if (!t_is_duplicated[i]) {
				count++;
			}
		}
		this->num_facet = count;

		delete[]t_is_duplicated;
		this->is_facet_extracted = true;
		std::cout<<"[OK]"<<std::endl;
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

	void VolumeMesh::selectAt(Vector3d _coord, int _type)
	{
		double dist;
		double min_dist = 1000;
		int t_index = -1;
		Vector3d Pc_o;
		if (this->is_enable.auto_trans)
			Pc_o = Tr.getUni2Model()*_coord;
		else
			Pc_o = Tr.getModel2World()*_coord;
		if (_type == SELECT_NODE) {
			for (int i = 0;i<this->num_node;i++) {
				dist = (this->node[i].position - Pc_o).abs();
				if (dist < min_dist && dist < 10) {
					min_dist = dist;
					t_index = i;
				}
			}
			if (0 <= t_index&&t_index<num_node)
				this->node[t_index].is_selected = !this->node[t_index].is_selected;
		}
		if (_type == SELECT_FACET) {
			for (int i = 0;i<this->num_facet;i++) {
				if (facet[i].is_enabled) {
					if (isProjectedPointOnFiniteFace(Pc_o, &facet[i])) {
						dist = getFacePointDistance(Pc_o, facet[i].normal[0], facet[i].position[0]);
						if (min_dist > fabs(dist)) {
							min_dist = fabs(dist);
							t_index = i;
						}
					}
				}
			}
			if (0 <= t_index&&t_index < num_facet)
				this->facet[t_index].is_selected = !this->facet[t_index].is_selected;
		}
	}

	void VolumeMesh::selectWithin(Vector3d *_coord, int _type)
	{
		if (!this->is_loaded)
			return;
		Vector3d select_normal;
		bool is_po_wi_tr;
		for (int i = 0;i<this->num_node;i++) {
			for (int j = 0;j<4;j++) {
				select_normal = (_coord[(j + 1) % 4] - _coord[j])
					% (_coord[j + 4] - _coord[j]);
				select_normal /= select_normal.abs();

				if (this->is_enable.auto_trans)
					is_po_wi_tr = EITS::isViewNodeWithinTriangle(10.0*(this->node[i].position - this->center) / this->scale, select_normal, _coord[j]);
				else
					is_po_wi_tr = EITS::isViewNodeWithinTriangle((this->node[i].position), select_normal, _coord[j]);
				if (!is_po_wi_tr) {
					if (j == 3) {
						if (_type == SELECT_NODE)
							this->node[i].is_selected = true;
						if (_type == SELECT_FACET) {
							for (int k = 0;k<num_facet;k++)
								for (int l = 0;l<3;l++)
									if (facet[k].index_node[l] == i)
										this->facet[k].is_selected = true;
						}
					}
				}
				else {
					break;
				}
			}
		}
	}

	void VolumeMesh::clearSelection(int _type)
	{
		if (_type == SELECT_FACET || _type == SELECT_NONE) {
			for (int i = 0;i < num_facet;i++) {
				facet[i].is_selected = false;
			}
		}
		if (_type == SELECT_NODE || _type == SELECT_NONE) {
			for (int i = 0;i < num_node;i++) {
				this->node[i].is_selected = false;
			}
		}
		if (_type == SELECT_LINE || _type == SELECT_NONE) {
			for (int i = 0;i < num_line;i++) {
				this->line[i].is_selected = false;
			}
		}
	}

	/****************Volume Mesh******************/

	bool isSharedLine(Line _new, Line *_old, int _num)
	{
		for(int i=0;i<_num;i++){
			if((_old[i].index_node[0]==_new.index_node[0]&&_old[i].index_node[1]==_new.index_node[1])
				||(_old[i].index_node[0] == _new.index_node[1] && _old[i].index_node[1] == _new.index_node[0]))
				return true;
		}
		return false;
	}

	int isSharedFacet(Facet _new, Facet *_old, int _num)
	{
		for(int i=0; i<_num; i++){
			if((_new.index_node[0]==_old[i].index_node[0] && _new.index_node[1] == _old[i].index_node[1] && _new.index_node[2] == _old[i].index_node[2]) ||
				(_new.index_node[0]==_old[i].index_node[1] && _new.index_node[1] == _old[i].index_node[2] && _new.index_node[2] == _old[i].index_node[0]) ||
				(_new.index_node[0]==_old[i].index_node[2] && _new.index_node[1] == _old[i].index_node[0] && _new.index_node[2] == _old[i].index_node[1]) ||
				(_new.index_node[0]==_old[i].index_node[0] && _new.index_node[1] == _old[i].index_node[2] && _new.index_node[2] == _old[i].index_node[1]) ||
				(_new.index_node[0]==_old[i].index_node[1] && _new.index_node[1] == _old[i].index_node[0] && _new.index_node[2] == _old[i].index_node[2]) ||
				(_new.index_node[0]==_old[i].index_node[2] && _new.index_node[1] == _old[i].index_node[1] && _new.index_node[2] == _old[i].index_node[0])){
					return i;
			}	
		}
		return -1;
	}

	bool convertTriangle2Quad(Facet *_tri1, Facet *_tri2)
	{
		if(_tri1->type!=GL_TRIANGLES||_tri2->type!=GL_TRIANGLES)
			return false;
		if(fabs(_tri1->normal[0]*_tri2->normal[0])<0.999)
			return false;

		Vector3d temp_position[4];
		int temp_index_node[4];
		int temp_shared_index[2];
		Vector3d temp_shared_position[2];
		int count=0;
		Facet temp_facet;
		temp_facet=*_tri1;

		//共有されている2点を抽出する
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				if(_tri1->index_node[i]==_tri2->index_node[j]){
					temp_shared_index[count]=_tri1->index_node[i];
					temp_shared_position[count]=_tri1->position[i];
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
					temp_position[0]=_tri1->position[i];
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
					temp_position[2]=_tri2->position[i];
				}
			}
		}
		//順番を考慮して一時変数に格納する
		temp_index_node[1]=temp_shared_index[0];
		temp_index_node[3]=temp_shared_index[1];
		temp_position[1]=temp_shared_position[0];
		temp_position[3]=temp_shared_position[1];

		//ポリゴン種類の変更
		_tri1->setFacetTypeAsPolygon();
		_tri2->setFacetTypeAsPolygon();	

		//値の継承
		_tri1->index_facet=temp_facet.index_facet;
		_tri1->normal_type=temp_facet.normal_type;
		_tri1->index_material=temp_facet.index_material;
		memcpy(_tri1->normal,temp_facet.normal,sizeof(Vector3d)*4);
		memcpy(_tri1->index_normal,temp_facet.index_normal,sizeof(int)*4);
		memcpy(_tri1->position,temp_position,sizeof(Vector3d)*4);
		memcpy(_tri1->index_node,temp_index_node,sizeof(int)*4);

		_tri2->index_facet=temp_facet.index_facet;
		_tri2->normal_type=temp_facet.normal_type;
		_tri2->index_material=temp_facet.index_material;
		memcpy(_tri2->normal,temp_facet.normal,sizeof(Vector3d)*4);
		memcpy(_tri2->index_normal,temp_facet.index_normal,sizeof(int)*4);
		_tri2->is_enabled=false;
		memcpy(_tri2->position,temp_position,sizeof(Vector3d)*4);
		memcpy(_tri2->index_node,temp_index_node,sizeof(int)*4);

		return true;
	}

}