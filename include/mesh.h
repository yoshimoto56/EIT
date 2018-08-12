#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <GL/freeglut.h>
#include <math.h>
#include <vectormatrix.h>
#include <transfermatrix.h>
#include <glutilities.h>
#include <pointlineface.h>

#define NUM_POSITION_MAX 6

//COLOR SETTINGS
#define GL_COLOR_SELECTED glColor3d(0.1, 1.0, 0.1)
#define COLOR_SELECTED Vector4f(0.1, 1.0, 0.1, 0.0)
#define COLOR_NODE_DEFAULT Vector3f(0.1, 0.5, 0.5)
#define COLOR_LINE_DEFAULT Vector3f(0.1, 0.1, 0.1)
#define COLOR_FACET_DEFAULT Vector3f(0.1, 0.1, 1.0)


namespace EITS{

	enum SELECT_MESH{
		SELECT_NONE,
		SELECT_NODE,
		SELECT_LINE,
		SELECT_FACET
	};

	enum NORMAL_TYPE{
		NORMAL_NONE,
		NORMAL_NODE,
		NORMAL_FACET,
	};

	enum COLOR_TYPE {
		COLOR_NONE,
		COLOR_FACET,
		COLOR_NODE
	};

	typedef struct{
		bool object;
		bool line;
		bool node;
		bool facet;
		bool element;
		bool label;
	}IsView;

	typedef struct {
		bool auto_trans;
		bool alpha_blend;
		bool display_list;
		bool vertex_color;
	}IsEnable;

	class Node{
	public:
		Node(){this->clear();}
		~Node(){}
		Node &operator = (const Node &_node){
			memcpy(&position, &_node.position, sizeof(Vector3d));
			memcpy(&normal, &_node.normal, sizeof(Vector3d));
			index = _node.index;
			state = _node.state;
			is_selected = _node.is_selected;
			return (*this);
		}
		//Node index
		int index;
		//State
		int state;
		//Label
		int label;
		//Coordinate
		Vector3d position;
		//Vertex normal
		Vector3d normal;
		//Vertex color;
		Vector3f color;
		//Selection interface
		bool is_selected;

		void clear(){
			this->index = -1;
			this->state = 0;
			this->label = -1;
			this->position = Vector3d(0,0,0);
			this->normal = Vector3d(0, 0, 0);
			this->color = COLOR_NODE_DEFAULT;
			this->is_selected = false;
		}
		void render() {
			if (this->is_selected)
				GL_COLOR_SELECTED;
			else{
				glColor3fv(color.X);
			}
			glBegin(GL_POINTS);
			glVertex3dv(position.X);
			glEnd();
		}
	};

	class Line
	{
	public:
		Line(){
			is_selected = false;
			color = COLOR_LINE_DEFAULT;
		}
		~Line(){}
		bool is_selected;
		Vector3f color;
		Vector3d position[2];
		int index_node[2];
		Line &operator = (const Line &_line){
			memcpy(position, _line.position, sizeof(Vector3d) * 2);
			memcpy(index_node, _line.index_node, sizeof(int) * 2);
			this->is_selected = _line.is_selected;
			this->color = _line.color;
			return (*this);
		}
		void render(){
			if (this->is_selected)
				GL_COLOR_SELECTED;
			else {
				glColor3fv(color.X);
			}
			glBegin(GL_LINE_STRIP);
			glVertex3dv(this->position[0].X);	
			glVertex3dv(this->position[1].X);	
			glEnd();		
		}
		double getLength(){return (position[0]-position[1]).abs();}
	};

	class Facet
	{
	public:
		GLenum type;
		int normal_type;
		int color_type;
		Line line[NUM_POSITION_MAX];
		int index_node[NUM_POSITION_MAX];
		int index_normal[NUM_POSITION_MAX];
		Vector3d position[NUM_POSITION_MAX];
		Vector2d position_local[NUM_POSITION_MAX];
		Vector3d normal[NUM_POSITION_MAX];
		Vector3f color[NUM_POSITION_MAX];
		Material material;
		bool is_selected;
		bool is_enabled;
		int index_facet;
		int index_material;
		int num_node;
		int num_normal;
		int num_line;
		int index_elem;
		Vector3d min;
		Vector3d max;
		double area;
		Facet &operator = (const Facet &_facet){
			type = _facet.type;
			this->is_selected = _facet.is_selected;
			this->is_enabled = _facet.is_enabled;
			index_facet = _facet.index_facet;
			index_material = _facet.index_material;
			num_node = _facet.num_node;
			num_line = _facet.num_line;
			num_normal = _facet.num_normal;
			normal_type = _facet.normal_type;
			color_type = _facet.color_type;
			index_elem = _facet.index_elem;
			area = _facet.area;
			memcpy(index_node, _facet.index_node, sizeof(int)*num_node);
			memcpy(position, _facet.position, sizeof(Vector3d)*num_node);
			memcpy(index_normal, _facet.index_normal, sizeof(int)*num_normal);
			memcpy(normal, _facet.normal, sizeof(Vector3d)*num_normal);
			memcpy(line, _facet.line, sizeof(Line)*num_line);
			memcpy(color, _facet.color, sizeof(Vector3f)*num_line);
			return (*this);
		}
		Facet(GLenum _type=GL_TRIANGLES, int _num_node=0,	int _normal_type=NORMAL_NODE) 
			: type(_type), num_node(_num_node), num_normal(_num_node), normal_type(_normal_type),is_selected(false)
			,color_type(COLOR_FACET),is_enabled(true),index_material(0),num_line(0){
				this->setFacetTypeAsTriangle();
				this->color[0] = COLOR_FACET_DEFAULT;
		}
		~Facet(){
		}
		void setFacetTypeAsTriangle(){
			type = GL_TRIANGLES;
			num_node = 3;
			num_line = 3;
			num_normal = 3;
			index_material = 0;
		}
		void setFacetTypeAsPolygon(int _num=4){
			if(_num<NUM_POSITION_MAX){
				type = GL_POLYGON;
				num_node = _num;
				num_line = _num;
				num_normal = _num;
				index_material = 0;
			}
		}
		Vector3d getMinPosition(){
			min = Vector3d(INT_MAX, INT_MAX, INT_MAX);
			for (int i = 0;i<num_node;i++) {
				if (min.x>position[i].x)min.x = position[i].x;
				if (min.y>position[i].y)min.y = position[i].y;
				if (min.z>position[i].z)min.z = position[i].z;
			}
			return min;
		}
		Vector3d getMaxPosition(){
			max = Vector3d(-INT_MAX, -INT_MAX, -INT_MAX);
			for (int i = 0;i<num_node;i++) {
				if (max.x<position[i].x)max.x = position[i].x;
				if (max.y<position[i].y)max.y = position[i].y;
				if (max.z<position[i].z)max.z = position[i].z;
			}
			return max;
		}
		void calPositionLocal(){
			Vector3d e1 = (this->position[1] - this->position[0])/(this->position[1] - this->position[0]).abs();
			Vector3d e2 = e1 % this->normal[0];
			e2 /=e2.abs();
			this->position_local[0] = Vector2d(0, 0);
			this->position_local[1] = Vector2d(e1*(this->position[1] - this->position[0]), e2*(this->position[1] - this->position[0]));
			this->position_local[2] = Vector2d(e1*(this->position[2] - this->position[0]), e2*(this->position[2] - this->position[0]));
		}
		Matrixd calArea(){
			this->calPositionLocal();
			Matrixd result(3, 3);
			Matrixd temp(3, 3);
			for(int i = 0;i < 3;i++){
				for(int j = 0;j < 3;j++){
					if(j == 0)temp.X[3 * i + j] = 1;
					else temp.X[3 * i + j] = this->position_local[i].X[ j - 1];
				}
			}
			this->area = fabs(inverseLU(result.X, temp.X, 3)) / 2.0;
			return result;
		}
		void render() {
			if (this->normal_type == NORMAL_FACET)
				glNormal3dv(this->normal[0].X);
			if (this->color_type == COLOR_FACET) {
				if (this->is_selected)
					GL_COLOR_SELECTED;
				else {
					glColor3fv(color[0].X);
				}
			}
			glBegin(this->type);
			for (int j = 0; j < this->num_node; j++) {
				if (this->normal_type == NORMAL_NODE)
					glNormal3dv(normal[j].X);
				if (this->color_type == COLOR_NODE) {
					if (this->is_selected)
						GL_COLOR_SELECTED;
					else {
						glColor3fv(color[j].X);
					}
				}
				glVertex3dv(this->position[j].X);
			}
			glEnd();
		}
	};

	class Tetrahedra{
	public:
		Tetrahedra(){this->position.malloc(12); this->clear();}
		~Tetrahedra(){this->position.free();}
		//ノードの座標
		VectorNd position;
		//材料番号
		int index_material;
		//ノード番号
		int index_node[4];
		//要素内のノード順序
		int order_node[16];
		//中心座標
		Vector3d center;
		//中心座標の計算
		void calCenter(){
			this->center = Vector3d(0,0,0);
			for(int i=0;i<4;i++)
				for(int j=0;j<3;j++)
					this->center.X[j] += this->position.X[3 * i + j] / 4;
		}
		//ノードの法線
		Vector3d normal[4];
		//法線の計算
		void calNormal(){
			Matrixd temp_coord(3,4);
			for(int i=0;i<4;i++){
				for(int j=0;j<3;j++)
					for(int k=0;k<4;k++)
						temp_coord.X[3*k+j] = position.X[3 * order_node[4*i+k]+j];
				for(int j=0; j<3; j++){
					normal[i].X[j] = temp_coord.X[3*1+(j+1)%3]*temp_coord.X[3*2+(j+2)%3]-temp_coord.X[3*1+(j+1)%3]*temp_coord.X[3*0+(j+2)%3]
									-temp_coord.X[3*0+(j+1)%3]*temp_coord.X[3*2+(j+2)%3]-temp_coord.X[3*1+(j+2)%3]*temp_coord.X[3*2+(j+1)%3]
									+temp_coord.X[3*1+(j+2)%3]*temp_coord.X[3*0+(j+1)%3]+temp_coord.X[3*0+(j+2)%3]*temp_coord.X[3*2+(j+1)%3];
				}
				if(normal[i].abs()!=0)
					normal[i]/=-normal[i].abs();	
				if(normal[i].X[0]*(temp_coord.X[3*3+0]-temp_coord.X[3*0+0])
					+normal[i].X[1]*(temp_coord.X[3*3+1]-temp_coord.X[3*0+1])
					+normal[i].X[2]*(temp_coord.X[3*3+2]-temp_coord.X[3*0+2])<0){
					normal[i]*=-1;
				}
			}
		}
		//体積
		double volume;
		//体積の計算
		Matrixd calVolume(){
			Matrixd result(4, 4);
			Matrixd temp(4, 4);
			for(int i = 0;i < 4;i++){
				for(int j = 0;j < 4;j++){
					if(j == 0)temp.X[4 * i + j] = 1;
					else temp.X[4 * i + j] = position.X[ 3 * i + j - 1];
				}
			}
			this->volume = fabs(inverseLU(result.X, temp.X, 4)) / 6.0;
			return result;
		}

		void clear(){
			this->index_material = 0;
			this->center = Vector3d(0,0,0);
			memset(this->index_node, -1, sizeof(int) * 4);
			this->volume = 0;
			order_node[0]=2;order_node[1]=1;order_node[2]=0;order_node[3]=3;
			order_node[4]=2;order_node[5]=3;order_node[6]=1;order_node[7]=0;
			order_node[8]=0;order_node[9]=3;order_node[10]=2;order_node[11]=1;
			order_node[12]=0;order_node[13]=1;order_node[14]=3;order_node[15]=2;
			memset(this->position.X, 0, sizeof(double) * this->position.n);
			for(int i=0;i<4;i++){
				this->normal[i]=Vector3d(0,0,0);
			}
		}
	};

	class SurfMesh
	{
	protected:
		int num_node;
		int num_facet;
		int num_line;
		int num_material;

		Facet *facet;
		Line *line;
		Node *node;
		Material *material;
		Vector3f *color;

		Material material_selected;
		bool is_selected;

		double scale;
		Vector3d size;
		Vector3d center;

		GLuint id;
		bool is_listed;
		bool is_loaded;
		bool is_line_extracted;

		void calFacetNormal();
		void calFacetVertex();
		void calLine();
		void calNode();
		void calCenter();
		void calScale();
		void calTransform();

		void addNode(Node &_node);
		void addFacet(Facet &_facet);
		void addMaterial(Material &_material);
	public:
		SurfMesh(void);
		~SurfMesh(void);

		localObject Tr;
		IsView is_view;
		IsEnable is_enable;

		bool load(const char *filename);
		bool save(const char *filename);

		void newMesh();
		void newMaterial();
		void deleteMesh();
		void deleteMaterial();
		void clearMesh();
		void clearView();
		void clearEnable();

		bool getIsLoaded(){return this->is_loaded;}

		Vector3d getSize(){return size;}
		double getScale(){return scale;}
		Vector3d getCenter(){return center;}
		Node getNode(int _index){return node[_index];}
		Node* getNodePointer(){return node;}
		Line getLine(int _index){return this->line[_index];}
		Vector3f getColor(int _index){return color[_index];}
		void setNode(int _index, Node _node){node[_index]=_node;}
		Facet getFacet(int _index){return facet[_index];}
		Facet* getFacetPointer(int _index){return &facet[_index];}
		void setFacet(int _index, Facet _facet){facet[_index]=_facet;}
		void setLine(int _index, Line _line){line[_index]=_line;}
		Material getMaterial(int _index){return this->material[_index];}
		Material* getMaterialPointerAt(int _index){return &this->material[_index];}

		std::ostream& getInfo(std::ostream &stream);
		void render();
		void renderList();
		GLuint makeDisplayList(int _id);

		void setColorAt(int _index, Vector3f _color){this->color[_index]=_color;}
		void setNumNode(int _num_node){this->num_node=_num_node;}
		int getNumNode() { return num_node; }
		void setNumLine(int _num_line) { this->num_line = _num_line; }
		int getNumLine() { return num_line; }
		void setNumFacet(int _num_facet){this->num_facet=_num_facet;}
		int getNumFacet() { return num_facet; }
		void setNumMaterial(int _num_material){this->num_material=_num_material;}
		int getNumMaterial() { return num_material; }

		Node* getNodeAt(int _index) { return &this->node[_index]; }
		Line* getLineAt(int _index) { return &this->line[_index]; }
		Facet* getFacetAt(int _index) { return &this->facet[_index]; }

		//Selection interface
		void setIsSelected(bool _is_selected) { this->is_selected = _is_selected; }
		bool getIsSelected() { return this->is_selected; }
		void selectAt(Vector3d _coord, int _type);
		void selectWithin(Vector3d *_coord, int _type);
		void clearSelection(int _type = SELECT_NONE);

		void quadrize();
		void centerize();
	};

	class VolumeMesh{

	protected:
		//Flag for data load
		bool is_loaded;
		//Node
		Node *node;
		//Number of nodes
		int num_node;
		//Element
		Tetrahedra *elem;
		//Number of element
		int num_elem;
		//Line
		Line *line;
		//Number of line
		int num_line;
		//Facet
		Facet *facet;
		//Number of facet
		int num_facet;

		//Center of object
		Vector3d center;
		//Size of object
		Vector3d size;
		//Scale of object
		double scale;

		bool is_line_extracted;
		bool is_facet_extracted;
		bool is_selected;

		void calCenter();
		//Calculate the scale
		void calScale();
		//Calculate Transformation matrix
		void calTransform();
		//Delete duplicated node
		void calNode();
		//Extract line from elements
		void calLine();
		//Extract facet from elements
		void calFacet();
	public:
		VolumeMesh();
		~VolumeMesh();

		IsEnable is_enable;
		IsView is_view;
		localObject Tr;

		void newMesh();
		void deleteMesh();

		//Clear FEM data
		void clear();
		void clearView();
		void clearEnable();
		//Load FEM data
		bool load(const char* _filename);
		//Save FEM data
		bool save(const char* _filename);
		//Calculate the scale and the center coordinate of the object

		void setup();

		//Information about drawing object
		std::ostream& getInfo(std::ostream &stream);
		void render();
		int makeDisplayList(int _id);

		//Accessor
		bool getIsLoaded() { return this->is_loaded; }
		void setIsLoaded(bool _is_loaded) { this->is_loaded = _is_loaded; }
		void setNumNode(int _num_node) { this->num_node = _num_node; }
		int getNumNode(){return this->num_node;}
		void setNumLine(int _num_line) { this->num_line = _num_line; }
		int getNumLine() { return this->num_line; }
		void setNumFacet(int _num_facet) { this->num_facet = _num_facet; }
		int getNumFacet() { return this->num_facet; }
		void setNumElem(int _num_elem) { this->num_elem = _num_elem; }
		int getNumElem() { return this->num_elem; }
		Vector3d getCenter(){return this->center;}
		Vector3d getSize(){return this->size;}
		Node* getNodeAt(int _index){return &this->node[_index];}
		Line* getLineAt(int _index) { return &this->line[_index]; }
		Facet* getFacetAt(int _index) { return &this->facet[_index]; }
		Tetrahedra* getElemAt(int _index){return &this->elem[_index];}

		//Selection interface
		void setIsSelected(bool _is_selected) { this->is_selected = _is_selected; }
		bool getIsSelected() { return this->is_selected; }
		void selectAt(Vector3d _coord, int _type);
		void selectWithin(Vector3d *_coord, int _type);
		void clearSelection(int _type = SELECT_NONE);

	};

	bool isSharedLine(Line _new, Line *_old, int _num);
	int isSharedFacet(Facet _new, Facet *_old, int _num);
	bool convertTriangle2Quad(Facet *_tri1, Facet *_tri2);
	//bool isCrossFacet(Facet *_facet1, Facet *_facet2, Vector3d *_Xend);
}