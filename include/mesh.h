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

#define NUM_VERTEX_MAX 6

namespace EITS{

	enum SELECT_MESH{
		SELECT_NONE,
		SELECT_FACET,
		SELECT_LINE,
		SELECT_POINT
	};

	enum NORMAL_TYPE{
		NORMAL_NONE,
		NORMAL_FACET,
		NORMAL_POINT
	};

	class Node{
	public:
		Node(){this->clear();}
		~Node(){}
		Node &operator = (const Node &_node){
			memcpy(&vertex, &_node.vertex, sizeof(Vector3d));
			memcpy(&normal, &_node.normal, sizeof(Vector3d));
			index = _node.index;
			state = _node.state;
			is_selected = _node.is_selected;
			return (*this);
		}
		//ノードの番号
		int index;
		//ノードの状態
		int state;
		//ノードのラベル
		int label;
		//座標
		Vector3d vertex;
		//法線
		Vector3d normal;
		//Selection interface
		bool is_selected;

		void clear(){
			this->index = -1;
			this->state = 0;
			this->label = -1;
			this->vertex = Vector3d(0,0,0);
			this->normal = Vector3d(0,0,0);
			this->is_selected = false;
		}
	};

	class Line
	{
	public:
		Line(){is_selected=false;}
		~Line(){}
		bool is_selected;
		Vector3d vertex[2];
		int index_node[2];
		Line &operator = (const Line &_line){
			memcpy(vertex,_line.vertex,sizeof(Vector3d)*2);
			memcpy(index_node,_line.index_node,sizeof(int)*2);
			is_selected=_line.is_selected;
			return (*this);
		}
		void render(){
			glBegin(GL_LINE_STRIP);
			glVertex3dv(this->vertex[0].X);	
			glVertex3dv(this->vertex[1].X);	
			glEnd();		
		}
		double getLength(){return (vertex[0]-vertex[1]).abs();}
	};

	class Facet
	{
	public:
		GLenum type;

		Line line[NUM_VERTEX_MAX];
		int index_node[NUM_VERTEX_MAX];
		int index_normal[NUM_VERTEX_MAX];
		Vector3d vertex[NUM_VERTEX_MAX];
		Vector2d vertex_local[NUM_VERTEX_MAX];
		Vector3d normal[NUM_VERTEX_MAX];

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
		//面積
		double area;

		int normal_type;
		Facet &operator = (const Facet &_facet){
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
			return (*this);
		}
		Facet(GLenum _type=GL_TRIANGLES, int _num_node=0,	int _normal_type=NORMAL_POINT) 
			: type(_type), num_node(_num_node), num_normal(_num_node), normal_type(_normal_type),is_selected(false)
			,is_enabled(true),index_material(0),num_line(0){
				this->setFacetTypeAsTriangle();
		}
		~Facet(){
		}
		void setFacetTypeAsTriangle(){
			type=GL_TRIANGLES;
			num_node=3;
			num_line=3;
			num_normal=3;
			index_material=0;
		}
		void setFacetTypeAsPolygon(int _num=4){
			if(_num<NUM_VERTEX_MAX){
				type=GL_POLYGON;
				num_node=_num;
				num_line=_num;
				num_normal=_num;
				index_material=0;
			}
		}
		Vector3d getMinVertex(){
			min=Vector3d(INT_MAX,INT_MAX,INT_MAX);
			for(int i=0;i<num_node;i++){
				if(min.x>vertex[i].x)min.x=vertex[i].x;
				if(min.y>vertex[i].y)min.y=vertex[i].y;
				if(min.z>vertex[i].z)min.z=vertex[i].z;
			}
			return min;
		}
		Vector3d getMaxVertex(){
			max=Vector3d(-INT_MAX,-INT_MAX,-INT_MAX);
			for(int i=0;i<num_node;i++){
				if(max.x<vertex[i].x)max.x=vertex[i].x;
				if(max.y<vertex[i].y)max.y=vertex[i].y;
				if(max.z<vertex[i].z)max.z=vertex[i].z;
			}
			return max;
		}
		bool getIsSelected(){return this->is_selected;}
		void setIsSelected(bool _is_selected){is_selected=_is_selected;}
		void calVertexLocal(){
			Vector3d e1 = (this->vertex[1] - this->vertex[0])/(this->vertex[1] - this->vertex[0]).abs();
			Vector3d e2 = e1 % this->normal[0];
			e2 /=e2.abs();
			this->vertex_local[0] = Vector2d(0, 0);
			this->vertex_local[1] = Vector2d(e1*(this->vertex[1] - this->vertex[0]), e2*(this->vertex[1] - this->vertex[0]));
			this->vertex_local[2] = Vector2d(e1*(this->vertex[2] - this->vertex[0]), e2*(this->vertex[2] - this->vertex[0]));
		}
		Matrixd calArea(){
			this->calVertexLocal();
			Matrixd result(3, 3);
			Matrixd temp(3, 3);
			for(int i = 0;i < 3;i++){
				for(int j = 0;j < 3;j++){
					if(j == 0)temp.X[3 * i + j] = 1;
					else temp.X[3 * i + j] = this->vertex_local[i].X[ j - 1];
				}
			}
			this->area = fabs(inverseLU(result.X, temp.X, 3)) / 2.0;
			return result;
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
		int num_normal;
		int num_facet;
		int num_line;
		int num_material;

		Vector3d *vertex;
		Vector3d *normal;
		Facet *facet;
		Line *line;
		Material *material;
		Material material_s;
		Material material_c;
		Vector3f *color;

		int *label_index;
		bool *is_selected;

		double scale;
		Vector3d size;
		Vector3d center;
		Vector3d selected_vertex;
		int selected_index;

		GLuint id;
		bool is_display_list;
		bool is_listed;
		bool is_loaded;
		bool is_view;
		bool is_view_line;
		bool is_view_node;
		bool is_view_facet;
		bool is_auto_scale;
		bool is_tri;
		bool is_identity;
		bool is_view_label;
		bool is_alpha_blend;
		bool is_view_line_extracted;
		bool is_vertex_color_enabled;

		void calFacetNormal();
		void calFacetVertex();
		void calLine();
		void calVertex();
		void calCenter();
		void calScale();
		void calTransform();

		void addVertex(Vector3d &_vertex);
		void addNormal(Vector3d &_normal);
		void addFacet(Facet &_facet);
		void addMaterial(Material &_material);
	public:
		SurfMesh(void);
		~SurfMesh(void);

		localObject Tr;

		bool load(const char *filename);
		bool save(const char *filename);

		void newMesh();
		void newMaterial();
		void deleteMesh();
		void deleteMaterial();
		void clearMesh();
		void clearView();

		bool getIsTri(){return this->is_tri;}
		bool getIsLoaded(){return this->is_loaded;}
		bool getIsView(){return this->is_view;}
		bool getIsAutoScale(){return this->is_auto_scale;}
		bool getIsIdentity(){return this->is_identity;}
		bool getIsSelected(){return this->is_selected;}
		bool getIsViewNode(){return this->is_view_node;}
		bool getIsViewLine(){return this->is_view_line;}
		bool getIsViewFacet(){return this->is_view_facet;}
		bool getIsDisplayList(){return this->is_display_list;}
		bool getIsVertex_color_enabled(){return this->is_vertex_color_enabled;}
		bool getIsViewLabel(){return this->is_view_label;}

		int getNumNode() { return num_node; }
		int getNumNormal() { return num_normal; }
		int getNumFacet() { return num_facet; }
		int getNumMaterial() { return num_material; }
		int getNumLine() {return num_line;}

		Vector3d getSize(){return size;}
		double getScale(){return scale;}
		Vector3d getCenter(){return center;}
		Vector3d getVertex(int _index){return vertex[_index];}
		Vector3d* getVertexPointer(){return vertex;}
		Vector3d getNormal(int _index){return normal[_index];}
		Line getLine(int _index){return this->line[_index];}
		Vector3f getColor(int _index){return color[_index];}
		void setVertex(int _index, Vector3d _vertex){vertex[_index]=_vertex;}
		Facet getFacet(int _index){return facet[_index];}
		Facet* getFacetPointer(int _index){return &facet[_index];}
		void setFacet(int _index, Facet _facet){facet[_index]=_facet;}
		void setNormal(int _index, Vector3d _normal){normal[_index]=_normal;}
		void setLine(int _index, Line _line){line[_index]=_line;}
		Material getMaterial(int _index){return this->material[_index];}
		Material* getMaterialPointerAt(int _index){return &this->material[_index];}
		int getLabelIndex(int _index){return this->label_index[_index];}

		std::ostream& getInfo(std::ostream &stream);
		void render();
		void renderList();
		GLuint makeDisplayList(int _id);

		void setColorAt(int _index, Vector3f _color){this->color[_index]=_color;}
		void setLabelIndexAt(int _index, int _label){this->label_index[_index] = _label;}
		void setNumNode(int _num_node){this->num_node=_num_node;}
		void setNumNormal(int _num_normal){this->num_normal=_num_normal;}
		void setNumFacet(int _num_facet){this->num_facet=_num_facet;}
		void setNumMaterial(int _num_material){this->num_material=_num_material;}
		void setNumLine(int _num_line){this->num_line=_num_line;}
		void setIsView(bool _is_view){this->is_view=_is_view;}
		void setIsViewLabel(bool _is_view_label){this->is_view_label=_is_view_label;}
		void setIsViewLine(bool _is_view_line){this->is_view_line=_is_view_line;}
		void setIsViewFacet(bool _is_view_facet){this->is_view_facet=_is_view_facet;}
		void setIsViewNode(bool _is_view_node){this->is_view_node=_is_view_node;}
		void setIsAutoScale(bool _is_auto_scale){this->is_auto_scale=_is_auto_scale;}
		void setIsTri(bool _is_tri){this->is_tri=_is_tri;}
		void setIsIdentity(bool _is_identity){this->is_identity=_is_identity;}
		void setIsAlpha_blend(bool _is_alpha_blend){this->is_alpha_blend=_is_alpha_blend;}
		void setIsVertex_color_enabled(bool _is_vertex_color_enabled){this->is_vertex_color_enabled=_is_vertex_color_enabled;}
		void setIsDisplayList(bool _is_display_list){this->is_display_list=_is_display_list;}

		//描画座標に最も近い頂点番号を取得する
		int getVertexIndexNear(Vector3d _pos);
		int getFacetIndexNear(Vector3d _pos);
		Vector3d getSelectedVertex(){return this->selected_vertex;}
		void setSelectedVertex(int _index){this->selected_vertex=this->vertex[_index];}
		void setSelectedIndex(int _selected_index){this->selected_index=_selected_index;}
		void clearSelection();

		void quadrize();
		void centerize();
	};

	class VolumeMesh{

	protected:
		//読み込みのフラグ
		bool is_loaded;
		//ノード
		Node *node;
		//ノード数
		int num_node;
		//要素
		Tetrahedra *elem;
		//要素数
		int num_elem;
		//線
		Line *line;
		//線数
		int num_line;
		//面
		Facet *facet;
		//面数
		int num_facet;

		localObject Tr;

		//物体のサイズを描画空間に合わせるかどうか
		bool is_auto_scale;
		//選択モード
		int selection_mode;
		//物体の中心
		Vector3d center;
		//サイズ
		Vector3d size;
		//スケール
		double scale;

		bool is_view_line_extracted;
		bool is_view_facet_extracted;

		bool is_view_info;
		bool is_view_node;
		bool is_view_line;
		bool is_view_facet;
		bool is_view_smooth;

		void calCenter();
		//Calculate the scale
		void calScale();
		//Calculate Transformation matrix
		void calTransform();
		//線の抽出
		void calLine();
		//面の抽出
		void calFacet();
	public:
		VolumeMesh();
		~VolumeMesh();

		void newMesh();
		void deleteMesh();

		//Clear FEM data
		void clear();
		void clearView();
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
		void setNumNode(int _num_node) { this->num_node = _num_node; }
		int getNumNode(){return this->num_node;}
		void setNumElem(int _num_elem) { this->num_elem = _num_elem; }
		int getNumElem() { return this->num_elem; }
		void setNumLine(int _num_line) { this->num_line = _num_line; }
		int getNumLine() { return this->num_line; }
		void setNumFacet(int _num_facet) { this->num_facet = _num_facet; }
		int getNumFacet() { return this->num_facet; }
		Vector3d getCenter(){return this->center;}
		Vector3d getSize(){return this->size;}
		bool getIsLoaded(){return this->is_loaded;}
		Node* getNodeAt(int _elem){return (&this->node[_elem]);}
		Tetrahedra* getElemAt(int _elem){return (&this->elem[_elem]);}
		void setIsViewNode(bool _is_view_node){this->is_view_node=_is_view_node;}
		void setIsViewFacet(bool _is_view_facet) { this->is_view_facet = _is_view_facet; }
		void setIsViewLine(bool _is_view_line){this->is_view_line = _is_view_line;}
		void setIsViewInfo(bool _is_view_info){this->is_view_info = _is_view_info;}
		void setIsAutoScale(bool _is_auto_scale){ is_auto_scale=_is_auto_scale;}
		void setIsLoaded(bool _is_loaded) { this->is_loaded = _is_loaded; }

	};

	bool isSharedLine(Line _new, Line *_old, int _num);
	bool isSharedFacet(Facet _new, Facet *_old, int _num);
	bool convertTriangle2Quad(Facet *_tri1, Facet *_tri2);
	//bool isCrossFacet(Facet *_facet1, Facet *_facet2, Vector3d *_Xend);
}