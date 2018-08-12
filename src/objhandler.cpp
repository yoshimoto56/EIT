#include "objhandler.h"

namespace EITS
{

	char * getDirectoryName(const char *filename, char *dest)
	{
		char *strings = NULL;
		char *dir;
		dir = new char [strlen(filename)+1];
		strcpy_s(dir, sizeof(char)*(strlen(filename) + 1), filename);//TODO error
		if ( strings = strrchr(dir, '/' ) ) strings[1] = '\0';
		else if ( strings = strrchr(dir, '\\') ) strings[1] ='\0';
		else dir[0] = '\0';
		strcpy_s(dest, sizeof(char)*(strlen(filename) + 1), dir);
		return dest;
	}
	bool ObjMesh::loadOBJFile(const char *filename)
	{
		std::ifstream file;
		int cvi = -1;
		int cfi = -1;
		int cmi =  0;
		int tmp_int;
		char tmp_char[256];
		char buf[256];
		char *pbuf;
		Node t_node;
		strcpy_s(filename_obj, filename);
		file.open(filename, std::ios::in);
		if ( !file.is_open() ){
			return false;
		}
		this->clearMesh();
		this->deleteMesh();
		file.getline(buf, sizeof(buf));
		while( !file.eof() ){
			file.getline(buf, sizeof(buf));
			if(strstr(buf, "of "));
			else if(strstr(buf, "v "))this->num_node++;
			else if(strstr(buf, "f "))this->num_facet++;
			else if(strstr(buf, "mtllib ")){
				sscanf_s(buf, "mtllib %s", &tmp_char, 256);
				sprintf_s(buf,"%s%s", directoryname, tmp_char);
				if ( !loadMTLFile(buf))
					return false;
			}
		}
		this->num_line = num_facet * 4;
		this->newMesh();
		file.clear();
		file.seekg(0, std::fstream::beg);
		while( !file.eof() ){
			file.getline(buf, sizeof(buf));
			if(strstr(buf, "v ")){
				cvi++;
				if ( sscanf_s(buf, "v %lf %lf %lf",&t_node.position.x, &t_node.position.y, &t_node.position.z) == 3 ){
					if ( sscanf_s(buf, "v %lf %lf %lf %f %f %f",&t_node.position.x, &t_node.position.y, &t_node.position.z, &color[cvi].x, &color[cvi].y, &color[cvi].z) == 6 ){
						this->is_enable.vertex_color = true;
						t_node.color = color[cvi];
					}
					this->node[cvi]= t_node;
				}
				else {
					return false;
				}
			}
			else if(strstr(buf, "f ")){
				cfi++;
				pbuf = buf;
				this->facet[cfi].num_node = 0;
				while ( *pbuf ){
					if ( *pbuf == ' ' ) this->facet[cfi].num_node++;
					pbuf++;
				}
				if ( this->facet[cfi].num_node < 3 ){
					return false;
				}
				else if ( this->facet[cfi].num_node ==3){
					this->facet[cfi].setFacetTypeAsTriangle();
				}
				else if ( this->facet[cfi].num_node >3){
					this->facet[cfi].setFacetTypeAsPolygon(this->facet[cfi].num_node);
				}
				else{
					return false;
				}
				pbuf = buf;
				for ( int i=0; i<this->facet[cfi].num_node; i++ ){
					pbuf = strchr(pbuf, ' ');
					pbuf++;
					if ( sscanf_s(pbuf, "%d/%d/%d", &this->facet[cfi].index_node[i], &tmp_int, &this->facet[cfi].index_normal[i] ) != 3 ){
						if ( sscanf_s(pbuf, "%d//%d", &this->facet[cfi].index_node[i], &this->facet[cfi].index_normal[i] ) != 2 ){
							if ( sscanf_s(pbuf, "%d/%d", &this->facet[cfi].index_node[i], &tmp_int ) != 2 ){
								sscanf_s(pbuf, "%d", &this->facet[cfi].index_node[i]);
								this->facet[cfi].normal_type = NORMAL_NONE;
							}
							else{
								this->facet[cfi].normal_type = NORMAL_NONE;
							}
						}
						else{
							this->facet[cfi].normal_type = NORMAL_FACET;
						}
					}
					else{
						this->facet[cfi].normal_type = NORMAL_NODE;
					}
					this->facet[cfi].index_node[i]--;
					this->facet[cfi].index_facet = cfi;
					if ( this->facet[cfi].normal_type!=NORMAL_NONE )
						this->facet[cfi].index_normal[i]--;
				}
				this->facet[cfi].index_material = cmi;
			}
			if(strstr(buf, "usemtl ")){
				sscanf_s(buf, "usemtl %s", &tmp_char, 256);
				for ( int i=0; i<num_material; i++ ){
					if ( _strcmpi(material[i].getName(), tmp_char) == 0 ) cmi = i;
				}
			}
		}
		file.close();
		return true;
	}

	bool ObjMesh::loadMTLFile(const char *filename)
	{
		std::ifstream file;
		int cmi = -1;
		char buf[256];
		int ti=4;
		float tmp_float=0.0f;
		Vector3f temp;
		Material _mat;
		strcpy_s(filename_mtl, filename); 
		file.open(filename, std::ios::in);
		if ( !file.is_open() ){
			return false;
		}
		this->num_material=0;
		this->deleteMaterial();
		while( !file.eof() ){
			file.getline(buf, sizeof(buf));
			if(strstr(buf, "newmtl "))this->num_material++;
		}
		this->newMaterial();
		file.clear();
		file.seekg(0, std::fstream::beg);
		while ( !file.eof() )
		{
			file.getline(buf, sizeof(buf));
			if(strstr(buf, "newmtl ")){
				cmi++;
				sscanf_s(buf, "newmtl %s", this->material[cmi].getName(), 256);
			}
			else if(strstr(buf, "Ka ")){
				sscanf_s(buf, "Ka %f %f %f", &temp.x, &temp.y, &temp.z);
				this->material[cmi].setAmbient(temp);
			}
			else if(strstr(buf, "Kd ")){
				sscanf_s(buf, "Kd %f %f %f", &temp.x, &temp.y, &temp.z);
				this->material[cmi].setColor(temp);
				this->material[cmi].setDiffuse(temp);
			}
			else if(strstr(buf, "Ks ")){
				sscanf_s(buf, "Ks %f %f %f", &temp.x, &temp.y, &temp.z);
				this->material[cmi].setSpecular(temp);
			}
			else if(strstr(buf, "Tr ")){
				if(sscanf_s(buf, "Tr %f %f %f", &temp.x, &temp.y, &temp.z)==3)
					this->material[cmi].setEmission(temp);
			}
			else if(strstr(buf, "Tf ")){
				if(sscanf_s(buf, "Tf %f %f %f", &temp.x, &temp.y, &temp.z)==3)
					this->material[cmi].setEmission(temp);
			}
			else if(strstr(buf, "Ni ")){
				if(sscanf_s(buf, "Ni %f", &temp.x))
					this->material[cmi].setShininess(temp.x);
			}
			else if(strstr(buf, "Ns ")){
				if(sscanf_s(buf, "Ns %f", &temp.x))
					this->material[cmi].setShininess(temp.x);
			}
			else if(strstr(buf, "illum ")){
				if(sscanf_s(buf, "illum %d", &ti))
					this->material[cmi].setIllum(ti);
			}
			else if(strstr(buf, "d ")){
				if(sscanf_s(buf, "d %f", &tmp_float)){
					this->material[cmi].setIllum(5);
					this->material[cmi].setColorA(tmp_float);
					this->material[cmi].setAmbientA(tmp_float);
					this->material[cmi].setDiffuseA(tmp_float);
					this->material[cmi].setSpecularA(tmp_float);
					this->material[cmi].setEmissionA(tmp_float);
				}
			}
		}
		file.close();
		return true;
	}

	bool ObjMesh::load(const char *_filename)
	{
		std::cout << "Loading .obj data ...";
		getDirectoryName(_filename, directoryname);
		if (!loadOBJFile(_filename)) {
			std::cout << "[FAIL]" << std::endl;
			return false;
		}
		this->setup();
		std::cout << "[OK]" << std::endl;
		return true;
	}
	void ObjMesh::setup()
	{
		this->calFacetVertex();
		this->calFacetNormal();
		this->calCenter();
		this->calScale();
		this->calTransform();
		this->calLine();
		this->is_loaded = true;
	}
}