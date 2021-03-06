#include "stlhandler.h"

namespace EITS
{

	bool StlMesh::load(const char *_filename)
	{
		this->is_loaded = false;
		std::cout <<"Loading .stl data ...";
		std::ifstream file;
		char buf[256];
		int count=0;
		if(_filename==NULL){
			std::cout << "[FAIL]" << std::endl;
			return false;
		}
		file.open(_filename, std::ios::in);
		if (!file.is_open()){
			std::cout << "[FAIL]" << std::endl;
			return false;
		}
		clearMesh();
		while(count<10){
			file.getline(buf, sizeof(buf));
			if(strstr(buf, "facet")){
				file.close();
				this->is_loaded = loadASCII(_filename);
			}
			else count++;
		}
		if(!this->is_loaded)
			this->is_loaded = loadBinary(_filename);
		if(is_loaded)
			std::cout <<"[OK]"<< std::endl;
		else std::cout <<"[FAIL]"<< std::endl;
		return is_loaded;
	}


	bool StlMesh::loadASCII(const char* _filename)
	{
		int temp_count_facet = 0;
		std::ifstream file;
		Facet t_facet;
		char buf[256];
		char dummy[256];
		char dummy2[256];
		file.open(_filename, std::ios::in);
		if (!file.is_open())return false;
		this->num_node=0;
		this->num_facet=0;
		while (!file.eof()){
			file.getline(buf, sizeof(buf));
			if(strstr(buf, "normal"))this->num_facet++;
		}
		this->num_node = this->num_facet * 3;
		this->num_material=1;
		this->deleteMesh();
		this->deleteMaterial();
		this->newMesh();
		this->newMaterial();
		t_facet.setFacetTypeAsTriangle();
		file.clear();
		file.seekg(0, std::fstream::beg);
		while (!file.eof()){
			file.getline(buf, sizeof(buf));
			if (strstr(buf, "normal")) {
				sscanf_s(buf, "%s %s %lf %lf %lf", dummy, 256, dummy2, 256, &t_facet.normal[0].x, &t_facet.normal[0].y, &t_facet.normal[0].z);
			}
			else if(strstr(buf, "vertex")){
				sscanf_s(buf, "%s %lf %lf %lf", dummy, 256, &t_facet.position[0].x, &t_facet.position[0].y, &t_facet.position[0].z);
				file.getline(buf, sizeof(buf));
				sscanf_s(buf, "%s %lf %lf %lf", dummy, 256, &t_facet.position[1].x, &t_facet.position[1].y, &t_facet.position[1].z);
				file.getline(buf, sizeof(buf));
				sscanf_s(buf, "%s %lf %lf %lf", dummy, 256, &t_facet.position[2].x, &t_facet.position[2].y, &t_facet.position[2].z);
				t_facet.index_facet = temp_count_facet;
				t_facet.normal_type = NORMAL_FACET;
				t_facet.normal[0] = temp_count_facet;
				this->facet[temp_count_facet] = t_facet;
				temp_count_facet++;
			}
		}
		file.close();
		this->setup();
		return true;
	}
	bool StlMesh::loadBinary(const char* _filename)
	{
		int count = 0;
		float *pf;
		int *pd;
		std::ifstream file;
		Facet t_facet;
		char buf[256];
		file.open(_filename, std::ios::in|std::ios::binary);
		if (!file.is_open())return false;
		file.read(buf,84);
		pd = (int *)&buf[80];
		this->num_facet = *pd;
		this->num_node = 3 * this->num_facet;
		this->num_material = 1;
		this->deleteMesh();
		this->deleteMaterial();
		this->newMesh();
		this->newMaterial();
		t_facet.setFacetTypeAsTriangle();
		while (file.read(buf,50)){
			pf = (float *)buf;
			t_facet.normal[0].x = *(pf + 0);
			t_facet.normal[0].y = *(pf + 1);
			t_facet.normal[0].z = *(pf + 2);
			t_facet.position[0].x = *(pf + 3);
			t_facet.position[0].y = *(pf + 4);
			t_facet.position[0].z = *(pf + 5);
			t_facet.position[1].x = *(pf + 6);
			t_facet.position[1].y = *(pf + 7);
			t_facet.position[1].z = *(pf + 8);
			t_facet.position[2].x = *(pf + 9);
			t_facet.position[2].y = *(pf + 10);
			t_facet.position[2].z = *(pf + 11);
			t_facet.index_facet = count;
			t_facet.index_normal[0] = count;
			facet[count] = t_facet;
			count++;
		}

		file.close();
		this->setup();
		return true;
	}
	bool StlMesh::save(const char *_filename, int _type)
	{
		bool ret =false;
		std::cout <<"Saving .stl data ...";
		if(_type==STL_ASCII)
			ret = saveASCII(_filename);
		if(_type==STL_BINARY)
			ret = saveBinary(_filename);

		if(ret)
			std::cout <<"[OK]"<< std::endl;
		else std::cout <<"[FAIL]"<< std::endl;
		return ret;
	}

	bool StlMesh::saveASCII(const char* _filename)
	{
		std::ofstream file;
		file.open(_filename, std::ios::out);
		if (!file.is_open())return false;
		file<<"solid "<<_filename<<" was written by ascii code."<<std::endl;
		for(int i=0;i<this->num_facet;i++){
			file<<"	facet normal "<<this->facet[i].normal[0].x<<" "<<this->facet[i].normal[0].y<<" "<<this->facet[i].normal[0].z<<std::endl;
			file<<"		outer loop"<<std::endl;		
			for(int j=0;j<3;j++)
				file<<"			vertex "<<this->facet[i].position[j].x<<" "<<this->facet[i].position[j].y<<" "<<this->facet[i].position[j].z<<std::endl;
			file<<"		endloop"<<std::endl;		
			file<<"	endfacet"<<std::endl;		
		}
		file<<"endsolid"<<std::endl;
		file.close();
		this->is_loaded=true;
		return true;
	}

	bool StlMesh::saveBinary(const char* _filename)
	{
		char *buf=new char[256];
		float temp_data;
		std::ofstream file;
		file.open(_filename, std::ios::out|std::ios::binary);
		if (!file.is_open())return false;
		for(int i=0;i<80;i++)buf[i]=NULL;
		sprintf_s(buf, 256, "Solid %s was written by binary code\n",_filename);
		file.write(buf,80);
		file.write((char*)&this->num_facet,4);
		for(int i=0;i<this->num_facet;i++){
			for(int j=0;j<3;j++){
				temp_data=this->facet[i].normal[0].X[j];
				file.write((char*)&temp_data,4);
			}
			for(int j=0;j<3;j++){
				for(int k=0;k<3;k++){
					temp_data=this->facet[i].position[j].X[k];
					file.write((char*)&temp_data,4);
				}
			}
			file.write("  ",2);
		}
		file.close();
		return true;
	}
	void StlMesh::setup()
	{
		this->calNode();
		this->calFacetNormal();
		this->calCenter();
		this->calScale();
		this->calTransform();
		this->calLine();
		this->is_loaded = true;
	}
};
