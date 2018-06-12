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
		int temp_countFacet=0;
		std::ifstream file;
		Facet sample;
		char buf[256];
		char dummy[256];
		char dummy2[256];
		file.open(_filename, std::ios::in);
		if (!file.is_open())return false;
		this->num_node=0;
		this->num_facet=0;
		while (!file.eof()){
			file.getline(buf, sizeof(buf));
//			if(strstr(buf, "vertex"))this->num_node++;
			if(strstr(buf, "normal"))this->num_facet++;
		}
		this->num_normal = this->num_facet;
		this->num_node = this->num_facet * 3;
		this->num_material=1;
		this->deleteMesh();
		this->deleteMaterial();
		this->newMesh();
		this->newMaterial();
		sample.setFacetTypeAsTriangle();
		file.clear();
		file.seekg(0, std::fstream::beg);
		while (!file.eof()){
			file.getline(buf, sizeof(buf));
			if (strstr(buf, "normal")) {
				sscanf_s(buf, "%s %s %lf %lf %lf", dummy, 256, dummy2, 256, &sample.normal[0].x, &sample.normal[0].y, &sample.normal[0].z);
			}
			else if(strstr(buf, "vertex")){
					sscanf_s(buf, "%s %lf %lf %lf", dummy, 256, &sample.vertex[0].x,&sample.vertex[0].y,&sample.vertex[0].z);
					file.getline(buf, sizeof(buf));
					sscanf_s(buf, "%s %lf %lf %lf", dummy, 256, &sample.vertex[1].x,&sample.vertex[1].y,&sample.vertex[1].z);
					file.getline(buf, sizeof(buf));
					sscanf_s(buf, "%s %lf %lf %lf", dummy, 256, &sample.vertex[2].x,&sample.vertex[2].y,&sample.vertex[2].z);
					sample.index_facet=temp_countFacet;
					sample.normalType=NORMAL_Facet;
					this->facet[temp_countFacet]=sample;
					sample.index_normal[0]=temp_countFacet;
					sample.index_normal[0]=temp_countFacet;
					sample.index_normal[0]=temp_countFacet;
					normal[temp_countFacet]=sample.normal[0];
					temp_countFacet++;
			}
		}
		file.close();
		this->setup();
		return true;
	}
	bool StlMesh::loadBinary(const char* _filename)
	{
		int count=0;
		float *pf;
		int *pd;
		std::ifstream file;
		Facet sample;
		char buf[256];
		file.open(_filename, std::ios::in|std::ios::binary);
		if (!file.is_open())return false;
		file.read(buf,84);
		pd = (int *)&buf[80];
		this->num_facet=*pd;
		this->num_normal=this->num_facet;
		this->num_node=3*this->num_facet;
		this->num_material = 1;
		this->deleteMesh();
		this->deleteMaterial();
		this->newMesh();
		this->newMaterial();
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
				file<<"			vertex "<<this->facet[i].vertex[j].x<<" "<<this->facet[i].vertex[j].y<<" "<<this->facet[i].vertex[j].z<<std::endl;
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
					temp_data=this->facet[i].vertex[j].X[k];
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
		this->calVertex();
		this->calFacetNormal();
		this->calCenter();
		this->calScale();
		this->calTransform();
		this->calLine();
		this->is_loaded = true;
	}
};
