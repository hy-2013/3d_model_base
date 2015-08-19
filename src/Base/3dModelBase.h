/*************************************************************************
    > File Name: ModelBase.h
    > Author: clzhang 
    > Mail: zcwtuibian@gmail.com 
    > Created Time: 2014年11月25日 星期二 05时28分42秒
 ************************************************************************/

#include<iostream>
#include<vector>
#include<cstring>
#include<algorithm>
#include<Eigen/Dense>

using namespace std;
/*#include <fstream>*/
/*#include <string.h>*/
/*#include <stdio.h>*/

class Face;
class Mesh;
class ModelBase;

// TYPE DEFINITION
class Vertex
{
	friend class Face;
	friend class Mesh;
	friend class SD;
	friend class ModelBase;
	private:
	float x,y,z;
	/*double x,y,z;*/ //if statement x,y,z to double type, and use the %f way to read .off, it will be a unknown error.
};
class Face //定义一个三角面片
{ 
	friend class Mesh;
	friend class SD;
	friend class ModelBase;
	private:
	Face(void): nverts(0),verts(0) {};//初始化
	int nverts;//一个面内有几个顶点；
	Vertex **verts;//一个顶点指针数组
};
class Mesh //一个对模型应一个mesh
{  
	friend class SD;
	friend class ModelBase;
  private:
	/*public:*/
	Mesh(void) : nverts(0), verts(0), nfaces(0), faces(0) {};//构造函数初始化
	int nverts; //顶点个数
	Vertex *verts; //顶点坐标
	int nfaces; //面片个数
	Face *faces; //面片顶点序号
};

/*class Area//中间结构体，用来存储每个面片的面积*/
/*{ */
	/*friend class ModelBase;*/
	/*private:*/
	/*Vertex *verts;*/
	/*float area; */
	/*Vertex gravity;  //三角面片重心*/
/*};*/

class ModelBase
{
 public:
	//ModelBase(string filename, int idx): _filename(filename), _idx(idx){}
	ModelBase(const char *filename): _filename(filename){}
	~ModelBase(){}

	void ReadOffFile();
    int getVertexNum();
    int getFaceNum();
	void off2pcd(); // convert vertex of .off to .pcd file directly.
     
    template<typename elemtype>
    double distVector(const vector<elemtype> v1, const vector<elemtype>  v2,  const char distType[]); //calculate the distance of the pair vector.
	/*Mesh *ReadOffFile();*/

 protected: 
	//string _filename;
	const char * _filename;
	Mesh * _mesh; // store the model's vertexes and face.
};

////////////////////////////////////////////////////////////
// OFF FILE READING CODE
////////////////////////////////////////////////////////////
void ModelBase::ReadOffFile()
{
	int i;
    // Open file
    FILE *fp;
	if (!(fp = fopen(_filename, "r"))) 
	{//"r"是以只读方式打开，此参数为const char* 
	fprintf(stderr, "Unable to open file %s\n", _filename);//stderr是标准错误流，一般把屏幕设为默认, 也可以输出到文件
        //return 0;
		exit(-1);
	}
   
    // Allocate _mesh structure
    _mesh = new Mesh();
    if (!_mesh)
	{
		fprintf(stderr, "Unable to allocate memory for file %s\n", _filename);
        fclose(fp);
        //return 0;
		exit(-1);
	}

    // Read file
    int nverts = 0;//顶点数
    int nfaces = 0;//面片数
    int nedges = 0;//边数
    int line_count = 0;//.off文件里有多少行
    char buffer[1024];
    while (fgets(buffer, 1023, fp))
	{
	// Increment line counter
		line_count++;
   
	// Skip white space
		char *bufferp = buffer;//buffer[1024]的头指针
		while (isspace(*bufferp)) bufferp++;//#include<ctype.h>里，若参数c为空格字符，则返回TRUE，否则返回NULL(0)。

	// Skip blank lines and comments
        if (*bufferp == '#') continue;
        if (*bufferp == '\0') continue;


		  // Check section
		if (nverts == 0) 
		{ //当顶点数为0时，也就是第一次读取时
            // Read header 
			if (!strstr(bufferp, "OFF"))
			{//off文件头有OFF字样，strstr返回指向第一次出现OFF位置的指针，如果没找到则返回NULL。
                // Read _mesh counts
				if ((sscanf(bufferp, "%d%d%d", &nverts, &nfaces, &nedges) != 3) || (nverts == 0)) 
				{
					//OFF下面有3个数，分别是点的个数，面片的个数，边的个数。
					fprintf(stderr, "Syntax error reading header on line %d in file %s\n", line_count, _filename);
					fclose(fp);
					//return 0;
					exit(-1);
				}
					// Allocate memory for _mesh
				_mesh->verts = new Vertex [nverts];//申请刚读入的nverts个点的数组
				assert(_mesh->verts);//在<assert.h>中,其作用是如果它的条件返回错误，则终止程序执行
	
				_mesh->faces = new Face [nfaces];//申请读入的nface个FACE数组
				assert(_mesh->faces);
			}
		 }
		 else if (_mesh->nverts < nverts) //当nvert不为0，也就是读取完OFF下面的三个点之后，当mesh里的nvert比读取进来的nvert小
		 {
			// Read vertex coordinates
			Vertex& vert = _mesh->verts[_mesh->nverts];//mesh->verts是一个数组，mesh->nverts从0开始，只要小于读入的nvert
			_mesh->nverts++;                          //mesh中的verts指针，指向这个vertex结构体
			if (sscanf(bufferp, "%f%f%f", &(vert.x), &(vert.y), &(vert.z)) != 3)//读入vert(vertex的结构体中的第_mesh->nvert个位置）
			{
				fprintf(stderr, "Syntax error with vertex coordinates on line %d in file %s\n", line_count, _filename);
				fclose(fp);
				//return 0; 
				exit(-1);
			}
		 }
         
		 else if (_mesh->nfaces < nfaces) 
		 {
			// Get next face
			Face& face = _mesh->faces[_mesh->nfaces];//新申请一个face,然后让mesh中的指针指向他，之后对face操作
			_mesh->nfaces++;

			// Read number of vertices in face 
			bufferp = strtok(bufferp, " \t");//strtok在bufferp中查找包含在'\t'的字符，返回之前的部分（此buffer[1024]被修改）
          
			if (bufferp) 
			{ 
				face.nverts = atoi(bufferp);
			}//把字符串转换成整型数,读取的是off文件中面片的N的数字。
            else
			{
				fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, _filename);
				fclose(fp);
				//return 0;
				exit(-1);
			}
			// Allocate memory for face vertices
			face.verts = new Vertex *[face.nverts];//为一个面中的n个顶点申请一个vertx*的数组，每一个vertex* 抓着一个vertx的结构体
			assert(face.verts);

			// Read vertex indices for face
			for (i = 0; i < face.nverts; i++)
			{
				bufferp = strtok(NULL, " \t");//接着读取N个面片个数之后的 顶点信息。
				if (bufferp) 
					face.verts[i] = &(_mesh->verts[atoi(bufferp)]);//buff读入的值，是mesh->verts数组的下标，代表第几个点，把点的坐标，放入verts[i];
                else 
				{
					fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, _filename);
					fclose(fp);
					//return 0;
					exit(-1);
				}
			 }
		}
	}//end of while

	// Check whether read all faces
	if (nfaces != _mesh->nfaces) 
	{
		fprintf(stderr, "Expected %d faces, but read only %d faces in file %s\n", nfaces, _mesh->nfaces, _filename);
	}

	// Close file
	fclose(fp);
	//cout << _mesh->verts[0].x <<endl;
	//return 0;
	/*return _mesh;		*/
}//end of readofffile  

int ModelBase::getVertexNum(){
    int vertexNum = _mesh->nverts;
    return vertexNum;
}

int ModelBase::getFaceNum(){
    int faceNum = _mesh->nfaces;
    return faceNum;
}

void ModelBase::off2pcd(){
	//int verts;
	//double *vertices;
	//vertices=new double[3*rand2];
	//vertices=sd.Row2Column(RAR1);
	//char filenamepcd[20];
	//string * filenamepcd = (string*)_filename;
	//string_replace(filenamepcd, "off", "pcd");
	//cout << __FILE__ <<endl;	
		
		char * filenamepcd = new char[100];
		for(int i=0;i<strlen(_filename);i++){
			if(_filename[i] != '.') filenamepcd[i] = _filename[i];
			else{ filenamepcd[i] = '\0'; strcat(filenamepcd, ".pcd"); break;}
		}
	//cout << filenamepcd <<endl;
		//ofstream saveFile("cup2.pcd"); //可以读取任意类型的文本文件，远不止txt等
		//ofstream saveFile("Airplane4_gravity.pcd");
		ofstream saveFile(filenamepcd);

		string str1="# .PCD v.7 - Point Cloud Data file format";
		string str2="VERSION .7";
		string str3="FIELDS x y z";
		string str4="SIZE 4 4 4";
		string str5="TYPE F F F";
		string str6="COUNT 1 1 1";
		string str7="WIDTH ";
		string str8="HEIGHT 1";
		string str9="VIEWPOINT 0 0 0 1 0 0 0";
		string str10="POINTS ";
		string str11="DATA ascii";
		
		saveFile << str1 << endl;
		saveFile<<str2<<endl;
		saveFile<<str3<<endl;
		saveFile<<str4<<endl;
		saveFile<<str5<<endl;
		saveFile<<str6<<endl;
		saveFile<<str7<<_mesh->nverts<<endl;
		saveFile<<str8<<endl;
		saveFile<<str9<<endl;
		saveFile<<str10<<_mesh->nverts<<endl;
		saveFile<<str11<<endl;

		//saveFile<<vertices[0][0]<<" "<<vertices[0][1]<<" "vertices[0][2]<<endl;
		for(int i=0;i<_mesh->nverts;i++){
			//for(int j=0;j<3;j++){
				//saveFile<< *(vertices+j+i*3)<<" ";
				saveFile<< _mesh->verts[i].x <<" ";
				saveFile<< _mesh->verts[i].y <<" ";
				saveFile<< _mesh->verts[i].z <<" ";
				//cout << _mesh->verts[i].x <<" ";
			//}
			saveFile<<endl;
		}
		
		saveFile.close();
		
		delete [] filenamepcd;
}


/**
Implement four similarity distance calculation,as
Euclidean, Manhaton, Chebyshev, Mahalanobis distance.
**/

// 实现欧氏距离的计算
template<typename elemtype>
double Oshi(const vector<elemtype> v1, const vector<elemtype>  v2)
{
	typename vector<elemtype>::const_iterator it1 = v1.begin(); //need ‘typename’ before ‘std::vector<elemtype>::const_iterator’ because ‘std::vector<elemtype>’ is a dependent scope
	/*vector<double>::iterator it1 = v1.begin();*/
	typename vector<elemtype>::const_iterator it2 = v2.begin();
	double ODistance=0;
	for(;it1 != v1.end();it1++, it2++)
	{
		ODistance=ODistance+(*it1 - *it2)*(*it1 - *it2);
		/*cout << ODistance << endl;*/
	}
	/*ODistance=ODistance/v1.size();*/
	ODistance=sqrt(ODistance);
	return ODistance;
}

// 实现曼哈顿距离的计算
template<typename elemtype>
double Manhaton(const vector<elemtype> &v1, const vector<elemtype> &v2)
{
	typename vector<elemtype>::const_iterator it1 = v1.begin(); //need ‘typename’ before ‘std::vector<elemtype>::const_iterator’ because ‘std::vector<elemtype>’ is a dependent scope
	typename vector<elemtype>::const_iterator it2 = v2.begin();
	double ODistance=0;
	for(;it1 != v1.end();it1++, it2++)
	{
		ODistance=ODistance+ abs(*it1 - *it2);
		/*cout << ODistance << endl;*/
	}
	return ODistance;
}

// 实现切比雪夫距离的计算
template<typename elemtype>
double Chebyshev(const vector<elemtype> &v1, const vector<elemtype> &v2)
{
	typename vector<elemtype>::const_iterator it1 = v1.begin(); //need ‘typename’ before ‘std::vector<elemtype>::const_iterator’ because ‘std::vector<elemtype>’ is a dependent scope
	typename vector<elemtype>::const_iterator it2 = v2.begin();
	double ODistance = 0;
	for(;it1 != v1.end();it1++, it2++)
	{
		/*cout << *it1 << "," << *it2 <<endl;*/
		/*cout << abs(*it1 - *it2) <<endl;*/
		ODistance = max(ODistance, abs(*it1-*it2)); 
		/*cout << ODistance << endl;*/
	}
	return ODistance;
}

// calculate the covariance of the v1 and v2 
void cov(const vector<double>& v1, const vector<double>& v2, vector<vector<double> > &ret)
{
    assert(v1.size() == v2.size() && v1.size() > 1);
     const int vs = v1.size();
     /*double ret[vs][vs];*/
     double va[vs];
	 for(int i=0;i<vs;i++){
		va[i] = (v1.at(i) + v2.at(i))/2;
	}
 
     for (vector<double>::size_type i = 0; i != v1.size(); ++i)
    {
		for (vector<double>::size_type j = 0; j != v1.size(); ++j)
			ret[i][j]=(v1[i]-va[i])*(v1[j]-va[j])+(v2[i]-va[i])*(v2[j]-va[j]);
            /*ret += (v1[i] - v1a) * (v2[i] - v2a);*/
    }

    /*return ret;*/
}

//实现马氏距离的计算
double Mashi(const vector<double>& v1, const vector<double>& v2)
{   
	const int elemnum = v1.size();
	/*double **cov = cov(dv3, dv4);*/
	vector<vector<double> > covdv; 
	covdv.resize(elemnum+1,vector<double>(elemnum+1)); //apply memory to more, as elemnum+1 (not elemnum), to avoid out of memory 
	cov(v1, v2, covdv);
	
	Eigen::MatrixXd mat(elemnum,elemnum);
	for(int i = 0;i<elemnum;i++){  
		for(int j = 0;j<elemnum;j++){  
			/*cout<< covdv[i][j] << " ";*/
			mat(i,j) = covdv[i][j];
		}
		/*cout << endl;*/
	}

	Eigen::VectorXd ev1(elemnum); 
	Eigen::VectorXd ev2(elemnum); 
	for(int i = 0;i<elemnum;i++){  
		ev1(i) = v1[i];	
		ev2(i) = v2[i];	
	}

	/*cout<< m <<endl;*/
	/*cout<< ev1 <<endl<<endl;*/
	/*cout<< ev2 <<endl<<endl;*/
	double ODistance = 0;
	if(!mat.determinant())
	cout << "the covariance matrix's determinant is zero." << endl;
	else{
		/*cout << (ev1 - ev2).transpose() <<endl<<endl;*/
		/*cout << mat.inverse() <<endl<<endl;*/
		ODistance = sqrt((ev1 - ev2).transpose() * mat.inverse() * (ev1 - ev2)); 
		/*ODistance = (ev1 - ev2).transpose() * m.inverse() * (ev1 - ev2); */
	}
	/*cout << ODistance <<endl;*/
	
	return ODistance;
} 

template<typename elemtype>
double ModelBase::distVector(const vector<elemtype> v1, const vector<elemtype>  v2,  const char distType[]){
/*double distVector(const vector<double> v1, const vector<double>  v2, char* distType){*/
	double Dist = 0;

	if(!strcmp(distType,"Oshi")) Dist = Oshi(v1, v2);
	/*if(distType == "Oshi") OshiDist = Oshi(v1, v2);*/
	else if(!strcmp(distType, "Manhaton")) Dist = Manhaton(v1, v2);
	else if(!strcmp(distType, "Chebyshev")) Dist = Chebyshev(v1, v2);
	/*else if(distType == (string)"Mashi") Dist = Mashi(v1, v2);*/
	else if(!strcmp(distType,"Mashi")) Dist = Mashi(v1, v2);
	/*else if(strcmp(distType, "Cosine")) Dist = Cosine(v1, v2);*/
	
	return Dist;
}
