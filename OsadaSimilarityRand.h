/*************************************************************************
    > File Name: OsadaSimilarityRand.h
    > Author: clzhang 
    > Mail: zcwtuibian@gmail.com 
    > Created Time: 2014年11月26日 星期三 05时28分42秒
 ************************************************************************/
 
#include<fstream>
#include<math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <ctype.h>
#include <cmath>
#include <cassert>
#include <vector>
#include <numeric>
#include <gsl/gsl_multifit.h>

#include "Base/3dModelBase.h"
//#include "Base/Distance_vector.h"

using std::string;

//#include "mclmcr.h"
//#include "mclmcrrt.h"
//#include "euclidean.h"

//#include "mclcppclass.h"

// Windows include files 
#ifdef _WIN32
#include <windows.h>
#endif
/*using namespace std;*/

/*const char *filename="m1701.off";
const char *filename1="m1701.off";
int rand1 = 46 * 80;
int rand2 = 46 * 80;*/
//nface = 328

class Area // store tri-face's area, gravity etc.
{ 
	friend class SD;
	friend int  getAreanum(Area *AR, int nface, float nrandarea); //add by levy @2013/10/07 for get the _mesh number according to the area
  private:
	Vertex *verts;
	float area; 
	Vertex gravity;  // store tri-face's gravity
};

class RArea //  store tri-face's random sampling points, and conrresponding vertex, area etc.
{
	friend class SD;
	friend int  cmp(const void *a,const void *b);
  private:
	Vertex *verts;
	float area;
	Vertex Sampling_points;
};

class SD:public ModelBase
{
 public:
	// SD(string filename, int idx): _filename(filename), _idx(idx){}
	// SD(char *filename, int rand, int idx): _filename(filename), _rand(rand), _idx(idx){}
	SD(const char *filename, int rand, int idx, int binnum, int randnum): ModelBase(filename), _rand(rand), _idx(idx), _binnum(binnum), _randnum(randnum){}
	
	~SD(){}

    //int cmp(const void *a,const void *b);
	/*int getAreanum(Area *AR, int nface, float nrandarea); //add by levy @2013/10/07 for get the mesh number according to the area*/
	/*Mesh *ReadOffFile();*/
	Area *CalArea(); // cal each tri-face area
	Area * ARGravity(Area *AR);  // cal each triangle's gravity
	Vertex * ModelGravity(Area *AR); // cal 3d model gravity
	Mesh * GravityNormalize(Vertex * MG); // move the model's gravity to origin
	Mesh * ScaleNormalize(); //ScaleNormalize the 3d model
	Area * CalArea_sort(); // Calculate Area of face, which has been scale normalized.
    void resetRandPointNum(); // set rand point num according to the num of mesh face
    RArea* calPseudoSampling(Area *AR); // Calculate Sampling random points by pseudo-random-sampling (not area weighted)
    void pseudoRandpoints2pcd(RArea *RAR); // output pseudo-randpoints by pcd file
	RArea* calsampling(Area *AR); // Calculate Sampling random points by random-sampling (area weighted, so more uniform)
	void Row2Column(RArea *RAR); // store rand sampling point to a matrix
    void scaledvertex2pcd(Mesh *mesh); // ouput the scaled vertexes to pcd file.
	void Randpoints2pcd(); // ouput the random points to pcd file.
	/*double * Row2Column(RArea *RAR);*/
	//Mesh * Normalize_nopose();
	//void Randsampling();
	void StatDist(); // statistics distance of each pair of sampling points
	void GetFeature(); // gain the coefficients of cubic polynomial and extract random points on the curve.

 private:
	int _rand;//the number of random points, which extracted on the surface of 3d mesh.
	const int _idx;//the index of .off file.
	const int _binnum;//the number of statistics bin.
	const int _randnum;//the number of the random points, which extracted on polyfit curve.
	vector<vector<double> > _mat; //save the coodinate of random points
	/*double _distnum[1000];//num in specific interval,the y to get the polyfit curve //dynamic array initialization in this way is an error. */
	vector<double> _distnum;//num in specific interval,the y to get the polyfit curve //dynamic array initialization in this way is an error. 
	/*double _medX[_binnum];//the x to get the polyfit curve.*/
	vector<double> _medX;//the x to get the polyfit curve.
	double _nsumArea;//the surface area of the model.
	/*double randx[_randnum] = {0.057896,0.084036,0.08974,0.11185,0.12164,0.13897,0.15259,0.18574,0.25386,0.27708,0.29458,0.3184,0.37944,0.38436,0.4369,0.5997,0.68802,0.72718,0.76114,0.79184,0.79644,0.83606,0.8777,0.88664,0.91984,0.94922,0.95262,0.99236,1.0217,1.0488,1.0529,1.1008,1.2062,1.2166,1.27,1.3054,1.3261,1.3324,1.4124,1.4155,1.4478,1.4587,1.4607,1.469,1.4932,1.5144,1.5263,1.5308,1.5675,1.597,1.6196,1.6279,1.6514,1.6612,1.6773,1.6924,1.7239,1.7759,1.8122,1.8317,1.8644,1.9232,1.9287,1.9922};*/
};

// cal each tri-face area
Area * SD::CalArea()
{
     Area *AR=new Area[_mesh->nfaces];
      for(int i = 0;i<_mesh->nfaces;i++)
	  {  
	     AR[i].verts=new Vertex[3];
          for(int j =0;j<3;j++)
			{
		      AR[i].verts[j].x=_mesh->faces[i].verts[j]->x;
		      AR[i].verts[j].y=_mesh->faces[i].verts[j]->y;
		      AR[i].verts[j].z=_mesh->faces[i].verts[j]->z;
             /*  cout<<AR[i].verts[j].x<<" ";
               cout<<AR[i].verts[j].y<<" ";
               cout<<AR[i].verts[j].z<<" "<<endl;*/
			}
	         //cout<<endl;
	  }
  //    ScaleNormalize(AR,_mesh->nfaces);

      for(int i = 0;i<_mesh->nfaces;i++)	//calculate the area of face	
	  {
	    float *temp=new float[3];
		for(int j = 0; j < 3; j++){temp[j]=0; }
		for(int j = 0; j < 3; j++){
	    temp[j] += pow(AR[i].verts[j].x-AR[i].verts[(j+1)%3].x,2);
		temp[j] += pow(AR[i].verts[j].y-AR[i].verts[(j+1)%3].y,2);
		temp[j] += pow(AR[i].verts[j].z-AR[i].verts[(j+1)%3].z,2);
		temp[j] = sqrt(temp[j]);
		}
        //Helen formula
        float p=(temp[0]+temp[1]+temp[2])/2;
        float area=sqrt(p*(p-temp[0])*(p-temp[1])*(p-temp[2])); 
		
		//a little number area is -nan, so that make the _nsumArea is -nan.
		if(area > 1.0E-12 && area == area) AR[i].area=area;
		else AR[i].area=0;
		//if(_isnan(area)) AR[i].area=0;
		 
		//cout<<"i="<<i+1<<":";
    	//cout<<AR[i].area<<endl;
	  }

    // qsort(AR,_mesh->nfaces,sizeof(Area),cmp);//ÅÅÐòAR
return AR;
}	

// cal each triangle's gravity
Area *  SD::ARGravity(Area *AR)  
{
	for (int i=0;i < _mesh->nfaces;i++)
	{
		/*AR[i].gravity.x=(AR[i].verts[j].x+
						AR[i].verts[++j].x+
						AR[i].verts[++j].x)/3;*/
		AR[i].gravity.x=(AR[i].verts[0].x+
						AR[i].verts[1].x+
						AR[i].verts[2].x)/3;
		AR[i].gravity.y=(AR[i].verts[0].y+
						AR[i].verts[1].y+
						AR[i].verts[2].y)/3;
		AR[i].gravity.z=(AR[i].verts[0].z+
						AR[i].verts[1].z+
						AR[i].verts[2].z)/3;
	}

	/*
	for (i=0;i < _mesh->nfaces;i++)
	{
		cout<<"AR["<<i<<"].gravity.x is"<<AR[i].gravity.x<<endl;
		cout<<"AR["<<i<<"].gravity.y is"<<AR[i].gravity.y<<endl;
		cout<<"AR["<<i<<"].gravity.z is"<<AR[i].gravity.z<<endl;
	}
	*/
	return AR;
}

// cal 3d model gravity
Vertex *  SD::ModelGravity(Area *AR) 
{
	Vertex * centergravity=new Vertex[1]; 
	centergravity->x=0;
	centergravity->y=0;
	centergravity->z=0;

	float sumarea=0;
	for (int i=0;i<_mesh->nfaces;i++)
	{
		sumarea+=AR[i].area;
	}

	//cout<<"sumarea:"<<sumarea<<endl;

	for (int i=0;i<_mesh->nfaces;i++)
	{
		centergravity->x+=AR[i].area*AR[i].gravity.x;
		centergravity->y+=AR[i].area*AR[i].gravity.y;
		centergravity->z+=AR[i].area*AR[i].gravity.z;
	}
	centergravity->x /=sumarea;
	centergravity->y /=sumarea;
	centergravity->z /=sumarea;
/*
	cout<<"centergravity.x is"<<centergravity->x<<endl;
	cout<<"centergravity.y is"<<centergravity->y<<endl;
	cout<<"centergravity.z is"<<centergravity->z<<endl;
*/
	return centergravity;
}

// move the model's gravity to origin
Mesh *  SD::GravityNormalize(Vertex * MG)
{
	for (int i=0;i<_mesh->nverts;i++)
	{
		_mesh->verts[i].x-=MG->x;
		_mesh->verts[i].y-=MG->y;
		_mesh->verts[i].z-=MG->z;
	}
	/*
	for (i=0;i<_mesh->nverts;i++)
	{
		cout<<"xŒõÈ¥Ä£ÐÍÖØÐÄ"<<_mesh->verts[i].x<<endl;
		cout<<"yŒõÈ¥Ä£ÐÍÖØÐÄ"<<_mesh->verts[i].y<<endl;
		cout<<"zŒõÈ¥Ä£ÐÍÖØÐÄ"<<_mesh->verts[i].z<<endl;
	}
	*/
	return _mesh;
}

//ScaleNormalize the 3d model
Mesh *  SD::ScaleNormalize()
{  
	float normal_length = 0.0;
	float max =0.0;
	for(int i = 0; i < _mesh->nverts; i++)
	{ 
		normal_length += _mesh->verts[i].x*_mesh->verts[i].x;
		normal_length += _mesh->verts[i].y*_mesh->verts[i].y;
		normal_length += _mesh->verts[i].z*_mesh->verts[i].z;
		if(max < normal_length) max = normal_length;
		normal_length = 0.0;		
	}
	//normal_length = sqrt(normal_length);
	max = sqrt(max);
	for (int i=0;i<_mesh->nverts;i++)
	{
		if (max > 1.0E-6) 
		{
			/*_mesh->verts[i].x  /= normal_length;
			_mesh->verts[i].y  /= normal_length;
			_mesh->verts[i].z  /= normal_length;*/
			_mesh->verts[i].x  /= max;
			_mesh->verts[i].y  /= max;
			_mesh->verts[i].z  /= max;
		}
	}
	//cout<< _mesh->verts[0].x<<" "<< mesh->verts[0].y<<" "<< mesh->verts[0].z<<" "<<endl;
/*
	for (i=0;i<_mesh->nverts;i++)
	{
		cout<< _mesh->verts[i].x<<" "<< mesh->verts[i].y<<" "<< mesh->verts[i].z<<" "<<endl;
	}
*/	 
	return _mesh;
}

// the compare style for qsort function
int  cmp(const void *a,const void *b)
{
	return (*(RArea*)a).area>(*(RArea*)b).area ? 1:-1;
}

// Calculate Area of face, which has been scale normalized.
Area * SD::CalArea_sort()
{
     Area *AR=new Area[_mesh->nfaces];
     // cout<<"ÎŽ±ê×Œ»¯"<<endl;
      for(int i = 0;i<_mesh->nfaces;i++)
	  {  
	     AR[i].verts=new Vertex[3];
          for(int j =0;j<3;j++)
			{
		      AR[i].verts[j].x=_mesh->faces[i].verts[j]->x;
		      AR[i].verts[j].y=_mesh->faces[i].verts[j]->y;
		      AR[i].verts[j].z=_mesh->faces[i].verts[j]->z;
             /*  cout<<AR[i].verts[j].x<<" ";
               cout<<AR[i].verts[j].y<<" ";
               cout<<AR[i].verts[j].z<<" "<<endl;*/
			}
	         //cout<<endl;
	  }
      //Normalize(AR,_mesh->nfaces);

	  _nsumArea = 0;//make the initialize.
      //ofstream ofArea("Area4debug.txt");
      for(int i = 0;i<_mesh->nfaces;i++)	//calculate the area of face	
	  {
	    float *temp=new float[3];
		for(int j = 0; j < 3; j++){temp[j]=0; }
		for(int j = 0; j < 3; j++){
	        temp[j] += pow(AR[i].verts[j].x-AR[i].verts[(j+1)%3].x,2);
		    temp[j] += pow(AR[i].verts[j].y-AR[i].verts[(j+1)%3].y,2);
		    temp[j] += pow(AR[i].verts[j].z-AR[i].verts[(j+1)%3].z,2);
		    temp[j] = sqrt(temp[j]);
		}
        //Helen formula
        double p=(temp[0]+temp[1]+temp[2])/2;
        double area=sqrt(p*(p-temp[0])*(p-temp[1])*(p-temp[2]));

		//a little number area is -nan, so that make the _nsumArea is -nan.
		if(area > 1.0E-12 && area == area) AR[i].area=area;
		else AR[i].area=0;
		//if(_isnan(area)) AR[i].area=0;

    	//cout<<"i="<<i<<":";
    	//cout<<AR[i].area<<" ";
        //ofArea << AR[i].area << endl;

        //_nsumArea += area;       //add by levy @2013/10/27 for cal the total area 
        _nsumArea += AR[i].area;       //add by levy @2013/10/27 for cal the total area 
	    //cout<< _nsumArea << " ";
	  }
      //cout << endl;
	  //cout<<"_nsumArea: "<< _nsumArea << endl;
      //ofArea.close();
      qsort(AR,_mesh->nfaces,sizeof(Area),cmp);//descend sort the area of tripatch.
	return AR;
}

// set rand point num according to the num of mesh face
void SD::resetRandPointNum(){
    if (_mesh->nfaces <= 0){
        cout << "Error: the mesh face num <= 0" <<endl;
        exit(-1);
    }
	if (_mesh->nfaces<1024){
        _rand = _mesh->nfaces;
	} 
	else if(_mesh->nfaces<2048){
	  _rand=1024+(_mesh->nfaces-1024)/2;
	}
	else if(_mesh->nfaces<4096){
	   _rand=1536+(_mesh->nfaces-2048)/4;
	}
	else if(_mesh->nfaces<8192){
		_rand=2048+(_mesh->nfaces-4096)/8;
	}
	else if(_mesh->nfaces<16384){
		_rand=2560+(_mesh->nfaces-8192)/16;
	}
	else if(_mesh->nfaces<32768){
		_rand=3072+(_mesh->nfaces-16384)/32;
	}
	else if(_mesh->nfaces<65536){
		_rand=3584+(_mesh->nfaces-32768)/64;
	}
	else{
		_rand = 4096;
	}		
    //cout << "randpointNum: " << _rand <<endl;
}

// Calculate Sampling random points by pseudo-sampling (not area weighted)
RArea*  SD::calPseudoSampling(Area *AR)//n=mehs-<nfaces
{  
	long i;
    long randidx; 
	RArea *RAR=new RArea[_rand];//return RAR
	srand((unsigned)time(0)); 
	
    //cout << "mesh->faces: " << _mesh->nfaces <<endl;
	for(i=0;i<_rand;i++)
	{
        randidx = rand()%(_mesh->nfaces); // get a random idx of triangle. 

		RAR[i].area=AR[randidx].area;
		RAR[i].verts=AR[randidx].verts;

        //cout << randidx << " " ;
        //ofstream ofRArea("RArea4debug.txt");
        //ofRArea << RAR[i].verts << endl;
        //ofRArea.close();

		float r1=0,r2=0;
		r1=rand()/(float)(RAND_MAX);
		r2=rand()/(float)(RAND_MAX);
		
		Vertex p;//sampling
		p.x=p.y=p.z=0;
		p.x =(1-sqrt(r1))*AR[randidx].verts[0].x;
		p.x += sqrt(r1)*(1-r2)*AR[randidx].verts[1].x;
		p.x += sqrt(r1)*r2*AR[randidx].verts[2].x;
		
		p.y =(1-sqrt(r1))*AR[randidx].verts[0].y;
		p.y += sqrt(r1)*(1-r2)*AR[randidx].verts[1].y;
		p.y += sqrt(r1)*r2*AR[randidx].verts[2].y;
		
		p.z =(1-sqrt(r1))*AR[randidx].verts[0].z;
		p.z += sqrt(r1)*(1-r2)*AR[randidx].verts[1].z;
		p.z += sqrt(r1)*r2*AR[randidx].verts[2].z;
		
		RAR[i].Sampling_points=p;
	}
    //cout << endl;
	/*cout << RAR[0].Sampling_points.x <<"," << RAR[0].Sampling_points.y << "," << RAR[0].Sampling_points.z <<endl;*/
//	qsort(RAR,_rand,sizeof(RArea),cmp);
	return RAR;
}

// output pseudo-randpoints by pcd file
void SD::pseudoRandpoints2pcd(RArea *RAR){
        /*
		char * filenamepcd = new char[30];
		for(int i=0;i<strlen(_filename);i++){
			if(_filename[i] != '.') filenamepcd[i] = _filename[i];
			else{ strcat(filenamepcd, ".pcd"); break;}
		}
        */

        string filenamepcd = _filename;
        //string filenametemp = _filename;

        filenamepcd.replace(strlen(_filename)-4,4,"_pseudoRand.pcd");

    	//cout << filenamepcd <<endl;

		ofstream saveFile(filenamepcd.c_str());

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
		saveFile<<str7<<_rand<<endl;
		saveFile<<str8<<endl;
		saveFile<<str9<<endl;
		saveFile<<str10<<_rand<<endl;
		saveFile<<str11<<endl;
        
        vector<double>::iterator vdi;
		//saveFile<<vertices[0][0]<<" "<<vertices[0][1]<<" "vertices[0][2]<<endl;
		for(int i=0;i<_rand;i++){
            saveFile << RAR[i].Sampling_points.x << " ";
            saveFile << RAR[i].Sampling_points.y << " ";
            saveFile << RAR[i].Sampling_points.z << endl;
		}
		
		saveFile.close();

		//delete [] filenamepcd;
}

// get the random mesh number according to the area
int  getAreanum(Area *AR, int nface, float nrandarea)         
{
    int i;
    float ntempArea = 0;
    for(i = 0;i < nface - 1 ;i++)
    {
        if(ntempArea > nrandarea)
            return i-1;
        else
            ntempArea += AR[i].area;
    }
    return i;
}

// Calculate Sampling random points by random-sampling (area weighted, so more uniform)
RArea*  SD::calsampling(Area *AR)
{  	
	long i,j;
    float nrandarea;
	RArea *RAR=new RArea[_rand];//return RAR
	srand((unsigned)time(0)); 
	
    //cout << "mesh->faces: " << _mesh->nfaces <<endl;
    cout << "randpointNum: " << _rand <<endl;
	for(i=0;i<_rand;i++)
	{
        /* modify by  levy @2013/10/07 for change the way to get a _mesh
		j=1+(int)((_mesh->nfaces-1)*rand()/(RAND_MAX+1.0));
        */

		nrandarea = rand()/(float)(RAND_MAX/_nsumArea);	//get a random area

		j = getAreanum(AR, _mesh->nfaces, nrandarea);

		RAR[i].area=AR[j].area;
		RAR[i].verts=AR[j].verts;

        //cout << j << " " ;
        //ofstream ofRArea("RArea4debug.txt");
        //ofRArea << RAR[i].verts << endl;
        //ofRArea.close();

		float r1=0,r2=0;  
		r1=rand()/(float)(RAND_MAX);
		r2=rand()/(float)(RAND_MAX);
		
		Vertex p;//sampling
		p.x=p.y=p.z=0;
		p.x =(1-sqrt(r1))*AR[j].verts[0].x;
		p.x += sqrt(r1)*(1-r2)*AR[j].verts[1].x;
		p.x += sqrt(r1)*r2*AR[j].verts[2].x;
		
		p.y =(1-sqrt(r1))*AR[j].verts[0].y;
		p.y += sqrt(r1)*(1-r2)*AR[j].verts[1].y;
		p.y += sqrt(r1)*r2*AR[j].verts[2].y;
		
		p.z =(1-sqrt(r1))*AR[j].verts[0].z;
		p.z += sqrt(r1)*(1-r2)*AR[j].verts[1].z;
		p.z += sqrt(r1)*r2*AR[j].verts[2].z;
		
		RAR[i].Sampling_points=p;
	}
        //cout << endl;
	//cout << "RAR[0].Sampling_points.x: " << RAR[0].Sampling_points.x <<"," << RAR[0].Sampling_points.y << "," << RAR[0].Sampling_points.z <<endl;
//	qsort(RAR,_rand,sizeof(RArea),cmp);
    return RAR;
}
	  
// store rand sampling point to a matrix
void SD::Row2Column(RArea *RAR)
{
	/*vector<vector<double> > _mat(_rand);*/
	//_mat.resize(_rand+1, vector<double>(0)); //if _mat.resize(_rand+1, vector<double>(4)); : the first 4 will be 0.(why???)
	_mat.resize(_rand+1); 
	for (int i=0;i<_rand;i++)   
	{
        //cout << _mat[i].size() << endl;
		_mat[i].push_back(RAR[i].Sampling_points.x);
		_mat[i].push_back(RAR[i].Sampling_points.y);
		_mat[i].push_back(RAR[i].Sampling_points.z);
        //cout << _mat[i].size() << endl;
	    //cout << _mat[i].at(0) <<"," << _mat[i].at(1) << "," << _mat[i].at(2) <<endl;
	}
	//cout << rar[0].sampling_points.x <<"," << rar[0].sampling_points.y << "," << rar[0].sampling_points.z <<endl;
	//cout << _mat[0].at(0) <<"," << _mat[0].at(1) << "," << _mat[0].at(2) <<endl;
    
    //vector<double>::iterator vdi = _mat[0].begin();
    //for(;vdi != _mat[0].end();vdi++){
        //cout << *vdi <<endl;
    //}
}

/*
double *  SD::Row2Column(RArea *RAR)
{
	double *TemPtr;
	TemPtr=new double[3*_rand];
	int j=0;
	
	for (int i=0;i<_rand;i++)   
	{
		if (j==0)
		{
			TemPtr[j]=RAR[i].Sampling_points.x;
			TemPtr[++j]=RAR[i].Sampling_points.y;
			TemPtr[++j]=RAR[i].Sampling_points.z;
			continue;
		}
		else
		{
			TemPtr[++j]=RAR[i].Sampling_points.x;
			TemPtr[++j]=RAR[i].Sampling_points.y;
			TemPtr[++j]=RAR[i].Sampling_points.z;
			continue;
		}
	}
	return TemPtr;
}
*/

// output randpoints by pcd file
void SD::Randpoints2pcd(){
        /*
		char * filenamepcd = new char[30];
		for(int i=0;i<strlen(_filename);i++){
			if(_filename[i] != '.') filenamepcd[i] = _filename[i];
			else{ strcat(filenamepcd, ".pcd"); break;}
		}
        */

        string filenamepcd = _filename;
        //string filenametemp = _filename;

        filenamepcd.replace(strlen(_filename)-4,4,"_Rand.pcd");

    	//cout << filenamepcd <<endl;

		ofstream saveFile(filenamepcd.c_str());

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
		saveFile<<str7<<_rand<<endl;
		saveFile<<str8<<endl;
		saveFile<<str9<<endl;
		saveFile<<str10<<_rand<<endl;
		saveFile<<str11<<endl;

        
    vector<double>::iterator vdi;
		//saveFile<<vertices[0][0]<<" "<<vertices[0][1]<<" "vertices[0][2]<<endl;
		for(int i=0;i<_rand;i++){
			//for(int j=0;j<3;j++){
				//saveFile<< *(vertices+j+i*3)<<" ";
            for(vdi=_mat[i].begin();vdi != _mat[i].end();vdi++){
				saveFile<< *vdi  <<" ";
                //cout << *vdi <<endl;
                /*
				saveFile<< _mat[i].at(0) <<" ";
				saveFile<< _mat[i].at(1) <<" ";
				saveFile<< _mat[i].at(2) <<" ";
				cout << _mat[i].at(0) <<" ";
                */
			}
			saveFile<<endl;
		}
		
		saveFile.close();

		//delete [] filenamepcd;
}

// statistics distance of each pair of sampling points
void SD::StatDist(){
		vector<double> distvec;//distance vector between each two different points
		for(int i=0;i<_rand-1;i++){
			for(int j=i+1;j<_rand;j++){
				/*cout << distVector(_mat[i], _mat[j], "Oshi") << " ";*/
				distvec.push_back(distVector(_mat[i], _mat[j], "Oshi"));
			}
		}
		cout << distvec.size() <<endl;
		cout << *max_element(distvec.begin(), distvec.end()) <<endl;//output the max distance between each two different points.*/

        /*double distnum[1000] = {0};*/
		/*int distnum[distvec.size()];*/
        /*_distnum.resize(_binnum);*/
		/*vector<double>::iterator itnum;*/
		/*for(int itnum=_distnum.begin();itnum != _distnum.end();itnum++) *itnum =0;//dynamic array initialization in this way is right.*/
		/*for(int i=0;i<_binnum;i++) _distnum[i] = 0;//dynamic array initialization in this way is right.*/
		for(int i=0;i<_binnum;i++) _distnum.push_back(0);//dynamic array initialization in this way is right.

        // stat distance num of bin by hash algorithm
		vector<double>::iterator iter;
		for(iter=distvec.begin();iter != distvec.end();iter++){
			/*cout << floor(*iter*1000/2) <<" ";*/
			/*_distnum[(int)(floor(*iter*1000/2))]++;*/
            /*cout << __FILE__ <<endl;*/
	        //double &tempnum = _distnum.at((int)(floor(*iter*1000/2)));
	        double &tempnum = _distnum.at((int)(floor(*iter*_binnum/2)));
            tempnum++;
		}
        /*_distnum(distnum, distnum+_binnum);*/
        /*_distnum(distnum, distnum+1000);*/
		/*cout << *max_element(_distnum, _distnum+_binnum) <<endl;*/
		//cout << _distnum.size() <<endl;
		/*cout << accumulate(_distnum.begin(), _distnum.end(), 0) <<endl;		*/
		/*cout << *max_element(distnum, distnum+binnum) <<endl;*/

        _medX.resize(_binnum);
		/*medX[0] = (2-0)/(binnum*2);*/ //why the medX is zero????????????????
		_medX[0] = 1; //to avoid the medX[0] is zero, I amplify 1000 time of the medX 
		float interval = 2;
		/*cout << _medX[0] <<endl;*/
		for(int i=1;i<_binnum;i++){
			_medX[i] = _medX[i-1] +interval;
		}
		/*cout << medX[0] <<" "<< medX[1] <<endl;*/
		/*ofstream wf("medX_distnum.txt");*/
}	


// gain the coefficients of cubic polynomial and extract random points on the curve.
void SD::GetFeature(){
  int dimension = 4;
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *c;

  X = gsl_matrix_alloc (_binnum, dimension);
  y = gsl_vector_alloc (_binnum);

  c = gsl_vector_alloc (dimension);
  cov = gsl_matrix_alloc (dimension, dimension);

  for (int i = 0; i < _binnum; i++)
    {
	  //cout << _medX[i] <<" "<<_distnum[i] <<endl;
	  //wf << _medX[i] <<" "<<_distnum[i] <<endl;
      
      gsl_matrix_set (X, i, 0, 1.0);
      gsl_matrix_set (X, i, 1, _medX[i]);
      gsl_matrix_set (X, i, 2, _medX[i]*_medX[i]);
      gsl_matrix_set (X, i, 3, _medX[i]*_medX[i]*_medX[i]);
      
      //gsl_vector_set (y, i, _distnum[i]);
      //gsl_vector_set (y, i, _distnum.at(i));
    }
	//wf.close();

    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (_binnum, dimension);
    gsl_multifit_linear (X, y, c, cov, &chisq, work); //polyfit function
    gsl_multifit_linear_free (work);

	#define C(i) (gsl_vector_get(c,(i)))
	#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  
    printf ("# best fit: Y = %g + %g X + %g X^2+ %g X^3\n", 
            C(0), C(1), C(2), C(3));

    //printf ("# covariance matrix:\n");
    //printf ("[ %+.5e, %+.5e, %+.5e  \n",
               //COV(0,0), COV(0,1), COV(0,2));
    //printf ("  %+.5e, %+.5e, %+.5e  \n", 
               //COV(1,0), COV(1,1), COV(1,2));
    //printf ("  %+.5e, %+.5e, %+.5e ]\n", 
               //COV(2,0), COV(2,1), COV(2,2));
    //cout << chisq <<endl;

	//double randx[_randnum] = {0.057896,0.084036,0.08974,0.11185,0.12164,0.13897,0.15259,0.18574,0.25386,0.27708,0.29458,0.3184,0.37944,0.38436,0.4369,0.5997,0.68802,0.72718,0.76114,0.79184,0.79644,0.83606,0.8777,0.88664,0.91984,0.94922,0.95262,0.99236,1.0217,1.0488,1.0529,1.1008,1.2062,1.2166,1.27,1.3054,1.3261,1.3324,1.4124,1.4155,1.4478,1.4587,1.4607,1.469,1.4932,1.5144,1.5263,1.5308,1.5675,1.597,1.6196,1.6279,1.6514,1.6612,1.6773,1.6924,1.7239,1.7759,1.8122,1.8317,1.8644,1.9232,1.9287,1.9922};
	double randx[] = {0.057896,0.084036,0.08974,0.11185,0.12164,0.13897,0.15259,0.18574,0.25386,0.27708,0.29458,0.3184,0.37944,0.38436,0.4369,0.5997,0.68802,0.72718,0.76114,0.79184,0.79644,0.83606,0.8777,0.88664,0.91984,0.94922,0.95262,0.99236,1.0217,1.0488,1.0529,1.1008,1.2062,1.2166,1.27,1.3054,1.3261,1.3324,1.4124,1.4155,1.4478,1.4587,1.4607,1.469,1.4932,1.5144,1.5263,1.5308,1.5675,1.597,1.6196,1.6279,1.6514,1.6612,1.6773,1.6924,1.7239,1.7759,1.8122,1.8317,1.8644,1.9232,1.9287,1.9922};
	//double randy[_randnum];
    vector<double> randy;

	ofstream wf("randy.txt", ios::app);
	//ofstream wf("randy.txt", ios::ate);
	for(int i=0;i<_randnum;i++){
		randy.push_back(C(0) + C(1)*randx[i] + C(2)*randx[i]*randx[i] + C(3)*randx[i]*randx[i]*randx[i]);
		//cout << randy[i] << " ";
		wf << randy.at(i) <<" ";
	}
	wf <<endl;
	wf.close();

  gsl_matrix_free (X);
  gsl_vector_free (y);
  gsl_vector_free (c);
  gsl_matrix_free (cov);
}

Mesh * Normalize_nopose(SD & sd){
	static Mesh *mesh=NULL; 
	static Area *AR=NULL;
	static Vertex *centergravity=NULL;

	//AR= SD::CalArea();   
	//AR=SD::ARGravity(AR);  
	//centergravity=SD::ModelGravity(AR); 
	//mesh=SD::GravityNormalize(centergravity);
	//mesh=SD::ScaleNormalize(); 
	AR= sd.CalArea();   // cal each triangle face area.
	AR=sd.ARGravity(AR);  //cal each tripatch gravity.
	centergravity=sd.ModelGravity(AR);  //cal model's gravity.
	mesh=sd.GravityNormalize(centergravity); //translate to origin.
	mesh=sd.ScaleNormalize();  // scaled the model to a unit sphere.
	
	delete(AR);
	delete(centergravity);

    return mesh;
}

// ouput scaled vertex by pcd file directly
void SD::scaledvertex2pcd(Mesh *mesh){
        /*
		char * filenamepcd = new char[30];
		for(int i=0;i<strlen(_filename);i++){
			if(_filename[i] != '.') filenamepcd[i] = _filename[i];
			else{ strcat(filenamepcd, ".pcd"); break;}
		}
        */

        string filenamepcd = _filename;
        //string filenametemp = _filename;

        filenamepcd.replace(strlen(_filename)-4,4,"_scaledvertex.pcd");

    	//cout << filenamepcd <<endl;

		ofstream saveFile(filenamepcd.c_str());

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
		saveFile<<str7<<mesh->nverts<<endl;
		saveFile<<str8<<endl;
		saveFile<<str9<<endl;
		saveFile<<str10<<mesh->nverts<<endl;
		saveFile<<str11<<endl;

        
    vector<double>::iterator vdi;
		//saveFile<<vertices[0][0]<<" "<<vertices[0][1]<<" "vertices[0][2]<<endl;
		for(int i=0;i<mesh->nverts;i++){
				saveFile<< mesh->verts[i].x <<" "<< mesh->verts[i].y <<" "<< mesh->verts[i].z <<endl;
		}
		
		saveFile.close();

		//delete [] filenamepcd;
}

