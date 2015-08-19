/*************************************************************************
    > File Name: OsadaSimilarityRand.cpp
    > Author: clzhang 
    > Mail: zcwtuibian@gmail.com 
    > Created Time: 2014年11月24日 星期一 06时57分01秒
 ************************************************************************/
//#include <stdio.h>
#include <iostream>
#include <ctime>
#include <fstream>
//#include <vector>
//#include <numeric>
//#include <gsl/gsl_multifit.h>
//extern "C"{
#include "OsadaSimilarityRand.h"
//}
using namespace std;


void Randsampling(SD & sd){
	static Area *AR=NULL;
	static RArea *RAR =NULL;

	//cout << __FILE__ <<endl;
	//AR=SD::CalArea_sort();   
	//RAR=SD::calsampling(AR);  
	//SD::Row2Column(RAR);
	AR=sd.CalArea_sort();   
	RAR=sd.calsampling(AR);  
	sd.Row2Column(RAR);

	delete(AR);
	delete(RAR);
}

//add by zcw,当上面的calsample的参数为三个时
int main(int argc, char *argv[])
{
        
    //add by levy @2013/10/30 for cal the time
    double time_start = (double)clock()/CLOCKS_PER_SEC;
   
    //FILE *fp;

    //char filename[20]="M000001.off";
    //char filename1[20]="airplane4.off";
    //char filename1[20],*filename2;
    
    //程序名 文件名 采样基数（一般是面片数） 采样点除以面片数（time）
    //strcpy(filename, argv[1]);
    //strcpy(filename1, argv[1]);
	//filename2 = argv[2];
     
		//char filename[100];
		//string filename;
		int randvertnum = 1024;//设置欲取三维模型表面的随机点数
		int idx;
		int binnum = 1024;//the number of statistic bin.
		int randcurvenum = 64;//the number of extracted random points on the polyfit curve.

/*
		string filename1;
		string filename2 = ".off";
		string filename;
		string filenameidx;
		string filepath = "RangeScans_QuerySet_2015_Final/";
		string filepath_name_temp;
        const char * filepath_name;
        //char * filepath_name;
        const int modelnum = 180;

		const int constVar1=10;
		const int constVar2=100;
		const int constVar3=1000;
		for(int i=1; i<=modelnum; i++){
			stringstream temps;
			temps << i; //将数值型变成字符串型
			filenameidx = temps.str();
			//cout << filenameidx <<endl;

			if(i< constVar1)
				filename1="RQ0000";
			else if(i< constVar2)
				filename1="RQ000";
			else if(i< constVar3)
				filename1="RQ00";
			else filename1="RQ0";

			filename=filename1 + filenameidx + filename2;
            cout << filename << endl;
			filepath_name_temp = filepath + filename;
*/

    //ofstream outf("facenum.pcd", ios::out);
    char filepath_name[100];
	for(int i=0;i<=1814;i++){
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/airplane/airplane%d.off",i); // the relative path is ok.
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/alarm_clock/alarm_clock%d.off",i); // the relative path is ok.
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/ant/ant%d.off",i); // the relative path is ok.
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/bed/bed%d.off",i); // the relative path is ok.
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/bee/bee%d.off",i); // the relative path is ok.
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/bench/bench%d.off",i); // the relative path is ok.
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/bus/bus%d.off",i); // the relative path is ok.
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/chair/chair%d.off",i); // the relative path is ok.
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/computer_monitor/computer_monitor%d.off",i); // the relative path is ok.
		//sprintf(filename,"SHREC14LSGTB_EXAMPLE_TARGET_MODELS/cup/cup%d.off",i); // the relative path is ok.
		sprintf(filepath_name,"psb_DB/m%d.off",i); // the relative path is ok.
		cout << filepath_name <<endl;
		//rand1 = (int)atoi(argv[2])*atof(argv[3]);
		//rand2 = (int)atoi(argv[2])*atof(argv[3]);

		static Mesh *mesh=NULL; 
		idx = i;	

        //filepath_name = filepath_name_temp.c_str();
		SD sd(filepath_name, randvertnum, idx, binnum, randcurvenum);

		//int n=5;//要拟合成的曲线的阶数

	    sd.ReadOffFile();  //读取了.off文件的所有信息—点坐标和三角面片顶点序号

        //int nface = sd.getFaceNum();
        //outf << nface << endl;
        //cout << nface << endl;

		//sd.off2pcd();  //将off文件转化为pcd文件（每个面片的顶点）

		mesh = Normalize_nopose(sd);
        
        //sd.scaledvertex2pcd(mesh); // convert scaled vertice to pcd file.

    	static Area *AR=NULL;
	    static RArea *RAR =NULL;

	    AR=sd.CalArea_sort();   // cal area of tripatch and sort them.
        sd.resetRandPointNum();

	    //RAR=sd.calPseudoSampling(AR);  
        //sd.pseudoRandpoints2pcd(RAR);

	    RAR=sd.calsampling(AR);  
	    sd.Row2Column(RAR);
		////Randsampling(sd);
        sd.Randpoints2pcd();
        
		//sd.StatDist();
		//sd.GetFeature();//存入到文件randy.txt中

        delete(mesh);
    }
    //outf.close();

/*
	static Mesh *mesh1=NULL; 
	static Area *AR1=NULL;
	static Vertex *centergravity1=NULL;
	static RArea *RAR1 =NULL;

	    sd.ReadOffFile();  //读取了.off文件的所有信息—点坐标和三角面片顶点序号
		AR1=sd.CalArea();   //得到三维模型的表面三角面片的面积AR1
		AR1=sd.ARGravity(AR1);  //计算每个面片的重心
		centergravity1=sd.ModelGravity(AR1);  //计算整个三维模型的重心
		mesh1=sd.GravityNormalize(centergravity1); //对整个三维模型顶点做标准化（平移到原点）
		mesh1=sd.ScaleNormalize();  //对整个三维模型的顶点做标准化（缩放）
		AR1=sd.CalArea_sort();   //得到三维模型的表面三角面片的面积AR1
		RAR1=sd.calsampling(AR1);  //在表面AR1的基础上得到Rand2个随机采样点
		sd.Row2Column(RAR1);//将mesh的顶点信息存入到一个vector的二维动态数组中

		sd.StatDist();
		//cout << __FILE__ <<endl;

		sd.GetFeature();
		 
		delete(mesh1);
		delete(centergravity1);
		//delete(AR);
		delete(AR1);
		//delete(RAR);
		delete(RAR1); 
*/
		double time_end = (double)clock()/CLOCKS_PER_SEC;
		double gap = time_end - time_start;
		//printf("time:%fs\n", gap);
		cout << gap <<endl;
		
		return 0;		
}
