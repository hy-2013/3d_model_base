/*************************************************************************
    > File Name: Distance_vector.cpp
    > Author: clzhang 
    > Mail: zcwtuibian@gmail.com 
    > Created Time: 2014年11月28日 星期五 22时26分32秒
 ************************************************************************/

/**
Implement four similarity distance calculation,as
Euclidean, Manhaton, Chebyshev, Mahalanobis distance.
**/

/*#include<iostream>*/
#include<vector>
/*#include<string>*/
/*#include<cmath>*/
#include<cstring>
/*#include<cassert>*/
#include<algorithm>
#include<Eigen/Dense>
using namespace std;

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
double distVector(const vector<elemtype> v1, const vector<elemtype>  v2,  const char distType[]){
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

/*
int main(void)
{
	const int elemnum = 3;
	double a[elemnum]={11.0,2.3,3.4};
	double b[elemnum]={2.1,15.8,9.54};
	vector<double> dv1(a, a+elemnum);
	vector<double> dv2(b, b+elemnum);
*/	
	/*string distname("Oshi");//the string cannot assignment as (distname= "Oshi"), since "Oshi" is a char[] type.*/
	/*string distname1("Manhaton");*/
	/*string distname2("Chebyshev");*/
	/*string distname3("Mashi");*/
	/*string distname3("Mashi");*/
	/*kk=Oshi(a,b, 3, 3);*/
	/*kk=distVector(&dv1,&dv2,(string)"Oshi");*/
    /*int	kk=distVector(dv1,dv2,distname);*/
	/*cout << "kk= " << kk <<endl;*/
/*
	cout<<"Oshi distance= "<<distVector(dv1,dv2,"Oshi")<<endl; // the method to call distance function.
	cout<<"Manhaton distance= "<<distVector(dv1,dv2,"Manhaton")<<endl;
	cout<<"Chebyshev distance= "<<distVector(dv1,dv2,"Chebyshev")<<endl;

	const int elemnum1 = 2;
	double his1[elemnum1]={11, 17};
	double his2[elemnum1]={2, 10};
	vector<double> dv3(his1, his1+elemnum1);
	vector<double> dv4(his2, his2+elemnum1);
	
	cout<<"Mashi distance= "<<distVector(dv3,dv4,"Mashi")<<endl;
*/
	/*Eigen::MatrixXd m(2,2);*/
	    /*m(0,0) = 3;*/
		    /*m(1,0) = 2;*/
			    /*m(0,1) = 2;*/
				    /*m(1,1) = 1;*/
	/*cout << m << endl<<endl;*/
	/*cout << m.determinant() << endl;*/

	/*return 0;*/
/*}*/


