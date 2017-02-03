// AlleleSurfing1_31_17.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <new>
#include <memory>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <complex>
#include <algorithm>

using namespace std;
ofstream out_Pars;
ofstream out_Data;
//Classes
class ind
{
public:
	bool**other;
	bool**hybrid;
	double W;
	ind();
	~ind();
	void ind::calcW(parsSet pars);
private:
};
class parsSet
{
public:
	int numPop,K,numOther,numHybrid;
	double*ES;
	parsSet();
	~parsSet();
private:
};
class pop
{
public:
	int numInd;
	double Wavg,**AFs;
	ind*indVec;
	pop::pop();
	pop::~pop();
	void pop::calcAF(parsSet pars);
	void pop::printAF(parsSet pars);
	void pop::calcWavg(parsSet pars);
private:
	
};
//Function Declarations
void input(parsSet*parsPtr);
void initalizeScinerio(parsSet*parsPtr);
void initalizePops(parsSet,pop**popVec1Ptr,pop**popVec2Ptr);
void initalizeAF(pop* popVec1,parsSet pars);
int _tmain(int argc, _TCHAR* argv[])
{
	//Variables
	pop* popVec1,*popVec2;
	parsSet pars;
	input(&pars);
	initalizePops(pars,&popVec1,&popVec2);

	getchar();
	return 0;
}
//Member functions
	//Population
pop::pop()
{
	indVec=NULL;
	numInd=0;
	Wavg=0;
}
pop::~pop()
{
	delete[] indVec;
}
void pop::calcAF(parsSet pars)
{
	for(int l=0;l<pars.numHybrid;l++)
	{
		AFs[0][l]=0;
		for(int i=0;i<numInd;i++)
		{
			for(int c=0;c<2;c++)
			{
					AFs[0][l]+=indVec[i].hybrid[c][l];
			}
		}
		if(numInd>0) AFs[0][l]/=2.0*numInd;
	}
	for(int l=0;l<pars.numOther;l++)
	{
		AFs[1][l]=0;
		for(int i=0;i<numInd;i++)
		{
			for(int c=0;c<2;c++)
			{
					AFs[1][l]+=indVec[i].other[c][l];
			}
		}
		if(numInd>0) AFs[1][l]/=2.0*numInd;
	}
}
void pop::printAF(parsSet pars)
{
	for(int l=0;l<pars.numHybrid;l++)
	{
		cout<<AFs[0][l]<<",";
	}
	cout<<",";
	for(int l=0;l<pars.numOther;l++)
	{
		cout<<AFs[1][l]<<",";
	}
	cout<<",";	
}
void pop::calcWavg(parsSet pars)
{
	Wavg=0.0;
	for(int i=0;i<numInd;i++)
	{
		Wavg+=indVec[i].W;
	}
	Wavg/=(double)numInd;
}
	//Pars
parsSet::parsSet()
{
}
parsSet::~parsSet()
{
}
	//Individual
ind::ind()
{
	other=new bool*[2];
	hybrid=new bool*[2];
	W=0.0;
}
ind::~ind()
{
	delete[] other;
	delete[] hybrid;
}
void ind::calcW(parsSet pars)
{
	W=1.0;
	//Scenerio 1.1//
	for(int l=0;l<pars.numOther;l++)
	{
		for(int c=0;c<2;c++)
		{
			W+=other[c][l]*pars.ES[l];
		}
	}
}
//Functions
void input(parsSet*parsPtr)
{
	ifstream inFile;
	inFile.open("InFile.txt");
	string myString,filename;
	myString="Nothing";
	while(myString != "name:" && inFile.good()){inFile>>myString;}
	inFile>>filename;
	cout<<"pars name: "<<filename<<endl;
	out_Pars.open(filename);
	myString="Nothing";
	while(myString != "name:" && inFile.good()){inFile>>myString;}
	inFile>>filename;
	cout<<"data name: "<<filename<<endl;
	out_Data.open(filename);
	while(myString != "populations:" && inFile.good()) {inFile>>myString;}
	inFile>>(*parsPtr).numPop;
	cout<<"number of populations: "<<(*parsPtr).numPop<<endl;
	out_Pars<<"# pops,"<<(*parsPtr).numPop<<endl;
	while(myString != "(K):" && inFile.good()) {inFile>>myString;}
	inFile>>(*parsPtr).K;
	cout<<"K: "<<(*parsPtr).K<<endl;
	out_Pars<<"K,"<<(*parsPtr).K<<endl;
	while(myString != "loci:" && inFile.good()) {inFile>>myString;}
	inFile>>(*parsPtr).numHybrid;
	cout<<"numHybrid: "<<(*parsPtr).numHybrid<<endl;
	out_Pars<<"numHybrid,"<<(*parsPtr).numHybrid<<endl;
	myString="Nothing";
	while(myString != "loci:" && inFile.good()) {inFile>>myString;}
	inFile>>(*parsPtr).numOther;
	cout<<"numOther: "<<(*parsPtr).numOther<<endl;
	out_Pars<<"numOther,"<<(*parsPtr).numOther<<endl;
	inFile.close();

	initalizeScinerio(parsPtr);
}
void initalizeScinerio(parsSet*parsPtr)
{
	//Sum of fitness effects must be less than -1.0"
//Set fitness effects:
	(*parsPtr).ES=new double[(*parsPtr).numOther];
	/*Scenerio 1.1*/
	for(int l=0;l<(*parsPtr).numOther;l=0)
	{
		(*parsPtr).ES[l]=0.0;
	}

	out_Pars<<"Fitness Effects: "<<endl;
	for(int l=0;l<(*parsPtr).numOther;l=0)
	{
		out_Pars<<(*parsPtr).ES[l]<<",";
	}
	out_Pars<<endl;
}
void initalizePops(parsSet pars,pop**popVec1Ptr,pop**popVec2Ptr)
{
	string myString;
	ifstream inFile;
	inFile.open("InFile.txt");
	while(myString != "Conditions:" && inFile.good()){inFile>>myString;}
	(*popVec1Ptr)=new pop[pars.numPop];
	(*popVec2Ptr)=new pop[pars.numPop];
	for(int p=0;p<pars.numPop;p++)
	{
		(*popVec1Ptr)[p].indVec=new ind[2*pars.K];
		(*popVec2Ptr)[p].indVec=new ind[2*pars.K];
		for(int i=0;i<2*pars.K;i++)
		{
			for(int c=0;c<2;c++)
			{
				(*popVec1Ptr)[p].indVec[i].hybrid[c]=new bool[pars.numHybrid];
				(*popVec1Ptr)[p].indVec[i].other[c]=new bool[pars.numOther];
				(*popVec2Ptr)[p].indVec[i].hybrid[c]=new bool[pars.numHybrid];
				(*popVec2Ptr)[p].indVec[i].other[c]=new bool[pars.numOther];
			}
		}
		(*popVec1Ptr)[p].AFs=new double*[2];
		(*popVec2Ptr)[p].AFs=new double*[2];
		for(int b=0;b<2;b++)/*element 1 is for hybrid loci and element 2 is for other loci*/
		{
			(*popVec1Ptr)[p].AFs[b]=new double[pars.numHybrid];
			(*popVec2Ptr)[p].AFs[b]=new double[pars.numOther];
		}
	}
	
	//Initalizing number of individuals
	cout<<"Inital number of individuals:"<<endl;
	while(myString != "individuals:" && inFile.good()){inFile>>myString;}
	for(int p=0;p<pars.numPop;p++)
	{
		inFile>>(*popVec1Ptr)[p].numInd;
		cout<<(*popVec1Ptr)[p].numInd<<",";
	}
	cout<<endl;
	cout<<"Inital allele Frequecies:"<<endl;
	for(int p=0;p<pars.numPop;p++)
	{
		for(int i=0;i<(*popVec1Ptr)[p].numInd;i++)
		{
			for(int c=0;c<2;c++)
			{
				for(int l=0;l<pars.numHybrid;l++)
				{
					(*popVec1Ptr)[p].indVec[i].hybrid[c][l]=0;
				}
				for(int l=0;l<pars.numOther;l++)
				{
					(*popVec1Ptr)[p].indVec[i].other[c][l]=0;
				}
			}
		}
	}
	initalizeAF((*popVec1Ptr),pars);
	cout<<endl;
	/**/
	//Calculate allele frequencies
	//Calculate individual fitnesses
	//Calculate population mean fitness
	//Test 
}
void initalizeAF(pop* popVec1,parsSet pars)
{
	/*Scinerio 1: (All populations) Hybrid: Fixed for 0, Other: Fixed for 0*/
	double*AFVec=new double[pars.numOther];
	for(int p=0;p<pars.numPop;p++)
	{
		for(int i=0;i<popVec1[p].numInd;i++)
		{
			for(int c=0;c<2;c++)
			{
				//Hybrid locous
				for(int l=0;l<pars.numHybrid;l++) popVec1[p].indVec[i].hybrid[c][l]=0;
				//Other locous
				for(int l=0;l<pars.numOther;l++) popVec1[p].indVec[i].other[c][l]=0;
			}
		}
		popVec1[p].calcAF(pars);
		cout<<"\n Population "<<p<<endl;
		//popVec1[p].printAF(pars);
	}
}
void selection1(parsSet pars,pop*popVecIn,pop*popVecOut)
{
	int cnt;
	for(int p=0;p<pars.numPop;p++)
	{
		cnt=0;
		for(int i=0;i<popVecIn[p].numInd;i++)
		{
			if(rand()/(double)RAND_MAX<popVecIn[p].indVec[i].W)//Survives
			{
				for(int c=0;c<2;c++)
				{
					for(int l=0;l<pars.numHybrid;l++)
					{
						popVecOut[p].indVec[cnt].hybrid[c][l]=popVecIn[p].indVec[i].hybrid[c][l];
					}
					for(int l=0;l<pars.numOther;l++)
					{
						popVecOut[p].indVec[cnt].other[c][l]=popVecIn[p].indVec[i].other[c][l];
					}
				}
				cnt++;
			}
		}
	}
}


//Add selection
//Add migration
//Add mating
