//============================================================================
// Name        : AlleleSurfing3_13.cpp
// Author      : Ailene
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

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
using namespace std;
char FileName[50];
ofstream out_Pars;
ofstream out_Dem;
ofstream out_Coal;
//classes
class parsSet;
class ind
{
	public:
		bool**hLoci,**oLoci;
		int***coal;
		ind();
		~ind();
		double calcW(parsSet pars);
		void printInd(parsSet pars);
		void initalizeInd(parsSet pars);

	private:
};
class parsSet
{
	public:
		int scene,nHybrid,nOther,nPop,*nInd1,*nInd2,K,coalLocus,nIndTot;
		double R,m,*r,mu,*Wavg;
		double*ES,**freqInit,**sMtrx,**AFs,*Savg;
		ind**pop1,**pop2;
		parsSet();
		~parsSet();
		void initalizePops();
		void inializeScene();
		void calcAFs(int pop);
		void switchPop();
		void migration(); //Deterministic migration
		void migration2();
		void migrationRand(); //Stochastic migration
		void mating(); //Mating with deterministic population sizes.
		void mating2();
		void matingRand();  //Mating with stochastic # of offspring
		void calcMeanW();
		void calcMeanS();
		void printCoal();
	private:
};
//Function declarations
void inPut(parsSet*parsPtr);
int poisson(double lambda);
void sort(int*array,int length);
double roundn(double in);
int g,gmax;
int main()
{
	srand((unsigned int)time(NULL));//Initializing the random seed.
	//Parameters
	parsSet pars;
	char name[50];
	inPut(&pars);
	out_Pars<<"Average number of pairwise differences"<<endl;
	for(int rep=0;rep<=0;rep++)
	{
		cout<<"Rep #: "<<rep<<endl;
		//Specifying output file for repatation.
		sprintf(name,"//Users//ailenemacpherson//Documents//workspace//AlleleSurfing3_13//data//coal%d_%d.csv",pars.scene,rep);
		out_Coal.open(name);
		//Initializing
		pars.initalizePops();
		for(int p=0;p<pars.nPop;p++) out_Dem<<pars.nInd1[p]<<",";
		for(g=0;g<gmax; g++)
		{
			if(g%50==0) cout<<"g: "<<g<<endl;
			out_Dem<<endl;
			pars.migration();
			pars.mating();
			for(int p=0;p<pars.nPop;p++) out_Dem<<pars.nInd1[p]<<",";
		}
		out_Dem<<endl;
		for(int p=0;p<pars.nPop;p++) out_Dem<<pars.nInd1[p]<<",";
		pars.printCoal();
		//number of pairwise differences.
		pars.calcMeanS();
		for(int p=0;p<pars.nPop;p++) out_Pars<<pars.Savg[p]<<",";
		out_Pars<<endl;
		//closing files
		out_Dem.close();
		out_Coal.close();
	}
	out_Pars.close();
	cout<<"DONE!"<<endl;
	return 0;
}
//member functions
parsSet::parsSet()
{
	nInd1=NULL;
	pop1=NULL;
	pop2=NULL;
	AFs=NULL;
	sMtrx=NULL;
	ES=NULL;
	r=NULL;
	freqInit=NULL;
	Savg=NULL;
	R=0;K=0;nPop=0;mu=0;coalLocus=0;nHybrid=0;nOther=0;nInd1=0;nInd2=0;Wavg=0;scene=0;m=0,nIndTot=0;
}
parsSet::~parsSet()
{
}
void parsSet::initalizePops()
{
	//Initializing hybrid loci
	nHybrid=2; //There are always 2 hybrid loci although only one of them may be used in a particular scenario
	sMtrx=new double*[3];
	for(int h1=0;h1<3;h1++)
	{
		sMtrx[h1]=new double[3];
	}
	AFs=new double*[2];
	AFs[0]=new double[nHybrid];
	AFs[1]=new double[nOther];
	Wavg=new double[nPop];
	Savg=new double[nPop];
	//Initializing populations and setting initial allele frequencies
	pop1=new ind*[nPop];
	pop2=new ind*[nPop];
	for(int p=0;p<nPop;p++)
	{
		pop1[p]=new ind[2*K];
		pop2[p]=new ind[2*K];
		for(int i=0;i<2*K;i++)
		{
			pop1[p][i].hLoci=new bool*[2];
			pop1[p][i].oLoci=new bool*[2];
			pop2[p][i].hLoci=new bool*[2];
			pop2[p][i].oLoci=new bool*[2];
			for(int c=0;c<2;c++)
			{
				pop1[p][i].hLoci[c]=new bool[nHybrid];
				pop1[p][i].oLoci[c]=new bool[nOther];
				pop2[p][i].hLoci[c]=new bool[nHybrid];
				pop2[p][i].oLoci[c]=new bool[nOther];
			}
			pop1[p][i].coal=new int**[gmax+1];
			pop2[p][i].coal=new int**[gmax+1];
			for(int gen=0;gen<gmax+1;gen++)
			{
				pop1[p][i].coal[gen]=new int*[2];
				pop2[p][i].coal[gen]=new int*[2];
				for(int c=0;c<2;c++)
				{
					pop1[p][i].coal[gen][c]=new int[3];
					pop2[p][i].coal[gen][c]=new int[3];
					for(int s=0;s<3;s++)
					{
						pop1[p][i].coal[gen][c][s]=0;
						pop2[p][i].coal[gen][c][s]=0;
					}
				}
			}
		}
	}
	inializeScene();
}
void parsSet::inializeScene()
{
	if(scene==0)//Single population with logistic growth
	{
		nPop=1;m=0;
		cout<<"Scenario 0: Single population under logistic growth- no load"<<endl;
		nInd1[0]=floor((double)K/50.0);
//		nInd1[0]=(double)K;
		//All wild type
		for(int p=0;p<nPop;p++)
		{
			//out_Pars<<p<<",";
			for(int i=0;i<nInd1[p];i++)
			{
				for(int c=0;c<2;c++)
				{
					for(int l=0;l<nHybrid;l++)	pop1[p][i].hLoci[c][l]=0;
					for(int l=0;l<nOther;l++)	pop1[p][i].oLoci[c][l]=0;
					pop1[p][i].coal[0][c][0]=p;pop1[p][i].coal[0][c][1]=i;pop1[p][i].coal[0][c][2]=c;
				}
			}
			calcAFs(p);
		}
	}
	else if(scene==1)//expansion only: no load and no hybrid effects.
	{
		cout<<"Scenario 1: Expansion- no load"<<endl;
		int nMove=roundn(m*K/2.0);
		if(nMove<1)
		{
			cout<<"nMove<1-> no population expansion will occur."<<endl<<"Is this okay? (Hit enter if yes!)"<<endl;
			getchar();
		}
		//initial # of individuals
		out_Pars<<"Initial population size(s): "<<endl;
		for(int p=0;p<nPop;p++)
		{
			if(p==0)
			{
				nInd1[p]=K; //Only population 0 begins with individuals in it.
			}
			else nInd1[p]=0;
		}
		out_Pars<<endl;
		//Initial allele frequency: all wild-type
		out_Pars<<"Initial Allele Frequencies:"<<endl<<"Pop,locus"<<endl;
		for(int p=0;p<nPop;p++)
		{
			out_Pars<<p<<",";
			for(int i=0;i<nInd1[p];i++)
			{
				for(int c=0;c<2;c++)
				{
					for(int l=0;l<nHybrid;l++)	pop1[p][i].hLoci[c][l]=0;
					for(int l=0;l<nOther;l++)	pop1[p][i].oLoci[c][l]=0;
					pop1[p][i].coal[0][c][0]=p;pop1[p][i].coal[0][c][1]=i;pop1[p][i].coal[0][c][2]=c;
				}
			}
			calcAFs(p);
		}

		//Selection: no fitness effects
		for(int l=0;l<nOther;l++) ES[l]=0.0;
//		out_Pars<<"Effect Sizes"<<endl;
//		for(int l=0;l<nOther;l++) out_Pars<<ES[l]<<",";
//		out_Pars<<endl;
		for(int h1=0;h1<3;h1++)
		{
			for(int h2=0;h2<3;h2++)
			{
				sMtrx[h1][h2]=0.0;
			}
		}
	}
	else if(scene==2)//expansion with load: load and no hybrid effects.
	{
		//initial # of individuals
		for(int p=0;p<nPop;p++)
		{
			if(p<5) nInd1[p]=K;
			else nInd1[p]=0;
		}
		//Initial allele frequency: all wild-type
		for(int p=0;p<nPop;p++)
		{
			for(int i=0;i<nInd1[p];i++)
			{
				for(int c=0;c<2;c++)
				{
					for(int l=0;l<nHybrid;l++)	pop1[p][i].hLoci[c][l]=0;
					for(int l=0;l<nOther;l++)	pop1[p][i].oLoci[c][l]=0;
				}
			}
		}
		//Selection: no fitness effects
		for(int l=0;l<nOther;l++) ES[l]=0.1;
//		out_Pars<<"Effect Sizes"<<endl;
//		for(int l=0;l<nOther;l++) out_Pars<<ES[l]<<",";
//		out_Pars<<endl;
		for(int h1=0;h1<3;h1++)
		{
			for(int h2=0;h2<3;h2++)
			{
				sMtrx[h1][h2]=rand()/(double)RAND_MAX*0.1;
			}
		}
	}
	else
	{
		cout<<"Error: enter different scenario"<<endl;
		getchar();
	}
}
void parsSet::calcAFs(int pop)
{
	//Calculates allele frequency of pop1[pop]
	for(int l=0;l<nHybrid;l++) AFs[0][l]=0.0;
	for(int l=0;l<nOther;l++) AFs[1][l]=0.0;
	for(int i=0;i<nInd1[pop];i++)
	{
		for(int l=0;l<nHybrid;l++) AFs[0][l]+=pop1[pop][i].hLoci[0][l]+pop1[pop][i].hLoci[1][l];
		for(int l=0;l<nOther;l++) AFs[1][l]+=pop1[pop][i].oLoci[0][l]+pop1[pop][i].oLoci[1][l];
	}
	if(nInd1[pop]>0)
	{
		for(int l=0;l<nHybrid;l++) AFs[0][l]/=2.0*nInd1[pop];
		for(int l=0;l<nOther;l++) AFs[1][l]/=2.0*nInd1[pop];
	}
//	for(int l=0;l<nHybrid;l++) out_Dem<<AFs[0][l]<<",";
//	out_Dem<<",";
//	for(int l=0;l<nOther;l++) out_Dem<<AFs[1][l]<<",";
//	out_Dem<<endl;

}
void parsSet::migrationRand()
{
	//Individuals move among populations, are copied from pop1 to pop2
	double random;
	for(int p=0;p<nPop;p++) nInd2[p]=0;
	for(int p=0;p<nPop;p++)
	{
		if(p==0)//left most population
		{
			for(int i=0;i<nInd1[p];i++)
			{
				if(rand()/(double)RAND_MAX<m/2.0)//Moves right
				{
					//cout<<"Individual moves 0->1"<<endl; getchar();
					for(int c=0;c<2;c++)
					{
						for(int l=0;l<nHybrid;l++) pop2[p+1][nInd2[p+1]].hLoci[c][l]=pop1[p][i].hLoci[c][l];
						for(int l=0;l<nOther;l++) pop2[p+1][nInd2[p+1]].oLoci[c][l]=pop1[p][i].oLoci[c][l];
						for(int gen=0;gen<=g;gen++)
						{
							pop2[p+1][nInd2[p+1]].coal[gen][c][0]=pop1[p][i].coal[gen][c][0];
							pop2[p+1][nInd2[p+1]].coal[gen][c][1]=pop1[p][i].coal[gen][c][1];
							pop2[p+1][nInd2[p+1]].coal[gen][c][2]=pop1[p][i].coal[gen][c][2];
						}
					}
					nInd2[p+1]++;
				}
				else//Stays
				{
					//cout<<"stays"<<endl; getchar();
					for(int c=0;c<2;c++)
					{
						for(int l=0;l<nHybrid;l++) pop2[p][nInd2[p]].hLoci[c][l]=pop1[p][i].hLoci[c][l];
						for(int l=0;l<nOther;l++) pop2[p][nInd2[p]].oLoci[c][l]=pop1[p][i].oLoci[c][l];
						for(int gen=0;gen<=g;gen++)
						{
							pop2[p][nInd2[p]].coal[gen][c][0]=pop1[p][i].coal[gen][c][0];
							pop2[p][nInd2[p]].coal[gen][c][1]=pop1[p][i].coal[gen][c][1];
							pop2[p][nInd2[p]].coal[gen][c][2]=pop1[p][i].coal[gen][c][2];
						}
					}
					nInd2[p]++;
				}
			}
		}
		else if(p==nPop-1) //right most population
		{
			for(int i=0;i<nInd1[p];i++)
			{
				if(rand()/(double)RAND_MAX<m/2.0)//Moves left
				{
					//cout<<"Individual moves  2->1"<<endl;
					for(int c=0;c<2;c++)
					{
						for(int l=0;l<nHybrid;l++) pop2[p-1][nInd2[p-1]].hLoci[c][l]=pop1[p][i].hLoci[c][l];
						for(int l=0;l<nOther;l++) pop2[p-1][nInd2[p-1]].oLoci[c][l]=pop1[p][i].oLoci[c][l];
						for(int gen=0;gen<=g;gen++)
						{
							pop2[p-1][nInd2[p-1]].coal[gen][c][0]=pop1[p][i].coal[gen][c][0];
							pop2[p-1][nInd2[p-1]].coal[gen][c][1]=pop1[p][i].coal[gen][c][1];
							pop2[p-1][nInd2[p-1]].coal[gen][c][2]=pop1[p][i].coal[gen][c][2];
						}
					}
					nInd2[p-1]++;
				}
				else//Stays
				{
					for(int c=0;c<2;c++)
					{
						for(int l=0;l<nHybrid;l++) pop2[p][nInd2[p]].hLoci[c][l]=pop1[p][i].hLoci[c][l];
						for(int l=0;l<nOther;l++) pop2[p][nInd2[p]].oLoci[c][l]=pop1[p][i].oLoci[c][l];
						for(int gen=0;gen<=g;gen++)
						{
							pop2[p][nInd2[p]].coal[gen][c][0]=pop1[p][i].coal[gen][c][0];
							pop2[p][nInd2[p]].coal[gen][c][1]=pop1[p][i].coal[gen][c][1];
							pop2[p][nInd2[p]].coal[gen][c][2]=pop1[p][i].coal[gen][c][2];
						}
					}
					nInd2[p]++;
				}
			}
		}
		else //migration both directions
		{
			for(int i=0;i<nInd1[p];i++)
			{
				random=rand()/(double)RAND_MAX;
				if(random<m/2.0)//Moves left
				{
					//cout<<"Individual moves  1->0"<<endl;
					for(int c=0;c<2;c++)
					{
						for(int l=0;l<nHybrid;l++) pop2[p-1][nInd2[p-1]].hLoci[c][l]=pop1[p][i].hLoci[c][l];
						for(int l=0;l<nOther;l++) pop2[p-1][nInd2[p-1]].oLoci[c][l]=pop1[p][i].oLoci[c][l];
						for(int gen=0;gen<=g;gen++)
						{
							pop2[p-1][nInd2[p-1]].coal[gen][c][0]=pop1[p][i].coal[gen][c][0];
							pop2[p-1][nInd2[p-1]].coal[gen][c][1]=pop1[p][i].coal[gen][c][1];
							pop2[p-1][nInd2[p-1]].coal[gen][c][2]=pop1[p][i].coal[gen][c][2];
						}
					}
					nInd2[p-1]++;
				}
				else if(random<m) //Moves right
				{
					//cout<<"Individual moves  1->2"<<endl;
					for(int c=0;c<2;c++)
					{
						for(int l=0;l<nHybrid;l++) pop2[p+1][nInd2[p+1]].hLoci[c][l]=pop1[p][i].hLoci[c][l];
						for(int l=0;l<nOther;l++) pop2[p+1][nInd2[p+1]].oLoci[c][l]=pop1[p][i].oLoci[c][l];
						for(int gen=0;gen<=g;gen++)
						{
							pop2[p+1][nInd2[p+1]].coal[gen][c][0]=pop1[p][i].coal[gen][c][0];
							pop2[p+1][nInd2[p+1]].coal[gen][c][1]=pop1[p][i].coal[gen][c][1];
							pop2[p+1][nInd2[p+1]].coal[gen][c][2]=pop1[p][i].coal[gen][c][2];
						}
					}
					nInd2[p+1]++;
				}
				else//Stays
				{
					for(int c=0;c<2;c++)
					{
						for(int l=0;l<nHybrid;l++) pop2[p][nInd2[p]].hLoci[c][l]=pop1[p][i].hLoci[c][l];
						for(int l=0;l<nOther;l++) pop2[p][nInd2[p]].oLoci[c][l]=pop1[p][i].oLoci[c][l];
						for(int gen=0;gen<=g;gen++)
						{
							pop2[p][nInd2[p]].coal[gen][c][0]=pop1[p][i].coal[gen][c][0];
							pop2[p][nInd2[p]].coal[gen][c][1]=pop1[p][i].coal[gen][c][1];
							pop2[p][nInd2[p]].coal[gen][c][2]=pop1[p][i].coal[gen][c][2];
						}
					}
					nInd2[p]++;
				}
			}
		}
	}
	/*
	out_Dem<<"After migration: ,";
	for(int p=0;p<nPop;p++) out_Dem<<nInd2[p]<<",";
	out_Dem<<endl;
	*/
}
void parsSet::migration()
{
	int moves,ct=0,nMove=roundn(m*K/2.0);
	int* move=new int[2*nMove],popTo;
	//This is a deterministic (from a population size standpoint) version of migration.  Exactly (m/2) N[p,t] individuals move left and (m/2)N[p,t] move right.
	//Important if m/2 K<1 no population expansion will occur.
	//Function moves individuals from population 1 to population 2.
	for(int p=0;p<nPop;p++) nInd2[p]=0; //Resetting individual counts
	for(int p=0;p<nPop;p++)
	{
		nMove=roundn(m*nInd1[p]/2.0);
		ct=0;
		for(int i=0;i<2*nMove;i++) move[i]=rand()%nInd1[p]; //Drawing random individuals to move, there is a very small probability I will draw the same individual twice but I won't worry about this.
		for(int i=0;i<nInd1[p];i++)
		{
			moves=0;
			for(int j=0;j<2*nMove;j++)
			{
				if(move[j]==i) moves=1;
			}
			if(moves==1 && ct<nMove)
			{
				popTo=max(p-1,0); //moves left
				ct++;
			}
			else if(moves==1)
			{
				popTo=min(p+1,nPop-1);//moves right
				ct++;
			}
			else popTo=p;
			for(int c=0;c<2;c++)
			{
				for(int l=0;l<nHybrid;l++) pop2[popTo][nInd2[popTo]].hLoci[c][l]=pop1[p][i].hLoci[c][l];
				for(int l=0;l<nOther;l++) pop2[popTo][nInd2[popTo]].oLoci[c][l]=pop1[p][i].oLoci[c][l];
				for(int gen=0;gen<=g;gen++)
				{
					pop2[popTo][nInd2[popTo]].coal[gen][c][0]=pop1[p][i].coal[gen][c][0];
					pop2[popTo][nInd2[popTo]].coal[gen][c][1]=pop1[p][i].coal[gen][c][1];
					pop2[popTo][nInd2[popTo]].coal[gen][c][2]=pop1[p][i].coal[gen][c][2];
				}
			}
			nInd2[popTo]++;
		}
//		cout<<"Population: "<<p<<" nInd2 "<<nInd2[p]<<endl;
	}
	delete[] move;
}
void parsSet::switchPop()
{
	//Moves individuals from pop1 to pop 2 without migration
	for(int p=0;p<nPop;p++)
	{
		nInd2[p]=nInd1[p];
		for(int i=0;i<nInd1[p];i++)
		{
			for(int c=0;c<2;c++)
			{
				for(int l=0;l<nHybrid;l++) pop2[p][i].hLoci[c][l]=pop1[p][i].hLoci[c][l];
				for(int l=0;l<nOther;l++) pop2[p][i].oLoci[c][l]=pop1[p][i].oLoci[c][l];
			}
		}
	}
}
void parsSet::matingRand()
{
	double lambda;
	int parent2,numOff,*remove;
	bool chromo;
	ind wrkOff;
	wrkOff.initalizeInd((*this));
	//Individuals in pop2 reproduce (each parent has a poisson number with lambda=1+R(1-N/K) offspring), surviving offspring are kept in pop1.
	for(int p=0;p<nPop;p++) nInd1[p]=0; //Reseting nInd counter.
	for(int p=0;p<nPop;p++)
	{
		lambda=1.0+R*(1.0-(double)nInd2[p]/K);
		for(int i=0;i<nInd2[p];i++)
		{
			numOff=poisson(lambda);
			for(int off=0;off<numOff;off++)
			{
				//generating parent 1 gamette
				chromo=rand()%2;
				for(int l=0;l<nHybrid;l++)
				{
					if(rand()/(double)RAND_MAX<r[0])
					{
						chromo=(chromo+1)%2;
					}
					wrkOff.hLoci[0][l]=pop2[p][i].hLoci[chromo][l];
				}
				if(rand()/(double)RAND_MAX<r[1])
				{
					chromo=(chromo+1)%2;
				}
				for(int l=0;l<nOther;l++)
				{
					if(rand()/(double)RAND_MAX<mu) wrkOff.oLoci[0][l]=(pop2[p][i].oLoci[chromo][l]+1)%2; //Mutation occurs
					else wrkOff.oLoci[0][l]=pop2[p][i].oLoci[chromo][l]; //No mutation
					if(l==coalLocus)
					{
						for(int gen=0;gen<=g;gen++)//coping history
						{
							wrkOff.coal[gen][0][0]=pop2[p][i].coal[gen][chromo][0];
							wrkOff.coal[gen][0][1]=pop2[p][i].coal[gen][chromo][1];
							wrkOff.coal[gen][0][2]=pop2[p][i].coal[gen][chromo][2];
						}
						//adding new entry
						wrkOff.coal[g+1][0][0]=p;wrkOff.coal[g+1][0][1]=i;wrkOff.coal[g+1][0][2]=chromo;
					}
					if(rand()/(double)RAND_MAX<r[2])
					{
						chromo=(chromo+1)%2;//Recombination between other loci
					}
				}
				//generating parent 2 gamete
				parent2=rand()%nInd2[p];//picking the other parent
				chromo=rand()%2;
				for(int l=0;l<nHybrid;l++)
				{
					if(rand()/(double)RAND_MAX<r[0]) chromo=(chromo+1)%2;//Recombination between hybrid loci
					wrkOff.hLoci[1][l]=pop2[p][parent2].hLoci[chromo][l];
				}
				if(rand()/(double)RAND_MAX<r[1]) chromo=(chromo+1)%2;
				for(int l=0;l<nOther;l++)
				{
					if(rand()/(double)RAND_MAX<mu) wrkOff.oLoci[1][l]=(pop2[p][parent2].oLoci[chromo][l]+1)%2; //Mutation occurs
					else wrkOff.oLoci[1][l]=pop2[p][parent2].oLoci[chromo][l]; //No mutation
					if(l==coalLocus)
					{
						for(int gen=0;gen<=g;gen++)//coping history
						{
							wrkOff.coal[gen][1][0]=pop2[p][parent2].coal[gen][chromo][0];
							wrkOff.coal[gen][1][1]=pop2[p][parent2].coal[gen][chromo][1];
							wrkOff.coal[gen][1][2]=pop2[p][parent2].coal[gen][chromo][2];
						}
						//adding new entry
						wrkOff.coal[g+1][1][0]=p;wrkOff.coal[g+1][1][1]=parent2;wrkOff.coal[g+1][1][2]=chromo;
					}
					if(rand()/(double)RAND_MAX<r[2]) chromo=(chromo+1)%2;//Recombination between other loci
				}
				if(rand()/(double)RAND_MAX<=wrkOff.calcW((*this)))//Survives
				{
					for(int c=0;c<2;c++)
					{
						for(int l=0;l<nHybrid;l++) pop1[p][nInd1[p]].hLoci[c][l]=wrkOff.hLoci[c][l];
						for(int l=0;l<nOther;l++) pop1[p][nInd1[p]].oLoci[c][l]=wrkOff.oLoci[c][l];
						for(int gen=0;gen<=g+1;gen++)
						{
							pop1[p][nInd1[p]].coal[gen][c][0]=wrkOff.coal[gen][c][0];
							pop1[p][nInd1[p]].coal[gen][c][1]=wrkOff.coal[gen][c][1];
							pop1[p][nInd1[p]].coal[gen][c][2]=wrkOff.coal[gen][c][2];
						}
					}
					nInd1[p]++;
				}
				else
				{
					//cout<<"DIES"<<endl;
					//getchar();
				}
			}
		}
		//Trimming population down to carrying capacity:
		int ct=0;
		if(nInd1[p]>K)
		{
//			cout<<"Above carrying capcity by: "<<nInd1[p]-K<<endl;
			remove=new int[nInd1[p]-K];
			for(int i=0;i<nInd1[p]-K;i++) remove[i]=rand()%(nInd1[p]);
//			for(int j=0;j<nInd1[p]-K;j++) cout<<remove[j]<<",";
//			cout<<endl;
			for(int i=0;i<nInd1[p];i++)
			{
				for(int j=0;j<nInd1[p]-K;j++)
				{
					if(i==remove[j])
					{
						ct++; //Remove individual i?
//						cout<<"Removeing individual: "<<i<<endl;
					}
				}
				if(ct>0)//If you must remove atleast one individual before individual i...
				{
					for(int c=0;c<2;c++)
					{
						for(int l=0;l<nHybrid;l++) pop1[p][max(i-ct,0)].hLoci[c][l]=pop1[p][i].hLoci[c][l];
						for(int l=0;l<nOther;l++) pop1[p][max(i-ct,0)].oLoci[c][l]=pop1[p][i].oLoci[c][l];
						for(int gen=0;gen<=g+1;gen++)
						{
							pop1[p][max(i-ct,0)].coal[gen][c][0]=pop1[p][i].coal[gen][c][0];
							pop1[p][max(i-ct,0)].coal[gen][c][1]=pop1[p][i].coal[gen][c][1];
							pop1[p][max(i-ct,0)].coal[gen][c][2]=pop1[p][i].coal[gen][c][2];
						}
					}
				}
			}
			delete[] remove;
			nInd1[p]=K;
		}
	}
}
void parsSet::mating()
{
	int par1,par2;
	bool chromo;
	ind wrkOff;
	wrkOff.initalizeInd((*this));
	for(int p=0;p<nPop;p++)
	{
		nInd1[p]=0;
//		cout<<"Population: "<<p<<" nOff: "<<roundn((double)nInd2[p]+R*(double)nInd2[p]*(1.0-(double)nInd2[p]/(double)K))<<endl;
//		getchar();
		for(int i1=0;i1<roundn((double)nInd2[p]+R*(double)nInd2[p]*(1.0-(double)nInd2[p]/(double)K));i1++)
		{
			//Draw parents at random
			par1=rand()%nInd2[p];par2=rand()%nInd2[p];
			//Produce offspring
			//generating parent 1 gamete
			chromo=rand()%2;
			for(int l=0;l<nHybrid;l++)
			{
				if(rand()/(double)RAND_MAX<r[0])
				{
					chromo=(chromo+1)%2;
				}
				wrkOff.hLoci[0][l]=pop2[p][par1].hLoci[chromo][l];
			}
			if(rand()/(double)RAND_MAX<r[1])
			{
				chromo=(chromo+1)%2;
			}
			for(int l=0;l<nOther;l++)
			{
				if(rand()/(double)RAND_MAX<mu) wrkOff.oLoci[0][l]=(pop2[p][par1].oLoci[chromo][l]+1)%2; //Mutation occurs
				else wrkOff.oLoci[0][l]=pop2[p][par1].oLoci[chromo][l]; //No mutation
				if(l==coalLocus)
				{
					for(int gen=0;gen<=g;gen++)//coping history
					{
						wrkOff.coal[gen][0][0]=pop2[p][par1].coal[gen][chromo][0];
						wrkOff.coal[gen][0][1]=pop2[p][par1].coal[gen][chromo][1];
						wrkOff.coal[gen][0][2]=pop2[p][par1].coal[gen][chromo][2];
					}
					//adding new entry
					wrkOff.coal[g+1][0][0]=p;wrkOff.coal[g+1][0][1]=par1;wrkOff.coal[g+1][0][2]=chromo;
				}
				if(rand()/(double)RAND_MAX<r[2])
				{
					chromo=(chromo+1)%2;//Recombination between other loci
				}
			}
			//generating parent 2 gamete
			chromo=rand()%2;
			for(int l=0;l<nHybrid;l++)
			{
				if(rand()/(double)RAND_MAX<r[0]) chromo=(chromo+1)%2;//Recombination between hybrid loci
				wrkOff.hLoci[1][l]=pop2[p][par2].hLoci[chromo][l];
			}
			if(rand()/(double)RAND_MAX<r[1]) chromo=(chromo+1)%2;
			for(int l=0;l<nOther;l++)
			{
				if(rand()/(double)RAND_MAX<mu) wrkOff.oLoci[1][l]=(pop2[p][par2].oLoci[chromo][l]+1)%2; //Mutation occurs
				else wrkOff.oLoci[1][l]=pop2[p][par2].oLoci[chromo][l]; //No mutation
				if(l==coalLocus)
				{
					for(int gen=0;gen<=g;gen++)//coping history
					{
						wrkOff.coal[gen][1][0]=pop2[p][par2].coal[gen][chromo][0];
						wrkOff.coal[gen][1][1]=pop2[p][par2].coal[gen][chromo][1];
						wrkOff.coal[gen][1][2]=pop2[p][par2].coal[gen][chromo][2];
					}
					//adding new entry
					wrkOff.coal[g+1][1][0]=p;wrkOff.coal[g+1][1][1]=par2;wrkOff.coal[g+1][1][2]=chromo;
				}
				if(rand()/(double)RAND_MAX<r[2]) chromo=(chromo+1)%2;//Recombination between other loci
			}
			if(rand()/(double)RAND_MAX<=wrkOff.calcW((*this)))//Survives
			{
				for(int c=0;c<2;c++)
				{
					for(int l=0;l<nHybrid;l++) pop1[p][nInd1[p]].hLoci[c][l]=wrkOff.hLoci[c][l];
					for(int l=0;l<nOther;l++) pop1[p][nInd1[p]].oLoci[c][l]=wrkOff.oLoci[c][l];
					for(int gen=0;gen<=g+1;gen++)
					{
						pop1[p][nInd1[p]].coal[gen][c][0]=wrkOff.coal[gen][c][0];
						pop1[p][nInd1[p]].coal[gen][c][1]=wrkOff.coal[gen][c][1];
						pop1[p][nInd1[p]].coal[gen][c][2]=wrkOff.coal[gen][c][2];
					}
				}
				nInd1[p]++;
			}
			else//Dies
			{
				//cout<<"Offspring dies"<<endl;
				i1--;
			}
		}
	}
}
void parsSet::calcMeanW()
{
	for(int p=0;p<nPop;p++)
	{
		Wavg[p]=0.0;
		for(int i=0;i<nInd1[p];i++)
		{
			Wavg[p]+=pop1[p][i].calcW((*this));
		}
		if(nInd1[p]>0) Wavg[p]/=(double)nInd1[p];
	}
}
void parsSet::printCoal()
{
	for(int gen=0;gen<=gmax;gen++)
	{
		for(int p=0;p<nPop;p++)
		{
			for(int i=0;i<nInd1[p];i++)
			{
				for(int c=0;c<2;c++)
				{
					out_Coal<<pop1[p][i].coal[gen][c][0]*(K*2)+pop1[p][i].coal[gen][c][1]*(2)+pop1[p][i].coal[gen][c][2]<<",";
				}
			}
		}
		out_Coal<<endl;
	}
	//Printing current state:
	for(int p=0;p<nPop;p++)
	{
		//cout<<"# of individuals: "<<nInd1[p]<<endl; getchar();
		for(int i=0;i<nInd1[p];i++)
		{
			for(int c=0;c<2;c++)
			{
				out_Coal<<p*(K*2)+i*(2)+c<<",";
			}
		}
	}
}
void parsSet::calcMeanS()
{
	double S,npairs=0;
	//function calculates the mean # of pairwise differences between all pairs of individuals in the population
	for(int p=0;p<nPop;p++)
	{
		Savg[p]=0;npairs=0;
		for(int i1=0;i1<nInd1[p];i1++)
		{
			for(int i2=i1+1;i2<nInd1[p];i2++)
			{
				/*for(int c=0;c<2;c++)
				{
					for(int l=0;l<nOther;l++) cout<<pop1[p][i1].oLoci[c][l]<<",";
					cout<<endl;
				}
				for(int c=0;c<2;c++)
				{
					for(int l=0;l<nOther;l++) cout<<pop1[p][i2].oLoci[c][l]<<",";
					cout<<endl;
				}*/
				S=0;npairs++;
				for(int c=0;c<2;c++)
				{
					for(int l=0;l<nOther;l++)
					{
						if(pop1[p][i1].oLoci[c][l]!=pop1[p][i2].oLoci[c][l]) S++;
					}
				}
				//cout<<"S is:"<<S<<endl;getchar();
				Savg[p]+=S;
			}
		}
		Savg[p]/=npairs;
	}
}
ind::ind()
{
	hLoci=NULL;
	oLoci=NULL;
	coal=NULL;
}
ind::~ind()
{
	for(int c=0;c<2;c++)
	{
		delete[] hLoci[c];
		delete[] oLoci[c];
	}
	delete[] hLoci;
	delete[] oLoci;
	for(int gen=0;gen<gmax;gen++)
	{
		for(int c=0;c<2;c++)
		{
			delete[] coal[gen][c];
		}
		delete[] coal[gen];
	}
	delete[] coal;
}
double ind::calcW(parsSet pars)
{
	double W=1.0;
	for(int c=0;c<2;c++)
	{
		for(int l=0;l<pars.nOther;l++) W*=pow(1.0-pars.ES[l],oLoci[c][l]);
	}
	W*=(1.0-pars.sMtrx[hLoci[0][0]+hLoci[1][0]][hLoci[0][1]+hLoci[1][1]]);
	//cout<<W<<endl;
	return W;
}
void ind::printInd(parsSet pars)
{
	for(int c=0;c<2;c++)
	{
		for(int l=0;l<pars.nHybrid;l++) cout<<hLoci[c][l]<<",";
		for(int l=0;l<pars.nOther;l++) cout<<oLoci[c][l]<<",";
		cout<<endl;
	}
}
void ind::initalizeInd(parsSet pars)
{
	hLoci=new bool*[2];
	oLoci=new bool*[2];
	for(int c=0;c<2;c++)
	{
		hLoci[c]=new bool[pars.nHybrid];
		oLoci[c]=new bool[pars.nOther];
	}
	coal=new int**[gmax+1];
	for(int gen=0;gen<gmax+1;gen++)
	{
		coal[gen]=new int*[2];
		for(int c=0;c<2;c++)
		{
			coal[gen][c]=new int[3];
			for(int s=0;s<3;s++)
			{
				coal[gen][c][s]=0;
			}
		}
	}
}
//functions
void inPut(parsSet*parsPtr)
{
	char name[50];
	ifstream inFile;
	inFile.open("inPut.txt");
	string myString;
	myString="Nothing";
	/*while(myString != "name:" && inFile.good())	{inFile>>myString;}
	inFile>>myString;
	cout<<"parameter file name: "<<myString<<endl;
	out_Pars.open(myString);
	while(myString != "name:" && inFile.good())	{inFile>>myString;}
	inFile>>myString;
	cout<<"data file name: "<<myString<<endl;
	out_Dem.open(myString);
	myString="Nothing";
	while(myString != "name:" && inFile.good())	{inFile>>myString;}
	inFile>>myString;
	cout<<"coalescent file name: "<<myString<<endl;
	out_Coal.open(myString);*/
	while(myString != "(scene):" && inFile.good())	{inFile>>myString;}
	inFile>>(*parsPtr).scene;
	sprintf(name,"pars%d.csv",(*parsPtr).scene);//Makes the parameter file name.
	out_Pars.open(name);
	sprintf(name,"data%d.csv",(*parsPtr).scene);//Makes the data file name.
	out_Dem.open(name);
	out_Pars<<"scenario:, "<<(*parsPtr).scene<<endl;
	while(myString != "(gmax):" && inFile.good())	{inFile>>myString;}
	inFile>>gmax;
	out_Pars<<"gmax:, "<<gmax<<endl;
	while(myString != "(nPop):" && inFile.good())	{inFile>>myString;}
	inFile>>(*parsPtr).nPop;
	out_Pars<<"# of populations:, "<<(*parsPtr).nPop<<endl;
	(*parsPtr).nInd1=new int[(*parsPtr).nPop];
	(*parsPtr).nInd2=new int[(*parsPtr).nPop];
	/*while(myString != "(nInd1):" && inFile.good())	{inFile>>myString;}
	out_Pars<<"inital # of individuals:, ";
	for(int p=0;p<(*parsPtr).nPop;p++)
	{
		inFile>>(*parsPtr).nInd1[p];
		out_Pars<<(*parsPtr).nInd1[p]<<",";
	}
	out_Pars<<endl;*/
	while(myString != "(nOther):" && inFile.good())	{inFile>>myString;}
	inFile>>(*parsPtr).nOther;
	out_Pars<<"# of other loci:, "<<(*parsPtr).nOther<<endl;
	for(int l=0;l<(*parsPtr).nOther;l++) (*parsPtr).ES=new double[(*parsPtr).nOther];
	while(myString != "(K):" && inFile.good())	{inFile>>myString;}
	inFile>>(*parsPtr).K;
	out_Pars<<"Carrying Capacity:, "<<(*parsPtr).K<<endl;
	while(myString != "(R):" && inFile.good())	{inFile>>myString;}
	inFile>>(*parsPtr).R;
	out_Pars<<"Max Growth rate:, "<<(*parsPtr).R<<endl;
	while(myString != "(m):" && inFile.good())	{inFile>>myString;}
	inFile>>(*parsPtr).m;
	out_Pars<<"migration rate:, "<<(*parsPtr).m<<endl;
	(*parsPtr).r=new double[3];
	while(myString != "(r):" && inFile.good())	{inFile>>myString;}
	out_Pars<<"recombination rate:, ";
	for(int c=0;c<3;c++)
	{
		inFile>>(*parsPtr).r[c];
		out_Pars<<(*parsPtr).r[c]<<",";
	}
	out_Pars<<endl;
	while(myString != "(mu):" && inFile.good())	{inFile>>myString;}
	inFile>>(*parsPtr).mu;
	out_Pars<<"mutation rate:, "<<(*parsPtr).mu<<endl;
	while(myString != "(coalLocus):" && inFile.good())	{inFile>>myString;}
	inFile>>(*parsPtr).coalLocus;
	out_Pars<<"focal coalescent lcous:, "<<(*parsPtr).coalLocus<<endl;
}
int poisson(double lambda)
{
	//Samples a random number from a poisson distribution with mean lambda.
	double L,p;
	int k;
	L=exp(-1.0*lambda);
	k=0;
	p=1;
	while(p>L)
	{
		k=k+1;
		p=p*rand()/(double)RAND_MAX;
	}
	return k-1;
}
void sort(int*array,int length)
{
	int temp;
	for(int i=0;i<length;i++)
	{
		for(int j=0;j<length-1;j++)
		{
			if(array[j]<array[i])//Swap them.
			{
				temp=array[i];
				array[i]=array[j];
				array[j]=temp;
			}
		}
	}
}
double roundn(double in)
{
	if((in-floor(in))<=0.5) return floor(in);
	else return ceil(in);
}
//Add a poisson sampling function
//Add/Finish mating function
