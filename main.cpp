#include <iostream>
#include "stdio.h"
#include <fstream>
#include "time.h" 

using namespace std; 
#define INFTY 2147483647

int BestScore(char **sDNA, int l, int t, char *bestCha)   //calculate the best score for each leaf
{
	int sA, sC, sG, sT, Tbest=0;
	int best = 0;
	
	
	
	for(int i=0; i<l; i++)
		
	{
		
		
		sA=0; sC=0; sG=0; sT=0;
		for(int j=0;j<t;j++)
		{
			
			if(sDNA[j][i]=='a')
			{
				sA++;
			}else if(sDNA[j][i]=='c')
			{
				
				sC++;
			}else if(sDNA[j][i]=='g')
			{
			    sG++;
			}else if(sDNA[j][i]=='t')
			{
				
			    sT++;
			}
			
			
		}
		
		if(sA>best)
		{
			best=sA;
			bestCha[i]='a';
		}
		if(sC>best)
		{
			best=sC;
			bestCha[i]='c';
		}
		if(sG>best)
		{
			best=sG;
			bestCha[i]='g';
		}
		if(sT>best)
		{
			best=sT;
			bestCha[i]='t';
		}
		
		Tbest=Tbest+best;
		best=0;
	}
	
	
	return Tbest;
}



char* BruteForceMotifSearchAgain(char **DNA, int t, int n, int l)
{
	int *a = new int[n-l+1];
	char **sDNA=new char *[t];
	char *bestCha=new char[l];
	char *mutif = new char[l];
	int k=0,Tbest=0,bScore=0,endflag=0;
	for (int i=0; i<n-l+1 ; i++) 
	{
		a[i]=0;
	}
	
	for (int i=0; i<t; i++) 
	{
		sDNA[i]=new char [l];
	}
	
	
	while (1)
	{
		
		for(int i=0;i<t;i++)
		{
			k=a[i];
			for(int j=0; j<l; j++ )
			{
				sDNA[i][j]= DNA[i][k+j];
			}
		}
		for(int i=t-1;i>=0;i--)      //calculate the next leaf
		{
			if(a[i]==(n-l))
			{
				a[i]=0;
				if (i==0) {
					endflag=1;
				}
				
			}else
			{
				a[i]++;
				break;
			}
			
		}
		
	    
		Tbest= BestScore(sDNA,l, t, bestCha);
		
		if(Tbest>bScore)    //find the leaf with a best score and put it in the mutif 
		{
			bScore=Tbest;
			for(int i=0; i<l; i++ )
			{
				mutif[i]=bestCha[i];
			}
			
		}
		if(endflag==1)   //if go over all the leaves, break loop
			break;
		
	}
	
	return mutif;    //return the best mutif
	
}



int Score(char **sDNA, int l, int ti, char *bestCha)   //calculate the best score for each leaf
{
	int sA, sC, sG, sT, Tbest=0;
	int best = 0;
	
	
	
	
	for(int i=0; i<l; i++)
		
	{
		sA=0; sC=0; sG=0; sT=0;
		for(int j=0;j<=ti;j++)
		{
			
			if(sDNA[j][i]=='a')
			{
				sA++;
			}else if(sDNA[j][i]=='c')
			{
				
				sC++;
			}else if(sDNA[j][i]=='g')
			{
			    sG++;
			}else if(sDNA[j][i]=='t')
			{
				
			    sT++;
			}
			
			
		}
		
		if(sA>best)
		{
			best=sA;
			bestCha[i]='a';
		}
		if(sC>best)
		{
			best=sC;
			bestCha[i]='c';
		}
		if(sG>best)
		{
			best=sG;
			bestCha[i]='g';
		}
		if(sT>best)
		{
			best=sT;
			bestCha[i]='t';
		}
		
		Tbest=Tbest+best;
		best=0;
	}
	
	//cout <<bestCha[0]<<bestCha[1]<<bestCha[2]<<bestCha[3]<<endl;
	//cout << Tbest<<endl;
	return Tbest;
}



char* BranchAndBoundMotifSearch(char **DNA, int t, int n, int l)
{
	int *a = new int[n-l+1];
	char **sDNA=new char *[t];
	char *bestCha=new char[l];
	char *mutif = new char[l];
	int k=0,iscore=0,tscore=0,bScore=0,endflag=0,optimisticScore;
	for (int i=0; i<n-l+1 ; i++) 
	{
		a[i]=0;
	}
	
	for (int i=0; i<t ; i++) 
	{
		sDNA[i]=new char [l];
	}
	
	
	while (1)
	{   
		for(int i=0;i<t;i++)
		{
			k=a[i];
			for(int j=0; j<l; j++ )
			{
				sDNA[i][j]= DNA[i][k+j];
				
			}
			
			iscore= Score(sDNA,l, i, bestCha);   
			optimisticScore= iscore+(t-i-1)*l;
			if(i!=(t-1))
			{
				if (optimisticScore<bScore) 
				{      //bypasss
					for(int j=i; j>=0; j--)
					{
						if(a[j]==(n-l))
						{
							a[j]=0;
							if (j==0)
							{
								endflag=1;
							}
						}else
						{
							a[j]++;
							break;
						}
						
						
					}
					for(int j=i+1; j<t; j++)
					{
						
						a[j]=0;
					}
				}
			}
		}
		
		for(int i=(t-1);i>=0;i--)      //calculate the next leaf
		{
			if(a[i]==(n-l))
			{
				a[i]=0;
				if (i==0)
				{
					endflag=1;
				}
				
			}else
			{
				a[i]++;
				break;
			}
			
		}
		
		tscore= Score(sDNA,l, t-1, bestCha);   
		
		if(iscore>bScore)
		{
			bScore=tscore;
			for(int i=0; i<l; i++ )
			{
				mutif[i]=bestCha[i];
			}
			
		}
		
		
		if(endflag==1)   //if go over all the leaves, break loop
			break;
		
	}
	
	return mutif;    //return the best mutif
	
}


int TotalDistance(char *s, char **DNA, int t, int n, int l)
{
    char *v = new char[l];
	int Distance, MinDistance=INFTY, SumDistance=0;
	
	for (int i=0; i< t; i++)
	{
		MinDistance=INFTY;	
		for (int j=0; j<(n-l+1); j++)
		{
			Distance=0;
			for(int k=0; k<l; k++)
			{
				v[k]=DNA[i][j+k];
				
				if(v[k]!=s[k])
				{
					Distance++;
					
				}
			}
			
			if(Distance<MinDistance)
			{
				MinDistance=Distance;
			}
			
			
		}
		
		
		SumDistance=SumDistance+MinDistance;
	}
	
	return SumDistance;
}


char* BruteForceMedianSearch(char **DNA, int t, int n, int l)
{
	int  *a = new int[l];
	char *s = new char[l];
	char *mutif = new char[l];
	int endflag=0,totaldistance=0,BestDistance=INFTY;
	
	for (int i=0; i<l; i++) {
		a[i]=0;
	}
	
	while (1)
	{
		
		//cout << "while";
		for (int i=0; i<l; i++)
		{
			if (a[i]==0) 
			{
				s[i]='a';
			}else if(a[i]==1){
				s[i]='c';
			}else if (a[i]==2) {
				s[i]='g';
			}else if(a[i]==3) {
				s[i]='t';
			}
			
			
		}
		
		totaldistance=TotalDistance(s,DNA,t,n,l);
		if(totaldistance<BestDistance)
		{
			BestDistance=totaldistance;
			for(int i=0; i<l; i++)
			{
				mutif[i]=s[i];
				//cout << mutif[i];
			}
			//cout << BestDistance<<endl;
			//cout << endl;
		}
		
		for(int i=(l-1);i>=0;i--)      //calculate the next leaf
		{
			if(a[i]==3)
			{
				a[i]=0;
				if (i==0) {
					endflag=1;
				}
				
			}else
			{
				a[i]++;
				break;
			}
			
		}
		
		
		if(endflag==1)   //if go over all the leaves, break loop
			break;
		
	}
	
	return mutif;    //return the best mutif
	
}


int TotalDistance2(char *s, char **DNA, int t, int n, int si)
{
    char *v = new char[si+1];
	int Distance, MinDistance=INFTY, SumDistance=0;
	
	for (int i=0; i< t; i++)
	{
		MinDistance=INFTY;	
		for (int j=0; j<(n-si); j++)
		{
			Distance=0;
			for(int k=0; k<(si+1); k++)
			{
				v[k]=DNA[i][j+k];
				
				if(v[k]!=s[k])
				{
					Distance++;
					
				}
			}
			
			if(Distance<MinDistance)
			{
				MinDistance=Distance;
			}
			
			
		}
		
		
		SumDistance=SumDistance+MinDistance;
	}
	
	return SumDistance;
}


char* BranchAndBoundMedianSearch(char **DNA, int t, int n, int l)
{
	int  *a = new int[l];
	char *s = new char[l];
	char *mutif = new char[l];
	int endflag=0,optimisticDistance=0,BestDistance=INFTY,si=0;
	
	for (int i=0; i<l; i++) {
		a[i]=0;
	}
	
	while (1)
	{
		
		
		
		for (int i=0; i<l; i++)
		{
			if (a[i]==0) 
			{
				s[i]='a';
			}else if(a[i]==1){
				s[i]='c';
			}else if (a[i]==2) {
				s[i]='g';
			}else if(a[i]==3) {
				s[i]='t';
			}
		}
		
		
			optimisticDistance=TotalDistance2(s,DNA,t,n,si);
			
			if(si!=(l-1))
			{
				if(optimisticDistance>BestDistance)
				{
					for(int i=si;i>=0;i--)      //calculate the next leaf
					{
						if(a[i]==3)
						{
							a[i]=0;
							if (i==0) {
								endflag=1;
							}
							
						}else
						{
							a[i]++;
							break;
						}
						
					}
					for(int j=si+1; j<l; j++)
					{
						a[j]=0;
					}
					
				}
				else {
				    si++;
					a[si]=0;
				     }


			}else
			{	
				if(optimisticDistance<BestDistance)
		       {
			
			     BestDistance=optimisticDistance;
			      for(int i=0; i<l; i++)
			      {
				     mutif[i]=s[i];
				//cout << mutif[i];
			      }
			
		         }
		
		     for(int i=(l-1);i>=0;i--)      //calculate the next leaf
		       {
			      if(a[i]==3)
			        {
				    si--;
				    //a[i]=0;
				     if (i==0) {
					endflag=1;
				     }
				
			       }else
			       {
				     a[i]++;
					 break;
			       }
			
		        }
		
			}
		if(endflag==1)   //if go over all the leaves, break loop
			break;
		
	}
	
	return mutif;    //return the best mutif
	
}







int main () {
    //***************************************
	//change the t, n, l ,ntimes here
	//****************************************
	int t=4, n=57, l=3, ntimes=1;
    char *mutif=new char [l];
	char** DNA = new char*[t];
	
	clock_t start, finish; 
	double duration; 
	
	for(int i=0; i<t; i++)
	{
	    DNA[i]=new char [n];
		
	}
	
	/*cout<<"please input the value of t:";
	cin>>t;
	cout<<endl;
	cout << "please input the value of l: ";
	cin>>l;
	cout<<endl;
	*/
	cout << "t="<<t<<",　l="<<l<<endl;
	
	ifstream fin("/data.txt");
	
	for(int i=0;i<t;i++)
	{
        for(int j=0;j<n;j++) 
		{
			fin>>DNA[i][j];
			//cout << DNA[i][j];
		}
		//cout<<endl;
	}
    fin.close();
	
	

	     //cout << "please input the times you want to run: ";
		 //cin>>ntimes;
		 //cout << endl;
		
		   
		    start   =   clock();
	        for(int i=0;i<ntimes;i++)
			mutif = BruteForceMotifSearchAgain(DNA, t, n, l);
			finish   =   clock();   
		   
	       duration = (double)(finish-start)/CLOCKS_PER_SEC;
	        duration=duration/ntimes;
	        cout<< "Average time for BruteForceMotifSearchAgain is "<<duration<<endl;

	        cout << "the motif is ";
	        cout<<"'";
	        for(int i=0; i<l; i++)
		    cout << mutif[i];
		    cout<<"'";
		    cout << endl;
		    cout << endl;

		   
		   start   =   clock();
	       for(int i=0;i<ntimes;i++)
		   mutif = BranchAndBoundMotifSearch(DNA, t, n, l);
	      finish   =   clock();   
	      
	      duration = (double)(finish-start)/CLOCKS_PER_SEC;
	      duration=duration/ntimes;
	       cout << "Average time for BranchAndBoundMotifSearch is "<<duration;
           cout << endl;
		   cout << "the motif is ";
		   cout<<"'";
		   for(int i=0; i<l; i++)
			   cout << mutif[i];
		   cout<<"'";
		   cout << endl;
		   cout << endl;

  
	  
		   
		   start   =   clock();
	      for(int i=0;i<ntimes;i++)
		  mutif = BruteForceMedianSearch(DNA, t, n, l);
	     finish   =   clock();   
	     
	     duration = (double)(finish-start)/CLOCKS_PER_SEC;
	      duration=duration/ntimes;
	       cout << "Average time for BruteForceMedianSearch is "<<duration;
           cout << endl;
	       cout << "the motif is ";

		   cout<<"'";
		   for(int i=0; i<l; i++)
			   cout << mutif[i];
		   cout << "'";
		   cout << endl;
		   cout << endl;

		   
		   start  =   clock();
	      for(int i=0;i<ntimes;i++)
		  mutif = BranchAndBoundMedianSearch(DNA, t, n, l);
	      finish   =   clock();   
          duration = (double)(finish-start)/CLOCKS_PER_SEC;
	      duration=duration/ntimes;
	       cout << "Average time for BranchAndBoundMedianSearch is "<<duration;
	       cout << endl;

		   cout << "the motif is ";
           cout<<"'";
		   for(int i=0; i<l; i++)
			   cout << mutif[i];
		   cout<<"'";
		   cout << endl;
		   cout << endl;

	
	
	return 0;
	
}
