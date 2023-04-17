
#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
using namespace std;



float GaussSeidel(int x,int y,float** Te2) //Function for GS iteration
{
float T;
float R=0;
  for(int i=1;i<y-1;i++)
    for(int j=1;j<x-1;j++)
    {
    T=(Te2[i+1][j]+Te2[i-1][j]+Te2[i][j+1]+Te2[i][j-1])/4;
    //if(abs(T-Te2[i][j])>R)
     R=R+abs(T-Te2[i][j]);
    Te2[i][j]=T;
    }

return R;
}

float PSOR(int x,int y,float** Te2,float omega) //Function for point GS iteration
{
float T;
float R=0;
    for(int i=1;i<y-1;i++)
    for(int j=1;j<x-1;j++)
    {
    T=(1-omega)*Te2[i][j]+omega*(Te2[i+1][j]+Te2[i-1][j]+Te2[i][j+1]+Te2[i][j-1])/4;
    //if(abs(T-Te2[i][j])>R)
     R=R+abs(T-Te2[i][j]);
    Te2[i][j]=T;
    }
return R;
}



float lineX(int x, int y, float** Te2) //Function for line GS iteration in X direction
{
    float* Tr=new float[x];
    float* a=new float[x];
    float* d=new float[x];
    float* b=new float[x];
    float* r=new float[x];

    float T;
    float R=0;

    for(int i=1;i<y-1;i++)
    {
     for(int j=0;j<x;j++)
     {


       if((j>0)&&(j<x-1))
        {
       d[j]=4;
       a[j]=-1;
       b[j]=-1;
       r[j]=Te2[i+1][j]+Te2[i-1][j];
        }
       else
        {

         if(j==0)
         {
         d[j]=1;
         b[j]=0;
         a[j]=0;
         r[j]=Te2[i][j];
         }
         else if(j==x-1)

         {
         d[j]=1;
         a[j]=0;
         b[j]=0;
         r[j]=Te2[i][j];

         }

       }

      }
     for(int j=1;j<x;j++)
      {
        d[j]=d[j]-a[j-1]*b[j]/d[j-1];
        r[j]=r[j]-r[j-1]*b[j]/d[j-1];
      }

      Te2[i][x-1]=r[x-1]/d[x-1];
      for(int j=x-2;j>0;j--)
        {
        T=(r[j]-(a[j]*Te2[i][j+1]))/d[j];
        //if(abs(T-Te2[i][j])>R)
        R=R+abs(T-Te2[i][j]);
        Te2[i][j]=T;

        }
     }

 return R;
}
float lineY(int x, int y, float** Te2)  //Function for line GS iteration in Y direction
{
    float* Tr=new float[y];
    float* a=new float[y];
    float* d=new float[y];
    float* b=new float[y];
    float* r=new float[y];

    float T;
    float R=0;
    for(int j=1;j<x-1;j++)
    {
     for(int i=0;i<y;i++)
     {

       if((i>0)&&(i<=y-2))
        {
        d[i]=4;
        a[i]=-1;
        b[i]=-1;
        r[i]=Te2[i][j+1]+Te2[i][j-1];
         }

       else
        {

         if(i==0)
         {
         d[i]=1;
         b[i]=0;
         a[i]=0;
         r[i]=Te2[i][j];
         }
         else if(i==y-1)
         {
         d[i]=1;
         a[i]=0;
         b[i]=-1;
         r[i]=Te2[i][j];
         }

       }

      }

     for(int i=1;i<y-1;i++)
      {
        d[i]=d[i]-a[i-1]*b[i]/d[i-1];
        r[i]=r[i]-r[i-1]*b[i]/d[i-1];

      }

      Te2[y-1][j]=r[y-1]/d[y-1];// Determination of the r[y-2] was the issue causing wrong results.

      for(int i=y-2;i>0;i--)
        {
        T=(r[i]-(a[i]*Te2[i+1][j]))/d[i];
        //if(abs(T-Te2[i][j])>R)
        R=R+abs(T-Te2[i][j]);
        Te2[i][j]=T;
        }
     }
 return R;
}
float lineSOR_X(int x, int y, float** Te2, float omega) //Function for line GS iteration in X direction
{
    float* Tr=new float[x];
    float* a=new float[x];
    float* d=new float[x];
    float* b=new float[x];
    float* r=new float[x];

    float T;
    float R=0;
    for(int i=1;i<y-1;i++)
    {
     for(int j=0;j<x;j++)
     {


       if((j>0)&&(j<x-1))
        {
       d[j]=4;
       a[j]=-omega;
       b[j]=-omega;
       r[j]=(1-omega)*4*Te2[i][j]+omega*(Te2[i+1][j]+Te2[i-1][j]);
        }
       else
        {

         if(j==0)
         {
         d[j]=1;
         b[j]=0;
         a[j]=0;
         r[j]=Te2[i][j];
         }
         else if(j==x-1)

         {
         d[j]=1;
         a[j]=0;
         b[j]=0;
         r[j]=Te2[i][j];

         }

       }

      }
     for(int j=1;j<x;j++)
      {
        d[j]=d[j]-a[j-1]*b[j]/d[j-1];
        r[j]=r[j]-r[j-1]*b[j]/d[j-1];
      }

      Te2[i][x-1]=r[x-1]/d[x-1];
      for(int j=x-2;j>0;j--)
        {
        T=(r[j]-(a[j]*Te2[i][j+1]))/d[j];
        //if(abs(T-Te2[i][j])>R)
        R=R+abs(T-Te2[i][j]);
        Te2[i][j]=T;
        }
     }
 return R;
}
float lineSOR_Y(int x, int y, float** Te2,float omega)  //Function for line GS iteration in Y direction
{
    float* Tr=new float[y];
    float* a=new float[y];
    float* d=new float[y];
    float* b=new float[y];
    float* r=new float[y];

    float T;
    float R=0;
    for(int j=1;j<x-1;j++)
    {
     for(int i=0;i<y;i++)
     {

       if((i>0)&&(i<=y-2))
        {
        d[i]=4;
        a[i]=-omega;
        b[i]=-omega;
        r[i]=(1-omega)*4*Te2[i][j]+omega*(Te2[i][j+1]+Te2[i][j-1]);
         }

       else
        {

         if(i==0)
         {
         d[i]=1;
         b[i]=0;
         a[i]=0;
         r[i]=Te2[i][j];
         }
         else if(i==y-1)
         {
         d[i]=1;
         a[i]=0;
         b[i]=-1;
         r[i]=Te2[i][j];
         }

       }

      }

     for(int i=1;i<y-1;i++)
      {
        d[i]=d[i]-a[i-1]*b[i]/d[i-1];
        r[i]=r[i]-r[i-1]*b[i]/d[i-1];

      }

      Te2[y-1][j]=r[y-1]/d[y-1];
      for(int i=y-2;i>0;i--)
        {
        T=(r[i]-(a[i]*Te2[i+1][j]))/d[i];
        //if(abs(T-Te2[i][j])>R)
        R=R+abs(T-Te2[i][j]);
        Te2[i][j]=T;
        }
     }
 return R;
}


float Unsteady(int x,int y,float** Te2) //Function for iteration
{

float R=0;
float** Te=new float*[y];
  for(int i=0;i<y;i++)
      Te[i]=new float[x];

  for(int i=1;i<y-1;i++)
    for(int j=1;j<x-1;j++)
    {
    Te[i][j]=(Te2[i+1][j]+Te2[i-1][j]+Te2[i][j+1]+Te2[i][j-1])/4;
     R=R+abs(Te[i][j]-Te2[i][j]);
    }

  for(int i=1;i<y-1;i++)
     for(int j=1;j<x-1;j++)
       Te2[i][j]=Te[i][j];
return R;
}

float UnsteadyW(int x,int y,float** Te2, float omega) //Function for iteration
{

float R=0;
float** Te=new float*[y];
  for(int i=0;i<y;i++)
      Te[i]=new float[x];

  for(int i=1;i<y-1;i++)
    for(int j=1;j<x-1;j++)
    {
    Te[i][j]=(1-omega)*Te2[i][j]+omega*(Te2[i+1][j]+Te2[i-1][j]+Te2[i][j+1]+Te2[i][j-1])/4;
     R=R+abs(Te[i][j]-Te2[i][j]);
    }

  for(int i=1;i<y-1;i++)
     for(int j=1;j<x-1;j++)
       Te2[i][j]=Te[i][j];
return R;
}



float** alocate(int x,int y)
{
    float** Te2=new float*[y];
    for(int i=0;i<y;i++)
    {
    Te2[i]=new float[x];
     }
     return Te2;
}




void dealocate(int x,float** Te2)
{
   for(int i=0;i<x;i++)
      {
      delete[] Te2[i];
      }
      delete[] Te2;

      Te2=NULL;
}

void boundaryFun(int x,int y, float** Te2,float L,float R,float T,float B)
{

for(int i=0;i<y;i++)
    {
     for(int j=0;j<x;j++)
        {
        if(i==0)
         Te2[i][j]=T;
        else if(i==y-1)
         Te2[i][j]=B;
             else if(j==0)
              Te2[i][j]=L;
                  else if(j==x-1)
                   Te2[i][j]=L;
                        else
                         Te2[i][j]=0;
        }
     }

}
void printResults(int x,int y, float** Te2)
{
 ofstream outfile;
 outfile.open("output.txt");
 for(int i=0;i<y;i++)
    {
     for(int j=0;j<x;j++)
        {
        cout<<Te2[i][j]<<"\t";
        outfile<<Te2[i][j]<<"\t";
        }
     cout<<endl;
     outfile<<endl;
     }
     outfile.close();
}



int main()
{
    ifstream infile;
    string word[11];
    float Num[11];
    infile.open("input.txt");
    int i=0;
    while(infile.good())
      {
       infile>>word[i]>>Num[i];
       cout<<word[i]<<" "<<Num[i]<<endl;
       i++;
      }
     infile.close();
     float L=Num[0];
     float H=Num[1];
     int Nx=Num[2];
     int Ny=Num[3];
     float Tl=Num[4],Tr=Num[5],Tt=Num[6],Tb=Num[7];
     float dx=L/(Nx-1);
     float dy=H/(Ny-1);
     float alpha=Num[8];
     float R=Num[9];
     float w;
     float** T2;
     T2=alocate(Nx,Ny); // Calling Function to allocate Memory
     boundaryFun(Nx,Ny,T2,Tl,Tr,Tt,Tb);// Calling Function to Initialize the boundary
     float residue=1;
     int itr=0;
     int x;
     cout<<"\nMenu\n1.Gauss Seidel\n2.Line Gauss Seidel\n3.ADI Method\n4.PSOR\n5.LSOR\n6.ADI with relaxation\n7.Parabolic Formulation\n\nChose solver(1,2,3,4,5,6,7) :";
     cin>>x;
     double start;
     double duration;
     float time=0;
     switch(x)
     {
      case 1:
             start=clock();
             while(residue>R)
             {
             residue=GaussSeidel(Nx,Ny,T2);
             itr++;
             }
             duration=clock()-start;
             break;
      case 2:
             start=clock();
             while(residue>R)
             {
              residue=lineX(Nx,Ny,T2);
              itr++;
             }
             duration=clock()-start;
              break;
      case 3:
             start=clock();
             while(residue>R)
             {
              lineX(Nx,Ny,T2);
              residue=lineY(Nx,Ny,T2);
              itr++;
             }
             duration=clock()-start;
              break;
      case 4:
             cout<<"\nEnter relaxation factor :";
             cin>>w;
             start=clock();
             while(residue>R)
             {
             residue=PSOR(Nx,Ny,T2,w);
             itr++;
             }
             duration=clock()-start;
             break;
      case 5:
             cout<<"\nEnter relaxation factor :";
             cin>>w;
             start=clock();
             while(residue>R)
             {
             residue=lineSOR_X(Nx,Ny,T2,w);
             itr++;
             }
             duration=clock()-start;
             break;
      case 6:
             cout<<"\nEnter relaxation factor :";
             cin>>w;
             start=clock();
             while(residue>R)
             {
             lineSOR_X(Nx,Ny,T2,w);
             residue=lineSOR_Y(Nx,Ny,T2,w);
             itr++;
             }
             duration=clock()-start;
             break;
      case 7:
             start=clock();
             while(residue>R)
             {
             residue=Unsteady(Nx,Ny,T2);
             itr++;
             }
             duration=clock()-start;
             break;
      case 8:
             cout<<"\nEnter relaxation factor :";
             cin>>w;
             start=clock();
             while(residue>R)
             {
             residue=UnsteadyW(Nx,Ny,T2,w);
             itr++;
             }
             duration=clock()-start;
             break;

      }
     cout<<endl;
     printResults(Nx,Ny,T2);
     if(x==7)
       {
       time=itr*0.222538;
       cout<<"\nTime to reach steadystate :"<<time<<" Seconds";
       }
     cout<<endl<<"Number of iterations :"<<itr;
     dealocate(Nx,T2);
     cout<<endl<<"Time taken by the machine to solve :"<<duration*1000000000<<" nano seconds";
     cout<<endl;
   return 0;
}
