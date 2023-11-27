#include "library.hpp"
#include<bits/stdc++.h>
#include<fstream>
#include<omp.h>
#define endl "\n"
#define int int32_t
using namespace std;
int n,m,k;
const int maxval=(1<<16)-1;
struct Block
{
    int posx;
    int posy;
    vector<int> block;
};
vector<Block> boutput,binput;

vector<Block> InputData(string fileName)
{
    vector<Block> v;
    ifstream f1(fileName);
    f1.read((char *)&n,sizeof(int));
    f1.read((char *)&m,sizeof(int));
    f1.read((char *)&k,sizeof(int));
    for(int i=0;i<k;i++)
    {
        Block b;
        f1.read((char *)&b.posx,sizeof(int));
        f1.read((char *)&b.posy,sizeof(int));
        b.block.resize(m*m);
        for(int x=0;x<m;x++)
        {
            for(int y=0;y<m;y++)
            {
                int8_t temp;
                f1.read((char *)&temp,sizeof(int8_t));
                b.block[x*m+y]=(int) temp;
            }
        }
        v.push_back(b);
    }
    f1.close();
    return v;
}
void OutputData(string fileName)
{
    ofstream f1(fileName,ios::binary);
    f1.write((char *)&n,sizeof(int32_t));
    f1.write((char *)&m,sizeof(int32_t));
    f1.write((char *)&k,sizeof(int32_t));
    int s=boutput.size();
    for(int i=0;i<s;i++)
    {
        f1.write((char *)&boutput[i].posx,sizeof(int32_t));
        f1.write((char *)&boutput[i].posy,sizeof(int32_t));
        for(int j=0;j<m;j++)
        {
            for(int k=0;k<m;k++)
            {
                int16_t temp=(int16_t)boutput[i].block[j*m+k];
                f1.write((char *)&temp,sizeof(int16_t));
            }
        }
    }
    f1.close();
}
vector<Block> transpose(vector<Block> b)
{
    int size=b.size();
    vector<Block> v;
    #pragma omp parallel
    {
        #pragma omp  for
        for(int i=0;i<size;i++)
        {
            Block b1;
            b1.block.resize(m*m);

            b1.posx=b[i].posy;
            b1.posy=b[i].posx;
            for(int j=0;j<m;j++)
            {
                for(int k=0;k<m;k++)
                {
                    b1.block[k*m+j]=b[i].block[j*m+k];
                }
            }
            #pragma omp critical
            v.push_back(b1);
        }
    }
    return v;
}
bool compare(Block b1,Block b2)
{
    if(b1.posx!=b2.posx)
    {
        return b1.posx<b2.posx;
    }
    else
    {
        return b1.posy<b2.posy;
    }
}
void multiply()
{
    vector<Block> bt=transpose(binput);
    vector<Block> binput1;
    #pragma omp parallel private(binput1) shared(bt) 
    {
        #pragma omp for
        for(int i=0;i<bt.size();i++)
        {
            if(bt[i].posx!=bt[i].posy)
            {
                binput1.push_back(bt[i]);
            }
        }

        #pragma omp critical
        {
            binput.insert(binput.end(),binput1.begin(),binput1.end());
        }
    }
    vector<Block> b1=binput;
    vector<Block> bmatrix;
    map<int,vector<Block>> m1;
    for(int i=0;i<binput.size();i++)
    {
        m1[binput[i].posx].push_back(binput[i]);
    }
    vector<bool> isSet;
    bmatrix.resize(n/m*n/m);
    isSet.resize(n/m*n/m);
    vector<int> zero(m*m,0);
    for(int i=0;i<n/m*n/m;i++)
    {
        isSet[i]=false;
        bmatrix[i].posx=i/(n/m);
        bmatrix[i].posy=i%(n/m);
    }
    #pragma omp parallel
    {
        #pragma omp master
        {
            for(int i=0;i<b1.size();i++)
            {
                vector<Block> b2=m1[b1[i].posy];
                for(int j=0;j<b2.size();j++)
                {
                    if(b1[i].posy==b2[j].posx)
                    {
                        int x1=b1[i].posx;
                        int y1=b2[j].posy;
                        if(x1>y1)
                        {
                            continue;
                        }
                        if(!isSet[x1*n/m+y1])
                        {
                            bmatrix[x1*n/m+y1].block=zero;
                        }
                        #pragma omp task depend(inout:bmatrix[x1*n/m+y1].block)
                        {
                            for(int a=0;a<m;a++)
                            {
                                for(int c=0;c<m;c++)
                                {
                                    for(int b=0;b<m;b++)
                                    {
                                        bmatrix[x1*n/m+y1].block[a*m+b]=min(maxval,Outer(bmatrix[x1*n/m+y1].block[a*m+b],Inner(b1[i].block[a*m+c],b2[j].block[c*m+b])));
                                    }
                                }
                            }
                        }
                        isSet[x1*n/m+y1]=true;
                    }
                }
            }
        }
    }
    boutput.clear();
    for(int i=0;i<n/m;i++)
    {
        for(int j=0;j<n/m;j++)
        {
            if(isSet[i*n/m+j] && i<=j)
            {
                boutput.push_back(bmatrix[i*n/m+j]);
            }
        }
    }
}
bool isZero(Block b)
{
    for(int i=0;i<m*m;i++)
    {
        if(b.block[i]!=0)
        {
            return false;
        }
    }
    return true;
}
                
vector<Block> nonZero(vector<Block> b)
{
    vector<Block> v;
    for(int i=0;i<b.size();i++)
    {
        if(!isZero(b[i]))
        {
            v.push_back(b[i]);
        }
    }
    return v;
}    

int main(int argc, char *argv[])
{
    omp_set_num_threads(omp_get_num_procs());
    string inputFile=argv[1];
    string outputFile=argv[2];
    binput=InputData(inputFile);
    binput=nonZero(binput);
    multiply();
    boutput=nonZero(boutput);
    k=boutput.size();
    OutputData(outputFile);
    return 0;
}