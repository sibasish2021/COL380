#include<bits/stdc++.h>
using namespace std;
#define int int32_t
#define endl "\n"
struct Block
{
    int posx;
    int posy;
    vector<vector<int>> block;
};

vector<Block> InputData(string fileName,int &n,int &m,int &k)
{
    vector<Block> v;
    ifstream f1(fileName);
    f1.read((char *)&n,sizeof(int32_t));
    f1.read((char *)&m,sizeof(int32_t));
    f1.read((char *)&k,sizeof(int32_t));
    for(int i=0;i<k;i++)
    {
        Block b;
        f1.read((char *)&b.posx,sizeof(int32_t));
        f1.read((char *)&b.posy,sizeof(int32_t));
        b.block.resize(m);
        for(int x=0;x<m;x++)
        {
            b.block[x].resize(m);
        }
        for(int x=0;x<m;x++)
        {
            for(int y=0;y<m;y++)
            {
                int16_t temp;
                f1.read((char *)&temp,sizeof(int16_t));
                b.block[x][y]=(int) temp;
            }
        }
        v.push_back(b);
    }
    f1.close();
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

void compareBinary(string file1,string file2)
{
    vector<Block> b1,b2;
    int n1,n2,m1,m2,k1,k2;
    b1=InputData(file1,n1,m1,k1);
    b2=InputData(file2,n2,m2,k2);
    sort(b1.begin(),b1.end(),compare);
    sort(b2.begin(),b2.end(),compare);
    if(n1!=n2)
    {
        cout<<"n mismatch"<<n1<<" "<<n2<<endl;
        return;
    }
    if(m1!=m2)
    {
        cout<<"m mismatch"<<m1<<" "<<m2<<endl;
        return;
    }
    if(n1!=n2)
    {
        cout<<"k mismatch"<<k1<<" "<<k2<<endl;
        return;
    }
    if(b1.size()!=b2.size())
    {
        cout<<"Block size mismatch"<<b1.size()<<" "<<b2.size()<<endl;
        return;
    }
    for(int i=0;i<b1.size();i++)
    {
        if(b1[i].posx!=b2[i].posx||b1[i].posy!=b2[i].posy)
        {
            cout<<"Index of blocks not matching"<<endl;
            cout<<"First index "<<b1[i].posx<<" "<<b1[i].posy<<endl<<"SecondIndex"<<b2[i].posx<<" "<<b2[i].posy<<endl;
            return;
        }
        for(int j=0;j<m1;j++)
        {
            for(int k=0;k<m1;k++)
            {
                if(b1[i].block[j][k]!=b2[i].block[j][k])
                {
                    cout<<"Elements do not match"<<endl;
                    cout<<"First element"<<b1[i].block[j][k]<<endl<<"Second element"<<b2[i].block[j][k]<<endl;
                    return;
                }
            }
        }
    }
    cout<<"Correct Matching"<<endl;
}

int main(int argc, char *argv[])
{
    compareBinary(argv[1],argv[2]);
}

