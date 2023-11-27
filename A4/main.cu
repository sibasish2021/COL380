#include<bits/stdc++.h>
#include<fstream>
#define endl "\n"
#define uint uint32_t
#define int int32_t
#define maxvalue UINT32_MAX
using namespace std;

__device__ __managed__ int n,m,k1,k2;
// __device__ __managed__ uint64_t maxvalue=((uint64_t)(1e32)-1);

void inputParameters(int* klocal,string fileName)
{
    ifstream f1(fileName);
    f1.read((char *)&n,sizeof(int));
    f1.read((char *)&m,sizeof(int));
    f1.read((char *)klocal,sizeof(int));
    f1.close();
}

void InputData(string fileName,int* matrix,uint* matrixdata,int* klocal)
{
    ifstream f1(fileName);
    f1.read((char *)&n,sizeof(int));
    f1.read((char *)&m,sizeof(int));
    f1.read((char *)klocal,sizeof(int));

    // cout<<"values of n,m,k are"<<n<<" "<<m<<" "<<*klocal<<endl;
    int k=*klocal;
    

    for(int i=0;i<n/m;i++)
    {
        for(int j=0;j<n/m;j++)
        {
            matrix[i*n/m+j]=-1;
        }
    }

    for(int i=0;i<k;i++)
    {
        int posx,posy;
        f1.read((char *)&posx,sizeof(int));
        f1.read((char *)&posy,sizeof(int));
        matrix[posx*n/m+posy]=i;
        // cout<<"Check3"<<endl;
        for(int x=0;x<m;x++)
        {
            for(int y=0;y<m;y++)
            {
                uint16_t temp;
                f1.read((char *)&temp,sizeof(uint16_t));
                matrixdata[i*m*m+x*m+y]=(uint) temp;
                // if(posx==11 && posy==19)
                //     cout<<"Value at position "<< x <<" "<<y<<" is "<<matrixdata[i*m*m+x*m+y]<<endl;
            }
        }
    }
    // cout<<"Check4"<<endl;
    f1.close();   
}

void printmatrix(int* matrix,uint* matrixdata,int* k)
{
    cout<<"n,m,k are"<<n<<" "<<m<<" "<<*k<<endl;

    for(int i=0;i<n/m;i++)
    {
        for(int j=0;j<n/m;j++)
        {
            cout<<matrix[i*n/m+j]<<" ";
        }
        cout<<endl;
    }
}


void outputData(string fileName,int* matrix,uint* matrixdata,int* klocal)
{
    // cout<<"values of n,m,k in output are "<<n<<" "<<m<<" "<<*klocal<<endl;
    ofstream f1(fileName,ios::binary);
    f1.write((char *)&n,sizeof(int32_t));
    f1.write((char *)&m,sizeof(int32_t));
    f1.write((char *)klocal,sizeof(int32_t));
    // int k=*klocal;

    for(int i=0;i<(n/m)*(n/m);i++)
    {
        if(matrix[i]<0)
            continue;
        int posx=i/(n/m);
        int posy=i%(n/m);    
        f1.write((char *)&posx,sizeof(int32_t));
        f1.write((char *)&posy,sizeof(int32_t));
        for(int j=0;j<m;j++)
        {
            for(int k=0;k<m;k++)
            {
                uint temp=(uint)(matrixdata[matrix[i]*m*m+j*m+k]);
                // if(posx==0 && posy==1)
                //     cout<<"Values are "<<j<<" "<<k<<" "<<temp<<endl; 
                f1.write((char *)&temp,sizeof(uint32_t));
            }
        }
    }
    f1.close();
}

__device__ void blockMult(uint* a,uint* b,uint* c)
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<m;j++)
        {
            uint64_t ans=c[i*m+j];
            for(int k=0;k<m;k++)
            {
                ans+=(uint64_t)(a[i*m+k]*b[k*m+j]);
                if (ans>(uint64_t)maxvalue) 
                {
                    ans=maxvalue;
                }
            }
            c[i*m+j]=ans;
        }
    }
}

__global__ void matrixMult(int* a,int* b,int* c,uint* aData,uint* bData,uint* cData)
{
    int blockNum=blockIdx.x*blockDim.x+threadIdx.x;
    if(blockNum >= (n/m)*(n/m)) 
    {
        return;
    }
    int posx= blockNum/(n/m);
    int posy= blockNum%(n/m);
    for (int i=0;i<n/m;i++)
    {
        if ((a[posx*(n/m)+i]<0)||(b[i*(n/m)+posy]<0)) 
        {
            continue;
        }
        blockMult(aData+(uint)a[posx*(n/m)+i]*m*m,bData+(uint)b[i*(n/m)+posy]*m*m, cData+(uint)c[blockNum]*m*m);
    }
    bool check = true;
    for (int i=0;i<m*m;i++) 
    {
        if(*(cData+(uint)c[blockNum]*m*m+i) != 0)
        {
            check=false;
            break;
        }
    }
    if (check) 
    {
        c[blockNum]=-1;
    }
}

signed main(int argc, char *argv[])
{
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);
    
    string inputFile1(argv[1]);
    string inputFile2(argv[2]);
    string outputFile(argv[3]);

    int* a;int* b;int* c;
    int* ga;int* gb;int* gc;
    uint* aData;uint* bData;uint* cData;
    uint* gaData;uint* gbData;uint* gcData;
    
    inputParameters(&k1,inputFile1);
    a=(int*)calloc(n/m*n/m,sizeof(int));
    aData=(uint*)calloc(m*m*k1,sizeof(uint));
    
    inputParameters(&k2,inputFile2);
    b=(int*)calloc(n/m*n/m,sizeof(int));
    bData=(uint*)calloc(m*m*k2,sizeof(uint));
    
    InputData(inputFile1,a,aData,&k1);
    InputData(inputFile2,b,bData,&k2);
    
    c=(int*)calloc((n/m)*(n/m),sizeof(int));
    cData=(uint*)calloc(n*n,sizeof(uint));
    for(int i=0;i<(n/m)*(n/m);i++)
    {
        c[i]=i;
    }

    cudaMalloc((void**)&ga, (n/m)*(n/m)* sizeof(int));
    cudaMalloc((void**)&gb, (n/m)*(n/m)* sizeof(int));
    cudaMalloc((void**)&gc, (n/m)*(n/m)* sizeof(int));
    
    cudaMalloc((void**)&gaData, k1*m*m* sizeof(uint));
    cudaMalloc((void**)&gbData, k2*m*m* sizeof(uint));
    cudaMalloc((void**)&gcData, n*n* sizeof(uint));

    int dev = -1;
    cudaGetDevice(&dev);
    cudaMemPrefetchAsync(&n, sizeof(int), dev);
    cudaMemPrefetchAsync(&m, sizeof(int), dev);

    cudaMemcpyAsync(ga, a, (n/m)*(n/m)*sizeof(int),cudaMemcpyHostToDevice,0);
    cudaMemcpyAsync(gb, b, (n/m)*(n/m)*sizeof(int),cudaMemcpyHostToDevice,0);
    cudaMemcpyAsync(gc, c, (n/m)*(n/m)*sizeof(int),cudaMemcpyHostToDevice,0);
    
    cudaStream_t stream1;
    cudaStreamCreate(&stream1);
    cudaStream_t stream2;
    cudaStreamCreate(&stream2);
    cudaStream_t stream3;
    cudaStreamCreate(&stream3);

    cudaMemcpyAsync(gaData,aData,k1*m*m*sizeof(uint),cudaMemcpyHostToDevice,stream1);
    cudaMemcpyAsync(gbData,bData,k2*m*m*sizeof(uint),cudaMemcpyHostToDevice,stream2);
    cudaMemcpyAsync(gcData,cData,n*n*sizeof(uint),cudaMemcpyHostToDevice,stream3);
    
    // cudaStreamSynchronize(0);
    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
    cudaStreamSynchronize(stream3);

    int threads_per_block=1024;
    int blocks_per_grid = ((n/m)*(n/m)+threads_per_block-1)/threads_per_block;
    
    dim3 dimGrid(blocks_per_grid,1,1);
    dim3 dimBlock(threads_per_block,1,1);
    
    cudaEvent_t	event1,	event2;
    cudaEventCreate(&event1);	
    cudaEventCreate(&event2);
    cudaEventRecord(event1,0);
    matrixMult<<<dimGrid,dimBlock,0>>>(ga,gb,gc,gaData,gbData,gcData);
    cudaEventRecord(event2,0);
    cudaEventSynchronize(event2);
    float ms;
    cudaEventElapsedTime(&ms,event1,event2);
    // cout<<"Time taken by the kernel is "<<ms<<endl;

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    { 
        printf("Error: %s\n", cudaGetErrorString(err));
    }
    
    // cudaDeviceSynchronize();

    cudaMemcpy(c,gc, (n/m)*(n/m)*sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(cData,gcData,n*n*sizeof(uint),cudaMemcpyDeviceToHost);
    // printBlocksCPU(c,cData,n,m);

    int k3=0;
    for(int i=0;i<n/m;i++)
    {
        for(int j=0;j<n/m;j++)
        {
            if(c[i*(n/m)+j]>=0)
            {
                k3++;
            }
        }
    }

    outputData(outputFile,c,cData,&k3);
    
    cudaFree(ga);
    cudaFree(gb);
    cudaFree(gc);
    cudaFree(gaData);
    cudaFree(gbData);
    cudaFree(gcData);
    
    free(a);
    free(b);
    free(c);
    free(aData);
    free(bData);
    free(cData);

    return 0;
}