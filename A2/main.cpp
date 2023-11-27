#include<bits/stdc++.h>
#include<mpi.h>
#include<fstream>
#include<chrono>
#include<iterator>
using namespace std;
using namespace std::chrono;
#define endl "\n"
#define  F first
#define  S second

map<int,set<int>> edges;
map<int,int> degree;
map<pair<int,int>,int> support;

char *inputpath, *headerpath, *outputpath;
int verbose;
int startk,endk;
int taskid;


int numproc;
int myrank;

int n,m;
string headerfile,inputfile;
int* isAvailable;
map<pair<int,int>,int> edge_available;

void inputData(string filename)
{
    ifstream f1(filename);
    f1.read((char *)&n,sizeof(int32_t));
    f1.read((char *)&m,sizeof(int32_t));

    for(int i=0;i<n;i++)
    {
        int tid,tdeg;
        f1.read((char *)&tid,sizeof(int32_t));
        f1.read((char *)&tdeg,sizeof(int32_t));
        degree[tid]=tdeg;
        for(int j=0;j<tdeg;j++)
        {
            int temp;
            f1.read((char *)&temp,sizeof(int32_t));
            // if(temp>tid)
            {
                edges[tid].insert(temp);
                edge_available.insert({{tid,temp},0});
            }
        }
    }
    f1.close();
}

void prefilter(int k)
{
    isAvailable=(int *)malloc(sizeof(int)*n);
    for(int i=0;i<n;i++)
    {
        isAvailable[i]=1;
    }
    stack<int> deletable;
    for(auto it:degree)
    {
        if(it.S>0 && it.S<k-1)
        {
            deletable.push(it.F);
            isAvailable[it.F]=0;
        }
    }
    while(!deletable.empty())
    {
        int node=deletable.top();
        deletable.pop();
        degree.erase(node);
        for(auto it1=edges[node].begin();it1!=edges[node].end();it1++)
        {
            int anode=*it1;
            support[{node,anode}]=-1;
            support[{anode,node}]=-1;
            {
                degree[anode]--;
                if(degree[anode]<k-1)
                {   
                    if(isAvailable[anode]==1){
                        deletable.push(anode);
                        isAvailable[anode]=0;
                    }
                }
            }
            edges[anode].erase(node);
        }
        edges.erase(node);
    }
    free(isAvailable);
}

stack<pair<int,int>> del;

void initialise(int k,int startk1)
{
    if(k==startk1+2)
    {
        
        vector<int> tempsup;
        vector<pair<int,int>> tempedges;
        tempedges.clear();
        
        for(auto it:edges)
        {
            for(auto it1=edges[it.F].begin();it1!=edges[it.F].end();it1++)
            {
                int tempn=*it1;
                if(it.F<tempn)
                {
                    tempedges.push_back({it.F,tempn});
                }
            }
        }
        int edgespercell=(int)ceil(tempedges.size()/numproc);
        for(int i=myrank;i<tempedges.size();i+=numproc)
        {
            int v1=tempedges[i].F;
            int v2=tempedges[i].S;
            int sup=0;
            for(auto it=edges[v2].begin();it!=edges[v2].end();it++)
            {
                int tempn=*it;
                if(edges[v1].count(tempn))
                {
                    sup++;
                }
            }
            tempsup.push_back(sup);
        }
        vector<int> tempsupall(tempedges.size(),0);
        if(myrank==0)
        {
            for(int i=1;i<numproc;i++)
            {
                vector<int> suprecv(edgespercell+numproc,0);
                MPI_Recv(suprecv.data(),edgespercell+numproc,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for(int j=i;j<tempedges.size();j+=numproc)
                {
                    int index=(int)j/numproc;
                    tempsupall[j]=suprecv[index];
                }
            }
            for(int j=0;j<tempedges.size();j+=numproc)
            {
                tempsupall[j]=tempsup[j/numproc];
            }
        }
        else
        {
            MPI_Send(tempsup.data(),tempsup.size(),MPI_INT,0,0,MPI_COMM_WORLD);
            return;
        }
        if(myrank==0)
        {
            for(int i=0;i<tempedges.size();i++)
            {
                int v1=tempedges[i].F;
                int v2=tempedges[i].S;
                support[{v1,v2}]=tempsupall[i];
                support[{v2,v1}]=tempsupall[i];
            }
        }
    }
    for(auto it:support)
    {
        if(it.S<k-2)
        {
            del.push(it.F);
            edge_available[it.F]=1;
        }
    }
}

void filteredges(int k)
{
    while(!del.empty())
    {
        pair<int,int> edge=del.top();
        del.pop();
        if(edge.F>edge.S)
        {
            continue;
        }
        int v1=edge.F;
        int v2=edge.S;
        for(auto it2=edges[v2].begin();it2!=edges[v2].end();it2++)
        {
            int n=*it2;
            if(edges[v1].count(n))
            {
                support[{v1,n}]--;
                support[{v2,n}]--;
                support[{n,v1}]--;
                support[{n,v2}]--;
                if(support[{v1,n}]<k-2 && !edge_available[{v1,n}])
                {
                    del.push({v1,n});
                    del.push({n,v1});
                    edge_available[{v1,n}]=1;
                    edge_available[{n,v1}]=1;
                }
                if(support[{v2,n}]<k-2 && !edge_available[{v2,n}])
                {
                    del.push({v2,n});
                    del.push({n,v2});
                    edge_available[{v2,n}]=1;
                    edge_available[{n,v2}]=1;
                }
            }
            edges[edge.F].erase(edge.S);
            edges[edge.S].erase(edge.F);
            support[{v1,v2}]=-1;
            support[{v2,v1}]=-1;
        }
        
    }
}

bool mintruss(int k)
{
    if(myrank!=0)
    return true;
    while(true)
    {
        bool check=true;//Checks if no more edges to be deleted
        for(auto it:support)
        {
            if(it.S==0)
            {
                edges[it.F.F].erase(it.F.S);
                edges[it.F.S].erase(it.F.F);
                it.S--;
            }
            if(it.S>0 && it.S<k-2 && it.F.F<it.F.S )//If there is an edge required to be deleted
            {
                check=false;
                int v1=it.F.F;
                int v2=it.F.S;
                degree[v1]--;
                degree[v2]--;
                edges[v1].erase(v2);
                edges[v2].erase(v1);
                
                for(auto it1=edges[v2].begin();it1!=edges[v2].end();it1++)
                {
                    int n=*it1;
                    if(edges[v1].count(n))
                    {
                        support[{v1,n}]--;
                        support[{v2,n}]--;
                        support[{n,v1}]--;
                        support[{n,v2}]--;
                    }
                }
                support[{v1,v2}]=0;
                support[{v2,v1}]=0;
            }
        }
        if(check)
        {
            break;
        }
    }
    int x1=0;
    int x2=0;
    x1=0;
    x2=0;
    for(auto it:support)
    {
        if(it.S>0)
        {
            x1++;
        }
    }
    for(auto it:edges)
    {
        degree[it.F]=it.S.size();
        x2+=it.S.size();
    }
    return (x1==0)&&(x2==0);
}

vector<vector<int>> ans;
vector<int> path;
vector<int> visited;

void dfa(int node)
{
    bool check =true;
    for(auto it=edges[node].begin();it!=edges[node].end();it++)
    {
        int nbr=*it;
        if(!visited[nbr])
        {
            visited[nbr]=1;
            path.push_back(nbr);
            dfa(nbr);
        }
    }
}

bool compare(vector<int> v1,vector<int> v2)
{
    return v1.size()>v2.size();
}

int size_support()
{
    int x1=0;
    for(auto it:support)
    {
        if(it.S>0)
        {
            x1++;
        }
    }
    return x1;
}

int size_edges()
{
    int x2=0;
    for(auto it:edges)
    {
        degree[it.F]=it.S.size();
        x2+=it.S.size();
    }
    return x2;
}

int main(int argc,char ** argv)
{
    taskid=atoi(argv[1]+9);
    inputpath=argv[2]+12;
    headerpath=argv[3]+13;
    outputpath=argv[4]+13;
    verbose=atoi(argv[5]+10);
    startk=atoi(argv[6]+9);
    endk=atoi(argv[7]+7);

    bool check=false;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    {
        inputData(inputpath);        
        ofstream outfile;
        if(myrank==0)
        {
            outfile.open(outputpath);
        }
        for(int index=startk;index<=endk;index++)
        {
            int x1,x2;

            if(check && myrank==0)
            {
                outfile<<0<<endl;
                continue;
            }
            int kvalue=index+2;

            prefilter(kvalue);
            initialise(kvalue,startk);
            filteredges(kvalue);
            check=mintruss(kvalue);

            outfile<<!check<<endl;

            if(!check && verbose==1 && myrank==0)
            {
                ans.clear();
                path.clear();
                visited.resize(n);
                for(int i=0;i<n;i++)
                visited[i]=0;
                for(int i=0;i<n;i++)
                {
                    if(!visited[i])
                    {
                        visited[i]=1;
                        path.push_back(i);
                        dfa(i);
                        if(path.size()>2)
                        ans.push_back(path);
                        path.clear();
                    }
                }
                outfile<<ans.size()<<endl;
                for(int i=0;i<ans.size();i++)
                {
                    sort(ans[i].begin(),ans[i].end());
                    for(int j=0;j<ans[i].size();j++)
                    {
                        outfile<<ans[i][j]<<" ";
                    }
                    outfile<<endl;
                }
            }
        }
        outfile.close();
    }
    MPI_Finalize();
    return 0;
}