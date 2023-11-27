#include<bits/stdc++.h>
#include<mpi.h>
#include<omp.h>
#include<vector>

using namespace std;
#define F first 
#define S second
#define endl "\n"

vector<vector<int>> adjacency_list;
vector<pair<int,int>> deg_ordered_vertices;
vector<int> id2owner,my_vertices;
map<pair<int,int>,set<int>> edge_triangles;
set<pair<int,int>> my_edges;
vector<int> offsets;
map<int,set<int>> edges;
ofstream outfile;
vector<int> parent;
map<int,set<int>> connected;

MPI_File gfile,hfile;
MPI_Offset gfsize;
MPI_Status status;
MPI_Info info;
MPI_Datatype packet,pair_type;

char *inputpath, *headerpath, *outputpath;
int verbose;
int startk,endk;
int taskid,p;

int numproc;
int myrank;
int n,m;

struct Packet
{
    int owner;
    int x;
    int y;
    int z;
};

void initialise()
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<adjacency_list[i].size();j++)
        {
            my_edges.insert({i,adjacency_list[i][j]});
        }
    }
    
    for(auto &it:my_edges)
    {
        int v1=it.F;
        int v2=it.S;
        int deg1 = (offsets[v1+1]-offsets[v1]-8)/4;
        int deg2 = (offsets[v2+1]-offsets[v2]-8)/4;
        vector<int> nbr1(deg1);
        vector<int> nbr2(deg2);
        MPI_File_read_at(gfile,offsets[v1]+8,nbr1.data(),deg1,MPI_INT,&status);
        MPI_File_read_at(gfile,offsets[v2]+8,nbr2.data(),deg2,MPI_INT,&status);
        int j = 0; int k =0;
        while(j<deg1 && k<deg2)
        {
            if(nbr1[j]==nbr2[k])
            {
                edge_triangles[{v1,v2}].insert(nbr1[j]);
                j++;
                k++;
            }
            else if(nbr1[j]<nbr2[k])
            {
                j++;
            }
            else
            {
                k++;
            }
        }
    }
}

void inputData()
{
    MPI_Info_create(&info);
    MPI_Info_set( info , "collective_buffering" , "true");
    MPI_File_open(MPI_COMM_WORLD, inputpath, MPI_MODE_RDONLY,info, &gfile);
    MPI_File_read_at(gfile, 0, &n, 1, MPI_INT, &status);
    MPI_File_read_at(gfile, 4, &m, 1, MPI_INT, &status);
    MPI_File_get_size(gfile, &gfsize);
    
    offsets.resize(n+1);
    parent.resize(n);
    
    MPI_File_open(MPI_COMM_SELF,headerpath, MPI_MODE_RDONLY, MPI_INFO_NULL, &hfile);
    MPI_File_read(hfile, offsets.data(), n, MPI_INT, &status);
    MPI_File_close(&hfile);
    
    offsets[n]=gfsize;
    for(int i=0;i<n;i++){
        int pp = (offsets[i+1]-offsets[i]-8)/4;
        deg_ordered_vertices.push_back({pp,i});
    }
    
    sort(deg_ordered_vertices.begin(),deg_ordered_vertices.end());
    id2owner.resize(n);
    for(int i=0;i<deg_ordered_vertices.size();i++)
    {
        id2owner[deg_ordered_vertices[i].S]=i%numproc;
    }
    adjacency_list.resize(n);

    #pragma omp for
    for(int i = myrank;i<n;i+=numproc)
    {
        my_vertices.push_back(deg_ordered_vertices[i].second);
        int temp_degree = deg_ordered_vertices[i].first;
        vector<int> tempnbr(temp_degree);
        MPI_File_read_at(gfile,offsets[deg_ordered_vertices[i].second]+8,tempnbr.data(),temp_degree,MPI_INT,&status);
        for(int j=0;j<temp_degree;j++)
        {
            if(deg_ordered_vertices[i].second<tempnbr[j])
            {
                adjacency_list[deg_ordered_vertices[i].second].push_back(tempnbr[j]);
            }
        }
    }
}

inline int get_owner(int v)
{
    return id2owner[v];
}

inline bool compare(Packet p1,Packet p2)
{
    return p1.owner<p2.owner;
}

bool mintruss(int k)
{
    
    stack<pair<int,int>> del;
    bool local_empty=false,global_empty=false;
    for(auto it:edge_triangles)
    {
        if(it.S.size()<k-2)
            del.push(it.F);
    }
    
    int v1,v2,owner;
    while(!global_empty)
    {
        vector<int> scount(numproc,0);
        vector<int> rcount(numproc,0);
        vector<Packet> sData;
        while(!del.empty())
        {
            pair<int,int> edge=del.top();
            del.pop();
            v1=edge.F;
            v2=edge.S;
            if(edge_triangles.find(make_pair(min(v1,v2),max(v1,v2)))==edge_triangles.end())
            {
                continue;
            }
            for(auto &it:edge_triangles[make_pair(min(v1,v2),max(v1,v2))])
            {
                owner=get_owner(min(v1,it));
                sData.push_back({owner,min(v1,it),max(v1,it),v2});
                scount[owner]++;
                owner=get_owner(min(v2,it));
                sData.push_back({owner,min(v2,it),max(v2,it),v1});
                scount[owner]++;
            }
            my_edges.erase(make_pair(min(v1,v2),max(v1,v2)));
            edge_triangles.erase({min(v1,v2),max(v1,v2)});
            adjacency_list[v1].erase(find(adjacency_list[v1].begin(),adjacency_list[v1].end(),v2));
        }
        sort(sData.begin(),sData.end(),compare);
        MPI_Alltoall(scount.data(),1,MPI_INT,rcount.data(),1,MPI_INT,MPI_COMM_WORLD);
        
        vector<int> sdisp(numproc,0);
        vector<int> rdisp(numproc,0);
        sdisp[0]=0;rdisp[0]=0;
        for(int i=1;i<numproc;i++)
        {
            sdisp[i]=sdisp[i-1]+scount[i-1];
            rdisp[i]=rdisp[i-1]+rcount[i-1];
        }
        
        int trecv=rdisp[numproc-1]+rcount[numproc-1];
        vector<Packet> rData;
        rData.resize(trecv+100);
        
        MPI_Alltoallv(sData.data(),scount.data(),sdisp.data(),packet,rData.data(),rcount.data(),rdisp.data(),packet,MPI_COMM_WORLD);
        
        for(int i=0;i<rData.size();i++)
        {
            int x=rData[i].x;
            int y=rData[i].y;
            int z=rData[i].z;
            if(edge_triangles.find(make_pair(x,y))==edge_triangles.end())
            {
                continue;
            }
            if(edge_triangles[{x,y}].find(z)!=edge_triangles[{x,y}].end())
            {   
                edge_triangles[{x,y}].erase(edge_triangles[{x,y}].find(z));    
            } 
            if(edge_triangles[{x,y}].size()<k-2)
            {
                del.push({x,y});
            }
        }
        local_empty=del.empty();
        MPI_Allreduce(&local_empty,&global_empty,1,MPI_C_BOOL,MPI_LAND,MPI_COMM_WORLD);
    }
    int x=0;
    for(auto it:edge_triangles)
    {
        x+=it.S.size();
    }
    return x!=0 && my_edges.size()!=0;
}

void createDatatype()
{
    MPI_Datatype types[4]={MPI_INT,MPI_INT,MPI_INT,MPI_INT};
    MPI_Datatype types1[2]={MPI_INT,MPI_INT};

    int blockcount[4]={1,1,1,1};
    int blockcount1[2]={1,1};
    
    MPI_Aint disp[4];
    MPI_Aint disp1[2];
    
    pair<int,int> e=make_pair(0,0);
    Packet arr[2];
    Packet p;
    p.owner=0;p.x=0;p.y=0;p.y=0;
    arr[0]=p;
    arr[1]=p;
    
    MPI_Get_address(arr,disp);
    MPI_Get_address(&arr[0].x,disp+1);
    MPI_Get_address(&arr[0].y,disp+2);
    MPI_Get_address(&arr[0].z,disp+3);
    
    MPI_Get_address(&e.F,disp1);
    MPI_Get_address(&e.S,disp1+1);
    
    disp1[1]-=disp1[0];
    disp1[0]-=disp1[0];
    for	(int i=3;i>=0;i--)	
        disp[i]	-=	disp[0];
    
    MPI_Type_create_struct(4,blockcount,disp,types,&packet);
    MPI_Type_commit(&packet);

    MPI_Type_create_struct(2,blockcount1,disp1,types1,&pair_type);
    MPI_Type_commit(&pair_type);   
}

int root(int a)
{
	if (a == parent[a]) {
		return a;
	}
	return parent[a] = root(parent[a]);
}

void connect(int a, int b)
{
	a = root(a);
	b = root(b);
	if(a != b)
    {
		parent[b] = a;
	}
}

int connectedComponents(int n)
{
    connected.clear();	
	for(int i = 0; i < n; i++) 
    {
        connected[root(parent[i])].insert(i);
	}
    int count=0;
    for(auto it:connected)
    {
        if(it.S.size()>1)
            count++;
    }
    return count;
}

void makedsu()
{
	for(int i = 0; i <= n; i++) 
    {
		parent[i] = i;
	}
    int rcounts[numproc];
    int scount=0;
    for(auto it:edge_triangles)
    {
        if(it.S.size()>0)
            scount++;
    }
    MPI_Gather(&scount,1,MPI_INT,rcounts,1,MPI_INT,0,MPI_COMM_WORLD);
    vector<pair<int,int>> sedges,redges;
    sedges.clear();
    for(auto it:edge_triangles)
    {
        if(it.S.size()>0)
            sedges.push_back(it.F);
    }
    if(myrank==0)
    {
        for(int i=0;i<sedges.size();i++)
        {
            connect(sedges[i].F,sedges[i].S);
        }

        #pragma omp for
        for(int i=1;i<numproc;i++)
        {
            redges.resize(rcounts[i]);
            MPI_Recv(redges.data(),rcounts[i],pair_type,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            for(int j=0;j<redges.size();j++)
            {
                connect(redges[j].F,redges[j].S);
            }
            redges.clear();
        }
    }
    else
    {
        MPI_Send(sedges.data(),scount,pair_type,0,0,MPI_COMM_WORLD);
        return;
    }
    int size=connectedComponents(n);
    if(taskid==1)
    {
        outfile<<size<<endl;
        for(auto it:connected)
        {
            if(it.S.size()>1)
            {
                for(auto &it1:it.S)
                {
                    outfile<<it1<<" ";
                }
                outfile<<endl;
            }
        }
    }
}

void task1_runner()
{
    if(myrank==0)
        outfile.open(outputpath);
    inputData();
    initialise();
    createDatatype();
    for(int kvalue=startk;kvalue<=endk;kvalue++)
    {
        bool lcheck,gcheck;
        lcheck=mintruss(kvalue+2);
        MPI_Allreduce(&lcheck,&gcheck,1,MPI_C_BOOL,MPI_LOR,MPI_COMM_WORLD);
        string s;
        if(verbose==1)
            s="\n";
        else
            s=" ";    
        if(myrank==0)
        {
            outfile<<gcheck<<s;
        }
        if(verbose==1 && gcheck)
        {
            makedsu();
        }        
    }
    if(myrank==0)
        outfile.close();
}

void calculate_ego()
{
    int gcount=1;
    map<int,int> groups;
    for(auto it:connected)
    {
        if(it.S.size()>1)
        {
            for(auto &it1:it.S)
            {
                groups[it1]=gcount;
            }
            gcount++;
        }
    }
    map<int,set<int>> has_group;
    for(int i=0;i<n;i++)
    {
        int degree=(offsets[i+1]-offsets[i]-8)/4;
        vector<int> nbr(degree);
        MPI_File_read_at(gfile , offsets[i]+8 , nbr.data() ,degree , MPI_INT ,&status);
        for(int j=0;j<degree;j++)
        {
            if(groups[nbr[j]]>0)
            {
                has_group[i].insert(groups[nbr[j]]);
            }
        }
    }
    int count =0;
    vector<int> ans1;
    ans1.clear();
    for(auto it:has_group)
    {
        if(it.S.size()>=p)
        {
            count++;
            ans1.push_back(it.F);
        }
    }
    ofstream outfile;
    outfile.open(outputpath);
    outfile<<ans1.size()<<endl;
    for(int i=0;i<ans1.size();i++)
    {
        outfile<<ans1[i]<<" ";
        if(verbose==1)
        {
            outfile<<endl;
            for(auto it:groups)
            {
                if(has_group[ans1[i]].count(it.S))
                outfile<<it.F<<" ";
            }
            outfile<<endl;
        }
    }
    outfile.close();
}

void task2_runner()
{
    inputData();
    initialise();
    createDatatype();
    bool lcheck,gcheck;
    lcheck=mintruss(endk+2);
    MPI_Allreduce(&lcheck,&gcheck,1,MPI_C_BOOL,MPI_LOR,MPI_COMM_WORLD);
    makedsu();
    if(myrank==0)
        calculate_ego();
}

int main(int argc,char ** argv)
{
    int provided;

    MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
    MPI_Comm_size(MPI_COMM_WORLD,&numproc);
    MPI_Comm_rank( MPI_COMM_WORLD ,&myrank);
    
    taskid=atoi(argv[1]+9);
    inputpath=argv[2]+12;
    headerpath=argv[3]+13;
    outputpath=argv[4]+13;
    verbose=atoi(argv[5]+10);
    startk=atoi(argv[6]+9);
    endk=atoi(argv[7]+7);

    if(taskid!=1)
        p=atoi(argv[8]+4);
    
    if(taskid==1)
    {
        task1_runner();
    }
    else if(taskid==2)
    {
        task2_runner();
    }

    MPI_File_close(&gfile);
    MPI_Finalize();
    return 0;
}