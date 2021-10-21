#ifndef _TOOL_H_
#define _TOOL_H_
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
using namespace std;
#define DATA_BASE_SIZE 100
#define MAX_LEN 150

const int inf=0x3f3f3f3f;
const int GAP=2;
const int MISMATCH=3;
const int EPOCH=400;
const int LIM_GEN=100;
const int LIM_SPACE_MAX=20;//20
const int LIM_SPACE_MIN=8;//8
const int ID_CHROM_NUM=65;
const int CHROM_BUM=130;
//const float P_CROSS=0.3;
//const float P_MUT=0.3;
//const float P_DBMUT=0.3;

bool cmp(int x,int y)
{
    return x<y;
}

struct Node2
{
    int x,y;
    Node2(int i=0,int j=0){x=i;y=j;}
};

struct Node3
{
    int x,y,z;
    Node3(int i=0,int j=0,int k=0){x=i;y=j;z=k;}
};

struct ASTAR_node2
{
    int x,y,g,h;
    Node2 pre;
    ASTAR_node2(int i=0,int j=0,int gi=inf,int hi=inf){x=i;y=j;g=gi;h=hi;pre=Node2(-1,-1);}
    bool operator < (ASTAR_node2 a)const
	{
		return a.g+a.h<g+h;
	}
};

struct ASTAR_node3
{
    int x,y,z,g,h;
    Node3 pre;
    ASTAR_node3(int i=0,int j=0,int k=0,int gi=inf,int hi=inf){x=i;y=j;z=k;g=gi;h=hi;pre=Node3(-1,-1,-1);}
    bool operator < (ASTAR_node3 a)const
	{
		return a.g+a.h<g+h;
	}
};

struct GA_node2
{
    int lenx,leny,final_len;
    int query[MAX_LEN+5];
    int comp[MAX_LEN+5];
    int vis[MAX_LEN*2+10];
    int fitness;
    GA_node2(int lx=0,int ly=0){lenx=lx;leny=ly;fitness=inf;}
    bool operator < (GA_node2 a)const
	{
		if(fitness==a.fitness)
        {
            if(lenx==a.lenx)
            {
                if(leny==a.leny)
                {
                    for(int i=1;i<=lenx;i++)
                        if(query[i]<a.query[i])
                            return 1;
                        else if(query[i]>a.query[i])
                            return 0;
                    for(int i=1;i<=leny;i++)
                        if(comp[i]<a.comp[i])
                            return 1;
                        else if(comp[i]>a.query[i])
                            return 0;
                    return 0;
                }
                return leny<a.leny;
            }
            else
                return lenx<a.lenx;
        }
        else
            return fitness<a.fitness;
	}
};

void visual_generate_2(Node2** pre_2,char* query_input,int len_input,char* query_comp,int len)
{
    vector<Node2> visual_path;
    Node2 cur_state=Node2(len_input,len);
    while(1)
    {
        if(cur_state.x==0&&cur_state.y==0)
            break;
        visual_path.push_back(cur_state);
        cur_state=pre_2[cur_state.x][cur_state.y];
    }

    int len_path=visual_path.size();
    int cur=0;
    printf("Output Comparsion:\n");
    for(int i=len_path-1;i>=0;i--)
    {
        cur_state=visual_path[i];
        if(cur_state.x==cur+1)
            printf("%c",query_input[cur_state.x]);
        else
            printf("-");
        cur=cur_state.x;
    }
    printf("\n");
    cur=0;
    for(int i=len_path-1;i>=0;i--)
    {
        cur_state=visual_path[i];
        if(cur_state.y==cur+1)
            printf("%c",query_comp[cur_state.y]);
        else
            printf("-");
        cur=cur_state.y;
    }
    printf("\n\n");
}

void visual_generate_3(Node3*** pre_3,char* query_input,int len_input,char* query_x,int lenx,char* query_y,int leny)
{
    vector<Node3> visual_path;
    Node3 cur_state=Node3(len_input,lenx,leny);
    while(1)
    {
        if(cur_state.x==0&&cur_state.y==0&&cur_state.z==0)
            break;
        visual_path.push_back(cur_state);
        cur_state=pre_3[cur_state.x][cur_state.y][cur_state.z];
    }

    int len_path=visual_path.size();
    int cur=0;
    printf("Output Comparsion:\n");
    for(int i=len_path-1;i>=0;i--)
    {
        cur_state=visual_path[i];
        if(cur_state.x==cur+1)
            printf("%c",query_input[cur_state.x]);
        else
            printf("-");
        cur=cur_state.x;
    }
    printf("\n");
    cur=0;
    for(int i=len_path-1;i>=0;i--)
    {
        cur_state=visual_path[i];
        if(cur_state.y==cur+1)
            printf("%c",query_x[cur_state.y]);
        else
            printf("-");
        cur=cur_state.y;
    }
    printf("\n");
    cur=0;
    for(int i=len_path-1;i>=0;i--)
    {
        cur_state=visual_path[i];
        if(cur_state.z==cur+1)
            printf("%c",query_y[cur_state.z]);
        else
            printf("-");
        cur=cur_state.z;
    }
    printf("\n\n");
}

#endif