#ifndef _ASTARWORKER_H_
#define _ASTARWORKER_H_
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <ctime>
#include <queue>
#include "Tool.h"
using namespace std;

class ASTAR_worker
{
    friend class MSA_worker;
private:
    int output_ans;
    int* len_data;
    int len_input;
    int data_size;
    int vis_2[MAX_LEN+5][MAX_LEN+5];
    int*** vis_3;
    char** MSA_data;
    char* query_input;
    Node2** pre_2;
    Node3*** pre_3;
    int dp_ix[MAX_LEN+5][MAX_LEN+5];
    int dp_iy[MAX_LEN+5][MAX_LEN+5];
    int dp_xy[MAX_LEN+5][MAX_LEN+5];

    void make_worker(char** input_data,int sz,int* ldata)
    {
        MSA_data=input_data;
        data_size=sz;
        len_data=ldata;
    }

    int get_2_H(int len,int x,int y)
    {
        int lx=len_input-x,ly=len-y;
        return (max(lx,ly)-min(lx,ly))*GAP;
    }

    int get_dp_H(int lenx,int leny,int x,int y,int z)
    {
        return dp_ix[x][y]+dp_iy[x][z]+dp_xy[y][z];
    }

    void dp_work(int dp_2[][MAX_LEN+5],char* query_i,int len_i,char* query_comp,int len)
    {
        dp_2[len_i][len]=0;
        for(int i=len;i>=1;i--)
            dp_2[len_i][i-1]=dp_2[len_i][i]+GAP;
        for(int i=len_i;i>=1;i--)
            dp_2[i-1][len]=dp_2[i][len]+GAP;
        for(int i=len_i-1;i>=0;i--)
            for(int j=len-1;j>=0;j--)
            {
                dp_2[i][j]=min(dp_2[i+1][j],dp_2[i][j+1])+GAP;
                if(query_i[i+1]==query_comp[j+1])
                    dp_2[i][j]=min(dp_2[i+1][j+1],dp_2[i][j]);
                else
                    dp_2[i][j]=min(dp_2[i+1][j+1]+MISMATCH,dp_2[i][j]);
            }
    }

    void Astar_dp_pre(char* query_x,int idx,char* query_y,int idy)
    {
        dp_work(dp_ix,query_input,len_input,query_x,len_data[idx]);
        dp_work(dp_iy,query_input,len_input,query_y,len_data[idy]);
        dp_work(dp_xy,query_x,len_data[idx],query_y,len_data[idy]);
        //printf("%d\n",dp_ix[len_input][len_data[idx]]);
        //printf("%d\n",dp_iy[len_input][len_data[idy]]);
    }

    void Astar_2_visual_step(char* query_comp,int id)
    {
        int len=len_data[id],re=inf;
        id=-1;
        priority_queue<ASTAR_node2> q;
        ASTAR_node2 st(0,0,0,0),xi,yi;
        st.h=get_2_H(len,0,0);
        q.push(st);
        while(!q.empty())
        {
            xi=q.top(); q.pop();
            if(vis_2[xi.x][xi.y]==id) continue;
            vis_2[xi.x][xi.y]=id;
            pre_2[xi.x][xi.y]=xi.pre;
            if(xi.x==len_input&&xi.y==len)
                break;

            if(xi.x<len_input)
            {
                yi=xi;
                yi.x=xi.x+1; yi.g+=GAP;
                yi.h=get_2_H(len,yi.x,yi.y);
                if(vis_2[yi.x][yi.y]!=id)
                {
                    yi.pre=Node2(xi.x,xi.y);
                    q.push(yi);
                }
            }

            if(xi.y<len)
            {
                yi=xi;
                yi.y=xi.y+1; yi.g+=GAP;
                yi.h=get_2_H(len,yi.x,yi.y);
                if(vis_2[yi.x][yi.y]!=id)
                {
                    yi.pre=Node2(xi.x,xi.y);
                    q.push(yi);
                }
            }

            if(xi.x<len_input&&xi.y<len)
            {
                yi=xi;
                yi.x=xi.x+1; yi.y=xi.y+1; 
                if(query_input[yi.x]!=query_comp[yi.y])
                    yi.g+=MISMATCH;
                yi.h=get_2_H(len,yi.x,yi.y);
                if(vis_2[yi.x][yi.y]!=id)
                {
                    yi.pre=Node2(xi.x,xi.y);
                    q.push(yi);
                }
            }
        }
        visual_generate_2(pre_2,query_input,len_input,query_comp,len);
    }

    int Astar_2_step(char* query_comp,int id)
    {
        int len=len_data[id];
        priority_queue<ASTAR_node2> q;
        ASTAR_node2 st(0,0,0,0),xi,yi;
        st.h=get_2_H(len,0,0);
        q.push(st);
        while(!q.empty())
        {
            xi=q.top(); q.pop();
            if(xi.g+xi.h>=output_ans)
                return inf;
            if(vis_2[xi.x][xi.y]==id) continue;
            vis_2[xi.x][xi.y]=id;
            if(xi.x==len_input&&xi.y==len)
                return xi.g;

            if(xi.x<len_input)
            {
                yi=xi;
                yi.x=xi.x+1; yi.g+=GAP;
                yi.h=get_2_H(len,yi.x,yi.y);
                if(vis_2[yi.x][yi.y]!=id)
                    q.push(yi);
            }

            if(xi.y<len)
            {
                yi=xi;
                yi.y=xi.y+1; yi.g+=GAP;
                yi.h=get_2_H(len,yi.x,yi.y);
                if(vis_2[yi.x][yi.y]!=id)
                    q.push(yi);
            }

            if(xi.x<len_input&&xi.y<len)
            {
                yi=xi;
                yi.x=xi.x+1; yi.y=xi.y+1; 
                if(query_input[yi.x]!=query_comp[yi.y])
                    yi.g+=MISMATCH;
                yi.h=get_2_H(len,yi.x,yi.y);
                if(vis_2[yi.x][yi.y]!=id)
                    q.push(yi);
            }
        }
        return inf;
    }

    int Astar_3_step(char* query_x,int idx,char* query_y,int idy,int ID)
    {
        Astar_dp_pre(query_x,idx,query_y,idy);
        int lenx=len_data[idx],leny=len_data[idy];
        int id=ID*MAX_LEN*MAX_LEN+idx*MAX_LEN+idy;
        priority_queue<ASTAR_node3> q;
        ASTAR_node3 st(0,0,0,0,0),xi,yi;
        st.h=get_dp_H(lenx,leny,0,0,0);
        //printf("%d\n",st.h);
        q.push(st);
        while(!q.empty())
        {
            xi=q.top(); q.pop();
            //printf("%d %d %d %d\n",xi.x,xi.y,xi.z,vis_3[xi.x][xi.y][xi.z]);
            if(xi.g+xi.h>=output_ans)
                return inf;
            if(vis_3[xi.x][xi.y][xi.z]==id) continue;
            vis_3[xi.x][xi.y][xi.z]=id;
            if(xi.x==len_input&&xi.y==lenx&&xi.z==leny)
                return xi.g;

            if(xi.x<len_input)
            {
                yi=xi;
                yi.x=xi.x+1; yi.g+=2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                    q.push(yi);
            }
            //printf("%d %d %d %d\n",xi.x,xi.y,xi.z,xi.g);
            if(xi.y<lenx)
            {
                yi=xi;
                yi.y=xi.y+1; yi.g+=2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                    q.push(yi);
            }

            if(xi.z<leny)
            {
                yi=xi;
                yi.z=xi.z+1; yi.g+=2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                    q.push(yi);
            }
            //printf("%d %d %d %d\n",xi.x,xi.y,xi.z,xi.g);
            if(xi.x<len_input&&xi.y<lenx)
            {
                yi=xi;
                yi.x=xi.x+1; yi.y=xi.y+1;
                if(query_input[yi.x]==query_x[yi.y])
                    yi.g+=2*GAP;
                else
                    yi.g+=MISMATCH+2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                    q.push(yi);
            }

            if(xi.x<len_input&&xi.z<leny)
            {
                yi=xi;
                yi.x=xi.x+1; yi.z=xi.z+1;
                if(query_input[yi.x]==query_y[yi.z])
                    yi.g+=2*GAP;
                else
                    yi.g+=MISMATCH+2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                    q.push(yi);
            }

            if(xi.y<lenx&&xi.z<leny)
            {
                yi=xi;
                yi.y=xi.y+1; yi.z=xi.z+1;
                if(query_x[yi.y]==query_y[yi.z])
                    yi.g+=2*GAP;
                else
                    yi.g+=MISMATCH+2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                    q.push(yi);
            }

            if(xi.x<len_input&&xi.y<lenx&&xi.z<leny)
            {
                yi=xi;
                yi.x=xi.x+1; yi.y=xi.y+1; yi.z=xi.z+1;
                if(query_input[yi.x]==query_x[yi.y]&&query_x[yi.y]==query_y[yi.z])
                    yi.g=yi.g;
                else if(query_input[yi.x]==query_x[yi.y]||query_x[yi.y]==query_y[yi.z]||query_input[yi.x]==query_y[yi.z])
                    yi.g+=2*MISMATCH;
                else
                    yi.g+=3*MISMATCH;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                    q.push(yi);
            }
        }
        return inf;
    }

    void Astar_3_visual_step(char* query_x,int idx,char* query_y,int idy)
    {
        Astar_dp_pre(query_x,idx,query_y,idy);
        int lenx=len_data[idx],leny=len_data[idy];
        int id=-1;
        priority_queue<ASTAR_node3> q;
        ASTAR_node3 st(0,0,0,0,0),xi,yi;
        st.h=get_dp_H(lenx,leny,0,0,0);
        //printf("%d\n",st.h);
        q.push(st);
        while(!q.empty())
        {
            xi=q.top(); q.pop();
            //printf("%d %d %d %d\n",xi.x,xi.y,xi.z,vis_3[xi.x][xi.y][xi.z]);
            if(vis_3[xi.x][xi.y][xi.z]==id) continue;
            vis_3[xi.x][xi.y][xi.z]=id;
            //printf("%d %d %d %d\n",xi.x,xi.y,xi.z,vis_3[xi.x][xi.y][xi.z]);
            //printf("%d\n",pre_3[xi.x][xi.y][xi.z].x);
            pre_3[xi.x][xi.y][xi.z]=xi.pre;
            //printf("%d %d %d %d\n",xi.x,xi.y,xi.z,vis_3[xi.x][xi.y][xi.z]);
            if(xi.x==len_input&&xi.y==lenx&&xi.z==leny)
                break;

            if(xi.x<len_input)
            {
                yi=xi;
                yi.x=xi.x+1; yi.g+=2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                {
                    yi.pre=Node3(xi.x,xi.y,xi.z);
                    q.push(yi);
                }
            }
            //printf("%d %d %d %d\n",xi.x,xi.y,xi.z,xi.g);
            if(xi.y<lenx)
            {
                yi=xi;
                yi.y=xi.y+1; yi.g+=2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                {
                    yi.pre=Node3(xi.x,xi.y,xi.z);
                    q.push(yi);
                }
            }

            if(xi.z<leny)
            {
                yi=xi;
                yi.z=xi.z+1; yi.g+=2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                {
                    yi.pre=Node3(xi.x,xi.y,xi.z);
                    q.push(yi);
                }
            }
            //printf("%d %d %d %d\n",xi.x,xi.y,xi.z,xi.g);
            if(xi.x<len_input&&xi.y<lenx)
            {
                yi=xi;
                yi.x=xi.x+1; yi.y=xi.y+1;
                if(query_input[yi.x]==query_x[yi.y])
                    yi.g+=2*GAP;
                else
                    yi.g+=MISMATCH+2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                {
                    yi.pre=Node3(xi.x,xi.y,xi.z);
                    q.push(yi);
                }
            }

            if(xi.x<len_input&&xi.z<leny)
            {
                yi=xi;
                yi.x=xi.x+1; yi.z=xi.z+1;
                if(query_input[yi.x]==query_y[yi.z])
                    yi.g+=2*GAP;
                else
                    yi.g+=MISMATCH+2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                {
                    yi.pre=Node3(xi.x,xi.y,xi.z);
                    q.push(yi);
                }
            }

            if(xi.y<lenx&&xi.z<leny)
            {
                yi=xi;
                yi.y=xi.y+1; yi.z=xi.z+1;
                if(query_x[yi.y]==query_y[yi.z])
                    yi.g+=2*GAP;
                else
                    yi.g+=MISMATCH+2*GAP;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                {
                    yi.pre=Node3(xi.x,xi.y,xi.z);
                    q.push(yi);
                }
            }

            if(xi.x<len_input&&xi.y<lenx&&xi.z<leny)
            {
                yi=xi;
                yi.x=xi.x+1; yi.y=xi.y+1; yi.z=xi.z+1;
                if(query_input[yi.x]==query_x[yi.y]&&query_x[yi.y]==query_y[yi.z])
                    yi.g=yi.g;
                else if(query_input[yi.x]==query_x[yi.y]||query_x[yi.y]==query_y[yi.z]||query_input[yi.x]==query_y[yi.z])
                    yi.g+=2*MISMATCH;
                else
                    yi.g+=3*MISMATCH;
                yi.h=get_dp_H(lenx,leny,yi.x,yi.y,yi.z);
                if(vis_3[yi.x][yi.y][yi.z]!=id)
                {
                    yi.pre=Node3(xi.x,xi.y,xi.z);
                    q.push(yi);
                }
            }
        }
        visual_generate_3(pre_3,query_input,len_input,query_x,lenx,query_y,leny);
    }

public:
    ASTAR_worker()
    {
        memset(vis_2,0,sizeof(vis_2));
        pre_2=new Node2*[MAX_LEN+5];
        for(int i=0;i<MAX_LEN+5;i++)
            pre_2[i]=new Node2[MAX_LEN+5];
        pre_3=new Node3**[MAX_LEN+5];
        vis_3=new int**[MAX_LEN+5];
        for(int i=0;i<MAX_LEN+5;i++)
        {
            pre_3[i]=new Node3*[MAX_LEN+5];
            vis_3[i]=new int*[MAX_LEN+5];
            for(int j=0;j<MAX_LEN+5;j++)
            {
                pre_3[i][j]=new Node3[MAX_LEN+5];
                vis_3[i][j]=new int[MAX_LEN+5];
            }
        }
    }

    ~ASTAR_worker()
    {
        for(int i=0;i<MAX_LEN+5;i++)
        {
            free(pre_2[i]);
            for(int j=0;j<MAX_LEN+5;j++)
                free(pre_3[i][j]),free(vis_3[i][j]);
            free(pre_3[i]);
            free(vis_3[i]);
        }
        free(pre_2);
        free(pre_3);
        free(vis_3);
    }

    void ASTAR_2(char* input)
    {
        memset(vis_2,0,sizeof(vis_2));
        query_input=input;
        output_ans=inf;
        int re;
        int ans_id=1;
        len_input=strlen(query_input+1);
        clock_t st,ed;
        st=clock();
        for(int i=1;i<=data_size;i++)
        {
            re=Astar_2_step(MSA_data[i],i);
            if(re<output_ans)
            {
                ans_id=i;
                output_ans=re;
            }
        }
        printf("Answer: %d\nComparsion ID: %d\n",output_ans,ans_id);
        ed=clock();
        printf("Time cost: %lfs\n",(double)(ed-st)/CLOCKS_PER_SEC);
        Astar_2_visual_step(MSA_data[ans_id],ans_id);
    }

    void ASTAR_3(char* input,int ID)
    {
        query_input=input;
        output_ans=inf;
        int re;
        int ans_idx=1,ans_idy=2;
        len_input=strlen(query_input+1);
        clock_t st,ed;
        st=clock();
        for(int i=1;i<=data_size;i++)
            for(int j=i+1;j<=data_size;j++)
            {
                //printf("%d %d\n",i,j);
                re=Astar_3_step(MSA_data[i],i,MSA_data[j],j,ID);
                if(re<output_ans)
                {
                    output_ans=re;
                    ans_idx=i; ans_idy=j;
                }
            }
        printf("Answer: %d\nComparsion ID: %d %d\n",output_ans,ans_idx,ans_idy);
        ed=clock();
        printf("Time cost: %lfs\n",(double)(ed-st)/CLOCKS_PER_SEC);
        Astar_3_visual_step(MSA_data[ans_idx],ans_idx,MSA_data[ans_idy],ans_idy);
    }
};

#endif
