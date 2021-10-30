#ifndef _DPWORKER_H_
#define _DPWORKER_H_
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include "Tool.h"
using namespace std;

class DP_worker
{
    friend class MSA_worker;
private:
    int** dp_2;
    Node2** pre_2;
    int*** dp_3;
    Node3*** pre_3;
    int* len_data;
    int len_input;
    int data_size;
    char** MSA_data;
    char* query_input;

    void make_worker(char** input_data,int sz,int* ldata)
    {
        MSA_data=input_data;
        data_size=sz;
        len_data=ldata;
    }

    int dp_2_step(char* query_comp,int id)
    {
        int len=len_data[id];
        dp_2[0][0]=0;
        for(int i=1;i<=len;i++)
            dp_2[0][i]=dp_2[0][i-1]+GAP;
        for(int i=1;i<=len_input;i++)
            dp_2[i][0]=dp_2[i-1][0]+GAP;
        for(int i=1;i<=len_input;i++)
            for(int j=1;j<=len;j++)
            {
                dp_2[i][j]=min(dp_2[i-1][j],dp_2[i][j-1])+GAP;
                if(query_input[i]==query_comp[j])
                    dp_2[i][j]=min(dp_2[i-1][j-1],dp_2[i][j]);
                else
                    dp_2[i][j]=min(dp_2[i-1][j-1]+MISMATCH,dp_2[i][j]);
            }
        return dp_2[len_input][len];
    }

    void dp_2_visual_step(char* query_comp,int id)
    {
        int len=len_data[id];
        dp_2[0][0]=0; pre_2[0][0]=Node2(-1,-1);
        for(int i=1;i<=len;i++)
        {
            dp_2[0][i]=dp_2[0][i-1]+GAP;
            pre_2[0][i]=Node2(0,i-1);
        }
        for(int i=1;i<=len_input;i++)
        {
            dp_2[i][0]=dp_2[i-1][0]+GAP;
            pre_2[0][i]=Node2(i-1,0);
        }
        for(int i=1;i<=len_input;i++)
            for(int j=1;j<=len;j++)
            {
                //dp_2[i][j]=min(dp_2[i-1][j],dp_2[i][j-1])+GAP;
                if(dp_2[i-1][j]>dp_2[i][j-1])
                {
                    dp_2[i][j]=dp_2[i][j-1]+GAP;
                    pre_2[i][j]=Node2(i,j-1);
                }
                else
                {
                    dp_2[i][j]=dp_2[i-1][j]+GAP;
                    pre_2[i][j]=Node2(i-1,j);
                }
                if(query_input[i]==query_comp[j])
                {
                    if(dp_2[i-1][j-1]<dp_2[i][j])
                    {
                        dp_2[i][j]=dp_2[i-1][j-1];
                        pre_2[i][j]=Node2(i-1,j-1);
                    }
                }
                else
                {
                    if(dp_2[i-1][j-1]+MISMATCH<dp_2[i][j])
                    {
                        dp_2[i][j]=dp_2[i-1][j-1]+MISMATCH;
                        pre_2[i][j]=Node2(i-1,j-1);
                    }
                }
            }
        visual_generate_2(pre_2,query_input,len_input,query_comp,len);
    }

    int dp_3_step(char* query_x,int idx,char* query_y,int idy)
    {
        int lenx=len_data[idx],leny=len_data[idy];
        //printf("%d %d | %d %d\n",idx,idy,lenx,leny);
        dp_3[0][0][0]=0;
        for(int i=1;i<=lenx;i++)
            dp_3[0][i][0]=dp_3[0][i-1][0]+2*GAP;
        for(int i=1;i<=leny;i++)
            dp_3[0][0][i]=dp_3[0][0][i-1]+2*GAP;
        for(int i=1;i<=len_input;i++)
            dp_3[i][0][0]=dp_3[i-1][0][0]+2*GAP;

        for(int i=1;i<=lenx;i++)
            for(int j=1;j<=leny;j++)
            {
                dp_3[0][i][j]=min(dp_3[0][i-1][j],dp_3[0][i][j-1])+2*GAP;
                if(query_x[i]==query_y[j])
                    dp_3[0][i][j]=min(dp_3[0][i-1][j-1]+2*GAP,dp_3[0][i][j]);
                else
                    dp_3[0][i][j]=min(dp_3[0][i-1][j-1]+MISMATCH+2*GAP,dp_3[0][i][j]);
            }
        
        for(int i=1;i<=len_input;i++)
            for(int j=1;j<=lenx;j++)
            {
                dp_3[i][j][0]=min(dp_3[i-1][j][0],dp_3[i][j-1][0])+2*GAP;
                if(query_input[i]==query_x[j])
                    dp_3[i][j][0]=min(dp_3[i-1][j-1][0]+2*GAP,dp_3[i][j][0]);
                else
                    dp_3[i][j][0]=min(dp_3[i-1][j-1][0]+MISMATCH+2*GAP,dp_3[i][j][0]);
            }
        
        for(int i=1;i<=len_input;i++)
            for(int j=1;j<=leny;j++)
            {
                dp_3[i][0][j]=min(dp_3[i-1][0][j],dp_3[i][0][j-1])+2*GAP;
                if(query_input[i]==query_y[j])
                    dp_3[i][0][j]=min(dp_3[i-1][0][j-1]+2*GAP,dp_3[i][0][j]);
                else
                    dp_3[i][0][j]=min(dp_3[i-1][0][j-1]+MISMATCH+2*GAP,dp_3[i][0][j]);
            }
        
        for(int i=1;i<=len_input;i++)
            for(int j=1;j<=lenx;j++)
                for(int k=1;k<=leny;k++)
                {
                    dp_3[i][j][k]=min(dp_3[i-1][j][k],min(dp_3[i][j-1][k],dp_3[i][j][k-1]))+2*GAP;
                    
                    if(query_input[i]==query_x[j])
                        dp_3[i][j][k]=min(dp_3[i-1][j-1][k]+2*GAP,dp_3[i][j][k]);
                    else
                        dp_3[i][j][k]=min(dp_3[i-1][j-1][k]+MISMATCH+2*GAP,dp_3[i][j][k]);
                    
                    if(query_input[i]==query_y[k])
                        dp_3[i][j][k]=min(dp_3[i-1][j][k-1]+2*GAP,dp_3[i][j][k]);
                    else
                        dp_3[i][j][k]=min(dp_3[i-1][j][k-1]+MISMATCH+2*GAP,dp_3[i][j][k]);
                    
                    if(query_x[j]==query_y[k])
                        dp_3[i][j][k]=min(dp_3[i][j-1][k-1]+2*GAP,dp_3[i][j][k]);
                    else
                        dp_3[i][j][k]=min(dp_3[i][j-1][k-1]+MISMATCH+2*GAP,dp_3[i][j][k]);

                    if(query_input[i]==query_x[j]&&query_input[i]==query_y[k])
                        dp_3[i][j][k]=min(dp_3[i-1][j-1][k-1],dp_3[i][j][k]);
                    else if(query_input[i]==query_x[j]||query_input[i]==query_y[k]||query_x[j]==query_y[k])
                        dp_3[i][j][k]=min(dp_3[i-1][j-1][k-1]+2*MISMATCH,dp_3[i][j][k]);
                    else
                        dp_3[i][j][k]=min(dp_3[i-1][j-1][k-1]+3*MISMATCH,dp_3[i][j][k]);
                }
        return dp_3[len_input][lenx][leny];
    }

    void dp_3_visual_step(char* query_x,int idx,char* query_y,int idy)
    {
        int lenx=len_data[idx],leny=len_data[idy];

        dp_3[0][0][0]=0; pre_3[0][0][0]=Node3(-1,-1,-1);
        for(int i=1;i<=lenx;i++)
        {
            dp_3[0][i][0]=dp_3[0][i-1][0]+2*GAP;
            pre_3[0][i][0]=Node3(0,i-1,0);
        }
        for(int i=1;i<=leny;i++)
        {
            dp_3[0][0][i]=dp_3[0][0][i-1]+2*GAP;
            pre_3[0][0][i]=Node3(0,0,i-1);
        }
        for(int i=1;i<=len_input;i++)
        {
            dp_3[i][0][0]=dp_3[i-1][0][0]+2*GAP;
            pre_3[i][0][0]=Node3(i-1,0,0);
        }

        for(int i=1;i<=lenx;i++)
            for(int j=1;j<=leny;j++)
            {   
                //dp_3[0][i][j]=min(dp_3[0][i-1][j],dp_3[0][i][j-1])+2*GAP;
                if(dp_3[0][i-1][j]>dp_3[0][i][j-1])
                {
                    dp_3[0][i][j]=dp_3[0][i][j-1]+2*GAP;
                    pre_3[0][i][j]=Node3(0,i,j-1);
                }
                else
                {
                    dp_3[0][i][j]=dp_3[0][i-1][j]+2*GAP;
                    pre_3[0][i][j]=Node3(0,i-1,j);
                }
                if(query_x[i]==query_y[j])
                {
                    if(dp_3[0][i-1][j-1]+2*GAP<dp_3[0][i][j])
                    {
                        dp_3[0][i][j]=dp_3[0][i-1][j-1]+2*GAP;
                        pre_3[0][i][j]=Node3(0,i-1,j-1);
                    }
                }
                else
                {
                    if(dp_3[0][i-1][j-1]+MISMATCH+2*GAP<dp_3[0][i][j])
                    {
                        dp_3[0][i][j]=dp_3[0][i-1][j-1]+MISMATCH+2*GAP;
                        pre_3[0][i][j]=Node3(0,i-1,j-1);
                    }
                }
            }
        
        for(int i=1;i<=len_input;i++)
            for(int j=1;j<=lenx;j++)
            {   
                //dp_3[i][j][0]=min(dp_3[i-1][j][0],dp_3[i][j-1][0])+2*GAP;
                if(dp_3[i-1][j][0]>dp_3[i][j-1][0])
                {
                    dp_3[i][j][0]=dp_3[i][j-1][0]+2*GAP;
                    pre_3[i][j][0]=Node3(i,j-1,0);
                }
                else
                {
                    dp_3[i][j][0]=dp_3[i-1][j][0]+2*GAP;
                    pre_3[i][j][0]=Node3(i-1,j,0);
                }
                if(query_input[i]==query_x[j])
                {
                    if(dp_3[i-1][j-1][0]+2*GAP<dp_3[i][j][0])
                    {
                        dp_3[i][j][0]=dp_3[i-1][j-1][0]+2*GAP;
                        pre_3[i][j][0]=Node3(i-1,j-1,0);
                    }
                }
                else
                {
                    if(dp_3[i-1][j-1][0]+MISMATCH+2*GAP<dp_3[i][j][0])
                    {
                        dp_3[i][j][0]=dp_3[i-1][j-1][0]+MISMATCH+2*GAP;
                        pre_3[i][j][0]=Node3(i-1,j-1,0);
                    }
                }
            }
        
        for(int i=1;i<=len_input;i++)
            for(int j=1;j<=leny;j++)
            {
                //dp_3[i][0][j]=min(dp_3[i-1][0][j],dp_3[i][0][j-1])+2*GAP;
                if(dp_3[i-1][0][j]>dp_3[i][0][j-1])
                {
                    dp_3[i][0][j]=dp_3[i][0][j-1]+2*GAP;
                    pre_3[i][0][j]=Node3(i,0,j-1);
                }
                else
                {
                    dp_3[i][0][j]=dp_3[i-1][0][j]+2*GAP;
                    pre_3[i][0][j]=Node3(i-1,0,j);
                }
                if(query_input[i]==query_y[j])
                {
                    if(dp_3[i-1][0][j-1]+2*GAP<dp_3[i][0][j])
                    {
                        dp_3[i][0][j]=dp_3[i-1][0][j-1]+2*GAP;
                        pre_3[i][0][j]=Node3(i-1,0,j-1);
                    }
                }
                else
                {
                    if(dp_3[i-1][0][j-1]+MISMATCH+2*GAP<dp_3[i][0][j])
                    {
                        dp_3[i][0][j]=dp_3[i-1][0][j-1]+MISMATCH+2*GAP;
                        pre_3[i][0][j]=Node3(i-1,0,j-1);
                    }
                }
            }
        
        for(int i=1;i<=len_input;i++)
            for(int j=1;j<=lenx;j++)
                for(int k=1;k<=leny;k++)
                {
                    //dp_3[i][j][k]=min(dp_3[i-1][j][k],min(dp_3[i][j-1][k],dp_3[i][j][k-1]))+2*GAP;
                    if(dp_3[i-1][j][k]<=dp_3[i][j-1][k]&&dp_3[i-1][j][k]<=dp_3[i][j][k-1])
                    {
                        dp_3[i][j][k]=dp_3[i-1][j][k]+2*GAP;
                        pre_3[i][j][k]=Node3(i-1,j,k);
                    }
                    else if(dp_3[i][j-1][k]<=dp_3[i-1][j][k]&&dp_3[i][j-1][k]<=dp_3[i][j][k-1])
                    {
                        dp_3[i][j][k]=dp_3[i][j-1][k]+2*GAP;
                        pre_3[i][j][k]=Node3(i,j-1,k);
                    }
                    else 
                    {
                        dp_3[i][j][k]=dp_3[i][j][k-1]+2*GAP;
                        pre_3[i][j][k]=Node3(i,j,k-1);
                    }
                    
                    if(query_input[i]==query_x[j])
                    {
                        if(dp_3[i-1][j-1][k]+2*GAP<dp_3[i][j][k])
                        {
                            dp_3[i][j][k]=dp_3[i-1][j-1][k]+2*GAP;
                            pre_3[i][j][k]=Node3(i-1,j-1,k);
                        }
                    }
                    else
                    {
                        if(dp_3[i-1][j-1][k]+MISMATCH+2*GAP<dp_3[i][j][k])
                        {
                            dp_3[i][j][k]=dp_3[i-1][j-1][k]+MISMATCH+2*GAP;
                            pre_3[i][j][k]=Node3(i-1,j-1,k);
                        }
                    }

                    if(query_input[i]==query_y[k])
                    {
                        if(dp_3[i-1][j][k-1]+2*GAP<dp_3[i][j][k])
                        {
                            dp_3[i][j][k]=dp_3[i-1][j][k-1]+2*GAP;
                            pre_3[i][j][k]=Node3(i-1,j,k-1);
                        }
                    }
                    else
                    {
                        if(dp_3[i-1][j][k-1]+MISMATCH+2*GAP<dp_3[i][j][k])
                        {
                            dp_3[i][j][k]=dp_3[i-1][j][k-1]+MISMATCH+2*GAP;
                            pre_3[i][j][k]=Node3(i-1,j,k-1);
                        }
                    }

                    if(query_x[j]==query_y[k])
                    {
                        if(dp_3[i][j-1][k-1]+2*GAP<dp_3[i][j][k])
                        {
                            dp_3[i][j][k]=dp_3[i][j-1][k-1]+2*GAP;
                            pre_3[i][j][k]=Node3(i,j-1,k-1);
                        }
                    }
                    else
                    {
                        if(dp_3[i][j-1][k-1]+MISMATCH+2*GAP<dp_3[i][j][k])
                        {
                            dp_3[i][j][k]=dp_3[i][j-1][k-1]+MISMATCH+2*GAP;
                            pre_3[i][j][k]=Node3(i,j-1,k-1);
                        }
                    }

                    if(query_input[i]==query_x[j]&&query_input[i]==query_y[k])
                    {
                        if(dp_3[i-1][j-1][k-1]<dp_3[i][j][k])
                        {
                            dp_3[i][j][k]=dp_3[i-1][j-1][k-1];
                            pre_3[i][j][k]=Node3(i-1,j-1,k-1);
                        }
                    }
                    else if(query_input[i]==query_x[j]||query_input[i]==query_y[k]||query_x[j]==query_y[k])
                    {
                        if(dp_3[i-1][j-1][k-1]+2*MISMATCH<dp_3[i][j][k])
                        {
                            dp_3[i][j][k]=dp_3[i-1][j-1][k-1]+2*MISMATCH;
                            pre_3[i][j][k]=Node3(i-1,j-1,k-1);
                        }
                    }
                    else
                    {
                        if(dp_3[i-1][j-1][k-1]+3*MISMATCH<dp_3[i][j][k])
                        {
                            dp_3[i][j][k]=dp_3[i-1][j-1][k-1]+3*MISMATCH;
                            pre_3[i][j][k]=Node3(i-1,j-1,k-1);
                        }
                    }
                }
        
        visual_generate_3(pre_3,query_input,len_input,query_x,lenx,query_y,leny);
    }

public:
    DP_worker()
    {
        dp_2=new int*[MAX_LEN+5];
        for(int i=0;i<MAX_LEN+5;i++)
            dp_2[i]=new int[MAX_LEN+5];
        dp_3=new int**[MAX_LEN+5];
        for(int i=0;i<MAX_LEN+5;i++)
        {
            dp_3[i]=new int*[MAX_LEN+5];
            for(int j=0;j<MAX_LEN+5;j++)
                dp_3[i][j]=new int[MAX_LEN+5];
        }
        pre_2=new Node2*[MAX_LEN+5];
        for(int i=0;i<MAX_LEN+5;i++)
            pre_2[i]=new Node2[MAX_LEN+5];
        pre_3=new Node3**[MAX_LEN+5];
        for(int i=0;i<MAX_LEN+5;i++)
        {
            pre_3[i]=new Node3*[MAX_LEN+5];
            for(int j=0;j<MAX_LEN+5;j++)
                pre_3[i][j]=new Node3[MAX_LEN+5];
        }
    }
    ~DP_worker()
    {
        for(int i=0;i<MAX_LEN+5;i++)
        {
            free(dp_2[i]);
            free(pre_2[i]);
            for(int j=0;j<MAX_LEN+5;j++)
                free(dp_3[i][j]),free(pre_3[i][j]);
            free(dp_3[i]);
            free(pre_3[i]);
        }
        free(dp_2);
        free(dp_3);
        free(pre_2);
        free(pre_3);
    }
    void DP_2(char* input)
    {
        query_input=input;
        int ans=inf,re;
        int ans_id=1;
        len_input=strlen(query_input+1);
        clock_t st,ed;
        st=clock();
        for(int i=1;i<=data_size;i++)
        {
            re=dp_2_step(MSA_data[i],i);
            if(re<ans)
            {
                ans_id=i;
                ans=re;
            }
        }
        printf("Answer: %d\nComparsion ID: %d\n",ans,ans_id);
        ed=clock();
        printf("Time cost: %lfs\n",(double)(ed-st)/CLOCKS_PER_SEC);
        dp_2_visual_step(MSA_data[ans_id],ans_id);
    }

    void DP_3(char* input)
    {
        query_input=input;
        int ans=inf,re;
        int ans_idx=1,ans_idy=2;
        len_input=strlen(query_input+1);
        //printf("%d\n",data_size);
        clock_t st,ed;
        st=clock();
        for(int i=1;i<=data_size;i++)
            for(int j=i+1;j<=data_size;j++)
            {
                //printf("%d %d\n",i,j);
                re=dp_3_step(MSA_data[i],i,MSA_data[j],j);
                if(re<ans)
                {
                    ans=re;
                    ans_idx=i; ans_idy=j;
                }
            }
        printf("Answer: %d\nComparsion ID: %d %d\n",ans,ans_idx,ans_idy);
        ed=clock();
        printf("Time cost: %lfs\n",(double)(ed-st)/CLOCKS_PER_SEC);
        dp_3_visual_step(MSA_data[ans_idx],ans_idx,MSA_data[ans_idy],ans_idy);
    }
};

#endif
