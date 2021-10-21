#ifndef _GAWORKER_H_
#define _GAWORKER_H_
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <set>
#include <ctime>
#include <algorithm>
#include "Tool.h"
using namespace std;

class GA_worker
{
    friend class MSA_worker;
private:
    int* len_data;
    int len_input;
    int len_x;
    int len_gap;
    int data_size;
    char** MSA_data;
    char* query_input;
    char* query_comp;
    set<GA_node2> chrom;
    vector<GA_node2> child_chrom;

    void make_worker(char** input_data,int sz,int* ldata)
    {
        MSA_data=input_data;
        data_size=sz;
        len_data=ldata;
    }

    void build_chrom()
    {
        int j=0,delta=LIM_SPACE_MAX-LIM_SPACE_MIN+1;
        int y;
        int ix,iy;
        GA_node2 xi;
        len_gap=max(len_x,len_input)-min(len_x,len_input);
        chrom.clear();
        for(int i=1;i<=ID_CHROM_NUM;i++,j=(j+1)%delta)
        {
            if(len_input<len_x)
            {
                xi.lenx=len_gap+LIM_SPACE_MIN+j;
                xi.leny=LIM_SPACE_MIN+j;
                xi.final_len=len_x+LIM_SPACE_MIN+j;
            }
            else
            {
                xi.lenx=LIM_SPACE_MIN+j;
                xi.leny=len_gap+LIM_SPACE_MIN+j;
                xi.final_len=len_input+LIM_SPACE_MIN+j;
            }
            if(xi.final_len>len_x+len_input-1)
            {
                cout<<"OVERFLOW:";
                cout<<len_x+len_input-1<<endl;
                cout<<xi.final_len<<endl;
                continue;
            }
            //cout<<3<<endl;
            for(int k=0;k<=xi.final_len+1;k++)
                xi.vis[k]=0;
            //cout<<4<<endl;
            for(int k=1;k<=xi.lenx;k++)
            {
                do
                {
                    y=1+rand()%xi.final_len;
                }while(xi.vis[y]==1);
                xi.vis[y]=1;
                xi.query[k]=y;
            }
            //sort(xi.query+1,xi.query+xi.lenx+1,cmp);
            for(int k=1;k<=xi.leny;k++)
            {
                do
                {
                    y=1+rand()%xi.final_len;
                }while(xi.vis[y]==1||xi.vis[y]==2);
                xi.vis[y]=2;
                xi.comp[k]=y;
            }
            //sort(xi.comp+1,xi.comp+xi.leny+1,cmp);
            xi.fitness=0; ix=1; iy=1;
            //cout<<"NUM: "<<i<<endl;
            //cout<<1<<endl;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]==1)
                {
                    xi.fitness+=GAP;
                    iy++;
                }
                else if(xi.vis[k]==2)
                {
                    xi.fitness+=GAP;
                    ix++;
                }
                else 
                {
                    if(query_input[ix]!=query_comp[iy])
                        xi.fitness+=MISMATCH;
                    ix++; iy++;
                }
            //cout<<2<<endl;
            /*debug ouput*/
            /*int l=1;
            for(int k=1;k<=xi.final_len;k++)
                if(vis[k]==2*i-1)
                    printf("-");
                else
                    printf("%c",query_input[l++]);
            printf("\n");
            l=1;
            //for(int k=1;k<=itr->leny;k++)
            //    vis[itr->comp[i]]=2*i;
            for(int k=1;k<=xi.final_len;k++)
            {
                if(vis[k]==2*i)
                    printf("-");
                else
                    printf("%c",query_comp[l++]);
            }
            printf("\nFitness: %d\n\n",xi.fitness);*/
            //cout<<xi.fitness<<endl;
            chrom.insert(xi);
            //cout<<5<<endl;
        }
    }

    void build_space_exc_chrom(int space_num)
    {
        int y;
        int ix,iy;
        GA_node2 xi;
        len_gap=max(len_x,len_input)-min(len_x,len_input);
        chrom.clear();
        for(int i=1;i<=CHROM_NUM;i++)
        {
            if(len_input<len_x)
            {
                xi.lenx=len_gap+LIM_SPACE_MIN+space_num;
                xi.leny=LIM_SPACE_MIN+space_num;
                xi.final_len=len_x+LIM_SPACE_MIN+space_num;
            }
            else
            {
                xi.lenx=LIM_SPACE_MIN+space_num;
                xi.leny=len_gap+LIM_SPACE_MIN+space_num;
                xi.final_len=len_input+LIM_SPACE_MIN+space_num;
            }

            if(xi.final_len>len_x+len_input-1)
            {
                cout<<"OVERFLOW:";
                cout<<len_x+len_input-1<<endl;
                cout<<xi.final_len<<endl;
                continue;
            }

            for(int k=0;k<=xi.final_len+1;k++)
                xi.vis[k]=0;

            for(int k=1;k<=xi.lenx;k++)
            {
                do
                {
                    y=1+rand()%xi.final_len;
                }while(xi.vis[y]==1);
                xi.vis[y]=1;
                xi.query[k]=y;
            }
 
            for(int k=1;k<=xi.leny;k++)
            {
                do
                {
                    y=1+rand()%xi.final_len;
                }while(xi.vis[y]==1||xi.vis[y]==2);
                xi.vis[y]=2;
                xi.comp[k]=y;
            }

            xi.fitness=0; ix=1; iy=1;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]==1)
                {
                    xi.fitness+=GAP;
                    iy++;
                }
                else if(xi.vis[k]==2)
                {
                    xi.fitness+=GAP;
                    ix++;
                }
                else 
                {
                    if(query_input[ix]!=query_comp[iy])
                        xi.fitness+=MISMATCH;
                    ix++; iy++;
                }
            chrom.insert(xi);
        }
    }
    
    void get_2_string_match(char* query_comp,set<GA_node2>::iterator& itr)
    {
        int j=1;
        for(int i=1;i<=itr->lenx;i++)
            printf("%d ",itr->query[i]);
        printf("\n");
        for(int i=1;i<=itr->leny;i++)
            printf("%d ",itr->comp[i]);
        printf("\n");
        for(int i=1;i<=itr->final_len;i++)
            printf("%d ",itr->vis[i]);
        printf("\n");
        printf("%d\n",itr->final_len);
        printf("%d %d\n",itr->lenx,itr->leny);
        for(int i=1;i<=itr->final_len;i++)
            if(itr->vis[i]==1)
                printf("-");
            else
                printf("%c",query_input[j++]);
        printf("\n");
        j=1;
        for(int i=1;i<=itr->final_len;i++)
            if(itr->vis[i]==2)
                printf("-");
            else
                printf("%c",query_comp[j++]);
        printf("\nFitness: %d\n\n",itr->fitness);
    }

    void build_debug()
    {
        set<GA_node2>::iterator itr;
        int ti=0;
        cout<<"Count:"<<chrom.size()<<endl;
        for(itr=chrom.begin();itr!=chrom.end();itr++)
        {
            printf("Num: %d\n",++ti);
            get_2_string_match(query_comp,itr);
        }
    }

    void get_2_mut_exc(set<GA_node2>::iterator& itr)
    {
        GA_node2 xi=*itr;
        int idx,idy;
        int chx,chy;
        int tox,toy;
        if(xi.lenx==0||xi.leny==0)
            return ;
        idx=1+rand()%xi.lenx;
        chx=xi.query[idx];
        idy=1+rand()%xi.leny;
        chy=xi.comp[idy];
        xi.vis[chx]=2;
        xi.vis[chy]=1;
        xi.query[idx]=chy;
        xi.comp[idy]=chx;
        int ix,iy;
        xi.fitness=0; ix=1; iy=1;
        for(int k=1;k<=xi.final_len;k++)
            if(xi.vis[k]==1)
            {
                xi.fitness+=GAP;
                iy++;
            }
            else if(xi.vis[k]==2)
            {
                xi.fitness+=GAP;
                ix++;
            }
            else 
            {
                if(query_input[ix]!=query_comp[iy])
                    xi.fitness+=MISMATCH;
                ix++; iy++;
            }
        child_chrom.push_back(xi);
    }

    void get_2_mut_eql(set<GA_node2>::iterator& itr)
    {
        GA_node2 xi=*itr;
        int idx,idy;
        int chx,chy;
        int tox,toy;
        int pi=rand()%100;
        if(xi.lenx==0&&xi.leny==0)
            return ;
        if(xi.leny==0||pi<=50)
        {
            idx=1+rand()%xi.lenx;
            chx=xi.query[idx];
            do
            {
                tox=1+rand()%xi.final_len;
            }while(xi.vis[tox]==1||xi.vis[tox]==2);
            xi.vis[chx]=0;
            xi.vis[tox]=1;
            xi.query[idx]=tox;
        }
        else
        {
            idy=1+rand()%xi.leny;
            chy=xi.comp[idy];
            do
            {
                toy=1+rand()%xi.final_len;
            }while(xi.vis[toy]==1||xi.vis[toy]==2);
            xi.vis[chy]=0;
            xi.vis[toy]=2;
            xi.comp[idy]=toy;
        }
        int ix,iy;
        xi.fitness=0; ix=1; iy=1;
        for(int k=1;k<=xi.final_len;k++)
            if(xi.vis[k]==1)
            {
                xi.fitness+=GAP;
                iy++;
            }
            else if(xi.vis[k]==2)
            {
                xi.fitness+=GAP;
                ix++;
            }
            else 
            {
                if(query_input[ix]!=query_comp[iy])
                    xi.fitness+=MISMATCH;
                ix++; iy++;
            }
        child_chrom.push_back(xi);
    }

    void select_good()
    {
        int cou=0;
        set<GA_node2>::iterator itr;
        for(itr=chrom.begin();itr!=chrom.end();itr++)
        {
            cou++;
            if(cou>ID_CHROM_NUM)
                break;
        }
        chrom.erase(itr,chrom.end());
        //cout<<chrom.size()<<endl;
    }

    void select_good_space_exc()
    {
        int cou=0;
        set<GA_node2>::iterator itr;
        for(itr=chrom.begin();itr!=chrom.end();itr++)
        {
            cou++;
            if(cou>CHROM_NUM)
                break;
        }
        chrom.erase(itr,chrom.end());
        //cout<<chrom.size()<<endl;
    }

    GA_node2 GA_2_iteration()
    {
        build_chrom();
        //build_debug();
        //system("pause");
        return *chrom.begin();
    }

    GA_node2 GA_2_step()
    {
        GA_node2 re,xi;
        int ans=inf;
        int delta=LIM_SPACE_MAX-LIM_SPACE_MIN+1;
        for(int ti=1;ti<=EPOCH;ti++)
            for(int j=0;j<=delta;j++)
            {
                build_space_exc_chrom(j);
                set<GA_node2>::iterator itr;
                for(int i=1;i<=LIM_GEN;i++)
                {
                    child_chrom.clear();
                    //cout<<"COUNT: "<<chrom.size()<<endl;
                    //build_debug();
                    //getchar();
                    //getchar();
                    for(itr=chrom.begin();itr!=chrom.end();itr++)
                    {
                        get_2_mut_eql(itr);
                        get_2_mut_exc(itr);
                    }
                    chrom.insert(child_chrom.begin(),child_chrom.end());
                    select_good_space_exc();
                }
                xi=*chrom.begin();
                if(xi.fitness<ans)
                {
                    ans=xi.fitness;
                    re=xi;
                }
                /*GA_2_visual_step(re);
                cout<<ans<<endl;
                getchar();
                getchar();*/
            }
        return re;
    }

    void GA_2_visual_step(GA_node2& xi)
    {
        int j=1;
        for(int i=1;i<=xi.final_len;i++)
            if(xi.vis[i]==1)
                printf("-");
            else
                printf("%c",query_input[j++]);
        printf("\n");
        j=1;
        for(int i=1;i<=xi.final_len;i++)
            if(xi.vis[i]==2)
                printf("-");
            else
                printf("%c",query_comp[j++]);
        printf("\n\n");
    }

public:
    GA_worker()
    {
        data_size=0;
    }

    ~GA_worker()
    {
        
    }

    void GA_2(char* input)
    {
        clock_t st,ed;
        st=clock();
        query_input=input;
        int ans=inf;
        GA_node2 re,ans_output;
        int ans_id=1;
        len_input=strlen(query_input+1);
        //pre-ID-GA
        for(int ti=1;ti<=EPOCH;ti++)
            for(int i=1;i<=data_size;i++)
            {
                query_comp=MSA_data[i];
                len_x=len_data[i];
                re=GA_2_iteration();
                //printf("NUM %d ANS: %d\n",i,re);
                if(re.fitness<ans)
                {
                    ans_id=i;
                    ans=re.fitness;
                    ans_output=re;
                }
            }
            
        //GA
        query_comp=MSA_data[ans_id];
        len_x=len_data[ans_id];
        re=GA_2_step();
        if(ans>re.fitness)
        {
            ans=re.fitness;
            ans_output=re;
        }
        printf("Answer: %d\nComparsion ID: %d\n",ans,ans_id);
        ed=clock();
        printf("Time cost: %lfs\n",(double)(ed-st)/CLOCKS_PER_SEC);
        GA_2_visual_step(ans_output);
    }
};

#endif
