#ifndef _IDGAWORKER_H_
#define _IDGAWORKER_H_
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <set>
#include <ctime>
#include <algorithm>
#include "Tool.h"
using namespace std;

class IDGA_worker
{
    friend class MSA_worker;
private:
    int* len_data;
    int len_input;
    int len_x;
    int len_y;//3-d
    int len_gap;
    int data_size;
    char** MSA_data;
    char* query_input;
    char* query_x;
    char* query_y;//3-d
    //2-d
    set<ID_GA_node2> chrom_2;
    vector<ID_GA_node2> child_chrom_2;
    //3-d
    set<ID_GA_node3> chrom_3;
    vector<ID_GA_node3> child_chrom_3;

    void make_worker(char** input_data,int sz,int* ldata)
    {
        MSA_data=input_data;
        data_size=sz;
        len_data=ldata;
    }

    void build_2_chrom()
    {
        int j=0,delta=LIM_SPACE_MAX-LIM_SPACE_MIN+1;
        int y;
        ID_GA_node2 xi;
        len_gap=max(len_x,len_input)-min(len_x,len_input);
        chrom_2.clear();
        for(int i=1;i<=ID_CHROM_NUM;i++,j=(j+1)%delta)
        {
            xi.final_len=max(len_x,len_input)+LIM_SPACE_MIN+j;
            xi.lenx=xi.final_len-len_input;
            xi.leny=xi.final_len-len_x;
            /*if(len_input<len_x)
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
            }*/
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
            
            /*xi.fitness=0; ix=1; iy=1;
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
                    if(query_input[ix]!=query_x[iy])
                        xi.fitness+=MISMATCH;
                    ix++; iy++;
                }*/
            fitness_update_2(xi,query_input,query_x);
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
                    printf("%c",query_x[l++]);
            }
            printf("\nFitness: %d\n\n",xi.fitness);*/
            //cout<<xi.fitness<<endl;
            chrom_2.insert(xi);
            //cout<<5<<endl;
        }
    }

    void build_2_space_exc_chrom(int space_num)
    {
        int y;
        ID_GA_node2 xi;
        len_gap=max(len_x,len_input)-min(len_x,len_input);
        chrom_2.clear();
        for(int i=1;i<=CHROM_NUM;i++)
        {
            xi.final_len=max(len_x,len_input)+LIM_SPACE_MIN+space_num;
            xi.lenx=xi.final_len-len_input;
            xi.leny=xi.final_len-len_x;

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

            /*xi.fitness=0; ix=1; iy=1;
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
                    if(query_input[ix]!=query_x[iy])
                        xi.fitness+=MISMATCH;
                    ix++; iy++;
                }*/
            fitness_update_2(xi,query_input,query_x);
            chrom_2.insert(xi);
        }
    }
    
    void build_3_chrom()
    {
        int j=0,delta=LIM_SPACE_MAX-LIM_SPACE_MIN+1;
        int y;
        int ix,iy,iz;
        ID_GA_node3 xi;
        len_gap=max(max(len_x,len_y),len_input)-min(max(len_x,len_y),len_input);
        chrom_3.clear();
        for(int i=1;i<=ID_CHROM_NUM_3;i++,j=(j+1)%delta)
        {
            xi.final_len=max(max(len_x,len_y),len_input)+LIM_SPACE_MIN+j;
            xi.lenx=xi.final_len-len_input;
            xi.leny=xi.final_len-len_x;
            xi.lenz=xi.final_len-len_y;
            if(xi.final_len>len_x+len_input+len_y-2)
            {
                cout<<"OVERFLOW:";
                cout<<len_x+len_input+len_y-2<<endl;
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
                }while((xi.vis[y]&2)==2);
                xi.vis[y]|=2;
                xi.compy[k]=y;
            }

            for(int k=1;k<=xi.lenz;k++)
            {
                do
                {
                    y=1+rand()%xi.final_len;
                }while(xi.vis[y]==3||((xi.vis[y]&4)==4));
                xi.vis[y]|=4;
                xi.compz[k]=y;
            }

            xi.fitness=0; ix=1; iy=1; iz=1;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]==1)
                {
                    xi.fitness+=2*GAP;
                    if(query_x[iy]!=query_y[iz])
                        xi.fitness+=MISMATCH;
                    iy++; iz++;
                }
                else if(xi.vis[k]==2)
                {
                    xi.fitness+=2*GAP;
                    if(query_input[ix]!=query_y[iz])
                        xi.fitness+=MISMATCH;
                    ix++; iz++;
                }
                else if(xi.vis[k]==4)
                {
                    xi.fitness+=2*GAP;
                    if(query_input[ix]!=query_x[iy])
                        xi.fitness+=MISMATCH;
                    ix++; iy++;
                }
                else if(xi.vis[k]==3)
                {
                    xi.fitness+=2*GAP;
                    iz++;
                }
                else if(xi.vis[k]==5)
                {
                    xi.fitness+=2*GAP;
                    iy++;
                }
                else if(xi.vis[k]==6)
                {
                    xi.fitness+=2*GAP;
                    ix++;
                }
                else
                {
                    if(query_input[ix]==query_x[iy]&&query_input[ix]==query_y[iz])
                        xi.fitness=xi.fitness;
                    else if(query_input[ix]==query_x[iy]||query_input[ix]==query_y[iz]||query_x[iy]==query_y[iz])
                        xi.fitness+=2*MISMATCH;
                    else
                        xi.fitness+=3*MISMATCH;
                    ix++; iy++; iz++;
                }
            //debug
            /*int l=1;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]&1)
                    printf("-");
                else
                    printf("%c",query_input[l++]);
            printf("\n");
            l=1;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]&2)
                    printf("-");
                else
                    printf("%c",query_x[l++]);
            printf("\n");
            l=1;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]&4)
                    printf("-");
                else
                    printf("%c",query_y[l++]);
            printf("\nFitness: %d\n",xi.fitness);
            getchar();
            getchar();*/
            chrom_3.insert(xi);
        }
    }

    void build_3_space_exc_chrom(int space_num)
    {
        int y;
        int ix,iy,iz;
        ID_GA_node3 xi;
        len_gap=max(max(len_x,len_y),len_input)-min(max(len_x,len_y),len_input);
        chrom_3.clear();
        for(int i=1;i<=ID_CHROM_NUM_3;i++)
        {
            xi.final_len=max(max(len_x,len_y),len_input)+LIM_SPACE_MIN+space_num;
            xi.lenx=xi.final_len-len_input;
            xi.leny=xi.final_len-len_x;
            xi.lenz=xi.final_len-len_y;
            if(xi.final_len>len_x+len_input+len_y-2)
            {
                cout<<"OVERFLOW:";
                cout<<len_x+len_input+len_y-2<<endl;
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
                }while((xi.vis[y]&2)==2);
                xi.vis[y]|=2;
                xi.compy[k]=y;
            }

            for(int k=1;k<=xi.lenz;k++)
            {
                do
                {
                    y=1+rand()%xi.final_len;
                }while(xi.vis[y]==3||((xi.vis[y]&4)==4));
                xi.vis[y]|=4;
                xi.compz[k]=y;
            }

            xi.fitness=0; ix=1; iy=1; iz=1;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]==1)
                {
                    xi.fitness+=2*GAP;
                    if(query_x[iy]!=query_y[iz])
                        xi.fitness+=MISMATCH;
                    iy++; iz++;
                }
                else if(xi.vis[k]==2)
                {
                    xi.fitness+=2*GAP;
                    if(query_input[ix]!=query_y[iz])
                        xi.fitness+=MISMATCH;
                    ix++; iz++;
                }
                else if(xi.vis[k]==4)
                {
                    xi.fitness+=2*GAP;
                    if(query_input[ix]!=query_x[iy])
                        xi.fitness+=MISMATCH;
                    ix++; iy++;
                }
                else if(xi.vis[k]==3)
                {
                    xi.fitness+=2*GAP;
                    iz++;
                }
                else if(xi.vis[k]==5)
                {
                    xi.fitness+=2*GAP;
                    iy++;
                }
                else if(xi.vis[k]==6)
                {
                    xi.fitness+=2*GAP;
                    ix++;
                }
                else
                {
                    if(query_input[ix]==query_x[iy]&&query_input[ix]==query_y[iz])
                        xi.fitness=xi.fitness;
                    else if(query_input[ix]==query_x[iy]||query_input[ix]==query_y[iz]||query_x[iy]==query_y[iz])
                        xi.fitness+=2*MISMATCH;
                    else
                        xi.fitness+=3*MISMATCH;
                    ix++; iy++; iz++;
                }
            //debug
            /*int l=1;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]&1)
                    printf("-");
                else
                    printf("%c",query_input[l++]);
            printf("\n");
            l=1;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]&2)
                    printf("-");
                else
                    printf("%c",query_x[l++]);
            printf("\n");
            l=1;
            for(int k=1;k<=xi.final_len;k++)
                if(xi.vis[k]&4)
                    printf("-");
                else
                    printf("%c",query_y[l++]);
            printf("\nFitness: %d\n",xi.fitness);
            getchar();
            getchar();*/
            chrom_3.insert(xi);
        }
    }


    void get_2_string_match(char* query_comp,set<ID_GA_node2>::iterator& itr)
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

    void build_2_debug()
    {
        set<ID_GA_node2>::iterator itr;
        int ti=0;
        cout<<"Count:"<<chrom_2.size()<<endl;
        for(itr=chrom_2.begin();itr!=chrom_2.end();itr++)
        {
            printf("Num: %d\n",++ti);
            get_2_string_match(query_x,itr);
        }
    }

    void get_2_mut_exc(set<ID_GA_node2>::iterator& itr)
    {
        ID_GA_node2 xi=*itr;
        int idx,idy;
        int chx,chy;
        //int tox,toy;
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
        /*int ix,iy;
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
                if(query_input[ix]!=query_x[iy])
                    xi.fitness+=MISMATCH;
                ix++; iy++;
            }*/
        fitness_update_2(xi,query_input,query_x);
        child_chrom_2.push_back(xi);
    }

    void get_3_mut_exc(set<ID_GA_node3>::iterator& itr)
    {
        ID_GA_node3 xi=*itr;
        int idx,idy,idz;
        int chx,chy,chz;
        //int tox,toy,toz;
        int pi=rand()%100;
        //cout<<1<<endl;
        if(pi<=33)
        {
            if(xi.lenx!=0&&xi.leny!=0)
            {
                idx=1+rand()%xi.lenx;
                chx=xi.query[idx];
                idy=1+rand()%xi.leny;
                chy=xi.compy[idy];
                if(((xi.vis[chy]&1)==1)||((xi.vis[chx]&2)==2))
                    return ;
                xi.vis[chx]-=1;
                xi.vis[chx]+=2;
                xi.vis[chy]-=2;
                xi.vis[chy]+=1;
                xi.query[idx]=chy;
                xi.compy[idy]=chx;
            }
        }
        else if(pi<=66)
        {
            if(xi.lenx!=0&&xi.lenz!=0)
            {
                idx=1+rand()%xi.lenx;
                chx=xi.query[idx];
                idz=1+rand()%xi.lenz;
                chz=xi.compz[idz];
                if(((xi.vis[chz]&1)==1)||((xi.vis[chx]&4)==4))
                    return ;
                xi.vis[chx]-=1;
                xi.vis[chx]+=4;
                xi.vis[chz]-=4;
                xi.vis[chz]+=1;
                xi.query[idx]=chz;
                xi.compz[idz]=chx;
            }
        }
        else
        {
            if(xi.leny!=0&&xi.lenz!=0)
            {
                idy=1+rand()%xi.leny;
                chy=xi.compy[idy];
                idz=1+rand()%xi.lenz;
                chz=xi.compz[idz];
                if(((xi.vis[chz]&2)==2)||((xi.vis[chy]&4)==4));
                    return ;
                xi.vis[chy]-=2;
                xi.vis[chy]+=4;
                xi.vis[chz]-=4;
                xi.vis[chz]+=2;
                xi.compy[idy]=chy;
                xi.compz[idz]=chx;
            }
        }
        
        int ix,iy,iz;
        xi.fitness=0; ix=1; iy=1; iz=1;
        for(int k=1;k<=xi.final_len;k++)
            if(xi.vis[k]==1)
            {
                xi.fitness+=2*GAP;
                if(query_x[iy]!=query_y[iz])
                    xi.fitness+=MISMATCH;
                iy++; iz++;
            }
            else if(xi.vis[k]==2)
            {
                xi.fitness+=2*GAP;
                if(query_input[ix]!=query_y[iz])
                    xi.fitness+=MISMATCH;
                ix++; iz++;
            }
            else if(xi.vis[k]==4)
            {
                xi.fitness+=2*GAP;
                if(query_input[ix]!=query_x[iy])
                    xi.fitness+=MISMATCH;
                ix++; iy++;
            }
            else if(xi.vis[k]==3)
            {
                xi.fitness+=2*GAP;
                iz++;
            }
            else if(xi.vis[k]==5)
            {
                xi.fitness+=2*GAP;
                iy++;
            }
            else if(xi.vis[k]==6)
            {
                xi.fitness+=2*GAP;
                ix++;
            }
            else
            {
                if(query_input[ix]==query_x[iy]&&query_input[ix]==query_y[iz])
                    xi.fitness=xi.fitness;
                else if(query_input[ix]==query_x[iy]||query_input[ix]==query_y[iz]||query_x[iy]==query_y[iz])
                    xi.fitness+=2*MISMATCH;
                else
                    xi.fitness+=3*MISMATCH;
                ix++; iy++; iz++;
            }
        /*cout<<xi.fitness<<endl;
        IDGA_3_visual_step(xi);
        getchar();
        getchar();*/
        child_chrom_3.push_back(xi);
    }

    void get_2_mut_eql(set<ID_GA_node2>::iterator& itr)
    {
        ID_GA_node2 xi=*itr;
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
        /*int ix,iy;
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
                if(query_input[ix]!=query_x[iy])
                    xi.fitness+=MISMATCH;
                ix++; iy++;
            }*/
        fitness_update_2(xi,query_input,query_x);
        child_chrom_2.push_back(xi);
    }

    void get_3_mut_eql(set<ID_GA_node3>::iterator& itr)
    {
        ID_GA_node3 xi=*itr;
        int idx,idy,idz;
        int chx,chy,chz;
        int tox,toy,toz;
        int pi=rand()%100;
        if(xi.lenx==0&&xi.leny==0&&xi.lenz==0)
            return ;
        if((xi.leny==0&&xi.lenz==0)||pi<=33)
        {
            idx=1+rand()%xi.lenx;
            chx=xi.query[idx];
            do
            {
                tox=1+rand()%xi.final_len;
            }while(((xi.vis[tox]&1)==1)||(xi.vis[tox]==6));
            xi.vis[chx]-=1;
            xi.vis[tox]|=1;
            xi.query[idx]=tox;
        }
        else if((xi.lenx==0&&xi.lenz==0)||(pi>33&&pi<=66))
        {
            idy=1+rand()%xi.leny;
            chy=xi.compy[idy];
            do
            {
                toy=1+rand()%xi.final_len;
            }while(((xi.vis[toy]&2)==2)||(xi.vis[toy]==5));
            xi.vis[chy]-=2;
            xi.vis[toy]|=2;
            xi.compy[idy]=toy;
        }
        else
        {
            idz=1+rand()%xi.lenz;
            chz=xi.compz[idz];
            do
            {
                toz=1+rand()%xi.final_len;
            }while(((xi.vis[toz]&4)==4)||(xi.vis[toz]==3));
            xi.vis[chz]-=4;
            xi.vis[toz]|=4;
            xi.compz[idz]=toz;
        }
        int ix,iy,iz;
        xi.fitness=0; ix=1; iy=1; iz=1;
        for(int k=1;k<=xi.final_len;k++)
            if(xi.vis[k]==1)
            {
                xi.fitness+=2*GAP;
                if(query_x[iy]!=query_y[iz])
                    xi.fitness+=MISMATCH;
                iy++; iz++;
            }
            else if(xi.vis[k]==2)
            {
                xi.fitness+=2*GAP;
                if(query_input[ix]!=query_y[iz])
                    xi.fitness+=MISMATCH;
                ix++; iz++;
            }
            else if(xi.vis[k]==4)
            {
                xi.fitness+=2*GAP;
                if(query_input[ix]!=query_x[iy])
                    xi.fitness+=MISMATCH;
                ix++; iy++;
            }
            else if(xi.vis[k]==3)
            {
                xi.fitness+=2*GAP;
                iz++;
            }
            else if(xi.vis[k]==5)
            {
                xi.fitness+=2*GAP;
                iy++;
            }
            else if(xi.vis[k]==6)
            {
                xi.fitness+=2*GAP;
                ix++;
            }
            else
            {
                if(query_input[ix]==query_x[iy]&&query_input[ix]==query_y[iz])
                    xi.fitness=xi.fitness;
                else if(query_input[ix]==query_x[iy]||query_input[ix]==query_y[iz]||query_x[iy]==query_y[iz])
                    xi.fitness+=2*MISMATCH;
                else
                    xi.fitness+=3*MISMATCH;
                ix++; iy++; iz++;
            }
        /*cout<<xi.fitness<<endl;
        IDGA_3_visual_step(xi);
        getchar();
        getchar();*/
        child_chrom_3.push_back(xi);
    }


    void select_2_good(const int& chrom_lim)
    {
        int cou=0;
        set<ID_GA_node2>::iterator itr;
        for(itr=chrom_2.begin();itr!=chrom_2.end();itr++)
        {
            cou++;
            if(cou>chrom_lim)
                break;
        }
        chrom_2.erase(itr,chrom_2.end());
        //cout<<chrom.size()<<endl;
    }

    void select_3_good(const int& chrom_lim)
    {
        int cou=0;
        set<ID_GA_node3>::iterator itr;
        for(itr=chrom_3.begin();itr!=chrom_3.end();itr++)
        {
            cou++;
            if(cou>chrom_lim)
                break;
        }
        chrom_3.erase(itr,chrom_3.end());
        //cout<<chrom.size()<<endl;
    }

    ID_GA_node2 IDGA_2_iteration()
    {
        build_2_chrom();
        //build_debug();
        //system("pause");
        return *chrom_2.begin();
    }

    ID_GA_node3 IDGA_3_iteration()
    {
        build_3_chrom();
        return *chrom_3.begin();
    }

    ID_GA_node2 IDGA_2_step()
    {
        ID_GA_node2 re,xi;
        int ans=inf;
        int delta=LIM_SPACE_MAX-LIM_SPACE_MIN+1;
        for(int ti=1;ti<=EPOCH;ti++)
            for(int j=0;j<=delta;j++)
            {
                build_2_space_exc_chrom(j);
                set<ID_GA_node2>::iterator itr;
                for(int i=1;i<=LIM_GEN;i++)
                {
                    child_chrom_2.clear();
                    //cout<<"COUNT: "<<chrom.size()<<endl;
                    //build_debug();
                    //getchar();
                    //getchar();
                    for(itr=chrom_2.begin();itr!=chrom_2.end();itr++)
                    {
                        get_2_mut_eql(itr);
                        get_2_mut_exc(itr);
                    }
                    chrom_2.insert(child_chrom_2.begin(),child_chrom_2.end());
                    select_2_good(CHROM_NUM);
                }
                xi=*chrom_2.begin();
                if(xi.fitness<ans)
                {
                    ans=xi.fitness;
                    re=xi;
                }
                /*IDGA_2_visual_step(re);
                cout<<ans<<endl;
                getchar();
                getchar();*/
            }
        return re;
    }

    ID_GA_node3 IDGA_3_step()
    {
        ID_GA_node3 re,xi;
        int ans=inf;
        int delta=LIM_SPACE_MAX-LIM_SPACE_MIN+1;
         for(int ti=1;ti<=EPOCH_3;ti++)
            for(int j=0;j<=delta;j++)
            {
                build_3_space_exc_chrom(j);
                set<ID_GA_node3>::iterator itr;
                for(int i=1;i<=LIM_GEN_3;i++)
                {
                    child_chrom_3.clear();
                    //cout<<"COUNT: "<<chrom.size()<<endl;
                    //build_debug();
                    //getchar();
                    //getchar();
                    for(itr=chrom_3.begin();itr!=chrom_3.end();itr++)
                    {
                        get_3_mut_eql(itr);
                        get_3_mut_exc(itr);
                    }
                    chrom_3.insert(child_chrom_3.begin(),child_chrom_3.end());
                    select_3_good(CHROM_NUM_3);
                }
                xi=*chrom_3.begin();
                if(xi.fitness<ans)
                {
                    ans=xi.fitness;
                    re=xi;
                }
                /*IDGA_3_visual_step(re);
                cout<<ans<<endl;
                getchar();
                getchar();*/
            }
        return re;
    }

    void IDGA_2_visual_step(ID_GA_node2& xi)
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
                printf("%c",query_x[j++]);
        printf("\n\n");
    }

    void IDGA_3_visual_step(ID_GA_node3& xi)
    {
        int l=1;
        for(int k=1;k<=xi.final_len;k++)
            if(xi.vis[k]&1)
                printf("-");
            else
                printf("%c",query_input[l++]);
        printf("\n");
        l=1;
        for(int k=1;k<=xi.final_len;k++)
            if(xi.vis[k]&2)
                printf("-");
            else
                printf("%c",query_x[l++]);
        printf("\n");
        l=1;
        for(int k=1;k<=xi.final_len;k++)
            if(xi.vis[k]&4)
                printf("-");
            else
                printf("%c",query_y[l++]);
        printf("\n\n");
    }

public:
    IDGA_worker()
    {
        data_size=0;
    }

    ~IDGA_worker(){}

    void IDGA_2(char* input)
    {
        clock_t st,ed;
        st=clock();
        query_input=input;
        int ans=inf;
        ID_GA_node2 re,ans_output;
        int ans_id=1;
        len_input=strlen(query_input+1);
        //pre-ID-GA
        for(int ti=1;ti<=ID_EPOCH;ti++)
            for(int i=1;i<=data_size;i++)
            {
                query_x=MSA_data[i];
                len_x=len_data[i];
                re=IDGA_2_iteration();
                //printf("NUM %d ANS: %d\n",i,re);
                if(re.fitness<ans)
                {
                    ans_id=i;
                    ans=re.fitness;
                    ans_output=re;
                }
            }
            
        //GA
        query_x=MSA_data[ans_id];
        len_x=len_data[ans_id];
        re=IDGA_2_step();
        if(ans>re.fitness)
        {
            ans=re.fitness;
            ans_output=re;
        }
        printf("Answer: %d\nComparsion ID: %d\n",ans,ans_id);
        ed=clock();
        printf("Time cost: %lfs\n",(double)(ed-st)/CLOCKS_PER_SEC);
        IDGA_2_visual_step(ans_output);
    }

    void IDGA_3(char* input)
    {
        clock_t st,ed;
        st=clock();
        query_input=input;
        int ans_idx=1,ans_idy=2;
        int ans=inf;
        ID_GA_node3 re,ans_output;
        ID_GA_ans3 ans3[3];
        len_input=strlen(query_input+1);
        for(int ti=1;ti<=ID_EPOCH_3;ti++)
            for(int i=1;i<=data_size;i++)
                for(int j=i+1;j<=data_size;j++)
                {
                    query_x=MSA_data[i];
                    len_x=len_data[i];
                    query_y=MSA_data[j];
                    len_y=len_data[j];
                    re=IDGA_3_iteration();
                    //printf("NUM %d ANS: %d\n",i,re);
                    if(re.fitness<ans3[0].ans)
                    {
                        if(ans3[0].idx!=i&&ans3[0].idy!=j)
                        {
                            if(ans3[0].ans!=inf)
                            {
                                ans3[2]=ans3[1];
                                ans3[1]=ans3[0];
                            }
                            ans3[0].ans_output=re;
                            ans3[0].ans=re.fitness;
                            ans3[0].idx=i; ans3[0].idy=j;
                        }
                        else
                        {
                            ans3[0].ans_output=re;
                            ans3[0].ans=re.fitness;
                        }
                    }
                    else if(re.fitness<ans3[1].ans)
                    {
                        if(ans3[1].idx!=i&&ans3[1].idy!=j)
                        {
                            if(ans3[1].ans!=inf)
                                ans3[2]=ans3[1];
                            ans3[1].ans_output=re;
                            ans3[1].ans=re.fitness;
                            ans3[1].idx=i; ans3[1].idy=j;
                        }
                        else
                        {
                            ans3[1].ans_output=re;
                            ans3[1].ans=re.fitness;
                        }
                    }
                    else if(re.fitness<ans3[2].ans)
                    {
                        if(ans3[2].idx!=i&&ans3[2].idy!=j)
                        {
                            ans3[2].ans_output=re;
                            ans3[2].ans=re.fitness;
                            ans3[2].idx=i; ans3[2].idy=j;
                        }
                        else
                        {
                            ans3[2].ans_output=re;
                            ans3[2].ans=re.fitness;
                        }
                    }
                }
        for(int i=0;i<3;i++)
        {
            query_x=MSA_data[ans3[i].idx];
            len_x=len_data[ans3[i].idx];
            query_y=MSA_data[ans3[i].idy];
            len_y=len_data[ans3[i].idy];
            re=IDGA_3_step();
            if(ans>re.fitness)
            {
                ans=re.fitness;
                ans_output=re;
                ans_idx=ans3[i].idx;
                ans_idy=ans3[i].idy;
            }
            query_x=MSA_data[ans_idx];
            len_x=len_data[ans_idx];
            query_y=MSA_data[ans_idy];
            len_y=len_data[ans_idy];
            /*cout<<ans<<endl;
            cout<<ans_idx<<' '<<ans_idy<<endl;
            IDGA_3_visual_step(ans_output);
            getchar(); getchar();*/
        }
        query_x=MSA_data[ans_idx];
        len_x=len_data[ans_idx];
        query_y=MSA_data[ans_idy];
        len_y=len_data[ans_idy];
        
        printf("Answer: %d\nComparsion ID: %d %d\n",ans,ans_idx,ans_idy);
        ed=clock();
        printf("Time cost: %lfs\n",(double)(ed-st)/CLOCKS_PER_SEC);
        IDGA_3_visual_step(ans_output);
    }
};

#endif