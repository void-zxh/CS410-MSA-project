#ifndef _MSAWORKER_H_
#define _MSAWORKER_H_
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <ctime>
#include "DP_worker.h"
#include "ASTAR_worker.h"
#include "IDGA_worker.h"
using namespace std;

class MSA_worker
{
private:
    int data_size;
    int len_data[DATA_BASE_SIZE+5];
    char* MSA_data[DATA_BASE_SIZE+5];
    int len_input;
    char query_input[MAX_LEN+5];

    int work_state;
    DP_worker dp;
    ASTAR_worker A_Star;
    IDGA_worker ga;

    void dimension_2(FILE* fp)
    {
        int ID=0;
        printf("Please choose the workstate (0: DP,1: A_Star, 2: GA): \n");
        scanf("%d",&work_state);
        while(~fscanf(fp,"%s",query_input+1))
        {
            printf("MSQ query 2_dimension ID: %d\nInput: %s\n",++ID,query_input+1);
            switch(work_state)
            {
                case 0:
                    printf("Dynamic programing: \n"); 
                    dp.DP_2(query_input);break;
                case 1:
                    printf("A_Star searching: \n"); 
                    A_Star.ASTAR_2(query_input);break;
                case 2:
                    srand(time(NULL));
                    printf("ID_GA algortihm: \n"); 
                    ga.IDGA_2(query_input);
                default:break;
            }
        }
    }

    void dimension_3(FILE* fp)
    {
        int ID=0;
        printf("Please choose the workstate (0: DP,1: A_Star, 2: GA): \n");
        scanf("%d",&work_state);
        while(~fscanf(fp,"%s",query_input+1))
        {
            printf("MSQ query 3_dimension ID: %d\nInput: %s\n",++ID,query_input+1);
            switch(work_state)
            {
                case 0:
                    printf("Dynamic programing: \n"); 
                    dp.DP_3(query_input);break;
                case 1:
                    printf("A_Star searching: \n"); 
                    A_Star.ASTAR_3(query_input,ID);break;
                case 2:
                    srand(time(NULL));
                    printf("GA algorithm: \n"); 
                    ga.IDGA_3(query_input);
                default:break;
            }
        }
    }

public:
    MSA_worker()
    {
        data_size=0;
        for(int i=0;i<DATA_BASE_SIZE+5;i++)
            MSA_data[i]=new char[MAX_LEN+5];
    }

    ~MSA_worker()
    {
        for(int i=0;i<DATA_BASE_SIZE+5;i++)
            free(MSA_data[i]);
    }

    void make_worker()
    {
        printf("Making worker ...\n");
        dp.make_worker(MSA_data,data_size,len_data); 
        A_Star.make_worker(MSA_data,data_size,len_data);     
        ga.make_worker(MSA_data,data_size,len_data);
    }

    void data_load(char* data_path)
    {
        printf("Data loading ...\n");
        FILE* fp=fopen(data_path,"r");
        while(~fscanf(fp,"%s",MSA_data[++data_size]+1))
        {
            len_data[data_size]=strlen(MSA_data[data_size]+1);
            printf("%s\n",MSA_data[data_size]+1);
        }
        fclose(fp);
        data_size--;
    }

    void query(char* query_path)
    {
        FILE* fp=fopen(query_path,"r");
        int query_d;
        fscanf(fp,"%d",&query_d);
        if(query_d==2)
            dimension_2(fp);
        else
            dimension_3(fp);
        fclose(fp);
    }
};

#endif
