#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <cstring>
using namespace std; 

const int GAP=2;
const int MISMATCH=3;

char si[100005];
char si2[100005];

int main()
{
	int i,len;
	scanf("%s",si+1);
	scanf("%s",si2+1);
	len=strlen(si+1);
	int cou1=0,cou2=0;
	for(i=1;i<=len;i++)
	{
    	if(si[i]=='-')
			cou1++;
		if(si2[i]=='-')
            cou2++;
    }
	printf("%d %d %d\n",cou1,cou2,len);	
	return 0;
}