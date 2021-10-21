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
	int ans=0;
	for(i=1;i<=len;i++)
		if(si[i]=='-'||si2[i]=='-')
			ans+=GAP;
		else if(si[i]!=si2[i])
			ans+=MISMATCH;
	printf("%d\n",ans);	
	return 0;
}