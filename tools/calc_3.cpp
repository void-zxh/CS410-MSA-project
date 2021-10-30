#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <cstring>
using namespace std; 

const int GAP=2;
const int MISMATCH=3;

char si[100005];
char si2[100005];
char si3[100005];

int main()
{
	int i,len;
	scanf("%s",si+1);
	scanf("%s",si2+1);
    scanf("%s",si3+1);
	len=strlen(si+1);
	int ans=0;
	for(i=1;i<=len;i++)
    {
        if((si[i]=='-'&&si2[i]=='-')||(si2[i]=='-'&&si3[i]=='-')||(si[i]=='-'&&si3[i]=='-'))
            ans+=2*GAP;
        else if(si[i]=='-')
        {
            if(si2[i]!=si3[i])
                ans+=MISMATCH+2*GAP;
            else
                ans+=2*GAP;
        }
        else if(si2[i]=='-')
        {
            if(si[i]!=si3[i])
                ans+=MISMATCH+2*GAP;
            else
                ans+=2*GAP;
        }
        else if(si3[i]=='-')
        {
            if(si[i]!=si2[i])
                ans+=MISMATCH+2*GAP;
            else
                ans+=2*GAP;
        }
        else
        {
            if(si[i]==si2[i]&&si[i]==si3[i])
                continue;
            else if(si[i]==si2[i]||si2[i]==si3[i]||si[i]==si3[i])
                ans+=2*MISMATCH;
            else
                ans+=3*MISMATCH;
        }
        //cout<<ans<<endl;
    }
	printf("%d\n",ans);	
	return 0;
}