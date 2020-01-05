#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
typedef struct//脱靶概率
{
  int spot;//脱靶位置
  int mis;//脱靶数量
  float po;//概率
  int spotnum[9]; //脱靶位点 
}offtarget;

offtarget p[100];
extern void getp_multi(int k2) ;

int getseed(int i,int pam,int c){
  int len;
  len=(i-pam)/c+1;
  return len;
}


int Index(char *S, char *T)                    //S为序列，T为seed
{
  int Slength = strlen(S);            //获得主串S的长度
  int Tlength = strlen(T);            //获得子串T的长度
  int i = 0;                                //记录主串S当前位置
  int j = 0;                                //记录子串T当前位置
  int k = 0;
  while(i < Slength){
    while( i < Slength&&j < Tlength)        //确保两个字符串的当前位置均小于其长度
      {
	if(S[i]=='\n')
	  i++;
	else if(S[i] == T[j])                    //判断主串S当前位置与子串T当前位置的字符是否相
	  {
	    i++;                            //主串S的当前位置加1（后移）
	    j++;                            //子串T的当前位置加1（后移）
	  }
	else                                //如果两字符串的当前位置字符不等
	  {
	    i = i - j + 1;                    //主串S的当前位置i回溯到j==0时i位置的下一位置
	    j = 0;                            //子串T的当前位置j归0
	  }
      }        //循环比较完毕
    if(j == Tlength){                        //判断位置j的数值是否与子串T的长度相等
      p[k].spot=(i - Tlength);                    //若是，说明搜索成功，返回T在S中出现的首位置
      i=i-j+1;
      j=0;
      k++;
    }
  }
  printf("Find %d off-target spot.\n",k);
  return k;
}

void getp(char *S,char *T,int num,int mode,float floor)//S为序列，T为seed后
{
  float pb;
  int misspot;//位点定位
  int spot1;
  int i,j=0;
  int mismatch=0;
  for( i=0;i<num;i++){
        printf("\rSites have been calculated: %2d/%2d\t", i+1,num); 
		fflush(stdout);
	    if(mode==1)//基因组测序链上的匹配		  
	      spot1=p[i].spot+11;
	    else if(mode==2)//另一条链上的匹配
	      spot1=p[i].spot-9;
	    mismatch=0;
	    j=0;
	    while(j<9){
	      if(S[spot1]==T[j]){
		  		 spot1++;
				 j++;
	      }
	      else{
		  	   mismatch++;
		  	   spot1++;
	  		   misspot=j+12;
		  	   p[i].spotnum[mismatch]=misspot;
	   		   j++;
	      }
	
	    }
	    if(mode==2)
	      misspot=31-misspot;
	    p[i].mis=mismatch;
	    if(mismatch==1){
	      pb=0.74/(1+pow(2.71828,-((misspot-11)*1.8)));
	      p[i].po=pb;
	    }
	    else if(mismatch==0){
	      p[i].po=0.74;
	    }
	    else{
	     getp_multi(i);
	    }
  }
  for(i=0;i<num;i++){
    printf("Sequence= ");
    if(mode==1&&p[i].po>=floor){
      for(spot1=p[i].spot;spot1<=p[i].spot+19;spot1++){
	printf("%c",S[spot1]);
      }
      printf(" Spot=%d,Mismatch number=%d,Possibility=%f\n",p[i].spot+1,p[i].mis,p[i].po);
    }
    else if(mode==2&&p[i].po>=floor){
      for(spot1=p[i].spot-9;spot1<=p[i].spot+10;spot1++){
	printf("%c",S[spot1]);
      }
      printf(" Spot=%d,Mismatch number=%d,Possibility=%f\n",p[i].spot-8,p[i].mis,p[i].po);
    }
  }
}

void getp_multi(int k2){// 
  float On_next=0.206;//匹配时后退概率
  float Off_next=0.997979;//错配时后退 概率
  double d=0.0;//0-1的随机数
  time_t t=time(NULL);//设定随机数种子     
  srand(t); 
  int k=1;//结合位点 
  int flag=0;//是否错配 
  int combine_num=0;//结合次数
  int sim_num;//模拟次数
  int sim_top=100000;//模拟上限 
  int step_top=1000;//步数上限 
  int i=0;
  int j=0;
  for(sim_num=0;sim_num<sim_top;sim_num++){
  	  if(sim_num%100==0);{
  	    printf(" Process:%2d%%\b\b\b\b\b\b\b\b\b\b\b\b", 100*sim_num/sim_top); 
        fflush(stdout);
  }
	  for(i=0;i<=step_top;i++){//步数计数 
	   	  if(k>0&&k<=20){
	         flag=0;
			 for(j=0;j<p[k2].mis;j++){
	            k==p[k2].spotnum[j]?(flag=1):(flag=flag);
	        }
			 d=((double)rand())/RAND_MAX;//随机数获取
			 if(flag==1)//不匹配 
			 	d-Off_next>0?(k=k+1):(k=k-1);
			 else if(flag==0)//匹配
			    d-On_next>0?(k=k+1):(k=k-1);		 		
	       }
	      else if (k>20){ 
	          combine_num++;
	          k=1;
			  break;
			  } 
		  else if (k<=0){ 
		      k=1;
		      break;
		      } 
	   }
	}   
   float multi_p=(float)combine_num/sim_top;
   p[k2].po=multi_p;
   }


int revbase(char *s){//之前有三个输入参数，用于确定开始位点以及反转长度 
  int x;
  int len=strlen(s); 
  int i=0;
  char rev[20];
  for(x=len;x>=0;x--){
    if (s[x]=='A')
      rev[i++]='T';
    else if (s[x]=='G')
      rev[i++]='C';
    else if (s[x]=='C')
      rev[i++]='G';
    else if (s[x]=='T')
      rev[i++]='A';		  	
  }
  rev[i]='\0';
  return rev;
}

int main(){
  // int delta_I=;
  // int delta_C=;
  // int delta_PAM=;
  // int delta_CLV=;
  // int seed=getseed(delta_I,delta_PAM,delta_C);
  printf("Welcome to off-target detector V1.5\nEnable the calculation of multi-off-target\nPlease name the name of fasta file as A.fa\n");
  char seed[12];
  char back[10];
  float floor=0.37;//检测下限 
  printf("please input the seed sequence(seed=11) \n");
  scanf("%s",seed);
  printf("please input the latter sequence \n");
  scanf("%s",back);
  printf("please set the floor of scanning.(Default 0.37)");
  scanf("%f",&floor);
  char seq[21];
  int x;
  int i=0;
  char *str = (char *)malloc(sizeof(char)*5000000);
  FILE *fa;
  fa = fopen("A.fa","r");
  int c;
  int num;
  i=0;
  while((c=fgetc(fa))!=EOF){
    if(c!='\n'){
      str[i]=c;
      i++;}
  }
  str[i]='\0';
  num=Index(str,seed);
  getp(str,back,num,1,floor);
  printf("Target=%s%s\n",seed,back);
  printf("Sequence on the other strand %s",revbase(back));
  printf("%s\n",revbase(seed));
  printf("Mismatch on another strand:\n");
  num=Index(str,revbase(seed));
  getp(str,revbase(back),num,2,floor);
  printf("Press q to quit; Press r to enter a new sequence");
  c=getch();
  if(c=='q')
    return 0;
  else if(c=='r'){
    printf("\n");
    main();
  }
}
       

