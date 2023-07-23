# include<iostream>
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>

using namespace std;
const int Q=9;          //D2Q9模型   //在c++语法中，由于这里的数字被写在了e[Q][2]作为数组的表达，所以只能作为常量 const ,而不能直接写成 int Q = 9; 
						//在 D2Q9 模型中，格子声速 Cs = C/sqrt(3)  
const int NX=256;       //x方向
const int NY=256;        //y方向
const double U=0.1;      //顶盖速度

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}; //构成9行2列的矩阵 ， e[9][2]; 
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};

//==================================================数组的定义================================================================================// 
double rho[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],F[NX+1][NY+1][Q],f[NX+1][NY+1][Q]; //u[NX+1][NY+1][0]表示流场内速度Ux,u[NX+1][NY+1][1]表示流场内速度Uy 
double M[NX+1][NY+1];	
double p[NX+1][NY+1];

//温度分布函数 
double T[NX+1][NY+1];
double FT[NX+1][NY+1][Q];
double ft[NX+1][NY+1][Q];	

//qx、qy热流线矢量

double qx[NX+1][NY+1],qy[NX+1][NY+1]; 

																				//[NX+1][NY+1]起到的是定位的作用，具体这个点对应哪个方向则是[0]或者[1]; 
//rho[NX+1][NY+1],标量值,流场中每个密度; 													//这里的u0用于判断收敛程度(更新后的 [U] -更新前的 [U0]) 
//u[NX+1][NY+1][2],矢量值，流场中每个速度值;
//f[NX+1][NY+1][Q], 每个点(x,y)对应的 Q(9)个方向的分布函数值; 


//=========================================================================================================================================//

int i,j,k,ip,jp,n; //定义了 int 类型的整数全局0变量 
double c,cs,Re,dx,dy,Lx,Ly,dt,rho0,T0,P0,tau_f,niu,error; //定义了 double 类型的全局变量 

void init();//对init()函数声明 

double feq(int k,double rho,double u[2]); //对feq()函数声明; 

double feqT(int k,double T,double u[2]); // 对feqT()函数声明; 

void evolution();

void output(int m);

void Error();



//============================================主程序======================================================// 
 int main()
 {
 	using namespace std;
 	init();//第一个执行函数 
 	for(n=0;;n++)
 	{
 		evolution();
 		if(n%1000==0)
 		{
 			Error();
 		cout<<"THE"<<n<<"th computation result:"<<endl<<"The u,v of point(NX/2,NY/2)is:"<<setprecision(6)<<u[NX/2][NY/2][0]<<","<<u[NX/2][NY/2][1]<<endl;
 		cout<<"The max relative error of uv is:"<<setiosflags(ios::scientific)<<error<<endl;
 		if(n>=1000)
 		{
 			if(n%10000==0)     output(n);
			if(error<10.e-9)  break;
			
	 		 }
		 }
	 }
	 return 0;
 }
 
//=======================================================================================================//  
 
 
//===========================================第一个执行函数================================================// 
 void init() //在该函数内，存在feq()函数 
 {
 	dx=1.0; //一个单元网格的长度为1 
 	dy=1.0;
 	Lx=dx*double(NX); //边长x方向的总长度 Lx = x方向的单位长度(dx)* x方向总的单元格数目(NX);
 	Ly=dy*double(NY);
 	dt=dx;
 	c=dx/dt;//1.0
 	cs = c/sqrt(3);
 	rho0=1.0;
 	T0 = 1.0;
 	Re=1000;
 	
 	//==================================粘度和无量纲松弛时间计算====================================================// 
 	niu=U*Lx/Re;//运动粘度系数 
 	tau_f=3.0*niu+0.5;
 	//======================================================================================================//
 
 	std::cout<<"tau_f="<<tau_f<<endl;
 	
 	for(i=0;i<=NX;i++) // 初始化流场内的每一个点 
	 for(j=0;j<=NY;j++) 
	 {
	 	u[i][j][0]=0;   //初始化流场内每个点的x方向的速度;
	 	u[i][j][1]=0;   //初始化流场内每个点的y方向的速度; 
	 	M[i][j] = sqrt(	u[i][j][0]*	u[i][j][0]+u[i][j][1]*u[i][j][1]); //计算矢量合成后的 U 
	 	rho[i][j]=rho0;
	 	p[i][j] = rho[i][j]*cs*cs;
	 	u[i][NY][0]=U;  //对顶盖的速度进行初始化; 
	 	M[i][NY] = U; 
	 	
	 	T[i][j] = T0;//对流场内温度初始化
		T[i][NY] = 330;  //对顶盖温度初始化:330K 
		
		//对热流线矢量初始化
		
		qx[i][j] = 0;
		qy[i][j] = 0; 
		
		
	 	for(k=0;k<Q;k++)
	 	{
	 		f[i][j][k]=feq(k,rho[i][j],u[i][j]);
   			ft[i][j][k] =feqT(k,T[i][j],u[i][j]); 		
			                                      //对流场内的每个点的每9个速度方向(f[i][j][k]进行初始化)，初始化的方法是调用了函数feq(int k, double rho, double u[2]) 
		}									     //从 double feq(int k, double rho, double u[2])函数来看，返回的是一个double类型，值为feq;	 
	  } 										//分布函数是一个小于1的概率(全局比例)，故 feq <= 1; 
  } 											//参数列表中的 u[i][j] 是一个地址，在feq()函数的参数列表中，double u[2]则表示 u[0]、u[1];两者合并就是u[i][j][2]; 
												//该形式的基本类型如下： 
//===========================================================// 
/* [该程序是一个可独立运作的程序]
#include<iostream>
using namespace std;
const int NX = 256;
const int NY = 256;
double u[NX+1][NX+1][2];
void feq(double u[2]);

void feq(double u[2])
{
	
	cout << u[1] << endl;
		
}

int main()
{
	
	
	for(int i=0;i<=NX;i++)
	{
		
		for(int j=0; j<=NY;j++)
		{
		
			u[i][j][0] = 0;
			
			u[i][j][1] = 1;
			
			feq(u[i][j]);		
		}
		
		
	}

	return 0;
	
}



*/
//===========================================================//



//==============================================================================================================//  

  
  
  
//====================================================init()函数内，调用了 feq()函数, 并且该函数是一个返回double数的函数===============================================// 
  
 double feq(int k,double rho,double u[2])//计算平衡态分布函数
 {
 	double eu,uv,feq;
 	//标量eu 
	 eu=(e[k][0]*u[0]+e[k][1]*u[1]);  //二维数组e[Q][2],e[k][0]表示第1列的 k 行(速度方向中的x方向); e[k][1]则对应表示为速度方向的 y 方向; 
	 								  //结合上面的分析，u[0]就是u[i][j][0];u[1]就是u[i][j][1]; 
 									  //eu表示的就是   x方向*x方向的速度值 + y方向*y方向的速度值; 
	 								  //就是每个点上的速度(ux,uy)在e[9][2]离散速度上的投影，对应就是：对应横/纵坐标相乘再相加 
	 //标量uv 
	 uv=(u[0]*u[0]+u[1]*u[1]);      //类似于矢量模的平方 |vector| = sqrt( ux*ux + uy*uy); 
 	
	 
	 feq=w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);  //w[k]是一个权重系数，即每个方向(Q=9,9个方向)的权重值; 
 	
	 return feq;
  } 
  //平衡态的温度函数 
 double feqT(int k,double T,double u[2])
 {
 	double eu,uv,feqT;
 	eu=(e[k][0]*u[0]+e[k][1]*u[1]);
 	uv=(u[0]*u[0]+u[1]*u[1]); 
 	feqT=w[k]*T*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
 	return feqT;
  } 

  //============================对于返回的feq值，则是全部小于1的值，下面是 返回 feq 完整独立程序=================================// 
  /*
  #include<iostream>
using namespace std;
const int NX = 256;
const int NY = 256;
const int Q = 9;
double u[NX+1][NX+1][2];
double rho[NX+1][NY+1];
int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}; 
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
void feq(int k, double rho, double u[2]);


void feq(int k, double rho, double u[2])
{
	double eu,uv,feq;
	
	eu = (e[k][0]*u[0]+e[k][1]*u[1]);
	uv = (u[0]*u[0]+u[1]*u[1]);
	
	feq = w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv); 
	
	cout << feq << endl;
		
}

int main()
{
	
	
	for(int i=0;i<=NX;i++)
	{
		
		for(int j=0; j<=NY;j++)
		{
		
			u[i][j][0] = 0;
			rho[i][j] = 1;
			u[i][j][1] = 1;
		for(int k=0;k<Q;k++)
	     	{
			
			feq(k,rho[i][j],u[i][j]);		
		
	        }
		}
		
		
	}

	return 0;
	
}
  
  
  */
 //===================================================================================================================================// 
  
  
  void evolution()
  {
  	for(i=1;i<NX;i++)//演化
	  for(j=1;j<NY;j++)
	    for(k=0;k<Q;k++)
	    {
	    	ip=i-e[k][0];
	    	jp=j-e[k][1]; //这一步的ip与jp是一定存在重复的
	    	              //在后续的更新分布函数中 F[i][j][k]，表示的是位置矢量由 r 变为 r+dr、t 变为 dt; ip、jp则是位置更新后的坐标；
						  /*
									Boltzman-BGK:f2(r+e*dt,t+dt)-f(r,t)    f2表示的是更新后(dt后)的分布函数，f1是更新前(dt前)的分布函数，e就是e[Q][2],对于 x方向则e[Q][0],对于 y 方向，则e[Q][1]；
									                                       其中e*dt,表示位移量;						  
 						  
						  */
//==========================验证ip与jp重复性的独立程序如下==================================//
/*
#include<iostream>
using namespace std;
const int NX = 256;
const int NY = 256;
const int Q = 9;
int ip,jp;
int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
int main()
{
	for(int i=1;i<NX;i++)
	{
		
		for(int j=1;j<NY;j++)
		{
			
			for(int k=0;k<Q;k++)
			{
				ip=i-e[k][0];
	         	jp=j-e[k][1];
				
				cout << ip <<'\t' << jp << endl;
				
			}
		}
		
	}
	
	return 0;
 } 

*/
//========================================================================================// 
	    
			F[i][j][k]=f[ip][jp][k]+(feq(k,rho[ip][jp],u[ip][jp])-f[ip][jp][k])/tau_f;//F[i][j][k]演化后的分布函数 
			FT[i][j][k] = ft[ip][jp][k]+(feqT(k,T[ip][jp],u[ip][jp])-ft[ip][jp][k])/tau_f;
//===========由上述的feq()函数，发现 feq(k,rho[ip][jp],u[ip][jp]) = f[ip][jp][k]，完整的证明程序如下====================================//		
		/*
		
# include<iostream>
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>

using namespace std;
const int Q=9;          //D2Q9模型   //在c++语法中，由于这里的数字被写在了e[Q][2]作为数组的表达，所以只能作为常量 const ,而不能直接写成 int Q = 9; 
						//在 D2Q9 模型中，格子声速 Cs = C/sqrt(3)  
const int NX=256;       //x方向
const int NY=256;        //y方向
const double U=0.1;      //顶盖速度

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}; //构成9行2列的矩阵 ， e[9][2]; 
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};

//==================================================数组的定义================================================================================// 
double rho[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],F[NX+1][NY+1][Q],f[NX+1][NY+1][Q]; //u[NX+1][NY+1][0]表示流场内速度Ux,u[NX+1][NY+1][1]表示流场内速度Uy 
																							//[NX+1][NY+1]起到的是定位的作用，具体这个点对应哪个方向则是[0]或者[1]; 
//rho[NX+1][NY+1],标量值,流场中每个密度; 													//这里的u0用于判断收敛程度(更新后的 [U] -更新前的 [U0]) 
//u[NX+1][NY+1][2],矢量值，流场中每个速度值;
//f[NX+1][NY+1][Q], 每个点(x,y)对应的 Q(9)个方向的分布函数值; 


//=========================================================================================================================================//

int i,j,k,ip,jp,n; //定义了 int 类型的整数全局0变量 
double c,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error; //定义了 double 类型的全局变量 

void init();//对init()函数声明 

double feq(int k,double rho,double u[2]); //对feq()函数声明; 

void evolution();

double feq(int k,double rho,double u[2])
{
	double eu,uv,feq;
	eu = (e[k][0]*u[0]+e[k][1]*u[1]);
	uv = (u[0]*u[0]+u[1]*u[1]);
	feq = w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
    return feq;
}

int main()
{
	for(int i=0;i<=NX;i++)
	{
		for(int j=0;j<=NY;j++)
		{
			u[i][j][0] = 0;
			u[i][j][1] = 0;
			rho[i][j] = 1.0;
			u[i][NY][0] = 50;
			for(int k=0;k<Q;k++)
			{
				
				f[i][j][k] = feq(k,rho[i][j],u[i][j]);
			
			}
		
		}
		
	
	 } 
	
	for(int i=1;i<NX;i++)
	{
	
	for(int j=1;j<NY;j++)
	{
	
	for(int k=0;k<Q;k++)
	{
		ip = i-e[k][0];
		jp = i-e[k][1];
		double pp;
		pp = feq(k,rho[ip][jp],u[ip][jp])-f[ip][jp][k];
		if(pp==0)
		
		cout << "相等" << endl;
		
	 }
}
}
	

	return 0;
}
				
		*/
		//============================================================================================================//
		}
	  for(i=1;i<NX;i++)//计算宏观量(内部流场) 
	    for(j=1;j<NY;j++)
		{
			
			//================================给予一个初始值================================// 
			u0[i][j][0]=u[i][j][0];
			u0[i][j][1]=u[i][j][1];
			rho[i][j]=0;
			u[i][j][0]=0;
			u[i][j][1]=0;
			M[i][j] = sqrt(	u[i][j][0]*	u[i][j][0]+u[i][j][1]*u[i][j][1]);
			p[i][j] = rho[i][j]*cs*cs;
			T[i][j] = 0;
			qx[i][j] = 0;
			qy[i][j] = 0;
		   //==============================================================================//
		   	
			for(k=0;k<Q;k++)
			{
				f[i][j][k]=F[i][j][k];
				ft[i][j][k] = FT[i][j][k];  //分布函数更新赋值 
				//===================================物理宏观量================================================// 
				rho[i][j]+=f[i][j][k];
				T[i][j]+=ft[i][j][k];//温度宏观量相加 
				//温度在x、y方向的一阶偏导数;
				 
				qx[i][j] = (T[i+1][j]-T[i-1][j])/(2*dx);
				qy[i][j] = (T[i][j+1]-T[i][j-1])/(2*dy);
				
				u[i][j][0]+=e[k][0]*f[i][j][k];
				u[i][j][1]+=e[k][1]*f[i][j][k]; 
			    
			}
			u[i][j][0]/=rho[i][j];
			u[i][j][1]/=rho[i][j];
			M[i][j] = sqrt(	u[i][j][0]*	u[i][j][0]+u[i][j][1]*u[i][j][1]);
			p[i][j] = rho[i][j]*cs*cs;
	   //============================================================================================//
	    } 
	 //边界处理 
	 
	 //=====================================左右边界=====================================================// 
	   for(j=1;j<NY;j++)
	     for(k=0;k<Q;k++)
		 {
		 	rho[NX][j]=rho[NX-1][j];
		 	T[NX][j] = T[NX-1][j]; 
		 
		 //右壁面设置为绝热壁面 
			 qx[NX][j] = 0;
			 qy[NX][j] = 0;
			  
		 	p[NX][j] = rho[NX][j]*cs*cs;
		 	f[NX][j][k]=feq(k,rho[NX][j],u[NX][j])+f[NX-1][j][k]-feq(k,rho[NX-1][j],u[NX-1][j]); //更新边界上的分布函数 
		 	ft[NX][j][k] = feqT(k,T[NX][j],u[NX][j])+ft[NX-1][j][k]-feqT(k,T[NX-1][j],u[NX-1][j]);
		 	
			rho[0][j]=rho[1][j];
			T[0][j] = T[1][j];
		
		//左壁面设置为绝热壁面 	
			qx[0][j] = 0;
			 qy[0][j] = 0;
			 
		 	p[0][j] = rho[0][j]*cs*cs;
		 	f[0][j][k]=feq(k,rho[0][j],u[0][j])+f[1][j][k]-feq(k,rho[1][j],u[1][j]);
		 	ft[0][j][k]=feqT(k,T[0][j],u[0][j])+ft[1][j][k]-feqT(k,T[1][j],u[1][j]);
		 	
	//====================================================================================================//	 	
		 }
		for(i=0;i<=NX;i++)//上下边界
		  for(k=0;k<Q;k++)
		  {
		  	rho[i][0]=rho[i][1];
		  	T[i][0] = T[i][1];
		  	
		  	//下壁面设置为绝热壁面
			  
			  qx[i][0] = 0;
			  qy[i][0] = 0;
		  	
		  	p[i][0] = rho[i][0]*cs*cs;
		  	f[i][0][k]=feq(k,rho[i][0],u[i][0])+f[i][1][k]-feq(k,rho[i][1],u[i][1]);
		  	ft[i][0][k]=feqT(k,T[i][0],u[i][0])+ft[i][1][k]-feqT(k,T[i][1],u[i][1]);
		  	
		  	rho[i][NY]=rho[i][NY-1];
		  //	T[i][NY]=T[i][NY-1];
		  		p[i][NY] = rho[i][NY]*cs*cs;
		  	  u[i][NY][0]=U;
		  	  T[i][NY] = 330;
		  	
			  //上壁面设置为绝热壁面  
		  	  qx[i][NY] = 0;
				qy[i][NY] = 0; 
		  	  
		  	  M[i][NY] = U; 
		  	  f[i][NY][k]=feq(k,rho[i][NY],u[i][NY])+f[i][NY-1][k]-feq(k,rho[i][NY-1],u[i][NY-1]);
		  	  ft[i][NY][k]=feqT(k,T[i][NY],u[i][NY])+ft[i][NY-1][k]-feqT(k,T[i][NY-1],u[i][NY-1]);
			}  
		  
		 
  }
  
  void output(int m)//输出
  {
  	ostringstream name;
	  name<<"cavity_"<<m<<".dat";
	  ofstream out(name.str().c_str());
	  out<< "Title= \"LBM Lid Driven Flow\"\n" << "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"M\",\"p\",\"T\",\"qx\",\"qy\",\"rho\"\n" << "ZONE T=\"BOX\",I=" << NX+1 << ",J=" << NY+1 << ",F=POINT" << endl;
	  for(j=0;j<=NY;j++)
	     for(i=0;i<=NX;i++)
	     {
	     	out<<double(i)/Lx<<" "<<double(j)/Ly<<" "<<u[i][j][0]<<" "<<u[i][j][1]<<" "<<M[i][j]<<" "<< p[i][j] << " "<<T[i][j]<<" "<<qx[i][j]<<" "<< qy[i][j]<<" "<<rho[i][j] <<endl;
		 }
   } 
   
   void Error()
   {
   	double temp1,temp2;
   	temp1=0;
   	temp2=0;
   	for(i=1;i<NX;i++)
   	   for(j=1;j<NY;j++)
   	   {
   	   	temp1 +=(
		(u[i][j][0]-u0[i][j][0])*(u[i][j][0]-u0[i][j][0])+(u[i][j][1]-u0[i][j][1])*(u[i][j][1]-u0[i][j][1]));
		  temp2 +=(u[i][j][0]*u[i][j][0]+u[i][j][1]*u[i][j][1]);
		  }
		  temp1=sqrt(temp1);
		  temp2=sqrt(temp2);
		    error=temp1/(temp2+1e-30);
   }
