#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>
#include<math.h>
#define N 10
#define T 100000000	

struct vec2 
{
       double x,y;
};

void out(double* _a)
{
    putchar('\n'); 
    for (int i=0; i<N; i++)
        printf("%f ",_a[i]);
    putchar('\n');
}

void out(double* _a, int n)
{
    putchar('\n'); 
    for (int i=0; i<n; i++)
        printf("%f ",_a[i]);
    putchar('\n');
}
void out(vec2 v)
{
	printf("%10lf %10lf\n",v.x,v.y);
}

/*double coeff_dif4(double* _a1, double * _a2)
{
    double res;
    for (int i=0; i<N; i++)
        printf("%f ",_a[i]);
    putchar('\n');putchar('\n');
}*/

double rand_1()
{
       return ((double)(rand()%5000000))/5000000.0;
}

double func(double* a, vec2 c)
{
       return a[0]*c.x*c.x*c.x + a[1]*c.y*c.y*c.y + a[2]*c.y*c.x*c.x + a[3]*c.y*c.y*c.x	+ a[4]*c.x*c.x + a[5]*c.y*c.y + a[6]*c.x*c.y + a[7]*c.x + a[8]*c.y + a[9];
		//a[3]*c.x*c.x + a[4]*c.y*c.y + a[5]*c.x*c.y + 
              	//             a[1]*c.x + a[2]*c.y + a[0];
}

double grad_x(double* a, vec2 c)
{
       return 3.0*a[0]*c.x*c.x + 2.0*a[2]*c.x*c.y + a[3]*c.y*c.y + 2.0*a[4]*c.x + a[6]*c.y + a[7];
	//2.0*a[3]*c.x + a[5]*c.y + a[1];
}

double grad_y(double* a, vec2 c)
{
       return 3.0*a[1]*c.y*c.y + a[2]*c.x*c.x + 2.0*a[3]*c.x*c.y + 2.0*a[5]*c.y + a[6]*c.x + a[8];
	//2.0*a[4]*c.y + a[5]*c.x + a[2];
}
double scalar(double* _v1, double* _v2)
{
	return _v1[0]*_v2[0]+_v1[1]*_v2[1];
}
double scalar(double _x0, double _y0, double _x1, double _y1)
{
	return _x0*_x1+_y0*_y1;
}
void intoRandomAxes(double* x, double* y, double *axes)
{
	double c[2]={	(*x)*axes[0]+(*y)*axes[2],
			(*x)*axes[1]+(*y)*axes[3]};
	*x=c[0];*y=c[1];
}
void intoRandomAxesGrad(double* x, double* y, double *axes)
{
	double c[2]={	(*x)*axes[0]+(*y)*axes[2],
			(*x)*axes[1]+(*y)*axes[3]};
	*x=c[0];*y=c[1];
}
void fromRandomAxes(double* x, double* y, double *a)
{
	double det = a[0]*a[3]-a[1]*a[2];
	double axesR[4]={   a[3]/det, - a[2]/det,
			  - a[1]/det,   a[0]/det};
	
//	printf("\n\nAR: \t%10lf  %10lf\n\t%10lf  %10lf\n\n",axesR[0],axesR[1],axesR[2],axesR[3]);
	intoRandomAxes(x,y,axesR);
}	
void fromRandomAxesGrad(double* x, double* y, double *a)
{
	double det = a[0]*a[3]-a[1]*a[2];
	double axesR[4]={   a[3]/det, - a[1]/det,
			  - a[2]/det,   a[0]/det};
//	printf("\n\nARG: \t%10lf  %10lf\n\t%10lf  %10lf\n\n",axesR[0],axesR[1],axesR[2],axesR[3]);
	intoRandomAxes(x,y,axesR);
}	
int main()
{
    printf("Interpolation polynom:\na1*x^2+a2*y^2+a3*z^2+a4*x*y*+a5*x*z+a6*y*z+a7*x+a8*y+a9*z+a10\n");
    double dif4_max=0.0, dif4_average, d;
    double coeff_exact[N];                                  
    srand(time(NULL));                                       

    for (int all=0; all<T; all++)
    {                                          
				//exact polynom coefficients         //                                    
    for (int i=0; i<N; i++)                                          //
        coeff_exact[i] = rand_1()*10;                                //
	coeff_exact[6] = 0.0;
        
    vec2 coord[3];                                                   //vertices
    double f[3], g_x[3], g_y[3];                             //exact function & gradient in vertices
    for (int i=0; i<3; i++)                                          //
    {                                                                //
        coord[i].x = rand_1();                                   //
        coord[i].y = rand_1();                                   //
        f[i]=func(coeff_exact,coord[i]);                             //
        g_x[i]=grad_x(coeff_exact,coord[i]);                         //
        g_y[i]=grad_y(coeff_exact,coord[i]);                         //
    }                                                                //
    
    
    double coeff_num[N];                                             //numerical coefficients
		//for (int i=0; i<4; i++) out(coord[i]);putchar('\n');
    int min=0,max=0;					//(y0-y3) appears in the following formulae,trying to 
    							//make the dif4erence as far from zero as possible
    vec2 temp; 
    double t;

    vec2 average={0,0};
    for (int i=0; i<3; i++)                                          //
    {                                                                //
	average.x += coord[i].x;
	average.y += coord[i].y;
    }                                                                //
    average.x /= 3.0;
    average.y /= 3.0;

	average.x += coord[1].x;
	average.y += coord[1].y;
	average.x /= 2.0;
	average.y /= 2.0;
	
    vec2 point = {average.x,average.y};
    vec2 pointN = {average.x,average.y};//=point;

    double   	x0=coord[0].x, x1=coord[1].x, x2=coord[2].x,  //
       		y0=coord[0].y, y1=coord[1].y, y2=coord[2].y;
        double 	f0=f[0],f1=f[1],f2=f[2],
		f0x=g_x[0],f1x=g_x[1],f2x=g_x[2],
		f0y=g_y[0],f1y=g_y[1],f2y=g_y[2]; //solve for 4x4
//	double znam = (y0-y1)*(y0-y2)*(y1-y2);
	double axes[4]={0.0};//,crd[2]={_crd[0],_crd[1]};
	
    	//if (fabs(znam) < 0.000001)  //coordinate transform

//		double 	l01 = scalar(x0-x1,y0-y1,x0-x1,y0-y1),
//			l02 = scalar(x0-x2,y0-y2,x0-x2,y0-y2),
//			l12 = scalar(x2-x1,y2-y1,x2-x1,y2-y1),
		//shift
		double	mid[2] = {x0,y0};
		x0 -= mid[0];		y0 -= mid[1];		
		x1 -= mid[0];		y1 -= mid[1];
		x2 -= mid[0];		y2 -= mid[1];
		pointN.x -= mid[0];	pointN.y -= mid[1];

		//to a simplex
		axes[0] = x1; axes[1] = x2;
		axes[2] = y1; axes[3] = y2;
		fromRandomAxes(&x0,&y0,axes);
		fromRandomAxes(&x1,&y1,axes);
		fromRandomAxes(&x2,&y2,axes);
		fromRandomAxes(&(pointN.x),&(pointN.y),axes);
		intoRandomAxesGrad(&f0x,&f0y,axes);
		intoRandomAxesGrad(&f1x,&f1y,axes);
		intoRandomAxesGrad(&f2x,&f2y,axes);

	vec2 new0,new1,new2,newN;
	new0.x=x0;new0.y=y0;
	new1.x=x1;new1.y=y1;
	new2.x=x2;new2.y=y2;
	double gr_x[3]={f0x,f1x,f2x}, gr_y[3]={f0y,f1y,f2y};

	coeff_num[9] = f0;
	coeff_num[8] = f0y;
	coeff_num[7] = f0x;
	coeff_num[6] = 0.0;//f1y - coeff_num[8];
	coeff_num[5] = 3.0*f2 - 2.0*coeff_num[8] - 3.0*coeff_num[9] - f2y;
	coeff_num[4] = 3.0*f1 - 2.0*coeff_num[7] - 3.0*coeff_num[9] - f1x;
	coeff_num[3] = f2x - coeff_num[7];//f2x - coeff_num[7];f2x - f1y + coeff_num[8] - coeff_num[7];//
	coeff_num[2] = f1y - coeff_num[8];
	coeff_num[1] = -2.0*f2 + coeff_num[8] + 2.0*coeff_num[9] + f2y;
	coeff_num[0] = -2.0*f1 + coeff_num[7] + 2.0*coeff_num[9] + f1x;
    
    
    double f_p_exact,f_p_numer;
        f_p_exact = func(coeff_exact,point);                      //exact function in point
        f_p_numer = func(coeff_num  ,pointN);                      //numer function in point
	double f_n[3]={func(coeff_num,new0),func(coeff_num,new1),func(coeff_num,new2)};
	double gx_n[3]={grad_x(coeff_num,new0),grad_x(coeff_num,new1),grad_x(coeff_num,new2)};
	double gy_n[3]={grad_y(coeff_num,new0),grad_y(coeff_num,new1),grad_y(coeff_num,new2)};
        d=fabs((f_p_exact-f_p_numer)/f_p_exact);
        if (d > dif4_max) dif4_max = d;
	dif4_average += d;
	if (d > 1.000001)
	{
		printf("%6d::  E: %10f  N: %10f  D: %10f\n",all,f_p_exact,f_p_numer,d);
		for (int i=0; i<3; i++) out(coord[i]);
		putchar('\n');
//		printf("vol: %lf",coord[1].x*coord[2].y+coord[0].x*coord[1].y+coord[0].y*coord[2].x-coord[1].x*coord[0].y-coord[2].x*coord[1].y-coord[0].x*coord[2].y);
//		out(coeff_exact); out(coeff_num);	
		out(new0);out(new1);out(new2);putchar('\n');out(point);out(pointN);
		out(f,3);out(f_n,3);putchar('\n');
		out(gr_x,3);out(gx_n,3);putchar('\n');//out(gr_x,3);putchar('\n');
		out(gr_y,3);out(gy_n,3);
		putchar('\n');putchar('\n');
	}
	if (!(all%100000)) printf("%10d|  Max: %10f  Average: %10f\n\n",all,dif4_max,dif4_average/T);
	//printf("%lf %lf %lf %lf\n",fabs(func(coeff_exact,coord[0])-func(coeff_num,coord[0])),fabs(func(coeff_exact,coord[1])-func(coeff_num,coord[1])),fabs(func(coeff_exact,coord[2])-func(coeff_num,coord[2])),fabs(func(coeff_exact,coord[3])-func(coeff_num,coord[3])));
    }
    ////double c_dif4=0.0;
    //for (int i=0; i<N; i++)
    //    c_dif4 += fabs(coeff_exact[i]-coeff_num[i])/fabs(coeff_exact[i]);
    dif4_average /= T;
    printf("\nMax: %10f  Average: %10f\n\n",dif4_max,dif4_average);
}
