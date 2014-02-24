#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>
#include<math.h>
#define N 6
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
    putchar('\n');putchar('\n');
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
       return a[3]*c.x*c.x + a[4]*c.y*c.y + a[5]*c.x*c.y + 
                           a[1]*c.x + a[2]*c.y + a[0];
}

double grad_x(double* a, vec2 c)
{
       return 2.0*a[3]*c.x + a[5]*c.y + a[1];
}

double grad_y(double* a, vec2 c)
{
       return 2.0*a[4]*c.y + a[5]*c.x + a[2];
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
	double c[2]={	(*x)*axes[0]+(*y)*axes[1],
			(*x)*axes[2]+(*y)*axes[3]};
	*x=c[0];*y=c[1];
}
void fromRandomAxes(double* x, double* y, double *a)
{
	double det = a[0]*a[3]-a[1]*a[2];
	double axesR[4]={   a[3]/det, - a[1]/det,
			   -a[2]/det,   a[0]/det};
	
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
 /*   for (int i=1; i<3; i++)
	if (coord[i].y < coord[min].y) min=i;
    temp=coord[min]; coord[min]=coord[0]; coord[0]=temp;
    t=f[0]; f[0]=f[min]; f[min]=t;
    t=g_x[0]; g_x[0]=g_x[min]; g_x[min]=t;
    t=g_y[0]; g_y[0]=g_y[min]; g_y[min]=t;
    t=g_z[0]; g_z[0]=g_z[min]; g_z[min]=t;
    for (int i=1; i<4; i++)
	if (coord[i].y > coord[max].y) max=i;
    temp=coord[max]; coord[max]=coord[3]; coord[3]=temp;
    t=f[3]; f[3]=f[max]; f[max]=t;
    t=g_x[3]; g_x[3]=g_x[max]; g_x[max]=t;
    t=g_y[3]; g_y[3]=g_y[max]; g_y[max]=t;
    t=g_z[3]; g_z[3]=g_z[max]; g_z[max]=t;*/
		//for (int i=0; i<4; i++) out(coord[i]);putchar('\n');putchar('\n');

    vec2 point,average={0,0};//,max=coord[0];                      // point inside the thetraedr
    for (int i=1; i<3; i++)                                          //
    {                                                                //
	average.x += coord[i].x;
	average.y += coord[i].y;
    }                                                                //
    average.x /= 4.0;
    average.y /= 4.0;
	
    point = average;
	vec2 pointN=point;

    double   	x0=coord[0].x, x1=coord[1].x, x2=coord[2].x,  //
       		y0=coord[0].y, y1=coord[1].y, y2=coord[2].y;
        double f0=g_x[0],f1=g_x[1],f2=g_x[2],
		f0y=g_y[0],f1y=g_y[1],f2y=g_y[2]; //solve for 4x4
	double znam = (y0-y1)*(y0-y2)*(y1-y2);
	double axes[4]={0.0};//,crd[2]={_crd[0],_crd[1]};
	
    	if (fabs(znam) < 0.000001)  //coordinate transform
	{
		double 	l01 = scalar(x0-x1,y0-y1,x0-x1,y0-y1),
			l02 = scalar(x0-x2,y0-y2,x0-x2,y0-y2),
			l12 = scalar(x2-x1,y2-y1,x2-x1,y2-y1),
		//shift
			mid[2] = {(x0+x1+x2)/3.0,(y0+y1+y2)/3.0};
		x0 -= mid[0];		y0 -= mid[1];		
		x1 -= mid[0];		y1 -= mid[1];
		x2 -= mid[0];		y2 -= mid[1];
		pointN.x -= mid[0];	pointN.y -= mid[1];
		//size
		double norm = 1.0;
		if (l01 >= l02 && l01 >= l12)
			 norm = l01/10.0; 
		if (l02 >= l01 && l02 >= l12)
			 norm = l02/10.0; 
		if (l12 >= l01 && l12 >= l02)	
			 norm = l12/10.0; 
		x0 /= norm;		y0 /= norm;
		x1 /= norm;		y1 /= norm;
		x2 /= norm;		y2 /= norm;
		pointN.x /= norm;		pointN.y /= norm;
		f0*=norm; f0y*=norm;
		f1*=norm; f1y*=norm;
		f2*=norm; f2y*=norm;

		//rotate
		if (l01 >= l02 && l01 >= l12)
			{ axes[0]=x2-(x0+x1)/2.0; axes[1]=y2-(y0+y1)/2.0; };
		if (l02 >= l01 && l02 >= l12)
			{ axes[0]=x1-(x0+x2)/2.0; axes[1]=y1-(y0+y2)/2.0; };
		if (l12 >= l01 && l12 >= l02)	
			{ axes[0]=x0-(x2+x1)/2.0; axes[1]=y0-(y2+y1)/2.0; };

		double 
		norma = sqrt(axes[0]*axes[0]+axes[1]*axes[1]);//scalar(axes+2,axes+2));	
		axes[0]/=norma; axes[1]/=norma; 	
		if (fabs(axes[0]) > 0.0)
		{
			axes[3] = 1.0;
			axes[2] = -axes[1]/axes[0];
			//printf("norm: %lf     %lf %lf\n",norma, axes[2], axes[3]);
		}
		else if (fabs(axes[1]) > 0.0)
		{
			axes[2] = 1.0;
			axes[3] = -axes[1]/axes[0];
		}
		else 
		{
			printf("Fail in determining axis direction");
			return 0.0;
		}
		norma = sqrt(scalar(axes,axes));
		axes[2]/=norma; axes[3]/=norma; 

	    	//if (fabs(znam) < 0.000001) 
		//printf("%4d %4d %4d before: %10lf %10lf %10lf\t%10lf %10lf %10lf\t%10lf\n",t->vert[0],t->vert[1],t->vert[2],y0,y1,y2,y0-y1,y0-y2,y1-y2, znam);
		intoRandomAxes(&x0,&y0,axes);
		intoRandomAxes(&x1,&y1,axes);
		intoRandomAxes(&x2,&y2,axes);
		intoRandomAxes(&(pointN.x),&(pointN.y),axes);
		intoRandomAxes(&f0,&f0y,axes);
		intoRandomAxes(&f1,&f1y,axes);
		intoRandomAxes(&f2,&f2y,axes);
	    	//if (fabs(znam) < 0.000001) 
		//printf("%4d %4d %4d after : %10lf %10lf %10lf\t%10lf %10lf %10lf\t%10lf\n\n",t->vert[0],t->vert[1],t->vert[2],y0,y1,y2,y0-y1,y0-y2,y1-y2, znam);

		znam = (y0-y1)*(y0-y2)*(y1-y2);
	    	if (fabs(znam) < 0.000000000001) 
		{
			printf("%10lf<- skipping by znam\n%10lf %10lf %10lf \n%10lf %10lf %10lf \n",znam, x0,x1,x2,y0,y1,y2); continue;		
		}
	}

	znam = x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2;
    	if (fabs(znam) < 0.00000001) 
	{
		printf("%10lf<- skipping by volume\n%10lf %10lf %10lf \n%10lf %10lf %10lf \n",znam, x0,x1,x2,y0,y1,y2); return 0.0;		
	}
        coeff_num[1] = 	(f2*x1*y0-f1*x2*y0-f2*x0*y1+f0*x2*y1+f1*x0*y2-f0*x1*y2)/znam;
        coeff_num[3] = 	(y0*(f1-f2)+y1*(f2-f0)+y2*(f0-f1))/znam/2.0;
        coeff_num[5] = 	-(f2*(x1-x0)+f1*(x0-x2)+f0*(x2-x1))/znam;
        double 	f3=f[0]-coeff_num[1]*x0-coeff_num[3]*x0*x0-coeff_num[5]*x0*y0,
		f4=f[1]-coeff_num[1]*x1-coeff_num[3]*x1*x1-coeff_num[5]*x1*y1,
		f5=f[2]-coeff_num[1]*x2-coeff_num[3]*x2*x2-coeff_num[5]*x2*y2; //solve for 4x4
	znam = (y0-y1)*(y0-y2)*(y1-y2);
        coeff_num[0] = 	(f5*y0*(y0-y1)*y1+y2*(f3*y1*(y1-y2)+f4*y0*(y2-y0)))/znam;
        coeff_num[2] = 	(f5*(y1+y0)*(y1-y0)+f4*(y0-y2)*(y0+y2)+f3*(y2-y1)*(y2+y1))/znam;
        coeff_num[4] = 	(f5*(y0-y1)+f3*(y1-y2)+f4*(y2-y0))/znam;
    
    
    double f_p_exact,f_p_numer;
        f_p_exact = func(coeff_exact,point);                      //exact function in point
        f_p_numer = func(coeff_num  ,pointN);                      //numer function in point
        d=fabs((f_p_exact-f_p_numer)/f_p_exact);
        if (d > dif4_max) dif4_max = d;
	dif4_average += d;
	if (d > 1.000001)
	{
		printf("%6d::  E: %10f  N: %10f  D: %10f\n",all,f_p_exact,f_p_numer,d);
		for (int i=0; i<3; i++) out(coord[i]);
		out(coeff_exact); out(coeff_num);
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
