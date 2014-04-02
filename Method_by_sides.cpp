	#include "Method.h"
#include "Triangle.h"
#include "Node.h"
#include "Mesh.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <stdio.h>
#define R_I 1000000
#define R_D 1000000.0

//double Method::axis[4]={0};
Method::Method()
{
}

Method::~Method()
{
}

int Method::init()
{
	srand(time(0));
	coeff[0] = coeff[1] = 0.2;
	coeff[1] = 0.53;
	order = 1;
	return 0;
}

int Method::init(double c0, double c1, int _order)
{
	srand(time(0));
	coeff[0] = c0; coeff[1] = c1;// coeff[] = 0.53;
	order = _order;
	return 0;
}

double Method::scalar(double* _v1, double* _v2)
{
	return _v1[0]*_v2[0]+_v1[1]*_v2[1];
}

double Method::scalar(double _x0, double _y0, double _x1, double _y1)
{
	return _x0*_x1+_y0*_y1;
}

void Method::randomizeAxis(double* axes)
{
	for (int i=0; i<4; i++)
		axes[i]=0.0;

	axes[0] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	axes[1] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	double norma = sqrt(scalar(axes,axes));
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
		axes[3] = -axes[0]/axes[1];
	}
	else 
	{
		printf("Fail in determining axis direction");
		return;
	}
	norma = sqrt(axes[2]*axes[2]+axes[3]*axes[3]);//scalar(axes+2,axes+2));	
	axes[2]/=norma; axes[3]/=norma; 
	//check!!!!!
	double check[3] = {scalar(axes,axes),scalar(axes+2,axes+2)
			,scalar(axes,axes+2)};
	if ((check[0] - 1.0)*(check[0] - 1.0) > 0.0001 || (check[1] - 1.0)*(check[1] - 1.0) > 0.0001 || check[2]*check[2] > 0.0001)
	{
		printf ("Warning! Recursion in axis random called! %lf %lf %lf   %lf %lf %lf %lf\n",check[0],check[1],check[2], axes[0], axes[1], axes[2], axes[3]);
		randomizeAxis(axes);
		return;
	}
}


void Method::calculateCoeff(double* _c, double* axes)
{
	_c[0] = coeff[0]*axes[0] + coeff[1]*axes[1];// + coeff[2]*axis[2];
	_c[1] = coeff[0]*axes[2] + coeff[1]*axes[3];// + coeff[2]*axis[5];
	//_c[2] = coeff[0]*axis[6] + coeff[1]*axis[7] + coeff[2]*axis[8];
}

double Method::interpolate_1_order(Triangle* t, double* _crd, int val, Mesh* mesh)
{	
	Node* nodes[3] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2])};
	/*int min[3]={0},max[3]={0},diff=0;
	for (int m=0; m<3; m++)
	{
		for (int i=1; i<4; i++)
			if (nodes[i]->coords[m] < nodes[min[m]]->coords[m]) min[m] = i;
	}
	if (min[1]) {Node* tmp = nodes[0]; nodes[0] = nodes[min[1]]; nodes[min[1]] = tmp;};
	for (int m=0; m<3; m++)
	{
		for (int i=1; i<4; i++)
			if (nodes[i]->coords[m] > nodes[max[m]]->coords[m]) max[m] = i;
	}
	if (max[1] != 3) {Node* tmp = nodes[3]; nodes[3] = nodes[max[1]]; nodes[max[1]] = tmp;};
*/
	double 	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], //
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1],   //
       		a,b,c;
        /*double axes[4], crd[2]={_crd[0],_crd[1]};
	randomizeAxis(axes);
	//printf("axes: %lf %lf %lf %lf\n",axes[0],axes[1],axes[2],axes[3]);
	intoRandomAxes(&x0,&y0,axes);
	intoRandomAxes(&x1,&y1,axes);
	intoRandomAxes(&x2,&y2,axes);
	intoRandomAxes(crd,crd+1,axes);*/
                                                   
        double f0=nodes[0]->u[val],f1=nodes[1]->u[val],f2=nodes[2]->u[val]; //solve for 4x4
	double znam = x1*y0-x2*y0-x0*y1+x2*y1+x0*y2-x1*y2;
	
    	if (fabs(znam) < 0.000001) 
	{
		printf("%10lf<- skipping\n%10lf %10lf %10lf \n%10lf %10lf %10lf \n",znam, x0,x1,x2,y0,y1,y2); return 0.0;		
	}
        a = 	(f2*x1*y0-f1*x2*y0-f2*x0*y1+f0*x2*y1+f1*x0*y2-f0*x1*y2)/znam;
        b = 	(y0*(f1-f2)+y1*(f2-f0)+y2*(f0-f1))/znam;
        c = 	-(f2*(x1-x0)+f1*(x0-x2)+f0*(x2-x1))/znam;
 
	return a+_crd[0]*b+_crd[1]*c;//+_crd[2]*d;
//	return 0.0;
}
double fsq(double* coeff_num, double* _crd)
{
	return coeff_num[3]*_crd[0]*_crd[0] + coeff_num[4]*_crd[1]*_crd[1] + coeff_num[5]*_crd[0]*_crd[1] + 
         		coeff_num[1]*_crd[0] + coeff_num[2]*_crd[1] + coeff_num[0];
}
double fsq(double* coeff_num, double x, double y)
{
	double _x[2]={x,y};
	return fsq(coeff_num,_x);
}

double fsqgx(double* coeff_num, double* _crd)
{
	return 2.0*coeff_num[3]*_crd[0] + coeff_num[5]*_crd[1] + coeff_num[1];
}
double fsqgx(double* coeff_num, double x, double y)
{
	double _x[2]={x,y};
	return fsqgx(coeff_num,_x);
}

double fsqgy(double* coeff_num, double* _crd)
{
	return 2.0*coeff_num[4]*_crd[1] + coeff_num[5]*_crd[0] + coeff_num[2];
}
double fsqgy(double* coeff_num, double x, double y)
{
	double _x[2]={x,y};
	return fsqgy(coeff_num,_x);
}
/*void Method::intoRandomAxes(double* x, double* y, int flag)
{
	if (flag)
		intoRandomAxes(x,x+1,y);
	else
	 	intoRandomAxes(x,y,axis);
}
void Method::intoRandomAxes(double* c)
{
	intoRandomAxes(c,c+1,0);
}*/

void Method::intoRandomAxes(double* x, double* y, double *axes)
{
	double c[2]={	(*x)*axes[0]+(*y)*axes[1],
			(*x)*axes[2]+(*y)*axes[3]};
	*x=c[0];*y=c[1];
}


void Method::intoRandomAxesGrad(double* x, double* y, double *axes)
{
	double c[2]={	(*x)*axes[0]+(*y)*axes[2],
			(*x)*axes[1]+(*y)*axes[3]};
	*x=c[0];*y=c[1];
}

void Method::fromRandomAxes(double* c, double *a)
{
	double det = a[0]*a[3]-a[1]*a[2];
	double axesR[4]={   a[3]/det, - a[1]/det,
			   -a[2]/det,   a[0]/det};
	
	intoRandomAxes(c,c+1,axesR);
}

void Method::fromRandomAxes(double* x, double* y, double *a)
{
	double det = a[0]*a[3]-a[1]*a[2];
	double axesR[4]={   a[3]/det, - a[1]/det,
			   -a[2]/det,   a[0]/det};
	
	intoRandomAxes(x,y,axesR);
}
void Method::fromRandomAxesGrad(double* x, double* y, double *a)
{
	double det = a[0]*a[3]-a[1]*a[2];
	double axesR[4]={   a[3]/det, - a[2]/det,
			   -a[1]/det,   a[0]/det};
	
	intoRandomAxes(x,y,axesR);
}
void swap(double* x, double* y)
{
	double tmp=*x; *x=*y; *y=tmp; 
}
double Method::interpolate_2_order(Triangle* t, double* _crd, Mesh* mesh, double* _res)//, Node* node)
{	
	Node* nodes[4] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2])};
	
	double 	coeff_num[6]={0.0},c[6]={
	 	nodes[0]->coords[0], nodes[1]->coords[0], nodes[2]->coords[0],
       		nodes[0]->coords[1], nodes[1]->coords[1], nodes[2]->coords[1]};
	double	x0=c[0], x1=c[1], x2=c[2], //
       		y0=c[3], y1=c[4], y2=c[5];
        double f0=nodes[0]->u[1],f1=nodes[1]->u[1],f2=nodes[2]->u[1],
		f0y=nodes[0]->vy,f1y=nodes[1]->vy,f2y=nodes[2]->vy; //solve for 4x4
	double axes[4]={1.0,0.0,0.0,1.0},crd[2]={_crd[0],_crd[1]};
	
		double 	l01 = sqrt(scalar(x0-x1,y0-y1,x0-x1,y0-y1)),
			l02 = sqrt(scalar(x0-x2,y0-y2,x0-x2,y0-y2)),
			l12 = sqrt(scalar(x2-x1,y2-y1,x2-x1,y2-y1)),
		//shift
			mid[2] = {(x0+x1+x2)/3.0,(y0+y1+y2)/3.0};
		x0 -= mid[0];		y0 -= mid[1];		
		x1 -= mid[0];		y1 -= mid[1];
		x2 -= mid[0];		y2 -= mid[1];
		crd[0] -= mid[0];	crd[1] -= mid[1];
				

		//rotate
		axes[2]=x2-(x0+x1)/2.0; axes[3]=y2-(y0+y1)/2.0; 
//		axes[2]*=10000000.0;
		double norma = sqrt(axes[2]*axes[2]+axes[3]*axes[3]);
		axes[2]/=norma; axes[3]/=norma; 	
		if (fabs(axes[2]) > 0.0)
		{
			axes[1] = 1.0;
			axes[0] = -axes[3]/axes[2];
		}
		else if (fabs(axes[3]) > 0.0)
		{
			axes[0] = 1.0;
			axes[1] = -axes[2]/axes[3];
		}
		else 
		{
			printf("Fail in determining axis direction");
			return 0.0;
		}
		norma = sqrt(axes[0]*axes[0]+axes[1]*axes[1]);
		axes[0]/=norma; axes[1]/=norma; 
		swap(axes,axes+2);swap(axes+1,axes+3);
printf("%lf %lf  %lf %lf  %lf %lf\n",x0,y0,x1,y1,x2,y2);
		intoRandomAxes(&x0,&y0,axes);
		intoRandomAxes(&x1,&y1,axes);
		intoRandomAxes(&x2,&y2,axes);
		intoRandomAxes(crd,crd+1,axes);
		fromRandomAxesGrad(&f0,&f0y,axes);
		fromRandomAxesGrad(&f1,&f1y,axes);
		fromRandomAxesGrad(&f2,&f2y,axes);
printf("%lf %lf  %lf %lf  %lf %lf\n\n",x0,y0,x1,y1,x2,y2);

		//size
		double norm = l01/1.0; 
		double l22 = sqrt(scalar(x2-(x0+x1)/2.0,y2-(y0+y1)/2.0,x2-(x0+x1)/2.0,y2-(y0+y1)/2.0));
		double norm_y = l22/1.0;
		x0 /= norm;		y0 /= norm_y;
		x1 /= norm;		y1 /= norm_y;
		x2 /= norm;		y2 /= norm_y;
		crd[0] /= norm;		crd[1] /= norm_y;
		f0*=norm; f0y*=norm_y;
		f1*=norm; f1y*=norm_y;
		f2*=norm; f2y*=norm_y;
	double znam = x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2;
    	if (fabs(znam) < 0.00000001) 
	{
		printf("%10lf<- skipping by volume\n%10lf %10lf %10lf \n%10lf %10lf %10lf \n",znam, x0,x1,x2,y0,y1,y2); return 0.0;		
	}
        coeff_num[1] = 	(f2*x1*y0-f1*x2*y0-f2*x0*y1+f0*x2*y1+f1*x0*y2-f0*x1*y2)/znam;
        coeff_num[3] = 	(y0*(f1-f2)+y1*(f2-f0)+y2*(f0-f1))/znam/2.0;
        coeff_num[5] = 	-(f2*(x1-x0)+f1*(x0-x2)+f0*(x2-x1))/znam;
        double 	f3=nodes[0]->u[0]-coeff_num[1]*x0-coeff_num[3]*x0*x0-coeff_num[5]*x0*y0,
		f4=nodes[1]->u[0]-coeff_num[1]*x1-coeff_num[3]*x1*x1-coeff_num[5]*x1*y1,
		f5=nodes[2]->u[0]-coeff_num[1]*x2-coeff_num[3]*x2*x2-coeff_num[5]*x2*y2; 
	znam = (y0-y1)*(y0-y2)*(y1-y2);
        coeff_num[0] = 	(f5*y0*(y0-y1)*y1+y2*(f3*y1*(y1-y2)+f4*y0*(y2-y0)))/znam;
        coeff_num[2] = 	(f5*(y1+y0)*(y1-y0)+f4*(y0-y2)*(y0+y2)+f3*(y2-y1)*(y2+y1))/znam;
        coeff_num[4] = 	(f5*(y0-y1)+f3*(y1-y2)+f4*(y2-y0))/znam;

	double coeff_num_0[6]={coeff_num[0],coeff_num[1],coeff_num[2],coeff_num[3],coeff_num[4],coeff_num[5]};
	double res_0 = fsq(coeff_num_0,crd),
		res_gx_0 = fsqgx(coeff_num_0,crd),
		res_gy_0 = fsqgy(coeff_num_0,crd);
	res_gx_0 /= norm; res_gy_0 /= norm;
	fromRandomAxes(&res_gx_0, &res_gy_0, axes);


		x0=c[0], x1=c[1], x2=c[2], //
       		y0=c[3], y1=c[4], y2=c[5];
        f0=nodes[0]->u[1],f1=nodes[1]->u[1],f2=nodes[2]->u[1],
		f0y=nodes[0]->vy,f1y=nodes[1]->vy,f2y=nodes[2]->vy; 
	crd[0]=_crd[0]; crd[1]=_crd[1];
	
		//shift
		x0 -= mid[0];		y0 -= mid[1];		
		x1 -= mid[0];		y1 -= mid[1];
		x2 -= mid[0];		y2 -= mid[1];
		crd[0] -= mid[0];	crd[1] -= mid[1];
				

		//rotate
		axes[2]=x1-(x0+x2)/2.0; axes[3]=y1-(y0+y2)/2.0; 
//		axes[2]*=10000000.0;
		norma = sqrt(axes[2]*axes[2]+axes[3]*axes[3]);
		axes[2]/=norma; axes[3]/=norma; 	
		if (fabs(axes[2]) > 0.0)
		{
			axes[1] = 1.0;
			axes[0] = -axes[3]/axes[2];
		}
		else if (fabs(axes[3]) > 0.0)
		{
			axes[0] = 1.0;
			axes[1] = -axes[2]/axes[3];
		}
		else 
		{
			printf("Fail in determining axis direction");
			return 0.0;
		}
		norma = sqrt(axes[0]*axes[0]+axes[1]*axes[1]);
		axes[0]/=norma; axes[1]/=norma; 
		swap(axes,axes+2);swap(axes+1,axes+3);
		intoRandomAxes(&x0,&y0,axes);
		intoRandomAxes(&x1,&y1,axes);
		intoRandomAxes(&x2,&y2,axes);
		intoRandomAxes(crd,crd+1,axes);
		fromRandomAxesGrad(&f0,&f0y,axes);
		fromRandomAxesGrad(&f1,&f1y,axes);
		fromRandomAxesGrad(&f2,&f2y,axes);

		//size
		double 	l11 = sqrt(scalar(x1-(x0+x2)/2.0,y1-(y0+y2)/2.0,x1-(x0+x2)/2.0,y1-(y0+y2)/2.0));
		norm_y = l11/1.0;
		norm = l02/1.0; 
		x0 /= norm;		y0 /= norm_y;
		x1 /= norm;		y1 /= norm_y;
		x2 /= norm;		y2 /= norm_y;
		crd[0] /= norm;		crd[1] /= norm_y;
		f0*=norm; f0y*=norm_y;
		f1*=norm; f1y*=norm_y;
		f2*=norm; f2y*=norm_y;
	znam = x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2;
    	if (fabs(znam) < 0.00000001) 
	{
		printf("%10lf<- skipping by volume\n%10lf %10lf %10lf \n%10lf %10lf %10lf \n",znam, x0,x1,x2,y0,y1,y2); return 0.0;		
	}
        coeff_num[1] = 	(f2*x1*y0-f1*x2*y0-f2*x0*y1+f0*x2*y1+f1*x0*y2-f0*x1*y2)/znam;
        coeff_num[3] = 	(y0*(f1-f2)+y1*(f2-f0)+y2*(f0-f1))/znam/2.0;
        coeff_num[5] = 	-(f2*(x1-x0)+f1*(x0-x2)+f0*(x2-x1))/znam;
        f3=nodes[0]->u[0]-coeff_num[1]*x0-coeff_num[3]*x0*x0-coeff_num[5]*x0*y0,
		f4=nodes[1]->u[0]-coeff_num[1]*x1-coeff_num[3]*x1*x1-coeff_num[5]*x1*y1,
		f5=nodes[2]->u[0]-coeff_num[1]*x2-coeff_num[3]*x2*x2-coeff_num[5]*x2*y2; 
	znam = (y0-y1)*(y0-y2)*(y1-y2);
        coeff_num[0] = 	(f5*y0*(y0-y1)*y1+y2*(f3*y1*(y1-y2)+f4*y0*(y2-y0)))/znam;
        coeff_num[2] = 	(f5*(y1+y0)*(y1-y0)+f4*(y0-y2)*(y0+y2)+f3*(y2-y1)*(y2+y1))/znam;
        coeff_num[4] = 	(f5*(y0-y1)+f3*(y1-y2)+f4*(y2-y0))/znam;

	double coeff_num_1[6]={coeff_num[0],coeff_num[1],coeff_num[2],coeff_num[3],coeff_num[4],coeff_num[5]};
	double res_1 = fsq(coeff_num_1,crd),
		res_gx_1 = fsqgx(coeff_num_1,crd),
		res_gy_1 = fsqgy(coeff_num_1,crd);
	res_gx_1 /= norm; res_gy_1 /= norm;
	fromRandomAxes(&res_gx_1, &res_gy_1, axes);


		x0=c[0], x1=c[1], x2=c[2], //
       		y0=c[3], y1=c[4], y2=c[5];
        f0=nodes[0]->u[1],f1=nodes[1]->u[1],f2=nodes[2]->u[1],
		f0y=nodes[0]->vy,f1y=nodes[1]->vy,f2y=nodes[2]->vy; 
	crd[0]=_crd[0]; crd[1]=_crd[1];
	
		//shift
		x0 -= mid[0];		y0 -= mid[1];		
		x1 -= mid[0];		y1 -= mid[1];
		x2 -= mid[0];		y2 -= mid[1];
		crd[0] -= mid[0];	crd[1] -= mid[1];
				

		//rotate
		axes[2]=x0-(x2+x1)/2.0; axes[3]=y0-(y2+y1)/2.0;  
//		axes[2]*=10000000.0;
		norma = sqrt(axes[2]*axes[2]+axes[3]*axes[3]);
		axes[2]/=norma; axes[3]/=norma; 	
		if (fabs(axes[2]) > 0.0)
		{
			axes[1] = 1.0;
			axes[0] = -axes[3]/axes[2];
		}
		else if (fabs(axes[3]) > 0.0)
		{
			axes[0] = 1.0;
			axes[1] = -axes[2]/axes[3];
		}
		else 
		{
			printf("Fail in determining axis direction");
			return 0.0;
		}
		norma = sqrt(axes[0]*axes[0]+axes[1]*axes[1]);
		axes[0]/=norma; axes[1]/=norma; 
		swap(axes,axes+2);swap(axes+1,axes+3);
		intoRandomAxes(&x0,&y0,axes);
		intoRandomAxes(&x1,&y1,axes);
		intoRandomAxes(&x2,&y2,axes);
		intoRandomAxes(crd,crd+1,axes);
		fromRandomAxesGrad(&f0,&f0y,axes);
		fromRandomAxesGrad(&f1,&f1y,axes);
		fromRandomAxesGrad(&f2,&f2y,axes);

		//size
		double 	l00 = sqrt(scalar(x0-(x1+x2)/2.0,y0-(y1+y2)/2.0,x0-(x1+x2)/2.0,y0-(y1+y2)/2.0));
		norm_y = l00/1.0;
		norm = l12/1.0; 
		x0 /= norm;		y0 /= norm_y;
		x1 /= norm;		y1 /= norm_y;
		x2 /= norm;		y2 /= norm_y;
		crd[0] /= norm;		crd[1] /= norm_y;
		f0*=norm; f0y*=norm_y;
		f1*=norm; f1y*=norm_y;
		f2*=norm; f2y*=norm_y;
	znam = x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2;
    	if (fabs(znam) < 0.00000001) 
	{
		printf("%10lf<- skipping by volume\n%10lf %10lf %10lf \n%10lf %10lf %10lf \n",znam, x0,x1,x2,y0,y1,y2); return 0.0;		
	}
        coeff_num[1] = 	(f2*x1*y0-f1*x2*y0-f2*x0*y1+f0*x2*y1+f1*x0*y2-f0*x1*y2)/znam;
        coeff_num[3] = 	(y0*(f1-f2)+y1*(f2-f0)+y2*(f0-f1))/znam/2.0;
        coeff_num[5] = 	-(f2*(x1-x0)+f1*(x0-x2)+f0*(x2-x1))/znam;
        f3=nodes[0]->u[0]-coeff_num[1]*x0-coeff_num[3]*x0*x0-coeff_num[5]*x0*y0,
		f4=nodes[1]->u[0]-coeff_num[1]*x1-coeff_num[3]*x1*x1-coeff_num[5]*x1*y1,
		f5=nodes[2]->u[0]-coeff_num[1]*x2-coeff_num[3]*x2*x2-coeff_num[5]*x2*y2; 
	znam = (y0-y1)*(y0-y2)*(y1-y2);
        coeff_num[0] = 	(f5*y0*(y0-y1)*y1+y2*(f3*y1*(y1-y2)+f4*y0*(y2-y0)))/znam;
        coeff_num[2] = 	(f5*(y1+y0)*(y1-y0)+f4*(y0-y2)*(y0+y2)+f3*(y2-y1)*(y2+y1))/znam;
        coeff_num[4] = 	(f5*(y0-y1)+f3*(y1-y2)+f4*(y2-y0))/znam;

	double coeff_num_2[6]={coeff_num[0],coeff_num[1],coeff_num[2],coeff_num[3],coeff_num[4],coeff_num[5]};
	double res_2 = fsq(coeff_num_2,crd),
		res_gx_2 = fsqgx(coeff_num_2,crd),
		res_gy_2 = fsqgy(coeff_num_2,crd);
	res_gx_2 /= norm; res_gy_2 /= norm;
	fromRandomAxes(&res_gx_2, &res_gy_2, axes);


	double minU = nodes[0]->u[0], maxU = nodes[0]->u[0];
	for (int i=1; i<3; i++)
	{
		if (nodes[i]->u[0] < minU) minU = nodes[i]->u[0];
		if (nodes[i]->u[0] > maxU) maxU = nodes[i]->u[0];
	}

	double res = 0.0, res_gx = (res_gx_0+res_gx_1+res_gx_2)/3.0, res_gy = (res_gy_0+res_gy_1+res_gy_2)/3.0;
	double diff_0=0,diff_1=0,diff_2=0;	
	int d0=0,d1=0,d2=0;
	if (res_0 < minU) { diff_0 += minU - res_0; d0=-1;};
	if (res_0 > maxU) { diff_0 += res_0 - maxU; d0=1;};
	if (res_1 < minU) { diff_1 += minU - res_1; d1=-1;};
	if (res_1 > maxU) { diff_1 += res_1 - maxU; d1=1;};
	if (res_2 < minU) { diff_2 += minU - res_2; d2=-1;};
	if (res_2 > maxU) { diff_2 += res_2 - maxU; d2=1;};
	
	int n=0;
	if (!d0) {res += res_0; n++;};
	if (!d1) {res += res_1; n++;};
	if (!d2) {res += res_2; n++;};
//	if (!d0 && !d1) {res = (res_0+res_1)/2.0; n++;}
//	else if (!d1 && !d2) {res = (res_1+res_2)/2.0; n++;}
//	else if (!d0 && !d2) {res = (res_0+res_2)/2.0; n++;}
//	if (n) res /= n;
//	else 
	res = (res_0+res_1+res_2)/3.0;
		//printf("Fall to 1st order\n");
		//res = interpolate_1_order(t, _crd, 0, mesh);
		res_gx = interpolate_1_order(t, _crd, 1, mesh);
		res_gy = interpolate_1_order(t, _crd, 2, mesh);
//	if (res < minU) res = minU;
//	if (res > maxU) res = maxU;
//	if (res!=res) res=0.0;
	_res[0] = res;
	_res[1] = res_gx;
	_res[2] = res_gy;
//	node->axis_method[0] = axes[0];
//	node->axis_method[1] = axes[1];
//	node->axis_method[2] = axes[2];
//	node->axis_method[3] = axes[3];
	//res = fsq(coeff_num,crd);
//res = (res_x + res_y)/2.0;
	
//	if (diff_x < diff_y)
//		res = res_x;
//	else res = res_y;
//	if (res != res) res = 0.0;
	//if (diff_x > 0.0001 && diff_y > 0.0001)
	//	printf("both axes failed %lf %lf\n",diff_x,diff_y);
//	if (d0 || d1 || d2)//(diff_x > 1.0 && diff_y > 1.0)//(res < minU || res > maxU)//res < minU-1.001)//(fabs(res) > fabs(10*(maxU+1))) 
	//if (res > 0.1)
//		printf(": res = %lf  res_0 = %lf  res_1 = %lf  res_2 = %lf   diff_0 = %lf  diff_1 = %lf  diff_2 = %lf\n"//%lf   u0 = %lf %lf %lf\nux = %lf %lf %lf\nuy = %lf %lf %lf\nf0 = %lf %lf %lf\nfx = %lf %lf %lf\nfy = %lf %lf %lf\nzn = %lf  vol = %lf\ncx = %lf %lf %lf \ncy = %lf %lf %lf \n\n",
//					,res, res_0,res_1, res_2, diff_0, diff_1,diff_2
		//			nodes[0]->u[0],nodes[1]->u[0],nodes[2]->u[0],
		//			nodes[0]->u[1],nodes[1]->u[1],nodes[2]->u[1],
		//			nodes[0]->u[2],nodes[1]->u[2],nodes[2]->u[2],
					//f0,f1,f2,f0y,f1y,f2y,
		//			fsq(coeff_num,x0,y0),fsq(coeff_num,x1,y1),fsq(coeff_num,x2,y2),
		//			f0,f1,f2,f0y,f1y,f2y,
					//fsqgx(coeff_num,x0,y0),fsqgx(coeff_num,x1,y1),fsqgx(coeff_num,x2,y2),
					//fsqgy(coeff_num,x0,y0),fsqgy(coeff_num,x1,y1),fsqgy(coeff_num,x2,y2),
				//coeff_num[0],coeff_num[1],coeff_num[2],coeff_num[3],coeff_num[4],coeff_num[5], 
		//		znam,x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2,
				//x0,x1,x2,y0,y1,y2);
		//		c[0],c[1],c[2],c[3],c[4],c[5]
//		);
	
//	if (res!=res) res=0.0;
	return res;
}

/*double fcb(double* a, double* c)
{
       return a[0] + a[1]*c[0] + a[2]*c[1] + a[3]*c[2] + a[4]*c[0]*c[0] + a[5]*c[1]*c[1] + a[6]*c[2]*c[2] + a[7]*c[0]*c[1] + a[8]*c[0]*c[2] + a[9]*c[1]*c[2] + a[10]*c[0]*c[0]*c[0] + a[11]*c[1]*c[1]*c[1] + a[12]*c[2]*c[2]*c[2] + a[13]*c[0]*c[0]*c[1] + a[14]*c[0]*c[0]*c[2] + a[15]*c[1]*c[1]*c[0] + a[16]*c[1]*c[1]*c[2] + a[17]*c[2]*c[2]*c[0] + a[18]*c[2]*c[2]*c[1] + a[19]*c[0]*c[1]*c[2];
                        
}

double fcb(double* a, double x, double y, double z)
{
       return a[0] + a[1]*x + a[2]*y + a[3]*z + a[4]*x*x + a[5]*y*y + a[6]*z*z + a[7]*x*y + a[8]*x*z + a[9]*y*z + a[10]*x*x*x + a[11]*y*y*y + a[12]*z*z*z + a[13]*x*x*y + a[14]*x*x*z + a[15]*y*y*x + a[16]*y*y*z + a[17]*z*z*x + a[18]*z*z*y + a[19]*x*y*z;
                        
}*/
double Method::interpolate_3_order(Triangle* t, double* _crd, Mesh* mesh)
{
	Node* nodes[4] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2])};
	//int min[3]={0},max[3]={0},diff=0;
	double 	c[6]={
	 	nodes[0]->coords[0], nodes[1]->coords[0], nodes[2]->coords[0],
       		nodes[0]->coords[1], nodes[1]->coords[1], nodes[2]->coords[1]};

		//printf("%10lf %10lf %10lf %10lf <- before rand\n%10lf %10lf %10lf %10lf\n%10lf %10lf %10lf %10lf\n",y0,y1,y2,y3, x0,x1,x2,x3,z0,z1,z2,z3); 
	/*double axes[9]; randomizeAxis(axes);
	intoRandomAxes(c,c+4,c+8,axes);
	intoRandomAxes(c+1,c+5,c+9,axes);
	intoRandomAxes(c+2,c+6,c+10,axes);
	intoRandomAxes(c+3,c+7,c+11,axes);
	int min=0;
	for (int i=1;i<4;i++) if (c[4+i]<c[4+min]) min=i;
	for (int i=0;i<3;i++) {double tmp=c[i*4+min]; c[i*4+min]=c[i*4]; c[i*4]=tmp;};
	int max=0;
	for (int i=1;i<4;i++) if (c[4+i]>c[4+max]) max=i;
	for (int i=0;i<3;i++) {double tmp=c[i*4+max]; c[i*4+max]=c[i*4+3]; c[i*4+3]=tmp;};
	//for (int i=0;i<12;i++) c[i]+=10.0;*/
	double	x0=c[0],x1=c[1],x2=c[2],x3=c[3],
		y0=c[4],y1=c[5],y2=c[6],y3=c[7],
		z0=c[8],z1=c[9],z2=c[10],z3=c[11];

       	double	coeff_num[20];
	return 0.0;
}


void Method::count(Mesh* mesh, Node* node, double timeStep, int ax)
{
	//randomizeAxis();
	//printf("count node %d axis %d start\n",node->local_num, ax);
	int v_n = 5;//number of value components 		sizeof(node->values.name)/sizeof(double);
	Node* next = new Node(node); 	//we store new time step node in here
	double nextValues[5]={0.0}; 	//(!)4=v_n new time step values
	double a_r_coeff[2]={0.0}; 	//coeffs in xi-eta-teta coords (axis[])
	double coord_char[2]={0.0}; 	//coordinates of point in old time step, where characteristic falls
	Triangle* t = 0;		//thetr for interpolation
	//for (int v_c=0; v_c<v_n; v_c++) 	//for each value component
	//{
	calculateCoeff(a_r_coeff,node->axis); 	//transform coefficients into random axes basis

		for (int i_crd=0; i_crd<2; i_crd++) //finding where characteristic for ax-th axis falls
			coord_char[i_crd] = node->coords[i_crd] - a_r_coeff[ax]*node->axis[2*ax+i_crd]*timeStep;
		//printf("COUNT:   ts %lf   c0 %lf   c1 %lf   ax %d %lf %lf\n",timeStep,a_r_coeff[0],a_r_coeff[1],ax,node->axis[2*ax],node->axis[2*ax+1]);
		//printf("COUNT:   ts %lf   nd: %lf %lf   cc: %lf %lf\n",timeStep,node->coords[0],node->coords[1],coord_char[0],coord_char[1]);
	//printf("count node %d axis %d tri to find at %lf %lf\n",node->local_num, ax,node->axis[2*ax],node->axis[2*ax+1]);
	//printf("axes randomized with %lf %lf %lf %lf\n",node->axis[0],node->axis[1],node->axis[2],node->axis[3]);
		t = mesh->findTriangle(coord_char,node);
	//printf("count node %d axis %d tri found\n",node->local_num, ax);
		if (!t) {printf("Fail! No thetr found for %lf %lf\n",coord_char[0],coord_char[1]); return;};
		if (order == 1)
		{
			nextValues[0] = interpolate_1_order(t, coord_char, 0, mesh); 
			nextValues[1] = interpolate_1_order(t, coord_char, 1, mesh); 
			nextValues[2] = interpolate_1_order(t, coord_char, 2, mesh);
		}
		else if (order == 2)
		{
			double res[3]={0.0};
//			nextValues[0] = 
			interpolate_2_order(t, coord_char, mesh, res); 
			nextValues[0] = res[0];
			nextValues[1] = interpolate_1_order(t, coord_char, 1, mesh); 
			nextValues[2] = interpolate_1_order(t, coord_char, 2, mesh);
		}
		else if (order == 3)
			nextValues[0] = interpolate_3_order(t, coord_char, mesh); 
	//printf("count node %d axis %d interpolated u with %lf\n",node->local_num, ax,nextValues[0]);
	//printf("count node %d axis %d interpolated ux with %lf\n",node->local_num, ax,nextValues[1]);
	//printf("count node %d axis %d interpolated uy with %lf\n",node->local_num, ax,nextValues[2]);
	next->setValues(nextValues);	//copy new time step values into new time step node
	//memcpy(next->coords, node->coords, sizeof(double)*3); //new time step coords, no mesh movement while
	node->nextStep = next;		//add link from old node to the new one
	//printf("count node %d axis %d finished\n",node->local_num, ax);
}
