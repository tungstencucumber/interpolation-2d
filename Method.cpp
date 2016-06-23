#include "Method.h"
#include "Triangle.h"
#include "Node.h"
#include "Mesh.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <stdio.h>
#include <iostream>
#define R_I 1000000
#define R_D 1000000.0
using std::cout;
using std::endl;

Method::Method()
{
    la = 1.0;
    mu = 0.7;
    rh = 1.0;
    c[0] = sqrt((la+2*mu)/rh);
    c[1] = sqrt(mu/rh);
    c[2] = 0;
    c[3] = - sqrt(mu/rh);
    c[4] = - sqrt((la+2*mu)/rh);
    om[0][0] = 1; 
    om[0][1] = 0;
    om[0][2] = 0;
    om[0][3] = 0;
    om[0][4] = - om[0][0];
    om[1][0] = 0;
    om[1][1] = 1;
    om[1][2] = 0;
    om[1][3] = 1;
    om[1][4] = 0;
    om[2][0] = sqrt(rh*(la+2*mu));
    om[2][1] = 0;
    om[2][2] = 0;
    om[2][3] = 0;
    om[2][4] = -om[2][0];
    om[3][0] = 0;
    om[3][1] = sqrt(rh*mu);
    om[3][2] = 0;
    om[3][3] = -om[3][1];
    om[3][4] = 0;
    om[4][0] = la*sqrt(rh/(la+2*mu));
    om[4][1] = 0;
    om[4][2] = 1;
    om[4][3] = 0;
    om[4][4] = -om[4][0];

    om_neg[0][0] = 0.5;
    om_neg[0][1] = 0;
    om_neg[0][2] = 1/(2*sqrt(rh*(la+2*mu)));
    om_neg[0][3] = 0;
    om_neg[0][4] = 0;
    om_neg[1][0] = 0;
    om_neg[1][1] = 0.5;
    om_neg[1][2] = 0;
    om_neg[1][3] = 1/(2*sqrt(rh*mu));
    om_neg[1][4] = 0;
    om_neg[2][0] = 0;
    om_neg[2][1] = 0;
    om_neg[2][2] = - la/(la+2*mu);
    om_neg[2][3] = 0;
    om_neg[2][4] = 1;
    om_neg[3][0] = 0;
    om_neg[3][1] = 0.5;
    om_neg[3][2] = 0;
    om_neg[3][3] = -1/(2*sqrt(rh*mu));
    om_neg[3][4] = 0;
    om_neg[4][0] = 0.5;
    om_neg[4][1] = 0;
    om_neg[4][2] = -1/(2*sqrt(rh*(la+2*mu)));
    om_neg[4][3] = 0;
    om_neg[4][4] = 0;
}

Method::~Method()
{
}

int Method::init()
{
	srand(time(0));
	//coeff[0] = coeff[1] = 0.2;
	//coeff[1] = 0.53;
	order = 1;
	return 0;
}

int Method::init(int _order)
{
	//srand(time(0));
	//coeff[0] = c0; coeff[1] = c1;
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
	norma = sqrt(axes[2]*axes[2]+axes[3]*axes[3]);
	axes[2]/=norma; axes[3]/=norma; 
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
	//_c[0] = coeff[0]*axes[0] + coeff[1]*axes[1];
	//_c[1] = coeff[0]*axes[2] + coeff[1]*axes[3];
}

double Method::interpolate_1_order(Triangle* t, double* _crd, int val, Mesh* mesh)
{	
	Node* nodes[3] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2])};
	double 	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], //
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1],   //
       		a,b,c;
    double f0 = 0, f1 = 0, f2 = 0;
    for (int i=0; i<5; i++)
    {
        f0 += om_neg[val][i]*nodes[0]->u[i];
        f1 += om_neg[val][i]*nodes[1]->u[i];
        f2 += om_neg[val][i]*nodes[2]->u[i];
    }

	double znam = x1*y0-x2*y0-x0*y1+x2*y1+x0*y2-x1*y2;
	
    	if (fabs(znam) < 0.000001) 
	{
		printf("%10lf<- skipping\n%10lf %10lf %10lf \n%10lf %10lf %10lf \n",znam, x0,x1,x2,y0,y1,y2); return 0.0;		
	}
        a = 	(f2*x1*y0-f1*x2*y0-f2*x0*y1+f0*x2*y1+f1*x0*y2-f0*x1*y2)/znam;
        b = 	(y0*(f1-f2)+y1*(f2-f0)+y2*(f0-f1))/znam;
        c = 	-(f2*(x1-x0)+f1*(x0-x2)+f0*(x2-x1))/znam;
 
	return a+_crd[0]*b+_crd[1]*c;
}
double fsq(double* a, double* c)
{
	return a[0] + a[1]*c[0] + a[2]*c[1] + a[3]*c[0]*c[1] + a[4]*c[0]*c[0] + a[5]*c[1]*c[1];
}
double fcb(double* a, double* c)
{
	return a[0] + a[1]*c[0] + a[2]*c[1] + a[3]*c[0]*c[1] + a[4]*c[0]*c[0] + a[5]*c[1]*c[1] + a[6]*c[0]*c[0]*c[1]+ a[7]*c[0]*c[1]*c[1]  +a[8]*c[0]*c[0]*c[0] + a[9]*c[1]*c[1]*c[1];
}
double fsq(double* coeff_num, double x, double y)
{
	double _x[2]={x,y};
	return fsq(coeff_num,_x);
}
double fcb(double* coeff_num, double x, double y)
{
	double _x[2]={x,y};
	return fcb(coeff_num,_x);
}
double fsqgx(double* a, double* c)
{
	return a[1] + a[3]*c[1] + 2.0*a[4]*c[0];
}
double fcbgx(double* a, double* c)
{
	return a[1] + a[3]*c[1] + 2.0*a[4]*c[0] + 2.0*a[6]*c[0]*c[1] + a[7]*c[1]*c[1] + 3.0*a[8]*c[0]*c[0];
}
double fsqgx(double* coeff_num, double x, double y)
{
	double _x[2]={x,y};
	return fsqgx(coeff_num,_x);
}

double fcbgx(double* coeff_num, double x, double y)
{
	double _x[2]={x,y};
	return fcbgx(coeff_num,_x);
}
double fsqgy(double* a, double* c)
{
	return a[2] + a[3]*c[0] + 2.0*a[5]*c[1];
}
double fcbgy(double* a, double* c)
{
	return a[2] + a[3]*c[0] + 2.0*a[5]*c[1] + a[6]*c[0]*c[0] + 2.0*a[7]*c[0]*c[1] + 3.0*a[9]*c[0]*c[0];
}

double fsqgy(double* coeff_num, double x, double y)
{
	double _x[2]={x,y};
	return fsqgy(coeff_num,_x);
}
double fcbgy(double* coeff_num, double x, double y)
{
	double _x[2]={x,y};
	return fcbgy(coeff_num,_x);
}

void Method::intoRandomAxes(double* x, double* y, double *axes)
{
	double c[2]={	(*x)*axes[0]+(*y)*axes[1],
			(*x)*axes[2]+(*y)*axes[3]};
	*x=c[0];*y=c[1];
}

void Method::intoAxes(Node* n, double* axis)
{
    if (!n) return;
    double vx = n->vx, vy = n->vy, sxx = n->sxx, sxy = n->sxy, syy = n->syy;
    n->vx = vx*axis[0] + vy*axis[1];
    n->vy = vx*axis[2] + vy*axis[3];
    n->sxx =    axis[0]*axis[0]*sxx + axis[0]*axis[2]*sxy + 
                axis[0]*axis[1]*sxy + axis[1]*axis[2]*syy;
    n->sxy =    axis[0]*axis[1]*sxx + axis[0]*axis[3]*sxy +
                axis[1]*axis[1]*sxy + axis[1]*axis[3]*syy;
    n->syy =    axis[1]*axis[2]*sxx + axis[2]*axis[3]*sxy +
                axis[1]*axis[3]*sxy + axis[3]*axis[3]*syy;
}


void Method::intoRandomAxesGrad(double* x, double* y, double *a)
{ 
	double _x = *x, _y = *y, det = a[2]*a[1] - a[0]*a[3];
	double 	new_x = _x*a[0] + _y*a[1],
		new_y = _x*a[2] + _y*a[3];
	*x = new_x; *y = new_y;
}


void Method::fromRandomAxes(double* c, double *a)
{	
	fromRandomAxes(c,c+1,a);
}

void Method::fromRandomAxes(double* x, double* y, double *a)
{
	double _x = *x, _y = *y, det = (a[0]*a[3] - a[2]*a[1]);
	double 	new_x = (a[3]*_x - _y*a[2])/det,
		new_y = (a[0]*_y - _x*a[1])/det;
	*x = new_x; *y = new_y;
}
void swap(double* x, double* y)
{
	double tmp=*x; *x=*y; *y=tmp; 
}
double Method::interpolate_2_order(Triangle* t, double* _crd, Mesh* mesh, double* _res)
{	
	Node* nodes[3] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2])};
	//no simplex, grads by OA, OB, OC
/*	double 	coeff_num[10]={0.0};
	double	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], x = _crd[0], 
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1], y = _crd[1];
        double 	f0 =nodes[0]->u[0],f1 =nodes[1]->u[0],f2 =nodes[2]->u[0],
		f0x=nodes[0]->u[1],f1x=nodes[1]->u[1],f2x=nodes[2]->u[1],
		f0y=nodes[0]->u[2],f1y=nodes[1]->u[2],f2y=nodes[2]->u[2];
	double	c01 = x0, c02 = y0, c03 = x0*y0, c04 = x0*x0, c05 = y0*y0,
		c11 = x1, c12 = y1, c13 = x1*y1, c14 = x1*x1, c15 = y1*y1,
		c21 = x2, c22 = y2, c23 = x2*y2, c24 = x2*x2, c25 = y2*y2,
		c31 = x-x0, c32 = y-y0, c33 = y0*(x-x0)+x0*(y-y0), c34 = 2*x0*(x-x0), c35 = 2*y0*(y-y0),
		c41 = x-x1, c42 = y-y1, c43 = y1*(x-x1)+x1*(y-y1), c44 = 2*x1*(x-x1), c45 = 2*y1*(y-y1),
		c51 = x-x2, c52 = y-y2, c53 = y2*(x-x2)+x2*(y-y2), c54 = 2*x2*(x-x2), c55 = 2*y2*(y-y2);
	double	b0 = f0, b1 = f1, b2 = f2, 	b3 = fx0*(x-x0)+fy0*(y-y0),
						b4 = fx1*(x-x1)+fy1*(y-y1),
						b5 = fx2*(x-x2)+fy2*(y-y2);
	coeff_num[0] = ;
	coeff_num[1] = ;
	coeff_num[2] = ;
	coeff_num[3] = ;
	coeff_num[4] = ;
	coeff_num[5] = ; */

	//no simplex, two lines
/*	double	xA=nodes[0]->coords[0], xB=nodes[1]->coords[0], xC=nodes[2]->coords[0], 
       		yA=nodes[0]->coords[1], yB=nodes[1]->coords[1], yC=nodes[2]->coords[1];
	double 	fA =nodes[0]->u[0],fB =nodes[1]->u[0],fC =nodes[2]->u[0],
		fAx=nodes[0]->u[1],fBx=nodes[1]->u[1],fCx=nodes[2]->u[1],
		fAy=nodes[0]->u[2],fBy=nodes[1]->u[2],fCy=nodes[2]->u[2];
	double x0 = _crd[0], y0 = _crd[1];

	double eps = 0.000001;
	if (fabs(x0 - xA) < eps && fabs(y0 - yA) < eps)
		return fA;
	else if (fabs(x0 - xB) < eps && fabs(y0 - yB) < eps)
		return fB;
	else if (fabs(x0 - xC) < eps && fabs(y0 - yC) < eps)
		return fC;
	double 	t2 = ((xB-xA)*(yC-yA)-(xC-xA)*(yB-yA))/((y0-yB)*(xC-xA)-(yC-yA)*(x0-xB));
	double 	t1 = 0.0;
	if (fabs(xC-xA) > eps) t1 = t2*(x0-xB)/(xC-xA) + (xB-xA)/(xC-xA);
	else t1 = t2*(y0-yB)/(yC-yA) + (yB-yA)/(yC-yA);
	double 	fAt = fAx*(xC-xA) + fAy*(yC-yA),
		fCt = fCx*(xC-xA) + fCy*(yC-yA);
	double a[4] = {0};
	a[0] = fA;
	a[1] = fAt;
	a[2] = 3*fC - 3*fA - 2*fAt - fCt;
	a[3] = fCt + fAt - 2*fC + 2*fA;
	double 	fM = a[0] + a[1]*t1 + a[2]*t1*t1 + a[3]*t1*t1*t1,
		fMx = (fCx - fAx)*t1 + fAx,
		fMy = (fCy - fAy)*t1 + fAy,
		xM = xA + (xC - xA)*t1,
		yM = yA + (yC - yA)*t1;
	double	fBt = fBx*(xM-xB) + fAy*(yM-yB),
		fMt = fMx*(xM-xB) + fMy*(yM-yB);
	a[0] = fB;
	a[1] = fBt;
	a[2] = 3*fM - 3*fB - 2*fBt - fMt;
	a[3] = fMt + fBt - 2*fM + 2*fB;
	double	tt = 0;
	if (fabs(xM-xB) > eps) tt = (x0 - xB)/(xM - xB);
	else tt = (y0 - yB)/(yM - yB);
	return a[0] + a[1]*tt + a[2]*tt*tt + a[3]*tt*tt*tt;*/


/*	//to a simplex
	double 	coeff_num[10]={0.0};
	double	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], 
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1];
        double 	f0 =nodes[0]->u[0],f1 =nodes[1]->u[0],f2 =nodes[2]->u[0],
		f0x=nodes[0]->u[1],f1x=nodes[1]->u[1],f2x=nodes[2]->u[1],
		f0y=nodes[0]->u[2],f1y=nodes[1]->u[2],f2y=nodes[2]->u[2];

	double axes[4]={1.0,0.0,0.0,1.0},crd[2]={_crd[0],_crd[1]};
	
		//shift
	double	mid[2] = {x0,y0};
	x0 -= mid[0];		y0 -= mid[1];		
	x1 -= mid[0];		y1 -= mid[1];
	x2 -= mid[0];		y2 -= mid[1];
	crd[0] -= mid[0];	crd[1] -= mid[1];

	axes[0] = x1; axes[1] = y1;
	axes[2] = x2; axes[3] = y2;
	fromRandomAxes(&x0,&y0,axes);
	fromRandomAxes(&x1,&y1,axes);
	fromRandomAxes(&x2,&y2,axes);
	fromRandomAxes(crd,crd+1,axes);
	intoRandomAxesGrad(&f0x,&f0y,axes);
	intoRandomAxesGrad(&f1x,&f1y,axes);
	intoRandomAxesGrad(&f2x,&f2y,axes);
//----------------quad (0x, 1x, 2x) polynom coeffs
	double	x = crd[0], y = crd[1];
	if (x == 0)
	{
		coeff_num[0] = f0;
		coeff_num[1] = f0y;
		coeff_num[2] = 3*f2 - 3*f0 - 2*f0y - f2y;
		coeff_num[3] = f2y + f0y - 2*f2 + 2*f0;
		return coeff_num[0] + coeff_num[1]*y + coeff_num[2]*y*y + coeff_num[3]*y*y*y;
	}
	if (y == 0)
	{
		coeff_num[0] = f0;
		coeff_num[1] = f0x;
		coeff_num[2] = 3*f1 - 3*f0 - 2*f0x - f1x;
		coeff_num[3] = f1x + f0x - 2*f1 + 2*f0;
		return coeff_num[0] + coeff_num[1]*x + coeff_num[2]*x*x + coeff_num[3]*x*x*x;
	}
	if (x + y - 1 == 0)
	{
		Node* nodes[3] = {mesh->getNode(t->vert[2]), mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1])};
		//to a simplex
		double 	coeff_num[10]={0.0};
		double	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], 
	       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1];
		double 	f0 =nodes[0]->u[0],f1 =nodes[1]->u[0],f2 =nodes[2]->u[0],
			f0x=nodes[0]->u[1],f1x=nodes[1]->u[1],f2x=nodes[2]->u[1],
			f0y=nodes[0]->u[2],f1y=nodes[1]->u[2],f2y=nodes[2]->u[2];

		double axes[4]={1.0,0.0,0.0,1.0},crd[2]={_crd[0],_crd[1]};
	
			//shift
		double	mid[2] = {x0,y0};
		x0 -= mid[0];		y0 -= mid[1];		
		x1 -= mid[0];		y1 -= mid[1];
		x2 -= mid[0];		y2 -= mid[1];
		crd[0] -= mid[0];	crd[1] -= mid[1];

		axes[0] = x1; axes[1] = y1;
		axes[2] = x2; axes[3] = y2;
		fromRandomAxes(&x0,&y0,axes);
		fromRandomAxes(&x1,&y1,axes);
		fromRandomAxes(&x2,&y2,axes);
		fromRandomAxes(crd,crd+1,axes);
		intoRandomAxesGrad(&f0x,&f0y,axes);
		intoRandomAxesGrad(&f1x,&f1y,axes);
		intoRandomAxesGrad(&f2x,&f2y,axes);
		double	x = crd[0], y = crd[1];
		coeff_num[0] = f0;
		coeff_num[1] = f0y;
		coeff_num[2] = 3*f2 - 3*f0 - 2*f0y - f2y;
		coeff_num[3] = f2y + f0y - 2*f2 + 2*f0;
		return coeff_num[0] + coeff_num[1]*y + coeff_num[2]*y*y + coeff_num[3]*y*y*y;
	}
	double	A = f1 - f0,
		B = f2 - f0,
		C = x*f0x + y*f0y,
		D = (x-1)*f1x + y*f1y,
		E = x*f2x + (y-1)*f2y;
	double	F = C - A*x - B*y,
		G = D - A*(x-1) - B*y,
		H = E - A*x - B*(y-1);
	double	I = -(F*(x+y-1) - H*y + G*x)/(2*y*(x+y-1));
	double	J = (H - F - 2*I*y + I)/x,
		K = -(F + I*y)/x;
	coeff_num[0] = f0;
	coeff_num[1] = A - K;
	coeff_num[2] = B - I;
	coeff_num[3] = J;
	coeff_num[4] = K;
	coeff_num[5] = I;
	return fsq(coeff_num,crd); */

//----------------quad polynom coeffs

/*	coeff_num[4] = (f1x-f0x)/2;
	coeff_num[5] = (f2y-f0y)/2;
	coeff_num[3] = ((f2x-f0x) + (f1y-f0y))/2;
	coeff_num[1] = f0x;
	coeff_num[2] = f0y;
	coeff_num[0] = (f0 + (f1-f0x-coeff_num[4]) + (f2-f0y-coeff_num[5]))/3;

	return fsq(coeff_num,crd); */

//----------------cubic polynom coeffs

/*	coeff_num[0] = f0;
	coeff_num[1] = f0x;
	coeff_num[2] = f0y;
	coeff_num[3] = f1y - f0y;
	coeff_num[4] = 3*f1 - 3*f0 - f1x - 2*f0x;
	coeff_num[5] = 3*f2 - 3*f0 - f2y - 2*f0y;
	coeff_num[6] = 0;
	coeff_num[7] = f2x - f0x - f1y + f0y;
	coeff_num[8] = 2*f0 - 2*f1 + f1x + f0x;
	coeff_num[9] = 2*f0 - 2*f2 + f2y + f0y; 

	return fcb(coeff_num,crd);
*/

/*--------------- minmax limiter
	double minU = nodes[0]->u[0], maxU = nodes[0]->u[0];
	for (int i=1; i<3; i++)
	{
		if (nodes[i]->u[0] < minU) minU = nodes[i]->u[0];
		if (nodes[i]->u[0] > maxU) maxU = nodes[i]->u[0];
	} */
	
//----------------Hermite polynom
	double	x1=nodes[0]->coords[0], x2=nodes[1]->coords[0], x3=nodes[2]->coords[0], 
       		y1=nodes[0]->coords[1], y2=nodes[1]->coords[1], y3=nodes[2]->coords[1];
        double 	f1 =nodes[0]->u[0],f2 =nodes[1]->u[0],f3 =nodes[2]->u[0],
		f1x=nodes[0]->u[1],f2x=nodes[1]->u[1],f3x=nodes[2]->u[1],
		f1y=nodes[0]->u[2],f2y=nodes[1]->u[2],f3y=nodes[2]->u[2],
		x = _crd[0], y = _crd[1];
	double	a = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)),
		b = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3)),
		c = sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2));
	double	p = (a + b + c)/2;
	double	s = sqrt(p*(p-a)*(p-b)*(p-c));
	double	l1 = fabs(((x2*y3 - x3*y2) + (y2-y3)*x + (x3-x2)*y)/(2*s)),
		l2 = fabs(((x3*y1 - x1*y3) + (y3-y1)*x + (x1-x3)*y)/(2*s)),
		l3 = fabs(((x1*y2 - x2*y1) + (y1-y2)*x + (x2-x1)*y)/(2*s));
	double	alpha1 = l1*l1*l1 + 3*l2*l1*l1 + 3*l3*l1*l1 + 3*l1*l2*l3,
		alpha2 = l2*l2*l2 + 3*l1*l2*l2 + 3*l3*l2*l2 + 3*l1*l2*l3,
		alpha3 = l3*l3*l3 + 3*l1*l3*l3 + 3*l2*l3*l3 + 3*l1*l2*l3;
	double	beta1 = (x2-x1)*(l1*l1*l2 + 0.5*l1*l2*l3) + (x3-x1)*(l1*l1*l3 + 0.5*l1*l2*l3),	
		beta2 = (x1-x2)*(l2*l2*l1 + 0.5*l1*l2*l3) + (x3-x2)*(l2*l2*l3 + 0.5*l1*l2*l3),
		beta3 = (x1-x3)*(l3*l3*l1 + 0.5*l1*l2*l3) + (x2-x3)*(l3*l3*l2 + 0.5*l1*l2*l3);
	double	gamma1 = (y2-y1)*(l1*l1*l2 + 0.5*l1*l2*l3) + (y3-y1)*(l1*l1*l3 + 0.5*l1*l2*l3),
		gamma2 = (y1-y2)*(l2*l2*l1 + 0.5*l1*l2*l3) + (y3-y2)*(l2*l2*l3 + 0.5*l1*l2*l3),
		gamma3 = (y1-y3)*(l3*l3*l1 + 0.5*l1*l2*l3) + (y2-y3)*(l3*l3*l2 + 0.5*l1*l2*l3); 
/*	double	alpha1 = l1*l1*l1 + 3*l1*l1*(l2+l3) - 7*l1*l2*l3,
		alpha2 = l2*l2*l2 + 3*l2*l2*(l1+l3) - 7*l1*l2*l3,
		alpha3 = l3*l3*l3 + 3*l3*l3*(l1+l2) - 7*l1*l2*l3;
	double	beta1 = (x2-x1)*(l1*l1*l2 - l1*l2*l3) + (x3-x1)*(l1*l1*l3 - l1*l2*l3),	
		beta2 = (x3-x2)*(l2*l2*l1 - l1*l2*l3) + (x1-x2)*(l2*l2*l3 - l1*l2*l3),
		beta3 = (x1-x3)*(l3*l3*l1 - l1*l2*l3) + (x2-x3)*(l3*l3*l2 - l1*l2*l3);
	double	gamma1 = (y2-y1)*(l1*l1*l2 - l1*l2*l3) + (y3-y1)*(l1*l1*l3 - l1*l2*l3),
		gamma2 = (y3-y2)*(l2*l2*l1 - l1*l2*l3) + (y1-y2)*(l2*l2*l3 - l1*l2*l3),
		gamma3 = (y1-y3)*(l3*l3*l1 - l1*l2*l3) + (y2-y3)*(l3*l3*l2 - l1*l2*l3); */
	double 	res = 	(alpha1*f1 + beta1*f1x + gamma1*f1y) +
			(alpha2*f2 + beta2*f2x + gamma2*f2y) + 
			(alpha3*f3 + beta3*f3x + gamma3*f3y);
//	printf("%lf %lf   ", f3, res);
	return res;
}

double Method::interpolate_3_order(Triangle* t, double* _crd, Mesh* mesh)
{
	return 0.0;
}


void Method::count_split(Mesh* mesh, Node* node, double timeStep)
{
	//randomizeAxis();
	int v_n = 5;
	Node* next = new Node(node); 	//we store new time step node in here
	//node->nextStep = next;
	double nextValues[5]={0.0}; 	//(!)4=v_n new time step values
	//double a_r_coeff[2]={0.0}; 	//coeffs in xi-eta-teta coords (axis[])
	double coord_char[2]={0.0}; 	//coordinates of point in old time step, where characteristic falls
    double riem;
	Triangle* t = 0;		//thetr for interpolation

    intoAxes(node, node->axis);

    for (int ci = 0; ci<5; ci++)
    {
		for (int i_crd=0; i_crd<2; i_crd++) 
			coord_char[i_crd] = node->coords[i_crd] - c[ci]*node->axis[i_crd]*timeStep;
		t = mesh->findTriangle(coord_char,node);
        if (!t) {printf("Fail! No thetr found for %lf %lf\n",coord_char[0],coord_char[1]); return;};

        for (int i=0; i<3; i++)
        {
            Node* tn = mesh->getNode(t->vert[i]);
            intoAxes(tn, node->axis);
        }
	
        riem = interpolate_1_order(t, coord_char, ci, mesh);
        //if (riem == riem ) cout <<" r "<<riem;
        for (int i=0; i<5; i++)
		    nextValues[i] += om[i][ci]*riem; 

        for (int i=0; i<3; i++)
        {
            Node* tn = mesh->getNode(t->vert[i]);
            intoAxes(tn, node->axis_neg);
        }
    }

    intoAxes(node, node->axis_neg);

	next->setValues(nextValues);	//copy new time step values into new time step node
    intoAxes(next, node->axis_neg);

	node->nextStep = next;		//add link from old node to the new one
}
