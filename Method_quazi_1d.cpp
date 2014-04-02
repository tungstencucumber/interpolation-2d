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
	coeff[0] = coeff[1] = 0;
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
	double res = 0.0;
{
	Node* nodes[3] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2])};
	double 	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], //
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1]   //
       		;
                                                   
        double f0=nodes[0]->u[val],f1=nodes[1]->u[val],f2=nodes[2]->u[val]; //solve for 4x4

	double Ax = x0, Ay = y0, Bx = x1, By = y1, Cx = x2, Cy = y2, Ox = _crd[0], Oy = _crd[1],
		ax = Cx - Ax, ay = Cy - Ay,
		bx = Ax, by = Ay,
		cx = Ox - Bx, cy = Oy - By,
		dx = Bx, dy = By, tMAC,
		fA = f0, fB = f1, fC = f2, fO, fM, xM, yM, tO;

	if (cx == 0.0)
	{
		tMAC = (dx-bx)/ax;
	}
	else if (cy == 0.0)
	{
		tMAC = (dy-by)/ay;
	}
	else
	{
		tMAC = (dx*cy - dy*cx - bx*cy + by*cx)/(ax*cy - ay*cx);
	}
	/*double znam = ay*cx - ax*cy;
	if (fabs(znam) <= 0.000000001)
	{
		printf ("ZERO znam ");
		printf("A: %lf %lf  B: %lf %lf  C: %lf %lf  O: %lf %lf  M: %lf %lf  tMAC: %lf\n",Ax,Ay,Bx,By,Cx,Cy,Ox,Oy,xM,yM,tMAC);
		return 0;		
	}
	if (fabs(cx) > 0.00001)
		tMAC = (cx*by-cx*dy+dx*cy)/(-znam);
	else if (fabs(cy) > 0.00001)
		tMAC = (bx*cy-dx*cy+cx*dy)/znam;
	else {
		//printf ("ZERO cx cy ");
		//printf("A: %lf %lf  B: %lf %lf  C: %lf %lf  O: %lf %lf  M: %lf %lf  tMAC: %lf\n",Ax,Ay,Bx,By,Cx,Cy,Ox,Oy,xM,yM,tMAC);
		return fB;		
	}*/

	xM = ax*tMAC + bx;
	yM = ay*tMAC + by;

	double a = fC - fA, b = fA;
	fM = a*tMAC + b;

	double ex = xM - Bx, ey = yM - By, fx = Bx, fy = By;
	if (fabs(ex) > 0.000000001)
		tO = (Ox - fx)/ex;
	else if (fabs(ey) > 0.000000001)
		tO = (Oy - fy)/ey;
	else 
	{
		printf ("ZERO ex ey ");
		return 0;
	}
	a = fM - fB; b = fB;
	fO = a*tO + b;
	res += f0;
//	printf("B: %lf %lf  M: %lf %lf  O: %lf %lf  tO: %lf\n",Bx,By,xM,yM,Ox,Oy,tO);
	printf("A: %10lf  B: %10lf  C: %10lf  tMAC: %lf  fM: %10lf  tO: %lf  fO: %10lf\n",fA,fB,fC,tMAC,fM,tO,fO);
//	printf("A: %lf %lf  B: %lf %lf  C: %lf %lf  O: %lf %lf  M: %lf %lf  tMAC: %lf  tO: %lf\n",Ax,Ay,Bx,By,Cx,Cy,Ox,Oy,xM,yM,tMAC,tO);
}
	return res;
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
void swap(double* x, double* y)
{
	double tmp=*x; *x=*y; *y=tmp; 
}
double Method::interpolate_2_order(Triangle* t, double* _crd, Mesh* mesh, double* _res, Node* node)
{
	_res[0] = _res[1] = _res[2] = 0.0;
	double minU, maxU;

{	
	Node* nodes[3] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2])};
	double 	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], //
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1]   //
       		;
                                                   
        double f0=nodes[0]->u[0],f1=nodes[1]->u[0],f2=nodes[2]->u[0]; //solve for 4x4

	double Ax = x0, Ay = y0, Bx = x1, By = y1, Cx = x2, Cy = y2, Ox = _crd[0], Oy = _crd[1],
		a0x = Cx - Ax, a0y = Cy - Ay,
		b0x = Ax, b0y = Ay,
		a1x = Ox - Bx, a1y = Oy - By,
		b1x = Bx, b1y = By, tM0, tM, xM,yM,tO,
		fA = f0, fB = f1, fC = f2, fO, fM,
		fxA=nodes[0]->u[1],fxB=nodes[1]->u[1],fxC=nodes[2]->u[1],
		fyA=nodes[0]->u[2],fyB=nodes[1]->u[2],fyC=nodes[2]->u[2], fxM, fyM, fxO, fyO, c, d;
	double znam = a0x*a1y - a0y*a1x;
	if ( fabs(a1x) <= 0 )
	{
		tM0 = (b1x - b0x)/a0x;
		tM = (Ox - b0x)/a0x;
	} else if ( fabs(a1y) <= 0 )
	{
		tM0 = (b1y - b0y)/a0y;
		tM = (Oy - b0y)/a0y;
	} else if ( fabs(znam) > 0 )
	{
		tM0 = ( a1x*(b0y-b1y) - a1y*(b0x-b1x) )/znam;
		tM = (Ox - b0x)/a0x;
	} else 
	{
		printf ("ZERO znam!!!\n");
		return 0;
	}

	xM = a0x*tM0 + b0x;
	yM = a0y*tM0 + b0y;

	double a = fxC - fxA, b = fxA;
	fxM = a*tM + b;
	a = fyC - fyA, b = fyA;
	fyM = a*tM + b;

	double dfA = scalar(fxA, fyA, Cx - Ax, Cy - Ay),
		dfC = scalar(fxC, fyC, Cx - Ax, Cy - Ay);

	a = -2.0*fC + dfA + 2.0*fA + dfC;
	b = 3.0*fC - 2.0*dfA - 3.0*fA - dfC;
	c = dfA;
	d = fA;
	fM = a*tM*tM*tM + b*tM*tM + c*tM + d;

	if (fabs(xM - Bx) > 0)
		tO = (Ox - Bx)/(xM - Bx);
	else if (fabs(yM - By) > 0)
		tO = (Oy - By)/(yM - By);
	else 
	{
		printf ("ZERO znam!!!\n");
		return 0;
	}

	a = fxM - fxB, b = fxB;
	fxO = a*tO + b;
	a = fyM - fyB, b = fyB;
	fyO = a*tO + b;

	double dfB = scalar(fxB, fyB, xM - Bx, yM - By),
		dfM = scalar(fxM, fyM, xM - Bx, yM - By);

	a = -2.0*fM + dfB + 2.0*fB + dfM;
	b = 3.0*fM - 2.0*dfB - 3.0*fB - dfM;
	c = dfB;
	d = fB;
	fO = a*tO*tO*tO + b*tO*tO + c*tO + d;
	_res[0] += fO; _res[1] += fxO; _res[2] += fyO;
}
if (false)
{
	Node* nodes[3] = {mesh->getNode(t->vert[2]), mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1])};
	double 	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], //
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1]   //
       		;
                                                   
        double f0=nodes[0]->u[0],f1=nodes[1]->u[0],f2=nodes[2]->u[0]; //solve for 4x4

	double Ax = x0, Ay = y0, Bx = x1, By = y1, Cx = x2, Cy = y2, Ox = _crd[0], Oy = _crd[1],
		a0x = Cx - Ax, a0y = Cy - Ay,
		b0x = Ax, b0y = Ay,
		a1x = Ox - Bx, a1y = Oy - By,
		b1x = Bx, b1y = By, tM0, tM, xM,yM,tO,
		fA = f0, fB = f1, fC = f2, fO, fM,
		fxA=nodes[0]->u[1],fxB=nodes[1]->u[1],fxC=nodes[2]->u[1],
		fyA=nodes[0]->u[2],fyB=nodes[1]->u[2],fyC=nodes[2]->u[2], fxM, fyM, fxO, fyO, c, d;
	double znam = a0x*a1y - a0y*a1x;
	if ( fabs(a1x) <= 0 )
	{
		tM0 = (b1x - b0x)/a0x;
		tM = (Ox - b0x)/a0x;
	} else if ( fabs(a1y) <= 0 )
	{
		tM0 = (b1y - b0y)/a0y;
		tM = (Oy - b0y)/a0y;
	} else if ( fabs(znam) > 0 )
	{
		tM0 = ( a1x*(b0y-b1y) - a1y*(b0x-b1x) )/znam;
		tM = (Ox - b0x)/a0x;
	} else 
	{
		printf ("ZERO znam!!!\n");
		return 0;
	}

	xM = a0x*tM0 + b0x;
	yM = a0y*tM0 + b0y;

	double a = fxC - fxA, b = fxA;
	fxM = a*tM + b;
	a = fyC - fyA, b = fyA;
	fyM = a*tM + b;

	double dfA = scalar(fxA, fyA, Cx - Ax, Cy - Ay),
		dfC = scalar(fxC, fyC, Cx - Ax, Cy - Ay);

	a = -2.0*fC + dfA + 2.0*fA + dfC;
	b = 3.0*fC - 2.0*dfA - 3.0*fA - dfC;
	c = dfA;
	d = fA;
	fM = a*tM*tM*tM + b*tM*tM + c*tM + d;

	if (fabs(xM - Bx) > 0)
		tO = (Ox - Bx)/(xM - Bx);
	else if (fabs(yM - By) > 0)
		tO = (Oy - By)/(yM - By);
	else 
	{
		printf ("ZERO znam!!!\n");
		return 0;
	}

	a = fxM - fxB, b = fxB;
	fxO = a*tO + b;
	a = fyM - fyB, b = fyB;
	fyO = a*tO + b;

	double dfB = scalar(fxB, fyB, xM - Bx, yM - By),
		dfM = scalar(fxM, fyM, xM - Bx, yM - By);

	a = -2.0*fM + dfB + 2.0*fB + dfM;
	b = 3.0*fM - 2.0*dfB - 3.0*fB - dfM;
	c = dfB;
	d = fB;
	fO = a*tO*tO*tO + b*tO*tO + c*tO + d;
	_res[0] += fO; _res[1] += fxO; _res[2] += fyO;
}
if (false)
{
	Node* nodes[3] = {mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2]), mesh->getNode(t->vert[0])};
	double 	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], //
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1]   //
       		;
	minU = nodes[0]->u[0], maxU = nodes[0]->u[0];
	for (int i=1; i<3; i++)
	{
		if (nodes[i]->u[0] < minU) minU = nodes[i]->u[0];
		if (nodes[i]->u[0] > maxU) maxU = nodes[i]->u[0];
	}
                                                   
        double f0=nodes[0]->u[0],f1=nodes[1]->u[0],f2=nodes[2]->u[0]; //solve for 4x4

	double Ax = x0, Ay = y0, Bx = x1, By = y1, Cx = x2, Cy = y2, Ox = _crd[0], Oy = _crd[1],
		a0x = Cx - Ax, a0y = Cy - Ay,
		b0x = Ax, b0y = Ay,
		a1x = Ox - Bx, a1y = Oy - By,
		b1x = Bx, b1y = By, tM0, tM, xM,yM,tO,
		fA = f0, fB = f1, fC = f2, fO, fM,
		fxA=nodes[0]->u[1],fxB=nodes[1]->u[1],fxC=nodes[2]->u[1],
		fyA=nodes[0]->u[2],fyB=nodes[1]->u[2],fyC=nodes[2]->u[2], fxM, fyM, fxO, fyO, c, d;
	double znam = a0x*a1y - a0y*a1x;
	if ( fabs(a1x) <= 0 )
	{
		tM0 = (b1x - b0x)/a0x;
		tM = (Ox - b0x)/a0x;
	} else if ( fabs(a1y) <= 0 )
	{
		tM0 = (b1y - b0y)/a0y;
		tM = (Oy - b0y)/a0y;
	} else if ( fabs(znam) > 0 )
	{
		tM0 = ( a1x*(b0y-b1y) - a1y*(b0x-b1x) )/znam;
		tM = (Ox - b0x)/a0x;
	} else 
	{
		printf ("ZERO znam!!!\n");
		return 0;
	}

	xM = a0x*tM0 + b0x;
	yM = a0y*tM0 + b0y;

	double a = fxC - fxA, b = fxA;
	fxM = a*tM + b;
	a = fyC - fyA, b = fyA;
	fyM = a*tM + b;

	double dfA = scalar(fxA, fyA, Cx - Ax, Cy - Ay),
		dfC = scalar(fxC, fyC, Cx - Ax, Cy - Ay);

	a = -2.0*fC + dfA + 2.0*fA + dfC;
	b = 3.0*fC - 2.0*dfA - 3.0*fA - dfC;
	c = dfA;
	d = fA;
	fM = a*tM*tM*tM + b*tM*tM + c*tM + d;

	if (fabs(xM - Bx) > 0)
		tO = (Ox - Bx)/(xM - Bx);
	else if (fabs(yM - By) > 0)
		tO = (Oy - By)/(yM - By);
	else 
	{
		printf ("ZERO znam!!!\n");
		return 0;
	}

	a = fxM - fxB, b = fxB;
	fxO = a*tO + b;
	a = fyM - fyB, b = fyB;
	fyO = a*tO + b;

	double dfB = scalar(fxB, fyB, xM - Bx, yM - By),
		dfM = scalar(fxM, fyM, xM - Bx, yM - By);

	a = -2.0*fM + dfB + 2.0*fB + dfM;
	b = 3.0*fM - 2.0*dfB - 3.0*fB - dfM;
	c = dfB;
	d = fB;
	fO = a*tO*tO*tO + b*tO*tO + c*tO + d;
	_res[0] += fO; _res[1] += fxO; _res[2] += fyO;
}
	_res[0] /= 1.0; _res[1] /= 1.0; _res[2] /= 1.0;
//	if (_res[0]!=_res[0]) _res[0] = 0.0;
//	if (_res[0] < minU) _res[0] = minU;
//	if (_res[0] > maxU) _res[0] = maxU;
	return _res[0];
//	double minU = nodes[0]->u[0], maxU = nodes[0]->u[0];
//	for (int i=1; i<3; i++)
//	{
//		if (nodes[i]->u[0] < minU) minU = nodes[i]->u[0];
//		if (nodes[i]->u[0] > maxU) maxU = nodes[i]->u[0];
//	}

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
			coord_char[i_crd] = node->coords[i_crd] - 2.0*a_r_coeff[ax]*node->axis[2*ax+i_crd]*timeStep;
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
			interpolate_2_order(t, coord_char, mesh, res, next); 
			nextValues[0] = res[0];
			nextValues[1] = res[1];
			nextValues[2] = res[2];
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
