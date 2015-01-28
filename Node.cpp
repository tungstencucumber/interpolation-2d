#include "Node.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <stdio.h>
#define R_I 1000000
#define R_D 1000000.0


Node::Node()
{
	nextStep = 0;
	trianglesNum = 0;
}

Node::Node(double* _c)
{
	for (int i=0; i<2; i++)
	{
		coords[i]=_c[i];
	}
	nextStep = 0;
	trianglesNum = 0;
}

Node::Node(Node* _n)
{
	for (int i=0; i<2; i++)
		coords[i]=_n->coords[i];
	setValues(_n->u);
	local_num = _n->local_num;
	nextStep = 0;
	trianglesNum = _n->trianglesNum;	
	for (int i=0; i<trianglesNum; i++)
		triangles[i]=_n->triangles[i];
	for (int i=0; i<4; i++)
		axis[i]=_n->axis[i];
	for (int i=0; i<4; i++)
		axis_method[i]=_n->axis_method[i];
//	if (local_num == 7) printf("C%d \n", trianglesNum);
}

Node::~Node()
{
	//do NOT delete nextStep!
}


double scalar(double* _v1, double* _v2)
{
	return _v1[0]*_v2[0]+_v1[1]*_v2[1];
}

double scalar(double _x0, double _y0, double _x1, double _y1)
{
	return _x0*_x1+_y0*_y1;
}

void Node::setValues(double* _v)
{
	if (!_v) return; 
	for (int i=0;i<5;i++) 
		u[i]=_v[i];
}

int Node::addTriangle(Triangle* t)
{
	if (!t) return 0;
	if (trianglesNum >= 30) return 0;
//	if (local_num == 7) printf("A%d \n", trianglesNum);
	triangles[trianglesNum] = t;
	trianglesNum++;
	if (nextStep)
		nextStep->addTriangle(t);
//	if (local_num == 7) printf("A%d \n", trianglesNum);
	return 1;
}

void Node::randomizeAxis()
{
	for (int i=0; i<4; i++)
		axis[i]= sqrt(2.0)/2.0;//0.0;//
	axis[3] = -sqrt(2.0)/2.0;
	axis[1] = axis[2] = 0.0; 
	axis[3] = axis[0] = 1.0; 
	return;
	axis[0] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	axis[1] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	double norma = sqrt(scalar(axis,axis));
	axis[0]/=norma; axis[1]/=norma; 
	
	if (fabs(axis[0]) > 0.0)
	{
		axis[3] = 1.0;
		axis[2] = -axis[1]/axis[0];
		//printf("norm: %lf     %lf %lf\n",norma, axis[2], axis[3]);
	}
	else if (fabs(axis[1]) > 0.0)
	{
		axis[2] = 1.0;
		axis[3] = -axis[0]/axis[1];
	}
	else 
	{
		printf("Fail in determining axis direction");
		return;
	}
	norma = sqrt(axis[2]*axis[2]+axis[3]*axis[3]);//scalar(axis+2,axis+2));	
	axis[2]/=norma; axis[3]/=norma; 
	//printf("axes randomized with %lf %lf %lf %lf\n",axis[0],axis[1],axis[2],axis[3]);
	//check!!!!!
	double check[3] = {scalar(axis,axis),scalar(axis+2,axis+2)
			,scalar(axis,axis+2)};
	//	printf ("axis random %.20lf %.20lf %.20lf\n",check[0],check[1],check[2]);
	if ((check[0] - 1.0)*(check[0] - 1.0) > 0.0001 || (check[1] - 1.0)*(check[1] - 1.0) > 0.0001 || check[2]*check[2] > 0.0001)
	{
		printf ("Warning! Recursion in axis random called! %lf %lf %lf   %lf %lf %lf %lf\n",check[0],check[1],check[2], axis[0], axis[1], axis[2], axis[3]);
		randomizeAxis();
		return;
	}
}
