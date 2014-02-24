#include "Triangle.h"
#include "Mesh.h"
#include "Node.h"
#include <math.h>
#include <stdio.h>


Triangle::Triangle()
{
}

Triangle::Triangle(int* _v, Mesh* mesh)
{
	for (int i=0; i<3; i++)
		vert[i]=_v[i];
	Node* nodes[3] = {mesh->getNode(vert[0]), mesh->getNode(vert[1]), mesh->getNode(vert[2])};
	aabb[0] = aabb[2] = nodes[0]->coords[0];
	aabb[1] = aabb[3] = nodes[0]->coords[1];
	for (int i=1; i<3; i++)
	{
		if (aabb[0] > nodes[i]->coords[0])
			aabb[0] = nodes[i]->coords[0];
		if (aabb[1] > nodes[i]->coords[1])
			aabb[1] = nodes[i]->coords[1];
		if (aabb[2] < nodes[i]->coords[0])
			aabb[2] = nodes[i]->coords[0];
		if (aabb[3] < nodes[i]->coords[1])
			aabb[3] = nodes[i]->coords[1];
	}
}

Triangle::~Triangle()
{
}

int Triangle::check(double* _crd, Mesh* mesh)
{
	//if (_crd[0] < aabb[0] || _crd[1] < aabb[1] || _crd[0] > aabb[2] || _crd[1] > aabb[3]) return 0; //fast cut of far thetrs

	Node* nodes[3] = {mesh->getNode(vert[0]), mesh->getNode(vert[1]), mesh->getNode(vert[2])};
	double axis[4];
	//check for a plane against node 0
	axis[0] = nodes[1]->coords[0] - nodes[2]->coords[0];
	axis[1] = nodes[1]->coords[1] - nodes[2]->coords[1];
	if (fabs(axis[0]) > 0.0)
	{
		axis[3] = 1.0;
		axis[2] = -axis[1]/axis[0];
	}
	else if (fabs(axis[1]) > 0.0)
	{
		axis[2] = 1.0;
		axis[3] = -axis[0]/axis[1];
	}
	else 
	{
		printf("Fail in determining axis direction");
		return 0;
	}
	if (   ((_crd[0]-nodes[2]->coords[0])*axis[2]+
		(_crd[1]-nodes[2]->coords[1])*axis[3])*
	       ((nodes[0]->coords[0]-nodes[2]->coords[0])*axis[2]+
		(nodes[0]->coords[1]-nodes[2]->coords[1])*axis[3])
									<-0.0000001)
	{
		//printf ("%lf %lf %lf node 0\n%lf %lf %lf node 1\n%lf %lf %lf node 2\n%lf %lf %lf node 3\n%lf %lf %lf crd\n",nodes[0]->coords[0],nodes[0]->coords[1],nodes[0]->coords[2],nodes[1]->coords[0],nodes[1]->coords[1],nodes[1]->coords[2],nodes[2]->coords[0],nodes[2]->coords[1],nodes[2]->coords[2],nodes[3]->coords[0],nodes[3]->coords[1],nodes[3]->coords[2],_crd[0],_crd[1],_crd[2]);

		return 0;
	}
	//check for a plane against node 1
	axis[0] = nodes[0]->coords[0] - nodes[2]->coords[0];
	axis[1] = nodes[0]->coords[1] - nodes[2]->coords[1];
	if (fabs(axis[0]) > 0.0)
	{
		axis[3] = 1.0;
		axis[2] = -axis[1]/axis[0];
	}
	else if (fabs(axis[1]) > 0.0)
	{
		axis[2] = 1.0;
		axis[3] = -axis[0]/axis[1];
	}
	else 
	{
		printf("Fail in determining axis direction");
		return 0;
	}
	if (   ((_crd[0]-nodes[2]->coords[0])*axis[2]+
		(_crd[1]-nodes[2]->coords[1])*axis[3])*
	       ((nodes[1]->coords[0]-nodes[2]->coords[0])*axis[2]+
		(nodes[1]->coords[1]-nodes[2]->coords[1])*axis[3])
									<-0.0000001)
	{
		//printf ("%lf %lf %lf node 0\n%lf %lf %lf node 1\n%lf %lf %lf node 2\n%lf %lf %lf node 3\n%lf %lf %lf crd\n",nodes[0]->coords[0],nodes[0]->coords[1],nodes[0]->coords[2],nodes[1]->coords[0],nodes[1]->coords[1],nodes[1]->coords[2],nodes[2]->coords[0],nodes[2]->coords[1],nodes[2]->coords[2],nodes[3]->coords[0],nodes[3]->coords[1],nodes[3]->coords[2],_crd[0],_crd[1],_crd[2]);

		return 0;
	}
	//check for a plane against node 2
	axis[0] = nodes[1]->coords[0] - nodes[0]->coords[0];
	axis[1] = nodes[1]->coords[1] - nodes[0]->coords[1];
	if (fabs(axis[0]) > 0.0)
	{
		axis[3] = 1.0;
		axis[2] = -axis[1]/axis[0];
	}
	else if (fabs(axis[1]) > 0.0)
	{
		axis[2] = 1.0;
		axis[3] = -axis[0]/axis[1];
	}
	else 
	{
		printf("Fail in determining axis direction");
		return 0;
	}
	if (   ((_crd[0]-nodes[0]->coords[0])*axis[2]+
		(_crd[1]-nodes[0]->coords[1])*axis[3])*
	       ((nodes[2]->coords[0]-nodes[0]->coords[0])*axis[2]+
		(nodes[2]->coords[1]-nodes[0]->coords[1])*axis[3])
									<-0.0000001)
	{
		//printf ("%lf %lf %lf node 0\n%lf %lf %lf node 1\n%lf %lf %lf node 2\n%lf %lf %lf node 3\n%lf %lf %lf crd\n",nodes[0]->coords[0],nodes[0]->coords[1],nodes[0]->coords[2],nodes[1]->coords[0],nodes[1]->coords[1],nodes[1]->coords[2],nodes[2]->coords[0],nodes[2]->coords[1],nodes[2]->coords[2],nodes[3]->coords[0],nodes[3]->coords[1],nodes[3]->coords[2],_crd[0],_crd[1],_crd[2]);

		return 0;
	}
	return 1;
}

/*
int Triangle::checkZborder( Mesh* mesh)
{	
	Node* nodes[4] = {mesh->getNode(vert[0]), mesh->getNode(vert[1]), mesh->getNode(vert[2]), mesh->getNode(vert[3])};
	double min=mesh->getMinZ(),
		max=mesh->getMaxZ();	
	if (nodes[0]->coords[2] == min && nodes[1]->coords[2] == min && nodes[2]->coords[2] == min) return 1;
	if (nodes[0]->coords[2] == min && nodes[1]->coords[2] == min && nodes[3]->coords[2] == min) return 1;
	if (nodes[0]->coords[2] == min && nodes[2]->coords[2] == min && nodes[3]->coords[2] == min) return 1;
	if (nodes[1]->coords[2] == min && nodes[2]->coords[2] == min && nodes[3]->coords[2] == min) return 1;
	if (nodes[0]->coords[2] == max && nodes[1]->coords[2] == max && nodes[2]->coords[2] == max) return 1;
	if (nodes[0]->coords[2] == max && nodes[1]->coords[2] == max && nodes[3]->coords[2] == max) return 1;
	if (nodes[0]->coords[2] == max && nodes[2]->coords[2] == max && nodes[3]->coords[2] == max) return 1;
	if (nodes[1]->coords[2] == max && nodes[2]->coords[2] == max && nodes[3]->coords[2] == max) return 1;
	return 0;
}*/
