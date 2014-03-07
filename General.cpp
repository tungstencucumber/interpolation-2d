#include "General.h"
#include "Mesh.h"
#include "Method.h"
#include "Node.h"
#include "VTKSnapshotWriter.h"
#include <math.h>
#define c0 1.0
#define c1 1.0
#define amp 10.0
#define wdt 0.2
#define PI 3.1415926

General::General()
{
	mesh = 0;
	method = 0;
}

General::~General()
{
	if (mesh) delete mesh;
	if (method) delete method;
}

int General::init()
{
	finalStep = 250;
	time = 0;
	//timeStep = 0.001;
	snapStep = 5;			//TODO:load parameters from file
	double courant = 0.08;		//lambda*tau/h
	double h = 0.0125;
	double min_c = c1;
	if (c0<c1) min_c = c0;
	timeStep = courant*h/min_c;
	if (mesh = new Mesh())
		{if (mesh->init(h)) return 1;} 		//TODO: add filepath
	else return 1;
	if (method = new Method())
		{if (method->init(c0,c1, 2)) return 1;}		//TODO: when added different methods, must be chosen from xml
	else return 1;
	if (sw = new VTKSnapshotWriter())
		{if (sw->init()) return 1;}		//TODO: when added different methods, must be chosen from xml
	else return 1;

	//mesh->setInitialConditionsStepY(amp, wdt);
	//setTimeStep();
	mesh->setInitialConditionsSin4(0.5,0.5,0.2,amp);
	return 0;
}

void General::setTimeStep()
{	
	//timeStep = 1;
}

int General::step(int currentStep)
{
	if (!currentStep)
		sw->dump_vtk(mesh,0);
	setTimeStep();
	time += timeStep;
	for (int i=0; i<mesh->getNodesNum(); i++)
	{
		//printf("randomizing %d\n",i);
		Node* n = mesh->getNode(i);
		n->randomizeAxis();
	}
	for (int axis=0; axis<2; axis++)
	{
		//printf("counting axis %d\n",axis);
		for (int i=0; i<mesh->getNodesNum(); i++)
		{
			Node* n = mesh->getNode(i);
			method->count(mesh, n, timeStep/2.0, axis);//axis);
		}
		mesh->transcend();
		//printf("counted axis %d\n",axis);
	}
	if (!(currentStep%snapStep))
		sw->dump_vtk(mesh,currentStep+1);
	if (!(currentStep%snapStep)) //L1 for stepY
	{
		double f=0, diffL1=0, diffL2=0, center_x=0, center_y = 0, top=0, dwn=0, dh=0, l,diffL0=0,maxdiff;
		for (int i=0; i<mesh->getNodesNum(); i++)
		{
			//printf("randomizing %d\n",i);
			Node* n = mesh->getNode(i);
			center_x = 0.5 + time*c0;
			center_y = 0.5 + time*c1;
			for (; center_x > 1.0 || center_x < 0.0; center_x -= 1.0*center_x/fabs(center_x));
			for (; center_y > 1.0 || center_y < 0.0; center_y -= 1.0*center_y/fabs(center_y));
			l = sqrt((n->coords[0]-center_x)*(n->coords[0]-center_x)+(n->coords[1]-center_y)*(n->coords[1]-center_y));
			if (l < wdt)
			{
				double arg = PI*l/wdt/2.0;
				l = sin(PI/2.0 + arg);
				f = amp*l*l*l*l;
			}
			else f = 0.0;

			maxdiff = fabs(f-n->u[0])/amp;
			if (diffL0 < maxdiff) diffL0 = maxdiff;
//			if (f>0) printf("f=%lf u=%lf  center:%lf %lf  l=%lf \n",f, n->u[0], center_x, center_y, l);
			diffL1 += fabs(f-n->u[0])/amp;
			diffL2 += fabs(f-n->u[0])*fabs(f-n->u[0])/amp/amp;
//			n->u[0] = f;
		}
		diffL1 /= (double)mesh->getNodesNum();
		diffL2 = sqrt(diffL2/(double)mesh->getNodesNum());
		printf("diff L0 = %lf\t\tdiff L1 = %lf\t\tdiff L2 = %lf\n", diffL0, diffL1, diffL2);
	}
		
}
