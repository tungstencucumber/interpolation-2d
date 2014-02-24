#include "General.h"
#include "Mesh.h"
#include "Method.h"
#include "Node.h"
#include "VTKSnapshotWriter.h"
#include <math.h>
#define c0 1.0
#define c1 0.0
#define amp 10.0
#define wdt 0.2

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
	timeStep = 0.005;
	snapStep = 5;			//TODO:load parameters from file
	double h = 0.025;
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
	//for (int axis=0; axis<2; axis++)
	//{
		//printf("counting axis %d\n",axis);
		for (int i=0; i<mesh->getNodesNum(); i++)
		{
			Node* n = mesh->getNode(i);
			method->count(mesh, n, timeStep, 0);//axis);
		}
		mesh->transcend();
		//printf("counted axis %d\n",axis);
	//}
	if (!(currentStep%snapStep))
		sw->dump_vtk(mesh,currentStep+1);
	/*if (!(currentStep%snapStep)) //L1 for stepY
	{
		double f=0, diffL1=0, diffL2=0, center=0, top=0, dwn=0, dh=0;
		for (int i=0; i<mesh->getNodesNum(); i++)
		{
			//printf("randomizing %d\n",i);
			Node* n = mesh->getNode(i);
			center = 5.0 + timeStep*currentStep*c1;
			for (; center > 10.0 || center < 0.0; center -= 10.0*center/fabs(center));
			dh = center - 5.0;
			center = 5.0;
			top = center + wdt/2.0;
			dwn = center - wdt/2.0;
			dh = n->coords[1] - dh;
			if (dh < top && dh > dwn) f=amp;
			else f=0.0;
			//printf("%lf %lf %lf %lf %lf \n",dwn, center, top, dh, f);
			diffL1 += fabs(f-n->u[0])/amp;
			diffL2 += fabs(f-n->u[0])*fabs(f-n->u[0])/amp/amp;
		}
		diffL1 /= (double)mesh->getNodesNum();
		diffL2 = sqrt(diffL2/(double)mesh->getNodesNum());
		printf("diff L1 = %lf\t\tdiff L2 = %lf\n", diffL1, diffL2);
	}*/
		
}
