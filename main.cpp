//#include <SnapWriter.h>
#include "General.h"
#include <iostream>
#include <stdio.h>
int main()
{
	//SnapWriter *sw = new SnapWriter();
	General *g = new General();
	if (g->init())
	{
		std::cout <<"Houston, we have a problem\n";
		return 0;
	}
	for (int i=0; i<g->getFinalStep(); i++)
	{
		printf("Step %d started\n",i);
		g->step(i);
		//if (g.writeSnap(i))
		//	sw.write();
	}
	delete g;
	return 0;
}
