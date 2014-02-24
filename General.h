#ifndef GENERAL
#define GENERAL 1

class Mesh;
class Method;
class VTKSnapshotWriter;
class General
{
public:
	General();
	~General();
	int init();
	int step(int currentStep);
	int writeSnap(int currentStep) {return finalStep%currentStep;};
	int getFinalStep() {return finalStep;};
private:	
	int finalStep;
	int snapStep;
	double time;
	double timeStep;
	double courant;
	int timeStepMode; //0 - fixed, 1 - courant condition 
	Mesh* mesh;	
	Method* method;
	VTKSnapshotWriter* sw;

	void setTimeStep();
};
#endif
