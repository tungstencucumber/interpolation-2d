#ifndef METHOD
#define METHOD 1

class Node;
class Triangle;
class Mesh;
class Method
{
public:
	Method();
	~Method();
	int init();
	int init(double c0, double c1, int _order);
	void count(Mesh* mesh, Node* n, double timeStep, int axis);

	static double scalar(double* _v1, double* _v2);
	static double scalar(double _x0, double _y0, double _x1, double _y1);
	static double interpolate_1_order(Triangle* t, double* _crd, int val, Mesh* mesh); 
	static double interpolate_2_order(Triangle* t, double* _crd, Mesh* mesh, double* _res, Node* node); 
	static double interpolate_3_order(Triangle* t, double* _crd, Mesh* mesh); 
private:
	//static double axis[4]; 	// current local randomised axis
	double coeff[2]; 	// dv/dt + coeff[0]*dv/dx + coeff[1]*dv/vy = 0
			 	// coeffs in global coordinates
	int order;
	static void randomizeAxis(double* axes);
	void calculateCoeff(double* _c, double* axis); //transform coordinates in accordance with axis[]

	//static void intoRandomAxes(double* x, double* y, int flag); //0 if x-y, 1 if coord-axis
	//static void intoRandomAxes(double* c);
	static void intoRandomAxes(double* x, double* y, double *axes);
	static void fromRandomAxes(double* x, double* y, double *a);
	static void intoRandomAxes(double* c, double *axes);
	static void fromRandomAxes(double* c, double *axes);
};
#endif
