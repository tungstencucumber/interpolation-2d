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
	int init(int _order);
	void count_split(Mesh* mesh, Node* n, double timeStep);
    double la, mu, rh, c[5];

	static double scalar(double* _v1, double* _v2);
	static double scalar(double _x0, double _y0, double _x1, double _y1);
	double interpolate_1_order(Triangle* t, double* _crd, int val, Mesh* mesh); 
	double interpolate_2_order(Triangle* t, double* _crd, Mesh* mesh, double* _res);//, Node* node); 
	double interpolate_3_order(Triangle* t, double* _crd, Mesh* mesh); 
private:

    double om_neg[5][5], om[5][5];
	int order;
	static void randomizeAxis(double* axes);
	void calculateCoeff(double* _c, double* axis); //transform coordinates in accordance with axis[]

	//static void intoRandomAxes(double* x, double* y, int flag); //0 if x-y, 1 if coord-axis
	//static void intoRandomAxes(double* c);
	static void intoRandomAxes(double* x, double* y, double *axes);
	static void fromRandomAxes(double* x, double* y, double *a);
//	static void fromRandomAxesGrad(double* x, double* y, double *a);
	static void intoRandomAxes(double* c, double *axes);
	static void intoRandomAxesGrad(double* x, double* y, double *axes);
//	static void intoRandomAxesGrad(double* c, double *axes);
	static void fromRandomAxes(double* c, double *axes);
    static void intoAxes(Node* n, double* axis);
};
#endif
