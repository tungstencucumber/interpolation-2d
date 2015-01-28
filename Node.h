#ifndef NODE
#define NODE 1

class Triangle;
class Node
{
public:
	Node();
	Node(double* _c);
	Node(Node* _n);
	~Node();
	/*union crd
	{
		double c[3];
		struct nm
		{
			double x;
			double y;
			double z;
		} name; 
	}*/
	double coords[2];
	union 
	{
		double u[5];
		struct 
		{
			double v;
			double vx;
			double vy;
			double vxx;
			double vyy;
		};
	};
	double axis[4];
	double axis_method[4];
	void randomizeAxis();
	int local_num;
	void setValues(double* _v);
	int addTriangle(Triangle* t);
	Node* nextStep;
	Triangle* triangles[30];//tethraeders
	int trianglesNum;//tethraeders
};
#endif
