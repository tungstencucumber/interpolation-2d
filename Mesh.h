#ifndef MESH
#define MESH 1

class Node;
class Triangle;
class Mesh
{
public:
	Mesh();
	~Mesh();
	int init(double _h); //TODO: add filepath
	int getNodesNum() {return nn;};
	int getTrianglesNum() {return nt;};
	Node* getNode(int num);
	Triangle* getTriangle(int num);
	void transcend();
	Triangle* findTriangle(double* _crd, Node* node);
	int load_msh_file(char* file_name);
	int load_smsh_file(char* file_name);
	void addNode(Node* newNode);
	void addTriangle(Triangle* newTriangle);
	void setInitialConditionsLinearX(double a);
	void setInitialConditionsLinearY(double a);
	//void setInitialConditionsLinearZ(double a);
	void setInitialConditionsLinear(double a);
	void setInitialConditionsStep(double a, double w);
	void setInitialConditionsStepX(double a, double w);
	void setInitialConditionsStepY(double a, double w);
	void setInitialConditionsSin4(double x, double y, double r, double amp);
	//void setInitialConditionsStepZ(double a, double w);
	void setInitialConditionsGradient();
	void setInitialConditionsGradientSecond();
	//double getMinZ() {return borders[2];};
	//double getMaxZ() {return borders[5];};
private:
	int nn,nt;
	Node** nodes;
	//Thetraeder** thtrs;
	Triangle** triangles;
	double borders[4];//0,1 - min values for x,y; 2,3 - max values
				//!!!Warning, borders only cyclic and for cubic mesh 
	double h; //space step
	//node list
	//tetr list
};
#endif
