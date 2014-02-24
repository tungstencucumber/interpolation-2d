#ifndef THETRAEDER
#define THETRAEDER 1
class Node;
class Mesh;
class Triangle
{
public:
	Triangle();
	Triangle(int* _v, Mesh* mesh);
	~Triangle();
	int check(double* _crd, Mesh* mesh);
	int vert[3];
	int local_num;
	//int checkZborder(Mesh* mesh);
private:
	double aabb[4];
};

#endif
