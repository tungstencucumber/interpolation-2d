#include "Mesh.h"
#include "Triangle.h"
#include "Node.h"
#include "Method.h"
#include "stdlib.h"
#include "stdio.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <math.h>
#define PI 3.1415926

using std::string;

Mesh::Mesh()
{
	nodes = 0;
	triangles = 0;
	nn = nt = 0;
}

Mesh::~Mesh()
{
	if (nodes) 
	{
		for (int i=0; i<nn; i++)
			delete nodes[i];
		free(nodes);
	}
	if (triangles) 
	{
		for (int i=0; i<nt; i++)
			delete triangles[i];
		free(triangles);
	}
}

int Mesh::init(double h)
{
	//h = _h;
	/*FILE* f=fopen("mesh.node","r");
	if (!f)
	{
		printf("FNF\n");
		return 1;
	}
	int _t[4],_n;
	double _c[3];
	int _ta[4];
	fscanf(f,"%d",&nn);
	printf("%d nodes:\n",nn);
	nodes = (Node**)malloc(sizeof(Node*)*nn);
	for (int i=0; i<nn; i++)
	{
		fscanf(f,"%d%lf%lf%lf",&_n,_c,_c+1,_c+2);
		if (_n!=i) { printf("Node file damaged\n"); return 1;};
		nodes[i] = new Node(_c);
		printf(" node %d: %lf %lf %lf\n",i,nodes[i]->coords[0],nodes[i]->coords[1],nodes[i]->coords[2]);
	}
	fscanf(f,"%d",&nt);
	printf("%d thtrs:\n",nt);
	thtrs = (Triangle**)malloc(sizeof(Triangle*)*nt);
	for (int i=0; i<nt; i++)
	{
		fscanf(f,"%d%d%d%d%d",&_n,_t,_t+1,_t+2,_t+3);
		for (int j=0; j<4; j++)
		{		
			_ta[j]=_t[j];
		}
		if (_n!=i) { printf("Thtr file damaged\n"); return 1;}; 
		thtrs[i] = new Triangle(_ta);
		printf(" thtr %d: %d %d %d %d\n",i,_t[0],_t[1],_t[2],_t[3]);
		for (int j=0; j<4; j++)
		{	
			nodes[_t[j]]->thetrs[nodes[_t[j]]->thetrsNum++]=thtrs[i];
		}
	}*/
	//load_msh_file("untitled_0.25.msh");
/*	if (h == 0.0125)
		load_smsh_file("untitled_0.0125.smsh");
	else if (h == 0.05)
		load_smsh_file("untitled_0.05.smsh");
	else if (h == 0.025)
		load_smsh_file("untitled_0.025.smsh");
	else if (h == 0.00625)
		load_smsh_file("untitled_0.00625.smsh");
	else if (h == 0.003125)
		load_smsh_file("untitled_0.003125.smsh");
	else if (h == 0.0015625)
		load_smsh_file("untitled_0.0015625.smsh");
	else return 1; */
/*	if (h == 0.0125)
		load_msh_file("untitled_0.0125.msh");
	else if (h == 0.05)
		load_msh_file("untitled_0.05.msh");
	else if (h == 0.025)
		load_msh_file("untitled_0.025.msh");
	else if (h == 0.00625)
		load_msh_file("untitled_0.00625.msh");
	else if (h == 0.003125)
		load_msh_file("untitled_0.003125.msh");
	else if (h == 0.0015625)
		load_msh_file("untitled_0.0015625.msh");
	else return 1;*/
	char fn[20];
	sprintf(fn, "untitled_%.2f.msh", h);
	//h = 0.025;
	//for (int i=0; i<n; i++)
	//	h *= 0.8;
	return load_msh_file(fn);
} //TODO: add filepath

//int Mesh::getNodesNum()
//{
//	return 0;
//}

Node* Mesh::getNode(int num)
{
	return nodes[num];
}

Triangle* Mesh::getTriangle(int num)
{
	return triangles[num];
}

void Mesh::transcend()
{
	Node* tmp;
	for (int i=0; i<nn; i++)
	{
		tmp = nodes[i]->nextStep;
		if (tmp == 0) printf("NULL Node called\n");
		delete nodes[i];
		nodes[i] = tmp;
	}
}

Triangle* Mesh::findTriangle(double* _crd, Node* node)
{
	//border correction
	while (_crd[0] < borders[0])
		_crd[0] += borders[2] - borders[0];
	while (_crd[1] < borders[1])
		_crd[1] += borders[3] - borders[1];
	//while (_crd[2] < borders[2])
	//	_crd[2] += borders[5] - borders[2];
	while (_crd[0] > borders[2])
		_crd[0] -= borders[2] - borders[0];
	while (_crd[1] > borders[3])
		_crd[1] -= borders[3] - borders[1];
	//while (_crd[2] > borders[5])
	//	_crd[2] -= borders[5] - borders[2];

	for (int i=0; i<node->trianglesNum; i++)
	{		
		if (node->triangles[i]->check(_crd,this))
		{
				return node->triangles[i];			
		}
	}
//	if (node->local_num == 7) printf("!!!Not found triangle in the neighbourhood crd: %lf %lf  node: %lf %lf\n",_crd[0],_crd[1],node->coords[0],node->coords[1]);
//	if (node->local_num == 7) printf("N%d ", node->local_num);
	for (int i=0; i<nt; i++)
		if (triangles[i])
			if (triangles[i]->check(_crd,this))
			{
//				if (node->local_num == 7) printf("!!!Thetr %i found for %lf %lf ----------------------\n",i,_crd[0],_crd[1]);
				if (!node->addTriangle(triangles[i]))
					printf("new triangle is not added!\n");
				return triangles[i];
			}
	return 0;
}


void Mesh::addNode(Node* newNode)
{
	//if (nn==2000) return 0;
	//printf("Adding node %5d\t",nn);
	nodes[nn] = newNode;
	nn++;
}

void Mesh::addTriangle(Triangle* newTriangle)
{
	//printf("Adding thetr %5d\t",nt);
	triangles[nt] = newTriangle;
	nt++;
	//printf("Added thetr %5d\t\n",nt);
}

int Mesh::load_smsh_file(const char* file_name)
{
	FILE* in = fopen(file_name,"r");
	if (!in) 
	{	
		printf("No smsh file\n");
		return 1;
	}
	printf("Reading file...");
	int Nn=0, Nt=0, n, t0=0, t1=0, t2=0;
	double x=0.0,y=0.0;
	Node* new_node=0;
	fscanf(in, "%d", &Nn);
	nodes = (Node**)malloc(sizeof(Node*)*Nn);
	for (int i=0; i<Nn; i++)
	{
		fscanf(in, "%d %lf %lf", &n, &x, &y);
		new_node = new Node();
		new_node->coords[0] = x;
		new_node->coords[1] = y;
		new_node->u[0] = new_node->u[1] = new_node->u[2] = 0;
		new_node->local_num = n;

		if (!i) 
		{
			borders[0] = new_node->coords[0];
			borders[1] = new_node->coords[1];
			borders[2] = new_node->coords[0];
			borders[3] = new_node->coords[1];
		}
		else
		{
			if (borders[0] > new_node->coords[0])
				borders[0] = new_node->coords[0];
			if (borders[1] > new_node->coords[1])
				borders[1] = new_node->coords[1];
			if (borders[2] < new_node->coords[0])
				borders[2] = new_node->coords[0];
			if (borders[3] < new_node->coords[1])
				borders[3] = new_node->coords[1];
		}

		addNode(new_node);
	}
	printf("Finished reading %d nodes\n",Nn);
	fscanf(in, "%d", &Nt);
	triangles = (Triangle**)malloc(sizeof(Triangle)*Nt);
	for (int i=0; i<Nt; i++)
	{
		fscanf(in, "%d   %d %d %d", &n, &t0, &t1, &t2);	
		int new_tetr_vert[3]={t0,t1,t2};
		Triangle *new_tetr = new Triangle(new_tetr_vert,this);
		addTriangle(new_tetr);
	}
	printf("Finished reading %d tetrs\n",Nt);
	fclose(in);
}

int Mesh::load_msh_file(const char* file_name)
{
	int fileVer;
	string str;
	int tmp_int;
	float tmp_float;
	int number_of_nodes;
	int number_of_elements;
	Node *new_node;

	std::ifstream infile;
	infile.open(file_name, std::ifstream::in);
	if(!infile.is_open())	
    {   
        printf("Can not open msh file %s\n", file_name);
        return 1;
    }
	printf("Reading file...");

	infile >> str;
	if(strcmp(str.c_str(),"$MeshFormat") != 0)
    {
	    printf("Wrong file format 0");
        return 1;
    }

	infile >> tmp_float >> tmp_int >> tmp_int;
	fileVer = (int)(tmp_float*10);

	infile >> str;
	if(strcmp(str.c_str(),"$EndMeshFormat") != 0)
    {
		printf("Wrong file format 1");
        return 1;
    }

	printf("INFO: Header Ok");
	infile >> str;
	if(strcmp(str.c_str(),"$Nodes") != 0)
    {
		printf("Wrong file format 2");
        return 1;
    }

	infile >> number_of_nodes;
	printf("the file contains %d nodes\n", number_of_nodes);
	nodes = (Node**)malloc(sizeof(Node*)*number_of_nodes);
	nn = 0;
	for(int i = 0; i < number_of_nodes; i++)
	{
		new_node = new Node();
		// Zero all values
		new_node->coords[0] = new_node->coords[1] = 0;
		new_node->u[0] = new_node->u[1] = new_node->u[2] = 0;
		//new_node->u[3] = 0;

		infile >> new_node->local_num;
		if(new_node->local_num > 0)
		{
			new_node->local_num--;
			infile >> new_node->coords[0] >> new_node->coords[1] >> new_node->coords[2];
			if (!i) 
			{
				borders[0] = new_node->coords[0];
				borders[1] = new_node->coords[1];
				borders[2] = new_node->coords[0];
				borders[3] = new_node->coords[1];
			}
			else
			{
				if (borders[0] > new_node->coords[0])
					borders[0] = new_node->coords[0];
				if (borders[1] > new_node->coords[1])
					borders[1] = new_node->coords[1];
				if (borders[2] < new_node->coords[0])
					borders[2] = new_node->coords[0];
				if (borders[3] < new_node->coords[1])
					borders[3] = new_node->coords[1];
		
	}
			new_node->local_num = i;
		}
		else
        {
			printf("Wrong file format 3");
            return 1;
        }
		addNode(new_node);
	}
	printf("Finished reading nodes");
	
	infile >> str;
	if(strcmp(str.c_str(),"$EndNodes") != 0)
    {
		printf("Wrong file format 4");
        return 1;
    }

	printf("INFO: Nodes Ok");
	
	infile >> str;
	if(strcmp(str.c_str(),"$Elements") != 0)
    {
		printf("Wrong file format 5");
        return 1;
    }

	infile >> number_of_elements;
	triangles = (Triangle**)malloc(sizeof(Triangle)*number_of_elements);
	nt = 0;
	printf("the file contains %d elements", number_of_elements);
	for(int i = 0; i < number_of_elements; i++)
	{
		int new_tetr_vert[3]={0}, local_num;
		infile >> tmp_int >> tmp_int;
		if(tmp_int != 2) {
			getline(infile, str);
		//printf("%d %d %d \n",new_tetr_vert[0],new_tetr_vert[1],new_tetr_vert[2]);
			continue;
		} else if (tmp_int == 2) {
			local_num = nt;
			if( fileVer == 22 ) {
				infile >> tmp_int >> tmp_int >> tmp_int 
					>> new_tetr_vert[0] >> new_tetr_vert[1] >> new_tetr_vert[2];
			} else {
				infile >> tmp_int >> tmp_int >> tmp_int >> tmp_int 
					>> new_tetr_vert[0] >> new_tetr_vert[1] >> new_tetr_vert[2];
			}

			if( (new_tetr_vert[0] <= 0) || (new_tetr_vert[1] <= 0) || (new_tetr_vert[2] <= 0))
            {
				printf("Wrong file format 6");
                return 1;
            }

			new_tetr_vert[0]--; new_tetr_vert[1]--; new_tetr_vert[2]--;
			Triangle *new_tetr = new Triangle(new_tetr_vert,this);
			addTriangle(new_tetr);
		}
	}
	printf("Finished reading elements");

	printf("INFO: Elements Ok");

	infile >> str;
	if(strcmp(str.c_str(),"$EndElements") != 0)
    {
		printf("Wrong file format 7");
        return 1;
    }

	printf("File successfully read.");

	infile.close();
	printf("\nMesh borders:\n%lf %lf  %lf %lf\n", borders[0], borders[1], borders[2], borders[3]);
	return 0;
};


void Mesh::setInitialConditionsStep(double a, double w)
{
	printf ("Step initial conditios:\n%lf %lf\n",a,w);
	double half[2] = {(borders[2]+borders[0])/2.0,(borders[3]+borders[1])/2.0};
	for (int i=0; i<nn; i++)
	{
		if (	nodes[i]->coords[0] < half[0] + w && nodes[i]->coords[0] > half[0] - w &&
			nodes[i]->coords[1] < half[1] + w && nodes[i]->coords[1] > half[1] - w)
		nodes[i]->vx = a;
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsStepX(double a, double w)
{
	printf ("Step initial conditios:\n%lf %lf\n",a,w);
	double half[2] = {(borders[2]+borders[0])/2.0,(borders[3]+borders[1])/2.0};
	for (int i=0; i<nn; i++)
	{
		if (	nodes[i]->coords[0] < half[0] + w && nodes[i]->coords[0] > half[0] - w 
		   )
		nodes[i]->vx = a;
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsStepY(double a, double w)
{
	printf ("Step initial conditios:\n%lf %lf\n",a,w);
	double half[2] = {(borders[2]+borders[0])/2.0,(borders[3]+borders[1])/2.0};
	for (int i=0; i<nn; i++)
	{
		if (	//nodes[i]->coords[0] < half[0] + w && nodes[i]->coords[0] > half[0] - w &&
			nodes[i]->coords[1] < half[1] + w && nodes[i]->coords[1] > half[1] - w 
			//nodes[i]->coords[2] < half[2] + w && nodes[i]->coords[2] > half[2] - w)
		   )
		nodes[i]->vx = a;
	}
	
	//printf ("Step initial conditios OK");
	setInitialConditionsGradient();
	//printf ("Step initial conditios Grad OK");
}

void Mesh::setInitialConditionsLinearY(double a)
{
	printf ("Linear Y initial conditios: %lf\n",a);
	double half = (borders[3]-borders[1])/2.0;
	for (int i=0; i<nn; i++)
	{
		nodes[i]->vx = a - a*(fabs(half-nodes[i]->coords[1]))/half;//fabs(half-nodes[i]->coords[0]) + fabs(half-nodes[i]->coords[1]) + 
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsLinearX(double a)
{
	printf ("Linear X initial conditios: %lf\n",a);
	double half = (borders[2]-borders[0])/2.0;
	for (int i=0; i<nn; i++)
	{
		nodes[i]->vx = a - a*(fabs(half-nodes[i]->coords[0]))/half;//fabs(half-nodes[i]->coords[0]) + fabs(half-nodes[i]->coords[1]) + 
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsLinear(double a)
{
	printf ("Linear initial conditios: %lf\n",a);
	double half = (borders[2]-borders[0])/2.0;
	for (int i=0; i<nn; i++)
	{
		nodes[i]->vx = a - a*(fabs(half-nodes[i]->coords[0]) + fabs(half-nodes[i]->coords[1]))/half;//
	}
	setInitialConditionsGradient();
}

void Mesh::setInitialConditionsSin4(double x, double y, double r, double amp)
{
	printf ("Sin^4 initial conditios: %lf\n",amp);
	double arg=0,l=0;
	for (int i=0; i<nn; i++)
	{
		l = sqrt((nodes[i]->coords[0]-x)*(nodes[i]->coords[0]-x)+(nodes[i]->coords[1]-y)*(nodes[i]->coords[1]-y));
		if (l < r)
		{
			arg = PI*l/r/2.0;
			l = sin(PI/2.0 + arg);
			nodes[i]->vx = amp*l*l*l*l;
            nodes[i]->sxy = nodes[i]->vx * sqrt(0.7);
            //nodes[i]->syy = amp*l*l*l*l;
		}
	}
	setInitialConditionsGradient();	
}

void Mesh::setInitialConditionsGradient()
{
/*
	double crd[2]={0.0},gP=0.0,gM=0.0;
	Triangle *tP=0, *tM=0;
	for (int i=0; i<nn; i++)
	{
		//h = ((nodes[i]->coords[0]-nodes[i]->coords[0])+(nodes[i]->coords[1]-nodes[i]->coords[1])+(nodes[i]->coords[2]-nodes[i]->coords[2]))/3.0;
		for (int ax=0; ax<2; ax++)
		{
			for (int c=0;c<2;c++) crd[c] = nodes[i]->coords[c];
			if (crd[ax]==borders[ax] || crd[ax]==borders[2+ax]) 
			{
				nodes[i]->u[1+ax] = 0.0;
				continue;
			}	
			crd[ax] += h/2.0;
			tP = findTriangle(crd,nodes[i]);
			gP = Method::interpolate_1_order(tP,crd,0,this);
			crd[ax] -= 1.0*h;
			tM=findTriangle(crd,nodes[i]);
			gM = Method::interpolate_1_order(tM,crd,0,this);
			nodes[i]->u[1+ax] = (gP-gM)/h;
			//printf("GRD: ax%d  u- %lf  u+ %lf  h %lf  grad %lf\n",ax,nodes[i]->u[0],gP,h,nodes[i]->u[1+ax]);
		}
	}
	//setInitialConditionsGradientSecond();
*/
}

void Mesh::setInitialConditionsGradientSecond()
{
/*
	double crd[2]={0.0},gP=0.0,gM=0.0;
	Triangle *tP=0, *tM=0;
	for (int i=0; i<nn; i++)
	{
		//h = ((nodes[i]->coords[0]-nodes[i]->coords[0])+(nodes[i]->coords[1]-nodes[i]->coords[1])+(nodes[i]->coords[2]-nodes[i]->coords[2]))/3.0;
		for (int ax=0; ax<2; ax++)
		{
			for (int c=0;c<2;c++) crd[c] = nodes[i]->coords[c];
			if (crd[ax]==borders[ax] || crd[ax]==borders[2+ax]) 
			{
				nodes[i]->u[3+ax] = 0.0;
				continue;
			}	
			crd[ax] += h/4.0;
			tP = findTriangle(crd,nodes[i]);
			gP = Method::interpolate_1_order(tP,crd,ax,this);
			crd[ax] -= h/2.0;
			tM=findTriangle(crd,nodes[i]);
			gM = Method::interpolate_1_order(tM,crd,ax,this);
			nodes[i]->u[3+ax] = 4.0*(gP-nodes[i]->u[ax])/h;
		}
	}
*/
}

