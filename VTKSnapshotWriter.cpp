#include "VTKSnapshotWriter.h"
#include "Node.h"
#include "Triangle.h"


VTKSnapshotWriter::VTKSnapshotWriter()
{
	snapshot_writer_type.assign("Generic snapshot writer");
	resultdir = "./";
	//if (!strcmp("@", param))
		fname = string ("snap_%n.vtu");
}
VTKSnapshotWriter::VTKSnapshotWriter(const char *param)
{
	snapshot_writer_type.assign("Generic snapshot writer");
	resultdir = "./";
	//if (!strcmp("@", param))
		fname = string ("snap_%n.vtu");
	//else
	//	fname = string (param);
};

VTKSnapshotWriter::~VTKSnapshotWriter () {}

string* VTKSnapshotWriter::get_snapshot_writer_type()
{
	return &snapshot_writer_type;
};

int VTKSnapshotWriter::dump_tetr_mesh(Mesh* tetr_mesh, int zone_num, int snap_num)
{
	//if(tetr_mesh == NULL)
	//	throw GCMException( GCMException::SNAP_EXCEPTION, "No mesh provided");

	//*logger < "WARN: SnapshotWriter::dump_tetr_mesh - not yet implemented!";	

	return 0;
};

/*int VTKSnapshotWriter::dump_vtk(int snap_num)
{
	TetrMeshSet *mesh_set = TetrMeshSet::getInstance();
	for(int i = 0; i < mesh_set->get_number_of_local_meshes(); i++) 
		if( dump_vtk( mesh_set->get_local_mesh(i), snap_num ) < 0 )
			return -1;

	for(int i = 0; i < (mesh_set->external_meshes).size(); i++)
		if( dump_external_mesh_vtk( (mesh_set->external_meshes)[i], snap_num ) < 0 )
			return -1;

	return 0;
};*/

string VTKSnapshotWriter::getFileName (int zoneNum, int snapNum)
{
	string filename = fname;
	Utils::replaceAll (filename, "%z", Utils::t_to_string (zoneNum));
	Utils::replaceAll (filename, "%n", Utils::t_to_string (snapNum));
	return filename;
}

// TODO - think about local, remote, unused, etc
int VTKSnapshotWriter::dump_vtk(Mesh* tetr_mesh, int snap_num)
{
	//int zone_num = tetr_mesh->zone_num;
	string filename = getFileName (0, snap_num);
	filename = resultdir+filename;

	return dump_vtk (filename, tetr_mesh, snap_num);
}

int VTKSnapshotWriter::dump_external_mesh_vtk(Mesh* tetr_mesh, int snap_num)
{
	//int zone_num = tetr_mesh->external_mesh_zone_id;
	string filename = "external_" + fname;
	Utils::replaceAll(filename, "%z", Utils::t_to_string(0));
	Utils::replaceAll(filename, "%n", Utils::t_to_string(snap_num));
	filename = resultdir+filename;

	return dump_vtk (filename, tetr_mesh, snap_num);
}

int VTKSnapshotWriter::dump_vtk (string filename, Mesh* tetr_mesh, int snap_num)
{
	vtkXMLUnstructuredGridWriter *xgw = vtkXMLUnstructuredGridWriter::New();
	vtkUnstructuredGrid *g = vtkUnstructuredGrid::New();

	Node *node;
	//Tetrahedron_1st_order tetr;
	Triangle *tetr;

	vtkPoints *pts = vtkPoints::New();
	vtkDoubleArray *vx = vtkDoubleArray::New();
	vtkDoubleArray *vy = vtkDoubleArray::New();
	vtkDoubleArray *sxx = vtkDoubleArray::New();
	vtkDoubleArray *sxy = vtkDoubleArray::New();
	vtkDoubleArray *syy = vtkDoubleArray::New();
	vtkDoubleArray *axis_0 = vtkDoubleArray::New();
	vtkDoubleArray *axis_1 = vtkDoubleArray::New();
	//vtkDoubleArray *ggv = vtkDoubleArray::New();
	axis_0->SetNumberOfComponents(3);
	axis_1->SetNumberOfComponents(3);
	axis_0->SetName("axis_0");
	axis_1->SetName("axis_1");
	//ggv->SetNumberOfComponents(3);
	//ggv->SetName("2-nd gradient");

	/*vtkDoubleArray *vel = vtkDoubleArray::New();
	vel->SetNumberOfComponents(3);
	vel->SetName("velocity");
	vtkDoubleArray *contact = vtkDoubleArray::New();
	vtkIntArray	   *nodeFlags = vtkIntArray::New ();
	vtkDoubleArray *sxx = vtkDoubleArray::New();
	vtkDoubleArray *sxy = vtkDoubleArray::New();
	vtkDoubleArray *sxz = vtkDoubleArray::New();
	vtkDoubleArray *syy = vtkDoubleArray::New();
	vtkDoubleArray *syz = vtkDoubleArray::New();
	vtkDoubleArray *szz = vtkDoubleArray::New();
	vtkDoubleArray *la = vtkDoubleArray::New();
	vtkDoubleArray *mu = vtkDoubleArray::New();
	vtkDoubleArray *rho = vtkDoubleArray::New();
	vtkDoubleArray *yield_limit = vtkDoubleArray::New();
	vtkDoubleArray *destruction = vtkDoubleArray::New();
	vtkDoubleArray *deformationLimit = vtkDoubleArray::New();
	vtkDoubleArray *meshQuality = vtkDoubleArray::New();

	vtkDoubleArray *maxCompression = vtkDoubleArray::New();
	vtkDoubleArray *maxTension = vtkDoubleArray::New();
	vtkDoubleArray *maxShear = vtkDoubleArray::New();
	vtkDoubleArray *maxDeviator = vtkDoubleArray::New();

	vtkDoubleArray *maxCompressionHist = vtkDoubleArray::New();
	vtkDoubleArray *maxTensionHist = vtkDoubleArray::New();
	vtkDoubleArray *maxShearHist = vtkDoubleArray::New();
	vtkDoubleArray *maxDeviatorHist = vtkDoubleArray::New();

	vtkDoubleArray *maxDeformation = vtkDoubleArray::New();*/

	float vec[3],vec0[3];

	for(int i = 0; i < tetr_mesh->getNodesNum(); i++) {
		node = tetr_mesh->getNode(i);
		//printf("%lf %lf %lf %lf\n",node->u[0],node->u[1],node->u[2],node->u[3]);
		pts->InsertNextPoint( node->coords[0], node->coords[1], 0.0 );
		vx->InsertNextValue( node->u[0] );
		vy->InsertNextValue( node->u[1] );
		sxx->InsertNextValue( node->u[2] );
		sxy->InsertNextValue( node->u[3] );
		syy->InsertNextValue( node->u[4] );
		vec[0] = node->axis[0];	vec[1] = node->axis[1];	vec[2] = 0.0;
		axis_0->InsertNextTuple(vec);
		vec0[0] = node->axis[2]; vec0[1] = node->axis[3]; vec0[2] = 0.0;
		axis_1->InsertNextTuple(vec0);

		/*sxx->InsertNextValue( node.values[3] );
		sxy->InsertNextValue( node.values[4] );
		sxz->InsertNextValue( node.values[5] );
		syy->InsertNextValue( node.values[6] );
		syz->InsertNextValue( node.values[7] );
		szz->InsertNextValue( node.values[8] );
		la->InsertNextValue( node.la );
		mu->InsertNextValue( node.mu );
		rho->InsertNextValue( node.rho );
		yield_limit->InsertNextValue( node.yield_limit );
		destruction->InsertNextValue( node.destruction );
		deformationLimit->InsertNextValue( node.deformation_limit );
		meshQuality->InsertNextValue( node.mesh_quality );
		contact->InsertNextValue( node.isInContact () ? 1: 0 );
		nodeFlags->InsertNextValue (node.getFlags ());

		maxCompression->InsertNextValue( node.max_compression );
		maxTension->InsertNextValue( node.max_tension );
		maxShear->InsertNextValue( node.max_shear );
		maxDeviator->InsertNextValue( node.max_deviator );
		maxCompressionHist->InsertNextValue( node.max_compression_history );
		maxTensionHist->InsertNextValue( node.max_tension_history );
		maxShearHist->InsertNextValue( node.max_shear_history );
		maxDeviatorHist->InsertNextValue( node.max_deviator_history );
		maxDeformation->InsertNextValue( node.max_deformation );*/
	}
	g->SetPoints(pts);

	vtkTriangle *tetra=vtkTriangle::New();
	for(int i = 0; i < tetr_mesh->getTrianglesNum(); i++) {
		tetr = tetr_mesh->getTriangle(i);
		tetra->GetPointIds()->SetId(0,tetr->vert[0]);
		tetra->GetPointIds()->SetId(1,tetr->vert[1]);
		tetra->GetPointIds()->SetId(2,tetr->vert[2]);
		g->InsertNextCell(tetra->GetCellType(),tetra->GetPointIds());
	}

	vx->SetName("vx");
	vy->SetName("vy");
	sxx->SetName("sxx");
	sxy->SetName("sxy");
	syy->SetName("syy");
	//gv->SetName("grad v");
	/*sxx->SetName("sxx");
	sxy->SetName("sxy");
	sxz->SetName("sxz");
	syy->SetName("syy");
	syz->SetName("syz");
	szz->SetName("szz");
	la->SetName("lambda");
	mu->SetName("mu");
	rho->SetName("rho");
	yield_limit->SetName("yieldLimit");
	destruction->SetName("destruction");
	deformationLimit->SetName("deformationLimit");
	meshQuality->SetName("meshQuality");
	contact->SetName("contact");
	nodeFlags->SetName ("flags");

	maxCompression->SetName("maxCompression");
	maxTension->SetName("maxTension");
	maxShear->SetName("maxShear");
	maxDeviator->SetName("maxDeviator");
	maxCompressionHist->SetName("maxCompressionHistory");
	maxTensionHist->SetName("maxTensionHistory");
	maxShearHist->SetName("maxShearHistory");
	maxDeviatorHist->SetName("maxDeviatorHistory");
	maxDeformation->SetName("maxDeformation");*/

	g->GetPointData()->SetVectors(axis_0);
	g->GetPointData()->AddArray(vx);
	g->GetPointData()->AddArray(vy);
	g->GetPointData()->AddArray(sxx);
	g->GetPointData()->AddArray(sxy);
	g->GetPointData()->AddArray(syy);
	//g->GetPointData()->SetVectors(axis_1);
	//g->GetPointData()->SetVectors(gv);
	/*g->GetPointData()->AddArray(sxx);
	g->GetPointData()->AddArray(sxy);
	g->GetPointData()->AddArray(sxz);
	g->GetPointData()->AddArray(syy);
	g->GetPointData()->AddArray(syz);
	g->GetPointData()->AddArray(szz);
	g->GetPointData()->AddArray(la);
	g->GetPointData()->AddArray(mu);
	g->GetPointData()->AddArray(rho);
	g->GetPointData()->AddArray(yield_limit);
	g->GetPointData()->AddArray(destruction);
	g->GetPointData()->AddArray(deformationLimit);
	g->GetPointData()->AddArray(meshQuality);
	g->GetPointData()->AddArray(contact);
	g->GetPointData ()->AddArray (nodeFlags);

	g->GetPointData()->AddArray(maxCompression);
	g->GetPointData()->AddArray(maxTension);
	g->GetPointData()->AddArray(maxShear);
	g->GetPointData()->AddArray(maxDeviator);
	g->GetPointData()->AddArray(maxCompressionHist);
	g->GetPointData()->AddArray(maxTensionHist);
	g->GetPointData()->AddArray(maxShearHist);
	g->GetPointData()->AddArray(maxDeviatorHist);
	g->GetPointData()->AddArray(maxDeformation);*/

	vx->Delete();
	vy->Delete();
	sxx->Delete();
	sxy->Delete();
	syy->Delete();
	axis_0->Delete();
	axis_1->Delete();
	//ggv->Delete();
	/*sxx->Delete();
	sxy->Delete();
	sxz->Delete();
	syy->Delete();
	syz->Delete();
	szz->Delete();
	la->Delete();
	mu->Delete();
	rho->Delete();
	yield_limit->Delete();
	destruction->Delete();
	deformationLimit->Delete();
	meshQuality->Delete();
	contact->Delete();
	nodeFlags->Delete ();

	maxCompression->Delete();
	maxTension->Delete();
	maxShear->Delete();
	maxDeviator->Delete();
	maxCompressionHist->Delete();
	maxTensionHist->Delete();
	maxShearHist->Delete();
	maxDeviatorHist->Delete();
	maxDeformation->Delete();*/

	//*logger << "Dumping VTK snapshot to file " < filename;

	xgw->SetInputData(g);
	xgw->SetFileName(filename.c_str());
	xgw->Update();

	xgw->Delete();
	g->Delete();
	pts->Delete();
	tetra->Delete();

	return 0;
};

//void VTKSnapshotWriter::dump(int snap_num)
//{
//	dump_vtk(snap_num);
//}

void VTKSnapshotWriter::parseArgs(int argc, char **argv)
{
	static struct option long_options[] =
	{
		{"output-dir", required_argument, 0, 'o'},
		{0           , 0                , 0, 0  }
	};
	
	int option_index = 0;

	int c;
	while ((c = getopt_long (argc, argv, "o:", long_options, &option_index)) != -1)
		switch (c)
		{
			case 'o':
				resultdir = optarg;
				if (resultdir[resultdir.length()-1] != '/')
					resultdir += '/';
				break;
			default:
				optind--;
				return;
		}
}

int VTKSnapshotWriter::init()
{
	return 0;
}
