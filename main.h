#ifndef MAIN_H
#define MAIN_H

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPoints2D.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData2DFS.h>

#include <string>
#include <iostream>  
#include <vector> 

#include <vtkButterflySubdivisionFilter.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkMetaImageWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPlane.h>
#include <vtkKdTree.h>
#include <vtkPolyDataNormals.h>
#include <vtkFloatArray.h>


#include <vtkPointPicker.h>
#include <vtkLookupTable.h>
#include <vtkCubeSource.h>


using namespace std;

#define CV_PI 3.141592653

#define CLRADIUS 1.5
#define  CIRCILEPOINTNUM 30
#define  R_INITIAL 0.25
#define  TH_F 5
#define  Hardness 2.0
#define POINTMASS 100.0
#define TIMESTEP 0.2
#define Jbarx 0.0
#define Jbary 0.0
#define Jbarz 0.0
#define CUBESTEP 1.0

const double BorderForce = 1 / (0.75 * 0.75 * 0.75 * 0.75 * 0.75 * 0.75) - 1 / 0.75;

const int range = 5;

#define CRADIUS 5.0
#define N 600

class CForce
{
public:
	double x;
	double y;
	double z;
public:
	CForce()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
};

void SaveVTKImage(vtkImageData *image, const char* fileName);

void SavePolyData(vtkPolyData *poly, const char* fileName);

int smoothvtkpolydata(vtkPolyData* Poly, int iternum, int TYPE = 1);

#define PutInAbsDecentOrder2(a, score)\
{\
	for (int i = 0; i < a.size(); i ++) \
{\
	for (int j = i + 1; j < a.size(); j ++) \
{\
	if(score.at(i) < score.at(j)) \
{\
	double tempscore = score.at(i);\
	score.at(i) = score.at(j);\
	score.at(j) = tempscore;\
	vtkIdType tempa = a.at(i);\
	a.at(i) = a.at(j);\
	a.at(j) = tempa;\
}\
}\
}\
}\


#define FindClosetDistance2CLi(clid, coord, dis)\
{\
	vtkIdType closestclpid = clpointLocator_eachcl[clid]->FindClosestPoint(coord);\
	double cloestclpoint[3];\
	clModel_eachcl[clid]->GetPoints()->GetPoint(closestclpid, cloestclpoint);\
	dis = sqrt(vtkMath::Distance2BetweenPoints(coord, cloestclpoint));\
}\

#define FindClosetDistanceandNorm2CLi(clid, coord, dis, dir)\
{\
	vtkIdType closestclpid = clpointLocator_eachcl[clid]->FindClosestPoint(coord);\
	double cloestclpoint[3];\
	clModel_eachcl[clid]->GetPoints()->GetPoint(closestclpid, cloestclpoint);\
	dis = sqrt(vtkMath::Distance2BetweenPoints(coord, cloestclpoint));\
	vtkMath::Subtract(coord, cloestclpoint, dir);\
	vtkMath::Normalize(dir);\
}\

#define PUTVECTOR3INORDER(a)\
{\
	for (int i = 0; i < 3; i ++)\
	{\
		for (int j = i + 1; j < 3; j ++)\
		{\
			if (a[i] > a[j])\
			{\
				double temp = a[i];\
				a[i] = a[j];\
				a[j] = temp;\
			}\
		}\
	}\
}\

#define INSERTNEWPOINT(coord)\
{\
	SamplePoly->GetPoints()->InsertNextPoint(coord);\
	double ndir[3] = {1.0, 0.0, 0.0};\
	NormalDir->InsertNextTuple(ndir);\
	pointflag->InsertNextValue(1);\
	R->InsertNextValue(R_INITIAL);	\
	R_Lumen->InsertNextValue(0.0); \
	F_sumabs->InsertNextValue(0.0);	\
	F_abssum->InsertNextValue(0.0);	\
	isSurfacePoint->InsertNextValue(0);\
	vector<vtkIdType> ThisAssociatedSegments;\
	ThisAssociatedSegments.push_back(CellId);\
	AssociatedSegments.push_back(ThisAssociatedSegments);\
	vtkAssociatedSegments->InsertNextTuple2(CellId, -1);\
}\

#define GENRANDOMPOINTNEARCOORD(coord, coord_new, R)\
{\
	double theta = vtkMath::Random(0*CV_PI, 1.0*CV_PI);\
	double phi = vtkMath::Random(0.0, 2.0*CV_PI);\
	coord_new[0] = coord[0] + R * sin(theta) * cos(phi);\
	coord_new[1] = coord[1] + R * sin(theta) * sin(phi);\
	coord_new[2] = coord[2] + R * cos(theta);\
}\

#define FINDMAXIDX(v, vlength, maxidx)\
{\
	double maxvalue = -10000.0;\
	for (int k = 0; k < 4; k ++)\
	{\
		if (v[k] > maxvalue)\
		{\
			maxvalue = v[k];\
			maxidx = k;\
		}\
	}\
}\

#define FINDMAXDIR(suggestdirs, givendir, maxidx, maxdir)\
{\
	double angle[4] = {0.0, 0.0, 0.0, 0.0};\
	for (int k = 0; k < 4; k ++)\
	{\
		for (int l = 0; l < 3; l ++)\
		{\
			angle[k] += givendir[l] * suggestdirs[k][l];\
		}\
	}\
	FINDMAXIDX(angle, 4, maxidx);\
	for (int l = 0; l < 3; l ++) maxdir[l] = suggestdirs[maxidx][l];\
}\


struct MyPoint
{
	int x;
	int y;
};
bool doIntersect(MyPoint p1, MyPoint q1, MyPoint p2, MyPoint q2);




#endif