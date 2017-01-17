#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPoints2D.h>
#include <vtkLine.h>
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
#include <vtkLine.h>
#include <vtkDataSetMapper.h>
#include <vtkTetra.h>
#include <vtkProperty.h>


#include <string>
#include <iostream>  
#include <vector> 

#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkMetaImageWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPlane.h>
#include <vtkKdTree.h>
#include <vtkPointLocator.h>
#include <vtkPolyLine.h>
#include <vtkPointLocator.h>
#include <vtkMatrix4x4.h>

#include "main.h"
 

using namespace std;


// Define interaction style
class MouseInteractorStyle: public vtkInteractorStyleTrackballCamera
{
public:

	static MouseInteractorStyle* New();
	vtkTypeMacro(MouseInteractorStyle, vtkInteractorStyleTrackballCamera);

	MouseInteractorStyle() 
	{
	}
	~MouseInteractorStyle()
	{
	}

	virtual void OnRightButtonDown()
	{
		ClickCount ++;

		SavePolyData(BoundaryPoly, "C:\\work\\smooth_deformation_3D\\testdata\\BoundaryPoly.vtp");
		vtkDoubleArray* BoundaryNormalArray = vtkDoubleArray::SafeDownCast(BoundaryPoly->GetPointData()->GetArray("Normals"));
		if (BoundaryNormalArray == NULL)
		{
			std::cerr << "cannot find BoundaryNormalArray" << std::endl;
			return;
		}

		vtkSmartPointer<vtkPoints> SamplePoints = SamplePoly->GetPoints();

		vector< vector<vtkIdType> > connections(SamplePoints->GetNumberOfPoints());

		for (vtkIdType i = 0; i < SamplePoints->GetNumberOfPoints(); i ++)
		{
			if (pointflag->GetValue(i) == 0) continue;
			double coordi[3];
			SamplePoints->GetPoint(i, coordi);
			double Ri = R->GetValue(i);

			double diri[3] = {0.0, 0.0, 0.0};
			for (int l = 0; l < 3; l ++) diri[l] = coordi[l];
			vtkMath::Normalize(diri);
	
			double Fi[3] = {0.0, 0.0, 0.0};
			double Fi_c[3] = {0.0, 0.0, 0.0};
			double Fi_p[3] = {0.0, 0.0, 0.0};

			// force between center and this point
			{
				vtkIdType nearestboundaryPID = boundarypointLocator->FindClosestPoint(coordi);
				double boundarycoord[3];
				BoundaryPoly->GetPoint(nearestboundaryPID, boundarycoord);
				double dir_p2boundary[3];
				vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
				vtkMath::Normalize(dir_p2boundary);
				double boundarynormal[3];
				BoundaryNormalArray->GetTuple(nearestboundaryPID, boundarynormal);

				double F_mode_c = 0.0;
				if (vtkMath::Dot(dir_p2boundary, boundarynormal) < 0)
					F_mode_c = 0.0;
				else
					F_mode_c = -10.0;

				for (int l = 0; l < 3; l ++) Fi_c[l] += F_mode_c * diri[l];
			}
		 
			// force between points
			double Fpi_abssum = 0.0;
			double Fpi_sumabs = 0.0;
			vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New(); 
			//pointLocator->FindPointsWithinRadius(4.0, coordi, NeighorpIds);
			pointLocator->FindClosestNPoints(20, coordi, NeighorpIds);

			for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj ++)
			{
				vtkIdType j = NeighorpIds->GetId(idxj);
				if (i == j) continue;

				double coordj[3];
				SamplePoints->GetPoint(j, coordj);

				double dir[3];
				vtkMath::Subtract(coordi, coordj, dir);
				double dis = vtkMath::Norm(dir);
				vtkMath::Normalize(dir);
				
				double Rj = R->GetValue(j);
				double R_total = Ri + Rj;

				if (dis < 0.55 * R_total)
				{	
					double temp = 0.5 / R_total * dis + 0.5;
					double Fij_p_mode = ( 1.0 / (temp*temp*temp*temp*temp*temp) - 1.0 / (temp));
					double Fij_p[3] = {0.0, 0.0, 0.0};
					for (int l = 0; l < 3; l ++) Fij_p[l] = Fij_p_mode * dir[l];
					for (int l = 0; l < 3; l ++) Fi_p[l] += Fij_p[l];
					Fpi_sumabs += vtkMath::Norm(Fij_p);

				//	connections[i].push_back(j);
				}
			}

			Fpi_abssum = vtkMath::Norm(Fi_p);
			
			// combine Fi_cl and Fi_p
			for (int l = 0; l < 3; l ++) Fi[l] = 100.0 * Fi_c[l] + 1.0 * Fi_p[l];

			// move points based on force
			if (vtkMath::Norm(Fi) > 10.0)
			{
				vtkMath::Normalize(Fi);						
				for (int l = 0; l < 3; l ++) Fi[l] = 10.0 * Fi[l];
			}
			double coordi_new[3];
			for (int l = 0; l < 3; l ++) coordi_new[l] = coordi[l] + Fi[l] / POINTMASS * TIMESTEP;
			//coordi_new[2] = 0.0; 

			SamplePoints->SetPoint(i, coordi_new);
			F_abssum->SetValue(i, Fpi_abssum);
			F_sumabs->SetValue(i, Fpi_sumabs);
		}

		// update radius Ri
		//double meanR = CRADIUS * CRADIUS / SamplePoints->GetNumberOfPoints();
		for (int i = 0; i < SamplePoints->GetNumberOfPoints(); i ++)
		{
			if (pointflag->GetValue(i) == 0)
				continue;

			double Ri = R->GetValue(i);
			double delta_Ri_Max = 1.0 * TIMESTEP;

			double delta_Ri = delta_Ri_Max;
			double CondenseForce = F_sumabs->GetValue(i) - F_abssum->GetValue(i) - 5.5 * BorderForce; // 6 means hexagon grid mesh
			//	delta_Ri = (delta_Ri_Max + 1.0 - 0.1 * exp(CondenseForce / Hardness)) * timestep;	
			double temp = 1.0 + 0.003 * CondenseForce;
			if (temp <= 0.05)
				delta_Ri = delta_Ri_Max;
			else
				delta_Ri = -log(temp) * TIMESTEP;	
			if (delta_Ri > delta_Ri_Max)
				delta_Ri = delta_Ri_Max;

			Ri += delta_Ri;	 
			Ri = Ri < 0.05? 0.05 : Ri;
			Ri = Ri > 10.0? 10.0 : Ri;

			R->SetValue(i, Ri);
		}		

		// draw the connections
		if (ClickCount > 00.0)
		{
			std::cout << "ClickCount > 00.0" << std::endl;
			connectionCellArray->Reset();

			for (vtkIdType i = 0; i < SamplePoints->GetNumberOfPoints(); i++)
			{
				if (pointflag->GetValue(i) == 0) continue;
				double coordi[3];
				SamplePoints->GetPoint(i, coordi);
				double Ri = R->GetValue(i);

				vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New();
				//pointLocator->FindPointsWithinRadius(4.0, coordi, NeighorpIds);
				pointLocator->FindClosestNPoints(20, coordi, NeighorpIds);

				for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj++)
				{
					vtkIdType j = NeighorpIds->GetId(idxj);
					if (i == j) continue;

					double coordj[3];
					SamplePoints->GetPoint(j, coordj);

					double dir[3];
					vtkMath::Subtract(coordi, coordj, dir);
					double dis = vtkMath::Norm(dir);

					double Rj = R->GetValue(j);
					double R_total = Ri + Rj;

					if (dis < 0.65 * R_total)
						connections[i].push_back(j);
				}
			}
			// refine connections
			for (vtkIdType i = 0; i < SamplePoints->GetNumberOfPoints(); i ++)
				for (int idj = 0; idj < connections[i].size(); idj ++)
			{
				vtkIdType j = connections[i].at(idj);
				bool flag_iinjfile = false;
				for (int idk = 0; idk < connections[j].size(); idk ++)
				{
					if (i == connections[j].at(idk))
					{
						flag_iinjfile = true;
						break;
					}
				}
				if (flag_iinjfile == false)
					connections[j].push_back(i);
			}
		
			for (vtkIdType i = 0; i < SamplePoints->GetNumberOfPoints(); i ++)
				for (int idj = 0; idj < connections[i].size(); idj ++)
			{
				vtkIdType j = connections[i].at(idj);
				vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
				line->GetPointIds()->SetId(0, i); 
				line->GetPointIds()->SetId(1, j);
				connectionCellArray->InsertNextCell(line);
			}

			connectionCellArray->Modified();
			connectionPolyData->Modified();
		}

		SamplePoints->Modified();
		SampleCell->Modified();
		SamplePoly->Modified();
		renderWindow->Render();

		std::cout << ClickCount << " calculation is done" << std::endl;

	//	SavePolyData(connectionPolyData, "C:\\work\\smooth_deformation_3D\\testdata\\connectionPolyData.vtp");

		// forward events
	//	vtkInteractorStyleTrackballCamera::OnRightButtonDown();
	}

public:

	// system
	int ClickCount;
	vtkRenderWindow* renderWindow;
	string dataname;
	
	// boundary
	vtkPolyData* BoundaryPoly;
	vtkPointLocator* boundarypointLocator;
	
	
	// p
	vtkPointLocator* pointLocator;
	vtkSmartPointer<vtkPolyData> SamplePoly;
	vtkCellArray* SampleCell;
	vtkIntArray* pointflag;
	vtkDoubleArray* R;	// radius
	vtkDoubleArray* F_sumabs;	// force with respect to other points
	vtkDoubleArray* F_abssum;	// force with respect to other points

	vtkPolyData* connectionPolyData;
	vtkCellArray* connectionCellArray;
};
vtkStandardNewMacro(MouseInteractorStyle);


int main(int argc, char *argv[])
{	
	vtkSmartPointer<vtkSphereSource> sphereSource =	vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(CRADIUS);
	sphereSource->SetPhiResolution(60);
	sphereSource->SetThetaResolution(60);
	sphereSource->Update();
	vtkSmartPointer<vtkPolyData> BoundaryPoly = vtkSmartPointer<vtkPolyData>::New();
	BoundaryPoly = sphereSource->GetOutput();

	vtkSmartPointer<vtkPolyDataNormals> BoundaryPolynormalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	BoundaryPolynormalGenerator->SetInputData(BoundaryPoly);
	BoundaryPolynormalGenerator->ComputePointNormalsOn();
	BoundaryPolynormalGenerator->ConsistencyOn();
	BoundaryPolynormalGenerator->AutoOrientNormalsOn();
	BoundaryPolynormalGenerator->Update();

	
	// generate initial SamplePoints
	vtkSmartPointer<vtkPoints> SamplePoints = vtkSmartPointer<vtkPoints>::New();
	for (int i = 0; i < N; i ++)
	{
		double coordi[3] = {0.0, 0.0, 0.0};
		coordi[0] = vtkMath::Random(-CRADIUS, CRADIUS);
		coordi[1] = vtkMath::Random(-CRADIUS, CRADIUS);
		coordi[2] = vtkMath::Random(-CRADIUS, CRADIUS);
		if (vtkMath::Norm(coordi) > CRADIUS) continue;

		SamplePoints->InsertNextPoint(coordi);
	}

	vtkSmartPointer<vtkCellArray> SampleCell = vtkSmartPointer<vtkCellArray>::New();
	for (vtkIdType i = 0; i < SamplePoints->GetNumberOfPoints(); i ++)
	{
		 SampleCell->InsertNextCell( 1, &i );
	}
	vtkSmartPointer<vtkPolyData> SamplePoly = vtkSmartPointer<vtkPolyData>::New();
	SamplePoly->SetPoints(SamplePoints);
	SamplePoly->SetVerts(SampleCell);
//	SavePolyData(SamplePoly, "C:\\work\\smooth_deformation_3D\\testdata\\SamplePoly.vtp");

	vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
	pointLocator->SetDataSet(SamplePoly);
	pointLocator->AutomaticOn();
	pointLocator->SetNumberOfPointsPerBucket(1);
	pointLocator->BuildLocator();

	vtkSmartPointer<vtkPointLocator> boundarypointLocator = vtkSmartPointer<vtkPointLocator>::New();
	boundarypointLocator->SetDataSet(BoundaryPoly);
	boundarypointLocator->AutomaticOn();
	boundarypointLocator->SetNumberOfPointsPerBucket(1);
	boundarypointLocator->BuildLocator();


	vtkSmartPointer<vtkIntArray> pointflag = vtkSmartPointer<vtkIntArray>::New();
	pointflag->SetName("pointflag");
	pointflag->SetNumberOfComponents(1);
	vtkSmartPointer<vtkDoubleArray> R = vtkSmartPointer<vtkDoubleArray>::New(); // radius
	R->SetName("radius");
	R->SetNumberOfComponents(1);
	vtkSmartPointer<vtkDoubleArray> F_abssum = vtkSmartPointer<vtkDoubleArray>::New(); // radius
	F_abssum->SetName("F_abssum");
	F_abssum->SetNumberOfComponents(1);
	vtkSmartPointer<vtkDoubleArray> F_sumabs = vtkSmartPointer<vtkDoubleArray>::New(); // radius
	F_sumabs->SetName("F_sumabs");
	F_sumabs->SetNumberOfComponents(1);

	for (int i = 0; i < SamplePoints->GetNumberOfPoints(); i ++)
	{
		pointflag->InsertNextValue(1);
		R->InsertNextValue(R_INITIAL);	// R is initialized with a small number
		F_abssum->InsertNextValue(0.0);
		F_sumabs->InsertNextValue(0.0);
	}

	SamplePoly->GetPointData()->AddArray(R);
	SamplePoly->GetPointData()->AddArray(F_abssum);
	SamplePoly->GetPointData()->AddArray(F_sumabs);
	SamplePoly->GetPointData()->AddArray(pointflag);

	std::cout << "SamplePoints.num = " << SamplePoints->GetNumberOfPoints() << std::endl;

	// 
	vtkSmartPointer<vtkCellArray> connectionCellArray =	vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPolyData> connectionPolyData = vtkSmartPointer<vtkPolyData>::New();
	connectionPolyData->SetPoints(SamplePoints);
	connectionPolyData->SetLines(connectionCellArray);
	
	
	// Render window
	vtkSmartPointer<vtkRenderWindow> renderWindow =	vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(200* 1,200); //(width, height)
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkRenderer> renderer =	vtkSmartPointer<vtkRenderer>::New();
	renderWindow->AddRenderer(renderer);
	renderer->SetViewport(static_cast<double>(0)/1,0,static_cast<double>(0+1)/1,1);
	vtkSmartPointer<vtkPolyDataMapper> mapper =	vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(SamplePoly); 
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
	actor->GetProperty()->SetPointSize(5.0);
	renderer->AddActor(actor);

	vtkSmartPointer<vtkPolyDataMapper> mapper1 =	vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper1->SetInputData(BoundaryPoly); 
	vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
	actor1->SetMapper(mapper1);
	actor1->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
	//renderer->AddActor(actor1);

	vtkSmartPointer<vtkPolyDataMapper> mapper2 =	vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper2->SetInputData(connectionPolyData); 
	vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
	actor2->SetMapper(mapper2);
	actor2->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
	renderer->AddActor(actor2);

	renderer->SetBackground(1.0, 1.0, 1.0);

	renderer->ResetCamera();
	renderWindow->Render();
	renderWindow->SetWindowName("Show Smoothed Points");

	vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
	style->BoundaryPoly = BoundaryPolynormalGenerator->GetOutput();
	style->boundarypointLocator = boundarypointLocator;

	style->pointLocator = pointLocator;

	style->SamplePoly = SamplePoly;
	style->SampleCell = SampleCell;
	style->pointflag = pointflag;
	style->R = R;
	style->F_abssum = F_abssum;
	style->F_sumabs = F_sumabs;

	style->connectionPolyData = connectionPolyData;
	style->connectionCellArray = connectionCellArray;


	style->renderWindow = renderWindow;
	style->ClickCount = 0;

	renderWindowInteractor->SetInteractorStyle( style );		
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}







void SaveVTKImage(vtkImageData *image, const char* fileName)
{
	vtkSmartPointer< vtkMetaImageWriter > writer = vtkSmartPointer< vtkMetaImageWriter >::New();
	writer->SetFileName(fileName);
	writer->SetInputData(image);
	try
	{
		writer->Write();
	}
	catch(...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}

void SavePolyData(vtkPolyData *poly, const char* fileName)
{
	if (!poly) return;
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(poly);
	writer->SetFileName(fileName);
	writer->SetDataModeToBinary();
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}

int smoothvtkpolydata(vtkPolyData* Poly, int iternum, int TYPE)
{
	vtkPoints* Points = Poly->GetPoints();
	vtkCellArray* Strips = NULL;
	if (TYPE == 1)
		Strips = Poly->GetPolys();
	else
		Strips = Poly->GetStrips();

	vtkIdType CellId = 0;
	vtkIdType npts = 0, *pts = NULL;

	int* adjcent = (int*)malloc(Points->GetNumberOfPoints() * 10 * sizeof(int) );
	int* num_adjcent = (int*)malloc(Points->GetNumberOfPoints() * sizeof(int));


/*	for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid ++)
	{
		num_adjcent[pid] = 0;
		for(CellId = 0, Strips->InitTraversal(); Strips->GetNextCell(npts, pts); CellId ++)
		{
			if (pts[0] != pid && pts[1] != pid && pts[2] != pid)
				continue;

			for (vtkIdType j = 0; j < npts; j ++)
			{
				if (pts[j] == pid)
					continue;

				bool find_ptsj_in_adjofpid = false;
				for (int k = 0; k < num_adjcent[pid]; k ++)
				{
					if (adjcent[pid * 10 + k] == pts[j])
					{
						find_ptsj_in_adjofpid = true;
						break;
					}
				}

				if (find_ptsj_in_adjofpid == false)
				{
					adjcent[pid * 10 + num_adjcent[pid]] = pts[j];
					num_adjcent[pid] ++;
				}
			}
		}
	}
*/
	for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid ++)
	{
		num_adjcent[pid] = 0;
	}
	for(CellId = 0, Strips->InitTraversal(); Strips->GetNextCell(npts, pts); CellId ++)
	{
		if (npts != 3)
		{
			std::cout << "not triangle, smooth cannot work!" << std::endl;
			return 0;
		}
		for (int i = 0; i < npts; i ++)
		{
			int p[2];
			int pidx = 0;
			for (int k = 0; k < npts; k ++)
			{
				if (k != i)
				{
					p[pidx] = k;
					pidx ++;
				}
			}
			for (int l = 0; l < 2; l ++)
			{
				bool find_pl_in_adjofptsi = false;
				for (int k = 0; k < num_adjcent[pts[i]]; k ++)
				{
					if (adjcent[pts[i] * 10 + k] == pts[p[l]])
					{
						find_pl_in_adjofptsi = true;
						break;
					}
				}

				if (find_pl_in_adjofptsi == false)
				{
					adjcent[pts[i] * 10 + num_adjcent[pts[i]]] = pts[p[l]];
					num_adjcent[pts[i]] ++;
				}
			}
		}
	}
	

	// the smooth algorithm
	{
		vtkSmartPointer<vtkPoints> Points_orig = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPoints> Points_last = vtkSmartPointer<vtkPoints>::New();
		Points_orig->DeepCopy(Points);

		double* b = (double*)malloc(Points->GetNumberOfPoints() * 3 * sizeof(double));
		double pi[3], qi[3], oi[3];
		const double alpha = 0.1, beta = 0.2;

		for (int iter = 0; iter < iternum; iter ++)
		{		
			Points_last->DeepCopy(Points);
			for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid ++)
			{
				if (num_adjcent[pid] > 0)
				{
					pi[0] = 0;
					pi[1] = 0;
					pi[2] = 0;
					for (int j = 0; j < num_adjcent[pid]; j ++)
					{
						Points_last->GetPoint(adjcent[pid * 10 + j], qi);
						pi[0] += qi[0];
						pi[1] += qi[1];
						pi[2] += qi[2];
					}
					pi[0] = pi[0] / num_adjcent[pid];
					pi[1] = pi[1] / num_adjcent[pid];
					pi[2] = pi[2] / num_adjcent[pid];
					Points->SetPoint(pid, pi);
					Points_orig->GetPoint(pid, oi);
					Points_last->GetPoint(pid, qi);

					b[pid * 3 + 0] = pi[0] - (alpha * oi[0] + (1.0 - alpha) * qi[0]);
					b[pid * 3 + 1] = pi[1] - (alpha * oi[1] + (1.0 - alpha) * qi[1]);
					b[pid * 3 + 2] = pi[2] - (alpha * oi[2] + (1.0 - alpha) * qi[2]);

				}
			}
			for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid ++)
			{
				if (num_adjcent[pid] > 0)
				{
					double sumbj[3];
					sumbj[0] = 0;
					sumbj[1] = 0;
					sumbj[2] = 0;

					for (int j = 0; j < num_adjcent[pid]; j ++)
					{					
						sumbj[0] += b[adjcent[pid * 10 + j] * 3 + 0];
						sumbj[1] += b[adjcent[pid * 10 + j] * 3 + 1];
						sumbj[2] += b[adjcent[pid * 10 + j] * 3 + 2];
					}
					sumbj[0] = sumbj[0] / num_adjcent[pid];
					sumbj[1] = sumbj[1] / num_adjcent[pid];
					sumbj[2] = sumbj[2] / num_adjcent[pid];

					Points->GetPoint(pid, pi);

					pi[0] = pi[0] - (beta * b[pid * 3 + 0] + (1.0 - beta) * sumbj[0]);
					pi[1] = pi[1] - (beta * b[pid * 3 + 1] + (1.0 - beta) * sumbj[1]);
					pi[2] = pi[2] - (beta * b[pid * 3 + 2] + (1.0 - beta) * sumbj[2]);

					Points->SetPoint(pid, pi);
				}
			}
		}
		Points_orig->Reset();
		Points_last->Reset();
		free(b);
	}

	Poly->SetPoints(Points);
	return 1;
}

int BuildConnectionMap(vtkPoints* points_in, vtkPointLocator* pointLocator, double suggestdirs[4][3])
{
	vector< vector<vtkIdType> > connectionmap (points_in->GetNumberOfPoints());
	
	for (int i = 0; i < points_in->GetNumberOfPoints(); i ++)
	{
		double coordi[3];
		points_in->GetPoint(i, coordi);

		vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New(); 
		pointLocator->FindClosestNPoints(20, coordi, NeighorpIds);

		for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj ++)
		{
			vtkIdType j = NeighorpIds->GetId(idxj);
			if (i == j) continue;

			double coordj[3];
			points_in->GetPoint(j, coordj);

			double dir[3];
			vtkMath::Subtract(coordj, coordi, dir);
			double dis = vtkMath::Norm(dir);
			vtkMath::Normalize(dir);

			int thissuggetid;
			double thissuggestdir[3];
			FINDMAXDIR(suggestdirs, dir, thissuggetid, thissuggestdir);
			 
			double angle = vtkMath::Dot(dir, thissuggestdir);
			if (angle < 0.85) continue;

			// score from dis and angle, if dis smaller angle larger, then score higher
			double scorej = -exp(0.1 * dis) + exp(1.0 * angle);

			

			
		}

	}
	
	return 1;
}



// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(MyPoint p, MyPoint q, MyPoint r)
{
	if (q.x <= vtkMath::Max(p.x, r.x) && q.x >= vtkMath::Min(p.x, r.x) &&
		q.y <= vtkMath::Max(p.y, r.y) && q.y >= vtkMath::Min(p.y, r.y))
		return true;

	return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(MyPoint p, MyPoint q, MyPoint r)
{
	// for details of below formula.
	int val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0;  // colinear

	return (val > 0)? 1: 2; // clock or counterclock wise
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(MyPoint p1, MyPoint q1, MyPoint p2, MyPoint q2)
{
	// Find the four orientations needed for general and
	// special cases
	int o1 = orientation(p1, q1, p2);
	int o2 = orientation(p1, q1, q2);
	int o3 = orientation(p2, q2, p1);
	int o4 = orientation(p2, q2, q1);

	// General case
	if (o1 != o2 && o3 != o4)
		return true;

	// Special Cases
	// p1, q1 and p2 are colinear and p2 lies on segment p1q1
	if (o1 == 0 && onSegment(p1, p2, q1)) return true;

	// p1, q1 and p2 are colinear and q2 lies on segment p1q1
	if (o2 == 0 && onSegment(p1, q2, q1)) return true;

	// p2, q2 and p1 are colinear and p1 lies on segment p2q2
	if (o3 == 0 && onSegment(p2, p1, q2)) return true;

	// p2, q2 and q1 are colinear and q1 lies on segment p2q2
	if (o4 == 0 && onSegment(p2, q1, q2)) return true;

	return false; // Doesn't fall in any of the above cases
}






