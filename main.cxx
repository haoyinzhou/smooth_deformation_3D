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
		SamplePoly = vtkSmartPointer<vtkPolyData>::New();
		isControlPoint = vtkSmartPointer<vtkIdTypeArray>::New();
		ControlPointCoord = vtkSmartPointer<vtkDoubleArray>::New();
		RelaxShape = vtkSmartPointer<vtkPoints>::New();
	}

	~MouseInteractorStyle(){	}

	bool InitialPointcloud()
	{
		if (BoundaryPolynormalGenerator == NULL)
		{
			std::cerr << "cannot find BoundaryPolynormalGenerator" << std::endl;
			return false;
		}
		vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPolynormalGenerator->GetOutput()->GetPointData()->GetArray("Normals"));
		if (BoundaryNormalArray == NULL)
		{
			std::cerr << "cannot find BoundaryNormalArray" << std::endl;
			return false;
		}

		connectionCellArray->Reset();
	//	relationships.clear();
	//	parameters.clear();

		vtkSmartPointer<vtkPoints> SamplePoints = SamplePoly->GetPoints();
		for (int i = 0; i < SamplePoints->GetNumberOfPoints();)
		{
			double coordi[3] = { 0.0, 0.0, 0.0 };
			coordi[0] = vtkMath::Random(-3.0 * CRADIUS, 3.0 * CRADIUS);
			coordi[1] = vtkMath::Random(-3.0 * CRADIUS, 3.0 * CRADIUS);
			coordi[2] = vtkMath::Random(-3.0 * CRADIUS, 3.0 * CRADIUS);

			vtkIdType nearestboundaryPID = boundarypointLocator->FindClosestPoint(coordi);
			double boundarycoord[3];
			BoundaryPoly->GetPoint(nearestboundaryPID, boundarycoord);
			double dir_p2boundary[3];
			vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
			vtkMath::Normalize(dir_p2boundary);
			double boundarynormal[3];
			BoundaryNormalArray->GetTuple(nearestboundaryPID, boundarynormal);
			if (vtkMath::Dot(dir_p2boundary, boundarynormal) < 0)
				continue;
			SamplePoints->SetPoint(i, coordi);
			i ++;
		}

		for (int i = 0; i < R->GetNumberOfTuples(); i ++)
		{
			R->SetValue(i, R_INITIAL);
		}

		return true;
	}

	bool UniformRedistribution(int iterationnumber)
	{
		if (BoundaryPolynormalGenerator == NULL)
		{
			std::cerr << "cannot find BoundaryPolynormalGenerator" << std::endl;
			return false;
		}
		vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPolynormalGenerator->GetOutput()->GetPointData()->GetArray("Normals"));
		if (BoundaryNormalArray == NULL)
		{
			std::cerr << "cannot find BoundaryNormalArray" << std::endl;
			return false;
		}

		for (int iter = 0; iter < iterationnumber; iter++)
		{
			for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				double coordi[3];
				SamplePoly->GetPoints()->GetPoint(i, coordi);
				double Ri = R->GetValue(i);

				double Fi[3] = { 0.0, 0.0, 0.0 };
				double Fi_c[3] = { 0.0, 0.0, 0.0 };
				double Fi_p[3] = { 0.0, 0.0, 0.0 };

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
					if (vtkMath::Dot(dir_p2boundary, boundarynormal) > 0)
						F_mode_c = 0.0;
					else
						F_mode_c = -1000.0;

					for (int l = 0; l < 3; l++) Fi_c[l] += F_mode_c * boundarynormal[l];
				}

				// force between points
				double Fpi_abssum = 0.0;
				double Fpi_sumabs = 0.0;
				vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New();
				//pointLocator->FindPointsWithinRadius(4.0, coordi, NeighorpIds);
				pointLocator->FindClosestNPoints(10, coordi, NeighorpIds);

				for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj++)
				{
					vtkIdType j = NeighorpIds->GetId(idxj);
					if (i == j) continue;

					double coordj[3];
					SamplePoly->GetPoint(j, coordj);

					double dir[3];
					vtkMath::Subtract(coordi, coordj, dir);
					double dis = vtkMath::Norm(dir);
					vtkMath::Normalize(dir);

					double Rj = R->GetValue(j);
					double R_total = Ri + Rj;

					if (dis < 0.55 * R_total)
					{
						double temp = 0.5 / R_total * dis + 0.5;
						double Fij_p_mode = (1.0 / (temp*temp*temp*temp*temp*temp) - 1.0 / (temp));
						double Fij_p[3] = { 0.0, 0.0, 0.0 };
						for (int l = 0; l < 3; l++) Fij_p[l] = Fij_p_mode * dir[l];
						for (int l = 0; l < 3; l++) Fi_p[l] += Fij_p[l];
						Fpi_sumabs += vtkMath::Norm(Fij_p);

						//	connections[i].push_back(j);
					}
				}

				Fpi_abssum = vtkMath::Norm(Fi_p);

				// combine Fi_cl and Fi_p
				for (int l = 0; l < 3; l++) Fi[l] = Fi_c[l] + 3.0 * Fi_p[l];

				// move points based on force
				if (vtkMath::Norm(Fi) > 20.0)
				{
					vtkMath::Normalize(Fi);
					for (int l = 0; l < 3; l++) Fi[l] = 20.0 * Fi[l];
				}

				{
					double coordi_new[3];
					for (int l = 0; l < 3; l++) coordi_new[l] = coordi[l] + Fi[l] / POINTMASS * TIMESTEP;
					SamplePoly->GetPoints()->SetPoint(i, coordi_new);
				}

				F_abssum->SetValue(i, Fpi_abssum);
				F_sumabs->SetValue(i, Fpi_sumabs);
			}

			// update radius Ri
			//double meanR = CRADIUS * CRADIUS / SamplePoints->GetNumberOfPoints();
			for (int i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				double Ri = R->GetValue(i);
				double delta_Ri_Max = 3.0 * TIMESTEP;

				double delta_Ri = delta_Ri_Max;
				double CondenseForce = F_sumabs->GetValue(i) - F_abssum->GetValue(i) - 5.5 * BorderForce; // 6 means hexagon grid mesh
				//	delta_Ri = (delta_Ri_Max + 1.0 - 0.1 * exp(CondenseForce / Hardness)) * timestep;	
				double temp = 1.0 + 0.003 * CondenseForce;
				if (temp <= 0.05)
					delta_Ri = delta_Ri_Max;
				else
					delta_Ri = -3.0 * log(temp) * TIMESTEP;
				if (delta_Ri > delta_Ri_Max)
					delta_Ri = delta_Ri_Max;

				Ri += delta_Ri;
				Ri = Ri < 0.05 ? 0.05 : Ri;
				Ri = Ri > 10.0 ? 10.0 : Ri;

				R->SetValue(i, Ri);
			}
		}

		return true;
	}

	bool BuildConnections()
	{
		this->Connections.clear();

		// draw the connections
		connectionCellArray->Reset();

		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
		{
			double coordi[3];
			SamplePoly->GetPoints()->GetPoint(i, coordi);
			double Ri = R->GetValue(i);

			vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New();
			//pointLocator->FindPointsWithinRadius(4.0, coordi, NeighorpIds);
			pointLocator->FindClosestNPoints(10, coordi, NeighorpIds);

			vector<vtkIdType> Connection_i;

			for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj++)
			{
				vtkIdType j = NeighorpIds->GetId(idxj);
				if (i == j) continue;

				double coordj[3];
				SamplePoly->GetPoints()->GetPoint(j, coordj);

				double dir[3];
				vtkMath::Subtract(coordi, coordj, dir);
				double dis = vtkMath::Norm(dir);

				double Rj = R->GetValue(j);
				double R_total = Ri + Rj;

				if (dis < 0.60 * R_total)
					Connection_i.push_back(j);
			}

			this->Connections.push_back(Connection_i);
		}

		// refine connections
		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			for (int idj = 0; idj < this->Connections[i].size(); idj++)
		{
			vtkIdType j = this->Connections[i].at(idj);
			bool flag_iinjfile = false;
			for (int idk = 0; idk < this->Connections[j].size(); idk++)
			{
				if (i == this->Connections[j].at(idk))
				{
					flag_iinjfile = true;
					break;
				}
			}
			if (flag_iinjfile == false)
				this->Connections[j].push_back(i);
		}

		// build connectionCellArray for display
		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			for (int idj = 0; idj < this->Connections[i].size(); idj++)
			{
				vtkIdType j = this->Connections[i].at(idj);
				if (j <= i)	continue;
				vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
				line->GetPointIds()->SetId(0, i);
				line->GetPointIds()->SetId(1, j);
				connectionCellArray->InsertNextCell(line);
			}

		return true;
	}

	//bool BuildRelationshipsParameters()
	//{
	//	if (relationships.size() == 0)
	//		return false;
	//	parameters.clear();
	//	for (int r = 0; r < relationships.size(); r ++)
	//	{
	//		double coord1[3], coord2[3];
	//		SamplePoly->GetPoint(relationships[r][0], coord1);
	//		SamplePoly->GetPoint(relationships[r][1], coord2);
	//		vector<double> thisparam;
	//		double dis = sqrt(vtkMath::Distance2BetweenPoints(coord1, coord2));
	//		thisparam.push_back(dis);
	//		double forcescale = 1.0; // should be different for different parts
	//		thisparam.push_back(forcescale);
	//		parameters.push_back(thisparam);
	//	}
	//	return true;
	//}

	bool FindSurfaceControlPoints()
	{
		RelatedBoundaryPids.clear();

		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
		{
			double coordi[3];
			SamplePoly->GetPoint(i, coordi);

			vtkIdType nearestboundaryPID = boundarypointLocator->FindClosestPoint(coordi);
			double boundarycoord[3];
			BoundaryPoly->GetPoint(nearestboundaryPID, boundarycoord);
			double dir_p2boundary[3];
			vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
			double dis_p2boundary = vtkMath::Norm(dir_p2boundary);

			if (dis_p2boundary < 5e-1)
			{
				isControlPoint->SetValue(i, 1);
				ControlPointCoord->SetTuple(i, boundarycoord);
				RelatedBoundaryPids.push_back(nearestboundaryPID);
			}
			else
			{
				double tempcoord[3] = { 0.0, 0.0, 0.0 };
				isControlPoint->SetValue(i, 0);
				ControlPointCoord->SetTuple(i, tempcoord);
				RelatedBoundaryPids.push_back(-1);
			}
		}

		return true;
	}

	bool DeformationMotion()
	{
		if (this->Connections.size() == 0)
			return false;
		if (isControlPoint == NULL)
			return false;
		if (RelaxShape == NULL)
			return false;

		if (BoundaryPolynormalGenerator == NULL)
		{
			std::cerr << "cannot find BoundaryPolynormalGenerator" << std::endl;
			return false;
		}
		vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPolynormalGenerator->GetOutput()->GetPointData()->GetArray("Normals"));
		if (BoundaryNormalArray == NULL)
		{
			std::cerr << "cannot find BoundaryNormalArray" << std::endl;
			return false;
		}
		
		{
			vtkSmartPointer<vtkPoints> SamplePoints = SamplePoly->GetPoints();
			std::memset(&forces.at(0), 0, forces.size() * sizeof(CForce));

			for (int i = 0; i < RelaxShape->GetNumberOfPoints(); i++)
			{
				double RelaxCenter[3];
				RelaxShape->GetPoint(i, RelaxCenter);
				double CurrentCenter[3];
				SamplePoints->GetPoint(i, CurrentCenter);

				double** RelaxShape_local = new double*[Connections[i].size()];
				double** CurrentShape_local = new double*[Connections[i].size()];
				for (int j = 0; j < Connections[i].size(); j ++)
				{
					double RelaxCoord[3], CurrentCoord[3];

					RelaxShape->GetPoint(Connections[i][j], RelaxCoord);
					SamplePoints->GetPoint(Connections[i][j], CurrentCoord);

					vtkMath::Subtract(RelaxCoord, RelaxCenter, RelaxCoord);
					vtkMath::Subtract(CurrentCoord, CurrentCenter, CurrentCoord);
					
					RelaxShape_local[j] = new double[3];
					CurrentShape_local[j] = new double[3];

					for (int l = 0; l < 3; l ++)
					{
						RelaxShape_local[j][l] = RelaxCoord[l];
						CurrentShape_local[j][l] = CurrentCoord[l];
					}
				}

				// add a syn rotation
				//double** rot = new double*[3];
				//for (int j = 0; j < 3; j++)
				//	rot[j] = new double[3];
				//rot[0][0] = 0.8754; rot[0][1] = 0.2346; rot[0][2] = -0.4226;
				//rot[1][0] = -0.2588; rot[1][1] = 0.9659; rot[1][2] = 0;
				//rot[2][0] = 0.4082; rot[2][1] = 0.1094; rot[2][2] = 0.9063;
				//vtkMath::MultiplyMatrix(RelaxShape_local, rot, Connections[i].size(), 3, 3, 3, CurrentShape_local);
				//
				//for (int j = 0; j < 3; j++)
				//	delete[] rot[j];
				//delete[] rot;

				// calculate rotation
				double** RelaxShape_local_T = new double*[3];
				for (int j = 0; j < 3; j++)
					RelaxShape_local_T[j] = new double[Connections[i].size()];

				for (int j = 0; j < 3; j++)
					for (int k = 0; k < Connections[i].size(); k++)
					{
						RelaxShape_local_T[j][k] = RelaxShape_local[k][j];
					}

				double** C = new double*[3];
				for (int j = 0; j < 3; j++)
				{
					C[j] = new double[3];
				}

				vtkMath::MultiplyMatrix(RelaxShape_local_T, CurrentShape_local, 3, Connections[i].size(), Connections[i].size(), 3, C);
				double A[3][3];
				memcpy(&A[0][0], &C[0][0], 3 * sizeof(double));
				memcpy(&A[1][0], &C[1][0], 3 * sizeof(double));
				memcpy(&A[2][0], &C[2][0], 3 * sizeof(double));

				double U[3][3], w[3], VT[3][3];
				vtkMath::SingularValueDecomposition3x3(A, U, w, VT);

				double UVT[3][3];
				vtkMath::Multiply3x3(U, VT, UVT);
				
				double** rot = new double*[3];
				for (int j = 0; j < 3; j++)
					rot[j] = new double[3];
				memcpy(&rot[0][0], &UVT[0][0], 3 * sizeof(double));
				memcpy(&rot[1][0], &UVT[1][0], 3 * sizeof(double));
				memcpy(&rot[2][0], &UVT[2][0], 3 * sizeof(double));

				double** RelaxShape_local_rot = new double*[Connections[i].size()];
				for (int j = 0; j < Connections[i].size(); j++)
				{
					RelaxShape_local_rot[j] = new double[3];
				}
				
				vtkMath::MultiplyMatrix(RelaxShape_local, rot, Connections[i].size(), 3, 3, 3, RelaxShape_local_rot);
				
				//if (i < 5)
				//{
				//	//std::cout << "i = " << i << std::endl;
				//	//std::cout << "CurrentShape_local = " << std::endl;
				//	//for (int k = 0; k < Connections[i].size(); k++)
				//	//{
				//	//	std::cout << CurrentShape_local[k][0] << ", " << CurrentShape_local[k][1] << ", " << CurrentShape_local[k][2] << std::endl;
				//	//}
				//	//std::cout << "RelaxShape_local_rot = " << std::endl;
				//	//for (int k = 0; k < Connections[i].size(); k++)
				//	//{
				//	//	std::cout << RelaxShape_local_rot[k][0] << ", " << RelaxShape_local_rot[k][1] << ", " << RelaxShape_local_rot[k][2] << std::endl;
				//	//}
				//	//std::cout << "CurrentShape_local_T = " << std::endl;
				//	//for (int j = 0; j < 3; j++)
				//	//{
				//	//	for (int k = 0; k < Connections[i].size(); k ++)
				//	//		std::cout << CurrentShape_local_T[j][k] << ", ";
				//	//	std::cout << std::endl;
				//	//}	
				//	//
				//	//std::cout << "RelaxShape_local = " << std::endl;
				//	//for (int k = 0; k < Connections[i].size(); k++)
				//	//{
				//	//	std::cout << RelaxShape_local[k][0] << ", " << RelaxShape_local[k][1] << ", " << RelaxShape_local[k][2] << std::endl;
				//	//}
				//	//std::cout << "RelaxShape_local_rot = " << std::endl;
				//	//for (int k = 0; k < Connections[i].size(); k++)
				//	//{
				//	//	std::cout << RelaxShape_local_rot[k][0] << ", " << RelaxShape_local_rot[k][1] << ", " << RelaxShape_local_rot[k][2] << std::endl;
				//	//}
				//	//std::cout << "C = " << std::endl;
				//	//for (int j = 0; j < 3; j++)
				//	//{
				//	//	for (int k = 0; k < 3; k++)
				//	//		std::cout << C[j][k] << ", ";
				//	//	std::cout << std::endl;
				//	//}
				//	//std::cout << "A = " << std::endl;
				//	//for (int j = 0; j < 3; j++)
				//	//{
				//	//	for (int k = 0; k < 3; k++)
				//	//		std::cout << A[j][k] << ", ";
				//	//	std::cout << std::endl;
				//	//}
				//	//std::cout << "UVT = " << std::endl;
				//	//for (int j = 0; j < 3; j++)
				//	//{
				//	//	std::cout << UVT[j][0] << ", " << UVT[j][1] << ", " << UVT[j][2] << std::endl;
				//	//}
				//	//std::cout << "w = " << w[0] << ", " << w[1] << ", " << w[2] << std::endl;
				//}

				// calculate force sent from i
				for (int j = 0; j < Connections[i].size(); j++)
				{
					double pid = Connections[i][j];
					double CurrentCoord[3], GoalCoord[3];
					memcpy(CurrentCoord, &CurrentShape_local[j][0], 3 * sizeof(double));
					memcpy(GoalCoord, &RelaxShape_local_rot[j][0], 3 * sizeof(double));

					double dir2goal[3];
					vtkMath::Subtract(GoalCoord, CurrentCoord, dir2goal);
					double dis2goal = vtkMath::Norm(dir2goal);
					vtkMath::Normalize(dir2goal);

					double force_i2pid[3];
					double force_i2pidnorm = 10 * dis2goal;
					//force_i2pidnorm = force_i2pidnorm > 10.0 ? 10.0 : force_i2pidnorm;
					for (int l = 0; l < 3; l++)	force_i2pid[l] = force_i2pidnorm * dir2goal[l];

					forces[pid].x = forces[pid].x + force_i2pid[0];
					forces[pid].y = forces[pid].y + force_i2pid[1];
					forces[pid].z = forces[pid].z + force_i2pid[2];
				}

				// delete memory
				for (int j = 0; j < Connections[i].size(); j++)
				{
					delete[] RelaxShape_local[j];
					delete[] CurrentShape_local[j];
					delete[] RelaxShape_local_rot[j];
				}
				delete[] RelaxShape_local;
				delete[] CurrentShape_local;
				delete[] RelaxShape_local_rot;

				for (int j = 0; j < 3; j++)
					delete[] RelaxShape_local_T[j];
				delete[] RelaxShape_local_T;
				for (int j = 0; j < 3; j++)
					delete[] C[j];
				delete[] C;
				for (int j = 0; j < 3; j++)
					delete[] rot[j];
				delete[] rot;

			}

		/*	for (int r = 0; r < relationships.size(); r++)
			{
				vtkIdType pid1 = relationships[r][0];
				vtkIdType pid2 = relationships[r][1];

				double coord1[3], coord2[3];
				SamplePoly->GetPoint(pid1, coord1);
				SamplePoly->GetPoint(pid2, coord2);

				double dir[3];
				vtkMath::Subtract(coord1, coord2, dir);
				double dis = vtkMath::Norm(dir);
				vtkMath::Normalize(dir);

				double zerodistance = parameters[r][0];
				//	zerodistance = 0.5 * zerodistance;
				double forcescale = parameters[r][1];

				double reletivedis = (dis + 1e-7) / (zerodistance + 1e-7);
				//double thisforcenorm = -log(reletivedis);
				double thisforcenorm = -20.0 * forcescale * (reletivedis - 1.0);
				thisforcenorm = thisforcenorm > 10.0 ? 10.0 : thisforcenorm;
				thisforcenorm = thisforcenorm < -10.0 ? -10.0 : thisforcenorm;

				double thisforce[3];
				for (int l = 0; l < 3; l++)
					thisforce[l] = thisforcenorm * dir[l];

				for (int l = 0; l < 3; l++)
				{
					forces[pid1][l] = forces[pid1][l] + thisforce[l];
					forces[pid2][l] = forces[pid2][l] - thisforce[l];
				}
			}
		*/
			//// add boundary force
			//for (int i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			//{
			//	if (isControlPoint->GetValue(i) == 1)
			//		continue;
			//	double coordi[3];
			//	SamplePoly->GetPoint(i, coordi);
			//	vtkIdType nearestboundaryPID = boundarypointLocator->FindClosestPoint(coordi);
			//	double boundarycoord[3];
			//	BoundaryPoly->GetPoint(nearestboundaryPID, boundarycoord);
			//	double dir_p2boundary[3];
			//	vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
			//	double dis_p2boundary = vtkMath::Norm(dir_p2boundary);
			//	vtkMath::Normalize(dir_p2boundary);
			//
			//	double boundarynormal[3];
			//	BoundaryNormalArray->GetTuple(nearestboundaryPID, boundarynormal);
			//
			//	double thisboundaryforce_norm = 0.0;
			//	if (vtkMath::Dot(dir_p2boundary, boundarynormal) > 0)
			//		thisboundaryforce_norm = 0.0;
			//	else
			//		thisboundaryforce_norm = 100.0;
			//
			//	for (int l = 0; l < 3; l++)
			//	{
			//		forces[i][l] = forces[i][l] + thisboundaryforce_norm * dir_p2boundary[l];
			//	}
			//}

			for (int i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				// move points according to forces
				if (isControlPoint->GetValue(i) != 0)
				{
					double GoalCoord[3];
					ControlPointCoord->GetTuple(i, GoalCoord);
					double CurrentCoord[3];
					SamplePoly->GetPoint(i, CurrentCoord);

					double dir2goal[3];
					vtkMath::Subtract(GoalCoord, CurrentCoord, dir2goal);
					double dis2goal = vtkMath::Norm(dir2goal);
					vtkMath::Normalize(dir2goal);

					double Fi[3];
					double force_i2pidnorm = 10 * dis2goal;
					force_i2pidnorm = force_i2pidnorm > 50.0 ? 50.0 : force_i2pidnorm;
					for (int l = 0; l < 3; l++)	Fi[l] = force_i2pidnorm * dir2goal[l];

					double coordnew[3];
					for (int l = 0; l < 3; l++) coordnew[l] = CurrentCoord[l] + Fi[l] / POINTMASS * TIMESTEP;
					SamplePoly->GetPoints()->SetPoint(i, coordnew);
				}
				else
				{
					double Fi[3];
					Fi[0] = forces[i].x;
					Fi[1] = forces[i].y;
					Fi[2] = forces[i].z;

					if (vtkMath::Norm(Fi) > 50.0)
					{
						vtkMath::Normalize(Fi);
						vtkMath::MultiplyScalar(Fi, 50.0);
					}

					double coordi[3];
					SamplePoly->GetPoint(i, coordi);
					for (int l = 0; l < 3; l++) coordi[l] = coordi[l] + Fi[l] / POINTMASS * TIMESTEP;
					SamplePoly->GetPoints()->SetPoint(i, coordi);
				}
			}
		
			//connectionCellArray->Modified();
			//connectionPolyData->Modified();
			//SamplePoly->Modified();
			//renderWindow->Render();
		}
	
		return true;
	}


	virtual void OnKeyPress()
	{
		std::string key = this->Interactor->GetKeySym();
	
		if (key == "s") // start to uniformly sampling
		{
		//	if (InitialPointcloud() == false)
		//		return;

			if (UniformRedistribution(10) == true)
			{
				BuildConnections();
			//	BuildRelationshipsParameters();
				FindSurfaceControlPoints();

				connectionCellArray->Modified();
				connectionPolyData->Modified();
				SamplePoly->Modified();
				renderWindow->Render();

				RelaxShape->DeepCopy(SamplePoly->GetPoints());
			//	SavePolyData(SamplePoly, "C:\\work\\smooth_deformation_3D\\testdata\\SamplePoly.vtp");
			}
		}

		else if (key == "Left")
		{
			if (ControlPointCoord == NULL)
				return;

			for (int i = 0; i < BoundaryPoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				double boundarycoordi[3];
				BoundaryPoly->GetPoint(i, boundarycoordi);
				for (int l = 0; l < 1; l++)
					boundarycoordi[l] = 0.99 * boundarycoordi[l];
				BoundaryPoly->GetPoints()->SetPoint(i, boundarycoordi);
			}
			BoundaryPoly->GetPoints()->Modified();
			BoundaryPoly->Modified();
			renderWindow->Render();

			BoundaryPolynormalGenerator->Update();

			for (int i = 0; i < ControlPointCoord->GetNumberOfTuples(); i++)
			{
				if (isControlPoint->GetValue(i) != 1)
					continue;

				double newcontrolcoord[3];
				BoundaryPoly->GetPoint(RelatedBoundaryPids[i], newcontrolcoord);
				ControlPointCoord->SetTuple(i, newcontrolcoord);
			}

		}
		else if (key == "Right")
		{
			if (ControlPointCoord == NULL)
				return;

			for (int i = 0; i < BoundaryPoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				double boundarycoordi[3];
				BoundaryPoly->GetPoint(i, boundarycoordi);
				for (int l = 0; l < 1; l++)
					boundarycoordi[l] = 1.01 * boundarycoordi[l];
				BoundaryPoly->GetPoints()->SetPoint(i, boundarycoordi);
			}
			BoundaryPoly->GetPoints()->Modified();
			BoundaryPoly->Modified();
			renderWindow->Render();

			BoundaryPolynormalGenerator->Update();

			for (int i = 0; i < ControlPointCoord->GetNumberOfTuples(); i++)
			{
				if (isControlPoint->GetValue(i) != 1)
					continue;

				double newcontrolcoord[3];
				BoundaryPoly->GetPoint(RelatedBoundaryPids[i], newcontrolcoord);
				ControlPointCoord->SetTuple(i, newcontrolcoord);
			}
		}
		else if (key == "g")
		{
			forces.resize(SamplePoly->GetPoints()->GetNumberOfPoints());
			for (int iter = 0; iter < 20; iter ++)
			{
				DeformationMotion();
				connectionCellArray->Modified();
				connectionPolyData->Modified();
				SamplePoly->Modified();
				renderWindow->Render();
			}
		
		//	BuildRelationshipsParameters();
		}	

	}

public:

	// system
	int ClickCount;
	vtkRenderWindow* renderWindow;
	string dataname;
	
	// boundary
	vtkPolyData* BoundaryPoly;
	vtkPointLocator* boundarypointLocator;	
	vtkSmartPointer<vtkPolyDataNormals> BoundaryPolynormalGenerator;

	// p
	vtkPointLocator* pointLocator;
	vtkSmartPointer<vtkPolyData> SamplePoly;
	vtkDoubleArray* R;	// radius
	vtkDoubleArray* F_sumabs;	// force with respect to other points
	vtkDoubleArray* F_abssum;	// force with respect to other points
	vtkSmartPointer<vtkPoints> RelaxShape;

	vtkSmartPointer<vtkIdTypeArray> isControlPoint;
	vtkSmartPointer<vtkDoubleArray> ControlPointCoord;
	vector<vtkIdType> RelatedBoundaryPids;

	vector< vector<vtkIdType> > Connections;
//	vector< vector<vtkIdType> > relationships; // each one has two points ids.
//	vector< vector<double> > parameters; // each one has some parameters for calculating forces
	vector< CForce > forces;		// forces used in phrase 2
	
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
	vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPolynormalGenerator->GetOutput()->GetPointData()->GetArray("Normals"));
	if (BoundaryNormalArray == NULL)
	{
		std::cerr << "cannot find BoundaryNormalArray" << std::endl;
		return false;
	}


	vtkSmartPointer<vtkPointLocator> boundarypointLocator = vtkSmartPointer<vtkPointLocator>::New();
	boundarypointLocator->SetDataSet(BoundaryPoly);
	boundarypointLocator->AutomaticOn();
	boundarypointLocator->SetNumberOfPointsPerBucket(1);
	boundarypointLocator->BuildLocator();
	
	// generate initial SamplePoints
	vtkSmartPointer<vtkPoints> SamplePoints = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkIdTypeArray> isControlPoint = vtkSmartPointer<vtkIdTypeArray>::New();
	isControlPoint->SetName("isControlPoint");
	isControlPoint->SetNumberOfComponents(1);
	isControlPoint->SetNumberOfValues(SamplePoints->GetNumberOfPoints());
	vtkSmartPointer<vtkDoubleArray> ControlPointCoord = vtkSmartPointer<vtkDoubleArray>::New();
	ControlPointCoord->SetName("ControlPointCoord");
	ControlPointCoord->SetNumberOfComponents(3);
	ControlPointCoord->SetNumberOfTuples(SamplePoints->GetNumberOfPoints());

	for (int i = 0; i < N; )
	{
		double coordi[3] = {0.0, 0.0, 0.0};
		coordi[0] = vtkMath::Random(-CRADIUS, CRADIUS);
		coordi[1] = vtkMath::Random(-CRADIUS, CRADIUS);
		coordi[2] = vtkMath::Random(-CRADIUS, CRADIUS);
		//	if (vtkMath::Norm(coordi) > CRADIUS) continue;
		
		vtkIdType nearestboundaryPID = boundarypointLocator->FindClosestPoint(coordi);
		double boundarycoord[3];
		BoundaryPoly->GetPoint(nearestboundaryPID, boundarycoord);
		double dir_p2boundary[3];
		vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
		vtkMath::Normalize(dir_p2boundary);
		double boundarynormal[3];
		BoundaryNormalArray->GetTuple(nearestboundaryPID, boundarynormal);
		if (vtkMath::Dot(dir_p2boundary, boundarynormal) < 0)
			continue;	

		SamplePoints->InsertNextPoint(coordi);
		isControlPoint->InsertNextValue(0);
		double temp[3] = { 0.0, 0.0, 0.0 };
		ControlPointCoord->InsertNextTuple(temp);
		i ++;
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
		R->InsertNextValue(R_INITIAL);	// R is initialized with a small number
		F_abssum->InsertNextValue(0.0);
		F_sumabs->InsertNextValue(0.0);
	}

	SamplePoly->GetPointData()->AddArray(R);
	SamplePoly->GetPointData()->AddArray(F_abssum);
	SamplePoly->GetPointData()->AddArray(F_sumabs);

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
	actor1->GetProperty()->SetOpacity(0.1);
	renderer->AddActor(actor1);

	vtkSmartPointer<vtkPolyDataMapper> mapper2 =	vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper2->SetInputData(connectionPolyData); 
	vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
	actor2->SetMapper(mapper2);
	actor2->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
//	actor2->GetProperty()->SetOpacity(0.1);
	renderer->AddActor(actor2);

	renderer->SetBackground(1.0, 1.0, 1.0);

	renderer->ResetCamera();
	renderWindow->Render();
	renderWindow->SetWindowName("Show Smoothed Points");

	vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
	style->BoundaryPoly = BoundaryPoly;
	style->boundarypointLocator = boundarypointLocator;
	style->BoundaryPolynormalGenerator = BoundaryPolynormalGenerator;

	style->pointLocator = pointLocator;

	style->SamplePoly = SamplePoly;
	style->R = R;
	style->F_abssum = F_abssum;
	style->F_sumabs = F_sumabs;

	style->connectionPolyData = connectionPolyData;
	style->connectionCellArray = connectionCellArray;
	
	style->isControlPoint = isControlPoint;
	style->ControlPointCoord = ControlPointCoord;

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






