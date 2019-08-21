//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//  Copyright 2018 UT-Battelle, LLC.
//  Copyright 2018 Los Alamos National Security.
//
//  Under the terms of Contract DE-NA0003525 with NTESS,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================
#ifndef vtk_m_exec_cellmetrics_RelativeSizeSquared_h
#define vtk_m_exec_cellmetrics_RelativeSizeSquared_h

/*
 * Mesh quality metric functions that compute the RelativeSizeSquared of mesh cells.
 *
 * These metric computations are adapted from the VTK implementation of the Verdict library,
 * which provides a set of mesh/cell metrics for evaluating the geometric qualities of regions 
 * of mesh spaces. 
 *
 * See: The Verdict Library Reference Manual (for per-cell-type metric formulae)
 * See: vtk/ThirdParty/verdict/vtkverdict (for VTK code implementation of this metric)
 */

#include "vtkm/CellShape.h"
#include "vtkm/CellTraits.h"
#include "vtkm/VecTraits.h"
#include "vtkm/VectorAnalysis.h"
#include "vtkm/exec/FunctorBase.h"

#define UNUSED(expr) (void)(expr);

namespace vtkm
{
namespace exec
{
namespace cellmetrics
{

using FloatType = vtkm::FloatDefault;

// ========================= Unsupported cells ==================================

// By default, cells have zero shape unless the shape type template is specialized below.
template <typename OutType, typename PointCoordVecType, typename CellShapeType>
VTKM_EXEC OutType ComputeRelativeSizeSquared(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              CellShapeType shape,
			      const vtkm::Id& numShapes,
                              const vtkm::exec::FunctorBase&)
{
  UNUSED(numPts);  
  UNUSED(pts);  
  UNUSED(shape);  
  UNUSED(numShapes);
  return OutType(0.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeRelativeSizeSquared(const vtkm::IdComponent& numPts,
                                 const PointCoordVecType& pts,
                                 vtkm::CellShapeTagPolygon,
				 const vtkm::Id& numPolygons,
                                 const vtkm::exec::FunctorBase& worklet)
{
  switch(numPts)
  {
    case 4:
	    return ComputeRelativeSizeSquared<OutType>(numPts, pts, vtkm::CellShapeTagQuad(), numPolygons, worklet);
    default:
	    break;
  }
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeRelativeSizeSquared(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagLine,
				 const vtkm::Id& numLines,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  UNUSED(numLines);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeRelativeSizeSquared(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagWedge,
				 const vtkm::Id& numWedges,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  UNUSED(numWedges);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeRelativeSizeSquared(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagPyramid,
				 const vtkm::Id& numPyramids,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  UNUSED(numPyramids);
  return OutType(-1.0);
}

// ========================= 2D cells ==================================
// Compute the RelativeSizeSquared of a triangle.
// Formula: (min(R, 1/R))^2
//     - R = A/~A
//     - ~A = average area of all triangles in the mesh
// value for unit triangle dependent on ~A
// Acceptable range: [0.25,1]
// Normal range: [0,1]
// Full range: [0,1]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeRelativeSizeSquared(const vtkm::IdComponent& numPts,
                                 const PointCoordVecType& pts,
                                 vtkm::CellShapeTagTriangle,
				 const vtkm::Id& numTriangles,
                                 const vtkm::exec::FunctorBase& worklet)
{
  if(numPts != 3)
  {
    worklet.RaiseError("RelativeSizeSquared metric(tri) requires 3 points");
    return OutType(0.0);
  }
  //auto shapeCountPortal = worklet.GetShapeCounts().GetPortalConstControl();
  //const int numTriangles = shapeCountPortal.Get(vtkm::CELL_SHAPE_TRIANGLE);
  OutType weight1_1, weight2_1, weight1_2, weight2_2, weightedDeterminant, weightedJacobian, norm, result; 
  const OutType root3 = OutType(sqrt(3.0));

  using Edge = typename PointCoordVecType::ComponentType;
  const Edge TriEdges[3] = {pts[1] - pts[0],
			     pts[2] - pts[1],
                             pts[0] - pts[2]
	      		    };
  weight1_1 = OutType(1);
  weight2_1 = OutType(0.0);
  weight1_2 = OutType(0.5);
  weight2_2 = OutType(0.5)*root3;
  OutType scale = sqrt((2.0f*numTriangles)/(weight1_1*weight2_2 - weight2_1*weight1_2));
  weight1_1 *= scale;
  weight1_2 *= scale;
  weight2_1 *= scale;
  weight2_2 *= scale;

  weightedDeterminant = (weight1_1*weight2_2 - weight2_1*weight1_2);
  auto cross = vtkm::Cross(TriEdges[0], TriEdges[2]);
  norm = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
  if(weightedDeterminant == 0.0 || norm == 0.0) 
  {
    return OutType(0.0);
  }
  weightedJacobian = norm/weightedDeterminant;
  result = vtkm::Min(weightedJacobian, 1/weightedJacobian);
  result *= result;

if (result > 0)
  {
      return vtkm::Min(result, vtkm::Infinity<OutType>()); //normal case
  }
  return vtkm::Max(result, vtkm::NegativeInfinity<OutType>());
}

// Compute the RelativeSizeSquared of a quadrilateral.
// Formula: (min(R, 1/R))^2
//     - R = A/~A
//     - ~A = average area of all triangles in the mesh
// value for unit triangle dependent on ~A
// Acceptable range: [0.3,1]
// Normal range:     [0,1]
// Full range:       [0,1]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeRelativeSizeSquared(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagQuad,
			      const vtkm::Id& numQuads,
                              const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("RelativeSizeSquared metric(quad) requires 4 points.");
    return OutType(0.0);
  }
  OutType negativeInfinity = vtkm::NegativeInfinity<OutType>();
  OutType area = 0;
  OutType weight1_1, weight2_1, weight1_2, weight2_2, avgArea, weightedJacobian, result; 
  using Edge = typename PointCoordVecType::ComponentType;
  const Edge QuadEdges[4] = {pts[1] - pts[0],
			     pts[2] - pts[1],
                             pts[3] - pts[2],
   		             pts[0] - pts[3] 
	      		    };
  //get area
  Edge principle1 = QuadEdges[0]+(pts[2]-pts[3]);
  Edge principle2 = QuadEdges[3]+QuadEdges[1];
  auto pAxis = vtkm::Cross(principle1, principle2);
  auto pNorm = pAxis/vtkm::MagnitudeSquared(pAxis);
  for(int i=0; i<3;i++)
  {
    auto currNorm = vtkm::Cross(QuadEdges[3+i%4], QuadEdges[i]);
    area += vtkm::Dot(pNorm, currNorm); 
  }
  area *= OutType(0.25);
  //get weights
  weight1_1 = weight2_2 = 1;
  weight2_1 = weight1_2 = 0;
  OutType scale = sqrt(numQuads/(weight1_1*weight2_2 - weight2_1*weight1_2));
  weight1_1 *= scale; 
  weight1_2 *= scale; 
  weight2_2 *= scale; 
  weight2_1 *= scale; 
  avgArea = weight1_1*weight2_2 - weight1_2*weight2_1;
  if(avgArea <= negativeInfinity) return OutType(0.0);
  if(area <= negativeInfinity) return OutType(0.0);
  weightedJacobian = area/avgArea;
  result = vtkm::Min(weightedJacobian, 1/weightedJacobian);
  result *= result;
if (result > 0)
  {
      return vtkm::Min(result, vtkm::Infinity<OutType>()); //normal case
  }
  return vtkm::Max(result, negativeInfinity);
}

// ============================= 3D Volume cells ==================================
// Compute the RelativeSizeSquared of a tetrahedron.
// Formula: q = min(R, 1/R)^2
//	- R = V/~V
//	- ~V = average volume of the tetrahedrons in the mesh
//	- if ~V < vtkm::NegativeInfinity, q = 0
//	- if R <= vtkm::NegativeInfinity, q = 0
// Value for equilateral tetra is unddefined(depends on ~V)
// Acceptable Range: [0.3, 1]
// Normal range: [0, 1]
// Full range: [0, 1]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeRelativeSizeSquared(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagTetra,
			      const vtkm::Id& numTets,
                              const vtkm::exec::FunctorBase& worklet)
{ 
  if (numPts != 4)
  {
    worklet.RaiseError("RelativeSizeSquared metric(tetrahedron) requires 4 points.");
    return OutType(0.0);
  }

  //The 6 edges of a tetrahedron
  const OutType negativeInfinity = vtkm::NegativeInfinity<OutType>();
  const float rt3 = sqrt(3.0f);
  const float rt2 = sqrt(2.0f);
  using Edge = typename PointCoordVecType::ComponentType;
  const Edge TetEdges[12] = {pts[1] - pts[0],   // 0
  	       		     pts[2] - pts[1], 
  	       		     pts[0] - pts[2],
  	       		     pts[3] - pts[0],   // 3
  	       		     pts[3] - pts[1],
  	       		     pts[3] - pts[2]
			    };
  // add weights to calculate the average volume
  Edge weight1 = {1, 0, 0};
  Edge weight2 = {0.5f, 0.5f*rt3, 0};
  Edge weight3 = {0.5f, rt3/6.0f, rt2/rt3};
  OutType avg_volume, volume, result, R;
  float scale = (float)pow(6*numTets/ vtkm::Dot(weight1, vtkm::Cross(weight2,weight3)), 0.333333333333333333);
  weight1 = weight1*scale;
  weight2 = weight2*scale;
  weight3 = weight3*scale;
  //calculate average volume and volume
  avg_volume = (vtkm::Dot(weight1, vtkm::Cross(weight2, weight3))/6);
  // calculate volume
  volume = vtkm::Dot(TetEdges[3], vtkm::Cross(TetEdges[2], TetEdges[0]))/6;
  if(avg_volume < negativeInfinity)
    return OutType(0.0);
  
  R = volume/avg_volume;
  if(R < negativeInfinity)
    return OutType(0.0);

  result = vtkm::Min(R, 1/R);
  result *= result;
  if (result > 0)
  {
      return vtkm::Min(result, vtkm::Infinity<OutType>()); //normal case
  }
  return vtkm::Max(result, negativeInfinity);
}
// Compute the RelativeSizeSquared of a hexahedron.
// Formula:
// Equals 0 for a unit cube
// Acceptable Range: [0, 0.5]
// Normal range: [0,FLOAT_MAX]
// Full range: [0,FLOAT_MAX]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeRelativeSizeSquared(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagHexahedron,
			      const vtkm::Id& numHexs,
                              const vtkm::exec::FunctorBase& worklet)
{ 
  if (numPts != 8)
  {
    worklet.RaiseError("RelativeSizeSquared metric(hexahedron) requires 8 points.");
    return OutType(0.0);
  }
  OutType avg_volume;
  const OutType negativeInfinity = vtkm::NegativeInfinity<OutType>();
  //The 12 edges of a hexahedron
  using Edge = typename PointCoordVecType::ComponentType;
  const Edge HexEdges[12] = {pts[1] - pts[0],   // 0
  	       		     pts[2] - pts[1], 
  	       		     pts[3] - pts[2],
  	       		     pts[3] - pts[0],   // 3
  	       		     pts[4] - pts[0],
  	       		     pts[5] - pts[1],
  	       		     pts[6] - pts[2],   // 6
  	       		     pts[7] - pts[3],
  	       		     pts[5] - pts[4],
  	       		     pts[6] - pts[5],   // 9
  	       		     pts[7] - pts[6],
  	       		     pts[7] - pts[4]    // 11
	          	    };

  Edge principleXAxis = HexEdges[0] + (pts[2] - pts[3]) + HexEdges[8] + (pts[6] - pts[7]);
  Edge principleYAxis = (pts[3] - pts[0]) + HexEdges[1] + (pts[7] - pts[4]) + HexEdges[9];
  Edge principleZAxis = HexEdges[4] + HexEdges[5] + HexEdges[6] + HexEdges[7];
  Edge hexJacobianMatrices[9][3] = { {HexEdges[0], HexEdges[3], HexEdges[4]},
                             {HexEdges[1], (-1*HexEdges[0]), HexEdges[5]},
                             {HexEdges[2], (-1*HexEdges[1]), HexEdges[6]},
                             {(-1*HexEdges[3]), (-1*HexEdges[2]), HexEdges[7]},
                             {HexEdges[11], HexEdges[8], (-1*HexEdges[4])},
                             {(-1*HexEdges[8]), HexEdges[9], (-1*HexEdges[5])},
                             {(-1*HexEdges[9]), HexEdges[10], (-1*HexEdges[6])},
                             {(-1*HexEdges[10]), (-1*HexEdges[11]), (-1*HexEdges[7])},
                             {principleXAxis, principleYAxis, principleZAxis}
                           };
  Edge weight1 = {1, 0, 0};
  Edge weight2 = {0, 1, 0};
  Edge weight3 = {0, 0, 1};
  
  float scale = (float)pow(numHexs/(vtkm::Dot(weight1, vtkm::Cross(weight2, weight3))), 0.3333333333333333);
  weight1 = weight1*scale;
  weight2 = weight2*scale;
  weight3 = weight3*scale;
  avg_volume = vtkm::Dot(weight1, vtkm::Cross(weight2, weight3));

  OutType tempMatrix1_1, tempMatrix1_2, tempMatrix1_3, tempMatrix2_2, tempMatrix2_3, tempMatrix3_3, determinant;
  OutType firstPartOfNumerator, secondPartOfNumerator, currentRelativeSizeSquared, maxRelativeSizeSquared = negativeInfinity;
  for (vtkm::IdComponent matrixNumber = 0; matrixNumber < 9; matrixNumber ++)
  {
/*
// these computations equal the value at X_Y for the matrix B,
// where B = matrix A multiplied by its transpose.
// the values at X_X are also used for the matrix multiplication AA.
// Note that the values 1_2 = 2_1, 1_3 = 3_1, and 2_3 = 3_2.
// This fact is used to optimize the computation
*/  
    tempMatrix1_1 = vtkm::Dot(hexJacobianMatrices[matrixNumber][0], hexJacobianMatrices[matrixNumber][0]);
    tempMatrix1_2 = vtkm::Dot(hexJacobianMatrices[matrixNumber][0], hexJacobianMatrices[matrixNumber][1]);
    tempMatrix1_3 = vtkm::Dot(hexJacobianMatrices[matrixNumber][0], hexJacobianMatrices[matrixNumber][2]);
    tempMatrix2_2 = vtkm::Dot(hexJacobianMatrices[matrixNumber][1], hexJacobianMatrices[matrixNumber][1]);
    tempMatrix2_3 = vtkm::Dot(hexJacobianMatrices[matrixNumber][1], hexJacobianMatrices[matrixNumber][2]);
    tempMatrix3_3 = vtkm::Dot(hexJacobianMatrices[matrixNumber][2], hexJacobianMatrices[matrixNumber][2]);
    determinant   = vtkm::Dot(hexJacobianMatrices[matrixNumber][0], vtkm::Cross(hexJacobianMatrices[matrixNumber][1], hexJacobianMatrices[matrixNumber][2]));
    if(determinant > vtkm::NegativeInfinity<OutType>())
    {
      firstPartOfNumerator =   (tempMatrix1_1 * tempMatrix1_1) +
			     2*(tempMatrix1_2 * tempMatrix1_2) + 
			     2*(tempMatrix1_3 * tempMatrix1_3) + 
			       (tempMatrix2_2 * tempMatrix2_2) +
			     2*(tempMatrix2_3 * tempMatrix2_3) + 
			       (tempMatrix3_3 * tempMatrix3_3);
      secondPartOfNumerator = tempMatrix1_1 + tempMatrix2_2 + tempMatrix3_3;
      secondPartOfNumerator *= secondPartOfNumerator;
      OutType temp_third = 0.0;
      secondPartOfNumerator *= temp_third;
      currentRelativeSizeSquared = (firstPartOfNumerator - secondPartOfNumerator)/(vtkm::Pow(determinant, OutType(4.0*temp_third)));
      maxRelativeSizeSquared = currentRelativeSizeSquared > maxRelativeSizeSquared ? currentRelativeSizeSquared : maxRelativeSizeSquared;
    }
    else {return vtkm::Infinity<OutType>(); }
  }
  if (maxRelativeSizeSquared > 0)
  {
      return vtkm::Min(maxRelativeSizeSquared, vtkm::Infinity<OutType>()); //normal case
  }
  return vtkm::Max(maxRelativeSizeSquared, negativeInfinity);
}

} // namespace cellmetrics
} // namespace exec
} // namespace vtkm

#endif // vtk_m_exec_cellmetrics_RelativeSizeSquared_h
