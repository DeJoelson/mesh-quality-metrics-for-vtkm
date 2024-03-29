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
#ifndef vtk_m_exec_cellmetrics_Oddy_h
#define vtk_m_exec_cellmetrics_Oddy_h

/*
 * Mesh quality metric functions that compute the Oddy of mesh cells.
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
VTKM_EXEC OutType ComputeOddy(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              CellShapeType shape,
                              const vtkm::exec::FunctorBase&)
{
  UNUSED(numPts);  
  UNUSED(pts);  
  UNUSED(shape);  
  return OutType(0.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeOddy(const vtkm::IdComponent& numPts,
                                 const PointCoordVecType& pts,
                                 vtkm::CellShapeTagPolygon,
                                 const vtkm::exec::FunctorBase& worklet)
{
  switch(numPts)
  {
    case 4:
	    return ComputeOddy<OutType>(numPts, pts, vtkm::CellShapeTagQuad(), worklet);
    default:
	    break;
  }
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeOddy(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagLine,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeOddy(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagTriangle,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeOddy(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagTetra,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeOddy(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagWedge,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeOddy(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagPyramid,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}
// ========================= 2D cells ==================================
// Compute the Oddy of a quadrilateral.
// Formula: for i 0 to 3: max{[(||Li||^2 - ||Li+1||^2)^2 + 4((Li * Li+1)^2)] / (2||Ni+1||^2)} 
//     - L4 = L0
//     - '*' symbolizes the dot product of two vectors
//     - Ni is the normal vector associated with each point
// Equals 0 for a unit square
// Acceptable range: [0,0.5]
// Normal range: [0,FLOAT_MAX]
// Full range: [0,FLOAT_MAX]
// Note, for loop avoided because L4 = L0
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeOddy(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagQuad,
                              const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("Oddy metric(quad) requires 4 points.");
    return OutType(0.0);
  }
 
  using Edge = typename PointCoordVecType::ComponentType;
  const Edge QuadEdges[4] = {pts[1] - pts[0],
			     pts[2] - pts[1],
                             pts[3] - pts[2],
   		             pts[0] - pts[3] 
	      		    };
  OutType negativeInfinity = vtkm::NegativeInfinity<OutType>();
  OutType oddyNumerator, oddyDenominator, currComputation, maxComputation = negativeInfinity;

  for (vtkm::IdComponent edgeIndex = 0; edgeIndex < 4; edgeIndex++)
  {
    oddyNumerator =(OutType) vtkm::MagnitudeSquared(vtkm::MagnitudeSquared(QuadEdges[edgeIndex]) - vtkm::MagnitudeSquared(QuadEdges[(edgeIndex+1%4)])) + 4*vtkm::MagnitudeSquared(vtkm::Dot(QuadEdges[edgeIndex], QuadEdges[(edgeIndex+1%4)]));
    oddyDenominator = (OutType) vtkm::MagnitudeSquared(vtkm::Cross(QuadEdges[(edgeIndex+3%4)], QuadEdges[edgeIndex]));
    if (oddyDenominator < vtkm::NegativeInfinity<OutType>()){return vtkm::Infinity<OutType>();}
    oddyDenominator *= OutType(2);
    currComputation = oddyNumerator/oddyDenominator;
    maxComputation = currComputation > maxComputation ? currComputation : maxComputation;
  }

if (maxComputation > 0)
  {
      return vtkm::Min(maxComputation, vtkm::Infinity<OutType>()); //normal case
  }
  return vtkm::Max(maxComputation, negativeInfinity);
}

// ============================= 3D Volume cells ==================================
// Compute the Oddy of a hexahedron.
// Formula:
// Equals 0 for a unit cube
// Acceptable Range: [0, 0.5]
// Normal range: [0,FLOAT_MAX]
// Full range: [0,FLOAT_MAX]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeOddy(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagHexahedron,
                              const vtkm::exec::FunctorBase& worklet)
{ 
  if (numPts != 8)
  {
    worklet.RaiseError("Oddy metric(hexahedron) requires 8 points.");
    return OutType(0.0);
  }

  //The 12 points of a hexahedron
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

  OutType third = (OutType) 1/3;
  OutType negativeInfinity = vtkm::NegativeInfinity<OutType>();
  OutType tempMatrix1_1, tempMatrix1_2, tempMatrix1_3, tempMatrix2_2, tempMatrix2_3, tempMatrix3_3, determinant;
  OutType firstPartOfNumerator, secondPartOfNumerator, currentOddy, maxOddy = negativeInfinity;
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
      secondPartOfNumerator *= third;
      currentOddy = (firstPartOfNumerator - secondPartOfNumerator)/(vtkm::Pow(determinant, OutType(4.0*third)));
      maxOddy = currentOddy > maxOddy ? currentOddy : maxOddy;
    }
    else {return vtkm::Infinity<OutType>(); }
  }
  if (maxOddy > 0)
  {
      return vtkm::Min(maxOddy, vtkm::Infinity<OutType>()); //normal case
  }
  return vtkm::Max(maxOddy, negativeInfinity);
}

} // namespace cellmetrics
} // namespace exec
} // namespace vtkm

#endif // vtk_m_exec_cellmetrics_Oddy_h
