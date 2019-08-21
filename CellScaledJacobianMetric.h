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
#ifndef vtk_m_exec_cellmetrics_ScaledJacobian_h
#define vtk_m_exec_cellmetrics_ScaledJacobian_h

/*
 * Mesh quality metric functions that computes the scaled jacobian of mesh cells.
 * The jacobian of a cell is defined as the determinant of the Jociabian matrix 
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
VTKM_EXEC OutType ComputeScaledJacobian(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              CellShapeType shape,
                              const vtkm::exec::FunctorBase&)
{
  UNUSED(numPts);  
  UNUSED(pts);  
  UNUSED(shape);  
  return OutType(0.0);
}

//TODO: Should polygons be supported? Maybe call Quad or Triangle function...
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeScaledJacobian(const vtkm::IdComponent& numPts,
                                 const PointCoordVecType& pts,
                                 vtkm::CellShapeTagPolygon,
                                 const vtkm::exec::FunctorBase& worklet)
{
  switch(numPts)
  {
    case 3:
            return ComputeScaledJacobian<OutType>(numPts, pts, vtkm::CellShapeTagTriangle(), worklet);
    case 4:
            return ComputeScaledJacobian<OutType>(numPts, pts, vtkm::CellShapeTagQuad(), worklet);
    default:
  	    return OutType(-2.0);
  }
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeScaledJacobian(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagLine,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-2.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeScaledJacobian(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagWedge,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-2.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeScaledJacobian(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagPyramid,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-2.0);
}
// ========================= 2D cells ==================================
//Compute the scaled jacobian of a triangle
//Formula: q = ((2*sqrt(3))/3) * (J/Lmax)
//	- J -> jacobian, if surface normal N is center of triangle and
//	     and N*L2*L1 < 0, then -jacobian 
//	- Lmax -> max{ |L0| * |L1|, |L0| * |L2|, |L1| * |L2| }
//Equals 1 for equilateral unit triangle
//Acceptable Range: [0.5, 2*sqrt(3)/3]
//Normal Range    : [ -(2*sqrt(3)/3) , 2*sqrt(3)/3]
//Full Range      : [-FLOAT_MAX, FLOAT_MAX]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeScaledJacobian(const vtkm::IdComponent& numPts,
                                 const PointCoordVecType& pts,
                                 vtkm::CellShapeTagTriangle,
                                 const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 3)
  {
    worklet.RaiseError("ScaledJacobian metric(tri) requires 3 points.");
    return OutType(0.0);
  }

  //The 3 edges of a triangle 
  using Edge = typename PointCoordVecType::ComponentType;
  Edge TriEdges[3] = {pts[2] - pts[1],
	              pts[0] - pts[2],
		      pts[1] - pts[0]
                     };
  Edge TriCross = vtkm::Cross((TriEdges[1] - TriEdges[0]), (TriEdges[2] - TriEdges[0]));
  OutType constant  = (OutType)(2.0/vtkm::Sqrt(3.0));
  OutType scaledJacobian = vtkm::Magnitude(TriCross);     
  scaledJacobian *= constant;

  OutType maxEdgeLengthProduct = vtkm::Max(vtkm::Max( vtkm::Magnitude(TriEdges[0])*vtkm::Magnitude(TriEdges[1]),
						      vtkm::Magnitude(TriEdges[1])*vtkm::Magnitude(TriEdges[2])
 						    ), 
					   vtkm::Magnitude(TriEdges[0])*vtkm::Magnitude(TriEdges[2])
 					  );
  scaledJacobian /= maxEdgeLengthProduct;

  //TODO: compute surface normal and adjust scaledjacobian accordingly
  auto surfaceNormalAtCenter = vtkm::Normal((TriEdges[0] + TriEdges[1] + TriEdges[2])/3);
  if (vtkm::Dot(surfaceNormalAtCenter, TriCross) < 0)
  {
    scaledJacobian *= -1;
  }

  if (scaledJacobian > 0)
    return vtkm::Min(scaledJacobian, vtkm::Infinity<OutType>()); //normal case

  return vtkm::Max(scaledJacobian, vtkm::NegativeInfinity<OutType>());
}

// Compute the scaled jacobian of a quadrilateral.
// Formula: min{J0/(L0*L3), J1/(L1*L0), J2/(L2*L1), J3/(L3*L2)}
//	-Ji -> Jacobian at corner i, the intersection of the edge vectors
//	       it is divided by
// Equals 1 for a unit square
//Acceptable Range: [0.3, 1]
//Normal Range    : [-1, 1]
// Full range     : [-1, 1]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeScaledJacobian(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagQuad,
                              const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("ScaledJacobian metric(quad) requires 4 points.");
    return OutType(0.0);
  }

  //The 4 edges of a quadrilateral
  using Edge = typename PointCoordVecType::ComponentType;
  Edge QuadEdges[4] = {pts[1] - pts[0],
	                     pts[2] - pts[1],
	                     pts[3] - pts[2],
	                     pts[0] - pts[3]
	              };

  OutType QuadEdgeLengths[4] = {vtkm::Magnitude(QuadEdges[0]),
			        vtkm::Magnitude(QuadEdges[1]),
			        vtkm::Magnitude(QuadEdges[2]),
			        vtkm::Magnitude(QuadEdges[3])
               	               };
  OutType negativeInfinity = vtkm::NegativeInfinity<OutType>();
  if (QuadEdgeLengths[0] < negativeInfinity || QuadEdgeLengths[1] < negativeInfinity || QuadEdgeLengths[2] < negativeInfinity || QuadEdgeLengths[2] < negativeInfinity)
  {
    return OutType(0);
  }
/*
3 * 0
0 * 1
1 * 2
2 * 3
*/
  
  //Calculate Jacobian & scaled jacobian of each corner, storing and returning the min scaled jacobian

  //cross product of principal axes
  Edge centerCalculation = vtkm::Cross(QuadEdges[0] - (pts[2] - pts[3]), QuadEdges[1] - (pts[3] - pts[0]));
  vtkm::Normalize(centerCalculation);
  
  OutType currJacobian, currScaledJacobian, minScaledJacobian = vtkm::Infinity<OutType>();

  currJacobian = vtkm::Dot(vtkm::Cross(QuadEdges[3], QuadEdges[0]), centerCalculation);
  currScaledJacobian = currJacobian / (QuadEdgeLengths[0] * QuadEdgeLengths[3]);
  minScaledJacobian = currScaledJacobian < minScaledJacobian ? currScaledJacobian : minScaledJacobian;
  
  currJacobian = vtkm::Dot(vtkm::Cross(QuadEdges[0], QuadEdges[1]), centerCalculation);
  currScaledJacobian = currJacobian / (QuadEdgeLengths[1] * QuadEdgeLengths[0]);
  minScaledJacobian = currScaledJacobian < minScaledJacobian ? currScaledJacobian : minScaledJacobian;
  
  currJacobian = vtkm::Dot(vtkm::Cross(QuadEdges[1], QuadEdges[2]), centerCalculation);
  currScaledJacobian = currJacobian / (QuadEdgeLengths[2] * QuadEdgeLengths[1]);
  minScaledJacobian = currScaledJacobian < minScaledJacobian ? currScaledJacobian : minScaledJacobian;
  
  currJacobian = vtkm::Dot(vtkm::Cross(QuadEdges[2], QuadEdges[3]), centerCalculation);
  currScaledJacobian = currJacobian / (QuadEdgeLengths[3] * QuadEdgeLengths[2]);
  minScaledJacobian = currScaledJacobian < minScaledJacobian ? currScaledJacobian : minScaledJacobian;
  
  if (minScaledJacobian > 0)
    return vtkm::Min(minScaledJacobian, vtkm::Infinity<OutType>()); //normal case

  return vtkm::Max(minScaledJacobian, negativeInfinity);
}

// ============================= 3D Volume cells ==================================
// Compute the scaled jacobian of a hexahedron.
// Formula: q = min{Ai}
//	Ai -> for i 1...8 (Jacobian determinant at respective corner, divided by corresponding edge lengths
// Equals 1 for a unit cube
// Acceptable Range: [0.5, 1]
// Normal Range    : [-1, 1]
// Full range      : [1,FLOAT_MAX]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeScaledJacobian(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagHexahedron,
                              const vtkm::exec::FunctorBase& worklet)
{ 
  if (numPts != 8)
  {
    worklet.RaiseError("ScaledJacobian metric(hexahedron) requires 8 points.");
    return OutType(0.0);
  }

  //The 12 edges of a hexahedron
  using Edge = typename PointCoordVecType::ComponentType;
  Edge HexEdges[12] = {pts[1] - pts[0],
  		       pts[2] - pts[1], 
       		       pts[3] - pts[2], 
       		       pts[3] - pts[0], 
       		       pts[4] - pts[0], 
       		       pts[5] - pts[1], 
       		       pts[6] - pts[2], 
       		       pts[7] - pts[3], 
       		       pts[5] - pts[4], 
       		       pts[6] - pts[5], 
       		       pts[7] - pts[6], 
       		       pts[7] - pts[4]
            	      };
  Edge principleXAxis = HexEdges[0] + (pts[2] - pts[3]) + HexEdges[8] + (pts[6] - pts[7]);
  Edge principleYAxis = (pts[3] - pts[0]) + HexEdges[1] + (pts[7] - pts[4]) + HexEdges[9];
  Edge principleZAxis = HexEdges[4] + HexEdges[5] + HexEdges[6] + HexEdges[7];
  Edge hexMatrices[9][3] = { {HexEdges[0], HexEdges[3], HexEdges[4]},
                             {HexEdges[1], (-1*HexEdges[0]), HexEdges[5]},
                             {HexEdges[2], (-1*HexEdges[1]), HexEdges[6]},
                             {(-1*HexEdges[3]), (-1*HexEdges[2]), HexEdges[7]},
                             {HexEdges[11], HexEdges[8], (-1*HexEdges[4])},
                             {(-1*HexEdges[8]), HexEdges[9], (-1*HexEdges[5])},
                             {(-1*HexEdges[9]), HexEdges[10], (-1*HexEdges[6])},
                             {(-1*HexEdges[10]), (-1*HexEdges[11]), (-1*HexEdges[7])},
                             {principleXAxis, principleYAxis, principleZAxis}
                           };
  OutType currDeterminant, minDeterminant = vtkm::Infinity<OutType>();
  FloatType lenSquared1, lenSquared2, lenSquared3, minLengthSquared = vtkm::Infinity<FloatType>();
  vtkm::IdComponent matrixIndex;
  for (matrixIndex = 0; matrixIndex < 9; matrixIndex++)
  {
    lenSquared1 = (FloatType)vtkm::MagnitudeSquared(hexMatrices[matrixIndex][0]);
    minLengthSquared = lenSquared1 < minLengthSquared ? lenSquared1 : minLengthSquared;
    lenSquared2 = (FloatType)vtkm::MagnitudeSquared(hexMatrices[matrixIndex][1]);
    minLengthSquared = lenSquared2 < minLengthSquared ? lenSquared2 : minLengthSquared;
    lenSquared3 = (FloatType)vtkm::MagnitudeSquared(hexMatrices[matrixIndex][2]);
    minLengthSquared = lenSquared3 < minLengthSquared ? lenSquared3 : minLengthSquared;

    vtkm::Normalize(hexMatrices[matrixIndex][0]);
    vtkm::Normalize(hexMatrices[matrixIndex][1]);
    vtkm::Normalize(hexMatrices[matrixIndex][2]);
    currDeterminant = (OutType)vtkm::Dot(hexMatrices[matrixIndex][0], vtkm::Cross(hexMatrices[matrixIndex][1], hexMatrices[matrixIndex][2]));
    if (currDeterminant < minDeterminant)
    {
	minDeterminant = currDeterminant;
    }
  }
  if (minLengthSquared < vtkm::NegativeInfinity<FloatType>())
  {
    return vtkm::Infinity<OutType>();
  }
  OutType toReturn = minDeterminant;
  if (toReturn > 0)
    return vtkm::Min(toReturn, vtkm::Infinity<OutType>()); //normal case

  return vtkm::Max(toReturn, vtkm::NegativeInfinity<OutType>());
}

// Compute the scaled jacobian of a tetrahedron.
// Formula: q = J*sqrt(2)/Lamda_max
//	J -> jacobian,(L2 * L0) * L3
//	Lamda_max -> max{ L0*L2*L3, L0*L1*L4, L1*L2*L5, L3*L4*L5}
// Equals Sqrt(2) / 2 for unit equilateral tetrahedron
// Acceptable Range: [0, FLOAT_MAX]
// Normal Range: [0, FLOAT_MAX]
// Full range: [FLOAT_MIN,FLOAT_MAX]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeScaledJacobian(const vtkm::IdComponent& numPts,
                                 const PointCoordVecType& pts,
                                 vtkm::CellShapeTagTetra,
                                 const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("ScaledJacobian metric requires 4 points");
    return OutType(0.0);
  }
  
  //the edge and side sets 
  using Edge = typename PointCoordVecType::ComponentType;
  Edge Edges[6] = {pts[1] - pts[0],
                         pts[2] - pts[1],
                         pts[0] - pts[2],
                         pts[3] - pts[0],
                         pts[3] - pts[1],
                         pts[3] - pts[2]
			};
  OutType EdgesSquared[6];
  OutType jacobian = vtkm::Dot(vtkm::Cross(Edges[2], Edges[0]), Edges[3]);
  // compute the scaled jacobian
  OutType currSide, maxSide = vtkm::NegativeInfinity<OutType>();
  vtkm::IdComponent edgeIndex, sideIndex;
  for (edgeIndex = 0; edgeIndex < 6; edgeIndex++)
  {
    EdgesSquared[edgeIndex] = vtkm::MagnitudeSquared(Edges[edgeIndex]);
  }
  OutType Sides[4] = {EdgesSquared[0]*EdgesSquared[2]*EdgesSquared[3],
                         EdgesSquared[0]*EdgesSquared[1]*EdgesSquared[4],
                         EdgesSquared[1]*EdgesSquared[2]*EdgesSquared[5],
                         EdgesSquared[3]*EdgesSquared[4]*EdgesSquared[5]
			};
  for (sideIndex = 0; sideIndex < 4; sideIndex++)
  {
    currSide = Sides[sideIndex];
    maxSide = currSide > maxSide ? currSide : maxSide;
  }
  maxSide = vtkm::Sqrt(maxSide);
  OutType toUseInCalculation = jacobian > maxSide ? jacobian : maxSide;
  if (toUseInCalculation < vtkm::NegativeInfinity<OutType>())
  {
    return vtkm::Infinity<OutType>();
  }  
  return (vtkm::Sqrt<OutType>(2)*jacobian)/toUseInCalculation;
}

} // namespace cellmetrics
} // namespace exec
} // namespace vtkm

#endif // vtk_m_exec_cellmetrics_CellEdgeRatioMetric_h
