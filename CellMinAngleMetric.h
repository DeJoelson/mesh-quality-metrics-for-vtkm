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
#ifndef vtk_m_exec_cellmetrics_CellMinAngleMetric_h
#define vtk_m_exec_cellmetrics_CellMinAngleMetric_h

/*
 * Mesh quality metric functions that compute the minimum angle of cell in a mesh.
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
VTKM_EXEC OutType CellMinAngleMetric(const vtkm::IdComponent& numPts,
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
VTKM_EXEC OutType CellMinAngleMetric(const vtkm::IdComponent& numPts,
                                 const PointCoordVecType& pts,
                                 vtkm::CellShapeTagPolygon,
                                 const vtkm::exec::FunctorBase& worklet)
{
  switch (numPts)
  {
    case 3:
            return CellMinAngleMetric<OutType>(numPts, pts, vtkm::CellShapeTagTriangle(), worklet);
    case 4:
            return CellMinAngleMetric<OutType>(numPts, pts, vtkm::CellShapeTagQuad(), worklet);
    default:
            break;
  } 
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMinAngleMetric(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagLine,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}


template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMinAngleMetric(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagHexahedron,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMinAngleMetric(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagWedge,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMinAngleMetric(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagPyramid,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}
// ========================= 2D cells ==================================
// Compute the minimum angle of a triangle.
// Formula: q = min( arccos((Ln dot Ln+1)/(||Ln|| * ||Ln+1||))(180º/π) for n 0,1, and 2 )
//     - L3 = L0
//     - if any edge has length 0, return q = 360º
//     - All angle measurements are in degrees
// q equals 60 for a unit triangle
// Acceptable range: [30º, 60º]
// Normal Range: [0º, 360º]
// Full range: [0º, 360º]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMinAngleMetric(const vtkm::IdComponent& numPts,
                                 const PointCoordVecType& pts,
                                 vtkm::CellShapeTagTriangle,
                                 const vtkm::exec::FunctorBase& worklet)
{
  if(numPts != 3)
  {
    worklet.RaiseError("Minimum angle metric(quad) requires 3 points.");
    return OutType(0.0);
  }
  using Edge = typename PointCoordVecType::ComponentType;
  const Edge TriEdges[3] = { pts[2] - pts[1],
			     pts[0] - pts[2],
			     pts[1] - pts[0]
			   };
  const OutType negativeInfinity = vtkm::NegativeInfinity<OutType>();
  const OutType infinity = vtkm::Infinity<OutType>();
  const double scalar = 57.2957795131;   // ~ 180/pi
  OutType result = infinity, tempResult;
  Edge edge1, edge2;
  int ctr, ctrPlusOne;
  for(ctr = 0; ctr < 3; ctr++)
  {
    if(TriEdges[ctr] == OutType(0.0))
    {
      return OutType(360);
    }
  }

  for(ctr = 0; ctr < 3; ctr++)
  {
    ctrPlusOne = (ctr+1)%4; // wrap around to first edge
    edge1 = TriEdges[ctr];
    edge2 = TriEdges[ctrPlusOne];
    tempResult = 1/vtkm::Cos(vtkm::Dot(edge1, edge2)/(vtkm::MagnitudeSquared(edge1)*vtkm::MagnitudeSquared(edge2)));
    tempResult *= scalar;
    result = tempResult < result ? tempResult : result;
  }
  
  if (result > 0)
    return vtkm::Min(result, infinity); //normal case

  return vtkm::Max(result, negativeInfinity);
}
// Compute the min angle of a quadrilateral.
// Formula: q = min( Ai for i 0,1,2, and 3 )
//     - L4 = L0
//     - Ai = -1^Si arccos(-1(Li dot Li+1)/(||Li||||Li+1||) )(180/π) + 360º*Si 
//     - if ||Li|| <= FLOAT_MIN or ||Li+1|| <= FLOAT_MIN, return q = 360º
// q = 90º for a unit square
// Acceptable range: [45º, 90º]
// Normal Range: [0º, 90º]
// Full range: [0º, 360º]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMinAngleMetric(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagQuad,
                              const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("Minimum angle metric(quad) requires 4 points.");
    return OutType(0.0);
  }
 
  using Edge = typename PointCoordVecType::ComponentType;
  const Edge QuadEdges[4] = { pts[1] - pts[0],
			      pts[2] - pts[1],
			      pts[3] - pts[2],
			      pts[0] - pts[3]
		    	    };
  const Edge principleXAxis = QuadEdges[0] + (pts[2] - pts[3]);
  const Edge principleYAxis = QuadEdges[1] + (pts[3] - pts[0]);
  Edge centerNormal = vtkm::Cross(principleXAxis, principleYAxis);
  vtkm::Normalize(centerNormal);
  const OutType negativeInfinity = vtkm::NegativeInfinity<OutType>();
  const OutType infinity = vtkm::Infinity<OutType>();
  const double scalar = 57.2957795131;   // ~ 180/pi
  OutType result = infinity, tempResult, tempNormal;
  Edge edge1, edge2;
  int ctr, ctrPlusOne;
  for(ctr = 0; ctr < 3; ctr++)
  {
    if(vtkm::MagnitudeSquared(QuadEdges[ctr]) <= negativeInfinity)
    {
      return OutType(360);
    }
  }
  Edge temp;
  for(ctr = 0; ctr < 3; ctr++)
  {
    ctrPlusOne = (ctr+1)%4; // wrap around to first edge
    edge1 = QuadEdges[ctr];
    edge2 = QuadEdges[ctrPlusOne];
    temp = vtkm::Cross(QuadEdges[(ctr+3)%4], edge1);
    temp = temp / vtkm::MagnitudeSquared(temp);
    tempNormal = vtkm::Dot(centerNormal, temp);
    tempResult = 1/vtkm::Cos(-1*vtkm::Dot(edge1, edge2)/(vtkm::MagnitudeSquared(edge1)*vtkm::MagnitudeSquared(edge2)));
    tempResult *= scalar;
    if(tempNormal < 0)
    {
      tempResult *= -1;
      tempResult += 360;
    }
    result = tempResult < result ? tempResult : result;
  }
  
  if (result > 0)
    return vtkm::Min(result, infinity); //normal case

  return vtkm::Max(result, negativeInfinity);
}

// ============================= 3D Volume cells ==================================
// Compute the max edge ratio of a tetrahedron.
// Formula: min{ Ai for i 0,1,2,3,4,5}
//     - Ai = 180º * arccos(ni1 dot ni2)
//     - ni1 & ni2 are the vector normal of the two tetrahedron faces adjacent to edge i
// For an equilateral tetrahedron, q = (180º/π) arccos(1/3) = ~70.528779
// Acceptable Range: [40º, (180º/π) arccos(1/3)]
// Normal Range: [0º, (180º/π) arccos(1/3)]
// Full range: [0º, 360º]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMinAngleMetric(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagTetra,
                              const vtkm::exec::FunctorBase& worklet)
{ 
  if (numPts != 4)
  {
    worklet.RaiseError("Max edge ratio metric(tetrahedron) requires 8 points.");
    return OutType(0.0);
  }


  //The 3 axes of a hexahedron
  using Edge = typename PointCoordVecType::ComponentType;
  const Edge TetEdges[6] = { pts[1]-pts[0],
			     pts[2]-pts[1],
			     pts[0]-pts[2],
			     pts[3]-pts[0],
			     pts[3]-pts[1],
			     pts[3]-pts[2]
	      		   };
  const OutType negativeInfinity = vtkm::NegativeInfinity<OutType>();
  const OutType infinity = vtkm::Infinity<OutType>();
  const double scalar = 57.2957795131;   // ~ 180/pi
  OutType result = infinity;
  Edge edge1 = vtkm::Cross(TetEdges[0], TetEdges[1]);
  Edge edge2 = vtkm::Cross(TetEdges[0], TetEdges[4]);
  Edge edge3 = vtkm::Cross(TetEdges[4], TetEdges[5]);
  Edge edge4 = vtkm::Cross(TetEdges[1], TetEdges[5]);
  OutType n1 = vtkm::MagnitudeSquared(edge1);
  OutType n2 = vtkm::MagnitudeSquared(edge2);
  OutType n3 = vtkm::MagnitudeSquared(edge3);
  OutType n4 = vtkm::MagnitudeSquared(edge4);

  OutType tempResult = 1/vtkm::Cos(vtkm::Dot(edge1, edge2)/(n1*n2)); 
  result = tempResult < result ? tempResult : result;

  tempResult = 1/vtkm::Cos(vtkm::Dot(edge1, edge3)/(n1*n3)); 
  result = tempResult < result ? tempResult : result;

  tempResult = 1/vtkm::Cos(vtkm::Dot(edge1, edge4)/(n1*n4)); 
  result = tempResult < result ? tempResult : result;

  tempResult = 1/vtkm::Cos(vtkm::Dot(edge2, edge3)/(n2*n3)); 
  result = tempResult < result ? tempResult : result;

  tempResult = 1/vtkm::Cos(vtkm::Dot(edge2, edge4)/(n2*n4)); 
  result = tempResult < result ? tempResult : result;

  tempResult = 1/vtkm::Cos(vtkm::Dot(edge3, edge4)/(n3*n4)); 
  result = tempResult < result ? tempResult : result;

  result *= scalar;

  if (result > 0)
    return vtkm::Min(result, infinity); //normal case

  return vtkm::Max(result, negativeInfinity);
}

} // namespace cellmetrics
} // namespace exec
} // namespace vtkm

#endif // vtk_m_exec_cellmetrics_CellEdgeRatioMetric_h
