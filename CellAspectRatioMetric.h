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
#ifndef vtk_m_exec_cellmetrics_CellAspectRatioMetric_h
#define vtk_m_exec_cellmetrics_CellAspectRatioMetric_h

/*
 * Mesh quality metric functions that compute the aspect ratio of mesh cells.
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
#include "TypeOfCellTriangle.h"
#include <iostream>

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
VTKM_EXEC OutType CellAspectRatioMetric(const vtkm::IdComponent& numPts,
                                          const PointCoordVecType& pts,
                                          CellShapeType shape,
                                          const vtkm::exec::FunctorBase&)
{
  std::cout << "An Unsupported Cell was passed; numPts: " << numPts << std::endl; //<< numPts << "; shape: " << shape << std::endl;
  UNUSED(numPts);
  UNUSED(pts);
  UNUSED(shape);
  return OutType(0);
}

// ========================= 2D cells ==================================
// Compute the diagonal ratio of a triangle.
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellAspectRatioMetric(const vtkm::IdComponent& numPts,
                                          const PointCoordVecType& pts,
                                          vtkm::CellShapeTagTriangle,
                                          const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 3)
  {
    worklet.RaiseError("Aspect ratio metric (triangle) requires 3 points.");
    return OutType(0.0);
  }

  /*
  q is the quantity of our metric for a triangle. It is given by:
  
  q = (lmax * (l0 + l1 + l2))/(4*sqrt(3)*A)
  
  Where
  
  lmax = max(l0, l1, l2)
  A = 0.5 * ||L0xL1||
  
  l0 = ||L0||
  l1 = ||L1||
  l2 = ||L2||
  
  ||vector|| implies the $\ell2$-norm (i.e., magnitude) of ``vector''.
  
  L0 = P2-P1
  L1 = P0-P2
  L2 = P1-P0
  
  
  Therefore:
  q = (lmax * (l0 + l1 + l2))/(4 * sqrt(3)* A)
  q = (lmax * (l0 + l1 + l2))/(4 * sqrt(3)* 0.5 * ||L0xL1||)
  q = (lmax * (l0 + l1 + l2))/(2 * sqrt(3) * ||L0xL1||)
  q = (lmax * (l0 + l1 + l2) * 0.5)/(sqrt(3) * ||L0xL1||)
  q = (lmax * (l0 + l1 + l2) * 0.5 * rsqrt(3))/(||L0xL1||)
  q = (lmax * (l0 + l1 + l2) * k)/(||L0xL1||); where k = 0.5 * rsqrt(3) = 0.28867513459
  q = (lmax * (l0 + l1 + l2) * k)/(c); where c = ||L0xL1||
  */

  using Scalar = OutType;
  using CollectionOfPoints = PointCoordVecType;
  using Vector = typename PointCoordVecType::ComponentType;
  const Vector TriangleEdges[3] = {pts[2] - pts[1],
                                 pts[0] - pts[2],
                                 pts[1] - pts[0]
                                };
  const Vector L0 = GetL0 <Scalar, Vector, CollectionOfPoints> (pts);
  const Vector L1 = TriangleEdges[1];
  const Vector L2 = TriangleEdges[2];  
  
  const OutType l0 = vtkm::Sqrt(vtkm::MagnitudeSquared(L0));
  const OutType l1 = vtkm::Sqrt(vtkm::MagnitudeSquared(L1));
  const OutType l2 = vtkm::Sqrt(vtkm::MagnitudeSquared(L2));
  const OutType lmax = vtkm::Max(l0, vtkm::Max(l1, l2));
  const OutType crossProductNormResult = vtkm::Sqrt(vtkm::MagnitudeSquared(vtkm::Cross(L0, L1)));
  const OutType k(0.28867513459);
  OutType q = (lmax * (l0 + l1 + l2) * k)/(crossProductNormResult);
  
  return q;
}

// Compute the diagonal ratio of a quad.
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellAspectRatioMetric(const vtkm::IdComponent& numPts,
                                          const PointCoordVecType& pts,
                                          vtkm::CellShapeTagQuad,
                                          const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("Aspect ratio metric (quad) requires 4 points.");
    return OutType(0.0);
  }
  
  /*
    q = (lmax * (l0+l1+l2+l3)) / (4 * A)
    
    Where
    
    lmax = max(l0, l1, l2, l3)
    
    l0 = ||L0||
    l1 = ||L1||
    l2 = ||L2||
    l3 = ||L2||
    
    ||vector|| implies the $\ell2$-norm (i.e., magnitude) of ``vector''.
  
    L0 = P1-P0
    L1 = P2-P1
    L2 = P3-P2
    L3 = P0-P3
    
    A = 0.25 * (a0 + a1 + a2 + a3)
    
    a0 = nc \cdot N0
    a1 = nc \cdot N1
    a2 = nc \cdot N2
    a3 = nc \cdot N3
    
    NC = (X1 \times X2)/||(X1 \times X2)||     
    
    
    X1 = (P1 - P0) + (P2 - P3)
    X2 = (P2 - P1) + (P3 - P0)
    
    N0 = L3 \times L0
    N1 = L0 \times L1
    N2 = L1 \times L2
    N3 = L2 \times L3
    
    
    Therefore:
    q = (lmax * (l0 + l1 + l2 + l3)) / (4 * A)
    q = (lmax * (l0 + l1 + l2 + l3)) / (4 * 0.25 * (a0 + a1 + a2 + a3))
    q = (lmax * (l0 + l1 + l2 + l3)) / (a0 + a1 + a2 + a3)
    
  
  */
    using Vector = typename PointCoordVecType::ComponentType;
    using Scalar = OutType;
    
  const Vector QuadEdges[4] = { pts[1] - pts[0],
			      pts[2] - pts[1],
			      pts[3] - pts[2],
			      pts[0] - pts[3]
		    	    };
  const Vector L0 = QuadEdges[0];
  const Vector L1 = QuadEdges[1];
  const Vector L2 = QuadEdges[2];
  const Vector L3 = QuadEdges[3];
  
  const Scalar l0 = vtkm::Sqrt(vtkm::MagnitudeSquared(L0));
  const Scalar l1 = vtkm::Sqrt(vtkm::MagnitudeSquared(L1));
  const Scalar l2 = vtkm::Sqrt(vtkm::MagnitudeSquared(L2));
  const Scalar l3 = vtkm::Sqrt(vtkm::MagnitudeSquared(L3));
  
  const Scalar lmax = vtkm::Max(l0, vtkm::Max(l1, vtkm::Max(l2, l3)));
  
  const Vector X1 = QuadEdges[0] + (pts[2] - pts[3]);
  const Vector X2 = QuadEdges[1] + (pts[3] - pts[0]);
  
  const Vector N0 = vtkm::Cross(L3, L0);
  const Vector N1 = vtkm::Cross(L0, L1);
  const Vector N2 = vtkm::Cross(L1, L2);
  const Vector N3 = vtkm::Cross(L2, L3);
  
  const Vector Xcross = vtkm::Cross(X1, X2);
  const Vector NC = Xcross/vtkm::Sqrt(vtkm::MagnitudeSquared(Xcross));
  
  const Scalar a0 = vtkm::Dot(NC, N0);
  const Scalar a1 = vtkm::Dot(NC, N1);
  const Scalar a2 = vtkm::Dot(NC, N2);
  const Scalar a3 = vtkm::Dot(NC, N3);
  
  const Scalar q = (lmax * (l0 + l1 + l2 + l3)) / (a0 + a1 + a2 + a3); 
  
  return q;
}

// ============================= 3D Volume cells ==================================
// Compute the aspect ratio of a tetrahedron.
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellAspectRatioMetric(const vtkm::IdComponent& numPts,
                                          const PointCoordVecType& pts,
                                          vtkm::CellShapeTagTetra,
                                          const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("Aspect ratio metric(tetrahedron) requires 4 points.");
    return OutType(0.0);
  }

  std::cout << "A Tetra Cell was passed" << std::endl;
  return OutType(13.31);
}

} // namespace cellmetrics
} // namespace exec
} // namespace vtkm

#endif // vtk_m_exec_cellmetrics_CellEdgeRatioMetric_h
