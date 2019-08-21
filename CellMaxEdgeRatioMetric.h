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
#ifndef vtk_m_exec_cellmetrics_CellMaxEdgeRatioMetric_h
#define vtk_m_exec_cellmetrics_CellMaxEdgeRatioMetric_h

/*
 * Mesh quality metric functions that compute the maximum edge ratio of mesh cells.
 * The maximum edge ratio of a cell is defined as the maximum of each pairwise combination 
 * division of principle axes.
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


template<typename OutType, typename VecType>
VTKM_EXEC inline OutType ComputeMaxEdgeRatio(const VecType& axes)
{
  const vtkm::Id numAxes = axes.GetNumberOfComponents();
  //Compare axis ratios to determine the maximum 

  FloatType  numerator, denominator;
  OutType currentRatio, maxRatio = vtkm::NegativeInfinity<OutType>();
  vtkm::IdComponent  numeratorIndex, denominatorIndex;
  
  for (numeratorIndex = 0; numeratorIndex < numAxes; numeratorIndex++)
  {
    numerator = (FloatType)vtkm::MagnitudeSquared(axes[numeratorIndex]);
    {
      for (denominatorIndex = 0; denominatorIndex < numAxes; denominatorIndex++)
      {
        if (numeratorIndex != denominatorIndex)
        {
          denominator = (FloatType)vtkm::MagnitudeSquared(axes[denominatorIndex]);
          currentRatio =(OutType)  vtkm::Sqrt(numerator / denominator);
	  maxRatio = maxRatio > currentRatio ?  maxRatio : currentRatio;
        }
      }
    } 
  } 


  if (maxRatio > 0)
    return vtkm::Min(maxRatio, vtkm::Infinity<OutType>()); //normal case

  return vtkm::Max(maxRatio, OutType(-1)*vtkm::Infinity<OutType>());
}  


// ========================= Unsupported cells ==================================

// By default, cells have zero shape unless the shape type template is specialized below.
template <typename OutType, typename PointCoordVecType, typename CellShapeType>
VTKM_EXEC OutType CellMaxEdgeRatioMetric(const vtkm::IdComponent& numPts,
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
VTKM_EXEC OutType CellMaxEdgeRatioMetric(const vtkm::IdComponent& numPts,
                                 const PointCoordVecType& pts,
                                 vtkm::CellShapeTagPolygon,
                                 const vtkm::exec::FunctorBase& worklet)
{
  switch (numPts)
  {
    case 4:
            return CellMaxEdgeRatioMetric<OutType>(numPts, pts, vtkm::CellShapeTagQuad(), worklet);
    default:
            break;
  } 
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMaxEdgeRatioMetric(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagLine,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMaxEdgeRatioMetric(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagTriangle,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMaxEdgeRatioMetric(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagTetra,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMaxEdgeRatioMetric(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagWedge,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMaxEdgeRatioMetric(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagPyramid,
                                 const vtkm::exec::FunctorBase& worklet)
{
  UNUSED(worklet);
  return OutType(-1.0);
}
// ========================= 2D cells ==================================
// Compute the max edge ratio of a quadrilateral.
// Formula: max{X1/X2, X2/X1}
//     - X<N> -> the principle axis of dimension N (X1 = x-axis, X2 = y-axis, ...)
// Equals 1 for a unit square
// Acceptable range: [1, 1.3]
// Normal Range: [1, DBL_MAX]
// Full range: [1,FLOAT_MAX]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMaxEdgeRatioMetric(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagQuad,
                              const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("Maximum edge ratio metric(quad) requires 4 points.");
    return OutType(0.0);
  }
 
  vtkm::IdComponent numAxes = 2; 

  //The 2 principle axes of a quadrilateral
  using AxesType = typename PointCoordVecType::ComponentType;
  const AxesType QuadAxes[2] = {(pts[1] - pts[0]) + (pts[2] - pts[3]),
				(pts[2] - pts[1]) + (pts[3] - pts[0]) 
	      		       };

  return vtkm::exec::cellmetrics::ComputeMaxEdgeRatio<OutType>(vtkm::make_VecC(QuadAxes, numAxes));
}

// ============================= 3D Volume cells ==================================
// Compute the max edge ratio of a hexahedron.
// Formula: max{A12, A13, A23}
//     - A<N1><N2> = max{XN1/XN2, XN2/XN1}
//     - X<N> -> the principle axis of dimension N (X1 = x-axis, X2 = y-axis, ...)
// Equals 1 for a unit cube
// Acceptable Range: [1 ,1.3] 
// Normal Range: [1, FLOAT_MAX]
// Full range: [1,FLOAT_MAX]
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellMaxEdgeRatioMetric(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagHexahedron,
                              const vtkm::exec::FunctorBase& worklet)
{ 
  if (numPts != 8)
  {
    worklet.RaiseError("Max edge ratio metric(hexahedron) requires 8 points.");
    return OutType(0.0);
  }

  vtkm::IdComponent numAxes = 3; 

  //The 3 axes of a hexahedron
  using AxesType = typename PointCoordVecType::ComponentType;
  const AxesType HexAxes[3] = {pts[1] - pts[0] + pts[2] - pts[3] + pts[5] - pts[4] + pts[6] - pts[7],
                           pts[3] - pts[0] + pts[2] - pts[1] + pts[7] - pts[4] + pts[6] - pts[5],
                           pts[4] - pts[0] + pts[5] - pts[1] + pts[6] - pts[2] + pts[7] - pts[3]
	      			       };

  return vtkm::exec::cellmetrics::ComputeMaxEdgeRatio<OutType>(vtkm::make_VecC(HexAxes, numAxes));
}

} // namespace cellmetrics
} // namespace exec
} // namespace vtkm

#endif // vtk_m_exec_cellmetrics_CellEdgeRatioMetric_h
