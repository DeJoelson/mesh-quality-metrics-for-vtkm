//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2014 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//  Copyright 2014 UT-Battelle, LLC.
//  Copyright 2014 Los Alamos National Security.
//
//  Under the terms of Contract DE-NA0003525 with NTESS,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================
#ifndef vtk_m_exec_CellConditionMetric_h
#define vtk_m_exec_CellConditionMetric_h
/*
 * Mesh quality metric functions that compute the shape, or weighted Jacobian, of mesh cells.
 * The Jacobian of a cell is weighted by the condition metric value of the cell.
 *
 * These metric computations are adapted from the VTK implementation of the Verdict library,
 * which provides a set of cell metrics for evaluating the geometric qualities of regions of mesh spaces.
 * 
 * See: The Verdict Library Reference Manual (for per-cell-type metric formulae)
 * See: vtk/ThirdParty/verdict/vtkverdict (for VTK code implementation of this metric)
 */

#include "vtkm/CellShape.h"
#include "vtkm/CellTraits.h"
#include "vtkm/VecTraits.h"
#include "vtkm/VectorAnalysis.h"
#include "vtkm/exec/FunctorBase.h"

namespace vtkm
{
namespace exec
{
namespace cellmetrics
{

using FloatType = vtkm::FloatDefault;

static constexpr FloatType rt3 = vtkm::Sqrt(3.0f);
static constexpr FloatType rt6 = vtkm::Sqrt(6.0f);
static constexpr FloatType FLOAT_MAX = vtkm::Infinity<FloatType>();
static constexpr FloatType FLOAT_MIN = vtkm::NegativeInfinity<FloatType>();

// By default, cells have zero shape unless the shape type template is specialized below.
template <typename OutType, typename PointCoordVecType, typename CellShapeType>
VTKM_EXEC OutType CellConditionMetric(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              CellShapeType shape,
                              const vtkm::exec::FunctorBase&)
{
  return OutType(0.0);
}

// ========================= Shapeless cells (?) ==================================

//TODO: What is the shape for a Line? (Magnitude/edge length may be only metric of a Line)
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellConditionMetric(const vtkm::IdComponent& numPts,
                                  const PointCoordVecType& pts,
                                  vtkm::CellShapeTagLine,
                                  const vtkm::exec::FunctorBase& worklet)
{
  if (numPts < 2)
  {
    worklet.RaiseError("Degenerate line has no shape.");
    return OutType(0.0);
  }
  OutType shape(0.0);
  return shape;
}


//TODO: What is the shape of a pyramid (sum of shape metrics of 2 Tets that make up pyramid? or just 0?)
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellConditionMetric(const vtkm::IdComponent& numPts,
                                  const PointCoordVecType& pts,
                                  vtkm::CellShapeTagPyramid,
                                  const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 5)
  {
    worklet.RaiseError("Shape metric(pyramid) requires 5 points.");
    return OutType(0.0);
  }
  OutType shape(0.0);
  return shape;
}


// =============================== Shape metric cells ==================================

// Compute the shape quality metric of a triangular cell.
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellConditionMetric(const vtkm::IdComponent& numPts,
                                      const PointCoordVecType& pts,
                                      vtkm::CellShapeTagTriangle,
                                      const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 3)
  {
    worklet.RaiseError("Condition metric(triangle) requires 3 points.");
    return OutType(0.0);
  }

  typename PointCoordVecType::ComponentType v1 = pts[1] - pts[0];
  typename PointCoordVecType::ComponentType v2 = pts[2] - pts[0];
  OutType area = OutType(0.5) * static_cast<OutType>(Magnitude(Cross(v1, v2)));
  return area;
}

/// Compute the area of a quadrilateral cell.
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellConditionMetric(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagQuad,
                              const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("Area(quad) requires 4 points.");
    return OutType(0.0);
  }

  typename PointCoordVecType::ComponentType edges[4] = {
    pts[1] - pts[0], pts[2] - pts[1], pts[3] - pts[2], pts[0] - pts[3],
  };

  typename PointCoordVecType::ComponentType cornerNormals[4] = {
    Cross(edges[3], edges[0]),
    Cross(edges[0], edges[1]),
    Cross(edges[1], edges[2]),
    Cross(edges[2], edges[3]),
  };

  // principal axes
  typename PointCoordVecType::ComponentType principalAxes[2] = {
    edges[0] - edges[2], edges[1] - edges[3],
  };

  // Unit normal at the quadrilateral center
  typename PointCoordVecType::ComponentType unitCenterNormal =
    Cross(principalAxes[0], principalAxes[1]);
  Normalize(unitCenterNormal);

  OutType area =
    (Dot(unitCenterNormal, cornerNormals[0]) + Dot(unitCenterNormal, cornerNormals[1]) +
     Dot(unitCenterNormal, cornerNormals[2]) + Dot(unitCenterNormal, cornerNormals[3])) *
    OutType(0.25);
  return area;
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType ComputeMeasure(const vtkm::IdComponent&,
                                 const PointCoordVecType&,
                                 vtkm::CellShapeTagPolygon,
                                 const vtkm::exec::FunctorBase& worklet)
{
  return OutType(0.0);
}


// ============================= 3D volumetric cells ==================================
/// Compute the condition metric of a tetrahedron.
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellConditionMetric(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagTetra,
                              const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 4)
  {
    worklet.RaiseError("Condition metric(tetrahedron) requires 4 points.");
    return OutType(0.0);
  }

  using FloatType = vtkm::FloatDefault;

  
  

  typename PointCoordVecType::ComponentType v1 = pts[1] - pts[0];
  typename PointCoordVecType::ComponentType v2 = pts[2] - pts[0];
  typename PointCoordVecType::ComponentType v3 = pts[3] - pts[0];
  OutType volume = Dot(Cross(v1, v2), v3) / OutType(6.0);
  return volume;
}

/// Compute the volume of a hexahedral cell (approximated via triple product of average edge along each parametric axis).
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellConditionMetric(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagHexahedron,
                              const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 8)
  {
    worklet.RaiseError("Volume(hexahedron) requires 8 points.");
    return OutType(0.0);
  }

  auto efg1 = pts[1];
  efg1 += pts[2];
  efg1 += pts[5];
  efg1 += pts[6];
  efg1 -= pts[0];
  efg1 -= pts[3];
  efg1 -= pts[4];
  efg1 -= pts[7];

  auto efg2 = pts[2];
  efg2 += pts[3];
  efg2 += pts[6];
  efg2 += pts[7];
  efg2 -= pts[0];
  efg2 -= pts[1];
  efg2 -= pts[4];
  efg2 -= pts[5];

  auto efg3 = pts[4];
  efg3 += pts[5];
  efg3 += pts[6];
  efg3 += pts[7];
  efg3 -= pts[0];
  efg3 -= pts[1];
  efg3 -= pts[2];
  efg3 -= pts[3];

  OutType volume = Dot(Cross(efg2, efg3), efg1) / OutType(64.0);
  return volume;
}

/// Compute the volume of a wedge cell (approximated as 3 tetrahedra).
template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellConditionMetric(const vtkm::IdComponent& numPts,
                              const PointCoordVecType& pts,
                              vtkm::CellShapeTagWedge,
                              const vtkm::exec::FunctorBase& worklet)
{
  if (numPts != 6)
  {
    worklet.RaiseError("Volume(wedge) requires 6 points.");
    return OutType(0.0);
  }

  typename PointCoordVecType::ComponentType v0 = pts[1] - pts[0];
  typename PointCoordVecType::ComponentType v1 = pts[2] - pts[0];
  typename PointCoordVecType::ComponentType v2 = pts[3] - pts[0];
  OutType volume = Dot(Cross(v0, v1), v2) / OutType(6.0);

  typename PointCoordVecType::ComponentType v3 = pts[4] - pts[1];
  typename PointCoordVecType::ComponentType v4 = pts[5] - pts[1];
  typename PointCoordVecType::ComponentType v5 = pts[3] - pts[1];
  volume += Dot(Cross(v3, v4), v5) / OutType(6.0);

  typename PointCoordVecType::ComponentType v6 = pts[5] - pts[1];
  typename PointCoordVecType::ComponentType v7 = pts[2] - pts[1];
  typename PointCoordVecType::ComponentType v8 = pts[3] - pts[1];
  volume += Dot(Cross(v6, v7), v8) / OutType(6.0);

  return volume;
}

} //end of namespace cellmetrics
}
}

#endif // vtk_m_exec_CellConditionMetric_h
