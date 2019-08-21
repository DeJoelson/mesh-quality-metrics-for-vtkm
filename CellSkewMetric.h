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
#ifndef vtk_m_exec_CellSkewMetric_h
#define vtk_m_exec_CellSkewMetric_h
/*
 */

#include "vtkm/CellShape.h"
#include "vtkm/CellTraits.h"
#include "vtkm/VecTraits.h"
#include "vtkm/VectorAnalysis.h"
#include "vtkm/exec/FunctorBase.h"
#include "vtkm/exec/cellmetrics/CellConditionMetric.h"

namespace vtkm
{
namespace exec
{
static constexpr FloatType FLOAT_MAX = vtkm::Infinity<FloatType>();
static constexpr FloatType FLOAT_MIN = vtkm::NegativeInfinity<FloatType>();

template <typename OutType, typename PointCoordVecType, typename CellShapeType>
VTKM_EXEC OutType CellStretchMetric(const vtkm::IdComponent& numPts,
                                    const PointCoordVecType& pts,
                                    CellShapeType shape,
                                    const vtkm::exec::FunctorBase& worklet)
{
    worklet.RaiseError("Shape type template must be specified to compute skew")
    return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellSkewMetric(const vtkm::IdComponents &numPts,
                                  const PointCoordVecType& pts,
                                  vtkm::CellShapeHex,
                                  const vtkm::exec::FunctorBase& worklet)
{
    FloatType X1 = (pts[1] - pts[0]) + (pts[2] - pts[3]) + (pts[5] - pts[4]) + (pts[6] - pts[7]);
    FloatType X1_mag = vtkm::Magnitude(X1);
    if (X1_mag < FLOAT_MIN)
        return FLOAT_MAX;
    FloatType x1 = X1/X1_mag;
    FloatType X2 = (pts[3] - pts[0]) + (pts[2] - pts[1]) + (pts[7] - pts[4]) + (pts[6] - pts[5]);
    FloatType X2_mag = vtkm::Magnitude(X2);
    if (X2_mag < FLOAT_MIN)
        return FLOAT_MAX;
    FloatType x2 = X2/X2_mag;
    FloatType X3 = (pts[4] - pts[0]) + (pts[5] - pts[1]) + (pts[6] - pts[2]) + (pts[7] - pts[3]);
    FloatType X3_mag = vtkm::Magnitude(X3);
    if (X3_mag < FLOAT_MIN)
        return FLOAT_MAX;
    FloatType x3 = X3/X3_mag;
    return vtkm::Max(vtkm::Dot(x1,x2),vtkm::Max(vtkm::Dot(x1,x3),vtkm::Dot(x2,x3)));
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellSkewMetric(const vtkm::IdComponents &numPts,
                                  const PointCoordVecType& pts,
                                  vtkm::CellShapeQuad,
                                  const vtkm::exec::FunctorBase& worklet)
{
    FloatType X1 = (pts[1] - pts[0]) + (pts[2] - pts[3]);
    FloatType X1_mag = vtkm::Magnitude(X1);
    if (X1_mag < FLOAT_MIN)
        return 0;
    FloatType x1 = X1/X1_mag;
    FloatType X2 = (pts[2] - pts[1]) + (pts[3] - pts[0]);
    FloatType X2_mag = vtkm::Magnitude(X2);
    if (X2_mag < FLOAT_MIN)
        return 0;
    FloatType x2 = X2/X2_mag;
    return vtkm::Abs(vtkm::Dot(x1,x2));

}

} // exec
} // vtkm

#endif vtk_m_exec_CellSkew_Metric_h

