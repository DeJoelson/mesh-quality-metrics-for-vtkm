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
#ifndef vtk_m_exec_CellWarpageMetric_h
#define vtk_m_exec_CellWarpageMetric_h
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
    worklet.RaiseError("Shape type template must be Quad to compute warpage")
    return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellWarpageMetric(const vtkm::IdComponents &numPts,
                                  const PointCoordVecType& pts,
                                  vtkm::CellShapeQuad,
                                  const vtkm::exec::FunctorBase& worklet)
{
    FloatType N0 = vtkm::Cross((pts[0]-pts[3]),(pts[1]-pts[0]));
    FloatType N1 = vtkm::Cross((pts[1]-pts[0]),(pts[2]-pts[1]));
    FloatTyoe N2 = vtkm::Cross((pts[2]-pts[1]),(pts[3]-pts[2]));
    FloatType N3 = vtkm::Cross((pts[3]-pts[2]),(pts[0]-pts[3]));
    FloatType N0_mag = vtkm::Magnitude(N0);
    FloatType N1_mag = vtkm::Magnitude(N1);
    FloatType N2_mag = vtkm::Magnitude(N2);
    FloatType N3_mag = vtkm::Magnitude(N3);
    if ((N0_mag < FLOAT_MIN) || (N1_mag < FLOAT_MIN) || (N2_mag < FLOAT_MIN) || (N3_mag < FLOAT_MIN))
        return FLOAT_MAX;
    FloatType n0 = N0/N0_mag;
    FloatType n1 = N1/N1_mag;
    FloatType n2 = N2/N2_mag;
    FloatType n3 = N3/N3_mag;
    return 1. - vtkm::Min(vtkm::Pow(vtkm::Dot(n0,n2),3),vtkm::Pow(vtkm::Dot(n1,n3),3));
}

} // exec
} // vtkm

#endif vtk_m_exec_CellWarpage_Metric_h

