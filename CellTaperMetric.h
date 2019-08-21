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
#ifndef vtk_m_exec_CellTaperMetric_h
#define vtk_m_exec_CellTaperMetric_h
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
    worklet.RaiseError("Shape type template must be specified to compute taper")
    return OutType(-1.0);
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellTaperMetric(const vtkm::IdComponents &numPts,
                                  const PointCoordVecType& pts,
                                  vtkm::CellShapeHex,
                                  const vtkm::exec::FunctorBase& worklet)
{
    FloatType X1 = (pts[1] - pts[0]) + (pts[2] - pts[3]) + (pts[5] - pts[4]) + (pts[6] - pts[7]);
    FloatType X2 = (pts[3] - pts[0]) + (pts[2] - pts[1]) + (pts[7] - pts[4]) + (pts[6] - pts[5]);
    FloatType X3 = (pts[4] - pts[0]) + (pts[5] - pts[1]) + (pts[6] - pts[2]) + (pts[7] - pts[3]);
    if ((X1 < FLOAT_MIN) || (X2 < FLOAT_MIN) || (X3 < FLOAT_MIN))
	return FLOAT_MAX;
    FloatType X12 = ((pts[2] - pts[3]) - (pts[1] - pts[0])) + ((pts[6] - pts[7]) - (pts[5] - pts[4]));
    FloatType X13 = ((pts[5] - pts[1]) - (pts[4] - pts[0])) + ((pts[6] - pts[2]) - (pts[7] - pts[3]));
    FloatType X23 = ((pts[7] - pts[4]) - (pts[3] - pts[0])) + ((pts[6] - pts[5]) - (pts[2] - pts[1]));
    FloatType T12 = X12 / vtkm::Min(X1,vtkm::Min(X2,X3));
    FloatType T13 = X13 / vtkm::Min(X1,vtkm::Min(X2,X3));
    FloatType T23 = X23 / vtkm::Min(X1,vtkm::Min(X2,X3));
    return vtkm::Max(T12, vtkm::Max(T13,T23));
}

template <typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellTaperMetric(const vtkm::IdComponents &numPts,
                                  const PointCoordVecType& pts,
                                  vtkm::CellShapeQuad,
                                  const vtkm::exec::FunctorBase& worklet)
{
    FloatType X1 = (pts[1] - pts[0]) + (pts[2] - pts[3]);
    FloatType X2 = (pts[0] - pts[3]) + (pts[2] - pts[1]);
    if ((X1 < FLOAT_MIN) || (X2 < FLOAT_MIN)) 
	return FLOAT_MAX;
    FloatType X12 = vtkm::Magnitude((pts[0] - pts[1]) + (pts[2] - pts[3]));
    return X12 / vtkm::Min(X1,X2);
}

} // exec
} // vtkm

#endif vtk_m_exec_CellTaper_Metric_h
