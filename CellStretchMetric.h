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
#ifndef vtk_m_exec_cellmetrics_CellStretchMetric_h
#define vtk_m_exec_cellmetrics_CellStretchMetric_h

/*
 * Mesh quality metric functions that compute the aspect frobenius of certain mesh cells.
 * The aspect frobenius metric generally measures the degree of regularity of a cell, with
 * a value of 1 representing a regular cell..
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

static constexpr FloatType FLOAT_MAX = vtkm::Infinity<FloatType>();
static constexpr FloatType FLOAT_MIN = vtkm::NegativeInfinity<FloatType>();

template <typename OutType, typename PointCoordVecType, typename CellShapeType>
VTKM_EXEC OutType CellStretchMetric(const vtkm::IdComponent& numPts,
                                    const PointCoordVecType& pts,
                                    CellShapeType shape,
                                    const vtkm::exec::FunctorBase& worklet)
{
    worklet.RaiseError("Shape type template must be specified to compute stretch")
    return OutType(-1.0);
}

template<typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellStretchMetric(const vtkm::IdComponent& numPts,
                                    const PointCoordVecType& pts,
                                    vtlm::CellShapeTagQuad,
                                    const vtkm::exec::FunctorBase& worklet)
{
    FloatType L0 = vtkm::MagnitudeSquared(pts[1] - pts[0]);
    FloatType L1 = vtkm::MagnitudeSquared(pts[2] - pts[1]);
    FloatType L2 = vtkm::MagnitudeSquared(pts[3] - pts[2]);
    FloatType L3 = vtkm::MagnitudeSquared(pts[0] - pts[3]);
    // Find the minimum length (use square of values to speed up)
    FloatType D0 = pts[2] - pts[0];
    FloatType D1 = pts[3] - pts[1];
    FloatType D_max = vtkm::Max(D0,D1);
    if(D_max < FLOAT_MIN) return FLOAT_MAX;

    return vtkm::Sqrt(2)/D_max * vtkm::Sqrt(L_min);
}

template<typename OutType, typename PointCoordVecType>
VTKM_EXEC OutType CellStretchMetric(const vtkm::IdComponent& numPts,
                                    const PointCoordVecType& pts,
                                    vtlm::CellShapeTagHex,
                                    const vtkm::exec::FunctorBase& worklet)
{
    FloatType L0 = vtkm::MagnitudeSquared(pts[1] - pts[0]);
    FloatType L1 = vtkm::MagnitudeSquared(pts[2] - pts[1]);
    FloatType L2 = vtkm::MagnitudeSquared(pts[3] - pts[2]);
    FloatType L3 = vtkm::MagnitudeSquared(pts[3] - pts[0]);
    FloatType L4 = vtkm::MagnitudeSquared(pts[4] - pts[0]);
    FloatType L5 = vtkm::MagnitudeSquared(pts[5] - pts[1]);
    FloatType L6 = vtkm::MagnitudeSquared(pts[6] - pts[2]);
    FloatType L7 = vtkm::MagnitudeSquared(pts[3] - pts[3]);
    FloatType L8 = vtkm::MagnitudeSquared(pts[5] - pts[4]);
    FloatType L9 = vtkm::MagnitudeSquared(pts[6] - pts[5]);
    FloatType L10 = vtkm::MagnitudeSquared(pts[7] - pts[6]);
    FloatType L11 = vtkm::MagnitudeSquared(pts[7] - pts[4]);
    FloatType D0 = pts[6] - pts[0];
    FloatType D1 = pts[7] - pts[1];
    FloatType D2 = pts[4] - pts[2];
    FloatType D3 = pts[5] - pts[3];
    FloatType D_max = vtkm::Max(D0,D1,D2,D3);
    if(D_max < FLOAT_MIN) return FLOAT_MAX;
    return vtkm::Sqrt(3) * vtkm::Sqrt(vtkm::Min(L0,L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11)) / D_max;
}


#endif  // vtk_m_exec_cellmetrics_CellStretchMetric_h
