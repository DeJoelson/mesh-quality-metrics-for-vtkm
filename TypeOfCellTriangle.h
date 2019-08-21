/**
 * The Verdict manual defines a set of commonly
 * used components of a triangle. For example,
 * area, side lengths, and so forth.
 *
 * This file contains a set of functions which
 * implement return the values of those commonly
 * used components for subsequent use in metrics.
 */


/**
 * Returns the L0 vector, as defined by the verdict manual.
 *
 *  \param [in] pts The three points which define the triangle.
 *  \return Returns the vector.
 */
template <typename Scalar, typename Vector, typename CollectionOfPoints>
VTKM_EXEC Vector GetL0(const CollectionOfPoints& pts)
{    
    const Vector L0(pts[2] - pts[1]);
    return L0;  
}

/**
 * Returns the L1 vector, as defined by the verdict manual.
 *
 *  \param [in] pts The three points which define the triangle.
 *  \return Returns the vector.
 */
template <typename Scalar, typename Vector, typename CollectionOfPoints>
VTKM_EXEC Vector GetL1(const CollectionOfPoints& pts)
{    
    const Vector L1(pts[0] - pts[2]);
    return L0;  
}

/**
 * Returns the L0 vector, as defined by the verdict manual.
 *
 *  \param [in] pts The three points which define the triangle.
 *  \return Returns the vector.
 */
template <typename Scalar, typename Vector, typename CollectionOfPoints>
VTKM_EXEC Vector GetL2(const CollectionOfPoints& pts)
{    
    const Vector L2(pts[1] - pts[0]);
    return L2;  
}


