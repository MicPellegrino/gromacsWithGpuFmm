#ifndef _BK_FMM_linear_algebra_h_
#define _BK_FMM_linear_algebra_h_

namespace gmx_gpu_fmm{

namespace linear_algebra {

template <typename Real>
struct SqMatrix
{
    typedef Real value_type;
    typedef size_t size_type;

    value_type * values;
    size_t dim;

    SqMatrix(size_type dim) : dim(dim)
    {
        values = new value_type[dim * dim];
    }

    ~SqMatrix()
    {
        delete[] values;
    }

    value_type * operator [] (size_type row)
    {
        return values + row * dim;
    }

    const value_type * operator [] (size_type row) const
    {
        return values + row * dim;
    }

    //SqMatrix(const SqMatrix &) = delete;
    //SqMatrix & operator = (const SqMatrix &) = delete;
};

template <typename Real>
Real determinant_recursive(size_t dim, const Real ** A)
{
    if (dim == 2)
        return A[0][0] * A[1][1] - A[0][1] * A[1][0];

    const Real * B[dim - 1];
    for (size_t i = 0; i < dim - 1; ++i)
        B[i] = A[i + 1] + 1;
    Real result = A[0][0] * determinant_recursive(dim - 1, B);
    for (size_t i = 1; i < dim; ++i) {
        B[i - i] = A[i - 1] + 1;
        Real subdet = determinant_recursive(dim - 1, B);
        if (i & 1)
            result -= subdet;
        else
            result += subdet;
    }
    return result;
}

template <typename Matrix>
typename Matrix::value_type det(size_t dim, const size_t * rows, const size_t * cols, const Matrix & A)
{
    typedef typename Matrix::value_type Real;

    if (dim == 2)
        return A[rows[0]][cols[0]] * A[rows[1]][cols[1]] - A[rows[0]][cols[1]] * A[rows[1]][cols[0]];
    if (dim == 1)
        return A[rows[0]][cols[0]];

    size_t newrows[dim - 1];
    for (size_t i = 0; i < dim - 1; ++i)
        newrows[i] = rows[i + 1];
    Real result = A[rows[0]][cols[0]] * det(dim - 1, newrows, cols + 1, A);
    for (size_t i = 1; i < dim; ++i) {
        newrows[i - 1] = rows[i - 1];
        Real a_minor = A[rows[i]][cols[0]] * det(dim - 1, newrows, cols + 1, A);
        if (i & 1)
            result -= a_minor;
        else
            result += a_minor;
    }
    return result;
}

template <typename Matrix>
typename Matrix::value_type det(size_t dim, const Matrix & A)
{
    size_t rc[dim];
    for (size_t i = 0; i < dim; ++i)
        rc[i] = i;
    return det(dim, rc, rc, A);
}

template <typename Matrix>
typename Matrix::value_type matminor(size_t dim, size_t row, size_t col, const size_t * rows, const size_t * cols, const Matrix & A)
{
    size_t newrows[dim - 1];
    size_t newcols[dim - 1];
    {
        size_t r;
        for (r = 0; r < row; ++r)
            newrows[r] = rows[r];
        for ( ; r < dim - 1; ++r)
            newrows[r] = rows[r + 1];
    }
    {
        size_t c;
        for (c = 0; c < col; ++c)
            newcols[c] = cols[c];
        for ( ; c < dim - 1; ++c)
            newcols[c] = cols[c + 1];
    }
    return det(dim - 1, newrows, newcols, A);
}


template <typename Matrix>
void adjugate(size_t dim, const Matrix & A, Matrix & Aadj)
{
    typedef typename Matrix::value_type Real;

    size_t rc[dim];
    for (size_t i = 0; i < dim; ++i)
        rc[i] = i;
    for (size_t j = 0; j < dim; ++j)
        for (size_t i = 0; i < dim; ++i) {
            // adjugate matrix = transpose of the cofactor matrix
            Real minor_ij = matminor(dim, i, j, rc, rc, A);
            if ((i + j) & 1)
                Aadj[j][i] = - minor_ij;
            else
                Aadj[j][i] = minor_ij;
        }
}

template <typename Matrix>
bool inverse(size_t dim, const Matrix & A, Matrix & Ainv)
{
    typedef typename Matrix::value_type Real;
    // A^-1 = 1 / det(A) * adj(A)
    Real detA = det(dim, A);
    if (detA == 0.)
        return false;
    Real scale = Real(1.) / det(dim, A);
    adjugate(dim, A, Ainv);
    for (size_t i = 0; i < dim; ++i)
        for (size_t j = 0; j < dim; ++j)
            Ainv[i][j] *= scale;
    return true;
}

// compute  b = Ax
template <typename Real, typename Matrix>
void multiply(size_t dim, const Matrix & A, const Real * x, Real * b, Real scale)
{
    for (size_t i = 0; i < dim; ++i) {
        Real r = 0;
        for (size_t j = 0; j < dim; ++j)
            r += A[i][j] * x[j];
        b[i] = scale * r;
    }
}

// solve  Ax = b  for x
template <typename Real, typename Matrix>
bool solve(size_t dim, const Matrix & A, Real * x, const Real * b)
{
    // Ax = b
    // A^-1 * A * x = A^-1 * b
    // x = 1 / det(A) * adj(A) * b
    Real detA = det(dim, A);
    if (detA == 0.)
        return false;
    Real scale = Real(1.) / detA;
    Matrix Aadj(dim);
    adjugate(dim, A, Aadj);
    multiply(dim, Aadj, b, x, scale);
    return true;
}

}  // namespace linear_algebra

}//namespace end

#endif
