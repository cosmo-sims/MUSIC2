#include <vector>
/*

 region_ellipsoid.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2010-13  Oliver Hahn

 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "region_generator.hh"

#include <array>
using vec3_t = std::array<double, 3>;

/***** Math helper functions ******/

//! return square of argument
template <typename X>
inline X sqr(X x)
{
    return x * x;
}

//! Determinant of 3x3 matrix
inline double Determinant_3x3(const double *data)
{
    double detS = data[0] * (data[4] * data[8] - data[7] * data[5]) - data[1] * (data[3] * data[8] - data[5] * data[6]) + data[2] * (data[3] * data[7] - data[4] * data[6]);

    return detS;
}

//! Inverse of 3x3 matrix
inline void Inverse_3x3(const double *data, double *m)
{
    double invdet = 1.0 / Determinant_3x3(data);

    m[0] = (data[4] * data[8] - data[7] * data[5]) * invdet;
    m[1] = -(data[1] * data[8] - data[2] * data[7]) * invdet;
    m[2] = (data[1] * data[5] - data[2] * data[4]) * invdet;
    m[3] = -(data[3] * data[8] - data[5] * data[6]) * invdet;
    m[4] = (data[0] * data[8] - data[2] * data[6]) * invdet;
    m[5] = -(data[0] * data[5] - data[2] * data[3]) * invdet;
    m[6] = (data[3] * data[7] - data[4] * data[6]) * invdet;
    m[7] = -(data[0] * data[7] - data[1] * data[6]) * invdet;
    m[8] = (data[0] * data[4] - data[1] * data[3]) * invdet;
}

void Inverse_4x4(double *mat)
{
    double tmp[12]; /* temp array for pairs */
    double src[16]; /* array of transpose source matrix */
    double det;     /* determinant */
    double dst[16];

    /* transpose matrix */
    for (int i = 0; i < 4; i++)
    {
        src[i] = mat[i * 4];
        src[i + 4] = mat[i * 4 + 1];
        src[i + 8] = mat[i * 4 + 2];
        src[i + 12] = mat[i * 4 + 3];
    }

    tmp[0] = src[10] * src[15];
    tmp[1] = src[11] * src[14];
    tmp[2] = src[9] * src[15];
    tmp[3] = src[11] * src[13];
    tmp[4] = src[9] * src[14];
    tmp[5] = src[10] * src[13];
    tmp[6] = src[8] * src[15];
    tmp[7] = src[11] * src[12];
    tmp[8] = src[8] * src[14];
    tmp[9] = src[10] * src[12];
    tmp[10] = src[8] * src[13];
    tmp[11] = src[9] * src[12];

    /* calculate first 8 elements (cofactors) */
    dst[0] = tmp[0] * src[5] + tmp[3] * src[6] + tmp[4] * src[7];
    dst[0] -= tmp[1] * src[5] + tmp[2] * src[6] + tmp[5] * src[7];
    dst[1] = tmp[1] * src[4] + tmp[6] * src[6] + tmp[9] * src[7];
    dst[1] -= tmp[0] * src[4] + tmp[7] * src[6] + tmp[8] * src[7];
    dst[2] = tmp[2] * src[4] + tmp[7] * src[5] + tmp[10] * src[7];
    dst[2] -= tmp[3] * src[4] + tmp[6] * src[5] + tmp[11] * src[7];
    dst[3] = tmp[5] * src[4] + tmp[8] * src[5] + tmp[11] * src[6];
    dst[3] -= tmp[4] * src[4] + tmp[9] * src[5] + tmp[10] * src[6];
    dst[4] = tmp[1] * src[1] + tmp[2] * src[2] + tmp[5] * src[3];
    dst[4] -= tmp[0] * src[1] + tmp[3] * src[2] + tmp[4] * src[3];
    dst[5] = tmp[0] * src[0] + tmp[7] * src[2] + tmp[8] * src[3];
    dst[5] -= tmp[1] * src[0] + tmp[6] * src[2] + tmp[9] * src[3];
    dst[6] = tmp[3] * src[0] + tmp[6] * src[1] + tmp[11] * src[3];
    dst[6] -= tmp[2] * src[0] + tmp[7] * src[1] + tmp[10] * src[3];
    dst[7] = tmp[4] * src[0] + tmp[9] * src[1] + tmp[10] * src[2];
    dst[7] -= tmp[5] * src[0] + tmp[8] * src[1] + tmp[11] * src[2];

    /* calculate pairs for second 8 elements (cofactors) */
    tmp[0] = src[2] * src[7];
    tmp[1] = src[3] * src[6];
    tmp[2] = src[1] * src[7];
    tmp[3] = src[3] * src[5];
    tmp[4] = src[1] * src[6];
    tmp[5] = src[2] * src[5];
    tmp[6] = src[0] * src[7];
    tmp[7] = src[3] * src[4];
    tmp[8] = src[0] * src[6];
    tmp[9] = src[2] * src[4];
    tmp[10] = src[0] * src[5];
    tmp[11] = src[1] * src[4];

    /* calculate second 8 elements (cofactors) */
    dst[8] = tmp[0] * src[13] + tmp[3] * src[14] + tmp[4] * src[15];
    dst[8] -= tmp[1] * src[13] + tmp[2] * src[14] + tmp[5] * src[15];
    dst[9] = tmp[1] * src[12] + tmp[6] * src[14] + tmp[9] * src[15];
    dst[9] -= tmp[0] * src[12] + tmp[7] * src[14] + tmp[8] * src[15];
    dst[10] = tmp[2] * src[12] + tmp[7] * src[13] + tmp[10] * src[15];
    dst[10] -= tmp[3] * src[12] + tmp[6] * src[13] + tmp[11] * src[15];
    dst[11] = tmp[5] * src[12] + tmp[8] * src[13] + tmp[11] * src[14];
    dst[11] -= tmp[4] * src[12] + tmp[9] * src[13] + tmp[10] * src[14];
    dst[12] = tmp[2] * src[10] + tmp[5] * src[11] + tmp[1] * src[9];
    dst[12] -= tmp[4] * src[11] + tmp[0] * src[9] + tmp[3] * src[10];
    dst[13] = tmp[8] * src[11] + tmp[0] * src[8] + tmp[7] * src[10];
    dst[13] -= tmp[6] * src[10] + tmp[9] * src[11] + tmp[1] * src[8];
    dst[14] = tmp[6] * src[9] + tmp[11] * src[11] + tmp[3] * src[8];
    dst[14] -= tmp[10] * src[11] + tmp[2] * src[8] + tmp[7] * src[9];
    dst[15] = tmp[10] * src[10] + tmp[4] * src[8] + tmp[9] * src[9];
    dst[15] -= tmp[8] * src[9] + tmp[11] * src[10] + tmp[5] * src[8];

    /* calculate determinant */
    det = src[0] * dst[0] + src[1] * dst[1] + src[2] * dst[2] + src[3] * dst[3];

    /* calculate matrix inverse */
    det = 1 / det;
    for (int j = 0; j < 16; j++)
    {
        dst[j] *= det;
        mat[j] = dst[j];
    }
}

/***** Minimum Volume Bounding Ellipsoid Implementation ******/
/*
 * Finds the minimum volume enclosing ellipsoid (MVEE) of a set of data
 * points stored in matrix P. The following optimization problem is solved:
 *
 *     minimize log(det(A))
 *     s.t. (P_i - c)'*A*(P_i - c)<= 1
 *
 * in variables A and c, where P_i is the i-th column of the matrix P.
 * The solver is based on Khachiyan Algorithm, and the final solution is
 * different from the optimal value by the pre-specified amount of 'tolerance'.
 *
 * The ellipsoid equation is given in the canonical form
 *     (x-c)' A (x-c) <= 1
 *
 * Code was adapted from the MATLAB version by Nima Moshtagh (nima@seas.upenn.edu)
 */
class min_ellipsoid
{
protected:
    size_t N;
    double X[16];
    double c[3];
    double A[9], Ainv[9];
    double *Q;
    double *u;

    double detA, detA13;

    double v1[3], v2[3], v3[3], r1, r2, r3;

    double V[9], mu[3];

    bool axes_computed, hold_point_data;

    void compute_axes(void)
    {
        gsl_vector *eval;
        gsl_matrix *evec;
        gsl_eigen_symmv_workspace *w;

        eval = gsl_vector_alloc(3);
        evec = gsl_matrix_alloc(3, 3);

        w = gsl_eigen_symmv_alloc(3);

        // promote to double, GSL wants double
        double dA[9];
        for (int i = 0; i < 9; ++i)
            dA[i] = (double)A[i];

        gsl_matrix_view m = gsl_matrix_view_array(dA, 3, 3);
        gsl_eigen_symmv(&m.matrix, eval, evec, w);

        gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

        gsl_vector_view evec_i;

        for (int i = 0; i < 3; ++i)
        {
            mu[i] = gsl_vector_get(eval, i);
            evec_i = gsl_matrix_column(evec, i);
            for (int j = 0; j < 3; ++j)
                V[3 * i + j] = gsl_vector_get(&evec_i.vector, j);
        }

        r1 = 1.0 / sqrt(gsl_vector_get(eval, 0));
        r2 = 1.0 / sqrt(gsl_vector_get(eval, 1));
        r3 = 1.0 / sqrt(gsl_vector_get(eval, 2));

        evec_i = gsl_matrix_column(evec, 0);
        v1[0] = gsl_vector_get(&evec_i.vector, 0);
        v1[1] = gsl_vector_get(&evec_i.vector, 1);
        v1[2] = gsl_vector_get(&evec_i.vector, 2);

        evec_i = gsl_matrix_column(evec, 1);
        v2[0] = gsl_vector_get(&evec_i.vector, 0);
        v2[1] = gsl_vector_get(&evec_i.vector, 1);
        v2[2] = gsl_vector_get(&evec_i.vector, 2);

        evec_i = gsl_matrix_column(evec, 2);
        v3[0] = gsl_vector_get(&evec_i.vector, 0);
        v3[1] = gsl_vector_get(&evec_i.vector, 1);
        v3[2] = gsl_vector_get(&evec_i.vector, 2);

        gsl_vector_free(eval);
        gsl_matrix_free(evec);
        gsl_eigen_symmv_free(w);

        axes_computed = true;
    }

    // use the Khachiyan Algorithm to find the minimum bounding ellipsoid
    void compute(double tol = 1e-4, int maxit = 10000)
    {
        double err = 10.0 * tol;
        double *unew = new double[N];
        int count = 0;
        double temp;

        while (err > tol && count < maxit)
        {
            for (int i = 0; i < 4; ++i)
                for (int j = 0, i4 = 4 * i; j < 4; ++j)
                {
                    const int k = i4 + j;
                    temp = 0.0;
                    for (size_t l = 0; l < N; ++l)
                        temp += (double)(Q[4 * l + i] * u[l] * Q[4 * l + j]);
                    X[k] = temp;
                }

            Inverse_4x4(X);

            int imax = 0;
            double Mmax = -1e30;
            double m;
            for (size_t i = 0; i < N; ++i)
            {
                m = 0.0;
                for (int k = 0; k < 4; ++k)
                    for (int l = 0; l < 4; ++l)
                        m += (double)(Q[4 * i + k] * X[4 * l + k] * Q[4 * i + l]);
                if (m > Mmax)
                {
                    imax = i;
                    Mmax = m;
                }
            }

            double step_size = (Mmax - 4.0) / (4.0 * (Mmax - 1.0)), step_size1 = 1.0 - step_size;
            for (size_t i = 0; i < N; ++i)
                unew[i] = u[i] * step_size1;
            unew[imax] += step_size;

            err = 0.0;
            for (size_t i = 0; i < N; ++i)
            {
                err += sqr(unew[i] - u[i]);
                u[i] = unew[i];
            }
            err = sqrt(err);
            ++count;
        }

        if (count >= maxit)
            music::elog.Print("No convergence in min_ellipsoid::compute: maximum number of iterations reached!");

        delete[] unew;
    }

public:
    min_ellipsoid(size_t N_, double *P)
        : N(N_), axes_computed(false), hold_point_data(true)
    {
        // --- initialize ---
        music::ilog.Print("computing minimum bounding ellipsoid from %lld points", N);

        Q = new double[4 * N];
        u = new double[N];

        // normalize coordinate frame
        double xcenter[3] = {0.0, 0.0, 0.0};
        double pmax[3] = {-1e30,-1e30,-1e30}, pmin[3] = {1e30,1e30,1e30};
        for (size_t i = 0; i < N; ++i)
        {
            int i3 = 3 * i;
            for (size_t j = 0; j < 3; ++j)
            {
                double d = P[i3 + j] - P[j];
                double p = (d > 0.5) ? P[i3 + j] - 1.0 : (d < -0.5) ? P[i3 + j] + 1.0
                                                                    : P[i3 + j];

                xcenter[j] += p;
                pmax[j] = p > pmax[j] ? p : pmax[j];
                pmin[j] = p < pmin[j] ? p : pmin[j];
            }
        }
        double L = 2.0 * std::min((pmax[0] - pmin[0]),std::min(pmax[1]-pmin[1],pmax[2]-pmin[2]));
        
        xcenter[0] = fmod(xcenter[0] / N + 1.0, 1.0);
        xcenter[1] = fmod(xcenter[1] / N + 1.0, 1.0);
        xcenter[2] = fmod(xcenter[2] / N + 1.0, 1.0);

        for (size_t i = 0; i < N; ++i)
        {
            int i3 = 3 * i;
            for (size_t j = 0; j < 3; ++j)
            {
                double d = P[i3 + j] - xcenter[j];
                d = (d > 0.5) ? d - 1.0 : (d < -0.5) ? d + 1.0
                                                     : d;

                P[i3 + j] = d / L + 0.5;
            }
        }

        for (size_t i = 0; i < N; ++i)
            u[i] = 1.0 / N;

        for (size_t i = 0; i < N; ++i)
        {
            int i4 = 4 * i, i3 = 3 * i;
            for (size_t j = 0; j < 3; ++j)
                Q[i4 + j] = P[i3 + j];
            Q[i4 + 3] = 1.0;
        }

        //--- compute the actual ellipsoid using the Khachiyan Algorithm ---
        compute();

        //--- determine the ellipsoid A matrix ---
        double Pu[3];
        for (int j = 0; j < 3; ++j)
        {
            Pu[j] = 0.0;
            for (size_t i = 0; i < N; ++i)
                Pu[j] += P[3 * i + j] * u[i];
        }

        // determine center
        c[0] = Pu[0];
        c[1] = Pu[1];
        c[2] = Pu[2];

        // need to do summation in double precision due to
        // possible catastrophic cancellation issues when
        // using many input points
        double Atmp[9];
        for (int i = 0; i < 3; ++i)
            for (int j = 0, i3 = 3 * i; j < 3; ++j)
            {
                const int k = i3 + j;
                Atmp[k] = 0.0;
                for (size_t l = 0; l < N; ++l)
                    Atmp[k] += P[3 * l + i] * u[l] * P[3 * l + j];
                Atmp[k] -= Pu[i] * Pu[j];
            }

        for (int i = 0; i < 9; ++i)
            Ainv[i] = Atmp[i] * L * L;

        c[0] = L * (c[0] - 0.5) + xcenter[0];
        c[1] = L * (c[1] - 0.5) + xcenter[1];
        c[2] = L * (c[2] - 0.5) + xcenter[2];

        Inverse_3x3(Ainv, A);
        for (size_t i = 0; i < 9; ++i)
        {
            A[i] /= 3.0;
            Ainv[i] *= 3.0;
        }

        detA = Determinant_3x3(A);
        detA13 = pow(detA, 1.0 / 3.0);
    }

    min_ellipsoid(const double *A_, const double *c_)
        : N(0), axes_computed(false), hold_point_data(false)
    {
        for (int i = 0; i < 9; ++i)
        {
            A[i] = A_[i];
            Ainv[i] = 0.0;
        }
        for (int i = 0; i < 3; ++i)
            c[i] = c_[i];
    }

    min_ellipsoid(const min_ellipsoid &e)
        : N(0), hold_point_data(false)
    {
        for (int i = 0; i < 16; ++i)
            X[i] = e.X[i];

        for (int i = 0; i < 3; ++i)
        {
            c[i] = e.c[i];
            v1[i] = e.v1[i];
            v2[i] = e.v2[i];
            v3[i] = e.v3[i];
            mu[i] = e.mu[i];
        }

        for (int i = 0; i < 9; ++i)
        {
            A[i] = e.A[i];
            Ainv[i] = e.Ainv[i];
            V[i] = e.V[i];
        }

        N = e.N;
        detA = e.detA;
        detA13 = e.detA13;
        axes_computed = e.axes_computed;
    }

    ~min_ellipsoid()
    {
        if (hold_point_data)
        {
            delete[] u;
            delete[] Q;
        }
    }

    bool check_point(const vec3_t &x, double dist = 0.0) const
    {
        dist = (dist + 1.0) * detA13;

        vec3_t q{{x[0] - c[0], x[1] - c[1], x[2] - c[2]}};

        double r = 0.0;
        for (int i = 0; i < 3; ++i)
            q[i] = (q[i] > 0.5) ? q[i] - 1.0 : (q[i] < -0.5) ? q[i] + 1.0
                                                             : q[i];

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                r += q[i] * A[3 * j + i] * q[j];

        return r <= 1.0;
    }

    void print(void)
    {
        std::cout << "A = \n";
        for (int i = 0; i < 9; ++i)
        {
            if (i % 3 == 0)
                std::cout << std::endl;
            std::cout << A[i] << "   ";
        }
        std::cout << std::endl;
        std::cout << "Ainv = \n";
        for (int i = 0; i < 9; ++i)
        {
            if (i % 3 == 0)
                std::cout << std::endl;
            std::cout << Ainv[i] << "   ";
        }
        std::cout << std::endl;
        std::cout << "c = (" << c[0] << ", " << c[1] << ", " << c[2] << ")\n";
    }

    void get_AABB(vec3_t &left, vec3_t &right)
    {
        for (int i = 0; i < 3; ++i)
        {
            left[i] = c[i] - sqrt(Ainv[3 * i + i]);
            right[i] = c[i] + sqrt(Ainv[3 * i + i]);

            if (left[i] < 0.0)
            {
                left[i] += 1.0;
                right[i] += 1.0;
            }
            if (left[i] >= 1.0)
            {
                left[i] -= 1.0;
                right[i] -= 1.0;
            }
        }
    }

    void get_center(vec3_t &xc) const
    {
        for (int i = 0; i < 3; ++i)
            xc[i] = c[i];
    }

    void get_matrix(double *AA)
    {
        for (int i = 0; i < 9; ++i)
            AA[i] = A[i];
    }

    double sgn(double x)
    {
        if (x < 0.0)
            return -1.0;
        return 1.0;
    }

    void expand_ellipsoid(double dr)
    {
        if (!axes_computed)
        {
            music::ulog.Print("computing ellipsoid axes.....");
            compute_axes();
        }
        // double muold[3] = {mu[0],mu[1],mu[2]};
        double munew[3];
        for (int i = 0; i < 3; ++i)
            munew[i] = sgn(mu[i]) / sqr(1.0 / sqrt(fabs(mu[i])) + dr);

        double Anew[9];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                Anew[3 * i + j] = 0.0;
                for (int k = 0; k < 3; ++k)
                    Anew[3 * i + j] += V[3 * k + i] * munew[k] * V[3 * k + j];
            }

        for (int i = 0; i < 9; ++i)
            A[i] = Anew[i];

        Inverse_3x3(A, Ainv);

        // music::ulog.Print("computing ellipsoid axes.....");
        compute_axes();

        // print();
    }
};

#include "point_file_reader.hh"

#if 0
#include "convex_hull.hh"
#endif

//! Minimum volume enclosing ellipsoid plugin
class region_ellipsoid_plugin : public region_generator_plugin
{
private:
    std::vector<min_ellipsoid *> pellip_;
    int shift[3], shift_level, padding_;
    double vfac_;
    bool do_extra_padding_;

    void apply_shift(size_t Np, double *p, int *shift, int levelmin)
    {
        double dx = 1.0 / (1 << levelmin);
        music::ilog.Print("unapplying shift of previous zoom region to region particles :\n"
                "\t [%d,%d,%d] = (%f,%f,%f)",
                shift[0], shift[1], shift[2], shift[0] * dx, shift[1] * dx, shift[2] * dx);

        for (size_t i = 0, i3 = 0; i < Np; i++, i3 += 3)
            for (size_t j = 0; j < 3; ++j)
                p[i3 + j] = p[i3 + j] - shift[j] * dx;
    }

public:
    explicit region_ellipsoid_plugin(config_file &cf)
        : region_generator_plugin(cf)
    {
        std::vector<double> pp;

        for (unsigned i = 0; i <= levelmax_; ++i)
            pellip_.push_back(NULL);

        // sanity check
        if (!cf.contains_key("setup", "region_point_file") &&
            !(cf.contains_key("setup", "region_ellipsoid_matrix[0]") &&
              cf.contains_key("setup", "region_ellipsoid_matrix[1]") &&
              cf.contains_key("setup", "region_ellipsoid_matrix[2]") &&
              cf.contains_key("setup", "region_ellipsoid_center")))
        {
            music::elog.Print("Insufficient data for region=ellipsoid\n Need to specify either \'region_point_file\' or the ellipsoid equation.");
            throw std::runtime_error("Insufficient data for region=ellipsoid");
        }
        //

        vfac_ = cf.get_value<double>("cosmology", "vfact");
        padding_ = cf.get_value<int>("setup", "padding");

        std::string point_file;

        if (cf.contains_key("setup", "region_point_file"))
        {
            point_file = cf.get_value<std::string>("setup", "region_point_file");

            point_reader pfr;
            pfr.read_points_from_file(point_file, vfac_, pp);

            // if file has more than three columns, just take first three
            // at the moment...
            if (pfr.num_columns > 3)
            {
                std::vector<double> xx;
                xx.reserve(3 * pp.size() / pfr.num_columns);

                for (size_t i = 0; i < pp.size() / pfr.num_columns; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        xx.push_back(pp[pfr.num_columns * i + j]);

                pp.swap(xx);
            }

            if (cf.contains_key("setup", "region_point_shift"))
            {
                std::string point_shift = cf.get_value<std::string>("setup", "region_point_shift");
                if (sscanf(point_shift.c_str(), "%d,%d,%d", &shift[0], &shift[1], &shift[2]) != 3)
                {
                    music::elog.Print("Error parsing triple for region_point_shift");
                    throw std::runtime_error("Error parsing triple for region_point_shift");
                }
                unsigned point_levelmin = cf.get_value<unsigned>("setup", "region_point_levelmin");

                apply_shift(pp.size() / 3, &pp[0], shift, point_levelmin);
                shift_level = point_levelmin;
            }

            pellip_[levelmax_] = new min_ellipsoid(pp.size() / 3, &pp[0]);
        }
        else
        {
            double A[9] = {0}, c[3] = {0};
            std::string strtmp;

            strtmp = cf.get_value<std::string>("setup", "region_ellipsoid_matrix[0]");
            if (sscanf(strtmp.c_str(), "%lf,%lf,%lf", &A[0], &A[1], &A[2]) != 3)
            {
                music::elog.Print("Error parsing triple for region_ellipsoid_matrix[0]");
                throw std::runtime_error("Error parsing triple for region_ellipsoid_matrix[0]");
            }
            strtmp = cf.get_value<std::string>("setup", "region_ellipsoid_matrix[1]");
            if (sscanf(strtmp.c_str(), "%lf,%lf,%lf", &A[3], &A[4], &A[5]) != 3)
            {
                music::elog.Print("Error parsing triple for region_ellipsoid_matrix[1]");
                throw std::runtime_error("Error parsing triple for region_ellipsoid_matrix[1]");
            }
            strtmp = cf.get_value<std::string>("setup", "region_ellipsoid_matrix[2]");
            if (sscanf(strtmp.c_str(), "%lf,%lf,%lf", &A[6], &A[7], &A[8]) != 3)
            {
                music::elog.Print("Error parsing triple for region_ellipsoid_matrix[2]");
                throw std::runtime_error("Error parsing triple for region_ellipsoid_matrix[2]");
            }
            strtmp = cf.get_value<std::string>("setup", "region_ellipsoid_center");
            if (sscanf(strtmp.c_str(), "%lf,%lf,%lf", &c[0], &c[1], &c[2]) != 3)
            {
                music::elog.Print("Error parsing triple for region_ellipsoid_center");
                throw std::runtime_error("Error parsing triple for region_ellipsoid_center");
            }

            pellip_[levelmax_] = new min_ellipsoid(A, c);
        }

#if 0
        if( false )
        {
            // compute convex hull and use only hull points to speed things up
            // BUT THIS NEEDS MORE TESTING BEFORE IT GOES IN THE REPO
            music::ilog.Print("Computing convex hull for %llu points", pp.size()/3 );
            convex_hull<double> ch( &pp[0], pp.size()/3 );
            std::set<int> unique;
            ch.get_defining_indices( unique );
            std::set<int>::iterator it = unique.begin();
            
            std::vector<double> pphull;
            pphull.reserve( unique.size()*3 );
            
            while( it != unique.end() )
            {
                int idx = *it;
                
                pphull.push_back( pp[3*idx+0] );
                pphull.push_back( pp[3*idx+1] );
                pphull.push_back( pp[3*idx+2] );
                
                ++it;
            }
            
            pphull.swap( pp );
        }
#endif

        // output the center
        vec3_t c;
        double A[9];
        pellip_[levelmax_]->get_center(c);
        pellip_[levelmax_]->get_matrix(A);

        music::ilog << "Region center for ellipsoid determined at" << std::endl
            <<"\t xc = ( " << c[0] << " " << c[1] << " " << c[2] << " )" << std::endl;
        music::ilog << "Ellipsoid matrix determined as" <<std ::endl 
            << "\t      / "  << std::setw(8) << A[0] << " " << std::setw(8) << A[1] << " " << std::setw(8) << A[2] << " \\" << std::endl
            << "\t  A = | "  << std::setw(8) << A[3] << " " << std::setw(8) << A[4] << " " << std::setw(8) << A[5] << " |" << std::endl
            << "\t      \\ " << std::setw(8) << A[6] << " " << std::setw(8) << A[7] << " " << std::setw(8) << A[8] << " /" << std::endl;

        // expand the ellipsoid by one grid cell

        unsigned levelmax = cf.get_value<unsigned>("setup", "levelmax");
        unsigned npad = cf.get_value<unsigned>("setup", "padding");
        double dx = 1.0 / (1ul << levelmax);
        pellip_[levelmax_]->expand_ellipsoid(dx);

        // generate the higher level ellipsoids
        for (int ilevel = levelmax_; ilevel > 0; --ilevel)
        {
            dx = 1.0 / (1ul << (ilevel - 1));
            pellip_[ilevel - 1] = new min_ellipsoid(*pellip_[ilevel]);
            pellip_[ilevel - 1]->expand_ellipsoid(npad * dx);
        }

        //-----------------------------------------------------------------
        // when querying the bounding box, do we need extra padding?
        do_extra_padding_ = false;

        // conditions should be added here
        {
            std::string output_plugin = cf.get_value<std::string>("output", "format");
            if (output_plugin == std::string("grafic2"))
                do_extra_padding_ = true;
        }
    }

    ~region_ellipsoid_plugin()
    {
        for (unsigned i = 0; i <= levelmax_; ++i)
            if (pellip_[i] != NULL)
                delete pellip_[i];
    }

    void get_AABB(vec3_t &left, vec3_t &right, unsigned level) const
    {
        if (level <= levelmin_)
        {
            left[0] = left[1] = left[2] = 0.0;
            right[0] = right[1] = right[2] = 1.0;
            return;
        }

        pellip_[level]->get_AABB(left, right);

        double dx = 1.0 / (1ul << level);
        double pad = (double)(padding_ + 1) * dx;

        if (!do_extra_padding_)
            pad = 0.0;

        double ext = sqrt(3) * dx + pad;

        for (int i = 0; i < 3; ++i)
        {
            left[i] -= ext;
            right[i] += ext;
        }
    }

    void update_AABB(vec3_t &left, vec3_t &right)
    {
        // we ignore this, the grid generator must have generated a grid that contains the ellipsoid
        // it might have enlarged it, but who cares...
    }

    bool query_point(const vec3_t &x, int level) const
    {
        return pellip_[level]->check_point(x);
    }

    bool is_grid_dim_forced(index3_t &ndims) const
    {
        return false;
    }

    void get_center(vec3_t &xc) const
    {
        pellip_[levelmax_]->get_center(xc);
    }

    void get_center_unshifted(vec3_t &xc) const
    {
        double dx = 1.0 / (1 << shift_level);
        pellip_[levelmax_]->get_center(xc);

        xc[0] += shift[0] * dx;
        xc[1] += shift[1] * dx;
        xc[2] += shift[2] * dx;
    }
};

namespace
{
    region_generator_plugin_creator_concrete<region_ellipsoid_plugin> creator("ellipsoid");
}
