#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "region_generator.hh"


/***** Math helper functions ******/

//! return square of argument
template <typename X>
inline X sqr( X x )
{ return x*x; }

//! Determinant of 3x3 matrix
inline double Determinant_3x3( const float *data )
{
    float detS = data[0]*(data[4]*data[8]-data[7]*data[5])
    - data[1]*(data[3]*data[8]-data[5]*data[6])
    + data[2]*(data[3]*data[7]-data[4]*data[6]);
    
    return detS;
}

//! Inverse of 3x3 matrix
inline void Inverse_3x3( const float *data, float *m )
{
    float invdet = 1.0f/Determinant_3x3( data );
    
    m[0] = (data[4]*data[8]-data[7]*data[5])*invdet;
    m[1] = -(data[1]*data[8]-data[2]*data[7])*invdet;
    m[2] = (data[1]*data[5]-data[2]*data[4])*invdet;
    m[3] = -(data[3]*data[8]-data[5]*data[6])*invdet;
    m[4] = (data[0]*data[8]-data[2]*data[6])*invdet;
    m[5] = -(data[0]*data[5]-data[2]*data[3])*invdet;
    m[6] = (data[3]*data[7]-data[4]*data[6])*invdet;
    m[7] = -(data[0]*data[7]-data[1]*data[6])*invdet;
    m[8] = (data[0]*data[4]-data[1]*data[3])*invdet;
}

void Inverse_4x4( float *mat )
{
    double tmp[12]; /* temp array for pairs */
    double src[16]; /* array of transpose source matrix */
    double det; /* determinant */
    double dst[16];
                
/* transpose matrix */
    for (int i = 0; i < 4; i++)
    {
        src[i] = mat[i*4];
        src[i + 4] = mat[i*4 + 1];
        src[i + 8] = mat[i*4 + 2];
        src[i + 12] = mat[i*4 + 3];
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
    dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
    dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
    dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
    dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
    dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
    dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
    dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
    dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
    dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
    dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
    dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
    dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
    dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
    dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
    dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
    dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
    
    /* calculate pairs for second 8 elements (cofactors) */
    tmp[0] = src[2]*src[7];
    tmp[1] = src[3]*src[6];
    tmp[2] = src[1]*src[7];
    tmp[3] = src[3]*src[5];
    tmp[4] = src[1]*src[6];
    tmp[5] = src[2]*src[5];
    tmp[6] = src[0]*src[7];
    tmp[7] = src[3]*src[4];
    tmp[8] = src[0]*src[6];
    tmp[9] = src[2]*src[4];
    tmp[10] = src[0]*src[5];
    tmp[11] = src[1]*src[4];
    
    /* calculate second 8 elements (cofactors) */
    dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
    dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
    dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
    dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
    dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
    dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
    dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
    dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
    dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
    dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
    dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
    dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
    dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
    dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
    dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
    dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
    
    /* calculate determinant */
    det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];
    
    /* calculate matrix inverse */
    det = 1/det;
    for (int j = 0; j < 16; j++)
    {    dst[j] *= det;
        mat[j] = dst[j];
    }
    
}

#include <xmmintrin.h>

//! Inversion of 4x4 matrix
//  Intel SIMD reference implementation
//  c.f. http://download.intel.com/design/PentiumIII/sml/24504301.pdf
inline void PIII_Inverse_4x4(float* src) {
    __m128 minor0, minor1, minor2, minor3;
    __m128 row0,   row1,   row2,   row3;
    __m128 det,    tmp1;
    tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src)), (__m64*)(src+ 4));
    row1 = _mm_loadh_pi(_mm_loadl_pi(row1, (__m64*)(src+8)), (__m64*)(src+12));
    row0 = _mm_shuffle_ps(tmp1, row1, 0x88);
    row1 = _mm_shuffle_ps(row1, tmp1, 0xDD);
    tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src+ 2)), (__m64*)(src+ 6));
    row3 = _mm_loadh_pi(_mm_loadl_pi(row3, (__m64*)(src+10)), (__m64*)(src+14));
    row2 = _mm_shuffle_ps(tmp1, row3, 0x88);
    row3 = _mm_shuffle_ps(row3, tmp1, 0xDD);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row2, row3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
    minor0 = _mm_mul_ps(row1, tmp1);
    minor1 = _mm_mul_ps(row0, tmp1);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
    minor0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), minor0);
    minor1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor1);
    minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row1, row2);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
    minor0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor0);
    minor3 = _mm_mul_ps(row0, tmp1);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
    minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row3, tmp1));
    minor3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor3);
    minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, 0x4E), row3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
    row2 = _mm_shuffle_ps(row2, row2, 0x4E);
    minor0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor0);
    minor2 = _mm_mul_ps(row0, tmp1);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
    minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row2, tmp1));
    minor2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor2);
    minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row0, row1);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
    minor2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor2);
    minor3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), minor3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
    minor2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), minor2);
    minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row2, tmp1));
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row0, row3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
    minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row2, tmp1));
    minor2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor2);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
    minor1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor1);
    minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row1, tmp1));
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row0, row2);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
    minor1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor1);
    minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row1, tmp1));
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
    minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row3, tmp1));
    minor3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor3);
    // -----------------------------------------------
    det = _mm_mul_ps(row0, minor0);
    det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
    det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
    tmp1 = _mm_rcp_ss(det);
    det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
    det = _mm_shuffle_ps(det, det, 0x00);
    minor0 = _mm_mul_ps(det, minor0);
    _mm_storel_pi((__m64*)(src), minor0);
    _mm_storeh_pi((__m64*)(src+2), minor0);
    minor1 = _mm_mul_ps(det, minor1);
    _mm_storel_pi((__m64*)(src+4), minor1);
    _mm_storeh_pi((__m64*)(src+6), minor1);
    minor2 = _mm_mul_ps(det, minor2);
    _mm_storel_pi((__m64*)(src+ 8), minor2);
    _mm_storeh_pi((__m64*)(src+10), minor2);
    minor3 = _mm_mul_ps(det, minor3);
    _mm_storel_pi((__m64*)(src+12), minor3);
    _mm_storeh_pi((__m64*)(src+14), minor3);
    
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
 */
class min_ellipsoid
{
protected:
    size_t N;
    float X[16];
    float c[3];
    float A[9], Ainv[9];
    float *Q;
    float *u;
    
    float v1[3],v2[3],v3[3],r1,r2,r3;
    
    float V[9], mu[3];
    
    bool axes_computed;
    
    void compute_axes( void )
    {
        gsl_vector     *eval;
        gsl_matrix     *evec;
        gsl_eigen_symmv_workspace *w;
        
        eval = gsl_vector_alloc(3);
        evec = gsl_matrix_alloc(3, 3);
        
        w = gsl_eigen_symmv_alloc(3);
        
        // promote to double, GSL wants double
        double dA[9];
        for( int i=0; i<9; ++i ) dA[i] = (double)A[i];
        
        gsl_matrix_view	m = gsl_matrix_view_array( dA, 3, 3);
        gsl_eigen_symmv( &m.matrix, eval, evec, w);
        
        gsl_eigen_symmv_sort( eval, evec, GSL_EIGEN_SORT_VAL_ASC);
        
        gsl_vector_view evec_i;
        
        for( int i=0; i<3; ++i )
        {
            mu[i] = gsl_vector_get(eval, i);
            evec_i = gsl_matrix_column (evec, i);
            for( int j=0; j<3; ++j )
                V[3*i+j] = gsl_vector_get(&evec_i.vector,j);
        }
        
        
        r1 = 1.0 / sqrt( gsl_vector_get(eval, 0) );
        r2 = 1.0 / sqrt( gsl_vector_get(eval, 1) );
        r3 = 1.0 / sqrt( gsl_vector_get(eval, 2) );
        
        evec_i = gsl_matrix_column (evec, 0);
        v1[0] = gsl_vector_get(&evec_i.vector,0);
        v1[1] = gsl_vector_get(&evec_i.vector,1);
        v1[2] = gsl_vector_get(&evec_i.vector,2);
        
        evec_i = gsl_matrix_column (evec, 1);
        v2[0] = gsl_vector_get(&evec_i.vector,0);
        v2[1] = gsl_vector_get(&evec_i.vector,1);
        v2[2] = gsl_vector_get(&evec_i.vector,2);
        
        evec_i = gsl_matrix_column (evec, 2);
        v3[0] = gsl_vector_get(&evec_i.vector,0);
        v3[1] = gsl_vector_get(&evec_i.vector,1);
        v3[2] = gsl_vector_get(&evec_i.vector,2);
        
        
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
        gsl_eigen_symmv_free (w);
        
        axes_computed = true;
    }
    
    void compute( double tol = 0.001, int maxit = 10000 )
    {
        double err = 10.0 * tol;
        float *unew = new float[N];
        int count = 0;
        double temp;
        
        while( err > tol && count < maxit )
        {
            for( int i=0; i<4; ++i )
                for( int j=0,i4=4*i; j<4; ++j )
                {
                    const int k = i4+j;
                    temp = 0.0;
                    for( size_t l=0; l<N; ++l )
                        temp += Q[4*l+i] * u[l] * Q[4*l+j];
                    X[k] = temp;
                }
            
            
            
            //PIII_Inverse_4x4(X);
            Inverse_4x4(X);
            
            
            int imax = 0; float Mmax = -1e30;
            double m;
            for( size_t i=0; i<N; ++i )
            {
                m = 0.0;
                for( int k=0; k<4; ++k )
                    for( int l=0; l<4; ++l )
                        m += Q[4*i+k] * X[4*l+k] * Q[4*i+l];
                if( m > Mmax )
                {
                    imax = i;
                    Mmax = m;
                }
            }
            
            float step_size = (Mmax-4.0f)/(4.0f*(Mmax-1.0f)), step_size1 = 1.0f-step_size;
            for( size_t i=0; i<N; ++i )
                unew[i] = u[i] * step_size1;
            unew[imax] += step_size;
            
            err = 0.0;
            for( size_t i=0; i<N; ++i )
            {
                err += sqr(unew[i]-u[i]);
                u[i] = unew[i];
            }
            err = sqrt(err);
            ++count;
        }
        
        if( count >= maxit )
            LOGERR("No convergence in min_ellipsoid::compute: maximum number of iterations reached!");
        
        delete[] unew;
    }
    
public:
    min_ellipsoid( size_t N_, float* P )
    : N( N_ ), axes_computed( false )
    {
        Q = new float[4*N];
        u = new float[N];
        
        for( size_t i=0; i<N; ++i )
            u[i] = 1.0/N;
        
        for( size_t i=0; i<N; ++i )
        {
            int i4=4*i, i3=3*i;
            for( size_t j=0; j<3; ++j )
                Q[i4+j] = P[i3+j];
            Q[i4+3] = 1.0f;
        }
        
        // compute the actual ellipsoid
        compute();
        
        // determine center
        for( int i=0; i<3; ++i )
        {
            c[i] = 0.0f;
            for( size_t j=0; j<N; ++j )
                c[i] += P[3*j+i] * u[j];
        }
        
        // determine A matrix
        float Pu[3];
        for( int j=0; j<3; ++j )
        {
            Pu[j] = 0.0;
            for( size_t i=0; i<N; ++i )
                Pu[j] += P[3*i+j] * u[i];
        }
        
        
        for( int i=0; i<3; ++i )
            for( int j=0,i3=3*i; j<3; ++j )
            {
                const int k = i3+j;
                Ainv[k] = 0.0f;
                for( size_t l=0; l<N; ++l )
                    Ainv[k] += P[3*l+i] * u[l] * P[3*l+j];
                Ainv[k] -= Pu[i]*Pu[j];
            }
        
        Inverse_3x3( Ainv, A );
        for( size_t i=0; i<9; ++i ){ A[i] *= 0.333333333; Ainv[i] *= 3.0; }
    }
    
    ~min_ellipsoid()
    {
        delete[] u;
        delete[] Q;
    }
    
    template<typename T>
    bool check_point( const T *x )
    {
        float q[3] = {x[0]-c[0],x[1]-c[1],x[2]-c[2]};
        
        double r = 0.0;
        for( int i=0; i<3; ++i )
            for( int j=0; j<3; ++j )
                r += q[i]*A[3*j+i]*q[j];
        
        return r <= 1.0;
    }
    
    void print( void )
    {
        std::cout << "A = \n";
        for( int i=0; i<9; ++i )
        {
            if( i%3==0 ) std::cout << std::endl;
            std::cout << A[i] << "   ";
        }
        std::cout << std::endl;
        std::cout << "Ainv = \n";
        for( int i=0; i<9; ++i )
        {
            if( i%3==0 ) std::cout << std::endl;
            std::cout << Ainv[i] << "   ";
        }
        std::cout << std::endl;
        std::cout << "c = (" << c[0] << ", " << c[1] << ", " << c[2] << ")\n";
    }
    
    template<typename T>
    void get_AABB( T *left, T *right )
    {
        for( int i=0; i<3; ++i )
        {
            left[i]  = c[i] - sqrt(Ainv[3*i+i]);
            right[i] = c[i] + sqrt(Ainv[3*i+i]);
        }
    }
    
    void get_center( float* xc )
    {
        for( int i=0; i<3; ++i ) xc[i] = c[i];
    }
    
    void expand_ellipsoid( float dr )
    {
        //print();
        
        if( !axes_computed )
            compute_axes();
        
        float munew[3];
        for( int i=0; i<3; ++i )
            munew[i] = 1.0/sqr(1.0/sqrt(mu[i])+dr);
        
        float Anew[9];
        for(int i=0; i<3; ++i )
            for( int j=0; j<3; ++j )
            {
                Anew[3*i+j] = 0.0;
                for( int k=0; k<3; ++k )
                    Anew[3*i+j] += V[3*k+i] * munew[k] * V[3*k+j];
            }
        
        for( int i=0; i<9; ++i )
            A[i] = Anew[i];
        
        Inverse_3x3( A, Ainv );
        
        //print();
    }
};




//! Minimum volume enclosing ellipsoid plugin
class region_ellipsoid_plugin : public region_generator_plugin{
private:
    
    min_ellipsoid *pellip_;
  int shift[3], shift_level;
    
    void read_points_from_file( std::string fname, std::vector<float>& p )
    {
        std::ifstream ifs(fname.c_str());
        if( !ifs )
        {
            LOGERR("region_ellipsoid_plugin::read_points_from_file : Could not open file \'%s\'",fname.c_str());
            throw std::runtime_error("region_ellipsoid_plugin::read_points_from_file : cannot open point file.");
        }
        while( ifs )
        {
            std::string s;
            if( !getline(ifs,s) )break;
            std::stringstream ss(s);
            while(ss)
            {
                if( !getline(ss,s,' ') ) break;
                p.push_back( strtod(s.c_str(),NULL) );
            }
        }
        
        if( p.size()%3 != 0 )
        {
            LOGERR("Region point file \'%s\' does not contain triplets",fname.c_str());
            throw std::runtime_error("region_ellipsoid_plugin::read_points_from_file : file does not contain triplets.");
        }
        
        double x0[3] = { p[0],p[1],p[2] }, dx;
        for( size_t i=3; i<p.size(); i+=3 )
        {
            for( size_t j=0; j<3; ++j )
            {
                dx = p[i+j]-x0[j];
                if( dx < -0.5 ) dx += 1.0;
                else if( dx > 0.5 ) dx -= 1.0;
                p[i+j] = x0[j] + dx;
            }
        }
    }
    
    void apply_shift( size_t Np, float *p, int *shift, int levelmin )
    {
        double dx = 1.0/(1<<levelmin);
        LOGINFO("unapplying previous shift to region particles :\n" \
                "\t [%d,%d,%d] = (%f,%f,%f)",shift[0],shift[1],shift[2],shift[0]*dx,shift[1]*dx,shift[2]*dx);
        
        for( size_t i=0,i3=0; i<Np; i++,i3+=3 )
            for( size_t j=0; j<3; ++j )
                p[i3+j] = p[i3+j]-shift[j]*dx;
    }
    
public:
    explicit region_ellipsoid_plugin( config_file& cf )
    : region_generator_plugin( cf )
    {
        std::vector<float> pp;
        
        
        std::string point_file = cf.getValue<std::string>("setup","region_point_file");
        read_points_from_file( point_file, pp );
        
        if( cf.containsKey("setup","region_point_shift") )
        {
            std::string point_shift = cf.getValue<std::string>("setup","region_point_shift");
            sscanf( point_shift.c_str(), "%d,%d,%d", &shift[0],&shift[1],&shift[2] );
            unsigned point_levelmin = cf.getValue<unsigned>("setup","region_point_levelmin");
        
            apply_shift( pp.size()/3, &pp[0], shift, point_levelmin );
            shift_level = point_levelmin;
        }
        
        
        pellip_ = new min_ellipsoid( pp.size()/3, &pp[0] );
        
        //expand the ellipsoid by one grid cell
        unsigned levelmax = cf.getValue<unsigned>("setup","levelmax");
        double dx = 1.0/(1<<levelmax);
        pellip_->expand_ellipsoid( dx );
        
        // output the center
        float c[3];
        pellip_->get_center( c );
        LOGINFO("Region center for ellipsoid determined at\n\t (%f,%f,%f)",c[0],c[1],c[2]);
    }
    
    ~region_ellipsoid_plugin()
    {
        delete pellip_;
    }
    
    void get_AABB( double *left, double *right, unsigned level )
    {
        pellip_->get_AABB( left, right );
        double dx = 1.0/(1ul<<level);
        
        for( int i=0;i<3;++i )
        {
            left[i] -= sqrt(3)*dx;
            right[i] += sqrt(3)*dx;
        }
        
    }
    
    bool query_point( double *x )
    {   return pellip_->check_point( x );   }
    
    bool is_grid_dim_forced( size_t* ndims )
    {   return false;   }
    
    void get_center( double *xc )
    {
        float c[3];
        pellip_->get_center( c );
        
        xc[0] = c[0];
        xc[1] = c[1];
        xc[2] = c[2];
    }

  void get_center_unshifted( double *xc )
  {
    double dx = 1.0/(1<<shift_level);
    float c[3];

    pellip_->get_center( c );        

    xc[0] = c[0]+shift[0]*dx;
    xc[1] = c[1]+shift[1]*dx;
    xc[2] = c[2]+shift[2]*dx;

  }

};

namespace{
    region_generator_plugin_creator_concrete< region_ellipsoid_plugin > creator("ellipsoid");
}
