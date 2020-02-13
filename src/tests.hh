/*
 
 tests.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
*/


#ifndef __TESTS_HH
#define __TESTS_HH

#include <math.h>


inline double CIC_interp_back( const MeshvarBnd<double>& A, double x, double y, double z )
{
	int 
		ix  = (int)x,
		iy  = (int)y,
		iz  = (int)z,
		ix1 = (ix+1),
		iy1 = (iy+1),
		iz1 = (iz+1);
	
	
    double
		dx = (double)(x - (double)ix),
		dy = (double)(y - (double)iy),
		dz = (double)(z - (double)iz),
		tx = 1.0-dx,
		ty = 1.0-dy,
		tz = 1.0-dz;
	
    double
		f_xyz = A(ix,iy,iz)*tx*ty*tz,
		f_Xyz = A(ix1,iy,iz)*dx*ty*tz,
		f_xYz = A(ix,iy1,iz)*tx*dy*tz,
		f_xyZ = A(ix,iy,iz1)*tx*ty*dz,
		f_XYz = A(ix1,iy1,iz)*dx*dy*tz,
		f_XyZ = A(ix1,iy,iz1)*dx*ty*dz,
		f_xYZ = A(ix,iy1,iz1)*tx*dy*dz,
		f_XYZ = A(ix1,iy1,iz1)*dx*dy*dz;
	
    return f_xyz + f_Xyz + f_xYz + f_xyZ + f_XYz + f_XyZ + f_xYZ + f_XYZ;
}

inline double TSC_interp_back( const MeshvarBnd<double>& A, double x, double y, double z )
{
	double val = 0.0;
    int xngp = (int)x, yngp = (int)y, zngp = (int)z;
        
	for( int xx = xngp-1; xx <= xngp+1; ++xx )
	{
		double weightx = 1.0;
		double dx = fabs(x-(double)xx);
		int axx(xx);
		
		if( xx==xngp )
			weightx *= 0.75-dx*dx;
		else{
			weightx *= 1.125 - 1.5*dx + 0.5*dx*dx;
		}
		
		for( int yy = yngp-1; yy <= yngp+1; ++yy )
		{
			double weighty = weightx;
			double dy = fabs(y-(double)yy);
			int ayy(yy);
			
			if( yy==yngp )
				weighty *= 0.75-dy*dy;
			else{
				weighty *= 1.125 - 1.5*dy + 0.5*dy*dy;
			}
			
			for( int zz = zngp-1; zz <= zngp+1; ++zz )
			{
				double weightz = weighty;
				double dz = fabs(z-(double)zz);
				int azz(zz);
				
				if( zz==zngp )
					weightz *= 0.75-dz*dz;
				else{
					weightz *= 1.125 - 1.5*dz + 0.5*dz*dz;
				}
				
				val += A(axx,ayy,azz) * weightz;
			}
		}
	}
	
	return val;
}

class TestProblem{
public:
	MeshvarBnd<double> m_rho, m_uana, m_ubnd, m_xgrad, m_ygrad, m_zgrad;
	int m_nb, m_nres;
	double m_h;
	
	TestProblem( int nb, int nres )
	: m_rho( nb, nres ), m_uana( nb, nres ), m_ubnd( nb, nres ),
	m_xgrad( nb, nres ), m_ygrad( nb, nres ), m_zgrad( nb, nres ),
	m_nb( nb ), m_nres( nres ), m_h( 1.0/((double)nres ) )//m_h( 1.0/((double)nres+1.0 ) )
	{ }
	
};

class TSC_Test : public TestProblem{
public:
	double m_q;
	
	class TSCcube{
	public:
		std::vector<double> m_data;
		
		
		TSCcube()
		{
			m_data.assign(27,0.0);
			
			//.. center
			(*this)( 0, 0, 0) = 27./64.;
			
			//.. faces
			(*this)(-1, 0, 0) = 
			(*this)(+1, 0, 0) = 
			(*this)( 0,-1, 0) = 
			(*this)( 0,+1, 0) = 
			(*this)( 0, 0,-1) = 
			(*this)( 0, 0,+1) = 9./128.;
			
			//.. edges
			(*this)(-1,-1, 0) =
			(*this)(-1,+1, 0) =
			(*this)(+1,-1, 0) =
			(*this)(+1,+1, 0) =
			(*this)(-1, 0,-1) =
			(*this)(-1, 0,+1) =
			(*this)(+1, 0,-1) =
			(*this)(+1, 0,+1) =
			(*this)( 0,-1,-1) =
			(*this)( 0,-1,+1) =
			(*this)( 0,+1,-1) =
			(*this)( 0,+1,+1) = 3./256.;
			
			//.. corners
			(*this)(-1,-1,-1) =
			(*this)(-1,+1,-1) =
			(*this)(-1,-1,+1) =
			(*this)(-1,+1,+1) =
			(*this)(+1,-1,-1) =
			(*this)(+1,+1,-1) =
			(*this)(+1,-1,+1) =
			(*this)(+1,+1,+1) = 1./512.;
			
		}
		
		double& operator()(int i, int j, int k)
		{ return m_data[ ((i+1)*3+(j+1))*3 +(k+1)]; }
		
		const double& operator()(int i, int j, int k) const
		{ return m_data[ ((i+1)*3+(j+1))*3 +(k+1)]; }
	};
	
	TSC_Test( int nb, int nres, double q=-1.0 )
	: TestProblem(nb, nres), m_q(q)
	{
		TSCcube c;
		int xm(nres/2-1), ym(nres/2-1), zm(nres/2-1);
		double xxm((double)xm*m_h), yym((double)ym*m_h), zzm((double)zm*m_h);
		
		double fourpi = 4.0*M_PI;
		
		m_uana.zero();
		m_ubnd.zero();
		m_xgrad.zero();
		m_ygrad.zero();
		m_zgrad.zero();
		
		for( int i=-nb; i<nres+nb; ++i )
			for( int j=-nb; j<nres+nb; ++j )
				for( int k=-nb; k<nres+nb; ++k )
				{
					
					//double xxm((double)xm), yym((double)ym), zzm((double)zm);
					double xx((double)i*m_h), yy((double)j*m_h), zz((double)k*m_h);
					
					
					for( int ix=-1; ix<=1; ++ix )
						for( int iy=-1; iy<=1; ++iy )
							for( int iz=-1; iz<=1; ++iz )
							{
								double dx(xx-(xxm+ix*m_h)), dy(yy-(yym+iy*m_h)), dz(zz-(zzm+iz*m_h));
								double d3 = pow(dx*dx+dy*dy+dz*dz,1.5);
								
								double dphi = m_q*c(ix,iy,iz)/sqrt(dx*dx+dy*dy+dz*dz);
								
								if( i==xm && j==ym && k==zm )
									m_rho(i+ix,j+iy,k+iz) = m_q*c(ix,iy,iz)/(m_h*m_h*m_h);
								
								if( d3 < 1e-10 )
									continue;
								
								m_uana(i,j,k) += dphi/fourpi;
								m_ubnd(i,j,k) += dphi/fourpi;
								
								m_xgrad(i,j,k) -= m_q*c(ix,iy,iz)*dx/d3/fourpi;
								m_ygrad(i,j,k) -= m_q*c(ix,iy,iz)*dy/d3/fourpi;
								m_zgrad(i,j,k) -= m_q*c(ix,iy,iz)*dz/d3/fourpi;
								
								
							}
				}
		
		
		//m_rho(xm,ym,zm) = 4.0*M_PI*m_q/(m_h*m_h*m_h);
		
	}
};


class PointMassTest : public TestProblem{
public:
	double m_q;
	
	PointMassTest( int nb, int nres, double q=-1.0 )
	: TestProblem(nb, nres), m_q( q )
	{
		//int xm(nres/2-1), ym(nres/2-1), zm(nres/2-1);
		int xm(nres/2), ym(nres/2), zm(nres/2);
		double xxm((double)xm*m_h), yym((double)ym*m_h), zzm((double)zm*m_h);
		m_rho.zero();
		
		double fourpi = 4.0*M_PI;
		
		for( int i=-nb; i<nres+nb; ++i )
			for( int j=-nb; j<nres+nb; ++j )
				for( int k=-nb; k<nres+nb; ++k )
				{
					
					//double xxm((double)xm), yym((double)ym), zzm((double)zm);
					double xx((double)i*m_h), yy((double)j*m_h), zz((double)k*m_h);
					
					double dx(xx-xxm), dy(yy-yym), dz(zz-zzm);
					m_uana(i,j,k) = m_q/sqrt(dx*dx+dy*dy+dz*dz)/fourpi;///2.0/M_PI;
					m_ubnd(i,j,k) = 0.0;//m_uana(i,j,k);
					
					double d3 = pow(dx*dx+dy*dy+dz*dz,1.5);
					m_xgrad(i,j,k) = -m_q*dx/d3/fourpi;
					m_ygrad(i,j,k) = -m_q*dy/d3/fourpi;
					m_zgrad(i,j,k) = -m_q*dz/d3/fourpi;
				}
		
		
		for( int iy=0; iy<nres; ++iy )
			for( int iz=0; iz<nres; ++iz )
			{
				double dx=0.5+0.5*m_h, dy=((double)iy+0.5)*m_h-0.5, dz=((double)iz+0.5)*m_h-0.5, d=sqrt(dx*dx+dy*dy+dz*dz);
				m_ubnd(-1,iy,iz) = m_q/d/fourpi;
				dx = 0.5-0.5*m_h;
				d = sqrt(dx*dx+dy*dy+dz*dz);
				m_ubnd(nres,iy,iz) = m_q/d/fourpi;
			}
		
		for( int ix=0; ix<nres; ++ix )
			for( int iz=0; iz<nres; ++iz )
			{
				double dx=((double)ix+0.5)*m_h-0.5, dy=0.5+0.5*m_h, dz=((double)iz+0.5)*m_h-0.5, d=sqrt(dx*dx+dy*dy+dz*dz);
				m_ubnd(ix,-1,iz) = m_q/d/fourpi;
				dy=0.5-0.5*m_h;
				d = sqrt(dx*dx+dy*dy+dz*dz);
				m_ubnd(ix,nres,iz) = m_q/d/fourpi;
			}
		
		for( int ix=0; ix<nres; ++ix )
			for( int iy=0; iy<nres; ++iy )
			{
				double dx=((double)ix+0.5)*m_h-0.5, dy=((double)iy+0.5)*m_h-0.5, dz=0.5+0.5*m_h, d=sqrt(dx*dx+dy*dy+dz*dz);
				m_ubnd(ix,iy,-1) = m_q/d/fourpi;
				dz=0.5-0.5*m_h;
				d = sqrt(dx*dx+dy*dy+dz*dz);
				m_ubnd(ix,iy,nres) = m_q/d/fourpi;
			}
		
		
		m_rho(xm,ym,zm) = m_q/(m_h*m_h*m_h);
	}	
};


#endif