/*
 *  solver.h
 *  GravitySolver
 *
 *  Created by Oliver Hahn on 1/20/10.
 *  Copyright 2010 KIPAC/Stanford University. All rights reserved.
 *
 */

#ifndef __SOLVER_HH
#define __SOLVER_HH

#include <cmath>
#include <iostream>
#include "mesh.hh"

#define BEGIN_MULTIGRID_NAMESPACE namespace multigrid {
#define END_MULTIGRID_NAMESPACE }

BEGIN_MULTIGRID_NAMESPACE
	
namespace opt {
	enum smtype { sm_jacobi, sm_gauss_seidel, sm_sor };
}

template< class S, class O, typename T=double >
class solver
{
public:
	typedef S scheme;
	typedef O mgop;

protected:
	scheme				m_scheme;
	mgop				m_gridop;
	unsigned			m_npresmooth, m_npostsmooth;
	opt::smtype			m_smoother;
	unsigned			m_ilevelmin;
	
	const static bool	m_bperiodic = true;

	GridHierarchy<T>	*m_pu, *m_pf, *m_pfsave;	
	GridHierarchy<bool> *m_pmask;
	const MeshvarBnd<T> *m_pubnd;
	
	double compute_error( const MeshvarBnd<T>& u, const MeshvarBnd<T>& unew );
	
	double compute_error( const GridHierarchy<T>& uh, const GridHierarchy<T>& uhnew, bool verbose );

protected:
	
	void Jacobi( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f );
	
	void GaussSeidel( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f );
	
	void SOR( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f );
	
	void twoGrid( unsigned ilevel );
	
	void interp_coarse_fine( unsigned ilevel, MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine, bool bcf=true );
	
	void setBC( unsigned ilevel );
	
	void make_periodic( MeshvarBnd<T> *u );
	
	void interp_cubic( MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine, int itop, int jtop, int ktop, int i, int j, int k );
	void interp_coarse_fine_cubic( unsigned ilevel, MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine, bool bcf );
	
public:
	solver( GridHierarchy<T>& f, //const MeshvarBnd<T>& uBC_top, 
				   opt::smtype smoother, unsigned npresmooth, unsigned npostsmooth );
	
	~solver()
	{ delete m_pmask; }
	
	double solve( GridHierarchy<T>& u, double accuracy, double h=-1.0, bool verbose=false );
	
	double solve( GridHierarchy<T>& u, double accuracy, bool verbose=false )
	{
		return this->solve ( u, accuracy, -1.0, verbose );
	}
	
	
	
};


template< class S, class O, typename T >
solver<S,O,T>::solver( GridHierarchy<T>& f, //const MeshvarBnd<T>& ubnd, 
					opt::smtype smoother, unsigned npresmooth, unsigned npostsmooth )
:	m_scheme(), m_gridop(), m_npresmooth( npresmooth ), m_npostsmooth( npostsmooth ), 
m_smoother( smoother ), m_ilevelmin( f.levelmin() ), m_pf( &f )//, m_pubnd( &ubnd )
{ 
	//... initialize the refinement mask
	m_pmask = new GridHierarchy<bool>( f.m_nbnd );
	m_pmask->create_base_hierarchy(f.levelmin());
			
	for( unsigned ilevel=f.levelmin()+1; ilevel<=f.levelmax(); ++ilevel )
	{
		meshvar_bnd* pf = f.get_grid(ilevel);
		m_pmask->add_patch( pf->offset(0), pf->offset(1), pf->offset(2), pf->size(0), pf->size(1), pf->size(2) );
	}
	
	m_pmask->zero();
	
	for( unsigned ilevel=0; ilevel<f.levelmin(); ++ilevel )
	{
		MeshvarBnd<T> *pf = f.get_grid(ilevel);
		for( int ix=0; ix < (int)pf->size(0); ++ix )
			for( int iy=0; iy < (int)pf->size(1); ++iy )
				for( int iz=0; iz < (int)pf->size(2); ++iz )
					(*m_pmask->get_grid(ilevel))(ix,iy,iz) = true;
	}
	
	for( unsigned ilevel=m_ilevelmin; ilevel<f.levelmax(); ++ilevel )
	{
		MeshvarBnd<T>* pf = f.get_grid(ilevel+1);//, *pfc = f.get_grid(ilevel);
		
		for( int ix=pf->offset(0); ix < (int)(pf->offset(0)+pf->size(0)/2); ++ix )
			for( int iy=pf->offset(1); iy < (int)(pf->offset(1)+pf->size(1)/2); ++iy )
				for( int iz=pf->offset(2); iz < (int)(pf->offset(2)+pf->size(2)/2); ++iz )
					(*m_pmask->get_grid(ilevel))(ix,iy,iz) = true;
	}
		
}


template< class S, class O, typename T >
void solver<S,O,T>::Jacobi( T h, MeshvarBnd<T> *u, const MeshvarBnd<T>* f )
{
	int
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
	
	double 
		c0 = -1.0/m_scheme.ccoeff(),
		h2 = h*h; 
	
	MeshvarBnd<T> uold(*u);
	
	double alpha = 0.95, ialpha = 1.0-alpha;
	
	#pragma omp parallel for
	for( int ix=0; ix<nx; ++ix )
		for( int iy=0; iy<ny; ++iy )
			for( int iz=0; iz<nz; ++iz )
				(*u)(ix,iy,iz) = ialpha * uold(ix,iy,iz) + alpha * (m_scheme.rhs( uold, ix, iy, iz ) + h2 * (*f)(ix,iy,iz))*c0;
	
	
	
	
}

template< class S, class O, typename T >
void solver<S,O,T>::SOR( T h, MeshvarBnd<T> *u, const MeshvarBnd<T>* f )
{
	int
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);

	double 
		c0 = -1.0/m_scheme.ccoeff(),
		h2 = h*h; 
		
	MeshvarBnd<T> uold(*u);
	
	double 
		alpha = 1.2, 
	//alpha = 2 / (1 + 4 * atan(1.0) / double(u->size(0)))-1.0,
		ialpha = 1.0-alpha;
	
	//std::cerr << "omega_opt = " << alpha << std::endl;
	
	#pragma omp parallel for
	for( int ix=0; ix<nx; ++ix )
		for( int iy=0; iy<ny; ++iy )
			for( int iz=0; iz<nz; ++iz )
				if( (ix+iy+iz)%2==0 )
					(*u)(ix,iy,iz) = ialpha * uold(ix,iy,iz) + alpha * (m_scheme.rhs( uold, ix, iy, iz ) + h2 * (*f)(ix,iy,iz))*c0;
	
	
	#pragma omp parallel for
	for( int ix=0; ix<nx; ++ix )
		for( int iy=0; iy<ny; ++iy )
			for( int iz=0; iz<nz; ++iz )
				if( (ix+iy+iz)%2!=0 )
					(*u)(ix,iy,iz) = ialpha * uold(ix,iy,iz) + alpha * (m_scheme.rhs( *u, ix, iy, iz ) + h2 * (*f)(ix,iy,iz))*c0;
	
	
	
}

template< class S, class O, typename T >
void solver<S,O,T>::GaussSeidel( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f )
{
	int 
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
	
	T
		c0 = -1.0/m_scheme.ccoeff(),
		h2 = h*h; 
	
	for( int color=0; color < 2; ++color )
		#pragma omp parallel for
		for( int ix=0; ix<nx; ++ix )
			for( int iy=0; iy<ny; ++iy )
				for( int iz=0; iz<nz; ++iz )
					if( (ix+iy+iz)%2 == color )
						(*u)(ix,iy,iz) = (m_scheme.rhs( *u, ix, iy, iz ) + h2 * (*f)(ix,iy,iz))*c0;
	
}


template< class S, class O, typename T >
void solver<S,O,T>::twoGrid( unsigned ilevel )
{
	MeshvarBnd<T> *uf, *uc, *ff, *fc;
	
	T 
		h = 1.0/(pow(2.0,ilevel)),
		c0 = -1.0/m_scheme.ccoeff(),
		h2 = h*h; 
	
	uf = m_pu->get_grid(ilevel);
	ff = m_pf->get_grid(ilevel);	
	
	uc = m_pu->get_grid(ilevel-1);
	fc = m_pf->get_grid(ilevel-1);	
	
	int 
		nx = uf->size(0), 
		ny = uf->size(1), 
		nz = uf->size(2);
	
	if( m_bperiodic && ilevel <= m_ilevelmin)
		make_periodic( uf );
	else if(!m_bperiodic)
		setBC( ilevel );
	
	//... do smoothing sweeps with specified solver
	for( unsigned i=0; i<m_npresmooth; ++i ){
		
		if( ilevel > m_ilevelmin )
			interp_coarse_fine(ilevel, *uc, *uf );
		
		if( m_smoother == opt::sm_gauss_seidel )
			GaussSeidel( h, uf, ff );
			
		else if( m_smoother == opt::sm_jacobi )
			Jacobi( h, uf, ff);		
			
		else if( m_smoother == opt::sm_sor )
			SOR( h, uf, ff );
		
		if( m_bperiodic && ilevel <= m_ilevelmin )
			make_periodic( uf );
	}
			
	
	m_gridop.restrict( *uf, *uc );
	
	//... essential!!
	if( m_bperiodic && ilevel <= m_ilevelmin )
		make_periodic( uc );
	else if( m_bperiodic )
		interp_coarse_fine(ilevel,*uc,*uf);
	
	meshvar_bnd Lu(*uf,false);
	Lu.zero();
	#pragma omp parallel for
	for( int ix=0; ix<nx; ++ix )
		for( int iy=0; iy<ny; ++iy )
			for( int iz=0; iz<nz; ++iz )
				Lu(ix,iy,iz) = m_scheme.apply( (*uf), ix, iy, iz )/h2;
	
	meshvar_bnd tLu(*uc,false);
	
	//... restrict Lu
	m_gridop.restrict( Lu, tLu );
	Lu.deallocate();
	
	//... restrict source term
	m_gridop.restrict( *ff, *fc );
	
	//... compute RHS tau-correction
	#pragma omp parallel for schedule(dynamic)
	for( int ix=0; ix<(int)uc->size(0); ++ix )
		for( int iy=0; iy<(int)uc->size(1); ++iy )
			for( int iz=0; iz<(int)uc->size(2); ++iz )
				if( (*m_pmask->get_grid(ilevel-1))(ix,iy,iz) == true )
					(*fc)(ix,iy,iz) += ((tLu( ix, iy, iz ) - (m_scheme.apply( *uc, ix, iy, iz )/(4.0*h2))));
				
					
	tLu.deallocate();
	
	meshvar_bnd ucsave(*uc,true);
						
	//... have we reached the end of the recursion or do we need to go up one level?
	if( ilevel == 1 )
		if( m_bperiodic )
			(*uc)(0,0,0) = 0.0;
		else 
			(*uc)(0,0,0) = (m_scheme.rhs( (*uc), 0, 0, 0 ) + 4.0 * h2 * (*fc)(0,0,0))*c0;
	else
		twoGrid( ilevel-1 );
	
	meshvar_bnd cc(*uc,false);
	
	//... compute correction on coarse grid
	#pragma omp parallel for
	for( int ix=0; ix<(int)cc.size(0); ++ix )
		for( int iy=0; iy<(int)cc.size(1); ++iy )
			for( int iz=0; iz<(int)cc.size(2); ++iz )
				cc(ix,iy,iz) = (*uc)(ix,iy,iz) - ucsave(ix,iy,iz);
		
	ucsave.deallocate();


	//... prolongate correction to fine grid	
	meshvar_bnd cf(*uf,false);
	m_gridop.prolong( cc, cf );
	
	cc.deallocate();
	
	
	#pragma omp parallel for 
	for( int ix=0; ix<nx; ++ix )
		for( int iy=0; iy<ny; ++iy )
			for( int iz=0; iz<nz; ++iz )
				(*uf)(ix,iy,iz) += cf(ix,iy,iz);
	

	cf.deallocate();
				
	//... interpolate and apply coarse-fine boundary conditions on fine level
	if( m_bperiodic && ilevel <= m_ilevelmin )
		make_periodic( uf );
	else if(!m_bperiodic)
		setBC( ilevel );
	
	//if( ilevel > m_ilevelmin )
	//	interp_coarse_fine(ilevel, *uc, *uf );

	//... do smoothing sweeps with specified solver
	for( unsigned i=0; i<m_npostsmooth; ++i ){
		
		if( ilevel > m_ilevelmin )
			interp_coarse_fine(ilevel, *uc, *uf );

		if( m_smoother == opt::sm_gauss_seidel )
			GaussSeidel( h, uf, ff );
		
		else if( m_smoother == opt::sm_jacobi )
			Jacobi( h, uf, ff);		
		
		else if( m_smoother == opt::sm_sor )
			SOR( h, uf, ff );
		
		if( m_bperiodic && ilevel <= m_ilevelmin )
			make_periodic( uf );

	}
}

template< class S, class O, typename T >
double solver<S,O,T>::compute_error( const MeshvarBnd<T>& u, const MeshvarBnd<T>& unew )
{
	int 
		nx = u.size(0), 
		ny = u.size(1), 
		nz = u.size(2);
	
	double err = 0.0;
	unsigned count = 0;
	
#pragma omp parallel for reduction(+:err,count)
	for( int ix=0; ix<nx; ++ix )
		for( int iy=0; iy<ny; ++iy )
			for( int iz=0; iz<nz; ++iz )
				if( fabs(unew(ix,iy,iz)) > 0.0 )//&& u(ix,iy,iz) != unew(ix,iy,iz) )
				{
					err += fabs(1.0 - u(ix,iy,iz)/unew(ix,iy,iz));
					++count;
				}
	
	if( count != 0 )
		err /= count;
	
	return err;
}

template< class S, class O, typename T >
double solver<S,O,T>::compute_error( const GridHierarchy<T>& uh, const GridHierarchy<T>& uhnew, bool verbose )
{
	double maxerr = 0.0;
	
	for( unsigned ilevel=uh.levelmin(); ilevel <= uh.levelmax(); ++ilevel )
	{
		double err = 0.0;
		err = compute_error( *uh.get_grid(ilevel), *uhnew.get_grid(ilevel) );
		
		if( verbose )
			std::cout << "    Level " << std::setw(6) << ilevel << ",   Error = " << err << std::endl;
		maxerr = std::max(maxerr,err);
		
	}
	return maxerr;
}

template< class S, class O, typename T >
double solver<S,O,T>::solve( GridHierarchy<T>& uh, double acc, double h, bool verbose )
{

	double err;
	
	GridHierarchy<T> uhnew(uh);//, fsave(*m_pf);
	m_pu = &uh;
	
    unsigned niter = 0;
	
	//... iterate ...//
	while (true)
	{
		
		
		twoGrid( uh.levelmax() );
		err = compute_error( *m_pu, uhnew, verbose );
		++niter;
		
		if( verbose ){
			std::cout << "--> Step No. " << std::setw(3) << niter << ", Max Err = " << err << std::endl;
			std::cout << "-------------------------------------------------------------\n";
		}
			
		if( (niter > 1) && ((err < acc) || (niter > 20)) )
			break;
		
		uhnew = *m_pu;
		//*m_pf = fsave;
	}		
	
	if( err > acc )
		std::cout << "Error : no convergence in Poisson solver" << std::endl;
	else if( verbose )
		std::cout << " - Converged in " << niter << " steps to req. acc. of " << acc << std::endl;

	
	//uh = uhnew;
	//*m_pf = fsave;
	return err;
}

inline double interp2( double x1, double x2, double x3, double f1, double f2, double f3, double x )
{
	double a,b,c;	
	a = (x1 * f3 - x3 * f1 - x2 * f3 - x1 * f2 + x2 * f1 + x3 * f2) / (x1 * x3 * x3 - x2 * x3 * x3 + x2 * x1 * x1 - x3 * x1 * x1 + x3 * x2 * x2 - x1 * x2 * x2);
	b = -(x1 * x1 * f3 - x1 * x1 * f2 - f1 * x3 * x3 + f2 * x3 * x3 - x2 * x2 * f3 + f1 * x2 * x2) / (x1 - x2) / (x1 * x2 - x1 * x3 + x3 * x3 - x2 * x3);
	c = (x1 * x1 * x2 * f3 - x1 * x1 * x3 * f2 - x2 * x2 * x1 * f3 + f2 * x1 * x3 * x3 + x2 * x2 * x3 * f1 - f1 * x2 * x3 * x3) / (x1 - x2) / (x1 * x2 - x1 * x3 + x3 * x3 - x2 * x3);
	
	return a*x*x+b*x+c;
}

inline double interp2( double fleft, double fcenter, double fright, double x )
{
	double a,b,c;
	a = 0.5*(fleft+fright)-fcenter;
	b = 0.5*(fright-fleft);
	c = fcenter;
	
	return a*x*x+b*x+c;
}


inline double interp2left( double fleft, double fcenter, double fright )
{
	double a,b,c;
	a = (6.0*fright-10.0*fcenter+4.0*fleft)/15.0;
	b = (-4.0*fleft+9.0*fright-5.0*fcenter)/15.0;
	c = fcenter;
	
	return a-b+c;
}

inline double interp2right( double fleft, double fcenter, double fright )
{
	double a,b,c;
	a = (6.0*fleft-10.0*fcenter+4.0*fright)/15.0;
	b = (4.0*fright-9.0*fleft+5.0*fcenter)/15.0;
	c = fcenter;
	
	return a+b+c;
}

template< class S, class O, typename T >
void solver<S,O,T>::interp_cubic( MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine, int i, int j, int k, int itop, int jtop, int ktop )
{
	MeshvarBnd<T> &u    = fine;
	MeshvarBnd<T> &utop = coarse;
	
	/*
	u(i+0,j+0,k+0) = ( -125.*utop(itop-2,jtop-2,ktop-2) +875.*utop(itop-2,jtop-2,ktop-1) +2625.*utop(itop-2,jtop-2,ktop) 
					  -175.*utop(itop-2,jtop-2,ktop+1) +875.*utop(itop-2,jtop-1,ktop-2) -6125.*utop(itop-2,jtop-1,ktop-1) 
					  -18375.*utop(itop-2,jtop-1,ktop) +1225.*utop(itop-2,jtop-1,ktop+1) +2625.*utop(itop-2,jtop,ktop-2) 
					  -18375.*utop(itop-2,jtop,ktop-1) -55125.*utop(itop-2,jtop,ktop) +3675.*utop(itop-2,jtop,ktop+1) 
					  -175.*utop(itop-2,jtop+1,ktop-2) +1225.*utop(itop-2,jtop+1,ktop-1) +3675.*utop(itop-2,jtop+1,ktop) 
					  -245.*utop(itop-2,jtop+1,ktop+1) +875.*utop(itop-1,jtop-2,ktop-2) -6125.*utop(itop-1,jtop-2,ktop-1) 
					  -18375.*utop(itop-1,jtop-2,ktop) +1225.*utop(itop-1,jtop-2,ktop+1) -6125.*utop(itop-1,jtop-1,ktop-2) 
					  +42875.*utop(itop-1,jtop-1,ktop-1) +128625.*utop(itop-1,jtop-1,ktop) -8575.*utop(itop-1,jtop-1,ktop+1) 
					  -18375.*utop(itop-1,jtop,ktop-2) +128625.*utop(itop-1,jtop,ktop-1) +385875.*utop(itop-1,jtop,ktop) 
					  -25725.*utop(itop-1,jtop,ktop+1) +1225.*utop(itop-1,jtop+1,ktop-2) -8575.*utop(itop-1,jtop+1,ktop-1) 
					  -25725.*utop(itop-1,jtop+1,ktop) +1715.*utop(itop-1,jtop+1,ktop+1) +2625.*utop(itop,jtop-2,ktop-2) 
					  -18375.*utop(itop,jtop-2,ktop-1) -55125.*utop(itop,jtop-2,ktop) +3675.*utop(itop,jtop-2,ktop+1) 
					  -18375.*utop(itop,jtop-1,ktop-2) +128625.*utop(itop,jtop-1,ktop-1) +385875.*utop(itop,jtop-1,ktop) 
					  -25725.*utop(itop,jtop-1,ktop+1) -55125.*utop(itop,jtop,ktop-2) +385875.*utop(itop,jtop,ktop-1) 
					  +1157625.*utop(itop,jtop,ktop) -77175.*utop(itop,jtop,ktop+1) +3675.*utop(itop,jtop+1,ktop-2) 
					  -25725.*utop(itop,jtop+1,ktop-1) -77175.*utop(itop,jtop+1,ktop) +5145.*utop(itop,jtop+1,ktop+1) 
					  -175.*utop(itop+1,jtop-2,ktop-2) +1225.*utop(itop+1,jtop-2,ktop-1) +3675.*utop(itop+1,jtop-2,ktop) 
					  -245.*utop(itop+1,jtop-2,ktop+1) +1225.*utop(itop+1,jtop-1,ktop-2) -8575.*utop(itop+1,jtop-1,ktop-1) 
					  -25725.*utop(itop+1,jtop-1,ktop) +1715.*utop(itop+1,jtop-1,ktop+1) +3675.*utop(itop+1,jtop,ktop-2) 
					  -25725.*utop(itop+1,jtop,ktop-1) -77175.*utop(itop+1,jtop,ktop) +5145.*utop(itop+1,jtop,ktop+1) 
					  -245.*utop(itop+1,jtop+1,ktop-2) +1715.*utop(itop+1,jtop+1,ktop-1) +5145.*utop(itop+1,jtop+1,ktop) 
					  -343.*utop(itop+1,jtop+1,ktop+1) )/2097152.;
	u(i+0,j+0,k+1) = ( -175.*utop(itop-2,jtop-2,ktop-1) +2625.*utop(itop-2,jtop-2,ktop) +875.*utop(itop-2,jtop-2,ktop+1) 
					  -125.*utop(itop-2,jtop-2,ktop+2) +1225.*utop(itop-2,jtop-1,ktop-1) -18375.*utop(itop-2,jtop-1,ktop) 
					  -6125.*utop(itop-2,jtop-1,ktop+1) +875.*utop(itop-2,jtop-1,ktop+2) +3675.*utop(itop-2,jtop,ktop-1) 
					  -55125.*utop(itop-2,jtop,ktop) -18375.*utop(itop-2,jtop,ktop+1) +2625.*utop(itop-2,jtop,ktop+2) 
					  -245.*utop(itop-2,jtop+1,ktop-1) +3675.*utop(itop-2,jtop+1,ktop) +1225.*utop(itop-2,jtop+1,ktop+1) 
					  -175.*utop(itop-2,jtop+1,ktop+2) +1225.*utop(itop-1,jtop-2,ktop-1) -18375.*utop(itop-1,jtop-2,ktop) 
					  -6125.*utop(itop-1,jtop-2,ktop+1) +875.*utop(itop-1,jtop-2,ktop+2) -8575.*utop(itop-1,jtop-1,ktop-1) 
					  +128625.*utop(itop-1,jtop-1,ktop) +42875.*utop(itop-1,jtop-1,ktop+1) -6125.*utop(itop-1,jtop-1,ktop+2) 
					  -25725.*utop(itop-1,jtop,ktop-1) +385875.*utop(itop-1,jtop,ktop) +128625.*utop(itop-1,jtop,ktop+1) 
					  -18375.*utop(itop-1,jtop,ktop+2) +1715.*utop(itop-1,jtop+1,ktop-1) -25725.*utop(itop-1,jtop+1,ktop) 
					  -8575.*utop(itop-1,jtop+1,ktop+1) +1225.*utop(itop-1,jtop+1,ktop+2) +3675.*utop(itop,jtop-2,ktop-1) 
					  -55125.*utop(itop,jtop-2,ktop) -18375.*utop(itop,jtop-2,ktop+1) +2625.*utop(itop,jtop-2,ktop+2) 
					  -25725.*utop(itop,jtop-1,ktop-1) +385875.*utop(itop,jtop-1,ktop) +128625.*utop(itop,jtop-1,ktop+1) 
					  -18375.*utop(itop,jtop-1,ktop+2) -77175.*utop(itop,jtop,ktop-1) +1157625.*utop(itop,jtop,ktop) 
					  +385875.*utop(itop,jtop,ktop+1) -55125.*utop(itop,jtop,ktop+2) +5145.*utop(itop,jtop+1,ktop-1) 
					  -77175.*utop(itop,jtop+1,ktop) -25725.*utop(itop,jtop+1,ktop+1) +3675.*utop(itop,jtop+1,ktop+2) 
					  -245.*utop(itop+1,jtop-2,ktop-1) +3675.*utop(itop+1,jtop-2,ktop) +1225.*utop(itop+1,jtop-2,ktop+1) 
					  -175.*utop(itop+1,jtop-2,ktop+2) +1715.*utop(itop+1,jtop-1,ktop-1) -25725.*utop(itop+1,jtop-1,ktop) 
					  -8575.*utop(itop+1,jtop-1,ktop+1) +1225.*utop(itop+1,jtop-1,ktop+2) +5145.*utop(itop+1,jtop,ktop-1) 
					  -77175.*utop(itop+1,jtop,ktop) -25725.*utop(itop+1,jtop,ktop+1) +3675.*utop(itop+1,jtop,ktop+2) 
					  -343.*utop(itop+1,jtop+1,ktop-1) +5145.*utop(itop+1,jtop+1,ktop) +1715.*utop(itop+1,jtop+1,ktop+1) 
					  -245.*utop(itop+1,jtop+1,ktop+2) )/2097152.;
	u(i+0,j+1,k+0) = ( -175.*utop(itop-2,jtop-1,ktop-2) +1225.*utop(itop-2,jtop-1,ktop-1) +3675.*utop(itop-2,jtop-1,ktop) 
					  -245.*utop(itop-2,jtop-1,ktop+1) +2625.*utop(itop-2,jtop,ktop-2) -18375.*utop(itop-2,jtop,ktop-1) 
					  -55125.*utop(itop-2,jtop,ktop) +3675.*utop(itop-2,jtop,ktop+1) +875.*utop(itop-2,jtop+1,ktop-2) 
					  -6125.*utop(itop-2,jtop+1,ktop-1) -18375.*utop(itop-2,jtop+1,ktop) +1225.*utop(itop-2,jtop+1,ktop+1) 
					  -125.*utop(itop-2,jtop+2,ktop-2) +875.*utop(itop-2,jtop+2,ktop-1) +2625.*utop(itop-2,jtop+2,ktop) 
					  -175.*utop(itop-2,jtop+2,ktop+1) +1225.*utop(itop-1,jtop-1,ktop-2) -8575.*utop(itop-1,jtop-1,ktop-1) 
					  -25725.*utop(itop-1,jtop-1,ktop) +1715.*utop(itop-1,jtop-1,ktop+1) -18375.*utop(itop-1,jtop,ktop-2) 
					  +128625.*utop(itop-1,jtop,ktop-1) +385875.*utop(itop-1,jtop,ktop) -25725.*utop(itop-1,jtop,ktop+1) 
					  -6125.*utop(itop-1,jtop+1,ktop-2) +42875.*utop(itop-1,jtop+1,ktop-1) +128625.*utop(itop-1,jtop+1,ktop) 
					  -8575.*utop(itop-1,jtop+1,ktop+1) +875.*utop(itop-1,jtop+2,ktop-2) -6125.*utop(itop-1,jtop+2,ktop-1) 
					  -18375.*utop(itop-1,jtop+2,ktop) +1225.*utop(itop-1,jtop+2,ktop+1) +3675.*utop(itop,jtop-1,ktop-2) 
					  -25725.*utop(itop,jtop-1,ktop-1) -77175.*utop(itop,jtop-1,ktop) +5145.*utop(itop,jtop-1,ktop+1) 
					  -55125.*utop(itop,jtop,ktop-2) +385875.*utop(itop,jtop,ktop-1) +1157625.*utop(itop,jtop,ktop) 
					  -77175.*utop(itop,jtop,ktop+1) -18375.*utop(itop,jtop+1,ktop-2) +128625.*utop(itop,jtop+1,ktop-1) 
					  +385875.*utop(itop,jtop+1,ktop) -25725.*utop(itop,jtop+1,ktop+1) +2625.*utop(itop,jtop+2,ktop-2) 
					  -18375.*utop(itop,jtop+2,ktop-1) -55125.*utop(itop,jtop+2,ktop) +3675.*utop(itop,jtop+2,ktop+1) 
					  -245.*utop(itop+1,jtop-1,ktop-2) +1715.*utop(itop+1,jtop-1,ktop-1) +5145.*utop(itop+1,jtop-1,ktop) 
					  -343.*utop(itop+1,jtop-1,ktop+1) +3675.*utop(itop+1,jtop,ktop-2) -25725.*utop(itop+1,jtop,ktop-1) 
					  -77175.*utop(itop+1,jtop,ktop) +5145.*utop(itop+1,jtop,ktop+1) +1225.*utop(itop+1,jtop+1,ktop-2) 
					  -8575.*utop(itop+1,jtop+1,ktop-1) -25725.*utop(itop+1,jtop+1,ktop) +1715.*utop(itop+1,jtop+1,ktop+1) 
					  -175.*utop(itop+1,jtop+2,ktop-2) +1225.*utop(itop+1,jtop+2,ktop-1) +3675.*utop(itop+1,jtop+2,ktop) 
					  -245.*utop(itop+1,jtop+2,ktop+1) )/2097152.;
	u(i+0,j+1,k+1) = ( -245.*utop(itop-2,jtop-1,ktop-1) +3675.*utop(itop-2,jtop-1,ktop) +1225.*utop(itop-2,jtop-1,ktop+1) 
					  -175.*utop(itop-2,jtop-1,ktop+2) +3675.*utop(itop-2,jtop,ktop-1) -55125.*utop(itop-2,jtop,ktop) 
					  -18375.*utop(itop-2,jtop,ktop+1) +2625.*utop(itop-2,jtop,ktop+2) +1225.*utop(itop-2,jtop+1,ktop-1) 
					  -18375.*utop(itop-2,jtop+1,ktop) -6125.*utop(itop-2,jtop+1,ktop+1) +875.*utop(itop-2,jtop+1,ktop+2) 
					  -175.*utop(itop-2,jtop+2,ktop-1) +2625.*utop(itop-2,jtop+2,ktop) +875.*utop(itop-2,jtop+2,ktop+1) 
					  -125.*utop(itop-2,jtop+2,ktop+2) +1715.*utop(itop-1,jtop-1,ktop-1) -25725.*utop(itop-1,jtop-1,ktop) 
					  -8575.*utop(itop-1,jtop-1,ktop+1) +1225.*utop(itop-1,jtop-1,ktop+2) -25725.*utop(itop-1,jtop,ktop-1) 
					  +385875.*utop(itop-1,jtop,ktop) +128625.*utop(itop-1,jtop,ktop+1) -18375.*utop(itop-1,jtop,ktop+2) 
					  -8575.*utop(itop-1,jtop+1,ktop-1) +128625.*utop(itop-1,jtop+1,ktop) +42875.*utop(itop-1,jtop+1,ktop+1) 
					  -6125.*utop(itop-1,jtop+1,ktop+2) +1225.*utop(itop-1,jtop+2,ktop-1) -18375.*utop(itop-1,jtop+2,ktop) 
					  -6125.*utop(itop-1,jtop+2,ktop+1) +875.*utop(itop-1,jtop+2,ktop+2) +5145.*utop(itop,jtop-1,ktop-1) 
					  -77175.*utop(itop,jtop-1,ktop) -25725.*utop(itop,jtop-1,ktop+1) +3675.*utop(itop,jtop-1,ktop+2) 
					  -77175.*utop(itop,jtop,ktop-1) +1157625.*utop(itop,jtop,ktop) +385875.*utop(itop,jtop,ktop+1) 
					  -55125.*utop(itop,jtop,ktop+2) -25725.*utop(itop,jtop+1,ktop-1) +385875.*utop(itop,jtop+1,ktop) 
					  +128625.*utop(itop,jtop+1,ktop+1) -18375.*utop(itop,jtop+1,ktop+2) +3675.*utop(itop,jtop+2,ktop-1) 
					  -55125.*utop(itop,jtop+2,ktop) -18375.*utop(itop,jtop+2,ktop+1) +2625.*utop(itop,jtop+2,ktop+2) 
					  -343.*utop(itop+1,jtop-1,ktop-1) +5145.*utop(itop+1,jtop-1,ktop) +1715.*utop(itop+1,jtop-1,ktop+1) 
					  -245.*utop(itop+1,jtop-1,ktop+2) +5145.*utop(itop+1,jtop,ktop-1) -77175.*utop(itop+1,jtop,ktop) 
					  -25725.*utop(itop+1,jtop,ktop+1) +3675.*utop(itop+1,jtop,ktop+2) +1715.*utop(itop+1,jtop+1,ktop-1) 
					  -25725.*utop(itop+1,jtop+1,ktop) -8575.*utop(itop+1,jtop+1,ktop+1) +1225.*utop(itop+1,jtop+1,ktop+2) 
					  -245.*utop(itop+1,jtop+2,ktop-1) +3675.*utop(itop+1,jtop+2,ktop) +1225.*utop(itop+1,jtop+2,ktop+1) 
					  -175.*utop(itop+1,jtop+2,ktop+2) )/2097152.;
	u(i+1,j+0,k+0) = ( -175.*utop(itop-1,jtop-2,ktop-2) +1225.*utop(itop-1,jtop-2,ktop-1) +3675.*utop(itop-1,jtop-2,ktop) 
					  -245.*utop(itop-1,jtop-2,ktop+1) +1225.*utop(itop-1,jtop-1,ktop-2) -8575.*utop(itop-1,jtop-1,ktop-1) 
					  -25725.*utop(itop-1,jtop-1,ktop) +1715.*utop(itop-1,jtop-1,ktop+1) +3675.*utop(itop-1,jtop,ktop-2) 
					  -25725.*utop(itop-1,jtop,ktop-1) -77175.*utop(itop-1,jtop,ktop) +5145.*utop(itop-1,jtop,ktop+1) 
					  -245.*utop(itop-1,jtop+1,ktop-2) +1715.*utop(itop-1,jtop+1,ktop-1) +5145.*utop(itop-1,jtop+1,ktop) 
					  -343.*utop(itop-1,jtop+1,ktop+1) +2625.*utop(itop,jtop-2,ktop-2) -18375.*utop(itop,jtop-2,ktop-1) 
					  -55125.*utop(itop,jtop-2,ktop) +3675.*utop(itop,jtop-2,ktop+1) -18375.*utop(itop,jtop-1,ktop-2) 
					  +128625.*utop(itop,jtop-1,ktop-1) +385875.*utop(itop,jtop-1,ktop) -25725.*utop(itop,jtop-1,ktop+1) 
					  -55125.*utop(itop,jtop,ktop-2) +385875.*utop(itop,jtop,ktop-1) +1157625.*utop(itop,jtop,ktop) 
					  -77175.*utop(itop,jtop,ktop+1) +3675.*utop(itop,jtop+1,ktop-2) -25725.*utop(itop,jtop+1,ktop-1) 
					  -77175.*utop(itop,jtop+1,ktop) +5145.*utop(itop,jtop+1,ktop+1) +875.*utop(itop+1,jtop-2,ktop-2) 
					  -6125.*utop(itop+1,jtop-2,ktop-1) -18375.*utop(itop+1,jtop-2,ktop) +1225.*utop(itop+1,jtop-2,ktop+1) 
					  -6125.*utop(itop+1,jtop-1,ktop-2) +42875.*utop(itop+1,jtop-1,ktop-1) +128625.*utop(itop+1,jtop-1,ktop) 
					  -8575.*utop(itop+1,jtop-1,ktop+1) -18375.*utop(itop+1,jtop,ktop-2) +128625.*utop(itop+1,jtop,ktop-1) 
					  +385875.*utop(itop+1,jtop,ktop) -25725.*utop(itop+1,jtop,ktop+1) +1225.*utop(itop+1,jtop+1,ktop-2) 
					  -8575.*utop(itop+1,jtop+1,ktop-1) -25725.*utop(itop+1,jtop+1,ktop) +1715.*utop(itop+1,jtop+1,ktop+1) 
					  -125.*utop(itop+2,jtop-2,ktop-2) +875.*utop(itop+2,jtop-2,ktop-1) +2625.*utop(itop+2,jtop-2,ktop) 
					  -175.*utop(itop+2,jtop-2,ktop+1) +875.*utop(itop+2,jtop-1,ktop-2) -6125.*utop(itop+2,jtop-1,ktop-1) 
					  -18375.*utop(itop+2,jtop-1,ktop) +1225.*utop(itop+2,jtop-1,ktop+1) +2625.*utop(itop+2,jtop,ktop-2) 
					  -18375.*utop(itop+2,jtop,ktop-1) -55125.*utop(itop+2,jtop,ktop) +3675.*utop(itop+2,jtop,ktop+1) 
					  -175.*utop(itop+2,jtop+1,ktop-2) +1225.*utop(itop+2,jtop+1,ktop-1) +3675.*utop(itop+2,jtop+1,ktop) 
					  -245.*utop(itop+2,jtop+1,ktop+1) )/2097152.;
	u(i+1,j+0,k+1) = ( -245.*utop(itop-1,jtop-2,ktop-1) +3675.*utop(itop-1,jtop-2,ktop) +1225.*utop(itop-1,jtop-2,ktop+1) -175.*utop(itop-1,jtop-2,ktop+2) +1715.*utop(itop-1,jtop-1,ktop-1) -25725.*utop(itop-1,jtop-1,ktop) -8575.*utop(itop-1,jtop-1,ktop+1) +1225.*utop(itop-1,jtop-1,ktop+2) +5145.*utop(itop-1,jtop,ktop-1) -77175.*utop(itop-1,jtop,ktop) -25725.*utop(itop-1,jtop,ktop+1) +3675.*utop(itop-1,jtop,ktop+2) -343.*utop(itop-1,jtop+1,ktop-1) +5145.*utop(itop-1,jtop+1,ktop) +1715.*utop(itop-1,jtop+1,ktop+1) -245.*utop(itop-1,jtop+1,ktop+2) +3675.*utop(itop,jtop-2,ktop-1) -55125.*utop(itop,jtop-2,ktop) -18375.*utop(itop,jtop-2,ktop+1) +2625.*utop(itop,jtop-2,ktop+2) -25725.*utop(itop,jtop-1,ktop-1) +385875.*utop(itop,jtop-1,ktop) +128625.*utop(itop,jtop-1,ktop+1) -18375.*utop(itop,jtop-1,ktop+2) -77175.*utop(itop,jtop,ktop-1) +1157625.*utop(itop,jtop,ktop) +385875.*utop(itop,jtop,ktop+1) -55125.*utop(itop,jtop,ktop+2) +5145.*utop(itop,jtop+1,ktop-1) -77175.*utop(itop,jtop+1,ktop) -25725.*utop(itop,jtop+1,ktop+1) +3675.*utop(itop,jtop+1,ktop+2) +1225.*utop(itop+1,jtop-2,ktop-1) -18375.*utop(itop+1,jtop-2,ktop) -6125.*utop(itop+1,jtop-2,ktop+1) +875.*utop(itop+1,jtop-2,ktop+2) -8575.*utop(itop+1,jtop-1,ktop-1) +128625.*utop(itop+1,jtop-1,ktop) +42875.*utop(itop+1,jtop-1,ktop+1) -6125.*utop(itop+1,jtop-1,ktop+2) -25725.*utop(itop+1,jtop,ktop-1) +385875.*utop(itop+1,jtop,ktop) +128625.*utop(itop+1,jtop,ktop+1) -18375.*utop(itop+1,jtop,ktop+2) +1715.*utop(itop+1,jtop+1,ktop-1) -25725.*utop(itop+1,jtop+1,ktop) -8575.*utop(itop+1,jtop+1,ktop+1) +1225.*utop(itop+1,jtop+1,ktop+2) -175.*utop(itop+2,jtop-2,ktop-1) +2625.*utop(itop+2,jtop-2,ktop) +875.*utop(itop+2,jtop-2,ktop+1) -125.*utop(itop+2,jtop-2,ktop+2) +1225.*utop(itop+2,jtop-1,ktop-1) -18375.*utop(itop+2,jtop-1,ktop) -6125.*utop(itop+2,jtop-1,ktop+1) +875.*utop(itop+2,jtop-1,ktop+2) +3675.*utop(itop+2,jtop,ktop-1) -55125.*utop(itop+2,jtop,ktop) -18375.*utop(itop+2,jtop,ktop+1) +2625.*utop(itop+2,jtop,ktop+2) -245.*utop(itop+2,jtop+1,ktop-1) +3675.*utop(itop+2,jtop+1,ktop) +1225.*utop(itop+2,jtop+1,ktop+1) -175.*utop(itop+2,jtop+1,ktop+2) )/2097152.;
	u(i+1,j+1,k+0) = ( -245.*utop(itop-1,jtop-1,ktop-2) +1715.*utop(itop-1,jtop-1,ktop-1) +5145.*utop(itop-1,jtop-1,ktop) -343.*utop(itop-1,jtop-1,ktop+1) +3675.*utop(itop-1,jtop,ktop-2) -25725.*utop(itop-1,jtop,ktop-1) -77175.*utop(itop-1,jtop,ktop) +5145.*utop(itop-1,jtop,ktop+1) +1225.*utop(itop-1,jtop+1,ktop-2) -8575.*utop(itop-1,jtop+1,ktop-1) -25725.*utop(itop-1,jtop+1,ktop) +1715.*utop(itop-1,jtop+1,ktop+1) -175.*utop(itop-1,jtop+2,ktop-2) +1225.*utop(itop-1,jtop+2,ktop-1) +3675.*utop(itop-1,jtop+2,ktop) -245.*utop(itop-1,jtop+2,ktop+1) +3675.*utop(itop,jtop-1,ktop-2) -25725.*utop(itop,jtop-1,ktop-1) -77175.*utop(itop,jtop-1,ktop) +5145.*utop(itop,jtop-1,ktop+1) -55125.*utop(itop,jtop,ktop-2) +385875.*utop(itop,jtop,ktop-1) +1157625.*utop(itop,jtop,ktop) -77175.*utop(itop,jtop,ktop+1) -18375.*utop(itop,jtop+1,ktop-2) +128625.*utop(itop,jtop+1,ktop-1) +385875.*utop(itop,jtop+1,ktop) -25725.*utop(itop,jtop+1,ktop+1) +2625.*utop(itop,jtop+2,ktop-2) -18375.*utop(itop,jtop+2,ktop-1) -55125.*utop(itop,jtop+2,ktop) +3675.*utop(itop,jtop+2,ktop+1) +1225.*utop(itop+1,jtop-1,ktop-2) -8575.*utop(itop+1,jtop-1,ktop-1) -25725.*utop(itop+1,jtop-1,ktop) +1715.*utop(itop+1,jtop-1,ktop+1) -18375.*utop(itop+1,jtop,ktop-2) +128625.*utop(itop+1,jtop,ktop-1) +385875.*utop(itop+1,jtop,ktop) -25725.*utop(itop+1,jtop,ktop+1) -6125.*utop(itop+1,jtop+1,ktop-2) +42875.*utop(itop+1,jtop+1,ktop-1) +128625.*utop(itop+1,jtop+1,ktop) -8575.*utop(itop+1,jtop+1,ktop+1) +875.*utop(itop+1,jtop+2,ktop-2) -6125.*utop(itop+1,jtop+2,ktop-1) -18375.*utop(itop+1,jtop+2,ktop) +1225.*utop(itop+1,jtop+2,ktop+1) -175.*utop(itop+2,jtop-1,ktop-2) +1225.*utop(itop+2,jtop-1,ktop-1) +3675.*utop(itop+2,jtop-1,ktop) -245.*utop(itop+2,jtop-1,ktop+1) +2625.*utop(itop+2,jtop,ktop-2) -18375.*utop(itop+2,jtop,ktop-1) -55125.*utop(itop+2,jtop,ktop) +3675.*utop(itop+2,jtop,ktop+1) +875.*utop(itop+2,jtop+1,ktop-2) -6125.*utop(itop+2,jtop+1,ktop-1) -18375.*utop(itop+2,jtop+1,ktop) +1225.*utop(itop+2,jtop+1,ktop+1) -125.*utop(itop+2,jtop+2,ktop-2) +875.*utop(itop+2,jtop+2,ktop-1) +2625.*utop(itop+2,jtop+2,ktop) -175.*utop(itop+2,jtop+2,ktop+1) )/2097152.;
	u(i+1,j+1,k+1) = ( -343.*utop(itop-1,jtop-1,ktop-1) +5145.*utop(itop-1,jtop-1,ktop) +1715.*utop(itop-1,jtop-1,ktop+1) -245.*utop(itop-1,jtop-1,ktop+2) +5145.*utop(itop-1,jtop,ktop-1) -77175.*utop(itop-1,jtop,ktop) -25725.*utop(itop-1,jtop,ktop+1) +3675.*utop(itop-1,jtop,ktop+2) +1715.*utop(itop-1,jtop+1,ktop-1) -25725.*utop(itop-1,jtop+1,ktop) -8575.*utop(itop-1,jtop+1,ktop+1) +1225.*utop(itop-1,jtop+1,ktop+2) -245.*utop(itop-1,jtop+2,ktop-1) +3675.*utop(itop-1,jtop+2,ktop) +1225.*utop(itop-1,jtop+2,ktop+1) -175.*utop(itop-1,jtop+2,ktop+2) +5145.*utop(itop,jtop-1,ktop-1) -77175.*utop(itop,jtop-1,ktop) -25725.*utop(itop,jtop-1,ktop+1) +3675.*utop(itop,jtop-1,ktop+2) -77175.*utop(itop,jtop,ktop-1) +1157625.*utop(itop,jtop,ktop) +385875.*utop(itop,jtop,ktop+1) -55125.*utop(itop,jtop,ktop+2) -25725.*utop(itop,jtop+1,ktop-1) +385875.*utop(itop,jtop+1,ktop) +128625.*utop(itop,jtop+1,ktop+1) -18375.*utop(itop,jtop+1,ktop+2) +3675.*utop(itop,jtop+2,ktop-1) -55125.*utop(itop,jtop+2,ktop) -18375.*utop(itop,jtop+2,ktop+1) +2625.*utop(itop,jtop+2,ktop+2) +1715.*utop(itop+1,jtop-1,ktop-1) -25725.*utop(itop+1,jtop-1,ktop) -8575.*utop(itop+1,jtop-1,ktop+1) +1225.*utop(itop+1,jtop-1,ktop+2) -25725.*utop(itop+1,jtop,ktop-1) +385875.*utop(itop+1,jtop,ktop) +128625.*utop(itop+1,jtop,ktop+1) -18375.*utop(itop+1,jtop,ktop+2) -8575.*utop(itop+1,jtop+1,ktop-1) +128625.*utop(itop+1,jtop+1,ktop) +42875.*utop(itop+1,jtop+1,ktop+1) -6125.*utop(itop+1,jtop+1,ktop+2) +1225.*utop(itop+1,jtop+2,ktop-1) -18375.*utop(itop+1,jtop+2,ktop) -6125.*utop(itop+1,jtop+2,ktop+1) +875.*utop(itop+1,jtop+2,ktop+2) -245.*utop(itop+2,jtop-1,ktop-1) +3675.*utop(itop+2,jtop-1,ktop) +1225.*utop(itop+2,jtop-1,ktop+1) -175.*utop(itop+2,jtop-1,ktop+2) +3675.*utop(itop+2,jtop,ktop-1) -55125.*utop(itop+2,jtop,ktop) -18375.*utop(itop+2,jtop,ktop+1) +2625.*utop(itop+2,jtop,ktop+2) +1225.*utop(itop+2,jtop+1,ktop-1) -18375.*utop(itop+2,jtop+1,ktop) -6125.*utop(itop+2,jtop+1,ktop+1) +875.*utop(itop+2,jtop+1,ktop+2) -175.*utop(itop+2,jtop+2,ktop-1) +2625.*utop(itop+2,jtop+2,ktop) +875.*utop(itop+2,jtop+2,ktop+1) -125.*utop(itop+2,jtop+2,ktop+2) )/2097152.;
	*/
	
	u(i+0,j+0,k+0) = ( -1.060835e-05*utop(itop-2,jtop-2,ktop-2) +9.901123e-05*utop(itop-2,jtop-2,ktop-1) +4.455505e-04*utop(itop-2,jtop-2,ktop) -5.940674e-05*utop(itop-2,jtop-2,ktop+1) +8.250936e-06*utop(itop-2,jtop-2,ktop+2) +9.901123e-05*utop(itop-2,jtop-1,ktop-2) -9.241048e-04*utop(itop-2,jtop-1,ktop-1) -4.158472e-03*utop(itop-2,jtop-1,ktop) +5.544629e-04*utop(itop-2,jtop-1,ktop+1) -7.700874e-05*utop(itop-2,jtop-1,ktop+2) +4.455505e-04*utop(itop-2,jtop,ktop-2) -4.158472e-03*utop(itop-2,jtop,ktop-1) -1.871312e-02*utop(itop-2,jtop,ktop) +2.495083e-03*utop(itop-2,jtop,ktop+1) -3.465393e-04*utop(itop-2,jtop,ktop+2) -5.940674e-05*utop(itop-2,jtop+1,ktop-2) +5.544629e-04*utop(itop-2,jtop+1,ktop-1) +2.495083e-03*utop(itop-2,jtop+1,ktop) -3.326777e-04*utop(itop-2,jtop+1,ktop+1) +4.620524e-05*utop(itop-2,jtop+1,ktop+2) +8.250936e-06*utop(itop-2,jtop+2,ktop-2) -7.700874e-05*utop(itop-2,jtop+2,ktop-1) -3.465393e-04*utop(itop-2,jtop+2,ktop) +4.620524e-05*utop(itop-2,jtop+2,ktop+1) -6.417395e-06*utop(itop-2,jtop+2,ktop+2) +9.901123e-05*utop(itop-1,jtop-2,ktop-2) -9.241048e-04*utop(itop-1,jtop-2,ktop-1) -4.158472e-03*utop(itop-1,jtop-2,ktop) +5.544629e-04*utop(itop-1,jtop-2,ktop+1) -7.700874e-05*utop(itop-1,jtop-2,ktop+2) -9.241048e-04*utop(itop-1,jtop-1,ktop-2) +8.624978e-03*utop(itop-1,jtop-1,ktop-1) +3.881240e-02*utop(itop-1,jtop-1,ktop) -5.174987e-03*utop(itop-1,jtop-1,ktop+1) +7.187482e-04*utop(itop-1,jtop-1,ktop+2) -4.158472e-03*utop(itop-1,jtop,ktop-2) +3.881240e-02*utop(itop-1,jtop,ktop-1) +1.746558e-01*utop(itop-1,jtop,ktop) -2.328744e-02*utop(itop-1,jtop,ktop+1) +3.234367e-03*utop(itop-1,jtop,ktop+2) +5.544629e-04*utop(itop-1,jtop+1,ktop-2) -5.174987e-03*utop(itop-1,jtop+1,ktop-1) -2.328744e-02*utop(itop-1,jtop+1,ktop) +3.104992e-03*utop(itop-1,jtop+1,ktop+1) -4.312489e-04*utop(itop-1,jtop+1,ktop+2) -7.700874e-05*utop(itop-1,jtop+2,ktop-2) +7.187482e-04*utop(itop-1,jtop+2,ktop-1) +3.234367e-03*utop(itop-1,jtop+2,ktop) -4.312489e-04*utop(itop-1,jtop+2,ktop+1) +5.989568e-05*utop(itop-1,jtop+2,ktop+2) +4.455505e-04*utop(itop,jtop-2,ktop-2) -4.158472e-03*utop(itop,jtop-2,ktop-1) -1.871312e-02*utop(itop,jtop-2,ktop) +2.495083e-03*utop(itop,jtop-2,ktop+1) -3.465393e-04*utop(itop,jtop-2,ktop+2) -4.158472e-03*utop(itop,jtop-1,ktop-2) +3.881240e-02*utop(itop,jtop-1,ktop-1) +1.746558e-01*utop(itop,jtop-1,ktop) -2.328744e-02*utop(itop,jtop-1,ktop+1) +3.234367e-03*utop(itop,jtop-1,ktop+2) -1.871312e-02*utop(itop,jtop,ktop-2) +1.746558e-01*utop(itop,jtop,ktop-1) +7.859512e-01*utop(itop,jtop,ktop) -1.047935e-01*utop(itop,jtop,ktop+1) +1.455465e-02*utop(itop,jtop,ktop+2) +2.495083e-03*utop(itop,jtop+1,ktop-2) -2.328744e-02*utop(itop,jtop+1,ktop-1) -1.047935e-01*utop(itop,jtop+1,ktop) +1.397246e-02*utop(itop,jtop+1,ktop+1) -1.940620e-03*utop(itop,jtop+1,ktop+2) -3.465393e-04*utop(itop,jtop+2,ktop-2) +3.234367e-03*utop(itop,jtop+2,ktop-1) +1.455465e-02*utop(itop,jtop+2,ktop) -1.940620e-03*utop(itop,jtop+2,ktop+1) +2.695306e-04*utop(itop,jtop+2,ktop+2) -5.940674e-05*utop(itop+1,jtop-2,ktop-2) +5.544629e-04*utop(itop+1,jtop-2,ktop-1) +2.495083e-03*utop(itop+1,jtop-2,ktop) -3.326777e-04*utop(itop+1,jtop-2,ktop+1) +4.620524e-05*utop(itop+1,jtop-2,ktop+2) +5.544629e-04*utop(itop+1,jtop-1,ktop-2) -5.174987e-03*utop(itop+1,jtop-1,ktop-1) -2.328744e-02*utop(itop+1,jtop-1,ktop) +3.104992e-03*utop(itop+1,jtop-1,ktop+1) -4.312489e-04*utop(itop+1,jtop-1,ktop+2) +2.495083e-03*utop(itop+1,jtop,ktop-2) -2.328744e-02*utop(itop+1,jtop,ktop-1) -1.047935e-01*utop(itop+1,jtop,ktop) +1.397246e-02*utop(itop+1,jtop,ktop+1) -1.940620e-03*utop(itop+1,jtop,ktop+2) -3.326777e-04*utop(itop+1,jtop+1,ktop-2) +3.104992e-03*utop(itop+1,jtop+1,ktop-1) +1.397246e-02*utop(itop+1,jtop+1,ktop) -1.862995e-03*utop(itop+1,jtop+1,ktop+1) +2.587494e-04*utop(itop+1,jtop+1,ktop+2) +4.620524e-05*utop(itop+1,jtop+2,ktop-2) -4.312489e-04*utop(itop+1,jtop+2,ktop-1) -1.940620e-03*utop(itop+1,jtop+2,ktop) +2.587494e-04*utop(itop+1,jtop+2,ktop+1) -3.593741e-05*utop(itop+1,jtop+2,ktop+2) +8.250936e-06*utop(itop+2,jtop-2,ktop-2) -7.700874e-05*utop(itop+2,jtop-2,ktop-1) -3.465393e-04*utop(itop+2,jtop-2,ktop) +4.620524e-05*utop(itop+2,jtop-2,ktop+1) -6.417395e-06*utop(itop+2,jtop-2,ktop+2) -7.700874e-05*utop(itop+2,jtop-1,ktop-2) +7.187482e-04*utop(itop+2,jtop-1,ktop-1) +3.234367e-03*utop(itop+2,jtop-1,ktop) -4.312489e-04*utop(itop+2,jtop-1,ktop+1) +5.989568e-05*utop(itop+2,jtop-1,ktop+2) -3.465393e-04*utop(itop+2,jtop,ktop-2) +3.234367e-03*utop(itop+2,jtop,ktop-1) +1.455465e-02*utop(itop+2,jtop,ktop) -1.940620e-03*utop(itop+2,jtop,ktop+1) +2.695306e-04*utop(itop+2,jtop,ktop+2) +4.620524e-05*utop(itop+2,jtop+1,ktop-2) -4.312489e-04*utop(itop+2,jtop+1,ktop-1) -1.940620e-03*utop(itop+2,jtop+1,ktop) +2.587494e-04*utop(itop+2,jtop+1,ktop+1) -3.593741e-05*utop(itop+2,jtop+1,ktop+2) -6.417395e-06*utop(itop+2,jtop+2,ktop-2) +5.989568e-05*utop(itop+2,jtop+2,ktop-1) +2.695306e-04*utop(itop+2,jtop+2,ktop) -3.593741e-05*utop(itop+2,jtop+2,ktop+1) +4.991307e-06*utop(itop+2,jtop+2,ktop+2));
	u(i+0,j+0,k+1) = ( +8.250936e-06*utop(itop-2,jtop-2,ktop-2) -5.940674e-05*utop(itop-2,jtop-2,ktop-1) +4.455505e-04*utop(itop-2,jtop-2,ktop) +9.901123e-05*utop(itop-2,jtop-2,ktop+1) -1.060835e-05*utop(itop-2,jtop-2,ktop+2) -7.700874e-05*utop(itop-2,jtop-1,ktop-2) +5.544629e-04*utop(itop-2,jtop-1,ktop-1) -4.158472e-03*utop(itop-2,jtop-1,ktop) -9.241048e-04*utop(itop-2,jtop-1,ktop+1) +9.901123e-05*utop(itop-2,jtop-1,ktop+2) -3.465393e-04*utop(itop-2,jtop,ktop-2) +2.495083e-03*utop(itop-2,jtop,ktop-1) -1.871312e-02*utop(itop-2,jtop,ktop) -4.158472e-03*utop(itop-2,jtop,ktop+1) +4.455505e-04*utop(itop-2,jtop,ktop+2) +4.620524e-05*utop(itop-2,jtop+1,ktop-2) -3.326777e-04*utop(itop-2,jtop+1,ktop-1) +2.495083e-03*utop(itop-2,jtop+1,ktop) +5.544629e-04*utop(itop-2,jtop+1,ktop+1) -5.940674e-05*utop(itop-2,jtop+1,ktop+2) -6.417395e-06*utop(itop-2,jtop+2,ktop-2) +4.620524e-05*utop(itop-2,jtop+2,ktop-1) -3.465393e-04*utop(itop-2,jtop+2,ktop) -7.700874e-05*utop(itop-2,jtop+2,ktop+1) +8.250936e-06*utop(itop-2,jtop+2,ktop+2) -7.700874e-05*utop(itop-1,jtop-2,ktop-2) +5.544629e-04*utop(itop-1,jtop-2,ktop-1) -4.158472e-03*utop(itop-1,jtop-2,ktop) -9.241048e-04*utop(itop-1,jtop-2,ktop+1) +9.901123e-05*utop(itop-1,jtop-2,ktop+2) +7.187482e-04*utop(itop-1,jtop-1,ktop-2) -5.174987e-03*utop(itop-1,jtop-1,ktop-1) +3.881240e-02*utop(itop-1,jtop-1,ktop) +8.624978e-03*utop(itop-1,jtop-1,ktop+1) -9.241048e-04*utop(itop-1,jtop-1,ktop+2) +3.234367e-03*utop(itop-1,jtop,ktop-2) -2.328744e-02*utop(itop-1,jtop,ktop-1) +1.746558e-01*utop(itop-1,jtop,ktop) +3.881240e-02*utop(itop-1,jtop,ktop+1) -4.158472e-03*utop(itop-1,jtop,ktop+2) -4.312489e-04*utop(itop-1,jtop+1,ktop-2) +3.104992e-03*utop(itop-1,jtop+1,ktop-1) -2.328744e-02*utop(itop-1,jtop+1,ktop) -5.174987e-03*utop(itop-1,jtop+1,ktop+1) +5.544629e-04*utop(itop-1,jtop+1,ktop+2) +5.989568e-05*utop(itop-1,jtop+2,ktop-2) -4.312489e-04*utop(itop-1,jtop+2,ktop-1) +3.234367e-03*utop(itop-1,jtop+2,ktop) +7.187482e-04*utop(itop-1,jtop+2,ktop+1) -7.700874e-05*utop(itop-1,jtop+2,ktop+2) -3.465393e-04*utop(itop,jtop-2,ktop-2) +2.495083e-03*utop(itop,jtop-2,ktop-1) -1.871312e-02*utop(itop,jtop-2,ktop) -4.158472e-03*utop(itop,jtop-2,ktop+1) +4.455505e-04*utop(itop,jtop-2,ktop+2) +3.234367e-03*utop(itop,jtop-1,ktop-2) -2.328744e-02*utop(itop,jtop-1,ktop-1) +1.746558e-01*utop(itop,jtop-1,ktop) +3.881240e-02*utop(itop,jtop-1,ktop+1) -4.158472e-03*utop(itop,jtop-1,ktop+2) +1.455465e-02*utop(itop,jtop,ktop-2) -1.047935e-01*utop(itop,jtop,ktop-1) +7.859512e-01*utop(itop,jtop,ktop) +1.746558e-01*utop(itop,jtop,ktop+1) -1.871312e-02*utop(itop,jtop,ktop+2) -1.940620e-03*utop(itop,jtop+1,ktop-2) +1.397246e-02*utop(itop,jtop+1,ktop-1) -1.047935e-01*utop(itop,jtop+1,ktop) -2.328744e-02*utop(itop,jtop+1,ktop+1) +2.495083e-03*utop(itop,jtop+1,ktop+2) +2.695306e-04*utop(itop,jtop+2,ktop-2) -1.940620e-03*utop(itop,jtop+2,ktop-1) +1.455465e-02*utop(itop,jtop+2,ktop) +3.234367e-03*utop(itop,jtop+2,ktop+1) -3.465393e-04*utop(itop,jtop+2,ktop+2) +4.620524e-05*utop(itop+1,jtop-2,ktop-2) -3.326777e-04*utop(itop+1,jtop-2,ktop-1) +2.495083e-03*utop(itop+1,jtop-2,ktop) +5.544629e-04*utop(itop+1,jtop-2,ktop+1) -5.940674e-05*utop(itop+1,jtop-2,ktop+2) -4.312489e-04*utop(itop+1,jtop-1,ktop-2) +3.104992e-03*utop(itop+1,jtop-1,ktop-1) -2.328744e-02*utop(itop+1,jtop-1,ktop) -5.174987e-03*utop(itop+1,jtop-1,ktop+1) +5.544629e-04*utop(itop+1,jtop-1,ktop+2) -1.940620e-03*utop(itop+1,jtop,ktop-2) +1.397246e-02*utop(itop+1,jtop,ktop-1) -1.047935e-01*utop(itop+1,jtop,ktop) -2.328744e-02*utop(itop+1,jtop,ktop+1) +2.495083e-03*utop(itop+1,jtop,ktop+2) +2.587494e-04*utop(itop+1,jtop+1,ktop-2) -1.862995e-03*utop(itop+1,jtop+1,ktop-1) +1.397246e-02*utop(itop+1,jtop+1,ktop) +3.104992e-03*utop(itop+1,jtop+1,ktop+1) -3.326777e-04*utop(itop+1,jtop+1,ktop+2) -3.593741e-05*utop(itop+1,jtop+2,ktop-2) +2.587494e-04*utop(itop+1,jtop+2,ktop-1) -1.940620e-03*utop(itop+1,jtop+2,ktop) -4.312489e-04*utop(itop+1,jtop+2,ktop+1) +4.620524e-05*utop(itop+1,jtop+2,ktop+2) -6.417395e-06*utop(itop+2,jtop-2,ktop-2) +4.620524e-05*utop(itop+2,jtop-2,ktop-1) -3.465393e-04*utop(itop+2,jtop-2,ktop) -7.700874e-05*utop(itop+2,jtop-2,ktop+1) +8.250936e-06*utop(itop+2,jtop-2,ktop+2) +5.989568e-05*utop(itop+2,jtop-1,ktop-2) -4.312489e-04*utop(itop+2,jtop-1,ktop-1) +3.234367e-03*utop(itop+2,jtop-1,ktop) +7.187482e-04*utop(itop+2,jtop-1,ktop+1) -7.700874e-05*utop(itop+2,jtop-1,ktop+2) +2.695306e-04*utop(itop+2,jtop,ktop-2) -1.940620e-03*utop(itop+2,jtop,ktop-1) +1.455465e-02*utop(itop+2,jtop,ktop) +3.234367e-03*utop(itop+2,jtop,ktop+1) -3.465393e-04*utop(itop+2,jtop,ktop+2) -3.593741e-05*utop(itop+2,jtop+1,ktop-2) +2.587494e-04*utop(itop+2,jtop+1,ktop-1) -1.940620e-03*utop(itop+2,jtop+1,ktop) -4.312489e-04*utop(itop+2,jtop+1,ktop+1) +4.620524e-05*utop(itop+2,jtop+1,ktop+2) +4.991307e-06*utop(itop+2,jtop+2,ktop-2) -3.593741e-05*utop(itop+2,jtop+2,ktop-1) +2.695306e-04*utop(itop+2,jtop+2,ktop) +5.989568e-05*utop(itop+2,jtop+2,ktop+1) -6.417395e-06*utop(itop+2,jtop+2,ktop+2));
	u(i+0,j+1,k+0) = ( +8.250936e-06*utop(itop-2,jtop-2,ktop-2) -7.700874e-05*utop(itop-2,jtop-2,ktop-1) -3.465393e-04*utop(itop-2,jtop-2,ktop) +4.620524e-05*utop(itop-2,jtop-2,ktop+1) -6.417395e-06*utop(itop-2,jtop-2,ktop+2) -5.940674e-05*utop(itop-2,jtop-1,ktop-2) +5.544629e-04*utop(itop-2,jtop-1,ktop-1) +2.495083e-03*utop(itop-2,jtop-1,ktop) -3.326777e-04*utop(itop-2,jtop-1,ktop+1) +4.620524e-05*utop(itop-2,jtop-1,ktop+2) +4.455505e-04*utop(itop-2,jtop,ktop-2) -4.158472e-03*utop(itop-2,jtop,ktop-1) -1.871312e-02*utop(itop-2,jtop,ktop) +2.495083e-03*utop(itop-2,jtop,ktop+1) -3.465393e-04*utop(itop-2,jtop,ktop+2) +9.901123e-05*utop(itop-2,jtop+1,ktop-2) -9.241048e-04*utop(itop-2,jtop+1,ktop-1) -4.158472e-03*utop(itop-2,jtop+1,ktop) +5.544629e-04*utop(itop-2,jtop+1,ktop+1) -7.700874e-05*utop(itop-2,jtop+1,ktop+2) -1.060835e-05*utop(itop-2,jtop+2,ktop-2) +9.901123e-05*utop(itop-2,jtop+2,ktop-1) +4.455505e-04*utop(itop-2,jtop+2,ktop) -5.940674e-05*utop(itop-2,jtop+2,ktop+1) +8.250936e-06*utop(itop-2,jtop+2,ktop+2) -7.700874e-05*utop(itop-1,jtop-2,ktop-2) +7.187482e-04*utop(itop-1,jtop-2,ktop-1) +3.234367e-03*utop(itop-1,jtop-2,ktop) -4.312489e-04*utop(itop-1,jtop-2,ktop+1) +5.989568e-05*utop(itop-1,jtop-2,ktop+2) +5.544629e-04*utop(itop-1,jtop-1,ktop-2) -5.174987e-03*utop(itop-1,jtop-1,ktop-1) -2.328744e-02*utop(itop-1,jtop-1,ktop) +3.104992e-03*utop(itop-1,jtop-1,ktop+1) -4.312489e-04*utop(itop-1,jtop-1,ktop+2) -4.158472e-03*utop(itop-1,jtop,ktop-2) +3.881240e-02*utop(itop-1,jtop,ktop-1) +1.746558e-01*utop(itop-1,jtop,ktop) -2.328744e-02*utop(itop-1,jtop,ktop+1) +3.234367e-03*utop(itop-1,jtop,ktop+2) -9.241048e-04*utop(itop-1,jtop+1,ktop-2) +8.624978e-03*utop(itop-1,jtop+1,ktop-1) +3.881240e-02*utop(itop-1,jtop+1,ktop) -5.174987e-03*utop(itop-1,jtop+1,ktop+1) +7.187482e-04*utop(itop-1,jtop+1,ktop+2) +9.901123e-05*utop(itop-1,jtop+2,ktop-2) -9.241048e-04*utop(itop-1,jtop+2,ktop-1) -4.158472e-03*utop(itop-1,jtop+2,ktop) +5.544629e-04*utop(itop-1,jtop+2,ktop+1) -7.700874e-05*utop(itop-1,jtop+2,ktop+2) -3.465393e-04*utop(itop,jtop-2,ktop-2) +3.234367e-03*utop(itop,jtop-2,ktop-1) +1.455465e-02*utop(itop,jtop-2,ktop) -1.940620e-03*utop(itop,jtop-2,ktop+1) +2.695306e-04*utop(itop,jtop-2,ktop+2) +2.495083e-03*utop(itop,jtop-1,ktop-2) -2.328744e-02*utop(itop,jtop-1,ktop-1) -1.047935e-01*utop(itop,jtop-1,ktop) +1.397246e-02*utop(itop,jtop-1,ktop+1) -1.940620e-03*utop(itop,jtop-1,ktop+2) -1.871312e-02*utop(itop,jtop,ktop-2) +1.746558e-01*utop(itop,jtop,ktop-1) +7.859512e-01*utop(itop,jtop,ktop) -1.047935e-01*utop(itop,jtop,ktop+1) +1.455465e-02*utop(itop,jtop,ktop+2) -4.158472e-03*utop(itop,jtop+1,ktop-2) +3.881240e-02*utop(itop,jtop+1,ktop-1) +1.746558e-01*utop(itop,jtop+1,ktop) -2.328744e-02*utop(itop,jtop+1,ktop+1) +3.234367e-03*utop(itop,jtop+1,ktop+2) +4.455505e-04*utop(itop,jtop+2,ktop-2) -4.158472e-03*utop(itop,jtop+2,ktop-1) -1.871312e-02*utop(itop,jtop+2,ktop) +2.495083e-03*utop(itop,jtop+2,ktop+1) -3.465393e-04*utop(itop,jtop+2,ktop+2) +4.620524e-05*utop(itop+1,jtop-2,ktop-2) -4.312489e-04*utop(itop+1,jtop-2,ktop-1) -1.940620e-03*utop(itop+1,jtop-2,ktop) +2.587494e-04*utop(itop+1,jtop-2,ktop+1) -3.593741e-05*utop(itop+1,jtop-2,ktop+2) -3.326777e-04*utop(itop+1,jtop-1,ktop-2) +3.104992e-03*utop(itop+1,jtop-1,ktop-1) +1.397246e-02*utop(itop+1,jtop-1,ktop) -1.862995e-03*utop(itop+1,jtop-1,ktop+1) +2.587494e-04*utop(itop+1,jtop-1,ktop+2) +2.495083e-03*utop(itop+1,jtop,ktop-2) -2.328744e-02*utop(itop+1,jtop,ktop-1) -1.047935e-01*utop(itop+1,jtop,ktop) +1.397246e-02*utop(itop+1,jtop,ktop+1) -1.940620e-03*utop(itop+1,jtop,ktop+2) +5.544629e-04*utop(itop+1,jtop+1,ktop-2) -5.174987e-03*utop(itop+1,jtop+1,ktop-1) -2.328744e-02*utop(itop+1,jtop+1,ktop) +3.104992e-03*utop(itop+1,jtop+1,ktop+1) -4.312489e-04*utop(itop+1,jtop+1,ktop+2) -5.940674e-05*utop(itop+1,jtop+2,ktop-2) +5.544629e-04*utop(itop+1,jtop+2,ktop-1) +2.495083e-03*utop(itop+1,jtop+2,ktop) -3.326777e-04*utop(itop+1,jtop+2,ktop+1) +4.620524e-05*utop(itop+1,jtop+2,ktop+2) -6.417395e-06*utop(itop+2,jtop-2,ktop-2) +5.989568e-05*utop(itop+2,jtop-2,ktop-1) +2.695306e-04*utop(itop+2,jtop-2,ktop) -3.593741e-05*utop(itop+2,jtop-2,ktop+1) +4.991307e-06*utop(itop+2,jtop-2,ktop+2) +4.620524e-05*utop(itop+2,jtop-1,ktop-2) -4.312489e-04*utop(itop+2,jtop-1,ktop-1) -1.940620e-03*utop(itop+2,jtop-1,ktop) +2.587494e-04*utop(itop+2,jtop-1,ktop+1) -3.593741e-05*utop(itop+2,jtop-1,ktop+2) -3.465393e-04*utop(itop+2,jtop,ktop-2) +3.234367e-03*utop(itop+2,jtop,ktop-1) +1.455465e-02*utop(itop+2,jtop,ktop) -1.940620e-03*utop(itop+2,jtop,ktop+1) +2.695306e-04*utop(itop+2,jtop,ktop+2) -7.700874e-05*utop(itop+2,jtop+1,ktop-2) +7.187482e-04*utop(itop+2,jtop+1,ktop-1) +3.234367e-03*utop(itop+2,jtop+1,ktop) -4.312489e-04*utop(itop+2,jtop+1,ktop+1) +5.989568e-05*utop(itop+2,jtop+1,ktop+2) +8.250936e-06*utop(itop+2,jtop+2,ktop-2) -7.700874e-05*utop(itop+2,jtop+2,ktop-1) -3.465393e-04*utop(itop+2,jtop+2,ktop) +4.620524e-05*utop(itop+2,jtop+2,ktop+1) -6.417395e-06*utop(itop+2,jtop+2,ktop+2));
	u(i+0,j+1,k+1) = ( -6.417395e-06*utop(itop-2,jtop-2,ktop-2) +4.620524e-05*utop(itop-2,jtop-2,ktop-1) -3.465393e-04*utop(itop-2,jtop-2,ktop) -7.700874e-05*utop(itop-2,jtop-2,ktop+1) +8.250936e-06*utop(itop-2,jtop-2,ktop+2) +4.620524e-05*utop(itop-2,jtop-1,ktop-2) -3.326777e-04*utop(itop-2,jtop-1,ktop-1) +2.495083e-03*utop(itop-2,jtop-1,ktop) +5.544629e-04*utop(itop-2,jtop-1,ktop+1) -5.940674e-05*utop(itop-2,jtop-1,ktop+2) -3.465393e-04*utop(itop-2,jtop,ktop-2) +2.495083e-03*utop(itop-2,jtop,ktop-1) -1.871312e-02*utop(itop-2,jtop,ktop) -4.158472e-03*utop(itop-2,jtop,ktop+1) +4.455505e-04*utop(itop-2,jtop,ktop+2) -7.700874e-05*utop(itop-2,jtop+1,ktop-2) +5.544629e-04*utop(itop-2,jtop+1,ktop-1) -4.158472e-03*utop(itop-2,jtop+1,ktop) -9.241048e-04*utop(itop-2,jtop+1,ktop+1) +9.901123e-05*utop(itop-2,jtop+1,ktop+2) +8.250936e-06*utop(itop-2,jtop+2,ktop-2) -5.940674e-05*utop(itop-2,jtop+2,ktop-1) +4.455505e-04*utop(itop-2,jtop+2,ktop) +9.901123e-05*utop(itop-2,jtop+2,ktop+1) -1.060835e-05*utop(itop-2,jtop+2,ktop+2) +5.989568e-05*utop(itop-1,jtop-2,ktop-2) -4.312489e-04*utop(itop-1,jtop-2,ktop-1) +3.234367e-03*utop(itop-1,jtop-2,ktop) +7.187482e-04*utop(itop-1,jtop-2,ktop+1) -7.700874e-05*utop(itop-1,jtop-2,ktop+2) -4.312489e-04*utop(itop-1,jtop-1,ktop-2) +3.104992e-03*utop(itop-1,jtop-1,ktop-1) -2.328744e-02*utop(itop-1,jtop-1,ktop) -5.174987e-03*utop(itop-1,jtop-1,ktop+1) +5.544629e-04*utop(itop-1,jtop-1,ktop+2) +3.234367e-03*utop(itop-1,jtop,ktop-2) -2.328744e-02*utop(itop-1,jtop,ktop-1) +1.746558e-01*utop(itop-1,jtop,ktop) +3.881240e-02*utop(itop-1,jtop,ktop+1) -4.158472e-03*utop(itop-1,jtop,ktop+2) +7.187482e-04*utop(itop-1,jtop+1,ktop-2) -5.174987e-03*utop(itop-1,jtop+1,ktop-1) +3.881240e-02*utop(itop-1,jtop+1,ktop) +8.624978e-03*utop(itop-1,jtop+1,ktop+1) -9.241048e-04*utop(itop-1,jtop+1,ktop+2) -7.700874e-05*utop(itop-1,jtop+2,ktop-2) +5.544629e-04*utop(itop-1,jtop+2,ktop-1) -4.158472e-03*utop(itop-1,jtop+2,ktop) -9.241048e-04*utop(itop-1,jtop+2,ktop+1) +9.901123e-05*utop(itop-1,jtop+2,ktop+2) +2.695306e-04*utop(itop,jtop-2,ktop-2) -1.940620e-03*utop(itop,jtop-2,ktop-1) +1.455465e-02*utop(itop,jtop-2,ktop) +3.234367e-03*utop(itop,jtop-2,ktop+1) -3.465393e-04*utop(itop,jtop-2,ktop+2) -1.940620e-03*utop(itop,jtop-1,ktop-2) +1.397246e-02*utop(itop,jtop-1,ktop-1) -1.047935e-01*utop(itop,jtop-1,ktop) -2.328744e-02*utop(itop,jtop-1,ktop+1) +2.495083e-03*utop(itop,jtop-1,ktop+2) +1.455465e-02*utop(itop,jtop,ktop-2) -1.047935e-01*utop(itop,jtop,ktop-1) +7.859512e-01*utop(itop,jtop,ktop) +1.746558e-01*utop(itop,jtop,ktop+1) -1.871312e-02*utop(itop,jtop,ktop+2) +3.234367e-03*utop(itop,jtop+1,ktop-2) -2.328744e-02*utop(itop,jtop+1,ktop-1) +1.746558e-01*utop(itop,jtop+1,ktop) +3.881240e-02*utop(itop,jtop+1,ktop+1) -4.158472e-03*utop(itop,jtop+1,ktop+2) -3.465393e-04*utop(itop,jtop+2,ktop-2) +2.495083e-03*utop(itop,jtop+2,ktop-1) -1.871312e-02*utop(itop,jtop+2,ktop) -4.158472e-03*utop(itop,jtop+2,ktop+1) +4.455505e-04*utop(itop,jtop+2,ktop+2) -3.593741e-05*utop(itop+1,jtop-2,ktop-2) +2.587494e-04*utop(itop+1,jtop-2,ktop-1) -1.940620e-03*utop(itop+1,jtop-2,ktop) -4.312489e-04*utop(itop+1,jtop-2,ktop+1) +4.620524e-05*utop(itop+1,jtop-2,ktop+2) +2.587494e-04*utop(itop+1,jtop-1,ktop-2) -1.862995e-03*utop(itop+1,jtop-1,ktop-1) +1.397246e-02*utop(itop+1,jtop-1,ktop) +3.104992e-03*utop(itop+1,jtop-1,ktop+1) -3.326777e-04*utop(itop+1,jtop-1,ktop+2) -1.940620e-03*utop(itop+1,jtop,ktop-2) +1.397246e-02*utop(itop+1,jtop,ktop-1) -1.047935e-01*utop(itop+1,jtop,ktop) -2.328744e-02*utop(itop+1,jtop,ktop+1) +2.495083e-03*utop(itop+1,jtop,ktop+2) -4.312489e-04*utop(itop+1,jtop+1,ktop-2) +3.104992e-03*utop(itop+1,jtop+1,ktop-1) -2.328744e-02*utop(itop+1,jtop+1,ktop) -5.174987e-03*utop(itop+1,jtop+1,ktop+1) +5.544629e-04*utop(itop+1,jtop+1,ktop+2) +4.620524e-05*utop(itop+1,jtop+2,ktop-2) -3.326777e-04*utop(itop+1,jtop+2,ktop-1) +2.495083e-03*utop(itop+1,jtop+2,ktop) +5.544629e-04*utop(itop+1,jtop+2,ktop+1) -5.940674e-05*utop(itop+1,jtop+2,ktop+2) +4.991307e-06*utop(itop+2,jtop-2,ktop-2) -3.593741e-05*utop(itop+2,jtop-2,ktop-1) +2.695306e-04*utop(itop+2,jtop-2,ktop) +5.989568e-05*utop(itop+2,jtop-2,ktop+1) -6.417395e-06*utop(itop+2,jtop-2,ktop+2) -3.593741e-05*utop(itop+2,jtop-1,ktop-2) +2.587494e-04*utop(itop+2,jtop-1,ktop-1) -1.940620e-03*utop(itop+2,jtop-1,ktop) -4.312489e-04*utop(itop+2,jtop-1,ktop+1) +4.620524e-05*utop(itop+2,jtop-1,ktop+2) +2.695306e-04*utop(itop+2,jtop,ktop-2) -1.940620e-03*utop(itop+2,jtop,ktop-1) +1.455465e-02*utop(itop+2,jtop,ktop) +3.234367e-03*utop(itop+2,jtop,ktop+1) -3.465393e-04*utop(itop+2,jtop,ktop+2) +5.989568e-05*utop(itop+2,jtop+1,ktop-2) -4.312489e-04*utop(itop+2,jtop+1,ktop-1) +3.234367e-03*utop(itop+2,jtop+1,ktop) +7.187482e-04*utop(itop+2,jtop+1,ktop+1) -7.700874e-05*utop(itop+2,jtop+1,ktop+2) -6.417395e-06*utop(itop+2,jtop+2,ktop-2) +4.620524e-05*utop(itop+2,jtop+2,ktop-1) -3.465393e-04*utop(itop+2,jtop+2,ktop) -7.700874e-05*utop(itop+2,jtop+2,ktop+1) +8.250936e-06*utop(itop+2,jtop+2,ktop+2));
	u(i+1,j+0,k+0) = ( +8.250936e-06*utop(itop-2,jtop-2,ktop-2) -7.700874e-05*utop(itop-2,jtop-2,ktop-1) -3.465393e-04*utop(itop-2,jtop-2,ktop) +4.620524e-05*utop(itop-2,jtop-2,ktop+1) -6.417395e-06*utop(itop-2,jtop-2,ktop+2) -7.700874e-05*utop(itop-2,jtop-1,ktop-2) +7.187482e-04*utop(itop-2,jtop-1,ktop-1) +3.234367e-03*utop(itop-2,jtop-1,ktop) -4.312489e-04*utop(itop-2,jtop-1,ktop+1) +5.989568e-05*utop(itop-2,jtop-1,ktop+2) -3.465393e-04*utop(itop-2,jtop,ktop-2) +3.234367e-03*utop(itop-2,jtop,ktop-1) +1.455465e-02*utop(itop-2,jtop,ktop) -1.940620e-03*utop(itop-2,jtop,ktop+1) +2.695306e-04*utop(itop-2,jtop,ktop+2) +4.620524e-05*utop(itop-2,jtop+1,ktop-2) -4.312489e-04*utop(itop-2,jtop+1,ktop-1) -1.940620e-03*utop(itop-2,jtop+1,ktop) +2.587494e-04*utop(itop-2,jtop+1,ktop+1) -3.593741e-05*utop(itop-2,jtop+1,ktop+2) -6.417395e-06*utop(itop-2,jtop+2,ktop-2) +5.989568e-05*utop(itop-2,jtop+2,ktop-1) +2.695306e-04*utop(itop-2,jtop+2,ktop) -3.593741e-05*utop(itop-2,jtop+2,ktop+1) +4.991307e-06*utop(itop-2,jtop+2,ktop+2) -5.940674e-05*utop(itop-1,jtop-2,ktop-2) +5.544629e-04*utop(itop-1,jtop-2,ktop-1) +2.495083e-03*utop(itop-1,jtop-2,ktop) -3.326777e-04*utop(itop-1,jtop-2,ktop+1) +4.620524e-05*utop(itop-1,jtop-2,ktop+2) +5.544629e-04*utop(itop-1,jtop-1,ktop-2) -5.174987e-03*utop(itop-1,jtop-1,ktop-1) -2.328744e-02*utop(itop-1,jtop-1,ktop) +3.104992e-03*utop(itop-1,jtop-1,ktop+1) -4.312489e-04*utop(itop-1,jtop-1,ktop+2) +2.495083e-03*utop(itop-1,jtop,ktop-2) -2.328744e-02*utop(itop-1,jtop,ktop-1) -1.047935e-01*utop(itop-1,jtop,ktop) +1.397246e-02*utop(itop-1,jtop,ktop+1) -1.940620e-03*utop(itop-1,jtop,ktop+2) -3.326777e-04*utop(itop-1,jtop+1,ktop-2) +3.104992e-03*utop(itop-1,jtop+1,ktop-1) +1.397246e-02*utop(itop-1,jtop+1,ktop) -1.862995e-03*utop(itop-1,jtop+1,ktop+1) +2.587494e-04*utop(itop-1,jtop+1,ktop+2) +4.620524e-05*utop(itop-1,jtop+2,ktop-2) -4.312489e-04*utop(itop-1,jtop+2,ktop-1) -1.940620e-03*utop(itop-1,jtop+2,ktop) +2.587494e-04*utop(itop-1,jtop+2,ktop+1) -3.593741e-05*utop(itop-1,jtop+2,ktop+2) +4.455505e-04*utop(itop,jtop-2,ktop-2) -4.158472e-03*utop(itop,jtop-2,ktop-1) -1.871312e-02*utop(itop,jtop-2,ktop) +2.495083e-03*utop(itop,jtop-2,ktop+1) -3.465393e-04*utop(itop,jtop-2,ktop+2) -4.158472e-03*utop(itop,jtop-1,ktop-2) +3.881240e-02*utop(itop,jtop-1,ktop-1) +1.746558e-01*utop(itop,jtop-1,ktop) -2.328744e-02*utop(itop,jtop-1,ktop+1) +3.234367e-03*utop(itop,jtop-1,ktop+2) -1.871312e-02*utop(itop,jtop,ktop-2) +1.746558e-01*utop(itop,jtop,ktop-1) +7.859512e-01*utop(itop,jtop,ktop) -1.047935e-01*utop(itop,jtop,ktop+1) +1.455465e-02*utop(itop,jtop,ktop+2) +2.495083e-03*utop(itop,jtop+1,ktop-2) -2.328744e-02*utop(itop,jtop+1,ktop-1) -1.047935e-01*utop(itop,jtop+1,ktop) +1.397246e-02*utop(itop,jtop+1,ktop+1) -1.940620e-03*utop(itop,jtop+1,ktop+2) -3.465393e-04*utop(itop,jtop+2,ktop-2) +3.234367e-03*utop(itop,jtop+2,ktop-1) +1.455465e-02*utop(itop,jtop+2,ktop) -1.940620e-03*utop(itop,jtop+2,ktop+1) +2.695306e-04*utop(itop,jtop+2,ktop+2) +9.901123e-05*utop(itop+1,jtop-2,ktop-2) -9.241048e-04*utop(itop+1,jtop-2,ktop-1) -4.158472e-03*utop(itop+1,jtop-2,ktop) +5.544629e-04*utop(itop+1,jtop-2,ktop+1) -7.700874e-05*utop(itop+1,jtop-2,ktop+2) -9.241048e-04*utop(itop+1,jtop-1,ktop-2) +8.624978e-03*utop(itop+1,jtop-1,ktop-1) +3.881240e-02*utop(itop+1,jtop-1,ktop) -5.174987e-03*utop(itop+1,jtop-1,ktop+1) +7.187482e-04*utop(itop+1,jtop-1,ktop+2) -4.158472e-03*utop(itop+1,jtop,ktop-2) +3.881240e-02*utop(itop+1,jtop,ktop-1) +1.746558e-01*utop(itop+1,jtop,ktop) -2.328744e-02*utop(itop+1,jtop,ktop+1) +3.234367e-03*utop(itop+1,jtop,ktop+2) +5.544629e-04*utop(itop+1,jtop+1,ktop-2) -5.174987e-03*utop(itop+1,jtop+1,ktop-1) -2.328744e-02*utop(itop+1,jtop+1,ktop) +3.104992e-03*utop(itop+1,jtop+1,ktop+1) -4.312489e-04*utop(itop+1,jtop+1,ktop+2) -7.700874e-05*utop(itop+1,jtop+2,ktop-2) +7.187482e-04*utop(itop+1,jtop+2,ktop-1) +3.234367e-03*utop(itop+1,jtop+2,ktop) -4.312489e-04*utop(itop+1,jtop+2,ktop+1) +5.989568e-05*utop(itop+1,jtop+2,ktop+2) -1.060835e-05*utop(itop+2,jtop-2,ktop-2) +9.901123e-05*utop(itop+2,jtop-2,ktop-1) +4.455505e-04*utop(itop+2,jtop-2,ktop) -5.940674e-05*utop(itop+2,jtop-2,ktop+1) +8.250936e-06*utop(itop+2,jtop-2,ktop+2) +9.901123e-05*utop(itop+2,jtop-1,ktop-2) -9.241048e-04*utop(itop+2,jtop-1,ktop-1) -4.158472e-03*utop(itop+2,jtop-1,ktop) +5.544629e-04*utop(itop+2,jtop-1,ktop+1) -7.700874e-05*utop(itop+2,jtop-1,ktop+2) +4.455505e-04*utop(itop+2,jtop,ktop-2) -4.158472e-03*utop(itop+2,jtop,ktop-1) -1.871312e-02*utop(itop+2,jtop,ktop) +2.495083e-03*utop(itop+2,jtop,ktop+1) -3.465393e-04*utop(itop+2,jtop,ktop+2) -5.940674e-05*utop(itop+2,jtop+1,ktop-2) +5.544629e-04*utop(itop+2,jtop+1,ktop-1) +2.495083e-03*utop(itop+2,jtop+1,ktop) -3.326777e-04*utop(itop+2,jtop+1,ktop+1) +4.620524e-05*utop(itop+2,jtop+1,ktop+2) +8.250936e-06*utop(itop+2,jtop+2,ktop-2) -7.700874e-05*utop(itop+2,jtop+2,ktop-1) -3.465393e-04*utop(itop+2,jtop+2,ktop) +4.620524e-05*utop(itop+2,jtop+2,ktop+1) -6.417395e-06*utop(itop+2,jtop+2,ktop+2));
	u(i+1,j+0,k+1) = ( -6.417395e-06*utop(itop-2,jtop-2,ktop-2) +4.620524e-05*utop(itop-2,jtop-2,ktop-1) -3.465393e-04*utop(itop-2,jtop-2,ktop) -7.700874e-05*utop(itop-2,jtop-2,ktop+1) +8.250936e-06*utop(itop-2,jtop-2,ktop+2) +5.989568e-05*utop(itop-2,jtop-1,ktop-2) -4.312489e-04*utop(itop-2,jtop-1,ktop-1) +3.234367e-03*utop(itop-2,jtop-1,ktop) +7.187482e-04*utop(itop-2,jtop-1,ktop+1) -7.700874e-05*utop(itop-2,jtop-1,ktop+2) +2.695306e-04*utop(itop-2,jtop,ktop-2) -1.940620e-03*utop(itop-2,jtop,ktop-1) +1.455465e-02*utop(itop-2,jtop,ktop) +3.234367e-03*utop(itop-2,jtop,ktop+1) -3.465393e-04*utop(itop-2,jtop,ktop+2) -3.593741e-05*utop(itop-2,jtop+1,ktop-2) +2.587494e-04*utop(itop-2,jtop+1,ktop-1) -1.940620e-03*utop(itop-2,jtop+1,ktop) -4.312489e-04*utop(itop-2,jtop+1,ktop+1) +4.620524e-05*utop(itop-2,jtop+1,ktop+2) +4.991307e-06*utop(itop-2,jtop+2,ktop-2) -3.593741e-05*utop(itop-2,jtop+2,ktop-1) +2.695306e-04*utop(itop-2,jtop+2,ktop) +5.989568e-05*utop(itop-2,jtop+2,ktop+1) -6.417395e-06*utop(itop-2,jtop+2,ktop+2) +4.620524e-05*utop(itop-1,jtop-2,ktop-2) -3.326777e-04*utop(itop-1,jtop-2,ktop-1) +2.495083e-03*utop(itop-1,jtop-2,ktop) +5.544629e-04*utop(itop-1,jtop-2,ktop+1) -5.940674e-05*utop(itop-1,jtop-2,ktop+2) -4.312489e-04*utop(itop-1,jtop-1,ktop-2) +3.104992e-03*utop(itop-1,jtop-1,ktop-1) -2.328744e-02*utop(itop-1,jtop-1,ktop) -5.174987e-03*utop(itop-1,jtop-1,ktop+1) +5.544629e-04*utop(itop-1,jtop-1,ktop+2) -1.940620e-03*utop(itop-1,jtop,ktop-2) +1.397246e-02*utop(itop-1,jtop,ktop-1) -1.047935e-01*utop(itop-1,jtop,ktop) -2.328744e-02*utop(itop-1,jtop,ktop+1) +2.495083e-03*utop(itop-1,jtop,ktop+2) +2.587494e-04*utop(itop-1,jtop+1,ktop-2) -1.862995e-03*utop(itop-1,jtop+1,ktop-1) +1.397246e-02*utop(itop-1,jtop+1,ktop) +3.104992e-03*utop(itop-1,jtop+1,ktop+1) -3.326777e-04*utop(itop-1,jtop+1,ktop+2) -3.593741e-05*utop(itop-1,jtop+2,ktop-2) +2.587494e-04*utop(itop-1,jtop+2,ktop-1) -1.940620e-03*utop(itop-1,jtop+2,ktop) -4.312489e-04*utop(itop-1,jtop+2,ktop+1) +4.620524e-05*utop(itop-1,jtop+2,ktop+2) -3.465393e-04*utop(itop,jtop-2,ktop-2) +2.495083e-03*utop(itop,jtop-2,ktop-1) -1.871312e-02*utop(itop,jtop-2,ktop) -4.158472e-03*utop(itop,jtop-2,ktop+1) +4.455505e-04*utop(itop,jtop-2,ktop+2) +3.234367e-03*utop(itop,jtop-1,ktop-2) -2.328744e-02*utop(itop,jtop-1,ktop-1) +1.746558e-01*utop(itop,jtop-1,ktop) +3.881240e-02*utop(itop,jtop-1,ktop+1) -4.158472e-03*utop(itop,jtop-1,ktop+2) +1.455465e-02*utop(itop,jtop,ktop-2) -1.047935e-01*utop(itop,jtop,ktop-1) +7.859512e-01*utop(itop,jtop,ktop) +1.746558e-01*utop(itop,jtop,ktop+1) -1.871312e-02*utop(itop,jtop,ktop+2) -1.940620e-03*utop(itop,jtop+1,ktop-2) +1.397246e-02*utop(itop,jtop+1,ktop-1) -1.047935e-01*utop(itop,jtop+1,ktop) -2.328744e-02*utop(itop,jtop+1,ktop+1) +2.495083e-03*utop(itop,jtop+1,ktop+2) +2.695306e-04*utop(itop,jtop+2,ktop-2) -1.940620e-03*utop(itop,jtop+2,ktop-1) +1.455465e-02*utop(itop,jtop+2,ktop) +3.234367e-03*utop(itop,jtop+2,ktop+1) -3.465393e-04*utop(itop,jtop+2,ktop+2) -7.700874e-05*utop(itop+1,jtop-2,ktop-2) +5.544629e-04*utop(itop+1,jtop-2,ktop-1) -4.158472e-03*utop(itop+1,jtop-2,ktop) -9.241048e-04*utop(itop+1,jtop-2,ktop+1) +9.901123e-05*utop(itop+1,jtop-2,ktop+2) +7.187482e-04*utop(itop+1,jtop-1,ktop-2) -5.174987e-03*utop(itop+1,jtop-1,ktop-1) +3.881240e-02*utop(itop+1,jtop-1,ktop) +8.624978e-03*utop(itop+1,jtop-1,ktop+1) -9.241048e-04*utop(itop+1,jtop-1,ktop+2) +3.234367e-03*utop(itop+1,jtop,ktop-2) -2.328744e-02*utop(itop+1,jtop,ktop-1) +1.746558e-01*utop(itop+1,jtop,ktop) +3.881240e-02*utop(itop+1,jtop,ktop+1) -4.158472e-03*utop(itop+1,jtop,ktop+2) -4.312489e-04*utop(itop+1,jtop+1,ktop-2) +3.104992e-03*utop(itop+1,jtop+1,ktop-1) -2.328744e-02*utop(itop+1,jtop+1,ktop) -5.174987e-03*utop(itop+1,jtop+1,ktop+1) +5.544629e-04*utop(itop+1,jtop+1,ktop+2) +5.989568e-05*utop(itop+1,jtop+2,ktop-2) -4.312489e-04*utop(itop+1,jtop+2,ktop-1) +3.234367e-03*utop(itop+1,jtop+2,ktop) +7.187482e-04*utop(itop+1,jtop+2,ktop+1) -7.700874e-05*utop(itop+1,jtop+2,ktop+2) +8.250936e-06*utop(itop+2,jtop-2,ktop-2) -5.940674e-05*utop(itop+2,jtop-2,ktop-1) +4.455505e-04*utop(itop+2,jtop-2,ktop) +9.901123e-05*utop(itop+2,jtop-2,ktop+1) -1.060835e-05*utop(itop+2,jtop-2,ktop+2) -7.700874e-05*utop(itop+2,jtop-1,ktop-2) +5.544629e-04*utop(itop+2,jtop-1,ktop-1) -4.158472e-03*utop(itop+2,jtop-1,ktop) -9.241048e-04*utop(itop+2,jtop-1,ktop+1) +9.901123e-05*utop(itop+2,jtop-1,ktop+2) -3.465393e-04*utop(itop+2,jtop,ktop-2) +2.495083e-03*utop(itop+2,jtop,ktop-1) -1.871312e-02*utop(itop+2,jtop,ktop) -4.158472e-03*utop(itop+2,jtop,ktop+1) +4.455505e-04*utop(itop+2,jtop,ktop+2) +4.620524e-05*utop(itop+2,jtop+1,ktop-2) -3.326777e-04*utop(itop+2,jtop+1,ktop-1) +2.495083e-03*utop(itop+2,jtop+1,ktop) +5.544629e-04*utop(itop+2,jtop+1,ktop+1) -5.940674e-05*utop(itop+2,jtop+1,ktop+2) -6.417395e-06*utop(itop+2,jtop+2,ktop-2) +4.620524e-05*utop(itop+2,jtop+2,ktop-1) -3.465393e-04*utop(itop+2,jtop+2,ktop) -7.700874e-05*utop(itop+2,jtop+2,ktop+1) +8.250936e-06*utop(itop+2,jtop+2,ktop+2));
	u(i+1,j+1,k+0) = ( -6.417395e-06*utop(itop-2,jtop-2,ktop-2) +5.989568e-05*utop(itop-2,jtop-2,ktop-1) +2.695306e-04*utop(itop-2,jtop-2,ktop) -3.593741e-05*utop(itop-2,jtop-2,ktop+1) +4.991307e-06*utop(itop-2,jtop-2,ktop+2) +4.620524e-05*utop(itop-2,jtop-1,ktop-2) -4.312489e-04*utop(itop-2,jtop-1,ktop-1) -1.940620e-03*utop(itop-2,jtop-1,ktop) +2.587494e-04*utop(itop-2,jtop-1,ktop+1) -3.593741e-05*utop(itop-2,jtop-1,ktop+2) -3.465393e-04*utop(itop-2,jtop,ktop-2) +3.234367e-03*utop(itop-2,jtop,ktop-1) +1.455465e-02*utop(itop-2,jtop,ktop) -1.940620e-03*utop(itop-2,jtop,ktop+1) +2.695306e-04*utop(itop-2,jtop,ktop+2) -7.700874e-05*utop(itop-2,jtop+1,ktop-2) +7.187482e-04*utop(itop-2,jtop+1,ktop-1) +3.234367e-03*utop(itop-2,jtop+1,ktop) -4.312489e-04*utop(itop-2,jtop+1,ktop+1) +5.989568e-05*utop(itop-2,jtop+1,ktop+2) +8.250936e-06*utop(itop-2,jtop+2,ktop-2) -7.700874e-05*utop(itop-2,jtop+2,ktop-1) -3.465393e-04*utop(itop-2,jtop+2,ktop) +4.620524e-05*utop(itop-2,jtop+2,ktop+1) -6.417395e-06*utop(itop-2,jtop+2,ktop+2) +4.620524e-05*utop(itop-1,jtop-2,ktop-2) -4.312489e-04*utop(itop-1,jtop-2,ktop-1) -1.940620e-03*utop(itop-1,jtop-2,ktop) +2.587494e-04*utop(itop-1,jtop-2,ktop+1) -3.593741e-05*utop(itop-1,jtop-2,ktop+2) -3.326777e-04*utop(itop-1,jtop-1,ktop-2) +3.104992e-03*utop(itop-1,jtop-1,ktop-1) +1.397246e-02*utop(itop-1,jtop-1,ktop) -1.862995e-03*utop(itop-1,jtop-1,ktop+1) +2.587494e-04*utop(itop-1,jtop-1,ktop+2) +2.495083e-03*utop(itop-1,jtop,ktop-2) -2.328744e-02*utop(itop-1,jtop,ktop-1) -1.047935e-01*utop(itop-1,jtop,ktop) +1.397246e-02*utop(itop-1,jtop,ktop+1) -1.940620e-03*utop(itop-1,jtop,ktop+2) +5.544629e-04*utop(itop-1,jtop+1,ktop-2) -5.174987e-03*utop(itop-1,jtop+1,ktop-1) -2.328744e-02*utop(itop-1,jtop+1,ktop) +3.104992e-03*utop(itop-1,jtop+1,ktop+1) -4.312489e-04*utop(itop-1,jtop+1,ktop+2) -5.940674e-05*utop(itop-1,jtop+2,ktop-2) +5.544629e-04*utop(itop-1,jtop+2,ktop-1) +2.495083e-03*utop(itop-1,jtop+2,ktop) -3.326777e-04*utop(itop-1,jtop+2,ktop+1) +4.620524e-05*utop(itop-1,jtop+2,ktop+2) -3.465393e-04*utop(itop,jtop-2,ktop-2) +3.234367e-03*utop(itop,jtop-2,ktop-1) +1.455465e-02*utop(itop,jtop-2,ktop) -1.940620e-03*utop(itop,jtop-2,ktop+1) +2.695306e-04*utop(itop,jtop-2,ktop+2) +2.495083e-03*utop(itop,jtop-1,ktop-2) -2.328744e-02*utop(itop,jtop-1,ktop-1) -1.047935e-01*utop(itop,jtop-1,ktop) +1.397246e-02*utop(itop,jtop-1,ktop+1) -1.940620e-03*utop(itop,jtop-1,ktop+2) -1.871312e-02*utop(itop,jtop,ktop-2) +1.746558e-01*utop(itop,jtop,ktop-1) +7.859512e-01*utop(itop,jtop,ktop) -1.047935e-01*utop(itop,jtop,ktop+1) +1.455465e-02*utop(itop,jtop,ktop+2) -4.158472e-03*utop(itop,jtop+1,ktop-2) +3.881240e-02*utop(itop,jtop+1,ktop-1) +1.746558e-01*utop(itop,jtop+1,ktop) -2.328744e-02*utop(itop,jtop+1,ktop+1) +3.234367e-03*utop(itop,jtop+1,ktop+2) +4.455505e-04*utop(itop,jtop+2,ktop-2) -4.158472e-03*utop(itop,jtop+2,ktop-1) -1.871312e-02*utop(itop,jtop+2,ktop) +2.495083e-03*utop(itop,jtop+2,ktop+1) -3.465393e-04*utop(itop,jtop+2,ktop+2) -7.700874e-05*utop(itop+1,jtop-2,ktop-2) +7.187482e-04*utop(itop+1,jtop-2,ktop-1) +3.234367e-03*utop(itop+1,jtop-2,ktop) -4.312489e-04*utop(itop+1,jtop-2,ktop+1) +5.989568e-05*utop(itop+1,jtop-2,ktop+2) +5.544629e-04*utop(itop+1,jtop-1,ktop-2) -5.174987e-03*utop(itop+1,jtop-1,ktop-1) -2.328744e-02*utop(itop+1,jtop-1,ktop) +3.104992e-03*utop(itop+1,jtop-1,ktop+1) -4.312489e-04*utop(itop+1,jtop-1,ktop+2) -4.158472e-03*utop(itop+1,jtop,ktop-2) +3.881240e-02*utop(itop+1,jtop,ktop-1) +1.746558e-01*utop(itop+1,jtop,ktop) -2.328744e-02*utop(itop+1,jtop,ktop+1) +3.234367e-03*utop(itop+1,jtop,ktop+2) -9.241048e-04*utop(itop+1,jtop+1,ktop-2) +8.624978e-03*utop(itop+1,jtop+1,ktop-1) +3.881240e-02*utop(itop+1,jtop+1,ktop) -5.174987e-03*utop(itop+1,jtop+1,ktop+1) +7.187482e-04*utop(itop+1,jtop+1,ktop+2) +9.901123e-05*utop(itop+1,jtop+2,ktop-2) -9.241048e-04*utop(itop+1,jtop+2,ktop-1) -4.158472e-03*utop(itop+1,jtop+2,ktop) +5.544629e-04*utop(itop+1,jtop+2,ktop+1) -7.700874e-05*utop(itop+1,jtop+2,ktop+2) +8.250936e-06*utop(itop+2,jtop-2,ktop-2) -7.700874e-05*utop(itop+2,jtop-2,ktop-1) -3.465393e-04*utop(itop+2,jtop-2,ktop) +4.620524e-05*utop(itop+2,jtop-2,ktop+1) -6.417395e-06*utop(itop+2,jtop-2,ktop+2) -5.940674e-05*utop(itop+2,jtop-1,ktop-2) +5.544629e-04*utop(itop+2,jtop-1,ktop-1) +2.495083e-03*utop(itop+2,jtop-1,ktop) -3.326777e-04*utop(itop+2,jtop-1,ktop+1) +4.620524e-05*utop(itop+2,jtop-1,ktop+2) +4.455505e-04*utop(itop+2,jtop,ktop-2) -4.158472e-03*utop(itop+2,jtop,ktop-1) -1.871312e-02*utop(itop+2,jtop,ktop) +2.495083e-03*utop(itop+2,jtop,ktop+1) -3.465393e-04*utop(itop+2,jtop,ktop+2) +9.901123e-05*utop(itop+2,jtop+1,ktop-2) -9.241048e-04*utop(itop+2,jtop+1,ktop-1) -4.158472e-03*utop(itop+2,jtop+1,ktop) +5.544629e-04*utop(itop+2,jtop+1,ktop+1) -7.700874e-05*utop(itop+2,jtop+1,ktop+2) -1.060835e-05*utop(itop+2,jtop+2,ktop-2) +9.901123e-05*utop(itop+2,jtop+2,ktop-1) +4.455505e-04*utop(itop+2,jtop+2,ktop) -5.940674e-05*utop(itop+2,jtop+2,ktop+1) +8.250936e-06*utop(itop+2,jtop+2,ktop+2));
	u(i+1,j+1,k+1) = ( +4.991307e-06*utop(itop-2,jtop-2,ktop-2) -3.593741e-05*utop(itop-2,jtop-2,ktop-1) +2.695306e-04*utop(itop-2,jtop-2,ktop) +5.989568e-05*utop(itop-2,jtop-2,ktop+1) -6.417395e-06*utop(itop-2,jtop-2,ktop+2) -3.593741e-05*utop(itop-2,jtop-1,ktop-2) +2.587494e-04*utop(itop-2,jtop-1,ktop-1) -1.940620e-03*utop(itop-2,jtop-1,ktop) -4.312489e-04*utop(itop-2,jtop-1,ktop+1) +4.620524e-05*utop(itop-2,jtop-1,ktop+2) +2.695306e-04*utop(itop-2,jtop,ktop-2) -1.940620e-03*utop(itop-2,jtop,ktop-1) +1.455465e-02*utop(itop-2,jtop,ktop) +3.234367e-03*utop(itop-2,jtop,ktop+1) -3.465393e-04*utop(itop-2,jtop,ktop+2) +5.989568e-05*utop(itop-2,jtop+1,ktop-2) -4.312489e-04*utop(itop-2,jtop+1,ktop-1) +3.234367e-03*utop(itop-2,jtop+1,ktop) +7.187482e-04*utop(itop-2,jtop+1,ktop+1) -7.700874e-05*utop(itop-2,jtop+1,ktop+2) -6.417395e-06*utop(itop-2,jtop+2,ktop-2) +4.620524e-05*utop(itop-2,jtop+2,ktop-1) -3.465393e-04*utop(itop-2,jtop+2,ktop) -7.700874e-05*utop(itop-2,jtop+2,ktop+1) +8.250936e-06*utop(itop-2,jtop+2,ktop+2) -3.593741e-05*utop(itop-1,jtop-2,ktop-2) +2.587494e-04*utop(itop-1,jtop-2,ktop-1) -1.940620e-03*utop(itop-1,jtop-2,ktop) -4.312489e-04*utop(itop-1,jtop-2,ktop+1) +4.620524e-05*utop(itop-1,jtop-2,ktop+2) +2.587494e-04*utop(itop-1,jtop-1,ktop-2) -1.862995e-03*utop(itop-1,jtop-1,ktop-1) +1.397246e-02*utop(itop-1,jtop-1,ktop) +3.104992e-03*utop(itop-1,jtop-1,ktop+1) -3.326777e-04*utop(itop-1,jtop-1,ktop+2) -1.940620e-03*utop(itop-1,jtop,ktop-2) +1.397246e-02*utop(itop-1,jtop,ktop-1) -1.047935e-01*utop(itop-1,jtop,ktop) -2.328744e-02*utop(itop-1,jtop,ktop+1) +2.495083e-03*utop(itop-1,jtop,ktop+2) -4.312489e-04*utop(itop-1,jtop+1,ktop-2) +3.104992e-03*utop(itop-1,jtop+1,ktop-1) -2.328744e-02*utop(itop-1,jtop+1,ktop) -5.174987e-03*utop(itop-1,jtop+1,ktop+1) +5.544629e-04*utop(itop-1,jtop+1,ktop+2) +4.620524e-05*utop(itop-1,jtop+2,ktop-2) -3.326777e-04*utop(itop-1,jtop+2,ktop-1) +2.495083e-03*utop(itop-1,jtop+2,ktop) +5.544629e-04*utop(itop-1,jtop+2,ktop+1) -5.940674e-05*utop(itop-1,jtop+2,ktop+2) +2.695306e-04*utop(itop,jtop-2,ktop-2) -1.940620e-03*utop(itop,jtop-2,ktop-1) +1.455465e-02*utop(itop,jtop-2,ktop) +3.234367e-03*utop(itop,jtop-2,ktop+1) -3.465393e-04*utop(itop,jtop-2,ktop+2) -1.940620e-03*utop(itop,jtop-1,ktop-2) +1.397246e-02*utop(itop,jtop-1,ktop-1) -1.047935e-01*utop(itop,jtop-1,ktop) -2.328744e-02*utop(itop,jtop-1,ktop+1) +2.495083e-03*utop(itop,jtop-1,ktop+2) +1.455465e-02*utop(itop,jtop,ktop-2) -1.047935e-01*utop(itop,jtop,ktop-1) +7.859512e-01*utop(itop,jtop,ktop) +1.746558e-01*utop(itop,jtop,ktop+1) -1.871312e-02*utop(itop,jtop,ktop+2) +3.234367e-03*utop(itop,jtop+1,ktop-2) -2.328744e-02*utop(itop,jtop+1,ktop-1) +1.746558e-01*utop(itop,jtop+1,ktop) +3.881240e-02*utop(itop,jtop+1,ktop+1) -4.158472e-03*utop(itop,jtop+1,ktop+2) -3.465393e-04*utop(itop,jtop+2,ktop-2) +2.495083e-03*utop(itop,jtop+2,ktop-1) -1.871312e-02*utop(itop,jtop+2,ktop) -4.158472e-03*utop(itop,jtop+2,ktop+1) +4.455505e-04*utop(itop,jtop+2,ktop+2) +5.989568e-05*utop(itop+1,jtop-2,ktop-2) -4.312489e-04*utop(itop+1,jtop-2,ktop-1) +3.234367e-03*utop(itop+1,jtop-2,ktop) +7.187482e-04*utop(itop+1,jtop-2,ktop+1) -7.700874e-05*utop(itop+1,jtop-2,ktop+2) -4.312489e-04*utop(itop+1,jtop-1,ktop-2) +3.104992e-03*utop(itop+1,jtop-1,ktop-1) -2.328744e-02*utop(itop+1,jtop-1,ktop) -5.174987e-03*utop(itop+1,jtop-1,ktop+1) +5.544629e-04*utop(itop+1,jtop-1,ktop+2) +3.234367e-03*utop(itop+1,jtop,ktop-2) -2.328744e-02*utop(itop+1,jtop,ktop-1) +1.746558e-01*utop(itop+1,jtop,ktop) +3.881240e-02*utop(itop+1,jtop,ktop+1) -4.158472e-03*utop(itop+1,jtop,ktop+2) +7.187482e-04*utop(itop+1,jtop+1,ktop-2) -5.174987e-03*utop(itop+1,jtop+1,ktop-1) +3.881240e-02*utop(itop+1,jtop+1,ktop) +8.624978e-03*utop(itop+1,jtop+1,ktop+1) -9.241048e-04*utop(itop+1,jtop+1,ktop+2) -7.700874e-05*utop(itop+1,jtop+2,ktop-2) +5.544629e-04*utop(itop+1,jtop+2,ktop-1) -4.158472e-03*utop(itop+1,jtop+2,ktop) -9.241048e-04*utop(itop+1,jtop+2,ktop+1) +9.901123e-05*utop(itop+1,jtop+2,ktop+2) -6.417395e-06*utop(itop+2,jtop-2,ktop-2) +4.620524e-05*utop(itop+2,jtop-2,ktop-1) -3.465393e-04*utop(itop+2,jtop-2,ktop) -7.700874e-05*utop(itop+2,jtop-2,ktop+1) +8.250936e-06*utop(itop+2,jtop-2,ktop+2) +4.620524e-05*utop(itop+2,jtop-1,ktop-2) -3.326777e-04*utop(itop+2,jtop-1,ktop-1) +2.495083e-03*utop(itop+2,jtop-1,ktop) +5.544629e-04*utop(itop+2,jtop-1,ktop+1) -5.940674e-05*utop(itop+2,jtop-1,ktop+2) -3.465393e-04*utop(itop+2,jtop,ktop-2) +2.495083e-03*utop(itop+2,jtop,ktop-1) -1.871312e-02*utop(itop+2,jtop,ktop) -4.158472e-03*utop(itop+2,jtop,ktop+1) +4.455505e-04*utop(itop+2,jtop,ktop+2) -7.700874e-05*utop(itop+2,jtop+1,ktop-2) +5.544629e-04*utop(itop+2,jtop+1,ktop-1) -4.158472e-03*utop(itop+2,jtop+1,ktop) -9.241048e-04*utop(itop+2,jtop+1,ktop+1) +9.901123e-05*utop(itop+2,jtop+1,ktop+2) +8.250936e-06*utop(itop+2,jtop+2,ktop-2) -5.940674e-05*utop(itop+2,jtop+2,ktop-1) +4.455505e-04*utop(itop+2,jtop+2,ktop) +9.901123e-05*utop(itop+2,jtop+2,ktop+1) -1.060835e-05*utop(itop+2,jtop+2,ktop+2));
	

}


template< class S, class O, typename T >
void solver<S,O,T>::interp_coarse_fine_cubic( unsigned ilevel, MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine, bool bcf=false )
{
	
	MeshvarBnd<T> *u    = &fine;
	MeshvarBnd<T> *utop = &coarse;
	
	
	bcf = true;
	
	int
	xoff = u->offset(0),
	yoff = u->offset(1),
	zoff = u->offset(2);
	
	//... don't do anything if we are not an additional refinement region
	if( ilevel <= m_ilevelmin )
		return;
	
	int
	nx = u->size(0), 
	ny = u->size(1), 
	nz = u->size(2);
	
	for( int j=0; j<ny; ++j )
		for( int k=0; k<nz; ++k )
		{
			int jtop = (int)(0.5*(double)(j))+yoff;
			int ktop = (int)(0.5*(double)(k))+zoff;
			
			interp_cubic( coarse, fine, -2, j, k, xoff-1,    jtop, ktop );
			interp_cubic( coarse, fine, nz, j, k, xoff+nz/2, jtop, ktop );
			
		}
	
	for( int i=0; i<nx; ++i )
		for( int k=0; k<nz; ++k )
		{
			int itop = (int)(0.5*(double)(i))+xoff;
			int ktop = (int)(0.5*(double)(k))+zoff;
			
			interp_cubic( coarse, fine, i, -2, k, itop, yoff-1,    ktop );
			interp_cubic( coarse, fine, i, ny, k, itop, yoff+ny/2, ktop );
			
		}
	
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )
		{
			int itop = (int)(0.5*(double)(i))+xoff;
			int jtop = (int)(0.5*(double)(j))+yoff;
			
			interp_cubic( coarse, fine, i, j, -2, itop, jtop, zoff-1 );
			interp_cubic( coarse, fine, i, j, nz, itop, jtop, zoff+nz/2 );
			
		}
}


template< class S, class O, typename T >
void solver<S,O,T>::interp_coarse_fine( unsigned ilevel, MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine, bool bcf )
{
	MeshvarBnd<T> *u    = &fine;
	MeshvarBnd<T> *utop = &coarse;
	
	
	bcf = true;;
	//bcf = false;
	
	int
		xoff = u->offset(0),
		yoff = u->offset(1),
		zoff = u->offset(2);

	//... don't do anything if we are not an additional refinement region
	if( xoff == 0 && yoff == 0 && zoff == 0 )
		return;
	
	int
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
	
	//... set boundary condition for fine grid
	
	#pragma omp parallel for schedule(dynamic)
	for( int ix=-1; ix<=nx; ++ix )
		for( int iy=-1; iy<=ny; ++iy )
			for( int iz=-1; iz<=nz; ++iz )
			{
				bool xbnd=(ix==-1||ix==nx),ybnd=(iy==-1||iy==ny),zbnd=(iz==-1||iz==nz);
				
				//if(ix==-1||ix==nx||iy==-1||iy==ny||iz==-1||iz==nz)
				if( xbnd || ybnd || zbnd )
				//if( xbnd ^ ybnd ^ zbnd )
				{
					
					//... only deal with proper ghostzones
					if( (xbnd&&ybnd) || (xbnd&&zbnd) || (ybnd&&zbnd) || (xbnd&&ybnd&&zbnd))
						continue;
					
					/*int ixtop = (int)(0.5*(double)(ix+2*xoff)+1e-3);
					int iytop = (int)(0.5*(double)(iy+2*yoff)+1e-3);
					int iztop = (int)(0.5*(double)(iz+2*zoff)+1e-3);*/
					
					int ixtop = (int)(0.5*(double)(ix))+xoff;
					int iytop = (int)(0.5*(double)(iy))+yoff;
					int iztop = (int)(0.5*(double)(iz))+zoff;
					
					if( ix==-1 ) ixtop=xoff-1;
					if( iy==-1 ) iytop=yoff-1;
					if( iz==-1 ) iztop=zoff-1;
					
					double ustar1, ustar2, ustar3, uhat;			
					double fac = 0.5;//0.25;
					double flux;;
				    if( ix == -1 && iy%2==0 && iz%2==0 )
					{
						flux = 0.0;
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
							{		
								ustar1 = interp2( (*utop)(ixtop,iytop-1,iztop-1),(*utop)(ixtop,iytop,iztop-1),(*utop)(ixtop,iytop+1,iztop-1), fac*((double)j-0.5) );
								ustar2 = interp2( (*utop)(ixtop,iytop-1,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop,iytop+1,iztop), fac*((double)j-0.5) );
								ustar3 = interp2( (*utop)(ixtop,iytop-1,iztop+1),(*utop)(ixtop,iytop,iztop+1),(*utop)(ixtop,iytop+1,iztop+1), fac*((double)j-0.5) );
									
								uhat   = interp2( /*-1.0, 0.0, 1.0, */ustar1, ustar2, ustar3, fac*((double)k-0.5) );
							
								//(*u)(ix,iy+j,iz+k) = 0.0;//(*utop)(ixtop,iytop,iztop);//interp2( -1.5, 0.0, 1.0, uhat, (*u)(ix+1,iy+j,iz+k), (*u)(ix+2,iy+j,iz+k), -1.0 );
								
								(*u)(ix,iy+j,iz+k) = interp2left( uhat, (*u)(ix+1,iy+j,iz+k), (*u)(ix+2,iy+j,iz+k) );
								
								flux += ((*u)(ix+1,iy+j,iz+k)-(*u)(ix,iy+j,iz+k));
							}
						
						
						
						flux /= 4.0;
						
						double dflux = ((*utop)(ixtop+1,iytop,iztop)-(*utop)(ixtop,iytop,iztop))/2.0 - flux;
						
						//dflux *= 2.0;
						
						if( bcf )
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
									(*u)(ix,iy+j,iz+k) -= dflux;
						else
							(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop+1,iytop,iztop) - 2.0*flux;
						
						
					}
					// right boundary
					if( ix == nx && iy%2==0 && iz%2==0 )
					{
						flux = 0.0;
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
							{		
								ustar1 = interp2( (*utop)(ixtop,iytop-1,iztop-1),(*utop)(ixtop,iytop,iztop-1),(*utop)(ixtop,iytop+1,iztop-1), fac*((double)j-0.5) );
								ustar2 = interp2( (*utop)(ixtop,iytop-1,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop,iytop+1,iztop), fac*((double)j-0.5) );
								ustar3 = interp2( (*utop)(ixtop,iytop-1,iztop+1),(*utop)(ixtop,iytop,iztop+1),(*utop)(ixtop,iytop+1,iztop+1), fac*((double)j-0.5) );
								
								uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
								
								//(*u)(ix,iy+j,iz+k) = 0.0;(*utop)(ixtop,iytop,iztop);//interp2( 1.5, 0.0, -1.0, uhat, (*u)(ix-1,iy+j,iz+k), (*u)(ix-2,iy+j,iz+k), 1.0 );
								(*u)(ix,iy+j,iz+k) = interp2right( (*u)(ix-2,iy+j,iz+k), (*u)(ix-1,iy+j,iz+k), uhat );
								flux += ((*u)(ix,iy+j,iz+k)-(*u)(ix-1,iy+j,iz+k));
							}
						flux /= 4.0;
						
						
						double dflux = ((*utop)(ixtop,iytop,iztop)-(*utop)(ixtop-1,iytop,iztop))/2.0 - flux;
						//dflux *= 2.0;
						
						if( bcf )
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
								(*u)(ix,iy+j,iz+k) += dflux;
						else
							(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop-1,iytop,iztop) + 2.0*flux;
						
					}
					// bottom boundary
					if( iy == -1 && ix%2==0 && iz%2==0 )
					{
						flux = 0.0;
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
							{
								ustar1 = interp2( (*utop)(ixtop-1,iytop,iztop-1),(*utop)(ixtop,iytop,iztop-1),(*utop)(ixtop+1,iytop,iztop-1), fac*(j-0.5) );
								ustar2 = interp2( (*utop)(ixtop-1,iytop,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop+1,iytop,iztop), fac*(j-0.5) );
								ustar3 = interp2( (*utop)(ixtop-1,iytop,iztop+1),(*utop)(ixtop,iytop,iztop+1),(*utop)(ixtop+1,iytop,iztop+1), fac*(j-0.5) );
								
								uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
								
								//(*u)(ix+j,iy,iz+k) = 0.0;(*utop)(ixtop,iytop,iztop);//interp2( -1.5, 0.0, 1.0, uhat, (*u)(ix+j,iy+1,iz+k), (*u)(ix+j,iy+2,iz+k), -1.0 );
								(*u)(ix+j,iy,iz+k) = interp2left( uhat, (*u)(ix+j,iy+1,iz+k), (*u)(ix+j,iy+2,iz+k) );
								
								flux += ((*u)(ix+j,iy+1,iz+k)-(*u)(ix+j,iy,iz+k));
							}
						flux /= 4.0;
						//(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop,iytop+1,iztop) - flux;
						double dflux = ((*utop)(ixtop,iytop+1,iztop)-(*utop)(ixtop,iytop,iztop))/2.0 - flux;
						//dflux *= 2.0;
						if( bcf )
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
								(*u)(ix+j,iy,iz+k) -= dflux;
						else
							(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop,iytop+1,iztop) - 2.0*flux;
						
					}
					// top boundary
					if( iy == ny && ix%2==0 && iz%2==0 )
					{
						flux = 0.0;
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
							{		
								ustar1 = interp2( (*utop)(ixtop-1,iytop,iztop-1),(*utop)(ixtop,iytop,iztop-1),(*utop)(ixtop+1,iytop,iztop-1), fac*(j-0.5) );
								ustar2 = interp2( (*utop)(ixtop-1,iytop,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop+1,iytop,iztop), fac*(j-0.5) );
								ustar3 = interp2( (*utop)(ixtop-1,iytop,iztop+1),(*utop)(ixtop,iytop,iztop+1),(*utop)(ixtop+1,iytop,iztop+1), fac*(j-0.5) );
								
								uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
								
								//(*u)(ix+j,iy,iz+k) = 0.0;(*utop)(ixtop,iytop,iztop);//interp2( 1.5, 0.0, -1.0, uhat, (*u)(ix+j,iy-1,iz+k), (*u)(ix+j,iy-2,iz+k), 1.0 );
								(*u)(ix+j,iy,iz+k) = interp2right( (*u)(ix+j,iy-2,iz+k), (*u)(ix+j,iy-1,iz+k), uhat  );
								
								flux += ((*u)(ix+j,iy,iz+k)-(*u)(ix+j,iy-1,iz+k));
							}
						flux /= 4.0;
						//(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop,iytop-1,iztop) + flux;
						double dflux = ((*utop)(ixtop,iytop,iztop)-(*utop)(ixtop,iytop-1,iztop))/2.0 - flux;
						//dflux *= 2.0;
						if( bcf )
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
								(*u)(ix+j,iy,iz+k) += dflux;
						else
							(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop,iytop-1,iztop) + 2.0*flux;
						
					}
					// front boundary
					if( iz == -1 && ix%2==0 && iy%2==0 )
					{
						flux = 0.0;
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
							{		
								ustar1 = interp2( (*utop)(ixtop-1,iytop-1,iztop),(*utop)(ixtop,iytop-1,iztop),(*utop)(ixtop+1,iytop-1,iztop), fac*(j-0.5) );
								ustar2 = interp2( (*utop)(ixtop-1,iytop,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop+1,iytop,iztop), fac*(j-0.5) );
								ustar3 = interp2( (*utop)(ixtop-1,iytop+1,iztop),(*utop)(ixtop,iytop+1,iztop),(*utop)(ixtop+1,iytop+1,iztop), fac*(j-0.5) );
								
								uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
								
								//(*u)(ix+j,iy+k,iz) = 0.0;(*utop)(ixtop,iytop,iztop);//interp2( -1.5, 0.0, 1.0, uhat, (*u)(ix+j,iy+k,iz+1), (*u)(ix+j,iy+k,iz+2), -1.0 );
								(*u)(ix+j,iy+k,iz) = interp2left( uhat, (*u)(ix+j,iy+k,iz+1), (*u)(ix+j,iy+k,iz+2) );
								
								flux += ((*u)(ix+j,iy+k,iz+1)-(*u)(ix+j,iy+k,iz));
							}
						flux /= 4.0;
						//(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop,iytop,iztop+1) - flux;
						double dflux = ((*utop)(ixtop,iytop,iztop+1)-(*utop)(ixtop,iytop,iztop))/2.0 - flux;
						//dflux *= 2.0;
						if( bcf )
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
								(*u)(ix+j,iy+k,iz) -= dflux;
						else
							(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop,iytop,iztop+1) - 2.0*flux;
						
					}
					// back boundary
					if( iz == nz && ix%2==0 && iy%2==0 )
					{
						flux = 0.0;
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
							{		
								ustar1 = interp2( (*utop)(ixtop-1,iytop-1,iztop),(*utop)(ixtop,iytop-1,iztop),(*utop)(ixtop+1,iytop-1,iztop), fac*(j-0.5) );
								ustar2 = interp2( (*utop)(ixtop-1,iytop,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop+1,iytop,iztop), fac*(j-0.5) );
								ustar3 = interp2( (*utop)(ixtop-1,iytop+1,iztop),(*utop)(ixtop,iytop+1,iztop),(*utop)(ixtop+1,iytop+1,iztop), fac*(j-0.5) );
								
								uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
								
								//(*u)(ix+j,iy+k,iz) = 0.0;(*utop)(ixtop,iytop,iztop);//interp2( 1.5, 0.0, -1.0, uhat, (*u)(ix+j,iy+k,iz-1), (*u)(ix+j,iy+k,iz-2), 1.0 );
								(*u)(ix+j,iy+k,iz) = interp2right( (*u)(ix+j,iy+k,iz-2), (*u)(ix+j,iy+k,iz-1), uhat );
								
								flux += ((*u)(ix+j,iy+k,iz)-(*u)(ix+j,iy+k,iz-1));
							}
						flux /= 4.0;
						//(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop,iytop,iztop-1) + flux;
						double dflux = ((*utop)(ixtop,iytop,iztop)-(*utop)(ixtop,iytop,iztop-1))/2.0 - flux;
						//dflux *= 2.0;
						if( bcf )
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
								(*u)(ix+j,iy+k,iz) += dflux;
						else
							(*utop)(ixtop,iytop,iztop) = (*utop)(ixtop,iytop,iztop-1) + 2.0*flux;
					}
					
				}
			}
	
}


#if 1
template< class S, class O, typename T >
void solver<S,O,T>::setBC( unsigned ilevel )
{
	//... set only on level before additional refinement starts
	//if( ilevel == m_ilevelmin )
	if( ilevel == m_ilevelmin )
	{
		MeshvarBnd<T> *u = m_pu->get_grid(ilevel);
		//int nbnd = u->m_nbnd,
		int
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
		
		/*for( int ix=-nbnd; ix<nx+nbnd; ++ix )
			for( int iy=-nbnd; iy<ny+nbnd; ++iy )
				for( int iz=-nbnd; iz<nz+nbnd; ++iz )
					if( ix<0||ix>=nx||iy<0||iy>=ny||iz<0||iz>=nz )
						(*u)(ix,iy,iz) = (*m_pubnd)(ix,iy,iz);*/
		
		
		for( int iy=0; iy<ny; ++iy )
			for( int iz=0; iz<nz; ++iz )
			{
				(*u)(-1,iy,iz) = 2.0*(*m_pubnd)(-1,iy,iz) - (*u)(0,iy,iz);
				(*u)(nx,iy,iz) = 2.0*(*m_pubnd)(nx,iy,iz) - (*u)(nx-1,iy,iz);;
			}
		
		for( int ix=0; ix<nx; ++ix )
			for( int iz=0; iz<nz; ++iz )
			{
				(*u)(ix,-1,iz) = 2.0*(*m_pubnd)(ix,-1,iz) - (*u)(ix,0,iz);
				(*u)(ix,ny,iz) = 2.0*(*m_pubnd)(ix,ny,iz) - (*u)(ix,ny-1,iz);
			}
		
		for( int ix=0; ix<nx; ++ix )
			for( int iy=0; iy<ny; ++iy )
			{
				(*u)(ix,iy,-1) = 2.0*(*m_pubnd)(ix,iy,-1) - (*u)(ix,iy,0);
				(*u)(ix,iy,nz) = 2.0*(*m_pubnd)(ix,iy,nz) - (*u)(ix,iy,nz-1);
			}		
		
		
		
	}/*else if( ilevel < m_ilevelmin ) {
		MeshvarBnd<T> *u = m_pu->get_grid(ilevel);
		int nbnd = u->m_nbnd,
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
		
		for( int ix=-nbnd; ix<nx+nbnd; ++ix )
			for( int iy=-nbnd; iy<ny+nbnd; ++iy )
				for( int iz=-nbnd; iz<nz+nbnd; ++iz )
					if( ix<0||ix>=nx||iy<0||iy>=ny||iz<0||iz>=nz )
						(*u)(ix,iy,iz) = 0.0;
	}*/


}

#else

//... enforce periodic boundary conditions
template< class S, class O, typename T >
void solver<S,O,T>::setBC( unsigned ilevel )
{
	MeshvarBnd<T> *u = m_pu->get_grid(ilevel);
	
	//... set only on level before additional refinement starts
	if( ilevel <= m_ilevelmin )
	{
		
		int nbnd = u->m_nbnd,
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
		
		//(*u)(0,0,0) = 0.0;
		
		
		double sum = 0.0;
		for( int ix=0; ix<nx; ++ix )
			for( int iy=0; iy<ny; ++iy )
				for( int iz=0; iz<nz; ++iz )
					sum += (*u)(ix,iy,iz);
		sum /= (nx*ny*nz);
		for( int ix=0; ix<nx; ++ix )
			for( int iy=0; iy<ny; ++iy )
				for( int iz=0; iz<nz; ++iz )
					(*u)(ix,iy,iz) -= sum;
		
		
		
		
		for( int iy=0; iy<ny; ++iy )
			for( int iz=0; iz<nz; ++iz )
			{
				(*u)(-1,iy,iz) = (*u)(nx-1,iy,iz);
				(*u)(nx,iy,iz) = (*u)(0,iy,iz);
			}
		
		for( int ix=0; ix<nx; ++ix )
			for( int iz=0; iz<nz; ++iz )
			{
				(*u)(ix,-1,iz) = (*u)(ix,ny-1,iz);
				(*u)(ix,ny,iz) = (*u)(ix,0,iz);
			}
					
		for( int ix=0; ix<nx; ++ix )
			for( int iy=0; iy<ny; ++iy )
			{
				(*u)(ix,iy,-1) = (*u)(ix,iy,nz-1);
				(*u)(ix,iy,nz) = (*u)(ix,iy,0);
			}	
		
		
		
	}
	
}
#endif


//... enforce periodic boundary conditions
template< class S, class O, typename T >
void solver<S,O,T>::make_periodic( MeshvarBnd<T> *u )
{

	int
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
		

	#pragma omp parallel
	{
		
		if( u->offset(0) == 0 )
			for( int iy=0; iy<ny; ++iy )
				for( int iz=0; iz<nz; ++iz )
				{
					(*u)(-1,iy,iz) = (*u)(nx-1,iy,iz);
					(*u)(nx,iy,iz) = (*u)(0,iy,iz);
				}
		
		if( u->offset(1) == 0 )
			for( int ix=0; ix<nx; ++ix )
				for( int iz=0; iz<nz; ++iz )
				{
					(*u)(ix,-1,iz) = (*u)(ix,ny-1,iz);
					(*u)(ix,ny,iz) = (*u)(ix,0,iz);
				}
		
		if( u->offset(2) == 0 )
			for( int ix=0; ix<nx; ++ix )
				for( int iy=0; iy<ny; ++iy )
				{
					(*u)(ix,iy,-1) = (*u)(ix,iy,nz-1);
					(*u)(ix,iy,nz) = (*u)(ix,iy,0);
				}									  
			
	}
}

END_MULTIGRID_NAMESPACE
 
#endif
