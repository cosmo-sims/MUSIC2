/*
 
 mg_solver.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
*/

#ifndef __MG_SOLVER_HH
#define __MG_SOLVER_HH

#include <cmath>
#include <iostream>

#include "mg_operators.hh"
#include "mg_interp.hh"

#include "mesh.hh"

#define BEGIN_MULTIGRID_NAMESPACE namespace multigrid {
#define END_MULTIGRID_NAMESPACE }

BEGIN_MULTIGRID_NAMESPACE
	
//! options for multigrid smoothing operation
namespace opt {
	enum smtype { sm_jacobi, sm_gauss_seidel, sm_sor };
}


//! actual implementation of FAS adaptive multigrid solver
template< class S, class I, class O, typename T=double >
class solver
{
public:
	typedef S scheme;
	typedef O mgop;
	typedef I interp;

protected:
	scheme				m_scheme;				//!< finite difference scheme
	mgop				m_gridop;				//!< grid prolongation and restriction operator
	unsigned			m_npresmooth,			//!< number of pre sweeps
						m_npostsmooth;			//!< number of post sweeps
	opt::smtype			m_smoother;				//!< smoothing method to be applied
	unsigned			m_ilevelmin;			//!< index of the top grid level
	
	const static bool	m_bperiodic = true;		//!< flag whether top grid is periodic
	
	std::vector<double> m_residu_ini;			//!< vector of initial residuals for each level
	bool m_is_ini;								//!< bool that is true for first iteration

	GridHierarchy<T>	*m_pu,					//!< pointer to GridHierarchy for solution u
						*m_pf,					//!< pointer to GridHierarchy for right-hand-side
						*m_pfsave;				//!< pointer to saved state of right-hand-side (unused)
	
	const MeshvarBnd<T> *m_pubnd;
	
	//! compute residual for a level
	double compute_error( const MeshvarBnd<T>& u, const MeshvarBnd<T>& unew );
	
	//! compute residuals for entire grid hierarchy
	double compute_error( const GridHierarchy<T>& uh, const GridHierarchy<T>& uhnew, bool verbose );
	
	//! compute residuals for entire grid hierarchy
	double compute_RMS_resid( const GridHierarchy<T>& uh, const GridHierarchy<T>& fh, bool verbose );

protected:
	
	//! Jacobi smoothing 
	void Jacobi( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f );
	
	//! Gauss-Seidel smoothing
	void GaussSeidel( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f );
	
	//! Successive-Overrelaxation smoothing
	void SOR( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f );
	
	//! main two-grid (V-cycle) for multi-grid iterations
	void twoGrid( unsigned ilevel );
	
	//! apply boundary conditions
	void setBC( unsigned ilevel );
	
	//! make top grid periodic boundary conditions
	void make_periodic( MeshvarBnd<T> *u );
	
	//void interp_coarse_fine_cubic( unsigned ilevel, MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine );
		
public:
	
	//! constructor
	solver( GridHierarchy<T>& f, opt::smtype smoother, unsigned npresmooth, unsigned npostsmooth );
	
	//! destructor
	~solver()
	{  }
	
	//! solve Poisson's equation 
	double solve( GridHierarchy<T>& u, double accuracy, double h=-1.0, bool verbose=false );
	
	//! solve Poisson's equation 
	double solve( GridHierarchy<T>& u, double accuracy, bool verbose=false )
	{
		return this->solve ( u, accuracy, -1.0, verbose );
	}
	
	
	
};


template< class S, class I, class O, typename T >
solver<S,I,O,T>::solver( GridHierarchy<T>& f, opt::smtype smoother, unsigned npresmooth, unsigned npostsmooth )
:	m_scheme(), m_gridop(), m_npresmooth( npresmooth ), m_npostsmooth( npostsmooth ), 
m_smoother( smoother ), m_ilevelmin( f.levelmin() ), m_is_ini( true ), m_pf( &f )
{ 
	m_is_ini = true;
	
	// TODO: maybe later : add more than one refinement region, then we need the mask
	//... initialize the refinement mask
	//m_pmask = new GridHierarchy<bool>( f.m_nbnd );
	//m_pmask->create_base_hierarchy(f.levelmin());
			
	/*for( unsigned ilevel=f.levelmin()+1; ilevel<=f.levelmax(); ++ilevel )
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
	*/	
}


template< class S, class I, class O, typename T >
void solver<S,I,O,T>::Jacobi( T h, MeshvarBnd<T> *u, const MeshvarBnd<T>* f )
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

template< class S, class I, class O, typename T >
void solver<S,I,O,T>::SOR( T h, MeshvarBnd<T> *u, const MeshvarBnd<T>* f )
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
	//alpha = 2 / (1 + 4 * atan(1.0) / double(u->size(0)))-1.0, //.. ideal alpha
		ialpha = 1.0-alpha;
	
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

template< class S, class I, class O, typename T >
void solver<S,I,O,T>::GaussSeidel( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f )
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


template< class S, class I, class O, typename T >
void solver<S,I,O,T>::twoGrid( unsigned ilevel )
{
	MeshvarBnd<T> *uf, *uc, *ff, *fc;
	
	T 
		h = 1.0/(1<<ilevel),
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
			interp().interp_coarse_fine(ilevel,*uc,*uf);
		
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
	else if( ilevel > m_ilevelmin )
		interp().interp_coarse_fine(ilevel,*uc,*uf);
		
	

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
	
	int oi, oj, ok;
	oi = ff->offset(0);
	oj = ff->offset(1);
	ok = ff->offset(2);
	
	#pragma omp parallel for
	for( int ix=oi; ix<oi+(int)ff->size(0)/2; ++ix )
		for( int iy=oj; iy<oj+(int)ff->size(1)/2; ++iy )
			for( int iz=ok; iz<ok+(int)ff->size(2)/2; ++iz )
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

	if( m_bperiodic && ilevel <= m_ilevelmin )
		make_periodic( &cc );

	m_gridop.prolong_add( cc, *uf );
	
	//... interpolate and apply coarse-fine boundary conditions on fine level
	if( m_bperiodic && ilevel <= m_ilevelmin )
		make_periodic( uf );
	else if(!m_bperiodic)
		setBC( ilevel );
	
	//... do smoothing sweeps with specified solver
	for( unsigned i=0; i<m_npostsmooth; ++i ){
		
		if( ilevel > m_ilevelmin )
			interp().interp_coarse_fine(ilevel,*uc,*uf);

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

template< class S, class I, class O, typename T >
double solver<S,I,O,T>::compute_error( const MeshvarBnd<T>& u, const MeshvarBnd<T>& unew )
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

template< class S, class I, class O, typename T >
double solver<S,I,O,T>::compute_error( const GridHierarchy<T>& uh, const GridHierarchy<T>& uhnew, bool verbose )
{
	double maxerr = 0.0;
	
	for( unsigned ilevel=uh.levelmin(); ilevel <= uh.levelmax(); ++ilevel )
	{
		double err = 0.0;
		err = compute_error( *uh.get_grid(ilevel), *uhnew.get_grid(ilevel) );
		
		if( verbose )
			std::cout << "      Level " << std::setw(6) << ilevel << ",   Error = " << err << std::endl;
		maxerr = std::max(maxerr,err);
		
	}
	return maxerr;
}

template< class S, class I, class O, typename T >
double solver<S,I,O,T>::compute_RMS_resid( const GridHierarchy<T>& uh, const GridHierarchy<T>& fh, bool verbose )
{
	if( m_is_ini )
		m_residu_ini.assign( uh.levelmax()+1, 0.0 );
	
	double maxerr=0.0;
	
	for( unsigned ilevel=uh.levelmin(); ilevel <= uh.levelmax(); ++ilevel )
	{
		int 
		nx = uh.get_grid(ilevel)->size(0), 
		ny = uh.get_grid(ilevel)->size(1), 
		nz = uh.get_grid(ilevel)->size(2);
		
		double h = 1.0/(1<<ilevel), h2=h*h, err;
		double sum = 0.0;
		unsigned count = 0;
		
		#pragma omp parallel for reduction(+:sum,count)
		for( int ix=0; ix<nx; ++ix )
			for( int iy=0; iy<ny; ++iy )
				for( int iz=0; iz<nz; ++iz )
				{
					double r = (m_scheme.apply( *uh.get_grid(ilevel), ix, iy, iz )/h2 + (*fh.get_grid(ilevel))(ix,iy,iz));
					sum += r*r;
					++count;
				}
		
		if( m_is_ini )
			m_residu_ini[ilevel] =  sqrt(sum)/count;
		
		err = sqrt(sum)/count/m_residu_ini[ilevel];
		
		if( verbose && !m_is_ini )
			std::cout << "      Level " << std::setw(6) << ilevel << ",   Error = " << err << std::endl;
		
		LOGDEBUG("[mg]      level %3d,  rms residual %g",ilevel,err);
		
		if( err > maxerr )
			maxerr = err;
		
	}
	
	if( m_is_ini )
		m_is_ini = false;
	
	return maxerr;
}


template< class S, class I, class O, typename T >
double solver<S,I,O,T>::solve( GridHierarchy<T>& uh, double acc, double h, bool verbose )
{

	double err, maxerr = 1e30;
	unsigned niter = 0;
	
	bool fullverbose = false;
	
	m_pu = &uh;
	
	//... iterate ...//
	while (true)
	{
		
		LOGUSER("Performing multi-grid V-cycle...");
		twoGrid( uh.levelmax() );
		
		err = compute_RMS_resid( *m_pu, *m_pf, fullverbose );
		++niter;
		
		if( fullverbose ){
			LOGUSER("  multigrid iteration %3d, maximum RMS residual = %g", niter, err );
			std::cout << "   - Step No. " << std::setw(3) << niter << ", Max Err = " << err << std::endl;
			std::cout << "     ---------------------------------------------------\n";
		}
		
		if( err < maxerr )
			maxerr = err;
			
		if( (niter > 1) && ((err < acc) || (niter > 20)) )
			break;
	}		
	
	if( err > acc )
	{	
		std::cout << "Error : no convergence in Poisson solver" << std::endl;
		LOGERR("No convergence in Poisson solver, final error: %g.",err);
	}
	else if( verbose )
	{	
		std::cout << " - Converged in " << niter << " steps to " << maxerr << std::endl;
		LOGUSER("Poisson solver converged to max. error of %g in %d steps.",err,niter);
	}

	
	//.. make sure that the RHS does not contain the FAS corrections any more
	for( int i=m_pf->levelmax(); i>0; --i )
		m_gridop.restrict( *m_pf->get_grid(i), *m_pf->get_grid(i-1) );
	
	
	return err;
}



//TODO: this only works for 2nd order! (but actually not needed)
template< class S, class I, class O, typename T >
void solver<S,I,O,T>::setBC( unsigned ilevel )
{
	//... set only on level before additional refinement starts
	if( ilevel == m_ilevelmin )
	{
		MeshvarBnd<T> *u = m_pu->get_grid(ilevel);
		int
			nx = u->size(0), 
			ny = u->size(1), 
			nz = u->size(2);
			
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
		
		
		
	}
}



//... enforce periodic boundary conditions
template< class S, class I, class O, typename T >
void solver<S,I,O,T>::make_periodic( MeshvarBnd<T> *u )
{
	

	int
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
	int nb = u->m_nbnd;
	
		
	//if( u->offset(0) == 0 )
		for( int iy=-nb; iy<ny+nb; ++iy )
			for( int iz=-nb; iz<nz+nb; ++iz )
			{
				int iiy( (iy+ny)%ny ), iiz( (iz+nz)%nz );
				
				for( int i=-nb; i<0; ++i )
				{
					(*u)(i,iy,iz) = (*u)(nx+i,iiy,iiz);
					(*u)(nx-1-i,iy,iz) = (*u)(-1-i,iiy,iiz);	
				}
				
			}
	
	//if( u->offset(1) == 0 )
		for( int ix=-nb; ix<nx+nb; ++ix )
			for( int iz=-nb; iz<nz+nb; ++iz )
			{
				int iix( (ix+nx)%nx ), iiz( (iz+nz)%nz );
				
				for( int i=-nb; i<0; ++i )
				{
					(*u)(ix,i,iz) = (*u)(iix,ny+i,iiz);
					(*u)(ix,ny-1-i,iz) = (*u)(iix,-1-i,iiz);
				}
			}
	
	//if( u->offset(2) == 0 )
		for( int ix=-nb; ix<nx+nb; ++ix )
			for( int iy=-nb; iy<ny+nb; ++iy )
			{
				int iix( (ix+nx)%nx ), iiy( (iy+ny)%ny );
				
				for( int i=-nb; i<0; ++i )
				{
					(*u)(ix,iy,i) = (*u)(iix,iiy,nz+i);
					(*u)(ix,iy,nz-1-i) = (*u)(iix,iiy,-1-i);
				}
			}
	
}


END_MULTIGRID_NAMESPACE
 
#endif
