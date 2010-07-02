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
	
namespace opt {
	enum smtype { sm_jacobi, sm_gauss_seidel, sm_sor };
}



template< class S, class I, class O, typename T=double >
class solver
{
public:
	typedef S scheme;
	typedef O mgop;
	typedef I interp;

protected:
	scheme				m_scheme;
	mgop				m_gridop;
	unsigned			m_npresmooth, m_npostsmooth;
	opt::smtype			m_smoother;
	unsigned			m_ilevelmin;
	
	const static bool	m_bperiodic = true;
	
	std::vector<double> m_residu_ini;
	bool m_is_ini;

	GridHierarchy<T>	*m_pu, *m_pf, *m_pfsave;	
	//GridHierarchy<bool> *m_pmask;
	const MeshvarBnd<T> *m_pubnd;
	
	double compute_error( const MeshvarBnd<T>& u, const MeshvarBnd<T>& unew );
	
	double compute_error( const GridHierarchy<T>& uh, const GridHierarchy<T>& uhnew, bool verbose );
	
	double compute_RMS_resid( const GridHierarchy<T>& uh, const GridHierarchy<T>& fh, bool verbose );

protected:
	
	void Jacobi( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f );
	
	void GaussSeidel( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f );
	
	void SOR( T h, MeshvarBnd<T>* u, const MeshvarBnd<T>* f );
	
	void twoGrid( unsigned ilevel );
	
	//void interp_coarse_fine( unsigned ilevel, MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine, bool bcf=true );
	
	void setBC( unsigned ilevel );
	
	void make_periodic( MeshvarBnd<T> *u );
	
	void interp_coarse_fine_cubic( unsigned ilevel, MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine );
		
public:
	solver( GridHierarchy<T>& f, //const MeshvarBnd<T>& uBC_top, 
				   opt::smtype smoother, unsigned npresmooth, unsigned npostsmooth );
	
	~solver()
	{ /*delete m_pmask;*/ }
	
	double solve( GridHierarchy<T>& u, double accuracy, double h=-1.0, bool verbose=false );
	
	double solve( GridHierarchy<T>& u, double accuracy, bool verbose=false )
	{
		return this->solve ( u, accuracy, -1.0, verbose );
	}
	
	
	
};


template< class S, class I, class O, typename T >
solver<S,I,O,T>::solver( GridHierarchy<T>& f, //const MeshvarBnd<T>& ubnd, 
					opt::smtype smoother, unsigned npresmooth, unsigned npostsmooth )
:	m_scheme(), m_gridop(), m_npresmooth( npresmooth ), m_npostsmooth( npostsmooth ), 
m_smoother( smoother ), m_ilevelmin( f.levelmin() ), m_is_ini( true ), m_pf( &f )
{ 
	m_is_ini = true;
	
	// TODO: maybe later : add more than one refinement region
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
	//alpha = 2 / (1 + 4 * atan(1.0) / double(u->size(0)))-1.0,
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
	//mg_straight().restrict( Lu, tLu );
	Lu.deallocate();
	
	//... restrict source term
	m_gridop.restrict( *ff, *fc );
	//mg_straight().restrict( *ff, *fc );
	
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
		
		double h = 1.0/pow(2,ilevel), h2=h*h, err;
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

	double err;
	
	//GridHierarchy<T> uhnew(uh);//, fsave(*m_pf);
	m_pu = &uh;
	
    unsigned niter = 0;
	
	//... iterate ...//
	while (true)
	{
		
		
		twoGrid( uh.levelmax() );
		//err = compute_error( *m_pu, uhnew, verbose );
		err = compute_RMS_resid( *m_pu, *m_pf, verbose );
		++niter;
		
		if( verbose ){
			std::cout << "  --> Step No. " << std::setw(3) << niter << ", Max Err = " << err << std::endl;
			std::cout << "      ---------------------------------------------------\n";
		}
			
		if( (niter > 1) && ((err < acc) || (niter > 8)) )
			break;
		
		//uhnew = *m_pu;
		//*m_pf = fsave;
	}		
	
	if( err > acc )
		std::cout << "Error : no convergence in Poisson solver" << std::endl;
	else if( verbose )
		std::cout << " - Converged in " << niter << " steps to req. acc. of " << acc << std::endl;

	
	//uh = uhnew;
	//*m_pf = fsave;
	
	//.. make sure that the RHS does not contain the FAS corrections any more
	for( int i=m_pf->levelmax(); i>0; --i )//(int)m_pf->levelmin(); --i )
		m_gridop.restrict( *m_pf->get_grid(i), *m_pf->get_grid(i-1) );
	
	
	return err;
}

template< class S, class I, class O, typename T >
void solver<S,I,O,T>::interp_coarse_fine_cubic( unsigned ilevel, MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine )
{
	MeshvarBnd<T> *u    = &fine;
	MeshvarBnd<T> *utop = &coarse;
	
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
	
	//.. perform cubic interpolation
	mg_cubic().prolong_bnd( *utop, *u );
	
	//.. match fluxes
	#pragma omp parallel for schedule(dynamic)
	for( int ix=-1; ix<=nx; ++ix )
		for( int iy=-1; iy<=ny; ++iy )
			for( int iz=-1; iz<=nz; ++iz )
			{
				bool xbnd=(ix==-1||ix==nx),ybnd=(iy==-1||iy==ny),zbnd=(iz==-1||iz==nz);
				
				if( xbnd || ybnd || zbnd )
				{
					
					//... only deal with proper ghostzones
					if( (xbnd&&ybnd) || (xbnd&&zbnd) || (ybnd&&zbnd) || (xbnd&&ybnd&&zbnd))
						continue;
					
					int ixtop = (int)(0.5*(double)(ix))+xoff;
					int iytop = (int)(0.5*(double)(iy))+yoff;
					int iztop = (int)(0.5*(double)(iz))+zoff;
					
					if( ix==-1 ) ixtop=xoff-1;
					if( iy==-1 ) iytop=yoff-1;
					if( iz==-1 ) iztop=zoff-1;
					
					double flux;;
				    if( ix == -1 && iy%2==0 && iz%2==0 )
					{					
						flux = 0.0;
						for( int j=0;j<=1;j++)
							for( int k=0;k<=1;k++)
							{		
								flux += ((*u)(ix+1,iy+j,iz+k)-(*u)(ix,iy+j,iz+k));
							}
						
						
						
						flux /= 4.0;
						
						double dflux = ((*utop)(ixtop+1,iytop,iztop)-(*utop)(ixtop,iytop,iztop))/2.0 - flux;
					}
				}
			}
	
}





template< class S, class I, class O, typename T >
void solver<S,I,O,T>::setBC( unsigned ilevel )
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
	
#if 0
	if( u->m_nbnd == 1 )
	{
			
		#pragma omp parallel
		{
			
			if( u->offset(0) == 0 )
				for( int iy=-1; iy<=ny; ++iy )
					for( int iz=-1; iz<=nz; ++iz )
					{
						int iiy( (iy+ny)%ny ), iiz( (iz+nz)%nz );

						(*u)(-1,iy,iz) = (*u)(nx-1,iiy,iiz);
						(*u)(nx,iy,iz) = (*u)(0,iiy,iiz);
					}
			
			if( u->offset(1) == 0 )
				for( int ix=-1; ix<=nx; ++ix )
					for( int iz=-1; iz<=nz; ++iz )
					{
						int iix( (ix+nx)%nx ), iiz( (iz+nz)%nz );

						(*u)(ix,-1,iz) = (*u)(iix,ny-1,iiz);
						(*u)(ix,ny,iz) = (*u)(iix,0,iiz);
					}
			
			if( u->offset(2) == 0 )
				for( int ix=-1; ix<=nx; ++ix )
					for( int iy=-1; iy<=ny; ++iy )
					{
						int iix( (ix+nx)%nx ), iiy( (iy+ny)%ny );

						(*u)(ix,iy,-1) = (*u)(iix,iiy,nz-1);
						(*u)(ix,iy,nz) = (*u)(iix,iiy,0);
					}									  
			
			
		}
	}else if( u->m_nbnd == 2 ){
		if( u->offset(0) == 0 )
			for( int iy=-2; iy<=ny+1; ++iy )
				for( int iz=-2; iz<=nz+1; ++iz )
				{
					int iiy( (iy+ny)%ny ), iiz( (iz+nz)%nz );
					
					(*u)(-1,iy,iz) = (*u)(nx-1,iiy,iiz);
					(*u)(-2,iy,iz) = (*u)(nx-2,iiy,iiz);
					(*u)(nx,iy,iz) = (*u)(0,iiy,iiz);
					(*u)(nx+1,iy,iz) = (*u)(1,iiy,iiz);
				}
		
		if( u->offset(1) == 0 )
			for( int ix=-2; ix<=nx+1; ++ix )
				for( int iz=-2; iz<=nz+1; ++iz )
				{
					int iix( (ix+nx)%nx ), iiz( (iz+nz)%nz );
					
					(*u)(ix,-1,iz) = (*u)(iix,ny-1,iiz);
					(*u)(ix,-2,iz) = (*u)(iix,ny-2,iiz);
					(*u)(ix,ny,iz) = (*u)(iix,0,iiz);
					(*u)(ix,ny+1,iz) = (*u)(iix,1,iiz);
				}
		
		if( u->offset(2) == 0 )
			for( int ix=-2; ix<=nx+1; ++ix )
				for( int iy=-2; iy<=ny+1; ++iy )
				{
					int iix( (ix+nx)%nx ), iiy( (iy+ny)%ny );
					
					(*u)(ix,iy,-1) = (*u)(iix,iiy,nz-1);
					(*u)(ix,iy,-2) = (*u)(iix,iiy,nz-2);
					(*u)(ix,iy,nz) = (*u)(iix,iiy,0);
					(*u)(ix,iy,nz+1) = (*u)(iix,iiy,1);
				}			
	}else
		throw std::runtime_error("solver<S,I,O,T>::make_periodic: don't know how to deal with boundary!");
#endif
}

#if 0
template< class S, class I, class O, typename T >
void solver<S,I,O,T>::interp_coarse_fine( unsigned ilevel, MeshvarBnd<T>& coarse, MeshvarBnd<T>& fine, bool bcf )
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
#endif

END_MULTIGRID_NAMESPACE
 
#endif
