/*
 
 constraints.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include "constraints.hh"

double dsigma2_tophat( double k, void *params );
double dsigma2_gauss( double k, void *params );
double find_coll_z( const std::vector<double>& z, const std::vector<double>& sigma, double nu );
void compute_sigma_tophat( config_file& cf, transfer_function *ptf, double R, std::vector<double>& z, std::vector<double>& sigma );
void compute_sigma_gauss( config_file& cf, transfer_function *ptf, double R, std::vector<double>& z, std::vector<double>& sigma );



double dsigma2_tophat( double k, void *pparams )
{
	if( k<=0.0 )
		return 0.0;
	
	char **params = reinterpret_cast<char**> (pparams);
	
	transfer_function *ptf = reinterpret_cast<transfer_function*>(params[0]);
	double x = k * *reinterpret_cast<double*>(params[1]);
	double nspect = *reinterpret_cast<double*>(params[2]);
	
	double w = 3.0*(sin(x)-x*cos(x))/(x*x*x);
	
	double tfk = ptf->compute(k,total);
	return k*k * w*w * pow(k,nspect) * tfk*tfk;
}

double dsigma2_gauss( double k, void *pparams )
{
	if( k<=0.0 )
		return 0.0;
	
	char **params = reinterpret_cast<char**> (pparams);
	
	transfer_function *ptf = reinterpret_cast<transfer_function*> (params[0]);
	double x = k * *reinterpret_cast<double*>(params[1]);
	double nspect = *reinterpret_cast<double*>(params[2]);
	
	double w = exp(-x*x*0.5);
	
	double tfk = ptf->compute(k,total);
	return k*k * w*w * pow(k,nspect) * tfk*tfk;
}

double find_coll_z( const std::vector<double>& z, const std::vector<double>& sigma, double nu )
{
	double dcoll = 1.686/nu;
	double zcoll = 0.0;
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, z.size());
	
	gsl_spline_init (spline, &sigma[0], &z[0], z.size() );
	
	zcoll = gsl_spline_eval(spline, dcoll, acc );
	
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
	
	return zcoll;
}


void compute_sigma_tophat( config_file& cf, transfer_function *ptf, double R, std::vector<double>& z, std::vector<double>& sigma )
{
	z.clear();
	sigma.clear();
	
	cosmology cosm( cf );
	CosmoCalc ccalc( cosm, ptf );
	
	double zmin = 0.0, zmax = 200.0;
	int nz = 100;
	for( int i=0; i <nz; ++i )
		z.push_back( zmax - i*(zmax-zmin)/(nz-1.0) );
	
	double D0 = ccalc.CalcGrowthFactor(1.0);

	double sigma8 = cf.getValue<double>("cosmology","sigma_8"); 
	double nspec = cf.getValue<double>("cosmology","nspec");
	
	double sigma0 = 0.0;
	{
		double eight=8.0;

		char *params[3];
		params[0] = reinterpret_cast<char*> (ptf);
		params[1] = reinterpret_cast<char*> (&eight);
		params[2] = reinterpret_cast<char*> (&nspec);
		
		sigma0 = sqrt(4.0*M_PI*integrate( &dsigma2_tophat, 1e-4, 1e4, reinterpret_cast<void*>(params) ));
	}
	
	for( int i=0; i <nz; ++i )
	{
		void *params[3];
		params[0] = reinterpret_cast<char*> (ptf);
		params[1] = reinterpret_cast<char*> (&R);
		params[2] = reinterpret_cast<char*> (&nspec);
		
		double sig = sqrt(4.0*M_PI*integrate( &dsigma2_tophat, 1e-4, 1e4, reinterpret_cast<void*>(params) ));
		double Dz  = ccalc.CalcGrowthFactor(1./(1.+z[i]));
		sigma.push_back( sig*sigma8/sigma0*Dz/D0 );
	}
}

void compute_sigma_gauss( config_file& cf, transfer_function *ptf, double R, std::vector<double>& z, std::vector<double>& sigma )
{
	z.clear();
	sigma.clear();
	
	cosmology cosm( cf );
	CosmoCalc ccalc( cosm, ptf );
	
	double zmin = 0.0, zmax = 200.0;
	int nz = 100;
	for( int i=0; i <nz; ++i )
		z.push_back( zmax - i*(zmax-zmin)/(nz-1.0) );
	
	double D0 = ccalc.CalcGrowthFactor(1.0);
	
	double sigma8 = cf.getValue<double>("cosmology","sigma_8"); 
	double nspec = cf.getValue<double>("cosmology","nspec");
	
	double sigma0 = 0.0;
	{
		double eight=8.0;
		
		char *params[3];
		params[0] = reinterpret_cast<char*> (ptf);
		params[1] = reinterpret_cast<char*> (&eight);
		params[2] = reinterpret_cast<char*> (&nspec);
		
		sigma0 = sqrt(4.0*M_PI*integrate( &dsigma2_tophat, 1e-4, 1e4, reinterpret_cast<void*>(params) ));
	}
	
	for( int i=0; i <nz; ++i )
	{
		void *params[3];
		params[0] = reinterpret_cast<char*> (ptf);
		params[1] = reinterpret_cast<char*> (&R);
		params[2] = reinterpret_cast<char*> (&nspec);
		
		double sig = sqrt(4.0*M_PI*integrate( &dsigma2_gauss, 1e-4, 1e4, reinterpret_cast<void*>(params) ));
		double Dz  = ccalc.CalcGrowthFactor(1./(1.+z[i]));
		
		//std::cerr << z[i] << "    " << sig << std::endl;
		sigma.push_back( sig*sigma8/sigma0*Dz/D0 );
	}
}


constraint_set::constraint_set( config_file& cf, transfer_function *ptf )
: pcf_( &cf ), ptf_( ptf )
{
	pcosmo_ = new Cosmology( cf );
	pccalc_ = new CosmoCalc( *pcosmo_, ptf_ );
	dplus0_ = 1.0;//pccalc_->CalcGrowthFactor( 1.0 );
	
	
	unsigned i=0;
	
	unsigned levelmin_TF = pcf_->getValue<unsigned>("setup","levelmin_TF");
	constr_level_ = pcf_->getValueSafe<unsigned>("constraints","level",levelmin_TF);
	
	constr_level_ = std::max(constr_level_,levelmin_TF);
	
	double omegam = pcf_->getValue<double>("cosmology","Omega_m");
	double rhom = omegam*2.77519737e11; //... mean matter density
	
	//... use EdS density for estimation
	//double rhom = 2.77519737e11;
	
	std::map< std::string, constr_type> constr_type_map;
	constr_type_map.insert( std::pair<std::string,constr_type>("halo",halo) );
	constr_type_map.insert( std::pair<std::string,constr_type>("peak",peak) );
	
	while(true)
	{
		char temp1[128];
		std::string temp2;
		sprintf(temp1,"constraint[%u].type",i);
		if( cf.containsKey( "constraints", temp1 ) )
		{
			std::string str_type = cf.getValue<std::string>( "constraints", temp1 );
			if( constr_type_map.find(str_type) == constr_type_map.end() )
				throw std::runtime_error("Unknown constraint type!\n");
			
			//... parse a new constraint
			constraint new_c;
			
			new_c.type = constr_type_map[ str_type ];
			
			//... read position of constraint
			sprintf(temp1,"constraint[%u].pos",i);
			temp2 = cf.getValue<std::string>( "constraints", temp1 );
			sscanf(temp2.c_str(), "%lf,%lf,%lf", &new_c.x, &new_c.y, &new_c.z);
			
			if( new_c.type == halo)
			{
				//.. halo type constraints take mass and collapse redshift
				sprintf(temp1,"constraint[%u].mass",i);
				double mass = cf.getValue<double>( "constraints", temp1 );
				
				sprintf(temp1,"constraint[%u].zform",i);
				double zcoll = cf.getValue<double>( "constraints", temp1 );
				new_c.Rg = pow((mass/pow(2.*M_PI,1.5)/rhom),1./3.);
				
				new_c.sigma = 1.686/(pccalc_->CalcGrowthFactor(1./(1.+zcoll))/pccalc_->CalcGrowthFactor(1.0));
				
				LOGINFO("Constraint %d : halo with %g h-1 M_o",i,pow(2.*M_PI,1.5)*rhom*pow(new_c.Rg,3));
			}
			else if( new_c.type == peak )
			{
				//... peak type constraints take a scale and a peak height
				//sprintf(temp1,"constraint[%u].Rg",i);
				//new_c.Rg = cf.getValue<double>( "constraints", temp1 );
				//double mass = pow(new_c.Rg,3.0)*rhom*pow(2.*M_PI,1.5);
				
				sprintf(temp1,"constraint[%u].mass",i);
				double mass = cf.getValue<double>( "constraints", temp1 );
				new_c.Rg = pow((mass/pow(2.*M_PI,1.5)/rhom),1./3.);
				double Rtophat = pow(mass/4.0*3.0/M_PI/rhom,1./3.);
				sprintf(temp1,"constraint[%u].nu",i);
				double nu = cf.getValue<double>( "constraints", temp1 );
				
				std::vector<double> z,sigma;
				compute_sigma_tophat( cf, ptf, Rtophat, z, sigma );
				double zcoll = find_coll_z( z, sigma, nu );
				
				//LOGINFO("Probable collapse redshift for constraint %d : z = %f @ M = %g", i, zcoll,mass );
				
				compute_sigma_gauss( cf, ptf, new_c.Rg, z, sigma );
				new_c.sigma = nu*sigma.back();
				
				//LOGINFO("Constraint %d : peak with Rg=%g h-1 Mpc and nu = %g",i,new_c.Rg,new_c.sigma);
				LOGINFO("Constraint %3d : peak",i);
				LOGINFO("   M = %g h-1 M_o, nu = %.2f sigma", mass, nu );
				LOGINFO("   estimated z_coll = %f, sigma = %f", zcoll, new_c.sigma );
				
			}
			
			new_c.Rg2 = new_c.Rg*new_c.Rg;
			
			cset_.push_back( new_c );
			
		}else
			break;
		
		++i;
	}
	
	LOGINFO("Found %d density constraint(s) to be obeyed.",cset_.size());
}


void constraint_set::wnoise_constr_corr( double dx, size_t nx, size_t ny, size_t nz, std::vector<double>& g0, matrix& cinv, fftw_complex* cw )
{
	double lsub = nx*dx;
	double dk = 2.0*M_PI/lsub, d3k=dk*dk*dk;
	
	double pnorm = pcf_->getValue<double>("cosmology","pnorm");
	double nspec = pcf_->getValue<double>("cosmology","nspec");
	pnorm *= dplus0_*dplus0_;
	
	size_t nconstr = cset_.size();
	size_t nzp=nz/2+1;
	
	
	
	/*for( size_t i=0; i<nconstr; ++i )
	 for( size_t j=0; j<nconstr; ++j )
	 {
	 std::cerr << "fact    = " << (cset_[j].sigma-g0[j])*cinv(i,j) << "\n";
	 std::cerr << "g(j)    = " << cset_[j].sigma << "\n";
	 std::cerr << "g0(j)   = " << g0[j] << "\n";
	 std::cerr << "qinv    = " << cinv(i,j) << "\n";
	 }
	 */
	
	
	double chisq = 0.0, chisq0 = 0.0;
	for( size_t i=0; i<nconstr; ++i )
		for( size_t j=0; j<nconstr; ++j )
		{
			chisq += cset_[i].sigma*cinv(i,j)*cset_[j].sigma;
			chisq0 += g0[i]*cinv(i,j)*g0[j];
		}
	LOGINFO("Chi squared for the constraints:\n       sampled = %f, desired = %f", chisq0, chisq );
	
	std::vector<double> sigma(nconstr,0.0);
	
	#pragma omp parallel 
	{
		std::vector<double> sigma_loc(nconstr,0.0);
		
		#pragma omp for 
		for( int ix=0; ix<(int)nx; ++ix )
		{	
			double iix(ix); if( iix > nx/2 ) iix-=nx;
			iix *= 2.0*M_PI/nx;
			
			for( size_t iy=0; iy<ny; ++iy )
			{	
				double iiy(iy); if( iiy > ny/2 ) iiy-=ny;
				iiy *= 2.0*M_PI/nx;
				for( size_t iz=0; iz<nzp; ++iz )
				{
					double iiz(iz);
					iiz *= 2.0*M_PI/nx;
					
					double k = sqrt(iix*iix+iiy*iiy+iiz*iiz)*(double)nx/lsub;
					
					double T = ptf_->compute(k,total);
					double Pk = pnorm*T*T*pow(k,nspec)*d3k;
					
					size_t q = ((size_t)ix*ny+(size_t)iy)*nzp+(size_t)iz;
					
					double fac = sqrt(Pk);
					
					for( unsigned i=0; i<nconstr; ++i )
						for( unsigned j=0; j<=i; ++j )
						{
							std::complex<double> 
							ci = eval_constr(i,iix,iiy,iiz),
							cj = eval_constr(j,iix,iiy,iiz);
							
							RE(cw[q]) += (cset_[j].sigma-g0[j])*cinv(i,j) * std::real(ci)*fac;
							IM(cw[q]) += (cset_[j].sigma-g0[j])*cinv(i,j) * std::imag(ci)*fac;
							
							if( i!=j )
							{
								RE(cw[q]) += (cset_[i].sigma-g0[i])*cinv(j,i) * std::real(cj)*fac;
								IM(cw[q]) += (cset_[i].sigma-g0[i])*cinv(j,i) * std::imag(cj)*fac;								
							}
							else
							{
								if( iz>0&&iz<nz/2 )
									sigma_loc[i] += 2.0*std::real(std::conj(ci)*std::complex<double>(RE(cw[q]),IM(cw[q])))*fac;
								else
									sigma_loc[i] += std::real(std::conj(ci)*std::complex<double>(RE(cw[q]),IM(cw[q])))*fac;
							}
						}
				}
				
			}
			
		}
		
		//.. 'critical' section for the global reduction
		#pragma omp critical
		{
			for(int i=0; i<(int)nconstr; ++i )
				sigma[i] += sigma_loc[i];
		}
	}
	
	for(int i=0; i<(int)nconstr; ++i )
		LOGINFO("Constraint %3d : sigma = %+6f (%+6f)",i,sigma[i],cset_[i].sigma);
}



void constraint_set::wnoise_constr_corr( double dx, fftw_complex* cw, size_t nx, size_t ny, size_t nz, std::vector<double>& g0 )
{
	size_t nconstr = cset_.size();
	size_t nzp=nz/2+1;
	
	g0.assign(nconstr,0.0);
	
	double pnorm = pcf_->getValue<double>("cosmology","pnorm");
	double nspec = pcf_->getValue<double>("cosmology","nspec");
	pnorm *= dplus0_*dplus0_;
	double lsub = nx*dx;
	double dk = 2.0*M_PI/lsub, d3k=dk*dk*dk;
	
	for( size_t i=0; i<nconstr; ++i )
	{
		double gg = 0.0;
		
		#pragma omp parallel for reduction(+:gg)
		for( int ix=0; ix<(int)nx; ++ix )
		{	
			double iix(ix); if( iix > nx/2 ) iix-=nx;
			iix *= 2.0*M_PI/nx;
			
			for( size_t iy=0; iy<ny; ++iy )
			{	
				double iiy(iy); if( iiy > ny/2 ) iiy-=ny;
				iiy *= 2.0*M_PI/nx;
				for( size_t iz=0; iz<nzp; ++iz )
				{
					double iiz(iz);
					iiz *= 2.0*M_PI/nx;
					
					double k = sqrt(iix*iix+iiy*iiy+iiz*iiz)*(double)nx/lsub;
					double T = ptf_->compute(k,total);
					
					std::complex<double> v(std::conj(eval_constr(i,iix,iiy,iiz)));
					
					v *= sqrt(pnorm*pow(k,nspec)*T*T*d3k);
					
					
					if( iz>0&&iz<nz/2)
						v*=2;
					
					size_t q = ((size_t)ix*ny+(size_t)iy)*nzp+(size_t)iz;

					std::complex<double> ccw(RE(cw[q]),IM(cw[q]));

					gg += std::real(v*ccw);
					
				}
			}
		}
		
		g0[i] = gg;
	}
}




void constraint_set::icov_constr( double dx, size_t nx, size_t ny, size_t nz, matrix& cij )
{
	size_t nconstr = cset_.size();
	size_t nzp=nz/2+1;
	
	double pnorm = pcf_->getValue<double>("cosmology","pnorm");
	double nspec = pcf_->getValue<double>("cosmology","nspec");
	pnorm *= dplus0_*dplus0_;
	
	cij		= matrix(nconstr,nconstr);
	
	double lsub = nx*dx;
	double dk = 2.0*M_PI/lsub, d3k=dk*dk*dk;
	
	//... compute lower triangle of covariance matrix
	//... and fill in upper triangle
	for( unsigned i=0; i<nconstr; ++i )
		for( unsigned j=0; j<=i; ++j )
		{
			
			float c1(0.0), c2(0.0);
			
#pragma omp parallel for reduction(+:c1,c2)
			for( int ix=0; ix<(int)nx; ++ix )
			{	
				double iix(ix); if( iix > nx/2 ) iix-=nx;
				iix *= 2.0*M_PI/nx;
				
				for( size_t iy=0; iy<ny; ++iy )
				{	
					double iiy(iy); if( iiy > ny/2 ) iiy-=ny;
					iiy *= 2.0*M_PI/nx;
					for( size_t iz=0; iz<nzp; ++iz )
					{
						double iiz(iz);
						iiz *= 2.0*M_PI/nx;
						
						double k = sqrt(iix*iix+iiy*iiy+iiz*iiz)*(double)nx/lsub;
						double T = ptf_->compute(k,total);
						std::complex<double> v(std::conj(eval_constr(i,iix,iiy,iiz)));
						v *= eval_constr(j,iix,iiy,iiz);
						v *= pnorm * pow(k,nspec) * T * T * d3k;
						
						if( iz>0&&iz<nz/2)
							v*=2;
						
						c1 += std::real(v);
						c2 += std::real(std::conj(v));
					}
				}
			}
			
			cij(i,j) = c1;
			cij(j,i) = c2;
		}
	
	//... invert convariance matrix
	cij.invert();
	
}

