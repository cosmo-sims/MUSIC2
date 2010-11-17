/*
 
 transfer_eisenstein.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 */

#include "transfer_function.hh"

//! Implementation of abstract base class TransferFunction for the Eisenstein & Hu transfer function 
/*!
 This class implements the analytical fit to the matter transfer
 function by Eisenstein & Hu (1999). In fact it is their code.
 */
class transfer_eisenstein_plugin : public transfer_function_plugin
{
protected:
	using transfer_function_plugin::cosmo_;
	
	//Cosmology m_Cosmology;
	double  m_h0;
	double	omhh,		/* Omega_matter*h^2 */
	obhh,		/* Omega_baryon*h^2 */
	theta_cmb,	/* Tcmb in units of 2.7 K */
	z_equality,	/* Redshift of matter-radiation equality, really 1+z */
	k_equality,	/* Scale of equality, in Mpc^-1 */
	z_drag,		/* Redshift of drag epoch */
	R_drag,		/* Photon-baryon ratio at drag epoch */
	R_equality,	/* Photon-baryon ratio at equality epoch */
	sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
	k_silk,		/* Silk damping scale, in Mpc^-1 */
	alpha_c,	/* CDM suppression */
	beta_c,		/* CDM log shift */
	alpha_b,	/* Baryon suppression */
	beta_b,		/* Baryon envelope shift */
	beta_node,	/* Sound horizon shift */
	k_peak,		/* Fit to wavenumber of first peak, in Mpc^-1 */
	sound_horizon_fit,	/* Fit to sound horizon, in Mpc */
	alpha_gamma;	/* Gamma suppression in approximate TF */
	
	//! private member function: sets internal quantities for Eisenstein & Hu fitting
	void TFset_parameters(double omega0hh, double f_baryon, double Tcmb)
	/* Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
	/* Input: omega0hh -- The density of CDM and baryons, in units of critical dens,
     multiplied by the square of the Hubble constant, in units
     of 100 km/s/Mpc */
	/* 	  f_baryon -- The fraction of baryons to CDM */
	/*        Tcmb -- The temperature of the CMB in Kelvin.  Tcmb<=0 forces use
	 of the COBE value of  2.728 K. */
	/* Output: Nothing, but set many global variables used in TFfit_onek(). 
     You can access them yourself, if you want. */
	/* Note: Units are always Mpc, never h^-1 Mpc. */
	{
		double z_drag_b1, z_drag_b2;
		double alpha_c_a1, alpha_c_a2, beta_c_b1, beta_c_b2, alpha_b_G, y;
		
		if (f_baryon<=0.0 || omega0hh<=0.0) {
			fprintf(stderr, "TFset_parameters(): Illegal input.\n");
			exit(1);
		}
		omhh = omega0hh;
		obhh = omhh*f_baryon;
		if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
		theta_cmb = Tcmb/2.7;
		
		z_equality = 2.50e4*omhh/POW4(theta_cmb);  /* Really 1+z */
		k_equality = 0.0746*omhh/SQR(theta_cmb);
		
		z_drag_b1 = 0.313*pow((double)omhh,-0.419)*(1+0.607*pow((double)omhh,0.674));
		z_drag_b2 = 0.238*pow((double)omhh,0.223);
		z_drag = 1291*pow(omhh,0.251)/(1+0.659*pow((double)omhh,0.828))*
		(1+z_drag_b1*pow((double)obhh,(double)z_drag_b2));
		
		R_drag = 31.5*obhh/POW4(theta_cmb)*(1000/(1+z_drag));
		R_equality = 31.5*obhh/POW4(theta_cmb)*(1000/z_equality);
		
		sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
	    log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));
		
		k_silk = 1.6*pow((double)obhh,0.52)*pow((double)omhh,0.73)*(1+pow((double)10.4*omhh,-0.95));
		
		alpha_c_a1 = pow((double)46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
		alpha_c_a2 = pow((double)12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
		alpha_c = pow(alpha_c_a1,-f_baryon)*
		pow(alpha_c_a2,-CUBE(f_baryon));
		
		beta_c_b1 = 0.944/(1+pow(458*omhh,-0.708));
		beta_c_b2 = pow(0.395*omhh, -0.0266);
		beta_c = 1.0/(1+beta_c_b1*(pow(1-f_baryon, beta_c_b2)-1));
		
		y = z_equality/(1+z_drag);
		alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
		alpha_b = 2.07*k_equality*sound_horizon*pow(1+R_drag,-0.75)*alpha_b_G;
		
		beta_node = 8.41*pow(omhh, 0.435);
		beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);
		
		k_peak = 2.5*3.14159*(1+0.217*omhh)/sound_horizon;
		sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1+10.0*pow(obhh,0.75));
		
		alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*
		SQR(f_baryon);
		
		return;
	}
	
	//! private member function: computes transfer function for mode k (k in Mpc)
	inline double TFfit_onek(double k, double *tf_baryon, double *tf_cdm)
	/* Input: k -- Wavenumber at which to calculate transfer function, in Mpc^-1.
	 *tf_baryon, *tf_cdm -- Input value not used; replaced on output if
	 the input was not NULL. */
	/* Output: Returns the value of the full transfer function fitting formula.
     This is the form given in Section 3 of Eisenstein & Hu (1997).
     *tf_baryon -- The baryonic contribution to the full fit.
     *tf_cdm -- The CDM contribution to the full fit. */
	/* Notes: Units are Mpc, not h^-1 Mpc. */
	{
		double T_c_ln_beta, T_c_ln_nobeta, T_c_C_alpha, T_c_C_noalpha;
		double q, xx, xx_tilde;//, q_eff;
		double T_c_f, T_c, s_tilde, T_b_T0, T_b, f_baryon, T_full;
		//double T_0_L0, T_0_C0, T_0, gamma_eff; 
		//double T_nowiggles_L0, T_nowiggles_C0, T_nowiggles;
		
		k = fabs(k);	/* Just define negative k as positive */
		if (k==0.0) {
			if (tf_baryon!=NULL) *tf_baryon = 1.0;
			if (tf_cdm!=NULL) *tf_cdm = 1.0;
			return 1.0;
		}
		
		q = k/13.41/k_equality;
		xx = k*sound_horizon;
		
		T_c_ln_beta = log(2.718282+1.8*beta_c*q);
		T_c_ln_nobeta = log(2.718282+1.8*q);
		T_c_C_alpha = 14.2/alpha_c + 386.0/(1+69.9*pow(q,1.08));
		T_c_C_noalpha = 14.2 + 386.0/(1+69.9*pow(q,1.08));
		
		T_c_f = 1.0/(1.0+POW4(xx/5.4));
		T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*SQR(q)) +
	    (1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*SQR(q));
		
		s_tilde = sound_horizon*pow(1.+CUBE(beta_node/xx),-1./3.);
		xx_tilde = k*s_tilde;
		
		T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*SQR(q));
		T_b = sin(xx_tilde)/(xx_tilde)*(T_b_T0/(1.+SQR(xx/5.2))+
										alpha_b/(1.+CUBE(beta_b/xx))*exp(-pow(k/k_silk,1.4)));
		
		f_baryon = obhh/omhh;
		T_full = f_baryon*T_b + (1-f_baryon)*T_c;
		
		/* Now to store these transfer functions */
		if (tf_baryon!=NULL) *tf_baryon = T_b;
		if (tf_cdm!=NULL) *tf_cdm = T_c;
		return T_full;
	}
	
public:
	//! Constructor for Eisenstein & Hu fitting for transfer function
	/*!
	 \param aCosm structure of type Cosmology carrying the cosmological parameters
	 \param Tcmb mean temperature of the CMB fluctuations (defaults to
	 Tcmb = 2.726 if not specified)
	 */
	transfer_eisenstein_plugin( config_file &cf )//Cosmology aCosm, double Tcmb = 2.726 )
    :  transfer_function_plugin(cf), m_h0( cosmo_.H0*0.01 )
	{
		double Tcmb = pcf_->getValueSafe("cosmology","Tcmb",2.726);
		TFset_parameters( (cosmo_.Omega_m)*cosmo_.H0*cosmo_.H0*(0.01*0.01), 
						 cosmo_.Omega_b/(cosmo_.Omega_m-cosmo_.Omega_b),//-aCosm.Omega_b), 
						 Tcmb);
		
		tf_distinct_ = false;
		tf_withvel_  = false;
	}
	
	//! Computes the transfer function for k in Mpc/h by calling TFfit_onek
	inline double compute( double k, tf_type type ){
		double tfb, tfcdm, fb, fc; 
		TFfit_onek( k*m_h0, &tfb, &tfcdm );
		
		fb = cosmo_.Omega_b/(cosmo_.Omega_m);
		fc = (cosmo_.Omega_m-cosmo_.Omega_b)/(cosmo_.Omega_m) ;
		
		return fb*tfb+fc*tfcdm;
	}
	
	
	inline double get_kmin( void ){
		return 1e-4;
	}
	
	inline double get_kmax( void ){
		return 1.e4;
	}
	
};



class transfer_eisenstein_wdm_plugin : public transfer_eisenstein_plugin
{
protected:
	using transfer_eisenstein_plugin::TFfit_onek;
	real_t m_WDMalpha; 
	double omegam_, wdmm_, wdmgx_, H0_, omegab_;
public:
	transfer_eisenstein_wdm_plugin( config_file &cf )
	: transfer_eisenstein_plugin( cf )
	{
		omegam_ = cf.getValue<double>("cosmology","Omega_m");
		omegab_ = cf.getValue<double>("cosmology","Omega_b");
		wdmm_   = cf.getValue<double>("cosmology","WDMmass");
		wdmgx_  = cf.getValue<double>("cosmology","WDMg_x");
		H0_     = cf.getValue<double>("cosmology","H0");
		m_WDMalpha = 0.05 * pow( omegam_/0.4,0.15)
					*pow(H0_/0.65,1.3)*pow(wdmm_,-1.15)
					*pow(1.5/wdmgx_,0.29);
	}
	
	inline real_t compute( real_t k, tf_type type )
	{
		double tfb, tfcdm, fb, fc;
		TFfit_onek( k*m_h0, &tfb, &tfcdm );
		
		fb = omegab_/omegam_;
		fc = (omegam_-omegab_)/omegam_ ;
		
		double tf = fb*tfb+fc*tfcdm;
		
		return tf*pow(1.0+(m_WDMalpha*k)*(m_WDMalpha*k),-5.0);
		
	}
	
};

#if 0

class TransferFunction_EisensteinNeutrino : public TransferFunction
{
	/* Fitting Formulae for CDM + Baryon + Massive Neutrino (MDM) cosmologies. */
	/* Daniel J. Eisenstein & Wayne Hu, Institute for Advanced Study */
	
	/* There are two primary routines here, one to set the cosmology, the
	 other to construct the transfer function for a single wavenumber k. 
	 You should call the former once (per cosmology) and the latter as 
	 many times as you want. */
	
	/* TFmdm_set_cosm() -- User passes all the cosmological parameters as
	 arguments; the routine sets up all of the scalar quantites needed 
	 computation of the fitting formula.  The input parameters are: 
	 1) omega_matter -- Density of CDM, baryons, and massive neutrinos,
	 in units of the critical density. 
	 2) omega_baryon -- Density of baryons, in units of critical. 
	 3) omega_hdm    -- Density of massive neutrinos, in units of critical 
	 4) degen_hdm    -- (Int) Number of degenerate massive neutrino species 
	 5) omega_lambda -- Cosmological constant 
	 6) hubble       -- Hubble constant, in units of 100 km/s/Mpc 
	 7) redshift     -- The redshift at which to evaluate */
	
	/* TFmdm_onek_mpc() -- User passes a single wavenumber, in units of Mpc^-1.
	 Routine returns the transfer function from the Eisenstein & Hu
	 fitting formula, based on the cosmology currently held in the 
	 internal variables.  The routine returns T_cb (the CDM+Baryon
	 density-weighted transfer function), although T_cbn (the CDM+
	 Baryon+Neutrino density-weighted transfer function) is stored
	 in the global variable tf_cbnu. */
	
	/* We also supply TFmdm_onek_hmpc(), which is identical to the previous
	 routine, but takes the wavenumber in units of h Mpc^-1. */
	
	/* We hold the internal scalar quantities in global variables, so that
	 the user may access them in an external program, via "extern" declarations. */
	
	/* Please note that all internal length scales are in Mpc, not h^-1 Mpc! */
	
	/* -------------------------- Prototypes ----------------------------- */
	
	/*	int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm,
	 int degen_hdm, float omega_lambda, float hubble, float redshift);
	 float TFmdm_onek_mpc(float kk);
	 float TFmdm_onek_hmpc(float kk);
	 */	
	
	/* ------------------------- Global Variables ------------------------ */
	
	/* The following are set in TFmdm_set_cosm() */
	float   alpha_gamma,	/* sqrt(alpha_nu) */
	alpha_nu,	/* The small-scale suppression */
	beta_c,		/* The correction to the log in the small-scale */
	num_degen_hdm,	/* Number of degenerate massive neutrino species */
	f_baryon,	/* Baryon fraction */
	f_bnu,		/* Baryon + Massive Neutrino fraction */
	f_cb,		/* Baryon + CDM fraction */
	f_cdm,		/* CDM fraction */
	f_hdm,		/* Massive Neutrino fraction */
	growth_k0,	/* D_1(z) -- the growth function as k->0 */
	growth_to_z0,	/* D_1(z)/D_1(0) -- the growth relative to z=0 */
	hhubble,	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
	k_equality,	/* The comoving wave number of the horizon at equality*/
	obhh,		/* Omega_baryon * hubble^2 */
	omega_curv,	/* = 1 - omega_matter - omega_lambda */
	omega_lambda_z, /* Omega_lambda at the given redshift */
	omega_matter_z,	/* Omega_matter at the given redshift */
	omhh,		/* Omega_matter * hubble^2 */
	onhh,		/* Omega_hdm * hubble^2 */
	p_c,		/* The correction to the exponent before drag epoch */
	p_cb,		/* The correction to the exponent after drag epoch */
	sound_horizon_fit,  /* The sound horizon at the drag epoch */
	theta_cmb,	/* The temperature of the CMB, in units of 2.7 K */
	y_drag,		/* Ratio of z_equality to z_drag */
	z_drag,		/* Redshift of the drag epoch */
	z_equality;	/* Redshift of matter-radiation equality */
	
	/* The following are set in TFmdm_onek_mpc() */
	float	gamma_eff,	/* Effective \Gamma */
	growth_cb,	/* Growth factor for CDM+Baryon perturbations */
	growth_cbnu,	/* Growth factor for CDM+Baryon+Neutrino pert. */
	max_fs_correction,  /* Correction near maximal free streaming */
	qq,		/* Wavenumber rescaled by \Gamma */
	qq_eff,		/* Wavenumber rescaled by effective Gamma */
	qq_nu,		/* Wavenumber compared to maximal free streaming */
	tf_master,	/* Master TF */
	tf_sup,		/* Suppressed TF */
	y_freestream; 	/* The epoch of free-streaming for a given scale */
	
	/* Finally, TFmdm_onek_mpc() and TFmdm_onek_hmpc() give their answers as */
	float   tf_cb,		/* The transfer function for density-weighted
						 CDM + Baryon perturbations. */
	tf_cbnu;	/* The transfer function for density-weighted
				 CDM + Baryon + Massive Neutrino perturbations. */
	
	/* By default, these functions return tf_cb */
	
	/* ------------------------- TFmdm_set_cosm() ------------------------ */
	int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm,
					   int degen_hdm, float omega_lambda, float hubble, float redshift)
	/* This routine takes cosmological parameters and a redshift and sets up
	 all the internal scalar quantities needed to compute the transfer function. */
	/* INPUT: omega_matter -- Density of CDM, baryons, and massive neutrinos,
	 in units of the critical density. */
	/* 	  omega_baryon -- Density of baryons, in units of critical. */
	/* 	  omega_hdm    -- Density of massive neutrinos, in units of critical */
	/* 	  degen_hdm    -- (Int) Number of degenerate massive neutrino species */
	/*        omega_lambda -- Cosmological constant */
	/* 	  hubble       -- Hubble constant, in units of 100 km/s/Mpc */
	/*        redshift     -- The redshift at which to evaluate */
	/* OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
	 sets many global variables for use in TFmdm_onek_mpc() */
	{
		float z_drag_b1, z_drag_b2, omega_denom;
		int qwarn;
		qwarn = 0;
		
		theta_cmb = 2.728/2.7;	/* Assuming T_cmb = 2.728 K */
		
		/* Look for strange input */
		if (omega_baryon<0.0) {
			fprintf(stderr,
					"TFmdm_set_cosm(): Negative omega_baryon set to trace amount.\n");
			qwarn = 1;
		}
		if (omega_hdm<0.0) {
			fprintf(stderr,
					"TFmdm_set_cosm(): Negative omega_hdm set to trace amount.\n");
			qwarn = 1;
		}
		if (hubble<=0.0) {
			fprintf(stderr,"TFmdm_set_cosm(): Negative Hubble constant illegal.\n");
			exit(1);  /* Can't recover */
		} else if (hubble>2.0) {
			fprintf(stderr,"TFmdm_set_cosm(): Hubble constant should be in units of 100 km/s/Mpc.\n");
			qwarn = 1;
		}
		if (redshift<=-1.0) {
			fprintf(stderr,"TFmdm_set_cosm(): Redshift < -1 is illegal.\n");
			exit(1);
		} else if (redshift>99.0) {
			fprintf(stderr,
					"TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate.\n");
			qwarn = 1;
		}
		if (degen_hdm<1) degen_hdm=1;
		num_degen_hdm = (float) degen_hdm;	
		/* Have to save this for TFmdm_onek_mpc() */
		/* This routine would crash if baryons or neutrinos were zero, 
		 so don't allow that */
		if (omega_baryon<=0) omega_baryon=1e-5;
		if (omega_hdm<=0) omega_hdm=1e-5;
		
		omega_curv = 1.0-omega_matter-omega_lambda;
		omhh = omega_matter*SQR(hubble);
		obhh = omega_baryon*SQR(hubble);
		onhh = omega_hdm*SQR(hubble);
		f_baryon = omega_baryon/omega_matter;
		f_hdm = omega_hdm/omega_matter;
		f_cdm = 1.0-f_baryon-f_hdm;
		f_cb = f_cdm+f_baryon;
		f_bnu = f_baryon+f_hdm;
		
		/* Compute the equality scale. */
		z_equality = 25000.0*omhh/SQR(SQR(theta_cmb));	/* Actually 1+z_eq */
		k_equality = 0.0746*omhh/SQR(theta_cmb);
		
		/* Compute the drag epoch and sound horizon */
		z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
		z_drag_b2 = 0.238*pow(omhh,0.223);
		z_drag = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*
		(1.0+z_drag_b1*pow(obhh,z_drag_b2));
		y_drag = z_equality/(1.0+z_drag);
		
		sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));
		
		/* Set up for the free-streaming & infall growth function */
		p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
		p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));
		
		omega_denom = omega_lambda+SQR(1.0+redshift)*(omega_curv+
													  omega_matter*(1.0+redshift));
		omega_lambda_z = omega_lambda/omega_denom;
		omega_matter_z = omega_matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
		growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/
	    (pow(omega_matter_z,4.0/7.0)-omega_lambda_z+
		 (1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
		growth_to_z0 = z_equality*2.5*omega_matter/(pow(omega_matter,4.0/7.0)
													-omega_lambda + (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0));
		growth_to_z0 = growth_k0/growth_to_z0;	
		
		/* Compute small-scale suppression */
		alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*
		pow(1+y_drag,p_cb-p_c)*
		(1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/
		(1-0.193*sqrt(f_hdm*num_degen_hdm)+0.169*f_hdm*pow(num_degen_hdm,0.2))*
		(1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
		alpha_gamma = sqrt(alpha_nu);
		beta_c = 1/(1-0.949*f_bnu);
		/* Done setting scalar variables */
		hhubble = hubble;	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
		return qwarn;
	}
	
	/* ---------------------------- TFmdm_onek_mpc() ---------------------- */
	
	float TFmdm_onek_mpc(float kk)
	/* Given a wavenumber in Mpc^-1, return the transfer function for the
	 cosmology held in the global variables. */
	/* Input: kk -- Wavenumber in Mpc^-1 */
	/* Output: The following are set as global variables:
	 growth_cb -- the transfer function for density-weighted
	 CDM + Baryon perturbations. 
	 growth_cbnu -- the transfer function for density-weighted
	 CDM + Baryon + Massive Neutrino perturbations. */
	/* The function returns growth_cb */
	{
		float tf_sup_L, tf_sup_C;
		float temp1, temp2;
		
		qq = kk/omhh*SQR(theta_cmb);
		
		/* Compute the scale-dependent growth functions */
		y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*
		SQR(num_degen_hdm*qq/f_hdm);
		temp1 = pow(growth_k0, 1.0-p_cb);
		temp2 = pow(growth_k0/(1+y_freestream),0.7);
		growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
		growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;
		
		/* Compute the master function */
		gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/
						 (1+SQR(SQR(kk*sound_horizon_fit*0.43))));
		qq_eff = qq*omhh/gamma_eff;
		
		tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
		tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
		tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff));
		
		qq_nu = 3.92*qq*sqrt(num_degen_hdm/f_hdm);
		max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(num_degen_hdm,0.3+0.6*f_hdm)/
		(pow(qq_nu,-1.6)+pow(qq_nu,0.8));
		tf_master = tf_sup*max_fs_correction;
		
		/* Now compute the CDM+HDM+baryon transfer functions */
		tf_cb = tf_master*growth_cb/growth_k0;
		tf_cbnu = tf_master*growth_cbnu/growth_k0;
		return tf_cb;
	}
	
	/* ---------------------------- TFmdm_onek_hmpc() ---------------------- */
	
	float TFmdm_onek_hmpc(float kk)
	/* Given a wavenumber in h Mpc^-1, return the transfer function for the
	 cosmology held in the global variables. */
	/* Input: kk -- Wavenumber in h Mpc^-1 */
	/* Output: The following are set as global variables:
	 growth_cb -- the transfer function for density-weighted
	 CDM + Baryon perturbations. 
	 growth_cbnu -- the transfer function for density-weighted
	 CDM + Baryon + Massive Neutrino perturbations. */
	/* The function returns growth_cb */
	{
		return TFmdm_onek_mpc(kk*hhubble);
	}
	
	
	
	TransferFunction_EisensteinNeutrino( Cosmology aCosm, real_t Omega_HDM, int degen_HDM, real_t Tcmb = 2.726 )
    :  TransferFunction(aCosm)//, m_h0( aCosm.H0*0.01 )
	{
		TFmdm_set_cosm( aCosm.Omega_m, aCosm.Omega_b, Omega_HDM, degen_HDM, aCosm.Omega_L, aCosm.H0, aCosm.astart );
		
	}
	
	//! Computes the transfer function for k in Mpc/h by calling TFfit_onek
	virtual inline real_t compute( real_t k ){
		return TFmdm_onek_hmpc( k );
		
		
		/*real_t tfb, tfcdm, fb, fc; //, tfull
		 TFfit_onek( k*m_h0, &tfb, &tfcdm );
		 
		 fb = m_Cosmology.Omega_b/(m_Cosmology.Omega_m);
		 fc = (m_Cosmology.Omega_m-m_Cosmology.Omega_b)/(m_Cosmology.Omega_m) ;
		 
		 return fb*tfb+fc*tfcdm;*/
		//return 1.0;
	}
	
	
};

#endif

namespace{
	transfer_function_plugin_creator_concrete< transfer_eisenstein_plugin > creator("eisenstein");
	transfer_function_plugin_creator_concrete< transfer_eisenstein_wdm_plugin > creator2("eisenstein_wdm");
}

