c=====================================================================================c
c        
c The code below was written by: Adrian Jenkins,                                            
c                                Institute for Computational Cosmology
c                                Department of Physics
c                                South Road
c                                Durham, DH1 3LE
c                                United Kingdom
c
c This file is part of the software made public in
c Jenkins and Booth 2013  - arXiv:1306.XXXX
c
c The software computes the Panphasia Gaussian white noise field
c realisation described in detail in Jenkins 2013 - arXiv:1306.XXXX
c 
c
c
c This software is free, subject to a agreeing licence conditions:
c
c
c (i)  you will publish the phase descriptors and reference Jenkins (13) 
c      for any new simulations that use Panphasia phases. You will pass on this 
c      condition to others for any software or data you make available publically 
c      or privately that makes use of Panphasia. 
c
c (ii) that you will ensure any publications using results derived from Panphasia 
c      will be submitted as a final version to arXiv prior to or coincident with
c      publication in a journal. 
c
c (iii) that you report any bugs in this software as soon as confirmed to 
c       A.R.Jenkins@durham.ac.uk 
c
c (iv)  that you understand that this software comes with no warranty and that is 
c       your responsibility to ensure that it is suitable for the purpose that 
c       you intend. 
c
c=====================================================================================c

c=====================================================================================
c       List of subroutines and arguments.  Each of these is documented in           c
c       arXiV/1306.XXXX                                                              c
c                                                                                    c
c       Adrian Jenkins, 24/6/2013.                                                   c
c-------------------------------------------------------------------------------------
c  Version 1.000
c===================================================================================

      module pan_state
      use Rand
      implicit none
      integer maxdim_, maxlev_, maxpow_
      parameter (maxdim_=60,maxlev_=50, maxpow_ = 3*maxdim_)
      integer nmulti_
      parameter (nmulti_=64)
      integer range_max
      parameter(range_max=10000)
      integer indmin,indmax
      parameter (indmin=-1, indmax=60)


      type state_data
      integer base_state(5), base_lev_start(5,0:maxdim_)
      TYPE(Rand_offset) :: poweroffset(0:maxpow_)
      TYPE(Rand_offset) :: superjump
      TYPE(Rand_state) :: current_state(-1:maxpow_)

      integer  layer_min,layer_max,indep_field

!  This module stores information needed to access the part of Panphasia
!  selected by a particular descriptor.
      integer*8 xorigin_store(0:1,0:1,0:1)
      integer*8 yorigin_store(0:1,0:1,0:1)
      integer*8 zorigin_store(0:1,0:1,0:1)

      integer*4 lev_common
      integer*4 layer_min_store,layer_max_store

      integer*8 ix_abs_store,iy_abs_store,iz_abs_store    
      integer*8 ix_per_store,iy_per_store,iz_per_store
      integer*8 ix_rel_store,iy_rel_store,iz_rel_store

      real*8 exp_coeffs(8,0:7,-1:maxdim_)
      integer*8 xcursor(0:maxdim_),ycursor(0:maxdim_),zcursor(0:maxdim_)

c    Local box parameters

      integer*4 ixshift(0:1,0:1,0:1)
      integer*4 iyshift(0:1,0:1,0:1)
      integer*4 izshift(0:1,0:1,0:1)


c     more state variables
      real*8 cell_data(9,0:7)
      integer*4 ixh_last,iyh_last,izh_last
      integer init

      integer return_cell_props_init
      integer reset_lecuyer_state_init
      integer*8 p_xcursor(indmin:indmax),p_ycursor(indmin:indmax),p_zcursor(indmin:indmax)



      end type state_data



c     Switch for enabling custom spherical function
c     Set isub_spherical_function = 1 to turn on the spherical function
      integer*4 isub_spherical_function
      parameter (isub_spherical_function=0)

      end module pan_state


c================================================================================
c       Begin white noise routines
c================================================================================
      recursive subroutine start_panphasia(ldata,descriptor,ngrid,VERBOSE)
      use pan_state
      implicit none
      type(state_data), intent(inout) :: ldata
      character*100 descriptor
      integer ngrid
      integer VERBOSE

      

      integer*4 wn_level_base,i_base,i_base_y,i_base_z
      integer*8 i_xorigin_base,i_yorigin_base,i_zorigin_base, check_rand
      character*20 name

      integer ratio
      integer lextra
      integer level_p

      
      integer*8 ix_abs,iy_abs,iz_abs
      integer*8 ix_per,iy_per,iz_per
      integer*8 ix_rel,iy_rel,iz_rel
      
      !integer  layer_min,layer_max,indep_field
      !common /oct_range/  layer_min,layer_max,indep_field

      call parse_descriptor(descriptor ,wn_level_base,i_xorigin_base,i_yorigin_base,
     &                      i_zorigin_base,i_base,i_base_y,i_base_z,check_rand,name)


      lextra = (log10(real(ngrid)/real(i_base))+0.001)/log10(2.0)
      ratio = 2**lextra

      if (ratio*i_base.ne.ngrid) 
     &stop 'Value of ngrid inconsistent with dim of region in Panphasia'

      level_p = wn_level_base + lextra

      ix_abs = ishft(i_xorigin_base,lextra)
      iy_abs = ishft(i_yorigin_base,lextra)
      iz_abs = ishft(i_zorigin_base,lextra)

      ix_per = i_base*ratio
      iy_per = i_base*ratio
      iz_per = i_base*ratio

c     Set the refinement position at the origin. 
   
      ix_rel = 0
      iy_rel = 0
      iz_rel = 0

      call set_phases_and_rel_origin(ldata,descriptor,level_p,ix_rel,iy_rel,iz_rel,VERBOSE)

c    Finally set the octree functions required for making cosmological
c    initial conditions.  These are passed using a common block.

      ldata%layer_min = 0
      ldata%layer_max = level_p
      ldata%indep_field  = 1

      end
c=================================================================================
      recursive subroutine set_phases_and_rel_origin(ldata,descriptor,lev,ix_rel,iy_rel,iz_rel,VERBOSE)
      use pan_state
      !use descriptor_phases
      implicit none
      type(state_data), intent(inout) :: ldata
      character*100 descriptor
      integer lev
      integer*8 ix_abs,iy_abs,iz_abs
      integer*8 ix_per,iy_per,iz_per
      integer*8 ix_rel,iy_rel,iz_rel
      integer*8 xorigin,yorigin,zorigin

      integer VERBOSE
      integer MYID
      integer*8 maxco
      integer i
      integer px,py,pz
      
      integer lnblnk
      integer*8 mconst
      parameter(mconst = 2147483647_Dint)

      integer*4 wn_level_base,i_base,i_base_y,i_base_z
      integer*8 i_xorigin_base,i_yorigin_base,i_zorigin_base, check_rand
      integer lextra,ratio
      character*20 phase_name

c-----------------------------------------------------------------------------------------------

      call initialise_panphasia(ldata)

      call validate_descriptor(ldata, descriptor,-1,check_rand)

      call parse_descriptor(descriptor ,wn_level_base,i_xorigin_base,i_yorigin_base,
     &                      i_zorigin_base,i_base,i_base_y,i_base_z,check_rand,phase_name)
      lextra = lev - wn_level_base
      ratio  = 2**lextra
      
      ix_abs = ishft(i_xorigin_base,lextra)
      iy_abs = ishft(i_yorigin_base,lextra)
      iz_abs = ishft(i_zorigin_base,lextra)

      ix_per = i_base*ratio
      iy_per = i_base*ratio
      iz_per = i_base*ratio

c-------------------------------------------------------------------------
c    Error checking
c-------------------------------------------------------------------------
      if ((lev.lt.0).or.(lev.gt.maxlev_)) stop 'Level out of range! (1)'
      

      maxco = 2_dint**lev

      if (ix_abs.lt.0) stop 'Error: ix_abs negative (1)'
      if (iy_abs.lt.0) stop 'Error: iy_abs negative (1)'
      if (iz_abs.lt.0) stop 'Error: iz_abs negative (1)'

      if (ix_rel.lt.0) stop 'Error: ix_rel negative (1)'
      if (iy_rel.lt.0) stop 'Error: iy_rel negative (1)'
      if (iz_rel.lt.0) stop 'Error: iz_rel negative (1)'


      if (ix_abs+ix_rel.ge.maxco)
     &   stop 'Error: ix_abs + ix_rel out of range. (1)'
      if (iy_abs+iy_rel.ge.maxco) 
     &   stop 'Error: iy_abs + iy_rel out of range. (1)'
      if (iz_abs+iz_rel.ge.maxco) 
     &   stop 'Error: iz_abs + iz_rel out of range. (1)'

c----------------------------------------------------------------------------------------
c  To allow the local box to wrap around, if needed, define a series of eight
c  'origins'.  For many purposes (ix,iy,iz) = (0,0,0) is the only origin needed.


      do px=0,1
       do py=0,1
        do pz=0,1

         xorigin = max(0,( ix_abs + ix_rel - px*ix_per )/2)
         yorigin = max(0,( iy_abs + iy_rel - py*iy_per )/2)
         zorigin = max(0,( iz_abs + iz_rel - pz*iz_per )/2)

         ldata%ixshift(px,py,pz) = max(0, ix_abs + ix_rel -px*ix_per) - 2*xorigin
         ldata%iyshift(px,py,pz) = max(0, iy_abs + iy_rel -py*iy_per) - 2*yorigin
         ldata%izshift(px,py,pz) = max(0, iz_abs + iz_rel -pz*iz_per) - 2*zorigin


c        Store box details:  store the positions at level lev-1
  

         ldata%xorigin_store(px,py,pz) = xorigin
         ldata%yorigin_store(px,py,pz) = yorigin
         ldata%zorigin_store(px,py,pz) = zorigin

        enddo
       enddo
      enddo

      ldata%lev_common = lev


      ldata%ix_abs_store = ix_abs
      ldata%iy_abs_store = iy_abs
      ldata%iz_abs_store = iz_abs

      ldata%ix_per_store = ix_per
      ldata%iy_per_store = iy_per
      ldata%iz_per_store = iz_per

      ldata%ix_rel_store = ix_rel
      ldata%iy_rel_store = iy_rel
      ldata%iz_rel_store = iz_rel

 
c  Reset all cursor values to negative numbers.

      do i=0,maxdim_
       ldata%xcursor(i) = -999
       ldata%ycursor(i) = -999
       ldata%zcursor(i) = -999
      enddo
      if (VERBOSE.gt.1) then
         if (MYID.lt.1) then
            print*,'----------------------------------------------------------'
            print*,'Successfully initialised Panphasia box at level ',lev
            write (6,105) ix_abs,iy_abs,iz_abs
            write (6,106) ix_rel,iy_rel,iz_rel
            write (6,107) ix_per,iy_per,iz_per
            write (6,*)  'Phases used: ',descriptor(1:lnblnk(descriptor))
            print*,'----------------------------------------------------------'
         endif
      endif
 105  format(' Abs origin: (',i12,',',i12,',',i12,')')
 106  format(' Rel origin: (',i12,',',i12,',',i12,')')
 107  format(' Periods   : (',i12,',',i12,',',i12,')') 
      end 
c================================================================================
      recursive subroutine initialise_panphasia( ldata )
      use Rand
      use pan_state
      implicit none

      type(state_data), intent(inout) :: ldata

      TYPE(Rand_state) :: state
      TYPE(Rand_offset) :: offset
      integer ninitialise
      parameter (ninitialise=218)
      integer i
      real*8 rand_num


      call Rand_seed(state,ninitialise)
      
      call Rand_save(ldata%base_state,state)

      call Rand_set_offset(offset,1)

c   Calculate offsets of powers of 2 times nmulti
c

      do i=0,maxpow_
        ldata%poweroffset(i) = Rand_mul_offset(offset,nmulti_)
        offset = Rand_mul_offset(offset,2)
      enddo


c   Compute the base state for each level. 

      call Rand_load(state,ldata%base_state)
      state = Rand_step(state,8)

      do i=0,maxdim_
       call Rand_save(ldata%base_lev_start(1,i),state)
       state = Rand_boost(state,ldata%poweroffset(3*i))
      enddo

c   Set superjump to value 2**137   - used occasionally in computing Gaussian variables
c   when the value of the returned random number is less an 10-6.

       call Rand_set_offset(ldata%superjump,1)

       do i=1,137
         ldata%superjump = Rand_mul_offset(ldata%superjump,2)
       enddo  


c   Run time test to see if one particular value can be recovered.
      
      call Rand_load(state,ldata%base_lev_start(1,34))
      call Rand_real(rand_num,state)

      if (abs(rand_num- 0.828481889948473d0).gt.1.e-14) then
        print*,'Error in initialisation!'
        print*,'Rand_num     = ',rand_num
        print*,'Target value = ', 0.828481889948473d0
        stop
      endif
      return
      end
c=================================================================================
      recursive subroutine panphasia_cell_properties(ldata,ixcell,iycell,izcell,cell_prop)
      use pan_state
      implicit none
      type(state_data), intent(inout) :: ldata
      !integer  layer_min,layer_max,indep_field
      !common /oct_range/  layer_min,layer_max,indep_field
      integer*4 ixcell,iycell,izcell
      real*8 cell_prop(9)

      call adv_panphasia_cell_properties(ldata,ixcell,iycell,izcell,ldata%layer_min,
     &                                           ldata%layer_max,ldata%indep_field,cell_prop)
      return
      end
c=================================================================================
      recursive subroutine adv_panphasia_cell_properties(ldata,ixcell,iycell,izcell,layer_min,
     &                                           layer_max,indep_field,cell_prop)
      use pan_state
      !use descriptor_phases
      implicit none

      type(state_data), intent(inout) :: ldata

      integer*4 lev
      integer*4 ixcell,iycell,izcell
      integer layer_min,layer_max,indep_field
      real*8 cell_prop(9)
c      real*8 cell_data(9,0:7)
      integer*4 j,l,lx,ly,lz
      integer*4 px,py,pz

c      integer*4 ixh_last,iyh_last,izh_last

c      integer init
c      data init/0/
c      save init,cell_data,ixh_last,iyh_last,izh_last  ! Keep internal state

      integer*4 ixh,iyh,izh

      lev = ldata%lev_common

c------- Error checking -----------------------------

       if (layer_min.gt.layer_max) then

              if (layer_min-layer_max.eq.1) then      ! Not necessarily bad. No octree basis functions
                   do j=1,9                           ! required at this level and position.
                   cell_prop(j) = 0.0d0               ! Set returned cell_prop data to zero.
                   enddo
                   return
              endif

              print*,'Warning: layer_min.gt.layer_max!'
              print*,'layer_min = ',layer_min
              print*,'layer_max = ',layer_max
              print*,'ixcell,iycell,izcell',ixcell,iycell,izcell

              call flush(6)
              stop 'Error: layer_min.gt.layer_max'
       endif

       if (layer_max.gt.ldata%lev_common) then
          print*,'lev_common = ',ldata%lev_common
          print*,'layer_min  = ',layer_min
          print*,'layer_max  = ',layer_max
          stop 'Error: layer_max.gt.lev_common'
       endif
       if ((indep_field.lt.-1).or.(indep_field.gt.1)) 
     & stop 'Error: indep_field out of range'

c----------------------------------------------------
c  Check which 'origin' to use.  

      px = 0
      py = 0
      pz = 0

      if (ldata%ix_rel_store+ixcell.ge.ldata%ix_per_store) px = 1  ! Crossed x-periodic bndy
      if (ldata%iy_rel_store+iycell.ge.ldata%iy_per_store) py = 1  ! Crossed y-periodic bndy
      if (ldata%iz_rel_store+izcell.ge.ldata%iz_per_store) pz = 1  ! Crossed z-periodic bndy
c----------------------------------------------------


      ixh = (ixcell+ldata%ixshift(px,py,pz) )/2
      iyh = (iycell+ldata%iyshift(px,py,pz) )/2
      izh = (izcell+ldata%izshift(px,py,pz) )/2

      lx  = mod(ixcell+ldata%ixshift(px,py,pz) ,2)
      ly  = mod(iycell+ldata%iyshift(px,py,pz) ,2)
      lz  = mod(izcell+ldata%izshift(px,py,pz) ,2)


      l = 4*lx + 2*ly + lz   ! Determine which cell is required

cc------------------   If no new evalation is needed skip assignment -----
      if ((ldata%init.eq.1).and.(ixh.eq.ldata%ixh_last).and.(iyh.eq.ldata%iyh_last).and.
     &   (izh.eq.ldata%izh_last).and.(layer_min.eq.ldata%layer_min_store).and.
     &   (layer_max.eq.ldata%layer_max_store)) goto 24
cc-----------------------------------------------------------------------------


       call   return_cell_props(ldata,lev,ixh,iyh,izh,px,py,pz,layer_min,
     &      layer_max,indep_field,ldata%cell_data)

c  Remember previous values.
 
       ldata%ixh_last = ixh
       ldata%iyh_last = iyh
       ldata%izh_last = izh


 24    continue

 
      do j=1,9
       cell_prop(j) = ldata%cell_data(j,l)  ! Copy the required data
      enddo
     
      if (ldata%init.eq.0) ldata%init=1

      return
      end
c=================================================================================
      recursive subroutine return_cell_props(ldata,lev_input,ix_half,iy_half,iz_half,
     &  px,py,pz,layer_min,layer_max,indep_field,cell_data)
      use Rand
      use pan_state
      !use descriptor_phases
      implicit none
      type(state_data), intent(inout) :: ldata
      integer lev_input,ix_half,iy_half,iz_half,px,py,pz
      integer layer_min,layer_max,indep_field 
      real*8 cell_data(9,0:7)

      real*8 garray(0:63)
      integer lev
      integer*8 xarray,yarray,zarray

      integer i,istart,icell_name


c      integer init
c      data init/0/
c      save init

 

c--------------------------------------------------------
c--------------------------- Initialise level -1 --------
c--------------------------------------------------------

      if (ldata%return_cell_props_init.eq.0) then                          ! First time called. Set up the Legendre coefficients    
      ldata%return_cell_props_init = 1                                     ! for the root cell.   This is the first term on the
      call Rand_load(ldata%current_state(-1),ldata%base_state) ! right hand side of the equation in appendix C of
      call return_gaussian_array(ldata,-1,8,garray)      ! Jenkins 2013 that defines PANPHASIA.
      ldata%exp_coeffs(1,0,-1) = garray(0)
      ldata%exp_coeffs(2,0,-1) = garray(1)
      ldata%exp_coeffs(3,0,-1) = garray(2)
      ldata%exp_coeffs(4,0,-1) = garray(3)
      ldata%exp_coeffs(5,0,-1) = garray(4)
      ldata%exp_coeffs(6,0,-1) = garray(5)
      ldata%exp_coeffs(7,0,-1) = garray(6)
      ldata%exp_coeffs(8,0,-1) = garray(7)

      ldata%layer_min_store = layer_min
      ldata%layer_max_store = layer_max
          
      endif

c--------------------------------------------------------
c---------------------------- Error checking ------------
c--------------------------------------------------------

      lev = lev_input-1

      if (lev_input.ne.ldata%lev_common) stop 'Box initialised at a different level !'
      if (ix_half.lt.0) then
          print*,'ix_half negative',ix_half
          stop 'ix_half out of range!'
      endif
      if (iy_half.lt.0) stop 'iy_half out of range!'
      if (iz_half.lt.0) then
          print*,'iz_half negative',iz_half
          stop 'iz_half out of range!' 
      endif


      xarray = ldata%xorigin_store(px,py,pz) + ix_half
      yarray = ldata%yorigin_store(px,py,pz) + iy_half
      zarray = ldata%zorigin_store(px,py,pz) + iz_half


c   If layer_max or layer_min have changed, rebuild from the start and reset the
c   recorded value of layer_max and layer_min

      if ((layer_max.ne.ldata%layer_max_store).or.(layer_min.ne.ldata%layer_min_store)) then

         if (layer_min.gt.layer_max) stop 'layer_min > layer_max : 2'

         istart = max(1,layer_min-1)

         ldata%layer_max_store = layer_max
         ldata%layer_min_store = layer_min

         goto 10

      endif


      if ((xarray.eq.ldata%xcursor(lev)).and.(yarray.eq.ldata%ycursor(lev)).and.(zarray.eq.ldata%zcursor(lev))) return ! Nothing to do.

c===========================================================================================================
c------------- First determine which levels need to be (re)computed
c===========================================================================================================

      istart = 0
      do i=lev-1,0,-1
        if ((ishft(xarray,i-lev).eq.ldata%xcursor(i)).and.(ishft(yarray,i-lev).eq.ldata%ycursor(i)).and.
     &         (ishft(zarray,i-lev).eq.ldata%zcursor(i))) then
            istart = i+1
            goto 10
        endif
      enddo

 10   continue


c====================================================================================
c------------- Now compute each level as required and update (x,y,z) cursor variables
c====================================================================================

      do i=istart,lev

       icell_name = 0

       ldata%xcursor(i) = ishft(xarray,i-lev)
       ldata%ycursor(i) = ishft(yarray,i-lev)
       ldata%zcursor(i) = ishft(zarray,i-lev)

       if (btest(ldata%xcursor(i),0)) icell_name = icell_name + 4
       if (btest(ldata%ycursor(i),0)) icell_name = icell_name + 2
       if (btest(ldata%zcursor(i),0)) icell_name = icell_name + 1

       call reset_lecuyer_state(ldata,i,ldata%xcursor(i),ldata%ycursor(i),ldata%zcursor(i))

       if (isub_spherical_function.ne.1) then
           call return_gaussian_array(ldata,i,64,garray)
       else
           call return_oct_sf_expansion(ldata,i,lev,ldata%xcursor(i),ldata%ycursor(i),ldata%zcursor(i),
     &                                    64,garray)
       endif


       call evaluate_panphasia(ldata,i,maxdim_,garray,layer_min,
     &    layer_max, indep_field, icell_name,cell_data,ldata%exp_coeffs)

      enddo
      return
      end
c=================================================================================
      recursive subroutine  evaluate_panphasia(ldata,nlev,maxdim,g,
     &   layer_min,layer_max,indep_field,icell_name,cell_data,leg_coeff)
      use pan_state
      implicit none
c---------------------------------------------------------------------------------
c    This subroutine calculates the Legendre block coefficients for the eight child
c    cells of an octree cell.
c
c----------------- Define subroutine arguments -----------------------------------
      type(state_data), intent(inout) :: ldata
      integer nlev,maxdim
      integer layer_min,layer_max,indep_field
      integer icell_name
      real*8 leg_coeff(0:7,0:7,-1:maxdim),cell_data(0:8,0:7)
      real*8 g(*)

c----------------- Define constants using notation from appendix A of Jenkins 2013
 
      real*8 a1,a2,b1,b2,b3,c1,c2,c3,c4

      parameter(a1 = 0.5d0*sqrt(3.0d0),      a2 = 0.5d0)

      parameter(b1 = 0.75d0,                 b2 = 0.25d0*sqrt(3.0d0))
      parameter(b3 = 0.25d0)

      parameter(c1 = sqrt(27.0d0/64.0d0),    c2 = 0.375d0)
      parameter(c3 = sqrt(3.0d0/64.0d0),     c4 = 0.125d0)

c----------------- Define octree variables --------------------------------

      real*8 coeff_p000, coeff_p001, coeff_p010, coeff_p011
      real*8 coeff_p100, coeff_p101, coeff_p110, coeff_p111

      real*8 positive_octant_lc(0:7,0:1,0:1,0:1),temp_value(0:7,0:7)
      integer i,j,ix,iy,iz
      integer icx,icy,icz
      integer iox,ioy,ioz
      real*8 parity,isig
      real*8 usually_rooteighth_factor
c--------------------------------------------------------------------------

c-------------  Set the Legendre block coefficients for the parent cell
c               itself. These are either inherited from the octree above
c               or set to zero depending on which levels of the octree
c               have been selected to be populated with the octree
c               basis functions.
c---------------------------------------------------------------------------
      if (nlev.ge.layer_min) then
             coeff_p000  = leg_coeff(0,icell_name,nlev-1)
             coeff_p001  = leg_coeff(1,icell_name,nlev-1)
             coeff_p010  = leg_coeff(2,icell_name,nlev-1)
             coeff_p011  = leg_coeff(3,icell_name,nlev-1)
             coeff_p100  = leg_coeff(4,icell_name,nlev-1)
             coeff_p101  = leg_coeff(5,icell_name,nlev-1)
             coeff_p110  = leg_coeff(6,icell_name,nlev-1)
             coeff_p111  = leg_coeff(7,icell_name,nlev-1)
      else
             coeff_p000  = 0.0d0
             coeff_p001  = 0.0d0
             coeff_p010  = 0.0d0
             coeff_p011  = 0.0d0
             coeff_p100  = 0.0d0
             coeff_p101  = 0.0d0
             coeff_p110  = 0.0d0
             coeff_p111  = 0.0d0 
      endif

c   Apply layer_max and indep_field inputs ---------------------------------

      if (indep_field.ne.-1) then
         usually_rooteighth_factor = sqrt(0.125d0)
      else
         usually_rooteighth_factor = 0.0d0  ! This option returns only the indep field.
      endif                                 ! For use in testing only.

      if (nlev.ge.layer_max) then
        do i=1,56
        g(i) = 0.0d0               ! Set octree coefficients to zero as not required.
        enddo
      endif

      if (indep_field.eq.0) then   ! Set the independent field to zero as not required.
        do i=57,64
        g(i) = 0.0d0
        enddo
      endif
c-----------------------------------------------------------------------------
c
c
c    The calculations immediately below evalute the eight Legendre block coefficients for the
c    child cell that is furthest from the absolute coordiate origin of the octree - we call
c    this the positive octant cell.
c
c    The coefficients are given by a set of matrix equations which combine the
c    coefficients of the Legendre basis functions of the parent cell itself, with
c    the coefficients from the octree basis functions that occupy the
c    parent cell.   
c
c    The Legendre basis function coefficients of the parent cell are stored in
c    the variables, coeff_p000 - coeff_p111 and are initialise above.
c
c    The coefficients of the octree basis functions are determined by the
c    first 56 entries of the array g, which is passed down into this
c    subroutine.
c
c    These two sources of information are combined using a set of linear equations.
c    The coefficients of these linear equations are taken from the inverses or
c    equivalently transposes of the matrices given in appendix A of Jenkins 2013.
c    The matrices in appendix A define the PANPHASIA octree basis functions
c    in terms of Legendre blocks.
c
c    All of the Legendre block functions of the parent cell, and the octree basis
c    functions of the parent cell share one of eight distinct symmetries with respect to
c    reflection about the x1=0,x2=0,x3=0 planes (where the origin is taken as the parent 
c    cell centre and x1,x2,x3 are parallel to the cell edges).
c
c    Each function has either purely reflectional symmetry (even parity) or
c    reflectional symmetry with a sign change (odd parity) about each of the three principal
c    planes through the cell centre. There are therefore 8 parity types. We can label each 
c    parity type with a binary triplet. So 000 is pure reflectional symmetry about 
c    all of the principal planes.
c  
c    In the code below the parent cell Legendre block functions, and octree functions are 
c    organised into eight groups each with eight members. Each group has a common
c    parity type.
c
c    We keep the contributions of each parity type to each of the eight Legendre basis
c    functions occupying the positive octant cell separate. Once they have all been
c    computed, we can apply the different symmetry operations and determine the
c    Legendre block basis functions for all eight child cells at the same time.
c---------------------------------------------------------------------------------------
c    000  parity

      positive_octant_lc(0, 0,0,0) =  1.0d0*coeff_p000
      positive_octant_lc(1, 0,0,0) = -1.0d0*g(1)
      positive_octant_lc(2, 0,0,0) = -1.0d0*g(2)
      positive_octant_lc(3, 0,0,0) =  1.0d0*g(3)
      positive_octant_lc(4, 0,0,0) = -1.0d0*g(4)
      positive_octant_lc(5, 0,0,0) =  1.0d0*g(5)
      positive_octant_lc(6, 0,0,0) =  1.0d0*g(6)
      positive_octant_lc(7, 0,0,0) = -1.0d0*g(7)

c    100 parity

      positive_octant_lc(0, 1,0,0) =  a1*coeff_p100  - a2*g(8)
      positive_octant_lc(1, 1,0,0) =  g(9)
      positive_octant_lc(2, 1,0,0) =  g(10)
      positive_octant_lc(3, 1,0,0) = -g(11)
      positive_octant_lc(4, 1,0,0) =  a2*coeff_p100  + a1*g(8)
      positive_octant_lc(5, 1,0,0) = -g(12) 
      positive_octant_lc(6, 1,0,0) = -g(13)
      positive_octant_lc(7, 1,0,0) =  g(14)

c     010 parity

      positive_octant_lc(0, 0,1,0) =  a1*coeff_p010 - a2*g(15)
      positive_octant_lc(1, 0,1,0) =  g(16) 
      positive_octant_lc(2, 0,1,0) =  a2*coeff_p010 + a1*g(15) 
      positive_octant_lc(3, 0,1,0) = -g(17)
      positive_octant_lc(4, 0,1,0) =  g(18)
      positive_octant_lc(5, 0,1,0) = -g(19)
      positive_octant_lc(6, 0,1,0) = -g(20)
      positive_octant_lc(7, 0,1,0) =  g(21)


c     001 parity

      positive_octant_lc(0, 0,0,1) =  a1*coeff_p001 - a2*g(22)
      positive_octant_lc(1, 0,0,1) =  a2*coeff_p001 + a1*g(22)
      positive_octant_lc(2, 0,0,1) =  g(23)
      positive_octant_lc(3, 0,0,1) = -g(24)
      positive_octant_lc(4, 0,0,1) =  g(25)
      positive_octant_lc(5, 0,0,1) = -g(26)
      positive_octant_lc(6, 0,0,1) = -g(27)
      positive_octant_lc(7, 0,0,1) =  g(28)

c    110 parity

      positive_octant_lc(0, 1,1,0) = b1*coeff_p110 - b2*g(29) + b3*g(30) - b2*g(31)
      positive_octant_lc(1, 1,1,0) = -g(32)
      positive_octant_lc(2, 1,1,0) = b2*coeff_p110 - b3*g(29) - b2*g(30) + b1*g(31)
      positive_octant_lc(3, 1,1,0) =  g(33)
      positive_octant_lc(4, 1,1,0) = b2*coeff_p110 + b1*g(29) + b2*g(30) + b3*g(31)
      positive_octant_lc(5, 1,1,0) =  g(34)
      positive_octant_lc(6, 1,1,0) = b3*coeff_p110 + b2*g(29) - b1*g(30) - b2*g(31)
      positive_octant_lc(7, 1,1,0) = -g(35) 


c     011 parity

      positive_octant_lc(0, 0,1,1) = b1*coeff_p011 - b2*g(36) + b3*g(37) - b2*g(38)
      positive_octant_lc(1, 0,1,1) = b2*coeff_p011 - b3*g(36) - b2*g(37) + b1*g(38)
      positive_octant_lc(2, 0,1,1) = b2*coeff_p011 + b1*g(36) + b2*g(37) + b3*g(38)
      positive_octant_lc(3, 0,1,1) = b3*coeff_p011 + b2*g(36) - b1*g(37) - b2*g(38)
      positive_octant_lc(4, 0,1,1) = -g(39) 
      positive_octant_lc(5, 0,1,1) =  g(40)
      positive_octant_lc(6, 0,1,1) =  g(41)
      positive_octant_lc(7, 0,1,1) = -g(42)

c     101 parity

      positive_octant_lc(0, 1,0,1) = b1*coeff_p101 - b2*g(43) + b3*g(44) - b2*g(45)
      positive_octant_lc(1, 1,0,1) = b2*coeff_p101 - b3*g(43) - b2*g(44) + b1*g(45) 
      positive_octant_lc(2, 1,0,1) = -g(46) 
      positive_octant_lc(3, 1,0,1) =  g(47)
      positive_octant_lc(4, 1,0,1) = b2*coeff_p101 + b1*g(43) + b2*g(44) + b3*g(45)
      positive_octant_lc(5, 1,0,1) = b3*coeff_p101 + b2*g(43) - b1*g(44) - b2*g(45)
      positive_octant_lc(6, 1,0,1) =  g(48)
      positive_octant_lc(7, 1,0,1) = -g(49)

c     111 parity

      positive_octant_lc(0, 1,1,1) = c1*coeff_p111 - c2*g(50) - c2*g(51) - c2*g(52) + c3*g(53) + c3*g(54) + c3*g(55) - c4*g(56)
      positive_octant_lc(1, 1,1,1) = c2*coeff_p111 + c1*g(50) - c2*g(51) + c2*g(52) - c3*g(53) + c3*g(54) + c4*g(55) + c3*g(56) 
      positive_octant_lc(2, 1,1,1) = c2*coeff_p111 + c2*g(50) + c1*g(51) - c2*g(52) - c3*g(53) - c4*g(54) + c3*g(55) - c3*g(56)
      positive_octant_lc(3, 1,1,1) = c3*coeff_p111 - c3*g(50) - c3*g(51) + c4*g(52) - c1*g(53) - c2*g(54) - c2*g(55) - c2*g(56) 
      positive_octant_lc(4, 1,1,1) = c2*coeff_p111 - c2*g(50) + c2*g(51) + c1*g(52) + c4*g(53) - c3*g(54) + c3*g(55) + c3*g(56)
      positive_octant_lc(5, 1,1,1) = c3*coeff_p111 + c3*g(50) - c4*g(51) - c3*g(52) + c2*g(53) - c1*g(54) - c2*g(55) + c2*g(56)
      positive_octant_lc(6, 1,1,1) = c3*coeff_p111 + c4*g(50) + c3*g(51) + c3*g(52) + c2*g(53) + c2*g(54) - c1*g(55) - c2*g(56)
      positive_octant_lc(7, 1,1,1) = c4*coeff_p111 - c3*g(50) + c3*g(51) - c3*g(52) - c2*g(53) + c2*g(54) - c2*g(55) + c1*g(56)
c--------------------------------------------------------------------------------------------
c 
c
c   We now calculate the Legendre basis coefficients for all eight child cells
c   by applying the appropriate reflectional parities to the coefficients 
c   calculated above for the positive octant child cell.
c
c   See equations A2 and A3 in appendix A of Jenkins 2013.
c
c   The reflectional parity is given by (ix,iy,iz) loops below.
c
c   The (icx,icy,icz) loops below, loop over the eight child cells.
c
c   The positive octant child cell is given below by  (icx=icy=icz=0) or i=7.
c
c   The combination ix*icx +iy*icy +iz*icz is either even or odd, depending
c   on whether the parity change is even or odd.
c
c   The variables iox,ioy,ioz are used to loop over the different
c   types of Legendre basis function.
c
c   The combination iox*icx + ioy*icy + ioz*icz is either even and odd
c   and identifies which coefficients keep or change sign respectively
c   due to a pure reflection about the principal planes.
c--------------------------------------------------------------------------------------------

      do iz=0,7
       do iy=0,7
        temp_value(iy,iz) = 0.0d0      ! Zero temporary sums
       enddo
      enddo
c--------------------------------------------------------------------------------------------
      do iz=0,1              ! Loop over z parity (0=keep sign, 1=change sign)
       do iy=0,1             ! Loop over y parity (0=keep sign, 1=change sign)
        do ix=0,1            ! Loop over x parity (0=keep sign, 1=change sign)


         do icx=0,1                      ! Loop over x-child cells
          do icy=0,1                     ! Loop over y-child cells
           do icz=0,1                    ! Loop over z-child cells

             if (mod(ix*icx+iy*icy+iz*icz,2).eq.0) then
                  parity = 1.0d0
             else
                  parity =-1.0d0
             endif

             i = 7 - 4*icx -2*icy - icz               ! Calculate which child cell this is.


             do iox=0,1                               ! Loop over Legendre basis function type                     
              do ioy=0,1                              ! Loop over Legendre basis function type
               do ioz=0,1                             ! Loop over Legendre basis function type

                  j = 4*iox + 2*ioy + ioz

                  if (mod(iox*icx + ioy*icy + ioz*icz,2).eq.0) then
                       isig =  parity
                  else
                       isig = -parity
                  endif

                  temp_value(j,i) = temp_value(j,i) + isig*positive_octant_lc(j,ix,iy,iz)

               enddo
              enddo
             enddo

           enddo   
          enddo
         enddo

        enddo
       enddo
      enddo


c   Assign values of the output variables

      do i=0,7
       do j=0,7
         leg_coeff(j,i,nlev) = temp_value(j,i)*usually_rooteighth_factor
         cell_data(j,i)      = leg_coeff(j,i,nlev)
       enddo
      enddo

c   Finally set the independent field values
    
      cell_data(8,0) = g(57)
      cell_data(8,1) = g(58)
      cell_data(8,2) = g(59)
      cell_data(8,3) = g(60)
      cell_data(8,4) = g(61)
      cell_data(8,5) = g(62)
      cell_data(8,6) = g(63)
      cell_data(8,7) = g(64)


      return
      end
c=================================================================================
      recursive subroutine reset_lecuyer_state(ldata,lev,xcursor,ycursor,zcursor)
      use pan_state
      implicit none
      
      type(state_data), intent(inout) :: ldata
      integer lev
      integer*8 xcursor,ycursor,zcursor

c      integer indmin,indmax
c      parameter (indmin=-1, indmax=60)
c      integer*8 p_xcursor(indmin:indmax),p_ycursor(indmin:indmax),p_zcursor(indmin:indmax)
c      save p_xcursor,p_ycursor,p_zcursor
      integer i
c      integer init
c      data init/0/
c      save init

      if (ldata%reset_lecuyer_state_init.eq.0) then       ! Initialise p_cursor variables with 
          ldata%reset_lecuyer_state_init = 1              ! negative values.
          do i=indmin,indmax
            ldata%p_xcursor(i) = -9999
            ldata%p_ycursor(i) = -9999
            ldata%p_zcursor(i) = -9999
          enddo
      endif

      if ( (xcursor.eq.ldata%p_xcursor(lev)).and.(ycursor.eq.ldata%p_ycursor(lev)).and.
     &      (zcursor.eq.ldata%p_zcursor(lev)+1)) then
          ldata%p_xcursor(lev) = xcursor
          ldata%p_ycursor(lev) = ycursor
          ldata%p_zcursor(lev) = zcursor
          return
      endif
      
      call advance_current_state(ldata,lev,xcursor,ycursor,zcursor)
          
      ldata%p_xcursor(lev) = xcursor
      ldata%p_ycursor(lev) = ycursor
      ldata%p_zcursor(lev) = zcursor
    

      return
      end
c=================================================================================
      recursive subroutine advance_current_state(ldata,lev,x,y,z)
      use Rand
      use pan_state
      !use descriptor_phases
      implicit none

      type(state_data), intent(inout) :: ldata

      integer lev
      integer*8 x,y,z

      integer*8 lev_range

      TYPE(Rand_offset) :: offset1,offset2
      TYPE(Rand_offset) :: offset_x,offset_y,offset_z,offset_total

      integer ndiv,nrem
      integer*8 ndiv8,nrem8
      integer nfactor
      parameter (nfactor=291071) ! Value unimportant except has to be > 262144


c-----   First some error checking ------------------------------------------
      if ((lev.lt.0).or.(lev.gt.maxlev_)) stop 'Level out of range! (2)'

      lev_range = 2_dint**lev


      if ((x.lt.0).or.(x.ge.lev_range)) then 
      print*,'x,lev,lev_range',x,lev,lev_range
      call flush(6)
      stop 'x out of range!'
      endif
      if ((y.lt.0).or.(y.ge.lev_range)) then
      print*,'y,lev,lev_range',y,lev,lev_range       
          stop 'y out of range!'
      endif
      if ((z.lt.0).or.(z.ge.lev_range)) stop 'z out of range!' 
c----------------------------------------------------------------------------          
c
c Note the Rand_set_offset subroutine takes an integer*4 value
c for the offset value. For this reason we need to use integer*4
c values - ndiv,nrem.  As a precaution an explicit check is made
c to be sure that these values are calculated correctly.
c---------------------------------------------------------------------------


      call Rand_load(ldata%current_state(lev),ldata%base_lev_start(1,lev))

      if (lev.eq.0) return

c     Calculate z-offset

      ndiv = z/nfactor
      nrem = z - ndiv*nfactor
      ndiv8 = ndiv
      nrem8 = nrem

      if (ndiv8*nfactor+nrem8.ne.z)  stop 'Error in z ndiv nrem'

      call Rand_set_offset(offset1,ndiv)
      offset1 = Rand_mul_offset(offset1,nfactor)
      call Rand_set_offset(offset2,nrem)
      offset2 = Rand_add_offset(offset1,offset2)
      offset_z = Rand_mul_offset(offset2,nmulti_)

c     Calculate y-offset

      ndiv = y/nfactor
      nrem = y - ndiv*nfactor
      ndiv8 = ndiv
      nrem8 = nrem

      if (ndiv8*nfactor+nrem8.ne.y) stop 'Error in y ndiv nrem'

      offset1 =  Rand_mul_offset(ldata%poweroffset(lev),ndiv)
      offset1 =  Rand_mul_offset(offset1,nfactor)
      offset2 =  Rand_mul_offset(ldata%poweroffset(lev),nrem)
      offset_y = Rand_add_offset(offset1,offset2)

c     Calculate x-offset

      ndiv = x/nfactor
      nrem = x - ndiv*nfactor
      ndiv8 = ndiv
      nrem8 = nrem

      if (ndiv8*nfactor+nrem8.ne.x) then
           print*,'ndiv,nfactor,nrem,x',ndiv,nfactor,nrem,x
           print*,'ndiv*nfactor+nrem',ndiv*nfactor+nrem
           print*,'x-ndiv*nfactor-nrem',x-ndiv*nfactor-nrem
           stop 'Error in x ndiv nrem'
      endif

      offset1 = Rand_mul_offset(ldata%poweroffset(2*lev),ndiv)
      offset1 = Rand_mul_offset(offset1,nfactor)
      offset2 = Rand_mul_offset(ldata%poweroffset(2*lev),nrem)
      offset_x = Rand_add_offset(offset1,offset2)

      offset1      = Rand_add_offset(offset_x,offset_y)
      offset_total = Rand_add_offset(offset1, offset_z)
 
      ldata%current_state(lev) = Rand_boost(ldata%current_state(lev),offset_total)
      
      return
      end
c=================================================================================
      recursive subroutine return_gaussian_array(ldata,lev,ngauss,garray)
      use Rand
      use pan_state
      implicit none 
      type(state_data), intent(inout) :: ldata
      integer lev,ngauss
      real*8 garray(0:*)
      TYPE(Rand_state) :: state
      real*8 PI
      parameter (PI=3.1415926535897932384d0)
      real*8 branch
      parameter (branch=1.d-6)
      integer iloop

      real*8 temp,mag,ang
      integer i

      if (mod(ngauss,2).ne.0) 
     & stop 'Error in return_gaussian_array - even pairs only'

c   First obtain a set of uniformly distributed pseudorandom numbers
c   between 0 and 1. The method used is described in detail in 
c   appendix B of Jenkins 2013.

      do i=0,ngauss-1
       call Rand_real(garray(i),ldata%current_state(lev))

       if (garray(i).lt.branch) then
          garray(i) = branch
          state = Rand_boost(ldata%current_state(lev),ldata%superjump)
          iloop = 0
 10       continue
          call Rand_real(temp,state)
          iloop = iloop+1
          if (temp.lt.branch) then
               garray(i) = garray(i)*branch
               state = Rand_boost(state,ldata%superjump)
               if (iloop.gt.100) then
               print*,'Too may iterations in return_gaussian_array!'
               call flush(6)
               stop
               endif               
               goto 10
          else
               garray(i) = garray(i)*temp
          endif
       endif
      enddo

c     Apply Box-Muller transformation to create pairs of Gaussian
c     pseudorandom numbers.

      do i=0,ngauss/2-1

       mag = sqrt(-2.0d0*log(garray(2*i)))
       ang = 2.0d0*PI*garray(2*i+1)

       garray(2*i)   = mag*cos(ang)
       garray(2*i+1) = mag*sin(ang)
 
      enddo
      end
c=================================================================================
      recursive subroutine parse_descriptor(string,l,ix,iy,iz,side1,side2,side3,check_int,name)
      implicit none
      integer nchar
      parameter(nchar=100)
      character*100  string
      integer*4 l,side1,side2,side3,ierror
      integer*8 ix,iy,iz
      integer*8 check_int
      character*20 name


      integer i,ip,iq,ir

      ierror = 0

      ip = 1
      do while (string(ip:ip).eq.' ')
       ip = ip + 1
      enddo

      if (string(ip:ip+7).ne.'[Panph1,') then
           ierror = 1
           print*,string(ip:ip+7)
           goto 10
      endif

      ip = ip+8
      if (string(ip:ip).ne.'L') then
          ierror = 2
          goto 10 
      endif

      ip = ip+1

      iq  = ip + scan( string(ip:nchar),',') -1

      if (ip.eq.iq) then
          ierror = 3
          goto 10
      endif


      read (string(ip:iq),*) l

      ip = iq+1
      
      if (string(ip:ip).ne.'(') then
         ierror = 4
         goto 10
      endif

      ip = ip+1

      iq = ip + scan( string(ip:nchar),')') -2

      read(string(ip:iq),*) ix,iy,iz

      ip = iq+2
      
      if (string(ip:ip).ne.',') then
         ierror = 5
         goto 10
      endif

      ip = ip+1
      if ((string(ip:ip).ne.'S').and.(string(ip:ip).ne.'D')) then
         ierror = 6
         goto 10
      endif

      if (string(ip:ip).eq.'S') then
        ip = ip + 1
        iq = ip + scan( string(ip:nchar),',') -2
        read (string(ip:iq),*) side1
        side2 = side1
        side3 = side1
        iq = iq+1
        if (string(iq:iq+2).ne.',CH') then
           print*,string(ip:iq),string(iq:iq+2)
           ierror = 6
           goto 10
        endif
      else
        ip = ip + 1
        if (string(ip:ip).ne.'(') then
           ierror = 7
           goto 10
        endif


        ip = ip + 1
        iq = ip + scan( string(ip:nchar),')') -2
        read (string(ip:iq),*) side1,side2,side3

        iq = iq + 1

        if (string(iq:iq).ne.')') then
           ierror = 8
           goto 10
        endif

        iq = iq + 1

         if (string(iq:iq+2).ne.',CH') then
            ierror = 9
            goto 10
        endif

      endif

      ip = iq + 3

      iq = ip + scan( string(ip:nchar),',') -2

      read (string(ip:iq),*) check_int

      ip = iq + 1

      if (string(ip:ip).ne.',') then
          ierror = 10
          goto 10
      endif

      ip = ip+1

      ir = ip + scan( string(ip:nchar),']') -2

      iq = min(ir,ip+19)

      do i=1,20
        name(i:i)=' '
      enddo

      do i=ip,iq
        name(i-ip+1:i-ip+1) = string(i:i)
      enddo

      iq = ir + 1

      if (string(iq:iq).ne.']') then
          ierror = 11
          goto 10
      endif


 10   continue

      if (ierror.eq.0) return

      print*,'Error reading panphasian descriptor. Error number:',ierror
      stop

      return
      end
c=================================================================================
      recursive subroutine compose_descriptor(l,ix,iy,iz,side,check_int,name,string)
      implicit none
      integer nchar
      parameter(nchar=100)
      character*100,intent(out)::string
      character*20 name
      integer*4 l,ltemp
      integer*8 side
      integer*8 ix,iy,iz
      integer*8 check_int

      character*50 temp1,temp2,temp3,temp4,temp5,temp6
      integer lnblnk

      integer ip1,ip2,ip3,ip4,ip5,ip6
      
      ltemp = l

 5    continue
      if ((mod(ix,2).eq.0).and.(mod(iy,2).eq.0).and.(mod(iz,2).eq.0).and.(mod(side,2).eq.0)) then
        ix = ix/2
        iy = iy/2
        iz = iz/2
        side = side/2
        ltemp = ltemp-1
        goto 5
      endif


      write (temp1,*) ltemp
      ip1= scan(temp1,'0123456789')
      write (temp2,*) ix
      ip2= scan(temp2,'0123456789')
      write (temp3,*) iy
      ip3= scan(temp3,'0123456789')
      write (temp4,*) iz
      ip4= scan(temp4,'0123456789')
      write (temp5,*) side
      ip5= scan(temp5,'0123456789')
      write (temp6,*) check_int
      ip6= scan(temp6,'-0123456789')


      string='[Panph1,L'//temp1(ip1:lnblnk(temp1))//',('//temp2(ip2:lnblnk(temp2))
     &   //','//temp3(ip3:lnblnk(temp3))//','//temp4(ip4:lnblnk(temp4))//'),S'
     &   // temp5(ip5:lnblnk(temp5))//',CH'//temp6(ip6:lnblnk(temp6))//
     &  ','//name(1:lnblnk(name))//']'
 
      return 

      end
c=================================================================================
      recursive subroutine validate_descriptor(ldata,string,MYID,check_number)
      use pan_state
      implicit none

      type(state_data), intent(inout) :: ldata
      character*100 string
      integer*8 check_number
      integer MYID

      character*20 phase_name
      integer*4 lev

      integer*8 ix_abs,iy_abs,iz_abs
      integer*4 ix_base,iy_base,iz_base
      

      integer*8 xval,yval,zval
      integer val_state(5)

      TYPE(Rand_state) :: state

      real*8 rand_num
      integer*8 mconst,check_total,check_rand
      parameter(mconst = 2147483647_Dint)
      integer ascii_list(0:255)
      integer*8 maxco
      integer i
      integer*8 ii
      integer lnblnk
      


      call parse_descriptor(string,lev,ix_abs,iy_abs,iz_abs,
     &                  ix_base,iy_base,iz_base,check_rand,phase_name)

c-------------------------------------------------------------------------
c    Some basic checking
c-------------------------------------------------------------------------
      if ((lev.lt.0).or.(lev.gt.maxlev_)) then
            print*,'lev,maxlev',lev,maxlev_
            call flush(6)
            stop 'Level out of range! (3)'
      endif

      if ((mod(ix_abs,2).eq.0).and.(mod(iy_abs,2).eq.0).and.(mod(iz_abs,2).eq.0).and.
     & (mod(ix_base,2).eq.0).and.(mod(iy_base,2).eq.0).and.(mod(iz_base,2).eq.0)) 
     &    stop 'Parameters not at lowest level'


      maxco = 2_dint**lev

      if (ix_abs.lt.0) stop 'Error: ix_abs negative (2)'
      if (iy_abs.lt.0) stop 'Error: iy_abs negative (2)'
      if (iz_abs.lt.0) stop 'Error: iz_abs negative (2)'


      if (ix_abs+ix_base.ge.maxco)
     &   stop 'Error: ix_abs + ix_per out of range.'
      if (iy_abs+iy_base.ge.maxco) 
     &   stop 'Error: iy_abs + iy_per out of range.'
      if (iz_abs+iz_base.ge.maxco) 
     &   stop 'Error: iz_abs + iz_per out of range.'

      check_total = 0

      call initialise_panphasia(ldata)
c    First corner
      xval = ix_abs + ix_base - 1
      yval = iy_abs
      zval = iz_abs
      call advance_current_state(ldata,lev,xval,yval,zval)
      call Rand_real(rand_num,ldata%current_state(lev))
      call Rand_save(val_state,ldata%current_state(lev))
      check_total = check_total + val_state(5)
      if (MYID.eq.0) print*,'--------------------------------------'
      if (MYID.eq.0)  print*,'X-corner rand = ',rand_num
      if (MYID.eq.0) print*,'State:',val_state
c    Second corner
      xval = ix_abs
      yval = iy_abs + iy_base - 1
      zval = iz_abs
      call advance_current_state(ldata,lev,xval,yval,zval)
      call Rand_real(rand_num,ldata%current_state(lev))
      call Rand_save(val_state,ldata%current_state(lev))
      check_total = check_total + val_state(5)
      if (MYID.eq.0)  print*,'Y-corner rand = ',rand_num
      if (MYID.eq.0) print*,'State:',val_state
c    Third corner
      xval = ix_abs
      yval = iy_abs
      zval = iz_abs + iz_base - 1
      call advance_current_state(ldata,lev,xval,yval,zval)
      call Rand_real(rand_num,ldata%current_state(lev))
      call Rand_save(val_state,ldata%current_state(lev))
      check_total = check_total + val_state(5)
      if (MYID.eq.0)  print*,'z-corner rand = ',rand_num
      if (MYID.eq.0) print*,'State:',val_state
      if (MYID.eq.0) print*,'--------------------------------------'

c     Now encode the name.  An integer for each ascii character is generated
c     starting from the state which gives r0 - the first random number in
c     Panphasia.   The integer is in the range 0 - m-1.  
c     After making the list, then loop over non-blank characters
c     in the name and take the ascii value, and sum the associated numbers.
c     To avoid simple anagrams giving the same score, weight the integer
c     by position in the string.  Finally take mod m  - to give the
c     check number.  

      call Rand_load(state,ldata%base_state)

      do i=0,255
      call Rand_real(rand_num,state)
      call Rand_save(val_state,state)
      ascii_list(i) = val_state(5)
      enddo



      do ii=1,lnblnk(phase_name)
       check_total = check_total + ii*ascii_list(iachar(phase_name(ii:ii)))
      enddo


      check_total = mod(check_total,mconst)
      if (check_rand.eq.-999) then         ! override the safety check number.
            check_number = check_total
            return
      else
          if (check_rand.ne.check_total) then
           print*,'Inconsistency in the input panphasia descriptor ',MYID
           print*,'Check_rand  = ',check_rand
           print*,'val_state(5) =',val_state(5)
           print*,'xval,yval,zval',xval,yval,zval
           print*,'lev_val =  ',lev
           call flush(6)
           stop
          endif
      endif


      return
      end
c=================================================================================
      recursive subroutine generate_random_descriptor(ldata,string)
      use Rand
      use pan_state
      implicit none
      type(state_data), intent(inout) :: ldata
      character*100  string
      character*100  instring
      character*20   name
      integer*4 unix_timestamp

      real*8 lbox
      real*8 lpanphasia
      parameter (lpanphasia = 25000000.0)  ! Units of Mpc/h
      integer level
      integer*8 cell_dim
      integer val_state(5)

      TYPE(Rand_state) :: state
      TYPE(Rand_offset) :: offset

      real*8 rand_num1,rand_num2
      integer*8 mconst,check_int
      parameter(mconst = 2147483647_Dint)
      integer*8 mfac,imajor,iminor
      parameter(mfac=33554332_Dint)
      integer ascii_list(0:255)
      integer i,lnblnk
      integer*8 ii
      integer mult

      integer*8 ixco,iyco,izco,irange

      print*,'___________________________________________________________'
      print*
      print*,'            Generate a random descriptor                   '
      print*
      print*,'The code uses the time (the unix timestamp) plus some extra '
      print*,'information personal to the user to choose a random region  '
      print*,'within PANPHASIA.  The user must also specify the side length'
      print*,'of the cosmological volume. The code assumes that the whole of'
      print*,'PANPHASIA is 25000 Gpc/h on a side and selects an appropriate '
      print*,'level in the octree for the descriptor.  '
      print*,'Assuming this scaling the small scale power is defined down '
      print*,'to a mass scale of around 10^{-12} solar masses.'
      print*      
      print*,'The user must also specify a human readable label for the '
      print*,'descriptor of less than 21 characters.'
      print*,'___________________________________________________________'
      print*
      print*,'Press return to continue '
      read (*,*) 
      print*
      print*,'___________________________________________________________'
      print*,'Enter the box side-length in Mpc/h units'
      read (*,*) lbox
      print*,'___________________________________________________________'
      print*
      print*
 5    continue
      print*,'Enter up to 20 character name to label the descriptor (no spaces)'
      read (*,'(a)') name
      if ((len_trim(instring).lt.21).or.(scan(name,' ').le.len_trim(name))) goto 5
      print*,'___________________________________________________________'
      print*
      print* 
      print*,'___________________________________________________________'
      print*,'The phases for the simulation are described by whole octree '
      print*,'cells. Enter an odd integer that defines the number of cells '
      print*,'you require in one dimension.  Choose this number carefully  '
      print*,'as it will limit the possible 1-D sizes of  the of the Fourier '
      print*,'transforms that can be used to make initial conditions to a product '
      print*,'of this integer times any power of two. In which case the only'
      print*,'choice is 1.)'
      print*,'(I would recommend 3 unless the initial condition code is'
      print*,'incapable of using grid sizes that are not purely powers of two.'
      print*,'___________________________________________________________'
      print*
 7    continue
      print*,'Enter number of octree cells on an edge (positive odd number only) '
      read (*,*) cell_dim
      if ((cell_dim.le.0).or.(mod(cell_dim,2).eq.0)) goto 7
      print*,'___________________________________________________________' 
      call system('date +%s>tempfile_42526037646')
      open(16,file='tempfile_42526037646',status='old')
      read (16,*) unix_timestamp
      close(16)
      call system('/bin/rm tempfile_42526037646') 

      print*,'Unix_timestamp determined. Value: ',unix_timestamp
      print*,'___________________________________________________________'
      print*
      print*
      print*
      print*,'___________________________________________________________'
      print*,'The code has just read the unix timestamp and will use this'
      print*,'to help choose a random region in PANPHASIA.  Although it is'
      print*,'perhaps unlikely that someone else is also running this code at '
      print*,'the same time to the nearest second, to make it more likely'
      print*,' still that the desciptor to be generated is unique'
      print*,'please enter your name or some other piece of information'
      print*,'below that you think is unlikely to be used by anyone else'
      print*,'___________________________________________________________'

      print*

 10   continue
      print*,'Please enter your name (a minimum of six characters)'
      read (*,'(a)') instring                         !'
      if (len_trim(instring).lt.6) goto 10

      level =  int(log10(dble(cell_dim)*lpanphasia/lbox)/log10(2.0d0))

      if (level.gt.50) stop 'level >50 '



c      'd' lines allow the generation of a large set of
c       descriptors. Use to check that they are randomly
c       positioned over the available volume.


c    First use the unix timestamp to initialises the
c    random generator.

      call Rand_seed(state,unix_timestamp)
      
      call Rand_save(ldata%base_state,state)
     

c   First generate an integer from the user data.
      call Rand_load(state,ldata%base_state)

      do i=0,255
      call Rand_real(rand_num1,state)
      call Rand_save(val_state,state)
      ascii_list(i) = val_state(5)
      enddo

      call Rand_set_offset(offset,1)

      do ii=1,lnblnk(instring)
       mult = mod(ii*ascii_list(iachar(instring(ii:ii))),mconst)
       offset =  Rand_mul_offset(offset,mult)
      enddo

      call Rand_load(state,ldata%base_state)
      state = Rand_boost(state,offset)          ! Starting point for choosing location.

 20   continue

      irange = 2_Dint**level
      imajor = irange/mfac
      iminor = mod(irange,mfac)

      call Rand_real(rand_num1,state)
      call Rand_real(rand_num2,state)

      ixco = int(rand_num1*imajor)*mfac + int(rand_num2*iminor)

      if (ixco+cell_dim.ge.irange) goto 20      ! Invalid descriptor

      call Rand_real(rand_num1,state)
      call Rand_real(rand_num2,state)

      iyco = int(rand_num1*imajor)*mfac + int(rand_num2*iminor)

      if (iyco+cell_dim.ge.irange) goto 20      ! Invalid descriptor

      call Rand_real(rand_num1,state)
      call Rand_real(rand_num2,state)

      izco = int(rand_num1*imajor)*mfac + int(rand_num2*iminor)

      if (izco+cell_dim.ge.irange) goto 20      ! Invalid descriptor


c     Value of the check digit is not known. Use validate_descriptor to compute it.

      check_int = -999  ! Special value required to make validate_descriptor 
                        ! return the check digit.

      call compose_descriptor(level,ixco,iyco,izco,cell_dim,check_int,name,string)

      call validate_descriptor(ldata,string,-1,check_int)

      call compose_descriptor(level,ixco,iyco,izco,cell_dim,check_int,name,string)


      return
      end
c=================================================================================
      recursive subroutine demo_basis_function_allocator

      implicit none
      integer nmax
      parameter (nmax=10)

      integer*4 wn_level(nmax)

      integer*8 ix_abs(nmax),iy_abs(nmax),iz_abs(nmax)
      integer*8 ix_per(nmax),iy_per(nmax),iz_per(nmax)
      integer*8 ix_rel(nmax),iy_rel(nmax),iz_rel(nmax)
      integer*8 ix_dim(nmax),iy_dim(nmax),iz_dim(nmax)

      integer ix,iy,iz,nref
      integer layer_min,layer_max,indep_field


      integer*8 itot_int,itot_ib

      integer inv_open

c      Assign some trial values

      nref = 3
      inv_open=9

      wn_level(1) = 22

      ix_abs(1) = 2000000
      iy_abs(1) = 1500032
      iz_abs(1) = 2500032

      ix_per(1) = 768
      iy_per(1) = 768
      iz_per(1) = 768

      ix_rel(1) = 0
      iy_rel(1) = 0
      iz_rel(1) = 0

      ix_dim(1) = 768
      iy_dim(1) = 768
      iz_dim(1) = 768


      wn_level(2) = 23

      ix_abs(2) = 4000000
      iy_abs(2) = 3000064
      iz_abs(2) = 5000064

      ix_per(2) = 1536
      iy_per(2) = 1536
      iz_per(2) = 1536

      ix_rel(2) = 256
      iy_rel(2) = 16
      iz_rel(2) = 720

      ix_dim(2) = 768
      iy_dim(2) = 768
      iz_dim(2) = 768


      wn_level(3) = 24

      ix_abs(3) = 8000000
      iy_abs(3) = 6000128
      iz_abs(3) = 10000128

      ix_per(3) = 3072
      iy_per(3) = 3072
      iz_per(3) = 3072

      ix_rel(3) = 896
      iy_rel(3) = 432
      iz_rel(3) = 1840

      ix_dim(3) = 768
      iy_dim(3) = 768
      iz_dim(3) = 768


      itot_int = 0
      itot_ib  = 0




      open(10,file='ascii_dump_r1',status='unknown')

      ix=320
      do iy=0,767
       do iz=0,767
        call layer_choice(ix,iy,iz,1,nref,ix_abs,iy_abs,iz_abs,
     &   ix_per,iy_per,iz_per,ix_rel,iy_rel,iz_rel,ix_dim,iy_dim,iz_dim,
     &   wn_level,inv_open,layer_min,layer_max,indep_field) 
        write(10,*) iy,iz,layer_min,layer_max,indep_field
       enddo
      enddo
      close(10)

      open(10,file='ascii_dump_r2',status='unknown')

      ix=384
      do iy=0,767
       do iz=0,767
        call layer_choice(ix,iy,iz,2,nref,ix_abs,iy_abs,iz_abs,
     &   ix_per,iy_per,iz_per,ix_rel,iy_rel,iz_rel,ix_dim,iy_dim,iz_dim,
     &   wn_level,inv_open,layer_min,layer_max,indep_field) 
        write(10,*) iy,iz,layer_min,layer_max,indep_field
       enddo
      enddo
      close(10)

      open(10,file='ascii_dump_r3',status='unknown')

      ix=384
      do iy=0,767
       do iz=0,767
        call layer_choice(ix,iy,iz,3,nref,ix_abs,iy_abs,iz_abs,
     &   ix_per,iy_per,iz_per,ix_rel,iy_rel,iz_rel,ix_dim,iy_dim,iz_dim,
     &   wn_level,inv_open,layer_min,layer_max,indep_field) 
        write(10,*) iy,iz,layer_min,layer_max,indep_field
       enddo
      enddo
      close(10)
      end
c=================================================================================    
      recursive subroutine layer_choice(ix0,iy0,iz0,iref,nref,
     &  ix_abs,iy_abs,iz_abs,ix_per,iy_per,iz_per,
     &  ix_rel,iy_rel,iz_rel,ix_dim,iy_dim,iz_dim,
     &  wn_level,x_fact,layer_min,layer_max,indep_field)
      implicit none

      integer ix0,iy0,iz0,iref,nref,isize,ibase
      integer ix,iy,iz,irefplus
      integer ione

      integer*8 ix_abs(nref),iy_abs(nref),iz_abs(nref)
      integer*8 ix_per(nref),iy_per(nref),iz_per(nref)
      integer*8 ix_rel(nref),iy_rel(nref),iz_rel(nref)
      integer*8 ix_dim(nref),iy_dim(nref),iz_dim(nref)

      integer wn_level(nref)
      integer layer_min,layer_max,indep_field,x_fact
      integer idebug


      integer interior,iboundary

      if (iref.eq.9999) then
         idebug = 1
      else
         idebug = 0
      endif

      ione =  1

      irefplus = min(iref+1,nref)

      if (nref.eq.1) then            ! Deal with simplest case
         layer_min = 0
         layer_max = wn_level(1)
         indep_field  = 1
         if (idebug.eq.1) print*,'return 1'
         return
      endif 

c-----------  Case of the top periodic refinement.  For this refinement layer_min=0 as
c-----------  all the larger basis functions must be included.  By default layer_max
c-----------  is set to wn_level(1) so all basis functions are included. A check is
c-----------  made to determine if the lowest basis function can be included in the
c-----------  next refinement. If it can the same process is repeated for the next
c-----------  largest basis function and this is repeated until a failure occurs.

      if ((iref.eq.1).and.(nref.gt.1)) then
         ibase = 1
 10      continue
         
         ix = ishft(ishft(ix_abs(iref)+ix_rel(iref)+ix0,-ibase),ibase)-ix_abs(iref)-ix_rel(iref)
         iy = ishft(ishft(iy_abs(iref)+iy_rel(iref)+iy0,-ibase),ibase)-iy_abs(iref)-iy_rel(iref)
         iz = ishft(ishft(iz_abs(iref)+iz_rel(iref)+iz0,-ibase),ibase)-iz_abs(iref)-iz_rel(iref)
         isize = ishft(ione,ibase)

         call inref(ix,iy,iz,isize,iref,irefplus,nref,wn_level,
     &   ix_abs,iy_abs,iz_abs,ix_per,iy_per,iz_per,
     &   ix_rel,iy_rel,iz_rel,ix_dim,iy_dim,iz_dim,x_fact,
     &   interior,iboundary) 

         if ((interior.eq.1).and.(iboundary.eq.1)) then
            ibase = ibase + 1
            goto 10
         endif

           layer_min = 0
           layer_max = wn_level(iref) - ibase + 1
         if (layer_max.ne.wn_level(iref)) then
           indep_field = 0
         else
           indep_field = 1
         endif

         if (idebug.eq.1) then
         print*,'iref,wn_level(iref)',iref,wn_level(iref)
         print*,'Return 2',layer_min,layer_max,indep_field
         endif

         return
      endif
c------------------------------------------------------------------------------------------
c------------------------------------------------------------------------------------------


c-----------  For second or higher refinement determine layer_min by reference 
c-----------  to itself.  In this case the loop continues until a basis function
c------------ is found which fits in a larger refinement
   
         ibase = 1

 20      continue


         ix = ishft(ishft(ix_abs(iref)+ix_rel(iref)+ix0,-ibase),ibase)-ix_abs(iref)-ix_rel(iref)
         iy = ishft(ishft(iy_abs(iref)+iy_rel(iref)+iy0,-ibase),ibase)-iy_abs(iref)-iy_rel(iref)
         iz = ishft(ishft(iz_abs(iref)+iz_rel(iref)+iz0,-ibase),ibase)-iz_abs(iref)-iz_rel(iref)
         isize = ishft(ione,ibase)

         call inref(ix,iy,iz,isize,iref,iref,nref,wn_level,
     &   ix_abs,iy_abs,iz_abs,ix_per,iy_per,iz_per,
     &   ix_rel,iy_rel,iz_rel,ix_dim,iy_dim,iz_dim,x_fact,
     &   interior,iboundary) 

         if ((interior.eq.1).and.(iboundary.eq.1)) then
            ibase = ibase + 1
            goto 20
         endif

         layer_min = wn_level(iref) - max(ibase-2,0)         ! Take last suitable refinement


c-----------  For an intermediate refinement define layer_max by reference to
c-----------  the next refinement

         if (iref.lt.nref) then
         ibase = 1

 30          continue

            ix = ishft(ishft(ix_abs(iref)+ix_rel(iref)+ix0,-ibase),ibase)-ix_abs(iref)-ix_rel(iref)
            iy = ishft(ishft(iy_abs(iref)+iy_rel(iref)+iy0,-ibase),ibase)-iy_abs(iref)-iy_rel(iref)
            iz = ishft(ishft(iz_abs(iref)+iz_rel(iref)+iz0,-ibase),ibase)-iz_abs(iref)-iz_rel(iref)
            isize = ishft(ione,ibase)

            call inref(ix,iy,iz,isize,iref,irefplus,nref,wn_level, 
     &      ix_abs,iy_abs,iz_abs,ix_per,iy_per,iz_per,
     &      ix_rel,iy_rel,iz_rel,ix_dim,iy_dim,iz_dim,x_fact,
     &      interior,iboundary) 

            if ((interior.eq.1).and.(iboundary.eq.1)) then
               ibase = ibase + 1
               goto 30
            endif

            layer_max = wn_level(iref) - ibase + 1

            if (layer_min.eq.wn_level(iref)) then
               indep_field = 1
            else
               indep_field = 0
            endif
         else
            layer_max = wn_level(iref)
            indep_field  = 1
         endif

         if (idebug.eq.1) then
           print*,'Return 3'
           print*,'layer_min,layer_max,indep_field',layer_min,layer_max,indep_field
           print*,'interior,iboundary',interior,iboundary
           print*,'ibase = ',ibase
           print*,'iref,nref,wn_level(iref)',iref,nref,wn_level(iref)
         endif


         return
     
      end




c   The function takes a given basis function specified by a corner ixc,iyc,izc
c   and a size isz at level wn_c in the oct-tree and returns two integer values.
c   (i)  interior:     
c                  Value 1 if the basis function is completely within the given
c                  refinement.
c
c                  Value 0 if the basis function is without the refinement, or
c                  overlaps the edges of the refinement, or the edges of the
c                  primary white noise patch.
c
c   (ii) iboundary:
c                  Value 1 if the basis function is sufficiently far from the
c                  refinement boundary.
c
c                  Value 0 otherwise.
c   The given refinement is defined at level wn_r in the oct-tree and by the variables
c   (ix_rel,iy_rel,iz_rel) which give the location of the refinement relative to
c   corner of the white noise patch, (ix_per,iy_per,iz_per) which define the
c   periodicity of the white noise patch, and (ix_dim,iy_dim,iz_dim) which
c   define the size of the refinement.
c
c
c
c=================================================================================
      recursive subroutine inref(ixc,iyc,izc,isz,ir1,ir2,nref,wn_level,
     &   ix_abs,iy_abs,iz_abs,ix_per,iy_per,iz_per,
     &   ix_rel,iy_rel,iz_rel,ix_dim,iy_dim,iz_dim,x_fact,
     &   interior,iboundary)
      implicit none

      integer nref
      integer ixc,iyc,izc,isz,ir1,ir2
      integer wn_level(nref)
      integer*8 ix_abs(nref),iy_abs(nref),iz_abs(nref)
      integer*8 ix_per(nref),iy_per(nref),iz_per(nref)
      integer*8 ix_rel(nref),iy_rel(nref),iz_rel(nref)
      integer*8 ix_dim(nref),iy_dim(nref),iz_dim(nref)
      integer interior, iboundary
      integer x_fact
     
      integer*8 ixco,iyco,izco,isize
      integer*8 ixref0,iyref0,izref0
      integer*8 ixref1,iyref1,izref1
      integer*8 idist

      integer delta_wn

c   Error checking
      if (ir2.lt.ir1) stop 'ir2<ir1'
      if ((ir1.lt.1).or.(ir2.gt.nref)) 
     &  stop 'Either/or ir1,ir2 out of range'

c  First copy coordinates to integer*8 variables

      ixco = ixc
      iyco = iyc
      izco = izc
      isize= isz

      delta_wn = wn_level(ir2)-wn_level(ir1)

c  Now translate coordinates from refinement ir1 to ir2 and express relative
c  to the origin of refinement 2.

      ixco =     ishft(ix_abs(ir1)+ix_rel(ir1)+ixco,delta_wn)-ix_abs(ir2)-ix_rel(ir2)
      iyco =     ishft(iy_abs(ir1)+iy_rel(ir1)+iyco,delta_wn)-iy_abs(ir2)-iy_rel(ir2)
      izco =     ishft(iz_abs(ir1)+iz_rel(ir1)+izco,delta_wn)-iz_abs(ir2)-iz_rel(ir2)
      isize=     ishft(isize,delta_wn)

      ixref0 = mod(ix_per(ir2) + ixco, ix_per(ir2))
      iyref0 = mod(iy_per(ir2) + iyco, iy_per(ir2))
      izref0 = mod(iz_per(ir2) + izco, iz_per(ir2))

      if ((ixref0.ge.ix_dim(ir2)).or.(iyref0.ge.iy_dim(ir2)).or.(izref0.ge.iz_dim(ir2))) then !
           interior =  0
          iboundary =  0
          return          ! The basis function is not inside the refinement
      endif

      ixref1 = mod(ix_per(ir2) + ixco + isize, ix_per(ir2))
      iyref1 = mod(iy_per(ir2) + iyco + isize, iy_per(ir2))
      izref1 = mod(iz_per(ir2) + izco + isize, iz_per(ir2))

      if ((ixref1.ge.ix_dim(ir2)).or.(iyref1.ge.iy_dim(ir2)).or.(izref1.ge.iz_dim(ir2))) then ! Location not in refinement
           interior =  0
          iboundary =  0
          return          ! The basis function is not inside the refinement
      endif

c     The basis function is within the refinement. Now calculate the 
c     minimum perpendicular distance of the basis function from the 
c     edge of the refinement.


      idist = min(ixref0,ix_dim(ir2)-ixref1,iyref0,iy_dim(ir2)-iyref1, izref0,iz_dim(ir2)-izref1)

      if (idist.gt.x_fact*isize) then  
          iboundary = 1                ! Sufficiently far from the boundary
      else
          iboundary = 0
      endif

c  Final check - does the basis function reside entirely in the white noise patch.

      ixref0 = mod(ix_rel(ir2)+ixco       ,ix_per(ir2))
      ixref1 = mod(ix_rel(ir2)+ixco +isize,ix_per(ir2))

      iyref0 = mod(iy_rel(ir2)+iyco       ,iy_per(ir2))
      iyref1 = mod(iy_rel(ir2)+iyco +isize,iy_per(ir2))

      izref0 = mod(iz_rel(ir2)+izco       ,iz_per(ir2))
      izref1 = mod(iz_rel(ir2)+izco +isize,iz_per(ir2))

      if ((ixref1.le.ixref0).or.(iyref1.le.iyref0).or.(izref1.le.izref0)) then
         interior  = 0
         iboundary = 0
         return          ! Basis function not completely in the refinement:
      endif              ! crosses white noise patch boundary.
      
      interior = 1    ! Basis function is completely within the refinement

      return   
      end
c==========================================================================================
      recursive subroutine set_local_box(ldata,lev,ix_abs,iy_abs,iz_abs,
     & ix_per,iy_per,iz_per, ix_rel,iy_rel,iz_rel,wn_level_base,check_rand,phase_name,MYID)
      use pan_state
      !use descriptor_phases
      implicit none

      type(state_data), intent(inout) :: ldata
      !integer  layer_min,layer_max,indep_field
      !common /oct_range/  layer_min,layer_max,indep_field


      integer lev
      integer*8 ix_abs,iy_abs,iz_abs
      integer*8 ix_per,iy_per,iz_per
      integer*8 ix_rel,iy_rel,iz_rel
      integer*8 xorigin,yorigin,zorigin
      integer wn_level_base
      integer*8 check_rand
      character*20 phase_name
      integer MYID
      integer*8 maxco
      integer i
      integer px,py,pz

      integer*8 xval,yval,zval,val_side
      integer lev_val
      character*100 outstring
      integer lnblnk
      integer*8 mconst
      parameter(mconst = 2147483647_Dint)
c-------------------------------------------------------------------------

      call initialise_panphasia(ldata)

c-------------------------------------------------------------------------
c    Error checking
c-------------------------------------------------------------------------
      if ((lev.lt.0).or.(lev.gt.maxlev_)) stop 'Level out of range! (4)'

      maxco = 2_dint**lev

      if (ix_abs.lt.0) stop 'Error: ix_abs negative (3)'
      if (iy_abs.lt.0) stop 'Error: iy_abs negative (3)'
      if (iz_abs.lt.0) stop 'Error: iz_abs negative (3)'

      if (ix_rel.lt.0) stop 'Error: ix_rel negative (2)'
      if (iy_rel.lt.0) stop 'Error: iy_rel negative (2)'
      if (iz_rel.lt.0) stop 'Error: iz_rel negative (2)'

      if (ix_abs+ix_rel.ge.maxco)
     &   stop 'Error: ix_abs + ix_rel out of range. (2)'
      if (iy_abs+iy_rel.ge.maxco) 
     &   stop 'Error: iy_abs + iy_rel out of range. (2)'
      if (iz_abs+iz_rel.ge.maxco) 
     &   stop 'Error: iz_abs + iz_rel out of range. (2)'
c-----------------------------------------------------------------------------
c  To allow the local box to wrap around, if needed, define a series of eight
c  'origins'.  For many purposes (ix,iy,iz) = (0,0,0) is the only origin needed.
c-----------------------------------------------------------------------------

      do px=0,1
       do py=0,1
        do pz=0,1

         xorigin = max(0,( ix_abs + ix_rel - px*ix_per )/2)
         yorigin = max(0,( iy_abs + iy_rel - py*iy_per )/2)
         zorigin = max(0,( iz_abs + iz_rel - pz*iz_per )/2)

         ldata%ixshift(px,py,pz) = max(0, ix_abs + ix_rel -px*ix_per) - 2*xorigin
         ldata%iyshift(px,py,pz) = max(0, iy_abs + iy_rel -py*iy_per) - 2*yorigin
         ldata%izshift(px,py,pz) = max(0, iz_abs + iz_rel -pz*iz_per) - 2*zorigin


c        Store box details:  store the positions at level lev-1
  

         ldata%xorigin_store(px,py,pz) = xorigin
         ldata%yorigin_store(px,py,pz) = yorigin
         ldata%zorigin_store(px,py,pz) = zorigin

        enddo
       enddo
      enddo

      ldata%lev_common = lev


      ldata%ix_abs_store = ix_abs
      ldata%iy_abs_store = iy_abs
      ldata%iz_abs_store = iz_abs

      ldata%ix_per_store = ix_per
      ldata%iy_per_store = iy_per
      ldata%iz_per_store = iz_per

      ldata%ix_rel_store = ix_rel
      ldata%iy_rel_store = iy_rel
      ldata%iz_rel_store = iz_rel

c------  Now validate the panphasian descriptor ---------------------------------
c-----   Use lowest level possible
      lev_val = wn_level_base
      xval = ix_abs/2_dint**(lev-lev_val)
      yval = iy_abs/2_dint**(lev-lev_val)
      zval = iz_abs/2_dint**(lev-lev_val)
      val_side = ix_per/2_dint**(lev-lev_val)
      call compose_descriptor(lev_val,xval,yval,zval,val_side,check_rand,phase_name,outstring)
      print*,'blabla: ',outstring
      call validate_descriptor(ldata,outstring,-1,check_rand)
c--------------------------------------------------------------------------------         
 
c  Reset all cursor values to negative numbers.

      do i=0,maxdim_
       ldata%xcursor(i) = -999
       ldata%ycursor(i) = -999
       ldata%zcursor(i) = -999
      enddo

      if (MYID.lt.1) then
         print*,'----------------------------------------------------------'
         print*,'Successfully initialised Panphasia box at level ',lev
         write (6,105) ix_abs,iy_abs,iz_abs
         write (6,106) ix_rel,iy_rel,iz_rel
         write (6,107) ix_per,iy_per,iz_per
         write (6,*)  'Phases used: ',outstring(1:lnblnk(outstring))
         print*,'----------------------------------------------------------'
      endif
 105  format(' Abs origin: (',i12,',',i12,',',i12,')')
 106  format(' Rel origin: (',i12,',',i12,',',i12,')')
 107  format(' Periods   : (',i12,',',i12,',',i12,')') 

c  Set default values 

      ldata%layer_min =  0
      ldata%layer_max =  lev
      ldata%indep_field= 1
      end 
c=================================================================================




c-------------------------------------------------------------------------------
c    The goal of this function is to replace the call to return_gaussian_array in
c    return_cell_props with an equivalent call that returns the Legendre
c    blocks for a special spherically symmetric function defined below.
c-------------------------------------------------------------------------------
      recursive subroutine return_oct_sf_expansion(ldata,ii,lev,x,y,z,ndim,garray)
      use pan_state
      implicit none

      type(state_data), intent(inout) :: ldata

      integer ii,jj
      integer lev,ndim
      real*8 garray(0:ndim-1)
      integer*8 x,y,z 
      real*8 xorig,yorig,zorig
      real*8 cell_data(0:8,0:7) 

c  Debugging variables ....      

      integer*8 xtemp,ytemp,ztemp
      integer ndimension
     
      integer*8 pstore,xstore,ystore,zstore


      integer i,j
c-------------------------------------------------------------------------------
      real*8 length,cube_centre(3),oct_cell_data(0:8,0:7)
      integer*8 lev_range
c-----   First some error checking ------------------------------------------
      if ((lev.lt.0).or.(lev.gt.maxlev_)) stop 'Level out of range! (2)'
      lev_range = 2_dint**lev
      if ((x.lt.0).or.(x.ge.lev_range)) then 
      print*,'x,lev,lev_range',x,lev,lev_range
      call flush(6)
      stop 'x out of range!'
      endif
      if ((y.lt.0).or.(y.ge.lev_range)) then
      print*,'y,lev,lev_range',y,lev,lev_range       
          stop 'y out of range!'
      endif
      if ((z.lt.0).or.(z.ge.lev_range)) stop 'z out of range!' 
c----------------------------------------------------------------------------
c     Define cell centre and cell size and get Legendre block coefficients
c     for an octree function expansion of a single layer of the octree
c----------------------------------------------------------------------------


      length = 1.0/dble(ldata%ix_per_store)*2.0d0**(1+lev-ii) 

      xorig = dble(ldata%ix_abs_store)/2.0d0**(1+lev-ii) 
      yorig = dble(ldata%iy_abs_store)/2.0d0**(1+lev-ii) 
      zorig = dble(ldata%iz_abs_store)/2.0d0**(1+lev-ii) 

      cube_centre(1) =  (dble(x)-xorig+0.5d0)*length
      cube_centre(2) =  (dble(y)-yorig+0.5d0)*length
      cube_centre(3) =  (dble(z)-zorig+0.5d0)*length


c----------------------------------------------------------------------------
       call octree_expansion(cube_centre,length,ndim,garray)
c--------------------------------------------------------------------------------

      return
      end



c-------------------------------------------------------------------------------
c   Expand function of interest in octree basis functions. The
c   result returned is the superposition of the octree functions
c   at a single octree level, expressed as Legendre block
c   functions
c-------------------------------------------------------------------------------    
      recursive subroutine octree_expansion(cube_centre,length,ndim,q)
      implicit none
      real*8 cube_centre(3),length,oct_cell_data(0:8,0:7)

      real*8 local_centre(3), small_data(0:8,0:7),small_len
      real*8 temp_data(0:8)

      real*8 moment(0:7)
      integer ndim


      real*8 p(0:7),q(ndim)


      integer ix,iy,iz,ind1,ind2
      integer i,i1,i2,i3
      integer isign

      small_len = 0.5d0*length

      do i1=0,1
       do i2=0,1
        do i3=0,1
          ind2 = 4*i1 + 2*i2 + i3
          local_centre(1) = cube_centre(1)+0.25d0*dble(2*i1-1)*length
          local_centre(2) = cube_centre(2)+0.25d0*dble(2*i2-1)*length
          local_centre(3) = cube_centre(3)+0.25d0*dble(2*i3-1)*length
          call spherical_perturbation(local_centre,small_len,temp_data)
           do i=0,8
            small_data(i,ind2) = temp_data(i)
          enddo
        enddo
       enddo
      enddo


      call expand_octree_coefficients(small_data,p,q)

      return
      end





      recursive subroutine spherical_perturbation(cube_centre,length,cell_data)
      implicit none
      real*8 cube_centre(3),length,cell_data(0:8)
      integer nfeature, nuse 

      parameter (nfeature=5,nuse=1)
      
      real*8 centre(3), amplitude(nfeature), sigma(nfeature)

      integer i,j
      real*8 cell_data_temp(0:8)
      real*8 pcentre(3),scaled_length
      real*8 CellVolume
      real*8 prefac0,prefac1,prefac2,prefac3

c   Set the parameters of the perturbation.  The periodic volume is
c   a cube of unit length, occupying the positive coordinate octant

      centre(1) = 0.60226666666d0
      centre(2) = 0.4025d0
      centre(3) = 0.5393d0
     
      amplitude(1) = 1.0d0
      sigma(1)     = 0.05d0

      amplitude(2) = 1.5d0
      sigma(2)     = 0.02d0

      amplitude(3) = 0.2d0
      sigma(3)     = 0.002d0

      amplitude(4) = 0.25d0
      sigma(4)     = 0.00024d0

      amplitude(5) = 0.3d0
      sigma(5)     = 0.00003d0

      do i=0,8
       cell_data(i) = 0.0d0
      enddo
      
      do j=1,nuse
       do i=1,3
        pcentre(i) = (cube_centre(i)-centre(i))/sigma(j)
       enddo
        scaled_length = length/sigma(j)

        CellVolume = length**3

 
       call evaluate_3d_integrals(pcentre,scaled_length,cell_data_temp)


c   Scaling factors for change of variables in 3-D integration from length to scaled_length

       prefac0 = amplitude(j)/scaled_length**3
       prefac1 = amplitude(j)/scaled_length**4
       prefac2 = amplitude(j)/scaled_length**5
       prefac3 = amplitude(j)/scaled_length**6

       cell_data(0) = cell_data(0)+prefac0*cell_data_temp(0)*sqrt(CellVolume)  !p000
       cell_data(1) = cell_data(1)+prefac1*cell_data_temp(1)*sqrt(CellVolume)  !p001
       cell_data(2) = cell_data(2)+prefac1*cell_data_temp(2)*sqrt(CellVolume)  !p010
       cell_data(3) = cell_data(3)+prefac2*cell_data_temp(3)*sqrt(CellVolume)  !p011
       cell_data(4) = cell_data(4)+prefac1*cell_data_temp(4)*sqrt(CellVolume)  !p100
       cell_data(5) = cell_data(5)+prefac2*cell_data_temp(5)*sqrt(CellVolume)  !p101
       cell_data(6) = cell_data(6)+prefac2*cell_data_temp(6)*sqrt(CellVolume)  !p110
       cell_data(7) = cell_data(7)+prefac3*cell_data_temp(7)*sqrt(CellVolume)  !p111
 
       cell_data(8) = 0.0d0

      enddo


      return
      end






c---------------------------------------------------------------------------
c  CODE WRITTEN BY ADRIAN JENKINS -   AUGUST 2015
c  ORCiD:    http://orcid.org/0000-0003-4389-2232
c---------------------------------------------------------------------------
      recursive subroutine evaluate_3d_integrals(pos_cen,len,cell_data)
c---------------------------------------------------------------------------
c  GOAL
c---------------------------------------------------------------------------
c
c
c  To expand the overdensity: 
c      
c            rho = (3-r^2)exp(-r^2/2) 
c
c   In terms of Legendre basis functions (Jenkins 2013)
c
c
c
c  This overdensity function is taken from Jenkins 2010.
c
c  It is easy to compute the Zeldovich and 2lpt displacements generated
c  by this function.  It can be used to test initial condition 
c  generator codes.
c
c  The coefficients are computed by a 3-d integral over a cell volume. 
c  However the integral can be written as a sum of products of 1-d integrals over x,y or z
c  coordinates, and the 1-d integrals can all be expressed in terms
c  of incomplete Gamma functions
c
c---------------------------------------------------------------------------
      implicit none
      real*8 pos_cen(3),len,cell_data(0:8)

      real*8 pos_min(3),pos_max(3)
      real*8 abs_u_min(3),abs_u_max(3)
      real*8 a(0:3),c(0:3),s(0:3),p_min,p_max
      real*8 stretch
      real*8 gammp

      real*8 si(0:3,3),ti(0:3,3)

      real*8 Coeff000,Coeff001,Coeff010,Coeff011
      real*8 Coeff100,Coeff101,Coeff110,Coeff111

      integer i,j,n

c---------------------------------------------------------------------------
c  SET UP ALL COEFFICIENTS
c---------------------------------------------------------------------------
c   Normalising coefficients for each integral
c
c   Each coefficient is the product of two numbers
c   (i) a coefficient from a Legendre block function (Jenkins 2013)
c       for c(0) and c(2) this is unity, for c(1) and c(3) this
c       is the sqrt(3)
c       
c   (ii) A coefficient of 2**( (n-1)/2)) \Gamma[ (n+1)/2]
c        for the integral  \int x^n exp(-x^2/2) dx  - 
c        which comes from the substitution u = x^2/2
c        and the definition of a incomplete Gamma function
c
c      P(a,x) = 1/\Gamma(a) * \int_0^t exp(-t) t^{a-1} dt
c
c
c---------------------------------------------------------------------------
c
 

      c(0) = sqrt(3.1415926535897932d0/2.0d0)
      c(1) = 1.0d0
      c(2) = c(0)
      c(3) = 2.0d0

c   Define a(n) = (n+1)/2

      a(0) = 0.5d0
      a(1) = 1.0d0
      a(2) = 1.5d0
      a(3) = 2.0d0

c   The substitution used to convert the desired integral into the incomplete
c   Gamma function (shown below) 'looses' the signs of the limits.  The n=0 and 2 the
c   desired integrand is symmetric about the origin, while for n=1 and 3,
c   it is antisymmetric.  The array 's' below encodes this information so
c   that the definite integral is evaluated correctly for positive or
c   negative values in pos_min and pos_max.  The fortran sign function
c   is used to carry the sign of the pos_min and pos_max arguments.

      s(0) = -1.0d0
      s(1) = 1.0d0
      s(2) = -1.0d0
      s(3) = 1.0d0

c----------------------------------------------------------------------------
c    Change of variable - from substitution in the integral above  (u=x^2/2)

      stretch = pos_cen(1)**2 + pos_cen(2)**2 + pos_cen(3)**2


      do i=1,3
        pos_min(i) = pos_cen(i) - 0.5d0*len
        pos_max(i) = pos_cen(i) + 0.5d0*len

        abs_u_min(i) = 0.5d0 * pos_min(i)**2
        abs_u_max(i) = 0.5d0 * pos_max(i)**2
      enddo
c----------------------------------------------------------------------------



c----------------------------------------------------------------------------
c    Create a 4x3 matrix of integrals  
c    First index n, second index coordinate (1=x,2=y,3=z)
c
c    si(n,j) = \int p^n exp(-p*p/2) dp     (for j=1,2,3, p = x,y,z)
c                                               n=0,1,2,3
c----------------------------------------------------------------------------

      do n=0,3  
       do j=1,3

        if (pos_min(j).lt.0.0) then
           p_min = s(n)
        else 
           p_min = 1.0d0
        endif

        if (pos_max(j).lt.0.0) then
           p_max = s(n)
        else 
           p_max = 1.0d0
        endif


       if (stretch.lt.200.0d0) then
     
       si(n,j) = c(n)*
     &  ( p_max*gammp(a(n),abs_u_max(j)) 
     &   -p_min*gammp(a(n),abs_u_min(j))) 
   
       else
         si(n,j) = 0.0d0
       endif



       ti(n,j) = si(n,j)
       enddo
      enddo
c----------------------------------------------------------------------------




c----------------------------------------------------------------------------
c     Compute integrals with respect to the Legendre block. Each block
c     factorises into x,y and z directions
c     p_j1j2j3    -  where j1,j2,j3 are either zero or 1
c
c     p_0(x)  = 1                   Zeroth moment
c     p_1(x)  = sqrt(12)(x-xcen)    First moment
c
c----------------------------------------------------------------------------c

      do n=1,3,2  ! Shift origins to cell centre for first moments
       do j=1,3
          ti(n,j) = sqrt(12.0d0) * (ti(n,j) - pos_cen(j)*ti(n-1,j))
       enddo
      enddo



c----------------------------------------------------------------------------
c   Combine the computed integrals to give the 8 Legendre block
c   expansion coefficients. 
c----------------------------------------------------------------------------


      Coeff000 = 3.0 * ti(0,1) * ti(0,2) * ti(0,3) 
     &  - ti(2,1) * ti(0,2) * ti(0,3)
     &  - ti(0,1) * ti(2,2) * ti(0,3)
     &  - ti(0,1) * ti(0,2) * ti(2,3)

 
      Coeff100 = 3.0 * ti(1,1) * ti(0,2) * ti(0,3) 
     &  - ti(3,1) * ti(0,2) * ti(0,3)
     &  - ti(1,1) * ti(2,2) * ti(0,3)
     &  - ti(1,1) * ti(0,2) * ti(2,3)

      Coeff010 = 3.0 * ti(0,1) * ti(1,2) * ti(0,3) 
     &  - ti(2,1) * ti(1,2) * ti(0,3)
     &  - ti(0,1) * ti(3,2) * ti(0,3)
     &  - ti(0,1) * ti(1,2) * ti(2,3)


      Coeff001 = 3.0 * ti(0,1) * ti(0,2) * ti(1,3) 
     &  - ti(2,1) * ti(0,2) * ti(1,3)
     &  - ti(0,1) * ti(2,2) * ti(1,3)
     &  - ti(0,1) * ti(0,2) * ti(3,3)


      Coeff110 = 3.0 * ti(1,1) * ti(1,2) * ti(0,3) 
     &  - ti(3,1) * ti(1,2) * ti(0,3)
     &  - ti(1,1) * ti(3,2) * ti(0,3)
     &  - ti(1,1) * ti(1,2) * ti(2,3)

      Coeff101 = 3.0 * ti(1,1) * ti(0,2) * ti(1,3) 
     &  - ti(3,1) * ti(0,2) * ti(1,3)
     &  - ti(1,1) * ti(2,2) * ti(1,3)
     &  - ti(1,1) * ti(0,2) * ti(3,3)


      Coeff011 = 3.0 * ti(0,1) * ti(1,2) * ti(1,3) 
     &  - ti(2,1) * ti(1,2) * ti(1,3)
     &  - ti(0,1) * ti(3,2) * ti(1,3)
     &  - ti(0,1) * ti(1,2) * ti(3,3)



      Coeff111 = 3.0 * ti(1,1) * ti(1,2) * ti(1,3) 
     &  - ti(3,1) * ti(1,2) * ti(1,3)
     &  - ti(1,1) * ti(3,2) * ti(1,3)
     &  - ti(1,1) * ti(1,2) * ti(3,3)     

c--------------------------------------------------------------------------
c     Copy into output structure - ordering matches the Panphasia code
c     (see Jenkins & Booth 2013)
c--------------------------------------------------------------------------


      cell_data(0) = Coeff000     ! Scales as len**1.5
      cell_data(1) = Coeff001     ! Scales as len**3.5
      cell_data(2) = Coeff010
      cell_data(3) = Coeff011     ! Scales as len**5.5
      cell_data(4) = Coeff100
      cell_data(5) = Coeff101
      cell_data(6) = Coeff110
      cell_data(7) = Coeff111     ! Scales as len**7.5

      cell_data(8) = 0.0d0      ! Set the 'independent' field to zero (J&B13)



      return
      end


c=========================================================================
c   NUMERICAL RECIPES ROUTINES BELOW - modified to make the
c   output double precision.  Routines taken from the Blue f77 Book
c   Value of EPS changed from 3e-7 to 3e-15, ITMAX increased from 100 to 200
c=========================================================================

      REAL*8 recursive FUNCTION gammp(a,x)
      REAL*8 a,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)stop 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END

      recursive SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=200,EPS=3.d-15,FPMIN=1.d-290)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      stop 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END

      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=200,EPS=3.d-15)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)stop 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      stop 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END


      REAL*8 recursive FUNCTION gammln(xx)
      REAL*8 xx
      INTEGER j
      REAL*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
c======================  END NR ================================================
c===============================================================================
      recursive subroutine expand_octree_coefficients(cell_data, p,q)
      implicit none
c----------------- Define subroutine arguments -----------------------------------
      real*8 cell_data(0:8,0:7)
      real*8 p(0:7),q(56)

c----------------- Define constants using notation from appendix A of Jenkins 2013
 
      real*8 a1,a2,b1,b2,b3,c1,c2,c3,c4,rooteighth_factor

      parameter(a1 = 0.5d0*sqrt(3.0d0),      a2 = 0.5d0)

      parameter(b1 = 0.75d0,                 b2 = 0.25d0*sqrt(3.0d0))
      parameter(b3 = 0.25d0)

      parameter(c1 = sqrt(27.0d0/64.0d0),    c2 = 0.375d0)
      parameter(c3 = sqrt(3.0d0/64.0d0),     c4 = 0.125d0)

      parameter(rooteighth_factor = sqrt(0.125d0))

c----------------- Define octree variables --------------------------------

      real*8 po(0:7,0:7),tsum(0:7,0:7)
      integer iparity
      integer i,j,ix,iy,iz
      integer icx,icy,icz
      integer iox,ioy,ioz
      real*8 parity,isig
c-----------------------------------------------------------------------------

c 
c
c   We now calculate the Legendre basis coefficients for all eight child cells
c   by applying the appropriate reflectional parities to the coefficients 
c   calculated above for the positive octant child cell.
c
c   See equations A2 and A3 in appendix A of Jenkins 2013.
c
c   The reflectional parity is given by (ix,iy,iz) loops below.
c
c   The (icx,icy,icz) loops below, loop over the eight child cells.
c
c   The positive octant child cell is given below by  (icx=icy=icz=0) or i=7.
c
c   The combination ix*icx +iy*icy +iz*icz is either even or odd, depending
c   on whether the parity change is even or odd.
c
c   The variables iox,ioy,ioz are used to loop over the different
c   types of Legendre basis function.
c
c   The combination iox*icx + ioy*icy + ioz*icz is either even and odd
c   and identifies which coefficients keep or change sign respectively
c   due to a pure reflection about the principal planes.
c--------------------------------------------------------------------------------------------
      do i=0,7
       p(i) = -9999.0d0
      enddo

      do i=1,56
       q(i) = -999.0d0
      enddo


      do iz=0,7
       do iy=0,7
        po(iy,iz) = 0.0d0      ! Set positive octant coefficients to zero
       enddo
      enddo
c--------------------------------------------------------------------------------------------
      do iz=0,1              ! Loop over z parity (0=keep sign, 1=change sign)
       do iy=0,1             ! Loop over y parity (0=keep sign, 1=change sign)
        do ix=0,1            ! Loop over x parity (0=keep sign, 1=change sign)
        iparity = 4*ix + 2*iy + iz

         do icx=0,1                      ! Loop over x-child cells
          do icy=0,1                     ! Loop over y-child cells
           do icz=0,1                    ! Loop over z-child cells

             if (mod(ix*icx+iy*icy+iz*icz,2).eq.0) then
                  parity = 1.0d0
             else
                  parity =-1.0d0
             endif

             i = 7 - 4*icx -2*icy - icz               ! Calculate which child cell this is.


             do iox=0,1                               ! Loop over Legendre basis function type                     
              do ioy=0,1                              ! Loop over Legendre basis function type
               do ioz=0,1                             ! Loop over Legendre basis function type

                  j = 4*iox + 2*ioy + ioz

                  if (mod(iox*icx + ioy*icy + ioz*icz,2).eq.0) then
                       isig =  parity
                  else
                       isig = -parity
                  endif

                  po(j,iparity) = po(j,iparity) + isig*cell_data(j,i)*rooteighth_factor

               enddo
              enddo
             enddo

           enddo   
          enddo
         enddo

        enddo
       enddo
      enddo


c
c    The calculations immediately below evalute the eight Legendre block coefficients for the
c    child cell that is furthest from the absolute coordiate origin of the octree - we call
c    this the positive octant cell.
c
c    The coefficients are given by a set of matrix equations which combine the
c    coefficients of the Legendre basis functions of the parent cell itself, with
c    the coefficients from the octree basis functions that occupy the
c    parent cell.   
c
c    The Legendre basis function coefficients of the parent cell are stored in
c    the variables, p(0) - p(7) and are initialise above.
c
c    The coefficients of the octree basis functions are determined by the
c    first 56 entries of the array g, which is passed down into this
c    subroutine.
c
c    These two sources of information are combined using a set of linear equations.
c    The coefficients of these linear equations are taken from the inverses or
c    equivalently transposes of the matrices given in appendix A of Jenkins 2013.
c    The matrices in appendix A define the PANPHASIA octree basis functions
c    in terms of Legendre blocks.
c
c    All of the Legendre block functions of the parent cell, and the octree basis
c    functions of the parent cell share one of eight distinct symmetries with respect to
c    reflection about the x1=0,x2=0,x3=0 planes (where the origin is taken as the parent 
c    cell centre and x1,x2,x3 are parallel to the cell edges).
c
c    Each function has either purely reflectional symmetry (even parity) or
c    reflectional symmetry with a sign change (odd parity) about each of the three principal
c    planes through the cell centre. There are therefore 8 parity types. We can label each 
c    parity type with a binary triplet. So 000 is pure reflectional symmetry about 
c    all of the principal planes.
c  
c    In the code below the parent cell Legendre block functions, and octree functions are 
c    organised into eight groups each with eight members. Each group has a common
c    parity type.
c
c    We keep the contributions of each parity type to each of the eight Legendre basis
c    functions occupying the positive octant cell separate. Once they have all been
c    computed, we can apply the different symmetry operations and determine the
c    Legendre block basis functions for all eight child cells at the same time.
c---------------------------------------------------------------------------------------
c    000/ 0-parity

      p(0) =  1.0d0*po(0,0)
      q(1) = -1.0d0*po(1,0)
      q(2) = -1.0d0*po(2,0)
      q(3) =  1.0d0*po(3,0)
      q(4) = -1.0d0*po(4,0)
      q(5) =  1.0d0*po(5,0)
      q(6) =  1.0d0*po(6,0)
      q(7) = -1.0d0*po(7,0)

c    100/ 4-parity

      p(4)  =  a1*po(0,4) +  a2*po(4,4)
      q(8)  = -a2*po(0,4)  + a1*po(4,4)
      q(9)  =  po(1,4)
      q(11) = -po(3,4)
      q(10) =  po(2,4)
      q(12) = -po(5,4)
      q(13) = -po(6,4)
      q(14) =  po(7,4)

c     010/ 2-parity

      p(2)  =  a1*po(0,2) + a2*po(2,2)
      q(15) = -a2*po(0,2) + a1*po(2,2) 
      q(16) =  po(1,2) 
      q(17) = -po(3,2)
      q(18) =  po(4,2)
      q(19) = -po(5,2)
      q(20) = -po(6,2)
      q(21) =  po(7,2)


c     001/ 1-parity

      p(1)  =  a1*po(0,1) + a2*po(1,1)
      q(22) = -a2*po(0,1) + a1*po(1,1)
      q(23) =  po(2,1)
      q(24) = -po(3,1)
      q(25) =  po(4,1)
      q(26) = -po(5,1)
      q(27) = -po(6,1)
      q(28) =  po(7,1)

c    110/ 6-parity

      p(6)  =  b1*po(0,6) + b2*po(2,6) + b2*po(4,6) + b3*po(6,6)
      q(29) = -b2*po(0,6) - b3*po(2,6) + b1*po(4,6) + b2*po(6,6)
      q(30) =  b3*po(0,6) - b2*po(2,6) + b2*po(4,6) - b1*po(6,6)
      q(31) = -b2*po(0,6) + b1*po(2,6) + b3*po(4,6) - b2*po(6,6)
      q(32) = -po(1,6)
      q(33) =  po(3,6)
      q(34) =  po(5,6)
      q(35) = -po(7,6)


c     011/ 3-parity

      p(3)  =  b1*po(0,3) + b2*po(1,3) + b2*po(2,3) + b3*po(3,3)
      q(36) = -b2*po(0,3) - b3*po(1,3) + b1*po(2,3) + b2*po(3,3)
      q(37) =  b3*po(0,3) - b2*po(1,3) + b2*po(2,3) - b1*po(3,3)
      q(38) = -b2*po(0,3) + b1*po(1,3) + b3*po(2,3) - b2*po(3,3)
      q(39) = -po(4,3)
      q(40) =  po(5,3)
      q(41) =  po(6,3)
      q(42) = -po(7,3)

c     101/ 5-parity


      p(5)  =  b1*po(0,5) + b2*po(1,5) + b2*po(4,5) + b3*po(5,5)
      q(43) = -b2*po(0,5) - b3*po(1,5) + b1*po(4,5) + b2*po(5,5)
      q(44) =  b3*po(0,5) - b2*po(1,5) + b2*po(4,5) - b1*po(5,5)
      q(45) = -b2*po(0,5) + b1*po(1,5) + b3*po(4,5) - b2*po(5,5)
      q(46) = -po(2,5)
      q(47) =  po(3,5)
      q(48) =  po(6,5)
      q(49) = -po(7,5)

c     111/ 7-parity

      p(7) = c1*po(0,7)+c2*po(1,7)+c2*po(2,7)+c3*po(3,7)+c2*po(4,7)+c3*po(5,7)+c3*po(6,7)+c4*po(7,7)
      q(50)=-c2*po(0,7)+c1*po(1,7)+c2*po(2,7)-c3*po(3,7)-c2*po(4,7)+c3*po(5,7)+c4*po(6,7)-c3*po(7,7)
      q(51)=-c2*po(0,7)-c2*po(1,7)+c1*po(2,7)-c3*po(3,7)+c2*po(4,7)-c4*po(5,7)+c3*po(6,7)+c3*po(7,7)
      q(52)=-c2*po(0,7)+c2*po(1,7)-c2*po(2,7)+c4*po(3,7)+c1*po(4,7)-c3*po(5,7)+c3*po(6,7)-c3*po(7,7)
      q(53)= c3*po(0,7)-c3*po(1,7)-c3*po(2,7)-c1*po(3,7)+c4*po(4,7)+c2*po(5,7)+c2*po(6,7)-c2*po(7,7)
      q(54)= c3*po(0,7)+c3*po(1,7)-c4*po(2,7)-c2*po(3,7)-c3*po(4,7)-c1*po(5,7)+c2*po(6,7)+c2*po(7,7)
      q(55)= c3*po(0,7)+c4*po(1,7)+c3*po(2,7)-c2*po(3,7)+c3*po(4,7)-c2*po(5,7)-c1*po(6,7)-c2*po(7,7)
      q(56)=-c4*po(0,7)+c3*po(1,7)-c3*po(2,7)-c2*po(3,7)+c3*po(4,7)+c2*po(5,7)-c2*po(6,7)+c1*po(7,7)


      return
      end
c===============================================================================


      recursive subroutine compound_octree_coefficients(p,q,cell_data)
      implicit none
c----------------- Define subroutine arguments -----------------------------------
      real*8 cell_data(0:8,0:7)
      real*8 p(0:7),q(56)

c----------------- Define constants using notation from appendix A of Jenkins 2013
 
      real*8 a1,a2,b1,b2,b3,c1,c2,c3,c4,rooteighth_factor

      parameter(a1 = 0.5d0*sqrt(3.0d0),      a2 = 0.5d0)

      parameter(b1 = 0.75d0,                 b2 = 0.25d0*sqrt(3.0d0))
      parameter(b3 = 0.25d0)

      parameter(c1 = sqrt(27.0d0/64.0d0),    c2 = 0.375d0)
      parameter(c3 = sqrt(3.0d0/64.0d0),     c4 = 0.125d0)

      parameter(rooteighth_factor = sqrt(0.125d0))

c----------------- Define octree variables --------------------------------

      real*8 po(0:7,0:7),tsum(0:7,0:7)
      integer iparity
      integer i,j,ix,iy,iz
      integer icx,icy,icz
      integer iox,ioy,ioz
      real*8 parity,isig
c-----------------------------------------------------------------------------
c
c
c    The calculations immediately below evalute the eight Legendre block coefficients for the
c    child cell that is furthest from the absolute coordiate origin of the octree - we call
c    this the positive octant cell.
c
c    The coefficients are given by a set of matrix equations which combine the
c    coefficients of the Legendre basis functions of the parent cell itself, with
c    the coefficients from the octree basis functions that occupy the
c    parent cell.   
c
c    The Legendre basis function coefficients of the parent cell are stored in
c    the variables, p(0) - p(7) and are initialise above.
c
c    The coefficients of the octree basis functions are determined by the
c    first 56 entries of the array g, which is passed down into this
c    subroutine.
c
c    These two sources of information are combined using a set of linear equations.
c    The coefficients of these linear equations are taken from the inverses or
c    equivalently transposes of the matrices given in appendix A of Jenkins 2013.
c    The matrices in appendix A define the PANPHASIA octree basis functions
c    in terms of Legendre blocks.
c
c    All of the Legendre block functions of the parent cell, and the octree basis
c    functions of the parent cell share one of eight distinct symmetries with respect to
c    reflection about the x1=0,x2=0,x3=0 planes (where the origin is taken as the parent 
c    cell centre and x1,x2,x3 are parallel to the cell edges).
c
c    Each function has either purely reflectional symmetry (even parity) or
c    reflectional symmetry with a sign change (odd parity) about each of the three principal
c    planes through the cell centre. There are therefore 8 parity types. We can label each 
c    parity type with a binary triplet. So 000 is pure reflectional symmetry about 
c    all of the principal planes.
c  
c    In the code below the parent cell Legendre block functions, and octree functions are 
c    organised into eight groups each with eight members. Each group has a common
c    parity type.
c
c    We keep the contributions of each parity type to each of the eight Legendre basis
c    functions occupying the positive octant cell separate. Once they have all been
c    computed, we can apply the different symmetry operations and determine the
c    Legendre block basis functions for all eight child cells at the same time.
c---------------------------------------------------------------------------------------
c    000/ 0-parity

      po(0,0) =  1.0d0*p(0)
      po(1,0) = -1.0d0*q(1)
      po(2,0) = -1.0d0*q(2)
      po(3,0) =  1.0d0*q(3)
      po(4,0) = -1.0d0*q(4)
      po(5,0) =  1.0d0*q(5)
      po(6,0) =  1.0d0*q(6)
      po(7,0) = -1.0d0*q(7)

c    100/ 4-parity

      po(0,4) =  a1*p(4)  - a2*q(8)
      po(1,4) =  q(9)
      po(2,4) =  q(10)
      po(3,4) = -q(11)
      po(4,4) =  a2*p(4)  + a1*q(8)
      po(5,4) = -q(12) 
      po(6,4) = -q(13)
      po(7,4) =  q(14)

c     010/ 2-parity

      po(0,2) =  a1*p(2) - a2*q(15)
      po(1,2) =  q(16) 
      po(2,2) =  a2*p(2) + a1*q(15) 
      po(3,2) = -q(17)
      po(4,2) =  q(18)
      po(5,2) = -q(19)
      po(6,2) = -q(20)
      po(7,2) =  q(21)


c     001/ 1-parity

      po(0,1) =  a1*p(1) - a2*q(22)
      po(1,1) =  a2*p(1) + a1*q(22)
      po(2,1) =  q(23)
      po(3,1) = -q(24)
      po(4,1) =  q(25)
      po(5,1) = -q(26)
      po(6,1) = -q(27)
      po(7,1) =  q(28)

c    110/ 6-parity

      po(0,6) = b1*p(6) - b2*q(29) + b3*q(30) - b2*q(31)
      po(1,6) = -q(32)
      po(2,6) = b2*p(6) - b3*q(29) - b2*q(30) + b1*q(31)
      po(3,6) =  q(33)
      po(4,6) = b2*p(6) + b1*q(29) + b2*q(30) + b3*q(31)
      po(5,6) =  q(34)
      po(6,6) = b3*p(6) + b2*q(29) - b1*q(30) - b2*q(31)
      po(7,6) = -q(35) 


c     011/ 3-parity

      po(0,3) = b1*p(3) - b2*q(36) + b3*q(37) - b2*q(38)
      po(1,3) = b2*p(3) - b3*q(36) - b2*q(37) + b1*q(38)
      po(2,3) = b2*p(3) + b1*q(36) + b2*q(37) + b3*q(38)
      po(3,3) = b3*p(3) + b2*q(36) - b1*q(37) - b2*q(38)
      po(4,3) = -q(39) 
      po(5,3) =  q(40)
      po(6,3) =  q(41)
      po(7,3) = -q(42)

c     101/ 5-parity

      po(0,5) = b1*p(5) - b2*q(43) + b3*q(44) - b2*q(45)
      po(1,5) = b2*p(5) - b3*q(43) - b2*q(44) + b1*q(45) 
      po(2,5) = -q(46) 
      po(3,5) =  q(47)
      po(4,5) = b2*p(5) + b1*q(43) + b2*q(44) + b3*q(45)
      po(5,5) = b3*p(5) + b2*q(43) - b1*q(44) - b2*q(45)
      po(6,5) =  q(48)
      po(7,5) = -q(49)

c     111/ 7-parity

      po(0,7) = c1*p(7) - c2*q(50) - c2*q(51) - c2*q(52) + c3*q(53) + c3*q(54) + c3*q(55) - c4*q(56)
      po(1,7) = c2*p(7) + c1*q(50) - c2*q(51) + c2*q(52) - c3*q(53) + c3*q(54) + c4*q(55) + c3*q(56) 
      po(2,7) = c2*p(7) + c2*q(50) + c1*q(51) - c2*q(52) - c3*q(53) - c4*q(54) + c3*q(55) - c3*q(56)
      po(3,7) = c3*p(7) - c3*q(50) - c3*q(51) + c4*q(52) - c1*q(53) - c2*q(54) - c2*q(55) - c2*q(56) 
      po(4,7) = c2*p(7) - c2*q(50) + c2*q(51) + c1*q(52) + c4*q(53) - c3*q(54) + c3*q(55) + c3*q(56)
      po(5,7) = c3*p(7) + c3*q(50) - c4*q(51) - c3*q(52) + c2*q(53) - c1*q(54) - c2*q(55) + c2*q(56)
      po(6,7) = c3*p(7) + c4*q(50) + c3*q(51) + c3*q(52) + c2*q(53) + c2*q(54) - c1*q(55) - c2*q(56)
      po(7,7) = c4*p(7) - c3*q(50) + c3*q(51) - c3*q(52) - c2*q(53) + c2*q(54) - c2*q(55) + c1*q(56)
c--------------------------------------------------------------------------------------------
c 
c
c   We now calculate the Legendre basis coefficients for all eight child cells
c   by applying the appropriate reflectional parities to the coefficients 
c   calculated above for the positive octant child cell.
c
c   See equations A2 and A3 in appendix A of Jenkins 2013.
c
c   The reflectional parity is given by (ix,iy,iz) loops below.
c
c   The (icx,icy,icz) loops below, loop over the eight child cells.
c
c   The positive octant child cell is given below by  (icx=icy=icz=0) or i=7.
c
c   The combination ix*icx +iy*icy +iz*icz is either even or odd, depending
c   on whether the parity change is even or odd.
c
c   The variables iox,ioy,ioz are used to loop over the different
c   types of Legendre basis function.
c
c   The combination iox*icx + ioy*icy + ioz*icz is either even and odd
c   and identifies which coefficients keep or change sign respectively
c   due to a pure reflection about the principal planes.
c--------------------------------------------------------------------------------------------

      do iz=0,7
       do iy=0,7
        tsum(iy,iz) = 0.0d0      ! Zero temporary sums
       enddo
      enddo
c--------------------------------------------------------------------------------------------
      do iz=0,1              ! Loop over z parity (0=keep sign, 1=change sign)
       do iy=0,1             ! Loop over y parity (0=keep sign, 1=change sign)
        do ix=0,1            ! Loop over x parity (0=keep sign, 1=change sign)
        iparity = 4*ix + 2*iy + iz

         do icx=0,1                      ! Loop over x-child cells
          do icy=0,1                     ! Loop over y-child cells
           do icz=0,1                    ! Loop over z-child cells

             if (mod(ix*icx+iy*icy+iz*icz,2).eq.0) then
                  parity = 1.0d0
             else
                  parity =-1.0d0
             endif

             i = 7 - 4*icx -2*icy - icz               ! Calculate which child cell this is.


             do iox=0,1                               ! Loop over Legendre basis function type                     
              do ioy=0,1                              ! Loop over Legendre basis function type
               do ioz=0,1                             ! Loop over Legendre basis function type

                  j = 4*iox + 2*ioy + ioz

                  if (mod(iox*icx + ioy*icy + ioz*icz,2).eq.0) then
                       isig =  parity
                  else
                       isig = -parity
                  endif

                  tsum(j,i) = tsum(j,i) + isig*po(j,iparity)

               enddo
              enddo
             enddo

           enddo   
          enddo
         enddo

        enddo
       enddo
      enddo


c   Assign values of the output variables an set independent field to zero

      do i=0,7
       do j=0,7
         cell_data(j,i) = tsum(j,i)*rooteighth_factor
        enddo
         cell_data(8,i) = 0.0d0
      enddo


      return
      end




