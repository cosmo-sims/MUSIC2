!=====================================================================================c
!        
! The code below was written by: Stephen Booth
!                                Edinburgh Parallel Computing Centre
!                                The University of Edinburgh
!                                JCMB
!                                Mayfield Road
!                                Edinburgh EH9 3JZ
!                                United Kingdom
!
! This file is part of the software made public in
! Jenkins and Booth 2013  - arXiv:1306.XXXX
!
! The software computes the Panphasia Gaussian white noise field
! realisation described in detail in Jenkins 2013 - arXiv:1306.XXXX
! 
!
!
! This software is free, subject to a agreeing licence conditions:
!
!
! (i)  you will publish the phase descriptors and reference Jenkins (13) 
!      for any new simulations that use Panphasia phases. You will pass on this 
!      condition to others for any software or data you make available publically 
!      or privately that makes use of Panphasia. 
!
! (ii) that you will ensure any publications using results derived from Panphasia 
!      will be submitted as a final version to arXiv prior to or coincident with
!      publication in a journal. 
!
!
! (iii) that you report any bugs in this software as soon as confirmed to 
!       A.R.Jenkins@durham.ac.uk 
!
! (iv)  that you understand that this software comes with no warranty and that is 
!       your responsibility to ensure that it is suitable for the purpose that 
!       you intend. 
!
!=====================================================================================c
!{{{Rand_base (define kind types) 
MODULE Rand_base
! This module just declares the base types 
! we may have to edit this to match to the target machine
! we really need a power of 2 selected int kind in fortran-95 we could
! do this with a PURE function I think.

!
! 10 decimal digits will hold 2^31
!

   INTEGER, PARAMETER :: Sint = SELECTED_INT_KIND(9)
!  INTEGER, PARAMETER :: Sint = SELECTED_INT_KIND(10)
!  INTEGER, PARAMETER :: Sint = 4

!
! 18-19 decimal digits will hold 2^63
! but all 19 digit numbers require 2^65 :-(
!

   INTEGER, PARAMETER :: Dint = SELECTED_INT_KIND(17)
!  INTEGER, PARAMETER :: Dint = SELECTED_INT_KIND(18)
!  INTEGER, PARAMETER :: Dint = 8

! type for index counters must hold Nstore
  INTEGER, PARAMETER :: Ctype = SELECTED_INT_KIND(3)
END MODULE Rand_base
!}}}

!{{{Rand_int (random integers mod 2^31-1) 

MODULE Rand_int
  USE Rand_base
  IMPLICIT NONE
! The general approach of this module is two have
! two types Sint and Dint 
! 
! Sint should have at least 31 bits
! dint shouldhave at least 63

!{{{constants

  INTEGER(KIND=Ctype), PARAMETER :: Nstate=5_Ctype
  INTEGER(KIND=Ctype), PRIVATE, PARAMETER :: Nbatch=128_Ctype
  INTEGER(KIND=Ctype), PRIVATE, PARAMETER :: Nstore=Nstate+Nbatch

  INTEGER(KIND=Sint), PRIVATE, PARAMETER  :: M = 2147483647_Sint
  INTEGER(KIND=Dint), PRIVATE, PARAMETER  :: Mask = 2147483647_Dint
  INTEGER(KIND=Dint), PRIVATE, PARAMETER  :: A1 = 107374182_Dint
  INTEGER(KIND=Dint), PRIVATE, PARAMETER  :: A5 = 104480_Dint
  LOGICAL, PARAMETER :: Can_step_int=.TRUE.
  LOGICAL, PARAMETER :: Can_reverse_int=.TRUE.

!}}}

!{{{Types
!
! This type holds the state of the generator
!
!{{{TYPE RAND_state

TYPE RAND_state
  PRIVATE
  INTEGER(KIND=Sint) :: state(Nstore) 
! do we need to re-fill state table this is reset when we initialise state.
  LOGICAL :: need_fill 
! position of the next state variable to output
  INTEGER(KIND=Ctype) :: pos
END TYPE RAND_state

!}}}

!
! This type defines the offset type used for stepping.
!
!{{{TYPE RAND_offset

TYPE RAND_offset
  PRIVATE
  INTEGER(KIND=Sint) :: poly(Nstate)
END TYPE RAND_offset

!}}}

!}}}

!{{{interface and overloads
!
! Allow automatic conversion between integers and offsets
!
INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE Rand_set_offset
  MODULE PROCEDURE Rand_load
  MODULE PROCEDURE Rand_save
  MODULE PROCEDURE Rand_seed
END INTERFACE
INTERFACE OPERATOR(+)
  MODULE PROCEDURE Rand_add_offset
END INTERFACE
INTERFACE OPERATOR(*)
  MODULE PROCEDURE Rand_mul_offset
END INTERFACE

!
! overload + as the boost/stepping operator
!
INTERFACE OPERATOR(+)
  MODULE PROCEDURE Rand_step
  MODULE PROCEDURE Rand_boost
END INTERFACE
!}}}


!{{{PUBLIC/PRIVATE 
  PRIVATE reduce,mod_saxpy,mod_sdot,p_saxpy,p_sdot,poly_mult
  PRIVATE poly_square, poly_power
  PRIVATE fill_state, repack_state

  PUBLIC Rand_sint, Rand_sint_vec

  PUBLIC Rand_save, Rand_load
  PUBLIC Rand_set_offset, Rand_add_offset, Rand_mul_offset
  PUBLIC Rand_step, Rand_boost, Rand_seed
!}}}

CONTAINS
  !{{{Internals
  !{{{RECURSIVE FUNCTION reduce(A)
  RECURSIVE FUNCTION reduce(A)
  !
  ! Take A Dint and reduce to Sint MOD M
  !
   INTEGER(KIND=Dint), INTENT(IN) :: A
   INTEGER(KIND=Sint) reduce
   INTEGER(KIND=Dint) tmp
  
    tmp = A  
    DO WHILE( ISHFT(tmp, -31) .GT. 0 )
      tmp = IAND(tmp,Mask) + ISHFT(tmp, -31)
    END DO
    IF( tmp .GE. M ) THEN
      reduce = tmp - M
    ELSE
      reduce = tmp
    END IF
  END FUNCTION reduce
  !}}}
  !{{{RECURSIVE SUBROUTINE fill_state(x)
  RECURSIVE SUBROUTINE fill_state(x)
  TYPE(RAND_state), INTENT(INOUT) ::  x
  INTEGER(KIND=Ctype) i
  INTRINSIC IAND, ISHFT
  INTEGER(KIND=Dint)  tmp
    DO i=Nstate+1,Nstore
      tmp = (x%state(i-5) * A5) + (x%state(i-1)*A1)
      !
      ! now reduce down to mod M efficiently
      ! really hope the compiler in-lines this
      !
      ! x%state(i) = reduce(tmp)
      DO WHILE( ISHFT(tmp, -31) .GT. 0 )
        tmp = IAND(tmp,Mask) + ISHFT(tmp, -31)
      END DO
      IF( tmp .GE. M ) THEN
        x%state(i) = tmp - M
      ELSE
        x%state(i) = tmp
      END IF
  
    END DO
    x%need_fill = .FALSE.
  END SUBROUTINE fill_state
  !}}}
  !{{{RECURSIVE SUBROUTINE repack_state(x)
  RECURSIVE SUBROUTINE repack_state(x)
  TYPE(RAND_state), INTENT(INOUT) ::  x
  INTEGER(KIND=Ctype) i
    DO i=1,Nstate
      x%state(i) = x%state(i+x%pos-(Nstate+1))
    END DO
    x%pos = Nstate + 1
    x%need_fill = .TRUE.  
  END SUBROUTINE repack_state
  !}}}
  !{{{RECURSIVE SUBROUTINE mod_saxpy(y,a,x)
  RECURSIVE SUBROUTINE mod_saxpy(y,a,x)
   INTEGER(KIND=Ctype) i
   INTEGER(KIND=Sint) y(Nstate)
   INTEGER(KIND=Sint) a
   INTEGER(KIND=Sint) x(Nstate)
   INTEGER(KIND=Dint) tx,ty,ta
  
     IF( a .EQ. 0_Sint ) RETURN
  
     ! We use KIND=Dint temporaries here to ensure
     ! that we don't overflow in the expression
  
     ta = a
     DO i=1,Nstate
       ty=y(i)
       tx=x(i)
       y(i) = reduce(ty + ta * tx)
     END DO
  
  END SUBROUTINE 
  !}}}
  !{{{RECURSIVE SUBROUTINE mod_sdot(res,x,y)
  RECURSIVE SUBROUTINE mod_sdot(res,x,y)
  INTEGER(KIND=Sint), INTENT(OUT) :: res
  INTEGER(KIND=Sint), INTENT(IN) :: x(Nstate) , y(Nstate)
  INTEGER(KIND=Dint) dx, dy, dtmp
  INTEGER(KIND=Sint) tmp
  INTEGER(KIND=Ctype) i
  
    tmp = 0
    DO i=1,Nstate
     dx = x(i)
     dy = y(i)
     dtmp = tmp
     tmp = reduce(dtmp + dx * dy)
    END DO
    res = tmp
  END SUBROUTINE
  !}}}
  !{{{RECURSIVE SUBROUTINE p_saxpy(y,a)
  RECURSIVE SUBROUTINE p_saxpy(y,a)
   ! Calculates mod_saxpy(y,a,P)
   INTEGER(KIND=Sint), INTENT(INOUT) :: y(Nstate)
   INTEGER(KIND=Sint), INTENT(IN) :: a
   INTEGER(KIND=Dint) tmp, dy, da
     dy = y(1)
     da = a
     tmp = dy + da*A5
     y(1) = reduce(tmp)
     dy = y(5)
     da = a
     tmp = dy + da*A1
     y(5) = reduce(tmp)
  
  END SUBROUTINE
  !}}}
  !{{{RECURSIVE SUBROUTINE p_sdot(res,n,x)
  RECURSIVE SUBROUTINE p_sdot(res,x)
  INTEGER(KIND=Sint), INTENT(OUT) :: res
  INTEGER(KIND=Sint), INTENT(IN) :: x(Nstate)
  INTEGER(KIND=Dint) dx1, dx5, dtmp
    dx1 = x(1)
    dx5 = x(5)
    
    dtmp = A1*dx5 + A5*dx1
    res = reduce(dtmp)
  END SUBROUTINE
  !}}}
  !{{{RECURSIVE SUBROUTINE poly_mult(a,b)
  RECURSIVE SUBROUTINE poly_mult(a,b)
    INTEGER(KIND=Sint), INTENT(INOUT) :: a(Nstate)
    INTEGER(KIND=Sint), INTENT(IN) :: b(Nstate)
    INTEGER(KIND=Sint) tmp((2*Nstate) - 1)
    INTEGER(KIND=Ctype) i
  
    tmp = 0_Sint
  
    DO i=1,Nstate
      CALL mod_saxpy(tmp(i:Nstate+i-1),a(i), b)
    END DO
    DO i=(2*Nstate)-1, Nstate+1, -1
      CALL P_SAXPY(tmp(i-Nstate:i-1),tmp(i))
    END DO
    a = tmp(1:Nstate)
  END SUBROUTINE
  !}}}
  !{{{RECURSIVE SUBROUTINE poly_square(a)
  RECURSIVE SUBROUTINE poly_square(a)
    INTEGER(KIND=Sint), INTENT(INOUT) :: a(Nstate)
    INTEGER(KIND=Sint) tmp((2*Nstate) - 1)
    INTEGER(KIND=Ctype) i
  
    tmp = 0_Sint
  
    DO i=1,Nstate
      CALL mod_saxpy(tmp(i:Nstate+i-1),a(i), a)
    END DO
    DO i=(2*Nstate)-1, Nstate+1, -1
      CALL P_SAXPY(tmp(i-Nstate:i-1),tmp(i))
    END DO
    a = tmp(1:Nstate)
  END SUBROUTINE
  !}}}
  !{{{RECURSIVE SUBROUTINE poly_power(poly,n)
  RECURSIVE SUBROUTINE poly_power(poly,n)
   INTEGER(KIND=Sint), INTENT(INOUT) :: poly(Nstate)
   INTEGER, INTENT(IN) :: n
   INTEGER nn
   INTEGER(KIND=Sint) x(Nstate), out(Nstate)
  
   IF( n .EQ. 0 )THEN
     poly = 0_Sint
     poly(1) = 1_Sint
     RETURN
   ELSE IF( n .LT. 0 )THEN
     poly = 0_Sint
     RETURN
   END IF
  
   out = 0_sint
   out(1) = 1_Sint
   x = poly
   nn = n
   DO WHILE( nn .GT. 0 )
     IF( MOD(nn,2) .EQ. 1 )THEN
       call poly_mult(out,x)
     END IF
     nn = nn/2
     IF( nn .GT. 0 )THEN
       call poly_square(x)
     END IF
   END DO 
   poly = out
  
  END SUBROUTINE poly_power
  !}}}
  !}}}

  !{{{RECURSIVE SUBROUTINE  Rand_seed( state, n )
  RECURSIVE SUBROUTINE  Rand_seed( state, n )
    TYPE(Rand_state), INTENT(OUT) :: state
    INTEGER, INTENT(IN) :: n
    ! initialise the genrator using a single integer
    ! fist initialise to an arbitrary state then boost by a multiple 
    ! of a long distance
    !
    ! state is moved forward by P^n steps
    ! we want this to be ok for seperating parallel sequences on MPP machines
    ! P is taken as a prime number as this should prevent strong correlations
    ! when the generators are operated in tight lockstep.
    ! equivalent points on different processors will also be related by a
    ! primative polynomial
    ! P is 2^48-59
    TYPE(Rand_state) tmp
    TYPE(Rand_offset), PARAMETER ::  P = &
         Rand_offset( (/ 1509238949_Sint ,2146167999_Sint ,1539340803_Sint , &
                     1041407428_Sint ,666274987_Sint /) )
  
    CALL Rand_load( tmp, (/ 5, 4, 3, 2, 1 /) )
    state = Rand_boost( tmp, Rand_mul_offset(P, n ))
  
  END SUBROUTINE Rand_seed
  !}}}
  !{{{RECURSIVE SUBROUTINE Rand_load( state, input )
  RECURSIVE SUBROUTINE Rand_load( state, input )
  TYPE(RAND_state), INTENT(OUT) :: state
  INTEGER, INTENT(IN) :: input(Nstate)
  
  INTEGER(KIND=Ctype) i
  
    state%state = 0_Sint
    DO i=1,Nstate
      state%state(i) = MOD(INT(input(i),KIND=Sint),M)
    END DO
    state%need_fill = .TRUE.
    state%pos = Nstate + 1
  END SUBROUTINE Rand_load
  !}}}
  !{{{RECURSIVE SUBROUTINE Rand_save( save_vec,state )
  RECURSIVE SUBROUTINE Rand_save( save_vec, x ) 
  INTEGER, INTENT(OUT) ::  save_vec(Nstate)
  TYPE(RAND_state), INTENT(IN) ::  x
  
  INTEGER(KIND=Ctype) i
    DO i=1,Nstate
      save_vec(i) = x%state(x%pos-(Nstate+1) + i)
    END DO
  END SUBROUTINE Rand_save
  !}}}

  !{{{RECURSIVE SUBROUTINE Rand_set_offset( offset, n )
  RECURSIVE SUBROUTINE Rand_set_offset( offset, n )
  TYPE(Rand_offset), INTENT(OUT) :: offset
  INTEGER, INTENT(IN) :: n
  
    offset%poly = 0_Sint
    IF ( n .GE. 0 ) THEN
      offset%poly(2) = 1_Sint
      call poly_power(offset%poly,n)
    ELSE
      !
      ! This is X^-1 
      !
      offset%poly(4) = 858869107_Sint
      offset%poly(5) = 1840344978_Sint    
      call poly_power(offset%poly,-n)
    END IF
  END SUBROUTINE Rand_set_offset
  !}}}
  !{{{TYPE(Rand_offset) RECURSIVE FUNCTION Rand_add_offset( a, b )
  TYPE(Rand_offset) RECURSIVE FUNCTION Rand_add_offset( a, b )
  TYPE(Rand_offset), INTENT(IN) :: a, b
  
    Rand_add_offset = a
    CALL poly_mult(Rand_add_offset%poly,b%poly)
    RETURN
  END FUNCTION Rand_add_offset
  !}}}
  !{{{TYPE(Rand_offset) RECURSIVE  FUNCTION Rand_mul_offset( a, n )
  TYPE(Rand_offset) RECURSIVE  FUNCTION Rand_mul_offset( a, n )
  TYPE(Rand_offset), INTENT(IN) :: a
  INTEGER, INTENT(IN) :: n
    Rand_mul_offset = a
    CALL poly_power(Rand_mul_offset%poly,n)
    RETURN
  END FUNCTION Rand_mul_offset
  !}}}
  !{{{RECURSIVE FUNCTION Rand_boost(x, offset)
  RECURSIVE FUNCTION Rand_boost(x, offset)
  TYPE(Rand_state) Rand_boost
  TYPE(Rand_state), INTENT(IN) ::  x
  TYPE(Rand_offset), INTENT(IN) :: offset
  INTEGER(KIND=Sint) tmp(2*Nstate-1), res(Nstate)
  INTEGER(KIND=Ctype) i
  
    DO i=1,Nstate
      tmp(i) = x%state(x%pos-(Nstate+1) + i)
    END DO
    tmp(Nstate+1:) = 0_Sint
  
    DO i=1,Nstate-1
      call P_SDOT(tmp(i+Nstate),tmp(i:Nstate+i-1))
    END DO
  
    DO i=1,Nstate
      call mod_sdot(res(i),offset%poly,tmp(i:Nstate+i-1))
    END DO
    Rand_boost%state = 0_Sint
    DO i=1,Nstate
      Rand_boost%state(i) = res(i)
    END DO
    Rand_boost%need_fill = .TRUE.
    Rand_boost%pos = Nstate + 1
  
  END FUNCTION Rand_boost
  !}}}
  !{{{RECURSIVE FUNCTION Rand_step(x, n)
  RECURSIVE FUNCTION Rand_step(x, n)
  TYPE(Rand_state) Rand_step
  TYPE(RAND_state), INTENT(IN) ::  x
  INTEGER, INTENT(IN) :: n
  TYPE(Rand_offset) tmp
  
    CALL Rand_set_offset(tmp,n)
    Rand_step=Rand_boost(x,tmp)
  
  END FUNCTION
  !}}}
  
  !{{{RECURSIVE FUNCTION Rand_sint(x)
  RECURSIVE FUNCTION Rand_sint(x)
    TYPE(RAND_state), INTENT(INOUT) :: x
    INTEGER(KIND=Sint)  Rand_sint
    IF( x%pos .GT. Nstore )THEN
      CALL repack_state(x)
    END IF
    IF( x%need_fill ) CALL fill_state(x)
    Rand_sint = x%state(x%pos)
    x%pos = x%pos + 1
    RETURN
  END FUNCTION Rand_sint
  !}}}
  !{{{RECURSIVE SUBROUTINE Rand_sint_vec(iv,x)
  RECURSIVE SUBROUTINE Rand_sint_vec(iv,x)
    INTEGER(KIND=Sint), INTENT(OUT)  :: iv(:)
    TYPE(RAND_state), INTENT(INOUT)  ::  x
    INTEGER left,start, chunk, i
  
    start=1
    left=SIZE(iv)
    DO WHILE( left .GT. 0 )
      IF( x%pos .GT. Nstore )THEN
        CALL repack_state(x)
      END IF
      IF( x%need_fill ) CALL fill_state(x)
  
      chunk = MIN(left,Nstore-x%pos+1)
      DO i=0,chunk-1
        iv(start+i) = x%state(x%pos+i)
      END DO
      start = start + chunk
      x%pos = x%pos + chunk
      left = left - chunk
    END DO
  
    RETURN
  END SUBROUTINE Rand_sint_vec
  !}}}


END MODULE Rand_int

!}}}

!{{{Rand (use Rand_int to make random reals)

MODULE Rand
  USE Rand_int
  IMPLICIT NONE

!{{{Parameters

  INTEGER, PARAMETER :: RAND_kind1 = SELECTED_REAL_KIND(10)
  INTEGER, PARAMETER :: RAND_kind2 = SELECTED_REAL_KIND(6)

  INTEGER, PARAMETER, PRIVATE :: Max_block=100
  INTEGER(KIND=Sint), PRIVATE, PARAMETER  :: M = 2147483647
  REAL(KIND=RAND_kind1), PRIVATE, PARAMETER :: INVMP1_1 = ( 1.0_RAND_kind1 / 2147483647.0_RAND_kind1 )
  REAL(KIND=RAND_kind2), PRIVATE, PARAMETER :: INVMP1_2 = ( 1.0_RAND_kind2 / 2147483647.0_RAND_kind2 )

  LOGICAL, PARAMETER :: Can_step = Can_step_int
  LOGICAL, PARAMETER :: Can_reverse = Can_reverse_int

!}}}
  PUBLIC Rand_real


INTERFACE Rand_real
  MODULE PROCEDURE Rand_real1
  MODULE PROCEDURE Rand_real2
  MODULE PROCEDURE Rand_real_vec1
  MODULE PROCEDURE Rand_real_vec2
END INTERFACE


CONTAINS

  !{{{RECURSIVE SUBROUTINE Rand_real1(y,x)
  RECURSIVE SUBROUTINE Rand_real1(y,x)
    REAL(KIND=RAND_kind1), INTENT(OUT) :: y
    TYPE(RAND_state), INTENT(INOUT) ::  x
    INTEGER(KIND=Sint) Z
  
    Z = Rand_sint(x)
    IF (Z .EQ. 0) Z = M
  
    y = ((Z-0.5d0)*INVMP1_1)
    RETURN
  END SUBROUTINE Rand_real1
  !}}}
  !{{{RECURSIVE SUBROUTINE Rand_real2(y,x)
  RECURSIVE SUBROUTINE Rand_real2(y,x)
    REAL(KIND=RAND_kind2), INTENT(OUT) :: y
    TYPE(RAND_state), INTENT(INOUT) ::  x
    INTEGER(KIND=Sint) Z
  
    Z = Rand_sint(x)
    IF (Z .EQ. 0) Z = M
  
    y = ((Z-0.5d0)*INVMP1_1)  ! generate in double and truncate.
    RETURN
  END SUBROUTINE Rand_real2
  !}}}

  !{{{RECURSIVE SUBROUTINE Rand_real_vec1(rv,x)
  RECURSIVE SUBROUTINE Rand_real_vec1(rv,x)
    TYPE(RAND_state), INTENT(INOUT) ::  x
    REAL(KIND=RAND_kind1)  rv(:)
    INTEGER left,start, chunk, i
    INTEGER(KIND=Sint) Z
    INTEGER(KIND=Sint) temp(MIN(SIZE(rv),Max_block))
  
    start=0
    left=SIZE(rv)
    DO WHILE( left .GT. 0 )
      chunk = MIN(left,Max_block)
      CALL Rand_sint_vec(temp(1:chunk),x)
      DO i=1,chunk
       Z = temp(i)
       IF (Z .EQ. 0) Z = M
       rv(start+i) = (Z-0.5d0)*INVMP1_1
      END DO 
      start = start + chunk
      left = left - chunk
    END DO
  
    RETURN
  END SUBROUTINE Rand_real_vec1
  !}}}
  !{{{RECURSIVE SUBROUTINE Rand_real_vec2(rv,x)
  RECURSIVE SUBROUTINE Rand_real_vec2(rv,x)
    TYPE(RAND_state), INTENT(INOUT) ::  x
    REAL(KIND=RAND_kind2)  rv(:)
    INTEGER left,start, chunk, i
    INTEGER(KIND=Sint) Z
    INTEGER(KIND=Sint) temp(MIN(SIZE(rv),Max_block))
  
    start=0
    left=SIZE(rv)
    DO WHILE( left .GT. 0 )
      chunk = MIN(left,Max_block)
      CALL Rand_sint_vec(temp(1:chunk),x)
      DO i=1,chunk
       Z = temp(i)
       IF (Z .EQ. 0) Z = M
       rv(start+i) = (Z-0.5d0)*INVMP1_2
      END DO 
      start = start + chunk
      left = left - chunk
    END DO
  
    RETURN
  END SUBROUTINE Rand_real_vec2
  !}}}
END MODULE Rand

!}}}

!{{{test program
! PROGRAM test_random
! use Rand
!     TYPE(RAND_state) x
!     REAL y
!      CALL Rand_load(x,(/5,4,3,2,1/)) 
!      DO I=0,10
!       CALL Rand_real(y,x)
!       WRITE(*,10) I,y
!      END DO
!
!10    FORMAT(I10,E25.16)
!
!     END

!         0   0.5024326127022505E-01
!         1   0.8260946767404675E-01
!         2   0.2123264316469431E-01
!         3   0.6926658791489899E+00
!         4   0.2076155943796039E+00
!         5   0.4327449947595596E-01
!         6   0.2204052871093154E-01
!         7   0.1288446951657534E+00
!         8   0.4859915426932275E+00
!         9   0.5721384193748236E-01
!        10   0.7996825082227588E+00
!


!}}}

