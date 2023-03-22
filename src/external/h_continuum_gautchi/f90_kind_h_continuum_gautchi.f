!  f90_kind_h_continuum_gautchi.f90
!
!  Module defining the KIND numbers for the NAGWare f90 Compiler.
!
      MODULE f90_kind_h_continuum_gautchi
!
      INTRINSIC kind,selected_int_kind,selected_real_kind  ! we use these,
!
! Indicator that the KIND= is not available for this compiler/host
      INTEGER, PARAMETER :: not_available = -1  
!
! Real and Complex numbers
!   Single precision
      INTEGER, PARAMETER :: single  = kind(0.0)
!   Double precision
      INTEGER, PARAMETER :: double  = kind(0.0d0)
!   Quadruple precision
      INTEGER, PARAMETER :: quad    = selected_real_kind(p=30)
!
! INTEGERs numbers:

!   Single byte INTEGER
      INTEGER, PARAMETER :: int8    = selected_int_kind(2)
!   Two byte INTEGER   
      INTEGER, PARAMETER :: int16   = selected_int_kind(4)
!   Four byte INTEGER
      INTEGER, PARAMETER :: int32   = selected_int_kind(9)
!   Eight byte INTEGER
      INTEGER, PARAMETER :: int64   = selected_int_kind(18)
!
! Logical values
!   Single byte logical
      INTEGER, PARAMETER :: byte    = 1
!   Two byte logical
      INTEGER, PARAMETER :: twobyte = 2
!   Four byte logical
      INTEGER, PARAMETER :: word    = kind(.TRUE.)
!
! Character type
!   Normal single byte character (ASCII sequence)
      INTEGER, PARAMETER :: ascii   = kind('x')
      
      END module f90_kind_h_continuum_gautchi
