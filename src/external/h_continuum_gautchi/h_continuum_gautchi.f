
      function fun(x)

        USE f90_kind_h_continuum_gautchi

 		
        implicit none 
	
        REAL(double) :: fun, x 
        
        REAL(double) ::  amybshec
      
        common/aPARAMETER/ amybshec  
      
	
        fun = x * DLOG(x) - amybshec 
        
      END  function fun

      MODULE h_continuum_gautchi

      CONTAINS

!=============================================================
!
! SUBROUTINE Hydrogenic_Continuum_Wave(kmom, rpt, ZE, Lmax,
!     &                                      Radial, DeltaL)! 
! 
! computes the radial wave function  Radial and the Coulomb
! phase sigma for given momentum kmom, grid point rpt,
! nuclear charge ZE, and the largest angurlar momentum Lmax,
! for all the L's: 0, 1, 2, 3, ...., Lmax 
! at a single call
!
! Authors: Liang-You Peng, Peking University,
! Dates: June, 2010 
!
!=============================================================            

      SUBROUTINE Hydrogenic_Continuum_Wave(kmom, rpt, ZE, Lmax,
     &                                      Radial)
      USE f90_kind_h_continuum_gautchi
      
      IMPLICIT NONE 

! The momentum kmom and grid point rpt, where radial wave function needed:
      
      REAL(double)  :: kmom, rpt
      
! The nuclear charge of the Hydrogenic atom:
      
      REAL(double)  :: ZE

! The largest angular momentum needed:
      
      INTEGER :: Lmax
  
! The returned radial wave function at rpt for all L's
! and the Coulomb phase 
     
      REAL(double) ::  Radial(0:Lmax)
      
!--------------------------------------------------------      
! Local variables:

      REAL(double) :: rho, eta 
      
      
      rho = rpt * kmom
      
      eta = - ZE / kmom 
      
! First calculate the radial wave function for all L's
! from L = 0, ..., Lmax :
      
       CALL Regular_Coulomb_Fun(eta, rho, Lmax, Radial)

! Due to the different convention, there is a factor of 2,
! which should be multipled in order to be consistent with 
! Landau's momentum normalized continuum states:

       Radial  =  2.D0 * Radial   
 
      END  SUBROUTINE Hydrogenic_Continuum_Wave   
      
      

      
      
!===============================================================
! This code is very accurate, the results can be 
! checked against those from Mathematica by using 
!
! 2^L Exp[ - Pi yita / 2.0] Abs[Gamma[L + 1.0 + I yita ]] / 
!  Gamma[ 2 L + 2.0] rho^(L + 1) Exp[-I rho]
!  Hypergeometric1F1[L + 1 - I yita, 2.0 * L + 2.0, 2 I rho]
!
!===============================================================
!   SUBROUTINE Regular_Coulomb_Fun(eta, ro, Lmax,  F) 
! returns an array F[0:Lmax] for the values of 
! regular Coulomb wave function at given arguments
! eta and ro. Accuracy: around 14 digits.
! [ALGOL  is a short for ALGOrithmic Language.]
!
! Author: Liang-You Peng, May 10, 2010
!
! This code is essentially a Fortran coding of the  
! original algorithm (in ALGOL 60) by  Walter Gautschi 
! " Algorithm 292: Regular Coulomb Wave Funcitons."
! References:
! (1) Communications of the ACM, V9, (1966) p793;
! (2) Communications of the ACM, V12, (1969) p280;
! (3) Communications of the ACM, V13, (1970) p573;
! (4) SIAM Review, V9, (1967) p24,  "Computational  
! Aspects of Three-Term Recurrence Relations", See Sect7   
! (5) Carl-Erik Fr\"oberg, Rev. Mod. Phys, V27 (1955), p399-411;
! "Numerical Treatment of Coulomb Wave Functions".
!================================================================
     
      SUBROUTINE Regular_Coulomb_Fun(eta, ro, Lmax,  F)

      USE f90_kind_h_continuum_gautchi
      
      IMPLICIT NONE 
                      
      INTEGER :: Lmax
     
      REAL(double)  :: eta, ro, t 
     
      REAL(double)  :: F(0:Lmax)
      
      REAL(double), PARAMETER  :: Ebase = 2.718281828459045235D0
      
      REAL(double), PARAMETER  :: EInv = 0.36787944117144232159D0

      REAL(double), PARAMETER  :: PI    = 3.141592653589793238D0
      
      INTEGER :: L, nu, nu1, mu, mu1, i,   d 
      
      REAL(double) :: scaling
      
      REAL(double)  :: epsilon, ro1, eta2, omega, d1, sum,  
     &           r, r1, s, t1, t2
      
      REAL(double) :: lmin(0:800000),Rr(0:Lmax),Rra(0:800000) 
    
      REAL(double) :: lam(0:800000),lambda(0:800000),Fapprox(0:800000)
      
      
      
            
      IF (ro<0) STOP 'Ro is negative!'
      
      IF (ro<1.D-16) THEN
         
         F = 0.D0
	 
	 RETURN
	 	 
      ENDIF 

! Scale for prevention of overflow:
      
      scaling = 1.D0**(-300)

! The precision for the calculation:
      
      d = 14
      
      epsilon =  1.0D-14
      
       

! The precision for the calculation.
      
      
      ro1 = 1.D0 / ro
      
      eta2 = eta**2
                        
         
      IF (eta>0) THEN
      
         t1 = 0.5D0 * ro / eta
      
      ELSE 
       
         t1 = 0.D0 
	 
      ENDIF 
              
	            
      IF (eta<1.D0) THEN
      
         omega = 0.0D0 
      
      ELSE IF(t1>=1.D0) THEN 
       
         omega = 0.5D0 * PI / t1
	 
      ELSE
      
         omega = ( 0.5D0*PI -  DATAN(DSQRT(1.D0/t1-1.D0)) +
     &           DSQRT(t1*(1.D0-t1)) ) / t1
      
      ENDIF


      
      lambda(0) = scaling  
  
      lmin(0)   = 1.D0 
       
      lambda(1) = scaling* ( omega - eta ) 
      
      sum       = ro * DEXP(omega*ro) * scaling


! Initialization:
    
      Fapprox = 0.D0 
      
!======================================
! log_e 10 =  2.30258509299404568D0
! 2 log_e 2 = 1.3862943611198906188D0
! e /2 =  2.7182818284590452353602875D0 /2  
!----------------------------                          
      
      d1 = 2.30258509299404568D0 * DBLE(d) +1.3862943611198906188D0
      
      t1 = 1.35914091422952261768D0 * ro

      IF ( Lmax<NINT(t1) ) THEN
      
         L = 1 + NINT(t1)
	 
      ELSE
      
         L = Lmax
	 
      ENDIF
      
      
      t1 = DEXP(PI*0.5D0*eta)
      
      s = DSQRT(1.D0+omega**2)
      
     
      IF (ABS(omega)<1.D-16) THEN
      
         t1 = t1 + 1/ t1
	 
      ELSE
      
         t1 = DEXP( -eta * DATAN(1.D0/omega) )
	 
      ENDIF 

      
      t2 = omega  + s 
      
      r  = 1.35914091422952261768D0 * ro * t2
      
      s  = ( d1 + DLOG(t1*DSQRT(t2/s)) -omega * ro ) / r

! The subroutine  inver_xlogx(y, x) calculates the value  x of 
! inverse  function of x Log x, for given value of y.
! The method is by numerical searching the root for x Log (x) - y = 0.      
    
      IF (s >=-EInv) THEN
      
        CALL inver_xlogx(s,  t) 
      
        nu = NINT( r * t )
      
      ELSE
      
        nu = NINT(r/Ebase) 
	
      ENDIF 

      
      CALL inver_xlogx(0.5D0*d1/DBLE(L),  t) 
                
     
      nu1 = NINT( DBLE(L) * t )
                
      IF (nu<nu1) nu = nu1
                            
      
      nu1 = 1 
      
      IF (ABS(omega)<1.D-14) THEN
      
         i = 1 ! L2
	 
      ELSE
      
         i = 2 ! L1;  i =3 ==> M1
	 
      ENDIF 

      
10    IF (i==1) THEN  ! L0, 10  
      
         GOTO 200 ! i = 1, L2
	 
      ELSE IF(i==2) THEN
      
         GOTO 100 ! i = 2, L1
	 
      ELSE  IF (i==3) THEN
      
         GOTO 300 ! i = 3, M1
	 
      ENDIF 
               
                                            
100   d1 = 2.D0 * eta / ( DEXP(2.D0*eta*DATAN(1.D0/omega)) - 1.D0 )
      

      
      i = 3 
      
300   mu1 = 0 ! 300, M1, i = 3
                 
      mu = NINT(1.25D0*nu)
      
      lam(mu1:nu) = 0.D0 
      

220   r = 0.D0  ! 220, M2
      
                     
      DO L = mu, mu1+ 1, -1 
      
         r = - ( DBLE(L**2)+ eta2 ) /  
     &     ( DBLE(L) * (  DBLE(2*L+1) * omega - DBLE(L+1) * r ) )
     
         IF (L<=nu) Rra(L-1) = r 
	 
      ENDDO
      
             
      
      DO L = mu1 + 1, nu
               
          lmin(L) = Rra(L-1) * lmin(L-1) 
	
	  IF (ABS(t1)<1.D-15) lmin(L) = 0.D0  	 	 	 	  
      
      ENDDO 
             
              
      
      DO L = mu1, nu
       
         IF ( ABS(lmin(L)-lam(L)) > epsilon * ABS(lmin(L)) ) THEN 
      
            lam(mu1:nu) = lmin(mu1:nu)
	    
	    mu = mu + 5 
	    
	    IF (mu<5*nu) THEN
	    
	       GOTO 220 ! M2
	       
	    ELSE
	    
	       STOP ' Convergence difficulty in Lambda[L]'
	       
	    ENDIF 
	    
	 
	 ENDIF 
      
      ENDDO 
              
	 
      r1     = 1.D-10	      
      
      lam(0) = - r1 * scaling 
      
      lam(1) = scaling  
      
      t1     = d1 / ( 1.D0 + r1**2 )
      
      
      DO L = 2, nu
       
        lam(L)   = (2.D0 - 1.d0/DBLE(L)) * omega * lam(L-1) 
     & + ( 1.d0-1.D0/DBLE(L) + eta2/DBLE(L-1)/DBLE(L) ) * lam(L-2)  
     
     
        lambda(L)= lmin(L) * scaling + t1 * (lam(L) 
     &            + scaling * r1 * lmin(L))        
     
      ENDDO 
      
     
      GOTO 330 ! L3
                       
    
200   DO L = nu1, nu-1
            
         lambda(L+1) = (2.D0+1.D0/DBLE(L)) / (1.D0+1.D0/DBLE(L))
     &   * omega * lambda(L) + ( 1.d0/(1.d0 + 1.d0/DBLE(L)) + 	 
     &   eta2/ DBLE(L)  / DBLE(L+1) ) * lambda(L-1)

      ENDDO
      
  
330   r = 0.D0 
      
      s = 0.D0 
      
       DO L = nu, 1, -1
       
         t1 = eta / DBLE(L+1)
	 
	 r  = 1.D0 / ( DBLE(2*L-1) * ( t1 / DBLE(L) + ro1 
     &       - (1.D0 + t1**2) * r / DBLE(2*L+3) ) ) 
     
         s  = r * ( Lambda(L) + s )  
	 
	 IF (L<=Lmax) Rr(L-1) = r 
	 
      ENDDO
      
 
      F(0) = sum / (scaling + s)
      
      DO L = 1, Lmax 
      
         F(L) = Rr(L-1) * F(L-1) 
	 
      ENDDO 
      
      
      DO L = 0, Lmax 
      
         IF ( ABS(F(L)-Fapprox(L)) > epsilon * ABS(F(L)) ) THEN 
	 
	    Fapprox(0:Lmax) = F(0:Lmax)
	  	  
	    mu1 = nu
	  
	    nu1 = mu1 
	  
	    nu = nu + 10
	  	  	  
	    IF (nu<800000) THEN
	  
	       GOTO 10 !L0
	      
	    ELSE
	   
	       STOP 'convergence difficulty in Coulomb!' 
	      
	    ENDIF
	  
	 ENDIF 
	   
      ENDDO 
      
    
      t1 = 2.D0*PI * eta
      
      IF (ABS(t1)<1) THEN 
       
         s  = 1.d0
	 
	 t2 = s
	 
	 L  = 1

400	 L = L + 1

         t2 = t1 * t2 / DBLE(L)
	 
	 s = s + t2
	 
	 IF ( ABS(t2)>epsilon*ABS(s) ) GOTO 400 ! L4
	 
	 s = DSQRT(1.D0/s)
      
      
      ELSE
      
         s = DEXP(-t1/4.D0) / DSQRT( ( DEXP(t1/2.D0)
     &      - DEXP(-t1/2.D0) ) / t1 )
	       
      ENDIF 
           
     
      F(0) = s * F(0) 
      
      
      DO L = 1, Lmax 
            
       s = (1.d0 - 1.d0/DBLE(2*L)) / (1.d0 + 1.d0/DBLE(2*L) ) 
     &         * DSQRT( 1.d0 + (eta/DBLE(L))**2 ) * s     
 
       F(L) = s * F(L) 
       
      ENDDO 
      
                                  
      END SUBROUTINE Regular_Coulomb_Fun ! L5








! The inverse of function y = t * log (t), for given y:

! Note, when x<1, there are two solutions, we take
! the larger one in this program

	
      SUBROUTINE inver_xlogx(y,  t) 
      
      USE f90_kind_h_continuum_gautchi
       
      IMPLICIT NONE 
      
      external fun 
      
      REAL(double) :: fun 
      
      REAL(double) :: t, y 
      
      REAL(double)::   EInv
      
      REAL(double) :: x0, x2, x1, eps 
      
      INTEGER :: imax
      
      REAL(double) ::  amybshec
      
      common/aPARAMETER/ amybshec  
      
      

    
      EInv = 0.36787944117144232159D0

! Set precision needed for the root finder:
     
      eps = 1.D-15
     
     
      amybshec = y 
      
 
      IF (amybshec < -EInv ) STOP 'y<-1/e in FUNCTION t(y)'
      

! Only one zero for the function for x*ln (x), at x = 1.d0
      
      IF ( ABS(amybshec) < 1.D-14) THEN
      
          t = 1.D0 	  	
	  
      ENDIF 
      
! If amybshec> 0, then the root is unique (x>1.0):
      
      IF ( amybshec > 0.D0 ) THEN
      
          x0 = 1.D0 + 1.D-10
	  
	  x2 = 1.D20

! Max iteration permitted:
	  
	  imax = 800000
	  	              
          CALL iterate(x0, x2, fun, eps, imax, x1)
	  
	  t = x1 
	  
      ENDIF 
      
! IF  amybshec < 0.D0  .AND.  amybshec >= -1/e, there are two solutions:
       
    
       IF ( amybshec < 0.D0  .AND.  amybshec >= -EInv ) THEN

! The smaller soluttion:
      
!!          x0 =   1.D-20
	  
!!	  x2 =   EInv
	  
!!	  imax = 8000000
           
!!         CALL iterate(x0, x2, fun, eps, imax, x1)
!!	  
!!	  t = x1 
	  	  	  
! The larger solution
      
          x0 =   EInv
	  
	  x2 =   1.D0 - 1.D-10
	  
	  imax = 800000
           
          CALL iterate(x0, x2, fun, eps, imax, x1)
	     
	  t = x1
	  
       ENDIF 
    
      END SUBROUTINE inver_xlogx 
	
      



! The world's best root finder alogrithm.
! http://www.embedded.com/story/OEG20030508S0030
! By Jack W. Crenshaw 
! Jack Crenshaw is a senior software engineer at Spectrum-Astro and 
! the author of Math Toolkit for Real-Time Programming, from CMP 
! Books. He holds a PhD in physics from Auburn University. 
! E-mail him at jcrens@earthlink.net. 

! Convert from C to Fortran by Liang-You Peng

! Make sure that there is ONLY one root between x0 and x2 
! Therefore, 
      
        subroutine iterate(x0, x2, fun, eps, imax, x1)

        USE f90_kind_h_continuum_gautchi

        implicit none 
  
        external fun
  
        REAL(double) :: fun 

        REAL(double) :: x1, y0, y1, y2, b, c, temp, y10, y20, y21,
     &            xm, ym, eps, x0,x2
 
        REAL(double) :: xmlast 
  
        INTEGER:: imax, i
  
        
	xmlast = x0 
  
        y0 = fun(x0)

! Is x0 a root?
  
        if ( abs(y0)<1.D-15 ) then
         
	   x1 = x0
   
           return
   
        ENDif 
  
        y2 = fun(x2)
	
! Is x2 a root?   
 
        if ( abs(y0)<1.D-15 ) then
  
           x1 = x2
  
           return 
  
        ENDif 

! Are y2 and y0 have the same sign?

        if (y2 * y0 > 0.0) THEN
  
           ! stop 'No zero OR at least two zeros between x0 and x2.'
           return
  
        ENDif

! Now do the iteration:


        do i = 0, imax 
  
           x1 = 0.5d0 * (x2 + x0)
           
	   y1 = fun(x1)
   
           if ( abs(y1)<1.D-16 .OR. abs(x1 - x0) < eps) then
	   
              return 
          
	   ENDif 
                 
           if (y1 * y0 > 0.d0) then 
    
               temp = x0
               x0   = x2
               x2   = temp
               temp = y0
               y0   = y2
               y2   = temp
      
           ENDif 
    
           y10 = y1 - y0
           y21 = y2 - y1
           y20 = y2 - y0
    
           
	   if (y2 * y20 <2.d0 * y1 * y10) then 
    
              x2 = x1
              y2 = y1
    
           else 
    
              b = (x1 - x0) / y10   
              c = (y10 -y21) / (y21 * y20) 
              xm = x0 - b * y0 * (1.d0-c * y1)
              ym = fun(xm)
      
              if (abs(ym)<1.D-16 .OR. abs (xm - xmlast) < eps) then
                 x1 = xm
             	 return 
              ENDif 
                 
              xmlast = xm
      
              if (ym * y0 < 0.d0) then 
      
                 x2 = xm
                 y2 = ym
       
              else
      
                 x0 = xm
                 y0 = ym
                 x2 = x1
                 y2 = y1
     
              ENDif  ! if (ym * y0 < 0.d0)
   
           ENDif !  if (y2 * y20 <2.d0 * y1 * y10)
    
         ENDdo ! do i = 0, imax 
            
       END subroutine iterate
         

!============================================
!
! Peng, 2010.05.20.
!
! Another way to calculate the Coulomb phase:
! 
! Given momentum xk, L, and   charge ZE,
! the code returns the Coulomb phase for all the
! L's: 0, 1, 2, ..., L,
! in the matrix DeltaL(0:L)
!
!=============================================

      SUBROUTINE sigma_arb(xk, l, Zar, DeltaL)
      
        USE f90_kind_h_continuum_gautchi
      
      IMPLICIT REAL(double) (a-h,o-z)
      
      REAL(double) :: DeltaL(0:L), sigma0
      
      INTEGER :: LL
      
      
!Computes sigma function appearing in coulomb phase for l=0
!according to Martins, JPB 1,154 (1968):

      c=100.d0-1.d0/(xk*xk)
      IF (c) 10,10,20
 10   xl=0.d0
      GOTO 30
 20   xl=idint(dsqrt(c))
 30   q=xl+1.d0
      z=Zar/xk
      a=-0.5d0*z*dlog(q*q+z*z)+z
      xx=z/q
      b=-(xl+0.5d0)*datan(xx)
      c=z/(12.d0*(q*q+z*z))
      d=z*(z*z-3.d0*q*q)/(360.d0*(q*q+z*z)**3)
      sigma0=a+b+c+d
 35   IF (xl) 50,40,50
 40   GOTO 549
 50   xy=z/xl
      sigma0=sigma0+datan(xy)
      xl=xl-1.d0
      GOTO 35
549   CONTINUE

      LL = NINT(xl) 
       
      DeltaL(0) = sigma0

351   IF (dble(l-LL)) 501,401,501
401   RETURN
501   LL=LL+1.d0
      xy=z/dble(LL)
      sigma0=sigma0-datan(xy)
      
      DeltaL(LL) = sigma0
      
      GOTO 351

      END SUBROUTINE sigma_arb

      END MODULE
