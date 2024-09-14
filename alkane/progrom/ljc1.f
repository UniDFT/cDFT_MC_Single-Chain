	PROGRAM SLIT_HSC_LJ

!     THIS PROGRAM IS FOR CALCULATING THE SEGMENTAL DENSITIES OF ONE-COMPONENT 
!     HARD-SPHERE CHAINS CONFINED IN A HARD SLIT-LIKE PORE
!	
!	PAPER: JCP,V117,2356 (2002)

!    Modified for Lennard-Jones chains, May 27, 2004
!    Lennard-Jones potential is represented by the double Yukawa potential

!	KEY PARAMETERS

!	DZ		step length in the units of hard-sphere radius  
!	ZH		half of the pore width in the units of hard-sphere radius 
!			(because of symmetry, only half of the density profile need to be calculated)
!	NZ		NZ=1 + half pore width/DZ or pore width=2*(NZ-1)*DZ  in units of radius 
!	F		mixing parameter in the Picard iteration
!	M		number of segments per chain 
!	RHOR	reduced number density of segments, RHO*SIGMA**3/8=rho*radius**3
!	ETA		packing fraction of segments  
!	RHO(I)	reduced density profile,rho*radius**3    
!	RE		file to save temporary density profile
!	RO(K,I)	the density profile of site K  
!	RD(I)   the difference in density profiles between two consective iterations
!	SUMP	average packing fraction in the pore
!		
	IMPLICIT NONE
	INTEGER NZ,M,ITER,ICODE,I,K
	DOUBLE PRECISION DZ,PI,F,ETA,Y,RHOR,SADD,SUMP,AMU,Rc0
	PARAMETER (NZ=1601,M=1)
	DOUBLE PRECISION RHO(NZ),RO(M,NZ),RHO1(NZ),RD(NZ)
	DOUBLE PRECISION FL(M,2*NZ),FR(M,2*NZ),FF(2*NZ)
	PARAMETER (PI=3.141592654)

	DOUBLE PRECISION SIGMA,RHOB,TF,DEFF,H,Z0,ZH   
	COMMON /HARDWALL/ Z0,ZH,TF,DEFF 
!  Initial Lennard-Jones parameters

	SIGMA=1.0	! Lennard-Jones diameter  
	RHOB=0.5	! reduced density, rho*sigma**3
	TF=1.35		! reduced temperature kT/episilon
	Rc0=4.0     ! cut off distance in units of LJ diameter
	DZ=0.01  	! step length in units of LJ diameter
	H=DZ*float(NZ-1)*2	! slit width in units of LJ diameter

	ICODE=0		! program control parameter, icode=1 restart
	F=0.1		! mixing parameter

!cccccccccccccccccccccccccccccccccccccccccccccccccc

!     effective hard-sphere diameter
	DEFF=SIGMA*(1.+0.2977*TF)/(1.+0.33163*TF+0.00104771*TF**2)  
	Z0=SIGMA/DEFF

	ZH=H/DEFF*2.	! pore width in units of hard-sphere radius
	DZ=DZ/DEFF*2.	! step length in units of hard-sphere radius

! hard-sphere packing fraction
	ETA=RHOB*DEFF**3*PI/6.

! use the effective hard-sphere radius as the unit length for calculations
	RHOR=RHOB*DEFF**3/8.	! RHOR=RHO*RADIUS_HS^3

! initialize the coefficients for the DCF of two-Yukawa potential 
	CALL C2INI(TF,SIGMA,DEFF,ETA,Rc0)

! start Picard iteration
	ITER=0

! initialize the density profile
	IF (ICODE.EQ.0) THEN
! use the bulk density as the initial guess
	DO I=1,NZ
	Y=FLOAT(I-1)*DZ
	IF(Y.GT.Z0) THEN
	RHO(I)=RHOR
	ELSE
        RHO(I)=0.D0
	ENDIF
	ENDDO

	ELSE

!  read the initial density profile from a file
	OPEN(2,FILE='hsc1.ini')
	DO I=1,NZ
	READ(2,*) Y,RHO(I)
	RHO(I)=RHO(I)*RHOR
	ENDDO
 	CLOSE(2)
	ENDIF

!  calculate the reduced excess chemical potential PER SEGMENT in the bulk phase, AMU
	CALL CPS(M,ETA,AMU)

!  calculate the recurrence functions and the effective external potential 
1	CALL FFLR(M,NZ,DZ,AMU,FL,FR,FF,RHO)

!  segmental densities,RO(K,I) 
	DO 10 K=1,M
	DO 10 I=1,NZ
  	RO(K,I)=FL(K,I)*FR(K,I)*FF(I)*RHOR/float(M)
10	CONTINUE
	
!  overall density at each point
        DO 20 I=1,NZ
	RHO1(I)=0.0
	DO K=1,M
	RHO1(I)=RHO1(I)+RO(K,I)
        ENDDO

!   density difference
	RD(I)=RHO1(I)-RHO(I)
20	CONTINUE

!   update the density profile
	SADD=0.0
	DO 30 I=1,NZ
	SADD=SADD+ABS(RD(I))												    
	RHO(I)=RHO(I)*(1.-F)+RHO1(I)*F
30	CONTINUE
	ITER=ITER+1

!  store intermediate results
	WRITE(*,*) ITER,SADD

	IF (MOD(ITER,100).EQ.0) THEN
	OPEN(2,FILE='hsc1.ini')
	DO 40 I=1,NZ
	Y=DZ*FLOAT(I-1)/2.0
	WRITE(2,41)Y,RHO(I)/RHOR
41	FORMAT(4X,F8.4,4X,F10.7)
40	CONTINUE
	CLOSE(2)
	ENDIF

! convergence criteria
	DO I=1,NZ
	IF(ABS(RD(I)).GE.1.E-5) GOTO 1
	ENDDO

c  the average total packing fraction inside the pore
	SUMP=(RHO(1)+RHO(NZ))/2
	DO I=2,NZ-1
	SUMP=SUMP+RHO(I)
	ENDDO
	SUMP=SUMP*8.*PI/6./FLOAT(NZ-1)

! output results
	OPEN(2,FILE='hsc1.ini')
	OPEN(3,FILE='hsc1.out')
	WRITE(3,*)'ETA_AV=',SUMP,'  H/sigma=',H
	DO 50 I=1,NZ
	Y=DZ*FLOAT(I-1)*DEFF/2.0 ! distance from the wall in units of LJ diameter
	RHO(I)=RHO(I)/RHOR
        WRITE(2,41)Y,RHO(I)
	IF (MOD(I-1,1).EQ.0) WRITE(3,51) Y,RHO(I)*RHOB
51	FORMAT(4X,F8.4,4X,F10.7)
50	CONTINUE
	CLOSE(3)
        CLOSE(2)

	STOP
	END

	SUBROUTINE FFLR(M,NZ,DZ,AMU,FL,FR,FF,RHO)
!	calculate the recurrence functions and 
!	the self-consistent field inside the hard slit pore
!       Unlike the density profile, FL,FR,FF must be calculated 
!	within the entire slit width=2*(NZ-1)*DZ

!	all other variables are DOUBLE PRECISION, why FL, FR, FF should be double precision?
	IMPLICIT NONE
	INTEGER NZ,M,I,K,IR,J,KIN,KIP
	DOUBLE PRECISION DZ,AMU,Z1,SUM1,S,DATT
	DOUBLE PRECISION FL(M,2*NZ),FR(M,2*NZ),FF(2*NZ)
	DOUBLE PRECISION RHO(NZ),PHI

!	IR: 	   hard-sphere diameter in units DZ, IR=2/DZ  
!	SUM1:      effective external potential 
!	ZI:        position in z-direction
!	FL:		   left Green function
!	FR:		   right Gren function
!	FF:		   Boltzmann term due to the self-consistent potential

	DO 10 I=1,2*NZ-1
	FL(1,I)=1.0
	FR(M,I)=1.0
	Z1=DZ*FLOAT(I-1)

! Excess local chemcial potential due to hard-sphere repulsion and chain connectivity
	CALL SFM1D(M,DZ,NZ,Z1,SUM1,RHO)

! Excess local chemical potential due to vdw attraction
	CALL DFEXATT(Z1,DZ,NZ,RHO,DATT)
! External potential
	CALL EXTP(Z1,PHI)
	IF (PHI.GT.30.d0) THEN
	FF(I)=0.0d0
	ELSE
	FF(I)=EXP(-PHI-SUM1+DATT+AMU)
	ENDIF

10    CONTINUE

	IF(M.LT.2)RETURN

	IR=INT(2./DZ)

	DO 30 J=2,M
	DO 20 I=1,2*NZ-1

! set integration boundaries
	KIN=I-IR
	KIP=I+IR
	IF(KIN.LE.0) KIN=1
	IF(KIP.GT.(2*NZ-1)) KIP=2*NZ-1

! integrate using the trapzoidal approximation
	S=(FL(J-1,KIN)*FF(KIN)+FL(J-1,KIP)*FF(KIP))/2.
	DO K=KIN+1,KIP-1
	S=S+FL(J-1,K)*FF(K)
	ENDDO

	FL(J,I)=S*DZ/4.0
	FR(M-J+1,I)=FL(J,I)
20	CONTINUE
30	CONTINUE
	RETURN
	END


	SUBROUTINE SFM1D(M,DZ,NZ,Z,SUM1,RHO)
! THE SELF-CONSISTENT-MEAN FIELD FROM THE HARD-SPHERE CHAIN REFERENCE SYSTEM
!	DZ:     step length in numerical integration 
!       Z1:	integration variable, from Z-R to Z+R 
!               (because of the definitions of the weighted densities) 

	IMPLICIT NONE
	INTEGER NZ,M,I,NP
	DOUBLE PRECISION PI,DZ,Z1,SUM1,SF21,SF31,SF21VI,Z,DF21,DF31,DF21VI

	DOUBLE PRECISION DFN(3),DFA(3),RHO(NZ)
	PARAMETER (PI=3.141592654)

!	the number of integration points 	
	NP=1+INT(2./DZ)
!	sum zero
	SF21=0.0
	SF31=0.0
	SF21VI=0.0

	DO 10 I=1,NP
	Z1=Z-1.+DZ*FLOAT(I-1)

! calculate the derivatives of the free energy density wrt weighted densities
!    n2, n3, nv2
	CALL FMN1D(M,DZ,NZ,Z1,DFN,DFA,RHO)
!	WRITE(*,*) DZ,Z1,DFN(1),DFA(1)

! integrate using the trapezoidal approximation
	DF21=DFN(1)+DFA(1)
	DF31=(DFN(2)+DFA(2))*(1.-(Z1-Z)**2)

	DF21VI=-(DFN(3)+DFA(3))*(Z-Z1)

	IF(I.EQ.1.OR.I.EQ.NP)THEN
	DF21=DF21*0.5
	DF31=DF31*0.5
	DF21VI=DF21VI*0.5
	ENDIF

	SF21=SF21+DF21
	SF31=SF31+DF31
	SF21VI=SF21VI+DF21VI
10	CONTINUE
	SF21=2.0*PI*SF21*DZ
	SF31=PI*SF31*DZ
	SF21VI=2.*PI*SF21VI*DZ
	SUM1=SF21+SF31+SF21VI
	RETURN
	END	
	
	SUBROUTINE FMN1D(M,DZ,NZ,Z1,DFN,DFA,RHO)
!  calculate the derivatives of the excess free energy density  
!  with respect to the weighted densities at position z1 given a density profile rho(z)
	IMPLICIT NONE
	INTEGER NZ,M,I,NP,IZ2
	DOUBLE PRECISION PI,PI2,PI4,DZ,Z1,Z2,XI1,XI2,DFN3,ARR
	DOUBLE PRECISION AN2,AN3,ANV2I,F21,F31,F21VI,GS11,DGN2,DGN3,DGNV2
	PARAMETER (PI=3.141592654,PI2=2.*PI,PI4=4.*PI)
	DOUBLE PRECISION DFN(3),DFA(3),RHO(NZ)

!	DZ:     step length  
!	AN2:    weighted density,n2  
!	AN3:    weighted density,n3    
!	ANV2I:  weighted density, nv2  
!	DFN(1): derivative of the hard-sphere Helmholtz energy density with respect to n2  
!	DFN(2): derivative of the hard-sphere Helmholtz energy density with respect to n3  
!	DFN(3): derivative of the hard-sphere Helmholtz energy density with respect to nv2  
!	DFA(1): derivative of the chain Helmholtz energy density with respect to n2  
!	DFA(2): derivative of the chain Helmholtz energy density with respect to n3  
!	DFA(3): derivative of the chain Helmholtz energy density with respect to nv2  
!	GS11:   radial distribution function of hard spheres at contact

!	the number of integration points from -sigma/2 to sigma/2 (-1 to 1 in units of radius)	
	NP=1+INT(2./DZ)
	
!	zero the sums for weighted densities n2, n3 and nv2
	AN2=0.0
	AN3=0.0
	ANV2I=0.0

!	Z2: integration variable, from Z1-1 to Z1+1	
!	see Eq.(13) in JCP, V117,10156(2002). 
	DO I=1,NP
	Z2=Z1-1.+DZ*FLOAT(I-1)
	IZ2=INT(Z2/DZ)+1

! impose the boundary conditions and system symmetry 
	IF ((IZ2.GT.(2*NZ-1)).OR.(IZ2.LT.1)) THEN
	F21=0.0
	ELSEIF (IZ2.LT.NZ) THEN
   	F21=RHO(IZ2)
	ELSE
	F21=RHO(2*NZ-IZ2)
	ENDIF

!	end points
	IF(I.EQ.1.OR.I.EQ.NP)THEN
	F21=F21*0.5
	ENDIF

	F31=F21*(1.-(Z1-Z2)**2)

!	nv2 is a vector
	F21VI=F21*(Z1-Z2)

	AN2=AN2+F21
	AN3=AN3+F31
	ANV2I=ANV2I+F21VI
	ENDDO
	AN2=2.*PI*DZ*AN2
	AN3=PI*DZ*AN3
	ANV2I=2.*PI*DZ*ANV2I
	
!	assume the system is ideal when the local overall density is small?
	IF(AN3.LE.1.E-5)THEN
	DFN(1)=0.
	DFN(2)=0.
	DFN(3)=0.
	DFA(1)=0.
	DFA(2)=0.
	DFA(3)=0.
	RETURN
	ENDIF

C	WRITE(*,*) AN2,AN3,ANV2I

!	calculate the excess free energy density and its derivatives wrt 
!       weighted densities using the MFMT
!	dfn(1) = dphi/dn2, dfn(2) = dphi/dn3, dfn(3) = dphi/dnv2  	

!  xi functions
	XI1=1.-ANV2I**2/AN2**2
	XI2=1.-3.*ANV2I**2/AN2**2
	DFN3=DLOG(1.-AN3)/(18.*PI*AN3**3)+1./(36.*PI*AN3**2*(1.-AN3))
     &    +(1.-3.*AN3)/(36.*PI*(1.-AN3)**3*AN3**2)

!	dfn(1), Eq.(?) in somewhere
	DFN(1)=-DLOG(1.-AN3)/PI4
     &       +AN2/PI2/(1.-AN3)
     & +(AN2**2-ANV2I**2)/(12.*PI*AN3)*(DLOG(1.-AN3)/AN3+1./(1.-AN3)**2)

!	dfn(2), Eq.(?) in somewhere
	DFN(2)=AN2/PI4/(1.-AN3)+AN2**2*XI1/PI4/(1.-AN3)**2-AN2**3*XI2*DFN3

!	dfn(3), Eq.(?) in somewhere
	DFN(3)=-ANV2I/PI2/(1.-AN3)
     &       -AN2*ANV2I/(6.*PI*AN3)*(DLOG(1.-AN3)/AN3+1./(1.-AN3)**2)

! derivatives of the excess free energy due to the hard-sphere chain wrt weighted densities
!   GS11, the contact value of the radial distribution function of uniform hard spheres
!	see Eq.(18) in JCP, V117,10156(2002). 

	GS11=1./(1.-AN3)
     &    +0.5*AN2*XI1/(1.-AN3)**2
     &    +AN2**2*XI1/(18.*(1.-AN3)**3)

!	derivatives of GS11 wrt n2, n3 and nv2
	DGN2=0.5*(2.-XI1)/(1.-AN3)**2+AN2/(1.-AN3)**3/9.0
	DGN3=1./(1.-AN3)**2+AN2*XI1/(1.-AN3)**3+AN2**2*XI1/(1.-AN3)**4/6.0
	DGNV2=-ANV2I/AN2/(1.-AN3)**2-ANV2I/(1.-AN3)**3/9.0

	ARR=(1./FLOAT(M)-1.)/PI4

	DFA(1)=ARR*(XI1*AN2/GS11*DGN2+(2.-XI1)*DLOG(GS11))
	DFA(2)=ARR*(XI1*AN2/GS11*DGN3) 
	DFA(3)=ARR*(XI1*AN2/GS11*DGNV2-2.*ANV2I/AN2*DLOG(GS11))

	RETURN
	END


	SUBROUTINE CPS(M,ETA,AMU)

!	CPS(M, RHOR, ETA, AMU) calculates the excess chemical potential PER SEGMENT of a one-component 
!		hard-sphere-chain fluid. Each chain consists of M identical segments
!		and the total reduced number density of segments is RHOR. 
!
! KEY PARAMETERS:
!	M      the number of segments per chain
!	ETA	   packing fraction of hard sphere segments
!	RHOR   the reduced density of segments, rho*sigma**3/8 
!	AMR	   the ideal gas chemical potential	
!	AMHS   the excess chemical potential due to hard sphere	collision
!	AMC    the excess chemical potential due to chain formation  
!	AMU    the total reduced chemical potential, mu/kBT/M  
!	GHS    contact value of the radial distribution function for hard spheres 
!	DGHS   derivative of LN(GHS) with respect to density
	IMPLICIT NONE
	INTEGER M
	DOUBLE PRECISION PI,ETA,PI4,RHOR,AMHS,GHS,DGHS,AMC,AMU

	PARAMETER (PI=3.141592654,PI4=4.0*PI)

	RHOR=ETA*6./PI/8.0

!	hard-sphere term from the Carnahan-Starling Equation state
	AMHS=ETA*(8.-9.*ETA+3.*ETA**2)/(1.-ETA)**3

!	chain connectivity
	GHS=(1.-0.5*ETA)/(1.-ETA)**3
	DGHS=3.*ETA/(1.-ETA)-ETA/(2.-ETA)
	AMC=(DLOG(GHS)+DGHS)*FLOAT(1-M)/FLOAT(M)

!	overall
	AMU=AMHS+AMC
	RETURN
	END


	SUBROUTINE DFEXATT(Z,DZ,NZ,RHO,DATT)
! Calculate the excess local chemical potential at point Z due to two-Yukawa attraction
	IMPLICIT NONE
	INTEGER NZ,K1
	DOUBLE PRECISION PI,DZ,RHO(NZ),DATT,Z,Z1
	DOUBLE PRECISION RCI,ZZ1,ZY(2),X(2,7),RhoBC,
     &                 cgt10,cgt11,cgt12,cgt13,Rc
	COMMON /C2DYUK/ ZY,X,RhoBC,cgt10,cgt11,cgt12,cgt13,Rc

	PARAMETER (PI=3.141592654)
	
	DATT=0.0

! First half slit
	DO 10 K1=1,NZ
	Z1=FLOAT(K1-1)*DZ
	ZZ1=ABS(Z-Z1)/2.
	DATT=DATT+RHO(K1)*RCI(ZZ1)
10      CONTINUE

! Second half slit (image)
	DO 20 K1=NZ,2*NZ-1
	Z1=FLOAT(K1-1)*DZ
	ZZ1=ABS(Z-Z1)/2.
	DATT=DATT+RHO(2*NZ-K1)*RCI(ZZ1)
20      CONTINUE

	DATT=DATT*DZ*8.*PI - RhoBC

	RETURN
	END


	FUNCTION RCI(ZZ1)
C Calcuulate the integration of rc(r) from |z-z'| to infinite	
	IMPLICIT NONE
	INTEGER J
	DOUBLE PRECISION RCI,ZZ1,ZY(2),X(2,7),TEMP1(2)
     &                ,RhoBC,cgt10,cgt11,cgt12,cgt13,Rc
	COMMON /C2DYUK/ ZY,X,RhoBC,cgt10,cgt11,cgt12,cgt13,Rc

	IF (ZZ1.GT.Rc) THEN	   ! Rc: Cut off distance for the LJ Potential

	RCI=0.0

	ELSEIF (ZZ1.GE.1.D0) THEN 

!	DO J=1,2
!	TEMP1(J)=X(J,1)*EXP(ZY(J)*(1.-ZZ1))/ZY(J)  ! Two Yukawa
!	ENDDO

	RCI=cgt11/ZZ1**10 + cgt12/ZZ1**4 - cgt13 ! Original LJ

	ELSE
	DO J=1,2
	TEMP1(J)=-(X(J,1)+X(J,2))*(1.-EXP(ZY(J)*(1.-ZZ1)))/ZY(J)
     &         +X(J,3)*(1.-EXP(ZY(J)*(ZZ1-1.0)))/ZY(J)
     &         +X(J,4)*(1.-ZZ1**5)/5.0
     &         +X(J,5)*(1.-ZZ1**3)/3.0
     &         +X(J,6)*(1.-ZZ1**2)/2.0
     &         +X(J,7)*(1.-ZZ1)
!     &         +X(J,1)/ZY(J)
	ENDDO
	RCI=TEMP1(1)-TEMP1(2)+ cgt11 + cgt12 - cgt13 
	ENDIF
	END

	SUBROUTINE C2INI(TF,SIGMA,DEFF,ETA,Rc0)
C calculate the coefficients of the direct correlation funciton 
C for a given temperature TF=KT/EPSILON and bulk density RHOB=RHO*SIGMA^3

	IMPLICIT NONE
	INTEGER I
	DOUBLE PRECISION PI,TF,SIGMA,DEFF,AK0,ETA,AQ,AL,AS,RhoBC
	DOUBLE PRECISION ZY(2),X(2,7),Z0(2),BETA(2),TEMP(2),Ck0(2),
     &       	 sid6,cgt10,cgt11,cgt12,cgt13,Rc,Rc0
	PARAMETER (PI=3.141592654)
	COMMON /C2DYUK/ ZY,X,RhoBC,cgt10,cgt11,cgt12,cgt13,Rc

  	AK0=2.1714
	Z0(1)=2.9637
	Z0(2)=14.0167	

! cut off distance
	Rc=Rc0*SIGMA/DEFF

	DO I=1,2
	ZY(I)=Z0(I)*DEFF/SIGMA
	BETA(I)=AK0*EXP(Z0(I)*(1.-DEFF/SIGMA))/(DEFF/SIGMA)/TF
	ENDDO

	sid6=(SIGMA/DEFF)**6
	cgt10=-4.*sid6**2/(9.*TF)*(1.-1./Rc**9)
     &     + 4.*sid6/(3.*TF)*(1.-1./Rc**3)
	cgt11=-0.4*sid6**2/TF
	cgt12=sid6/TF
	cgt13=cgt11/Rc**10 + cgt12/RC**4  ! Original LJ

	DO I=1,2
	AL=(1.+0.5*ETA)*ZY(I)+1.+2.*ETA
	AS=(1.-ETA)**2*ZY(I)**3+6.*ETA*(1.-ETA)*ZY(I)**2
     &   +18.*ETA**2*ZY(I)-12.*ETA*(1.+2.*ETA)
	AQ=(AS+12.*ETA*AL*EXP(-ZY(I)))/((1.-ETA)**2*ZY(I)**3)
	TEMP(I)=-BETA(I)/((1.-ETA)**4*ZY(I)**6*AQ**2)
	X(I,1)=BETA(I)
	X(I,2)=TEMP(I)*AS**2
	X(I,3)=TEMP(I)*144.*ETA**2*AL**2
	X(I,4)=-TEMP(I)*12.*ETA**2
     &       *((1.+2.*ETA)**2*ZY(I)**4+(1.-ETA)*(1.+2.*ETA)*ZY(I)**5)
	X(I,5)=TEMP(I)*12.*ETA
     &       *(AS*AL*ZY(I)**2-(1.-ETA)**2*(1.+0.5*ETA)*ZY(I)**6)
	X(I,6)=-TEMP(I)*24.*ETA
     &       *((1.+2.*ETA)**2*ZY(I)**4
     &       +(1.-ETA)*(1.+2.*ETA)*ZY(I)**5)
	X(I,7)=TEMP(I)*24.*ETA*AS*AL

! Calculate the integration for the bulk term	
      Ck0(I)=(X(I,1)+X(I,2))*(exp(ZY(i))-ZY(i)-1.d0)/ZY(i)**2
     &      + X(I,3)*(exp(-ZY(i))+ZY(i)-1.d0)/ZY(i)**2
     &      + X(I,4)/6.+X(I,5)/4.+X(I,6)/3.+X(I,7)/2.
!     &     + X(I,1)*(ZY(i)+1)/ZY(i)**2    ! Original Two Yukawa potential

	ENDDO
	RhoBC=24.*ETA*(CK0(1)-CK0(2)+cgt10)
	RETURN
	END

	SUBROUTINE EXTP(Z,PHI)
	IMPLICIT NONE
        DOUBLE PRECISION Z,PHI,Z0,ZH,TF,DEFF
        COMMON /HARDWALL/ Z0,ZH,TF,DEFF
	IF ((Z.LT.Z0).OR.(Z.GT.ZH-Z0)) THEN
	PHI=1.D10
	ELSE
	PHI=0.0
	ENDIF
	END

        SUBROUTINE EXTPAT(Z,PHI)
        IMPLICIT NONE
        DOUBLE PRECISION Z0,ZH,Z,Z1,DELTA,
     &         TF,SIGW,EPSW,TEMP1,TEMP2,PHI,DEFF
        Parameter (DELTA=0.7071,SIGW=1.d0,epsw=6.283)

        COMMON /HARDWALL/ Z0,ZH,TF,DEFF

	IF(Z.LE.0.5d0) THEN
	PHI=1.d10
	ELSE
        Z1=Z*DEFF/2.0
        TEMP1=0.4*(SIGW/Z1)**10-(SIGW/Z1)**4
     &      -SIGW**4/(3.*DELTA*(0.61*DELTA+Z1)**3)

	Z1=(ZH-Z)*DEFF/2.0
        TEMP2=0.4*(SIGW/Z1)**10-(SIGW/Z1)**4
     &      -SIGW**4/(3.*DELTA*(0.61*DELTA+Z1)**3)

        PHI=epsw*(TEMP1+TEMP2)/TF
	ENDIF
        END


