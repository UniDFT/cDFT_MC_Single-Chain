C********************************************************
C      DFT for Alkanes confined in slit pore
C
C
C********************************************************

	PROGRAM DIBLK_SLIT

	IMPLICIT NONE
	INTEGER I,J,K,M,MA,NZ,ICODE,ITER
	REAL PI,DZ,XB(2),F,ZH,ETA,Y,SADD,diam
	PARAMETER (PI=3.141592654)

	PARAMETER (ICODE=0,M=4,MA=2,NZ=501,DZ=0.02,F=0.01)
	REAL RHO(2,NZ),ROA,ROB,RO(M,NZ),RHO1(2,NZ),RD(NZ),SUMP(2)
	REAL FL(M,NZ),FR(M,NZ),FF(M,NZ)

	ETA=0.3  ! overall packing fraction of segemnts A and B in the bulk
	Diam=1.0
	XB(1)=FLOAT(MA)/FLOAT(M) ! mole fractions of A,B segments
	XB(2)=FLOAT(M-MA)/FLOAT(M)

! overall reduced number densities for segments A and B, rho*radius**3 
	ROA=XB(1)*ETA*6./PI/diam**3
	ROB=XB(2)*ETA*6./PI/diam**3

	ZH=2.*(NZ-1)*DZ	! pore width in units of radius 20 radius = 10 diameter

	ITER=1  ! start Picard iteration

	IF (ICODE.EQ.0) THEN
! initialize the density profiles
	DO I=1,NZ
	RHO(1,I)=ROA			! use the average density as the initial guess
	RHO(2,I)=ROB			
	ENDDO
	ELSE
!  read the initial density profileS from a file
	OPEN(2,FILE='blk.ini')
	DO I=1,NZ
	READ(2,*) Y,RHO(1,I),RHO(2,I)	!Y in units of diameter
	RHO(1,I)=RHO(1,I)*ROA
	RHO(2,I)=RHO(2,I)*ROB		
	ENDDO
 	CLOSE(2)
	ENDIF


!  calculate the recurrence functions and the effective external potential 
1	CALL FFLR(M,MA,ETA,NZ,DZ,FL,FR,FF,RHO,diam)

!  segmental densities,RO(K,I) 
	DO 10 K=1,M
	DO 10 I=1,NZ
!left, right Green fucntions and self-consistent field
  	RO(K,I)=FL(K,I)*FR(K,I)*FF(K,I)		   
!	WRITE(*,*) RO(K,I),FL(K,I),FR(K,I),FF(2,I)
10	CONTINUE


	
!  overall densities of A and B segments at each point
      DO 20 I=1,NZ

	RHO1(1,I)=0.0
	DO K=1,MA
	RHO1(1,I)=RHO1(1,I)+RO(K,I)
      ENDDO

	RHO1(2,I)=0.0
	DO K=MA+1,M
	RHO1(2,I)=RHO1(2,I)+RO(K,I)
      ENDDO

!   density difference
	RD(I)=0
	DO J=1,2
	RD(I)=RD(I)+ABS(RHO1(J,I)-RHO(J,I))
	ENDDO
20	CONTINUE

!   update the density profiles
	SADD=0.0
	DO 30 I=1,NZ
	SADD=SADD+ABS(RD(I))												    
	DO J=1,2
	RHO(J,I)=RHO(J,I)*(1.-F)+RHO1(J,I)*F
	ENDDO
30	CONTINUE

	ITER=ITER+1
!  display intermediate results
	WRITE(*,*) ITER,SADD,RHO(1,1),RHO1(1,1),RD(1)
C	PAUSE

	IF (MOD(ITER,10).EQ.0) THEN
	OPEN(2,FILE='blk.ini')
	DO 40 I=1,NZ
	Y=DZ*FLOAT(I-1)/diam		! in units of diameter (2)
	WRITE(2,41) Y,RHO(1,I)/ROA,RHO(2,I)/ROB	
41	FORMAT(4X,F8.4,4X,F10.7,4X,F10.7)
40	CONTINUE
	CLOSE(2)
	ENDIF

! convergence criteria
	DO I=1,NZ
	IF(ABS(RD(I)).GE.1.E-3) GOTO 1
	ENDDO

c  the average total packing fraction inside the pore
	DO 50 J=1,2
	SUMP(J)=(RHO(J,1)+RHO(J,NZ))/2
	DO I=2,NZ-1
	SUMP(J)=SUMP(J)+RHO(J,I)
	ENDDO
	SUMP(J)=SUMP(J)*(diam**3)*PI/6./FLOAT(NZ-1)
50	CONTINUE

! output results
	OPEN(2,FILE='blk.ini')
	OPEN(3,FILE='blk.out')
	WRITE(3,*)'ETA_AV=',SUMP(1),SUMP(2)

	DO 60 I=1,NZ
	IF (I.LT.50.OR.MOD(I,10).EQ.0) THEN
	Y=DZ*FLOAT(I-1)/diam		! Y is in the units of DIAMETER!
	RHO(1,I)=RHO(1,I)/ROA	! REDUCED density
	RHO(2,I)=RHO(2,I)/ROB
	WRITE(2,41) Y,RHO(1,I),RHO(2,I)	
	IF (MOD(I,5).EQ.0) WRITE(3,61)Y,RHO(1,I),RHO(2,I)
	ELSE
	ENDIF
61	FORMAT(4X,F8.4,4X,F10.7,4X,F10.7)
60	CONTINUE
	CLOSE(3)
	CLOSE(2)

	STOP
	END

	SUBROUTINE FFLR(M,MA,ETA,NZ,DZ,FL,FR,FF,RHO,diam)
!	calculate the recurrence functions and 
!	the self-consistent field inside the hard slit pore
!     Unlike the density profile, FL,FR,FF must be calculated 
!	within the entire slit width=2*(NZ-1)*DZ

!	IR: 	   hard-sphere diameter in units DZ, IR=2/DZ  
!	SUM1:      effective external potential 
!	Z1:        position in z-direction
!	FL:		   left Green function
!	FR:		   right Gren function
!	FF:		   Boltzmann term due to the self-consistent potential

! SUBROUTINES USED: SFM1D(M,MA,ETA,DZ,NZ,Z1,SUM1,RHO), SELF-CONSISTENT MEAN-FIELD AT Z1

	IMPLICIT NONE
	INTEGER I,J,IR,K,M,MA,NZ,KIN,KIP
	REAL DZ,RHO(2,NZ),Z1,SUM1(M),S,ETA,diam
	REAL FL(M,NZ),FR(M,NZ),FF(M,NZ)

	DO 10 I=1,NZ
	FL(1,I)=1.0
	FR(M,I)=1.0
	Z1=DZ*FLOAT(I-1)
	CALL SFM1D(M,MA,ETA,DZ,NZ,Z1,SUM1,RHO, diam)
	DO J=1,M
	FF(J,I)=EXP(-SUM1(J))
	ENDDO
10    CONTINUE

! calculate Green functions

	IR=INT(diam/DZ)

! left Green functions
	DO 30 J=2,M
	DO 20 I=1,NZ

! set integration boundaries
	KIN=I-IR
	KIP=I+IR
	IF(KIN.LE.0) KIN=1
	IF(KIP.GT.(2*NZ-1)) KIP=2*NZ-1

! integrate using the trapzoidal approximation
	IF (KIP.LT.NZ) THEN
	S=(FL(J-1,KIN)*FF(J-1,KIN)
     &  +FL(J-1,KIP)*FF(J-1,KIP))/2.
	ELSE
	S=(FL(J-1,KIN)*FF(J-1,KIN)
     &  +FL(J-1,2*NZ-KIP)*FF(J-1,2*NZ-KIP))/2.
	ENDIF

	DO K=KIN+1,KIP-1
	IF (K.LT.NZ) THEN
	S=S+FL(J-1,K)*FF(J-1,K)
	ELSE
	S=S+FL(J-1,2*NZ-K)*FF(J-1,2*NZ-K)
	ENDIF
	ENDDO

	FL(J,I)=S*DZ/(2.0*diam)
20	CONTINUE
30	CONTINUE


! right Green functions
	DO 50 J=M-1,1,-1
	DO 40 I=1,NZ

! set integration boundaries
	KIN=I-IR
	KIP=I+IR
	IF(KIN.LE.0) KIN=1
	IF(KIP.GT.(2*NZ-1)) KIP=2*NZ-1

! integrate using the trapzoidal approximation
	IF (KIP.LT.NZ) THEN
	S=(FR(J+1,KIN)*FF(J+1,KIN)
     &  +FR(J+1,KIP)*FF(J+1,KIP))/2.
	ELSE
	S=(FR(J+1,KIN)*FF(J+1,KIN)
     &  +FR(J+1,2*NZ-KIP)*FF(J+1,2*NZ-KIP))/2.
	ENDIF

	DO K=KIN+1,KIP-1
	IF (K.LT.NZ) THEN
	S=S+FR(J+1,K)*FF(J+1,K)
	ELSE
	S=S+FR(J+1,2*NZ-K)*FF(J+1,2*NZ-K)
	ENDIF
	ENDDO

	FR(J,I)=S*DZ/(2.0*diam)
40	CONTINUE
50	CONTINUE

	RETURN
	END



	SUBROUTINE SFM1D(M,MA,ETA,DZ,NZ,Z,SUM1,RHO,diam)
! THE SELF-CONSISTENT-MEAN FIELDS OF M INDIVIDUAL SEGMENTS 
!  AT POSITION Z FROM THE SQUARE-WELL TANGENT CHAIN SYSTEM
!	DZ:     step length in numerical integration 
!     Z1:		integration variable, from Z-R to Z+R 
!      (because of the definitions of the weighted densities) 

!	SUBROUTNES USED:
!		CPS(M,ETA,AMU), bulk chemical potential
!		PHI_E(M,MA,Z,PHI),  external potential
!		FMN1D(M,DZ,NZ,Z1,DFN,DFA,RHOT), MFMT functions
!		CPEVDW(M,MA,Z,DZ,NZ,RHO,CPVDW), van der Waals excess chemical potential
	IMPLICIT NONE
	INTEGER M,MA,NZ,J
	REAL DZ,ETA,AMU,Z,CPHSC,diam
	REAL RHO(2,NZ),SUM1(M),CPVDW(M),PHI(M)

!  calculate the reduced chemical potential PER SEGMENT in the bulk phase, AMU
	CALL CPS(M,MA,ETA,AMU,diam)
!	WRITE(*,*) AMU*FLOAT(M)

!  calculate external potential
	CALL PHI_E(M,MA,Z,PHI)

!  calculate the excess chemical potential due to chain connectivity and hard-sphere interactions 
	CALL CPEHSC(M,DZ,NZ,Z,RHO,CPHSC,diam)


! calculate the excess chemical potential due to van der Waals attractions
	CALL CPEVDW(M,MA,Z,DZ,NZ,RHO,CPVDW)

! overall density field
	DO J=1,MA
	SUM1(J)=CPHSC+CPVDW(J)+PHI(J)-AMU
	ENDDO

	DO J=MA+1,M
	SUM1(J)=CPHSC+CPVDW(J)+PHI(J)-AMU
	ENDDO

	RETURN
	END	
	
	SUBROUTINE CPEHSC(M,DZ,NZ,Z,RHO,CPHSC, diam) 
	IMPLICIT NONE
	INTEGER M,NP,NZ,I
	REAL PI,SF21,SF31,SF21VI,DF21,DF31,DF21VI,Z,Z1,DZ,CPHSC
	REAL DFN(3),DFA(3),RHO(2,NZ),RHOT(NZ),diam
	PARAMETER (PI=3.141592654)
!  calculate the excess chemical potential due to chain connectivity and hard-sphere interactions 
!  overall density
!	SUBROUTNES USED:
!		FMN1D(M,DZ,NZ,Z1,DFN,DFA,RHOT), MFMT functions
	DO I=1,NZ
	RHOT(I)=RHO(1,I)+RHO(2,I)
	ENDDO

!	the number of integration points 	
	NP=1+INT(diam/DZ)

!	sum zero
	SF21=0.0
	SF31=0.0
	SF21VI=0.0

	DO 10 I=1,NP
	Z1=Z-diam/2.+DZ*FLOAT(I-1)

! calculate the derivatives of the free energy density wrt weighted densities
!    n2, n3, nv2
	CALL FMN1D(M,DZ,NZ,Z1,DFN,DFA,RHOT,diam)
!	WRITE(*,*) DZ,Z1,DFN(1),DFA(1)

! integrate using the trapezoidal approximation
	DF21=DFN(1)+DFA(1)
	DF31=(DFN(2)+DFA(2))*((diam/2.)**2-(Z1-Z)**2)

!  explain why Z-Z1 instead of Z1-Z?
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
	SF21=diam*PI*SF21*DZ
	SF31=PI*SF31*DZ
	SF21VI=diam*PI*SF21VI*DZ

	CPHSC=SF21+SF31+SF21VI
	RETURN
	END


	SUBROUTINE FMN1D(M,DZ,NZ,Z1,DFN,DFA,RHOT,diam)
!  calculate the derivatives of the excess free energy density from MFMT 
!  with respect to the weighted densities at position z1 given a density profile rho(z)

	IMPLICIT NONE
	INTEGER I,M,IZ2,NP,NZ
	REAL PI,PI2,PI4,Z1,Z2,DZ,GS11,AN2,AN3,ANV2I,F21,DGN3,DGNV2,ARR
	REAL XI2,XI1,F31,F21VI,DFN3,DGN2,diam
	PARAMETER (PI=3.141592654,PI2=2.*PI,PI4=4.*PI)
	REAL DFN(3),DFA(3),RHOT(NZ)

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
	NP=1+INT(diam/DZ)
	
!	zero the sums for weighted densities n2, n3 and nv2
	AN2=0.0
	AN3=0.0
	ANV2I=0.0

!	Z2: integration variable, from Z1-1 to Z1+1	
!	see Eq.(13) in JCP, V117,10156(2002). 
	DO I=1,NP
	Z2=Z1-diam/2.0+DZ*FLOAT(I-1)
	IZ2=INT(Z2/DZ)+1

! impose the boundary conditions and system symmetry 
	IF ((IZ2.GT.(2*NZ-1)).OR.(IZ2.LT.1)) THEN
	F21=0.0
	ELSEIF (IZ2.LT.NZ) THEN
   	F21=RHOT(IZ2)
	ELSE
	F21=RHOT(2*NZ-IZ2)
	ENDIF

!	end points
	IF(I.EQ.1.OR.I.EQ.NP)THEN
	F21=F21*0.5
	ENDIF

	F31=F21*((diam/2.)**2-(Z1-Z2)**2)

!	nv2 is a vector, explain why Z1-Z2 instead of Z2-Z1
	F21VI=F21*(Z1-Z2)

	AN2=AN2+F21
	AN3=AN3+F31
	ANV2I=ANV2I+F21VI
	ENDDO
	AN2=diam*PI*DZ*AN2
	AN3=PI*DZ*AN3
	ANV2I=2.*PI*DZ*ANV2I
	
!	assume the system is ideal when the local overall density is small
	IF(AN3.LE.1.E-5)THEN
	DFN(1)=0.
	DFN(2)=0.
	DFN(3)=0.
	DFA(1)=0.
	DFA(2)=0.
	DFA(3)=0.
	RETURN
	ELSEIF (AN3.GT.1.0) THEN
	AN3=0.74
	ENDIF

C	WRITE(*,*) AN2,AN3,ANV2I

!	calculate the excess free energy density and its derivatives wrt 
!     weighted densities using the MFMT
!	dfn(1) = dphi/dn2, dfn(2) = dphi/dn3, dfn(3) = dphi/dnv2  	

!  xi functions
	XI1=1.-ANV2I**2/AN2**2
	XI2=1.-3.*ANV2I**2/AN2**2

	DFN3=ALOG(1.-AN3)/(18.*PI*AN3**3)+1./(36.*PI*AN3**2*(1.-AN3))
     &    +(1.-3.*AN3)/(36.*PI*(1.-AN3)**3*AN3**2)

	DFN(1)=-ALOG(1.-AN3)/PI/diam**2
     &       +AN2/PI/diam/(1.-AN3)
     & +(AN2**2-ANV2I**2)/(12.*PI*AN3)*(ALOG(1.-AN3)/AN3+1./(1.-AN3)**2)

	DFN(2)=AN2/PI/diam**2/(1.-AN3)+AN2**2*XI1/PI/(2*diam)/(1.-AN3)**2
     &       -AN2**3*XI2*DFN3

	DFN(3)=-ANV2I/PI/diam/(1.-AN3)
     &       -AN2*ANV2I/(6.*PI*AN3)*(ALOG(1.-AN3)/AN3+1./(1.-AN3)**2)

! derivatives of the excess free energy due to the hard-sphere chain wrt weighted densities
!   GS11, the contact value of the radial distribution function of uniform hard spheres
!	see Eq.(18) in JCP, V117,10156(2002). 

	GS11=1./(1.-AN3)
     &    +0.25*diam*AN2*XI1/(1.-AN3)**2
     &    +AN2**2*diam**2*XI1/(72.*(1.-AN3)**3)

!	derivatives of GS11 wrt n2, n3 and nv2
	DGN2=0.25*diam*(2.-XI1)/(1.-AN3)**2+AN2*diam**2/(1.-AN3)**3/36.0
	DGN3=1./(1.-AN3)**2+0.5*diam*AN2*XI1/(1.-AN3)**3
     &       +AN2**2*diam**2*XI1/(1.-AN3)**4/24.0
	DGNV2=-0.5*diam*ANV2I/AN2/(1.-AN3)**2
     &      -diam**2*ANV2I/(1.-AN3)**3/36.0

	ARR=(1./FLOAT(M)-1.)/PI/diam**2

	DFA(1)=ARR*(XI1*AN2/GS11*DGN2+(2.-XI1)*ALOG(GS11))
	DFA(2)=ARR*(XI1*AN2/GS11*DGN3) 
	DFA(3)=ARR*(XI1*AN2/GS11*DGNV2-2.*ANV2I/AN2*ALOG(GS11))

C	WRITE(*,*) DFN,DFA
	RETURN
	END

!	calculate the van der Waals self-consistent-mean-field potentials
!     at position Z for a binary mixture with	the
!	square-well potential for attraction between segments (1,2)


	SUBROUTINE CPEVDW(M,MA,Z,DZ,NZ,RHO,CPVDW)
	IMPLICIT NONE
	REAL Z,DZ,XH,XL,ZK,FVDW,PI,FRHO
	INTEGER	M,MA,IZ,IH,IL,I,J,K,NZ
	REAL LAMDA(2,2),SIGMA(2,2),EPS(2,2),CPMF(2),RHO(2,NZ),CPVDW(M)
	PARAMETER (PI=3.141592654)
	COMMON /VDW/LAMDA,EPS,SIGMA
	EXTERNAL FRHO

	IZ=INT(Z/DZ)

	DO 20 I=1,2

	CPMF(I)=0.0

	DO 10 J=1,2

! integration from z-sigma(i,j) to z+sigma(i,j)
	IL=IZ-INT(SIGMA(I,J)/DZ)
	IH=IZ+INT(SIGMA(I,J)/DZ)

	FVDW=(FRHO(J,IL,NZ,RHO)+FRHO(J,IH,NZ,RHO))/2.
	DO K=IL+1,IH-1
	FVDW=FVDW+FRHO(J,K,NZ,RHO)
	ENDDO
	CPMF(I)=CPMF(I)
     &     +EPS(I,J)*SIGMA(I,J)**2*(LAMDA(I,J)**2-1.)*FVDW

! integration from sigma(i,j) to lamda(i,j)*sigma(i,j)
	XL=SIGMA(I,J)
	XH=LAMDA(I,J)*SIGMA(I,J)
	IL=INT(XL/DZ)
	IH=INT(XH/DZ)

	   FVDW=FRHO(J,IL+IZ,NZ,RHO)*(XH**2-XL**2)/2.
    	DO K=IL+1,IH-1
	  ZK=FLOAT(K)*DZ
	  FVDW=FVDW+FRHO(J,K+IZ,NZ,RHO)*(XH**2-ZK**2)
	  ENDDO

! integration from -lamda(i,j)*sigma(i,j) to -sigma(i,j) 
	  FVDW=FVDW+FRHO(J,-IL+IZ,NZ,RHO)*(XH**2-XL**2)/2.
	  DO K=-IH+1,-IL-1
  	  ZK=FLOAT(K)*DZ
	  FVDW=FVDW+FRHO(J,K+IZ,NZ,RHO)*(XH**2-ZK**2)
	  ENDDO

	  CPMF(I)=CPMF(I)+EPS(I,J)*FVDW

10	  CONTINUE

	  CPMF(I)=-PI*CPMF(I)*DZ

20	  CONTINUE

	  DO I=1,MA
	  CPVDW(I)=CPMF(1)
	  ENDDO

	  DO I=MA+1,M
	  CPVDW(I)=CPMF(2)
	  ENDDO

	  RETURN
	  END

	  FUNCTION FRHO(I,IZ,NZ,RHO)
	  IMPLICIT NONE
	  INTEGER I,IZ,NZ
	  REAL FRHO,RHO(2,NZ)
	  IF ((IZ.GT.(2*NZ-1)).OR.(IZ.LT.1)) THEN
	  FRHO=0.0
	  ELSEIF (IZ.LT.NZ) THEN
   	  FRHO=RHO(I,IZ)
	  ELSE
	  FRHO=RHO(I,2*NZ-IZ)
	  ENDIF
	  END

	  SUBROUTINE CPS(M,MA,ETA,AMU,diam)

!CPS(M, ROSS, ETA, AMU) calculates the chemical potential PER SEGMENT of a one-component 
!hard-sphere-chain fluid. Each chain consists of M identical segments
!	and the total reduced number density of segments is ROSS. 
!
! KEY PARAMETERS:
!	M      the number of segments per chain
!	ETA	   packing fraction of hard sphere segments
!	ROSS   the reduced density of segments, rho*sigma**3/8 
!	AMR	   the ideal gas chemical potential	
!	AMHS   the excess chemical potential due to hard sphere	collision
!	AMC    the excess chemical potential due to chain formation  
!	AMU    the total reduced chemical potential, mu/kBT/M  
!	GHS    contact value of the radial distribution function for hard spheres 
!	DGHS   derivative of LN(GHS) with respect to density

!	lamda,sigma,eps are square-well parameters

	IMPLICIT NONE
	REAL VDWA,PI,PI4,ROSS,ETA,AMR,AMHS,GHS,DGHS,AMC,AMVDW,AMU,diam
	INTEGER	I,J,M,MA
	REAL LAMDA(2,2),SIGMA(2,2),EPS(2,2),XB(2)
	PARAMETER (PI=3.141592654,PI4=4.0*PI)
	COMMON /VDW/LAMDA,EPS,SIGMA

	XB(1)=FLOAT(MA)/FLOAT(M) ! mole fractions of A,B segments
	XB(2)=FLOAT(M-MA)/FLOAT(M)

	ROSS=ETA*6./PI/diam**3

!	ideal-gas term
	AMR=ALOG(ROSS/FLOAT(M))/FLOAT(M)

!	hard-sphere term from the Carnahan-Starling equation state
	AMHS=ETA*(8.-9.*ETA+3.*ETA**2)/(1.-ETA)**3

!	chain connectivity
	GHS=(1.-0.5*ETA)/(1.-ETA)**3
	DGHS=3.*ETA/(1.-ETA)-ETA/(2.-ETA)
	AMC=(ALOG(GHS)+DGHS)*FLOAT(1-M)/FLOAT(M)

!	van der Waals mean-field term
!	VDWA: van der Waals energy parameter
	VDWA=0.
	DO I=1,2
	DO J=1,2
	VDWA=VDWA + XB(I)*XB(J)* 
     &  (-4.*PI/3.)*EPS(I,J)*SIGMA(I,J)**3*(LAMDA(I,J)**3-1.)
	ENDDO
	ENDDO

	AMVDW=ROSS*VDWA

!	overall
	AMU=AMR+AMHS+AMC+AMVDW
	RETURN
	END

	SUBROUTINE PHI_E(M,MA,Z,PHI)
	IMPLICIT NONE
	INTEGER M,MA,I
	REAL Z,PHI(M)
	REAL WALL(2),WLAM(2)
! wall potential, wall: enregy,kT; wlam: attraction distance, in units of radius 
	DATA WALL/-0.,0./
	DATA WLAM/0.5,0.5/

	DO I=1,MA
	IF ((Z.LT.0.0).OR.(Z.GT.WLAM(1))) THEN
	PHI(I)=0.0
	ELSE
	PHI(I)=WALL(1)
	ENDIF
	ENDDO

	DO I=MA+1,M
	IF ((Z.LT.0.0).OR.(Z.GT.WLAM(1))) THEN
	PHI(I)=0.0
	ELSE
	PHI(I)=WALL(2)
	ENDIF
	ENDDO
	RETURN
	END

	BLOCKDATA
	REAL LAMDA(2,2),SIGMA(2,2),EPS(2,2)
	COMMON /VDW/LAMDA,EPS,SIGMA
	DATA LAMDA/1.2,1.2,1.2,1.2/
	DATA EPS/0.0,0.0,0.0,0.0/
	DATA SIGMA/1.0,1.0,1.0,1.0/
	END
