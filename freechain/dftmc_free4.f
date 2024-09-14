!   ****************************************************************
!    * A combination of the DFT and MC for claculation of         *
!    * microstructure of fluids (Free chian)                      *
!!   *                           Author:Dapeng Cao                *
!    *                            April 19,2004                   *
!   ****************************************************************
	

      PROGRAM SLIT_DFT_MC

!	THIS PROGRAM IS FOR CALCULATING THE SEGMENTAL DENSITIES OF ONE-COMPONENT 
!     HARD-SPHERE CHAINS CONFINED IN A HARD SLIT-LIKE PORE
!	
!	PAPER: JCP,V117,2356 (2002)

!	KEY PARAMETERS

!	DZ		step length in the units of hard-sphere radius  
!	ZH		half of the pore width in the units of hard-sphere radius 
!			(because of symmetry, only half of the density profile need to be calculated)
!	NZ		NZ=1 + half pore width/DZ or pore width=2*(NZ-1)*DZ  
!	F		mixing parameter in the Picard iteration
!	M		number of segments per chain 
!	ROSS	reduced number density of segments, RHO*SIGMA**3/8=rho*radius**3
!	ETA		packing fraction of segments  
!	RHO(I)	reduced density profile,rho*radius**3    
!	RE		file to save temporary density profile
!	RO(K,I)	the density profile of site K  
!	RD(I)   the difference in density profiles between two consective iterations
!	SUMP	average packing fraction in the pore
!		
	PARAMETER (PI=3.141592654)
	PARAMETER (M=3,NZ=501,numb=50000)
	REAL RHO(NZ),RHO1(NZ),RD(NZ),ros(m,2*nz) 
c	DOUBLE PRECISION FL(M,2*NZ),FR(M,2*NZ),FF(2*NZ)
      real  confg(numb,m,3)
     	
	data idum/-1/
	          
	
!     MC simulations

        open(unit=14,file='mcstep.dat',status='old')
        read(14,*)wh,boxl,dmove,ncycle,nsave
        close(14)

        open(unit=6,file='con')

  !      pi=4.0d0*datan(1.0d0)
  !      twopi=2.0d0*pi
         
	  write(*,*) '...MC start for generating configuration....'
	   
	   do 70 ik=1,numb
	     	  
	   confg(ik,1,1)=0.0
	   confg(ik,1,2)=0.0
	   confg(ik,1,3)=0.0
	
	    Do 80 i=2,m

81        sita=2.0 * pi *ran2(idum)
	    fai=pi*ran2(idum)
	    xk=2.0*ran2(idum)-1.0 	
!		if(i-2.le.0) then 
		        
	       confg(ik,i,1)=confg(ik,i-1,1)+cos(sita)*sqrt(1.-xk*xk)
	       confg(ik,i,2)=confg(ik,i-1,2)+sin(sita)*sqrt(1.-xk*xk)
	       confg(ik,i,3)=confg(ik,i-1,3)+xk
	       
!	  	else 

!             xix=confg(ik,i-1,1)+cos(sita)*sqrt(1.-xk*xk)
!	       yiy=confg(ik,i-1,2)+sin(sita)*sqrt(1.-xk*xk)
!	       ziz=confg(ik,i-1,3)+xk
	    
!		   do j=1,i-2
!		     xii=xix-confg(ik,j,1)
!		     yii=yiy-confg(ik,j,2)
!		     zii=ziz-confg(ik,j,3)
!		     rii=xii*xii+yii*yii+zii*zii
!		     if(rii.lt.1.0) goto 81
!		   enddo
		   		 
!		   confg(ik,i,1)=xix
!	       confg(ik,i,2)=yiy
!	       confg(ik,i,3)=ziz
		   
!		 endif 
	
80      continue
70      continue
		      
810	  open(unit=9,file='config0.dat',status='unknown')
        
	    Do 610 jj=1,numb
	    do 610 kk=1,m
		  write(9,999) confg(jj,kk,1), confg(jj,kk,2), confg(jj,kk,3)
610	    continue   
999     format(3f12.6)
	  close(9)

!	 MC codes for generating configurations end here and start DFT codes 

! model parameters
	ETA=0.1
	ROSS=ETA*6./PI/8.0

! numerical parameters
	ICODE=0
	DZ=0.02
	F=0.1    

! pore width in units of radius
	ZH=2.*(NZ-1)*DZ
!	WRITE(*,*) ZH,ROSS

! start Picard iteration
	ITER=0

! initialize the density profile
	IF (ICODE.EQ.0) THEN
! use the bulk density as the initial guess
	DO I=1,NZ
	RHO(I)=ROSS
	ENDDO
	ELSE

!  read the initial density profile from a file
	OPEN(2,FILE='RE')
	DO I=1,NZ
	READ(2,*) Y,RHO(I)
	ENDDO
 	CLOSE(2)
	ENDIF

!  calculate the reduced chemical potential PER SEGMENT in the bulk phase, AMU
	CALL CPS(M,ETA,AMU)
!	WRITE(*,*) AMU*FLOAT(M)
 	  
	  write(*,*) '.....start DFT iteration.......'
	  write(*,*)'eta=',eta, 'ross=',ross,'zh=',zh

!  calculate the recurrence functions and the effective external potential 

1	  call aver(m,nz,dz,amu,ros,rho,confg)

c1      CALL FFLR(M,NZ,DZ,AMU,FL,FR,FF,RHO)

!  segmental densities,RO(K,I) 

!  overall density at each point
      DO 20 I=1,NZ
	RHO1(I)=ros(1,i)+ros(2,i)+ros(3,i)
		      
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
	WRITE(*,*) ITER,SADD,RHO(1),RHO1(1),RD(1)

	IF (MOD(ITER,10).EQ.0) THEN
	OPEN(2,FILE='RE')
	DO 40 I=1,NZ
	Y=DZ*FLOAT(I-1)/2.0
	WRITE(2,41)Y,RHO(I)
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
	roav=sump/float(nz-1)
	SUMP=SUMP*8.*PI/6./FLOAT(NZ-1)

! output results
	OPEN(3,FILE='RES1')
	WRITE(3,*)'ETA_AV=',SUMP
	DO 50 I=1,NZ
! Y is in the units of DIAMETER!
	Y=DZ*FLOAT(I-1)/2.0
!	RHO(I)=RHO(I)/ROSS
	WRITE(3,51)Y,RHO(I)/ross, rho(i)/roav
	open(4,file='ros.dat')
	write(4,52)y,3.*ros(1,i)/roav,3.*ros(2,i)/roav,3.*ros(3,i)/roav
51	FORMAT(1X,F8.4,4X,2F10.6)
52    format(1x,f8.4,4x,3f10.6)	 
50	CONTINUE
	CLOSE(3)

	STOP
	END

*   ================= subroutine program =======================

	subroutine aver(m,nz,dz,amu,ros,rho,confg)
	 parameter(numb=50000)
	 real confg(numb,m,3),ros(m,2*nz-1),ff(2*nz-1),rho(nz),zz(m)
	 double precision rosum(m,2*nz-1)
	 integer ict(m,2*nz-1), iz(m)
	 
	 
	 Do 10 i=1,2*nz-2
	  z1=dz*float(i-1)
	  call sfm1d(m,dz,nz,z1,sum1,rho)
	  ff(i)=exp(-sum1+amu)
!	  write(*,*) i, '  ', ff(i)
10     continue

	 do ii=1,2*nz-2
	   do jj=1,m
	     ict(jj,ii)=0
	     rosum(jj,ii)=0.0
	   enddo
	 enddo
     	 	 
	 do 200 i=1,2*nz-2
	  
	    zrz=float(i)*dz

	 do 100 j=1,numb
	  
	   do 50 im=1,m

	      do kj=1,m
	  	   zz(kj)=(confg(j,kj,3)-confg(j,im,3))*2.0+zrz
	       iz(kj)=int(zz(kj)/dz)+1
	   
	       if(kj.ne.im) then
		     if(iz(kj).le.0.or.iz(kj).ge.2*nz-2) goto 50   
	  	   endif
	      enddo
	    
		  fsum=1.0
	      do kk=1,m
	        fsum=fsum*ff(iz(kk))
	      enddo
	    rosum(im,i)=rosum(im,i)+fsum
50	 continue

100    continue
        
	   do km=1,m
	     ros(km,i)=rosum(km,i)/float(numb)
	   enddo

200     continue
	 
	return
	end


	SUBROUTINE FFLR(M,NZ,DZ,AMU,FL,FR,FF,RHO)
!	calculate the recurrence functions and 
!	the self-consistent field inside the hard slit pore
!     Unlike the density profile, FL,FR,FF must be calculated 
!	within the entire slit width=2*(NZ-1)*DZ

!	all other variables are real, why FL, FR, FF should be double precision?
	DOUBLE PRECISION FL(M,2*NZ),FR(M,2*NZ),FF(2*NZ)
	REAL RHO(NZ)

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
	CALL SFM1D(M,DZ,NZ,Z1,SUM1,RHO)
	FF(I)=EXP(-SUM1+AMU)
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
!	WRITE(*,*) I,J,FL(J,I),FF(I)
20	CONTINUE
30	CONTINUE
	RETURN
	END


	SUBROUTINE SFM1D(M,DZ,NZ,Z,SUM1,RHO)
! THE SELF-CONSISTENT-MEAN FIELD FROM THE HARD-SPHERE CHAIN REFERENCE SYSTEM
!	DZ:     step length in numerical integration 
!     Z1:		integration variable, from Z-R to Z+R (because of the definitions of the weighted densities) 

	REAL DFN(3),DFA(3),RHO(NZ)
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
	SF21=2.0*PI*SF21*DZ
	SF31=PI*SF31*DZ
	SF21VI=2.*PI*SF21VI*DZ

	SUM1=SF21+SF31+SF21VI
!	WRITE(*,*) Z,SUM1

	RETURN
	END	
	
	SUBROUTINE FMN1D(M,DZ,NZ,Z1,DFN,DFA,RHO)
!  calculate the derivatives of the excess free energy density  
!  with respect to the weighted densities at position z1 given a density profile rho(z)


	PARAMETER (PI=3.141592654,PI2=2.*PI,PI4=4.*PI)
	REAL DFN(3),DFA(3),RHO(NZ)

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

!	nv2 is a vector, explain why Z1-Z2 instead of Z2-Z1
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
!     weighted densities using the MFMT
!	dfn(1) = dphi/dn2, dfn(2) = dphi/dn3, dfn(3) = dphi/dnv2  	

!  xi functions
	XI1=1.-ANV2I**2/AN2**2
	XI2=1.-3.*ANV2I**2/AN2**2
	DFN3=ALOG(1.-AN3)/(18.*PI*AN3**3)+1./(36.*PI*AN3**2*(1.-AN3))
     &    +(1.-3.*AN3)/(36.*PI*(1.-AN3)**3*AN3**2)

!	dfn(1), Eq.(?) in somewhere
	DFN(1)=-ALOG(1.-AN3)/PI4
     &       +AN2/PI2/(1.-AN3)
     & +(AN2**2-ANV2I**2)/(12.*PI*AN3)*(ALOG(1.-AN3)/AN3+1./(1.-AN3)**2)

!	dfn(2), Eq.(?) in somewhere
	DFN(2)=AN2/PI4/(1.-AN3)+AN2**2*XI1/PI4/(1.-AN3)**2-AN2**3*XI2*DFN3

!	dfn(3), Eq.(?) in somewhere
	DFN(3)=-ANV2I/PI2/(1.-AN3)
     &       -AN2*ANV2I/(6.*PI*AN3)*(ALOG(1.-AN3)/AN3+1./(1.-AN3)**2)

! derivatives of the excess free energy due to the hard-sphere chain wrt weighted densities
!   GS11, the contact value of the radial distribution function of uniform hard spheres
!	see Eq.(18) in JCP, V117,10156(2002). 

	GS11=1./(1.-AN3)
     &    +0.5*AN2*XI1/(1.-AN3)**2
     &    +AN2**2*XI1/(18.*(1.-AN3)**3)

C	WRITE(*,*) GS11
!	derivatives of GS11 wrt n2, n3 and nv2
	DGN2=0.5*(2.-XI1)/(1.-AN3)**2+AN2/(1.-AN3)**3/9.0
	DGN3=1./(1.-AN3)**2+AN2*XI1/(1.-AN3)**3+AN2**2*XI1/(1.-AN3)**4/6.0
	DGNV2=-ANV2I/AN2/(1.-AN3)**2-ANV2I/(1.-AN3)**3/9.0

	ARR=(1./FLOAT(M)-1.)/PI4

	DFA(1)=ARR*(XI1*AN2/GS11*DGN2+(2.-XI1)*ALOG(GS11))
	DFA(2)=ARR*(XI1*AN2/GS11*DGN3) 
	DFA(3)=ARR*(XI1*AN2/GS11*DGNV2-2.*ANV2I/AN2*ALOG(GS11))

C	WRITE(*,*) DFN,DFA
	RETURN
	END


	SUBROUTINE CPS(M,ETA,AMU)

!	CPS(M, ROSS, ETA, AMU) calculates the chemical potential PER SEGMENT of a one-component 
!		hard-sphere-chain fluid. Each chain consists of M identical segments
!		and the total reduced number density of segments is ROSS. 
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

	PARAMETER (PI=3.141592654,PI4=4.0*PI)

	ROSS=ETA*6./PI/8.0

!	ideal-gas term
	AMR=ALOG(ROSS/FLOAT(M))/FLOAT(M)

!	hard-sphere term from the Carnahan-Starling Equation state
	AMHS=ETA*(8.-9.*ETA+3.*ETA**2)/(1.-ETA)**3

!	chain connectivity
	GHS=(1.-0.5*ETA)/(1.-ETA)**3
	DGHS=3.*ETA/(1.-ETA)-ETA/(2.-ETA)
	AMC=(ALOG(GHS)+DGHS)*FLOAT(1-M)/FLOAT(M)

!	overall
	AMU=AMR+AMHS+AMC
	RETURN
	END



	real FUNCTION ran2(idum)

      implicit real*8 (a-h,o-z)

      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.-EPS)
      INTEGER iv(NTAB)
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=dmin1(AM*iy,RNMX)
      return
      END

