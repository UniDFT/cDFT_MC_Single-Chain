!$DEBUG
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
	PARAMETER (M=3,NZ=501,numb=100000)
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
       	  	  
	   confg(1,1,1)=0.0
	   confg(1,1,2)=0.0
	   confg(1,1,3)=wh/2.0
	
	    Do 80 i=2,m

55        sita=2.0 * pi *ran2(idum)
	    fai=pi*ran2(idum)
	  
	    xr=confg(1,i-1,1)+cos(sita)*sin(fai)
	    yr=confg(1,i-1,2)+sin(sita)*sin(fai)
	    zr=confg(1,i-1,3)+cos(fai)
	  
!	   if(zr.gt.wh. or.zr.lt.0)	goto 55
	    confg(1,i,1)=xr-boxl*int(2.0*xr /boxl)
	    confg(1,i,2)=yr-boxl*int(2.0*yr /boxl)
	    confg(1,i,3)=zr

80      continue
	           
        write(*,*)   'Chain move starts (NVT MC)......'
	  icy=2
      Do 800 ii=2,ncycle

        idice=int(5.0*ran2(idum))+1
        if(idice.eq.1) then
          goto 200
        else if(idice.eq.2) then
          goto 300
	  else if(idice .eq.3) then
	    goto 400
	  else if(idice .eq.4) then
	    goto 500
        endif

c       Self_avoid walk or moves of chains
c    xi,yi, zi are old coordinates; xr; yr; zr are the generated coordinates
C       nmm=nmm+1

100     xi=confg(icy-1,m,1)
        yi=confg(icy-1,m,2)
	  zi=confg(icy-1,m,3)

        sita=2.0*pi*ran2(idum)
	  fai=pi*ran2(idum)

	  xr=xi+cos(sita)*sin(fai)
	  yr=yi+sin(sita)*sin(fai)
	  zr=zi+cos(fai)
	  
	  if(zr.gt.0.and.zr.lt.wh) then
	    goto 105
	  else
	    do ik=2,m
		 xtst=confg(icy-1,ik,3)
		 if(xtst.gt.0.and.xtst.lt.wh) goto 105
		enddo
		goto 600   
	  endif

105	  xr=xr-boxl*int(2.0*xr /boxl)
	  yr=yr-boxl*int(2.0*yr /boxl)
	  
	  do j=2,m
	    confg(icy,j-1,1)=confg(icy-1,j,1)
	    confg(icy,j-1,2)=confg(icy-1,j,2)
	    confg(icy,j-1,3)=confg(icy-1,j,3)
	   enddo

	    confg(icy,m,1)=xr
	    confg(icy,m,2)=yr
	    confg(icy,m,3)=zr
	    icy=icy+1	  
        goto 600

C	  inmm=inmm+1
        
200     xi=confg(icy-1,1,1)
        yi=confg(icy-1,1,2)
	  zi=confg(icy-1,1,3)
	          
        sita=2.0*pi*ran2(idum)
	  fai=pi*ran2(idum)

	  xr=xi+cos(sita)*sin(fai)
	  yr=yi+sin(sita)*sin(fai)
	  zr=zi+cos(fai)

       if(zr.gt.0.and.zr.lt.wh) then
	    goto 205
	  else
	    do ik=2,m
		 xtst=confg(icy-1,ik,3)
		 if(xtst.gt.0.and.xtst.lt.wh) goto 205
		enddo
		goto 600   
	  endif
	  

205	  xr=xr-boxl*int(2.0*xr /boxl)
	  yr=yr-boxl*int(2.0*yr /boxl)
	 	  	  
	   do j=m,2,-1
	    confg(icy,j,1)=confg(icy-1,j-1,1)
	    confg(icy,j,2)=confg(icy-1,j-1,2)
	    confg(icy,j,3)=confg(icy-1,j-1,3)
	   enddo
	
          confg(icy,1,1)=xr
	    confg(icy,1,2)=yr
	    confg(icy,1,3)=zr
	  	icy=icy+1
		  
        goto 600

C        plane moves

 
300	 Do jj=1,m
	  confg(icy,jj,1)=confg(icy-1,jj,1)+dmove*(2.0*ran2(idum)-1.0)
	  confg(icy,jj,2)=confg(icy-1,jj,2)
	  confg(icy,jj,3)=confg(icy-1,jj,3)
	 Enddo

	  icy=icy+1
        goto 600

400    Do jj=1,m
        confg(icy,jj,2)=confg(icy-1,jj,2)+dmove*(2.0*ran2(idum)-1.0)
	  confg(icy,jj,1)=confg(icy-1,jj,1)
	  confg(icy,jj,3)=confg(icy-1,jj,3)
	 Enddo
	 
	  icy=icy+1
	  goto 600

500     Do jj=1,m
         confg(icy,jj,3)=confg(icy-1,jj,3)+dmove*(2.0*ran2(idum)-1.0)
	   confg(icy,jj,1)=confg(icy-1,jj,1)
	   confg(icy,jj,2)=confg(icy-1,jj,2)
	  Enddo
	  
	  do ik=1,m
	     xtt=confg(icy,ik,3)
	    if(xtt.lt.wh.and.xtt.gt.0) goto 501
	  enddo
	  goto 600

501	   icy=icy+1 

c      if(icy.gt.nequl) then
          
c          ngd=ngd+1
c          call sumation(u1,uff,ufw,rostar,ucm,ucf,ucw,nk)
c          ngl=ngl+1
c          if(ngl.eq.nglong) then
c            call radial
c            ncc=ncc+1
c            call sumzero
c            ngl=0
c          endif
c        endif

600     rrm1=float(icy)/float(ii)
        if(mod(icy,nsave).eq.0) then
	  write(*,*) 'icy=',icy,'ii=',ii,'ratio=',rrm1
	  endif
	if(icy.gt.numb) goto 810

800     continue
       
810	  open(unit=9,file='config0.dat',status='unknown')
        
	    Do 610 jj=1,numb
	    do 610 kk=1,m
		  write(9,999) confg(jj,kk,1), confg(jj,kk,2), confg(jj,kk,3)
610	    continue   
999     format(3f12.6)
	  close(9)



! model parameters
	ETA=0.1
	ROSS=ETA*6./PI/8.0

! numerical parameters
	ICODE=0
	DZ=0.02
	F=0.05    

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
	  write(*,*)'eta=',eta, 'ross=',ross,'zh=',zh,'\'

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
	write(4,51)y,3.0*ros(1,i)/roav, 3.0*ros(2,i)/roav
51	FORMAT(1X,F8.4,4X,2F10.6)
50	CONTINUE
	CLOSE(3)

	STOP
	END

*   ================= subroutine program =======================

	subroutine aver(m,nz,dz,amu,ros,rho,confg)
	 parameter(numb=100000)
	 real confg(numb,m,3),ros(m,2*nz-1),ff(2*nz-1),rho(nz)
	 double precision rosum1(2*nz-1),rosum2(2*nz-1),rosum3(2*nz-1)
	 integer ict1(2*nz-1),ict2(2*nz-1),ict3(2*nz-1)
	 
	 
	 Do 10 i=1,2*nz-2
	  z1=dz*float(i-1)
	  call sfm1d(m,dz,nz,z1,sum1,rho)
	  ff(i)=exp(-sum1+amu)
!	  write(*,*) i, '  ', ff(i)
10     continue

	 do ii=1,2*nz-2
	 ict1(ii)=0
	 ict2(ii)=0
	 ict3(ii)=0
	 
	 rosum1(ii)=0.0
	 rosum2(ii)=0.0
	 rosum3(ii)=0.0
	 enddo
     
	 ic1=0
	 ic2=0
	 ic3=0

	 do 200 i=1,2*nz-2

	 do 100 j=1,numb
	  	  
	  zz1=confg(j,1,3)*2.0
	  zz2=confg(j,2,3)*2.0
	  zz3=confg(j,3,3)*2.0
	  iz1=int(zz1/dz)+1
	  iz2=int(zz2/dz)+1
	  iz3=int(zz3/dz)+1 

	  if(iz1.eq.i) then
	   if(iz2.le.0.or.iz2.ge.2*nz-2) goto 202
	   if(iz3.le.0.or.iz3.ge.2*nz-2) goto 202
	
         rosum1(i)=rosum1(i)+ff(iz1)*ff(iz2)*ff(iz3)
202	   ict1(i)=ict1(i)+1
	  endif

	  if(iz2.eq.i)then
          if(iz1.le.0.or.iz1.ge.2*nz-2) goto 204
	    if(iz3.le.0.or.iz3.ge.2*nz-2) goto 204
		
	 	rosum2(i)=rosum2(i)+ff(iz1)*ff(iz2)*ff(iz3)
204	    ict2(i)=ict2(i)+1
         endif
	 
	  if(iz3.eq.i) then
	    if(iz1.le.0.or.iz1.ge.2*nz-2) goto 206
	    if(iz2.le.0.or.iz2.ge.2*nz-2) goto 206

	    rosum3(i)=rosum3(i)+ff(iz1)*ff(iz2)*ff(iz3)
206	    ict3(i)=ict3(i)+1
	  endif

100    continue
	  ros(1,i)=rosum1(i)/float(ict1(i))
	  ros(2,i)=rosum2(i)/float(ict2(i))
	  ros(3,i)=rosum3(i)/float(ict3(i))

	 ic1=ic1+ict1(i)
	 ic2=ic2+ict2(i)
	 ic3=ic3+ict3(i)
!	 write(*,*)i, ' ','ff9i0=',ff(i)
200     continue
	 
!	 if(ic1.ne.numb) then
!	 write(*,*) 'ic1=',ic1,'ic2=',ic2,'ic3=',ic3,ic1+ic2+ic3
!	 stop 'ic1 wrong'
!	 endif

!	 if(ic2.ne.numb) stop 'ic2 wrong'
!	 if(ic3.ne.numb) stop 'ic3 wrong' 

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

