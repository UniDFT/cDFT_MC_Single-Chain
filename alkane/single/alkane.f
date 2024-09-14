C***************************************************************
C      *  DFT for Alkanes (LJ Copolymer) confined in slit pore *
C	 * The original code by J.Wu is for SW copolymer         *
C      * Them modified code combines the single chain          *
C      * simulation to DFT    for LJ Copolymer chain           *
C      *      .....Version 2...            By Dapeng CAO       *  
C	 *			                     June 15, 2004	     *
C***************************************************************

	PROGRAM DIBLK_SLIT_LJ

	IMPLICIT NONE
	INTEGER I,J,K,M,MA,NZ,ICODE,ITER,numb,id
	REAL PI,DZ,XB(2),F,ZH,ETA,Y,SADD,diam,cmsum,av_cm
	PARAMETER (PI=3.141592654,numb=10000)
	PARAMETER (ICODE=1,M=4,MA=2,NZ=551,DZ=0.02,F=0.02)
	REAL RHO(2,NZ),ROA,ROB,RHO1(2,NZ),RD(NZ),SUMP(2)
	REAL confg(numb,m,3),ros(m,2*nz-1),rocm(2*nz-1)
!	real ang_num(2*nz-1), rend(2*nz-1)

	ETA=0.1  ! overall packing fraction of segemnts A and B in the bulk
	Diam=3.9
			      
      XB(1)=FLOAT(MA)/FLOAT(M) ! mole fractions of A,B segments
	XB(2)=FLOAT(M-MA)/FLOAT(M)

! overall reduced number densities for segments A and B, rho*radius**3 
	ROA=XB(1)*ETA*6./PI/diam**3
	ROB=XB(2)*ETA*6./PI/diam**3

	ZH=2.*(NZ-1)*DZ	! pore width in units of radius 20 radius = 10 diameter

	ITER=1  ! start Picard iteration
	write(*,*)'zh=',zh,'A;  Diam=',diam

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
!      Generating the configurations from single chain simulation     	
	call simulation(m,numb,confg)
!  calculate the recurrence functions and the effective external potential 
!1	CALL FFLR(M,MA,ETA,NZ,DZ,FL,FR,FF,RHO,diam)
1	Call aver(m,ma,eta,nz,dz,ros,rho,rocm,confg,diam)

!  segmental densities,RO(K,I) 
	DO 10 K=1,M
	DO 10 I=1,NZ
!left, right Green fucntions and self-consistent field
!  	RO(K,I)=FL(K,I)*FR(K,I)*FF(K,I)		   
!	WRITE(*,*) RO(K,I),FL(K,I),FR(K,I),FF(2,I)
10	CONTINUE
 	
!  overall densities of A and B segments at each point
      DO 20 I=1,NZ
!	RHO1(1,I)=0.0
!	DO K=1,MA
	RHO1(1,I)=ROS(1,I)+ROS(M,I)
!      ENDDO

	RHO1(2,I)=0.0
	DO K=2,M-1
	RHO1(2,I)=RHO1(2,I)+ROS(K,I)
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
        I=int(diam/dz)
!  display intermediate results
	WRITE(*,*) ITER,SADD,RHO(1,i),RHO1(1,i),RD(i)
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
	
      cmsum=(rocm(1)+rocm(m))/2.0
	do i=2,nz-1
	cmsum=cmsum+rocm(i)
	enddo
	av_cm=cmsum/float(nz-1)

! output results
!	OPEN(2,FILE='blk.ini')
	OPEN(3,FILE='result.dat')
	open(5,file='cm.dat')
	WRITE(3,*)'ETA_AV=',sump(1)+sump(2),SUMP(1),SUMP(2)
	
      id=diam/dz

	DO 60 I=1,2*(NZ-1)-1
	IF (I.LT.id.OR.MOD(I,5).EQ.0) THEN
	Y=DZ*FLOAT(I-1)/diam		! Y is in the units of DIAMETER!
!	RHO(1,I)=RHO(1,I)/ROA	! REDUCED density
!	RHO(2,I)=RHO(2,I)/ROB
!	WRITE(2,41) Y,RHO(1,I)/roa,RHO(2,I)/rob
	
      if (i.le.nz) then
	write(5,62) y, rocm(i)/av_cm

      IF (MOD(I,1).EQ.0) WRITE(3,61)Y,RHO(1,I)/roa,RHO(2,I)/rob,
     &                    (rho(1,I)+rho(2,I))/(roa+rob)
	else
	write(3,61) Y,rho(1,2*(nz-1)-I)/roa, rho(2,2*(nz-1)-i)/rob,
     &            (rho(1,2*(nz-1)-I)+rho(2,2*(nz-1)-I))/(roa+rob)
	endif

        
	ENDIF
61	FORMAT(4X,F8.4,4X,F10.7,4X,F10.7,4x,f10.7)
62    format(4x,f8.4,4x,f10.7)
60	CONTINUE
	CLOSE(3)
!	CLOSE(2)
 	STOP
	END
	  
	 subroutine simulation(m,numb,confg)
      parameter (pi=3.1415926535)
	real confg(numb,m,3)
		
  !    open(unit=6,file='con')

  !      pi=4.0d0*datan(1.0d0)
 
	  write(*,*) '...MC start for generating configuration....'
	   
         eslg_star=5.0
	   irigid=0

!	   irigid and eslg_star are the same, but with the different types (integral and real)
!	   cct=acos(1/eslg_star-1.)
!	 write(*,*) 'e*=',eslg_star, 'cct=',cct
!	   cta0=0.0
!	   inumb=0

	   do 70 ik=1,numb
	     	  
	   confg(ik,1,1)=0.0
	   confg(ik,1,2)=0.0
	   confg(ik,1,3)=0.0
	
	    Do 80 i=2,m

	 	if(i-2.le.0.or.irigid.eq.0) then 
            
		   sita=2.0 * pi *ran2(idum)
	       fai=pi*ran2(idum)
	       xk=2.0*ran2(idum)-1.0 		        
	    
		   confg(ik,i,1)=confg(ik,i-1,1)+cos(sita)*sqrt(1.-xk*xk)
	       confg(ik,i,2)=confg(ik,i-1,2)+sin(sita)*sqrt(1.-xk*xk)
	       confg(ik,i,3)=confg(ik,i-1,3)+xk
	       
	  	else 
  
!		  randm1=ran2(idum)
!	      randm2=ran2(idum)
!	      cta0=pi*randm1/2.0
!121	      engr1=exp(-eslg_star*(1.0+cos(pi-cta0)))
		  	          
!		  if(engr1.lt.randm2) then
	        
!			if(randm2.gt.0.95) then
!	         cta0=0
!	         goto 122
!	        endif 

!              cta0=cta0*ran2(idum)
!	        goto 121
!	      endif

	      cta0=pi/3.0

122		  if (cta0.gt.pi/2) then
		    isign=-1
			cta0=pi-cta0
		  else
		    isign=1
		  endif			  	 

	 	  xvi=confg(ik,i-1,1)-confg(ik,i-2,1)
	      yvi=confg(ik,i-1,2)-confg(ik,i-2,2)
	      zvi=confg(ik,i-1,3)-confg(ik,i-2,3)
		  
88	      idice=int(3.0*ran2(idum))+1
	      if(idice.eq.1) then
	        goto 90
	      elseif(idice.eq.2) then
	        goto 95
	      endif

89	       if((1.0-xvi**2).lt.sin(cta0/2)**2) goto 90

	       coscit=1.0-2.0*sin(cta0/2)**2/(1.0-xvi*xvi)
	       sincit=sqrt(1.0-coscit**2)

	       xplus=xvi
	       yplus=yvi*coscit+zvi*sincit
	       zplus=-yvi*sincit+zvi*coscit

	       goto 100

90          if((1.0-yvi**2).lt.sin(cta0/2)**2) goto 95
		   coscit=1.0-2.0*sin(cta0/2)**2/(1.0-yvi*yvi)
	       sincit=sqrt(1.0-coscit**2)

             xplus=coscit*xvi-sincit*zvi
             yplus=yvi
		   zplus=sincit*xvi+coscit*zvi
		   	 
		 goto 100


95		  if((1.0-zvi**2).lt.sin(cta0/2)**2) goto 89

	      coscit=1.0-2.0*sin(cta0/2)**2/(1.0-zvi*zvi)
	      sincit=sqrt(1.0-coscit**2)
			  
             xplus=xvi*coscit+yvi*sincit
             yplus=yvi*coscit-xvi*sincit
		     zplus=zvi
 		  
100    	  if(isign.eq.-1) then
             xplus=xplus-confg(ik,i-1,1)
	       yplus=yplus-confg(ik,i-1,2)
	       zplus=zplus-confg(ik,i-1,3)
		  endif
		    
            confg(ik,i,1)=confg(ik,i-1,1)+xplus
            confg(ik,i,2)=confg(ik,i-1,2)+yplus
	      confg(ik,i,3)=confg(ik,i-1,3)+zplus
	    endif

	
80      continue
70      continue

		igg=0      
        	if(igg.eq.0) goto 99
!	    Do 610 jj=1,m
	    do 610 kk=1,100
	      dis1=confg(kk,1,1)-confg(kk,m,1)
              dis2=confg(kk,1,2)-confg(kk,m,2)
              dis3=confg(kk,1,3)-confg(kk,m,3)   
              dr1=sqrt(dis1**2+dis2**2+dis3**2)
             
               dd2=confg(kk,2,1)-confg(kk,m,1)
              dd22=confg(kk,2,2)-confg(kk,m,2)
              dd23=confg(kk,2,3)-confg(kk,m,3)
              ddr1=sqrt(dd2**2+dd22**2+dd23**2)
              
              dd11=confg(kk,1,1)-confg(kk,2,1)
              dd12=confg(kk,1,2)-confg(kk,2,2)
             dd13=confg(kk,1,3)-confg(kk,2,3)
             ddr2=sqrt(dd11**2+dd12**2+dd13**2) 
   
          write(*,*) dr1,ddr1,ddr2 
610	    continue   
!999     format(3f12.6)

99	 call MAK_RAS_SAD(confg,m)
	 write(*,*)'..'
	 write(*,*) '....Single chain simulation done......'
	 return
	 end
	 

	 subroutine aver(m,ma,eta,nz,dz,ros,rho,rocm,confg,diam)
	 parameter(numb=10000)
	 real confg(numb,m,3),ros(m,2*nz-1),rocm(2*nz-1),diam,sum1(m)
!	 real ang_num(2*nz-1),rend(2*nz-1)
	 real ff(M,2*nz-1),rho(nz),zz(m),cmsum(2*nz-1),cz(m)
	 double precision rosum(m,2*nz-1)
	 integer iz(m),icz(m)
	 
	 
	 Do 10 i=1,nz
	  z1=dz*float(i-1)
	  call sfm1d(m,ma,eta,dz,nz,z1,sum1,rho,diam)
	  
	  do j=1,m
	  	  ff(j,i)=exp(-sum1(j))
	  enddo
!	  write(*,*) i, '  ', ff(1,i),ff(2,i),ff(m,i)
10     continue

       do i=nz+1,2*nz-1
	   do j=1,m
	     ff(j,i)=ff(j,2*nz-i)
	   enddo
	 enddo    

	 do ii=1,2*(nz-1)
	   do jj=1,m
	    rosum(jj,ii)=0.0
	   enddo
           cmsum(ii)=0.0
!           rend(ii)=0.0
!           ang_num(ii)=0.0
	 enddo
     	 	 
	 do 200 i=1,2*(nz-1)-1
!	   write(*,*) i,'....test...'
	    zrz=float(i)*dz
          icount=0
	 do 100 j=1,numb
	  	  cm=0.0
	   do 50 im=1,m
	      
	      do kj=1,m
	       if(im.eq.1) then
		   cm=cm+confg(j,kj,3)/3.0
	       endif
	  	   
		   zz(kj)=(confg(j,kj,3)-confg(j,im,3))*diam+zrz
	       iz(kj)=int(zz(kj)/dz)+1
	   
	       if(kj.ne.im) then
		     if(iz(kj).le.0.or.iz(kj).ge.2*nz-1) goto 50   
	  	   endif
	      enddo
	    
		  fsum=1.0
	      do kk=1,m
	        fsum=fsum*ff(kk,iz(kk))
	      enddo
	    rosum(im,i)=rosum(im,i)+fsum
50	 continue
	  
	    do kj=1,m
	    cz(kj)=(confg(j,kj,3)-cm)*2.0+zrz
	    icz(kj)=int(cz(kj)/dz)+1

	    if(icz(kj).le.0.or.icz(kj).ge.2*nz-1) goto 100
	    enddo

	    cfsum=1.0
	    do kk=1,m
	    cfsum=cfsum*ff(kk,icz(kk))
	    enddo
	    cmsum(i)=cmsum(i)+cfsum
!            icount=icount+1
!            cmx=diam*(confg(j,m,1)-confg(j,1,1))
!            cmy=diam*(confg(j,m,2)-confg(j,1,2))
!            cmz=diam*(confg(j,m,3)-confg(j,1,3))
!            rcm2=cmx*cmx+cmy*cmy+cmz*cmz
!            rend(i)=rend(i)+sqrt(rcm2)
!           ang=abs(cmz)/sqrt(rcm2)
!            ang_num(i)=ang_num(i)+ang 

100    continue
        
	   do km=1,m
	     ros(km,i)=rosum(km,i)/float(numb)
	   enddo

	     rocm(i)=cmsum(i)/float(numb)
!           rend(i)=rend(i)/float(icount)
!           ang_num(i)=ang_num(i)/float(icount)
200     continue
!	 write(*,*) 'average over'
	return
	end
	 
	SUBROUTINE FFLR(M,MA,ETA,NZ,DZ,FL,FR,FF,RHO,diam)
!	calculate the recurrence functions and 
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

!	write(*,*)i,'  ff(j,i)=',ff(j,i)
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
	REAL DZ,ETA,AMU,Z,CPHSC,diam,eps(2),tmp
	REAL RHO(2,NZ),SUM1(M),CPVDW(M),PHI(M)
	common eps,tmp
!	 write(*,*) 'start test.....',diam
!  calculate the reduced chemical potential PER SEGMENT in the bulk phase, AMU
	CALL CPS(M,MA,ETA,AMU,diam)
!	 write(*,*) 'test1-------'
!  calculate external potential
	CALL PHI_E(M,MA,nz,dz,Z,PHI,diam)
!	 write(*,*) 'test2--------'
!  calculate the excess chemical potential due to chain connectivity and hard-sphere interactions 
	CALL CPEHSC(M,DZ,NZ,Z,RHO,CPHSC,diam)
!	 write(*,*) 'test3------'
! calculate the excess chemical potential due to van der Waals attractions
	CALL CPEVDW(M,MA,Z,DZ,NZ,RHO,CPVDW,diam)
! overall density field
!	DO J=1,MA
	SUM1(1)=CPHSC+CPVDW(1)+PHI(1)-AMU
	sum1(m)=cphsc+cpvdw(m)+phi(m)-amu
!	ENDDO
	DO J=2,M-1
	SUM1(J)=CPHSC+CPVDW(J)+PHI(J)-AMU
!	sum1(j)=cphsc-amu
	ENDDO
!      write(*,*) 'sum1(1)=',sum1(1),sum1(2)

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


	SUBROUTINE CPEVDW(M,MA,Z,DZ,NZ,RHO,CPVDW,diam)
	IMPLICIT NONE
	REAL Z,DZ,PI,FRHO,diam
	INTEGER	M,MA,I,J,K,NZ,iz
	REAL Evdw(2,2),CPMF(2),RHO(2,NZ),CPVDW(M)
	Real eps(2),tmp, amidd,datt, z1
	PARAMETER (PI=3.141592654)
	COMMON EPS,tmp
	EXTERNAL FRHO

	Do i=1,2
	do j=1,2
	  Evdw(i,j)=sqrt(eps(i)*eps(j))
	enddo
	enddo
     	iz=int(diam/dz/2.0)
C	write(*,*) 'test9---',tmp,z,iz
	DO 20 I=1,2

	CPMF(I)=0.0
	
	DO 10 J=1,2
	datt=0.0
	Do k=iz,2*(nz-1)-iz
	 z1=float(k-1)*dz
	 amidd=0.1*diam**12/max(diam,abs(z-z1))**10
     &       -0.25*diam**6/max(diam,abs(z-z1))**4
	 datt=datt+Frho(J,k,nz,rho,diam)*amidd
	 enddo

	  CPMF(i)=cpmf(i)+datt*dz*8.0*pi*evdw(i,j)/tmp 

10	  CONTINUE
20	  CONTINUE

!	  DO I=1,MA
        
	  if (ma.eq.2) then
	  CPVDW(1)=CPMF(1)
	  cpvdw(m)=cpmf(1)
	  endif
!	write(*,*) cpmf(1),cpmf(2)
!	  ENDDO

	  DO I=2,M-1
	  CPVDW(I)=CPMF(2)
	  ENDDO
	  RETURN
	  END

	  FUNCTION FRHO(I,IZ,NZ,RHO,diam)
	  IMPLICIT NONE
	  INTEGER I,IZ,NZ
	  REAL FRHO,RHO(2,NZ),diam
	  IF ((IZ.GT.(2*NZ-1)).OR.(IZ.LT.diam/2.)) THEN
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
	    REAL XB(2),eps(2),tmp, evdw(2,2)
	    PARAMETER (PI=3.141592654,PI4=4.0*PI)
    	    COMMON EPS,tmp

		tmp=298.
	    do i=1,2
	    do j=1,2
	    evdw(i,j)=sqrt(eps(i)*eps(j))
	    enddo
	    enddo

	    XB(1)=FLOAT(MA)/FLOAT(M) ! mole fractions of A,B segments
	    XB(2)=FLOAT(M-MA)/FLOAT(M)

	    ROSS=ETA*6./PI/diam**3

!   	ideal-gas term
    	AMR=ALOG(ROSS/FLOAT(M))/FLOAT(M)

!   	hard-sphere term from the Carnahan-Starling equation state
    	AMHS=ETA*(8.-9.*ETA+3.*ETA**2)/(1.-ETA)**3

!   	chain connectivity
	    GHS=(1.-0.5*ETA)/(1.-ETA)**3
	    DGHS=3.*ETA/(1.-ETA)-ETA/(2.-ETA)
    	    AMC=(ALOG(GHS)+DGHS)*FLOAT(1-M)/FLOAT(M)

!    	van der Waals mean-field term
!    	VDWA: van der Waals energy parameter
	    VDWA=0.
	    DO I=1,2
	    DO J=1,2
       	VDWA=VDWA + XB(I)*XB(J)* 
     &      (-32.*PI/9.)*Evdw(I,J)*ross*diam**3/tmp
	 
	    ENDDO
	    ENDDO

	    AMVDW=ROSS*VDWA
	
!   	    overall
    	    AMU=AMR+AMHS+AMC+AMVDW
    	    RETURN
    	    END

	    SUBROUTINE PHI_E(M,MA,nz,dz,Z,PHI,diam)
	    IMPLICIT NONE
!	    parameter (rowall=114.0,sigw=3.4,eswall=28.,delt=3.35)
!   	    parameter (pi=3.1415926)
	
         INTEGER M,MA,I,nz
    	   REAL Z,PHI(M),z2,dz,Atemp1,Atemp2,sigfw
    	   REAL eu12(2),eps(2),tmp,eswall,delt

	    Real diam,rowall,pi,sigw
	    common eps,tmp

!	write(*,*) 'test9------',z,diam/2.0
	
    	    rowall=0.1140
	    sigw=3.4
	    eswall=28.0
	    delt=3.35
	    pi=3.1415926
	    tmp=298.0
	    eps(1)=114.0
	    eps(2)=47.0
        	
         if (z.le.diam/4.0) then
	       do I=1,M
	        phi(i)=1.e30
	       enddo
	   else 
	   sigfw=0.5*(diam+sigw)
	   z2=2.*(nz-1)*dz-z
	   eu12(1)=sqrt(eps(1)*eswall)
	   eu12(2)=sqrt(eps(2)*eswall)

	   eu12(1)=eu12(1)*2.0*pi*rowall*sigfw**2*delt/tmp
	   eu12(2)=eu12(2)*2.0*pi*rowall*sigfw**2*delt/tmp
	 
	   If(Ma.eq.2) then
	    ATemp1=0.4*(sigfw/z)**10-(sigfw/z)**4
     &       -sigfw**4/(3.*delt*(0.61*delt+z)**3)
	    ATemp2=0.4*(sigfw/z2)**10-(sigfw/z2)**4
     &       -sigfw**4/(3.*delt*(0.61*delt+z2)**3)
	    Phi(1)=eu12(1)*(Atemp1+Atemp2)
	    Phi(m)=eu12(1)*(Atemp1+Atemp2) 	  
    	   endif 	
		
    	   DO I=2,M-1
	   ATemp1=0.4*(sigfw/z)**10-(sigfw/z)**4
     &       -sigfw**4/(3.*delt*(0.61*delt+z)**3)
	   ATemp2=0.4*(sigfw/z2)**10-(sigfw/z2)**4
     &       -sigfw**4/(3.*delt*(0.61*delt+z2)**3)
	   Phi(i)=eu12(2)*(Atemp1+Atemp2)
    	   ENDDO
	  endif

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



	 SUBROUTINE  MAK_RAS_SAD(confg,m)
	  
	 implicit real*8 (a-h, o-z)
	        	   
       parameter(numb=10000)
       
	 real    confg(numb,m,3)
	          
	 REAL        SIGMA0
       INTEGER     I
	 REAL        R1, R2

       CHARACTER   PFILE*30
       CHARACTER   C1*4, C2*2, C3*4

        DATA        SIGMA0/5.0/
        DATA        C1/'ATOM'/
        DATA        C2/'C '/
        DATA        C3/'1CRN'/
        DATA        R1/1.00/
        DATA        R2/10.00/

	    PFILE = 'sad.pdb'
	  
	    OPEN ( UNIT = 36, FILE = PFILE, STATUS = 'UNKNOWN')

	 do jj=1,10
       DO 530 I = 1, m
		 	
		 rrx=confg(jj,i,1) * sigma0
		 rry=confg(jj,i,2) * sigma0	        
		 rrz=confg(jj,i,3) * sigma0

	  WRITE(36,10005) C1, I, C2, rrx,rry,rrz, R1, R2, C3, I

530        CONTINUE
	 enddo

         CLOSE(36)

10005   FORMAT(A4, 4X, I3, 2X, A2, 15X, 3F8.3, 2F6.2, 4X, A4, I5)

        RETURN
        END

