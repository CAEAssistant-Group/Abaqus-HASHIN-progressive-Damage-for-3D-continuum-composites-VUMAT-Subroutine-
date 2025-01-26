C                                                 ######################################################################
C                                                 #################      CAE Assistant Company          ################
C                                                 ##############         www CAEassistant com              #############
C                                                 ###########   Copy right © 2021 by CAE Assistant Company    ##########
C                                                 ###################################################################### 
C                                                                         CAE Assisitant Services: 
C                                                 Toturial Packages,Consultancy,Articles,Q&A,Video Gallery,Online Course
C                                                 ######################################################################
C                                                                      Need help with your project? 
C                                                   You can get initial free consultation from (Support@CAEassistant com)
C                                                 ######################################################################  
C                                                                      CAE Assistant  © All Rights Reserved
      
      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
C  LOCAL ARRAY VARIABLES
C
      real  C(6,6), W(NBLOCK,9), M(6,6), CW(6,6), F(6,6,6), 
     *      G1(6,6), G(6,6), E_RATE(6), Q(6,6), 
     *      EPSO(NBLOCK,6),wnew(nblock,9),
     *      EPSN(NBLOCK,6), EFT(NBLOCK,7),STRESSNEWE(NBLOCK,6)
         
     
C  
C  LOCAL REAL AND INTEGER VARIABLES
C   
      real D, E1, E2, E3, NU12, NU21, NU23, NU32, NU31,NU13,
     *     GS12, GS23, GS31, X11T, X22T, X33T, X1M, X2M, X3M,
     *     S12, S23, S13, S12M, S23M, S31M, INEL, ED, X11P, 
     *     X22P, X33P, SMT, SMC, GFT, GFC, GMT, GMC, E11FT,
     *     E11FC, E22MT, E22MC, EFC, EMT, EMC,GSHEAR12,E12SH,
     *     E31SH,E23SH, GSHEAR13,GSHEAR23  

	INTEGER I, K, J, Z, L, P, R, N, V

	
	E1 = PROPS(1)
	E2 = PROPS(2)
	E3 = PROPS(3)
	NU12 = PROPS(4)	
	NU23 = PROPS(5)	
	NU13 = PROPS(6)
	GS12 = PROPS(7)
	GS23 = PROPS(8)
	GS31 = PROPS(9)
C------------TENSILE STRENGTH---------------------------------------
	X11T = PROPS(10)
	X22T = PROPS(11)
	X33T = PROPS(12)
C--------------SHEAR STRENGTH---------------------------------------
	S12 = PROPS(13)
	S23 = PROPS(14)
	S13 = PROPS(15)
C----------COMPRESSION STRENGTH-------------------------------------
	X11P = PROPS(16)
	X22P = PROPS(17)
	X33P = PROPS(18)
C SHEAR STIFFNESS LOSS DUE TO MATRIX FAILURE UNDER LOADING-SMT & SMC-
	SMT = PROPS(19)
	SMC = PROPS(20)
C-----INTRALAMINAR TOUGHNESS PROPERTIES-----------------------------
      GFT = PROPS(21)
      GFC = PROPS(22)
      GMT = PROPS(23)
      GMC = PROPS(24) 
C----------COMPRESSION STRENGTH OF MATRIX ---------------------------
      MCS = PROPS(25)
C----------SHEAR TOUGHNESS PROPERTIES--------------------------------
      GSHEAR12 = PROPS(26)
      GSHEAR13 = PROPS(27)
      GSHEAR23 = PROPS(28)
c-----------------Statev----------------------------------------------
c       statev 1=w1
c       statev 2-w2
c       statev 3=w3
c       statev 4=w4
c       statev 5=w5
c       statev 6=w6
c       statev 7=w7
c       statev 8=w8
c       statev 9=w9
c       statev10=Eps1
c       statev11=eps2
c       statev12=eps3
c       statev13=eps12
c       statev14=eps23
c       statev15=eps13
c       statev16=element deletion
c       statev17=eft1
c       statev18=eft2
c       statev19=eft3
c       statev20=eft4
c       statev21=eft5
c       statev22=eft6
c       statev23=eft7
c       statev24-29=elastic stress
          	
C   INITIALISING VARIABLE FOR ELEMENT DELETION (ED=1->ELEMENTS ACTIVE; ED=0->ELEMENTS INACTIVE)
	ED = 1
C------COMPUTE MINOR POISON RATIO---------------------
      NU21=NU12*E2/E1
      NU31=NU13*E3/E1
      NU32=NU23*E3/E2
C-----COMPUTE FAILURE STRAINS-------------------------
      E11FT=(X11T/E1)
      E11FC=(X11P/E1)
      E33MT=(X33T/E3)
      E33MC=(X33P/E3)
      E22MT=(X22t/E2)
      E22MC=(X22P/E2)
      E12SH=S12/(2.d0*GS12)
      E23SH=S23/(2.d0*GS23)
      E31SH=S13/(2.d0*GS31)
C   CREATE A ZERO 6*6 STIFFNESS MATRIX
         DO I=1,6
            DO J=1,6
                C(I,J)=0.d0
           END DO
         END DO

C
C   UNDAMAGED STIFFNESS MATRIX FOR ORTHOTROPIC LINEAR ELASTIC MATERIAL

	D = (1.d0-NU12*NU21-NU23*NU32-NU31*NU13-2.d0*NU21*NU32*NU13)/(E1*E2*E3)

		C(1,1) = (1.d0-NU23*NU32)/(E2*E3*D)
		C(2,1) = (NU21+NU23*NU31)/(E2*E3*D)
		C(3,1) = (NU31+NU21*NU32)/(E3*E2*D)
		C(4,1) = 0.
		C(5,1) = 0.
		C(6,1) = 0.

		C(1,2) = (NU21 + NU23*NU31)/(E2*E3*D)
		C(2,2) = (1.d0- NU13*NU31)/(E1*E3*D)
		C(3,2) = (NU32 + NU12*NU31)/(E1*E3*D)
		C(4,2) = 0.
		C(5,2) = 0.
		C(6,2) = 0.

		C(1,3) = (NU31+NU21*NU32)/(E3*E2*D)
		C(2,3) = (NU32+NU12*NU31)/(E1*E3*D)
		C(3,3) = (1.d0-NU12*NU21)/(E1*E2*D)
		C(4,3) = 0.
		C(5,3) = 0.
		C(6,3) = 0.

		C(1,4) = 0.
		C(2,4) = 0.
		C(3,4) = 0.
		C(4,4) = GS12
		C(5,4) = 0.
		C(6,4) = 0.

		C(1,5) = 0.
		C(2,5) = 0.
		C(3,5) = 0.
		C(4,5) = 0.
		C(5,5) = GS23
		C(6,5) = 0.

		C(1,6) = 0.
		C(2,6) = 0.
		C(3,6) = 0.
		C(4,6) = 0. 
		C(5,6) = 0.
		C(6,6) = GS31

C
C LOOP OVER BLOCK OF ELEMENTS
C
	DO 1 K = 1, NBLOCK
	
        
            DO I=1,9
               W(K,I) = STATEOLD(K,I)               
            END DO
			DO I = 1,6
				EPSO(K,I) = STATEOLD(K,I+9)
			END DO
			DO I = 1,3
				EPSN(K,I) =STRAININC(K,I) + EPSO(K,I)
			END DO
			DO I = 4,6
				EPSN(K,I) =STRAININC(K,I) + EPSO(K,I)
			END DO
c			print*,"epsn(k,6)=",epsn(k,6)
			DO I=1,7
			    EFT(K,I)=STATEOLD(K,I+16)
			END DO
			DO I=1,6
			    STRESSNEWE(K,I)=STATEOLD(K,I+23)
              END DO
              
C			DO I=1,4
C			    EFT(K,I) = STATEOLD(K,I+12)	
C			END DO	

C     CREATE DAMAGE STIFFNESS MATRIX (CW)
          DO I=1,6
           DO J=1,6
                CW(I,J)=0.
            END DO
          END DO
      
              CW(1,1)=(1.-W(K,5))* C(1,1)
              CW(2,2)=!Hidden
              CW(3,3)=(1.-W(K,5))*(1.-W(K,6))*C(3,3)
              CW(1,2)=(1.-W(K,5))*(1.-W(K,6))*C(1,2)
			!Hidden
              CW(1,3)=(1.-W(K,5))*(1.-W(K,6))*C(1,3)
              CW(3,1)=(1.-W(K,5))*(1.-W(K,6))*C(3,1)
              CW(4,4)=!Hidden
     *            C(4,4)
              CW(5,5)=!Hidden
     *            
              CW(6,6)=!Hidden
C      COMPUTE STRESSNEW = !Hidden
c============================================================================
		CALL MATVEC(STRESSNEW(K,:),CW,EPSN(K,:),6,6)
         
c============================================================================		 
C---------------MATRIX CRUSHING MODELLING-----------------------------------
C            IF (STRESSNEW(K,1) .LT. 0)THEN
C                IF (STRESSNEW(K,1) .EQ. 0)THEN
C                STRESSNEW(K,1) = -MCS
C                END IF
C            END IF
C---------------------------------------------------------------------------
     
      IF (STRESSNEW(K,1) .GT. 0.) THEN
            
        IF (W(K,1) .GT. 0.) THEN
                CALL ELASTRE(STRESSNEWE(K,:),C,EPSN(K,:),6,6)
			!Hidden                  

C----------------------------------------------------------------------------
C----------------------------------------------------------------------------
C----------------------------------------------------------------------------
C----------------------------------------------------------------------------
               
        ELSE IF (W(K,1) .EQ. 0.)THEN                 
            EFT(K,1) = (STRESSNEW(K,1)/X11T)**2
     1                      +(STRESSNEW(K,4)/S12)**2.
				!Hidden(about 1000 lines)
               end if

!C----------------------------------------------------------------------------

!C----------------------------------------------------------------------------
                    IF ((STRESSNEW(K,2) + STRESSNEW(K,3)) .GT. 0.)THEN
                        IF (W(K,3) .GT. 0.)THEN
                  CALL ELASTRE (STRESSNEWE(K,:),C,EPSN(K,:),6,6)
                  
                        EFT(K,3) = ((STRESSNEWE(K,2) + 
     *                   STRESSNEWE(K,3))/X22T)**2+
     *                  ((STRESSNEWE(K,5)**2 - STRESSNEWE(K,2)*
     *                STRESSNEWE(K,3))/S23**2)+(STRESSNEWE(K,4)/S12)**2+
     *                  (STRESSNEWE(K,6)/S13)**2
                        if(eft(k,3) .gt. 1)then
                         W(K,3)= 1.d0-(1.d0/EFT(K,3))*exp(-C(2,2)*
     *                       E22MT*E22MT*
     *                       (EFT(K,3)-1.d0)*CHARLENGTH(K)/GMT)
			!Hidden
                         end if                   
                            
     			!Hidden(about 500 lines)

      
           
      end if
c-----------------mode: fiberkinking-------------------------------  
	  IF(STRESSNEW(K, 1) .LT. 0.)THEN
c			!Hidden(about 300 lines)

	  END IF
C----------------SHEAR MODE---------Maximum stress failure Criteria----------------------------------- 
c====================================================================
c===================shear mode 12==================================== 
c        print*,"e12=",epsn(k,4)
	  IF (ABS(STRESSNEW (K,4)) .GT. 0.)THEN
			!Hidden(about 50 lines)

        END IF
c========================shear mode 23===================================        
        IF (ABS(STRESSNEW (K,5)) .GT. 0.)THEN
		!Hidden(about 100 lines)
        END IF
c   =======================shear mode 13===================================
c        print*,"STRESSNEW (K,6)=",STRESSNEW (K,6)
			!Hidden(about 50 lines)
c   =======================shear mode 13===================================
	       

     
			DO I =1,9
				wnew(K,I) = W(K,I)
			END DO

			DO I = 1,6
				STATENEW(K,I+9) = EPSN(K,I)
			END DO
			DO I = 1,7
			    STATENEW(K,I+16) = EFT(K,I)
			END DO
			DO I=1,6
			    STATENEW(K,I+23)=STRESSNEWE(K,I)
              END DO
             
C   STATE VARIABLE CONTROLLING ELEMENT DELETION			    
			    

			
			do i=1,9
			   if (w(k,i) .lt. 0.)then
					w(k,i)=0.
			   end if
			end do
			
			DO I = 1,9
			!Hidden
			END DO
				
			!Hidden
C			IF (W(K,3) .GE. 0.99)THEN
C			   !Hidden
C			END IF
			
			STATENEW(K,16) = ED
c		END IF

1	CONTINUE
	RETURN
	END SUBROUTINE VUMAT

C ---------------------------------------------------------
      SUBROUTINE MATVEC(V,A,U,NA,NB)
C     PERFORMS MATRIXXVECTOR PRODUCT V=A*U
C     A(NA,NB)
C     U(NB)
      IMPLICIT NONE
      INTEGER I,J,NA,NB
      REAL V(NA),A(NA,NB),U(NB)

      DO I=1,NA
		V(I)=0
        DO J=1,NB
          V(I)=V(I)+A(I,J)*U(J)
        ENDDO
      ENDDO
      END SUBROUTINE MATVEC
      
C---------------------------------------------------------

      SUBROUTINE ELASTRE(V,A,U,NA,NB)
      IMPLICIT NONE
      INTEGER I,J,NA,NB
      REAL V(NA),A(NA,NB),U(NB)

      DO I=1,NA
		V(I)=0
        DO J=1,NB
          V(I)=V(I)+A(I,J)*U(J)
        ENDDO
      ENDDO
      END SUBROUTINE ELASTRE