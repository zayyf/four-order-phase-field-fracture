! ======================================================================
! User SUBROUTINE UEL for Abaqus PHASE-Field implementation:
! All rights of reproduction or distribution in any form are reserved.
! ======================================================================
! N_ELEM : Numder of elements.
! Material properties to be given through the input file, where
! PROPS(1) = Young's Modulus
! PROPS(2) = Poisson's Ratio
! PROPS(3) = Failure stress
! PROPS(4) = Length Scale Parameter, Characteristic Crack Width
! PROPS(5) = Strain Energy Release Rate
! PROPS(6) = Softening Parameter 
! PROPS(7) = Monolithic (kms=0) and staggered (kms=1)
! ======================================================================
! Used variables 
!
! NSVARS/3 -  solution dependent variables for the phase field element 
!             and displacement field emement
!            (1/phase, 1/gradient, 1/engery)  
!     sdv(1)=phase
! ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
!     ==================================================================
      INCLUDE 'ABA_PARAM.INC'
!     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 SIX=6.D0,TOLER=1.0D-6,RP25=0.25D0,HALF=0.5D0,RP7=0.7D0,
     2 PI=3.1415926535897932384626433832D0,NSTV=4,N_ELEM=5527)
!     ==================================================================
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
!
       REAL*8 AINTW(8),XII(8,MCRD),XI(MCRD),dNdxi(NNODE,MCRD),  
     1 dNdx(NNODE,MCRD),VJACOB(MCRD,MCRD),VJABOBINV(MCRD,MCRD),  
     2 BB(6,NDOFEL),EPS(6),STRESS(6),CMAT(6,6),DDSDDE(6,6), 
     3 SDV(NSVARS),AN(NNODE),BP(MCRD,NNODE),DP(MCRD),VNI(MCRD,NDOFEL),
     4 ULOC(MCRD),B_B(NNODE,NNODE),EIGV(3),ALPHAI(3),VECTI(3),
     5 DCMAT(6,6),EPSC(6),BG(MCRD,NNODE),DG(MCRD)
	
!
       REAL*8 DTM,EG,EG2,EG3,EBULK,EBULK3,ELAM,EMOD,ENU,FTPAR
      
       REAL*8 SAVG,SDIF,SDEV,SMAX,SMIN,EAVG,EDIF,EDEV,EMAX,EMIN
	 
       REAL*8 ENG0,ENG,ENGN,PHASE,PARK,CLPAR,GCPAR,PHASE0,DPHASE,HIST,
     1 HISTN,ALPHAE,THCK,VM_STRESS,GRAD,ENGM,ENGP,MK,PARMS,HISTM
	 
       REAL*8 ALPHA,DALPHA,DDALPHA,OMEGA,DOMEGA,DDOMEGA,P,A1,A2,A3,C0,
     1 PHASE_SOURCE,DPHASE_SOURCE,ELAMEG,ELAMEL,ENGKG,ENGEG,ENGDG
!
       INTEGER I,J,K,L,K1,K2,K3,K4,NM,NN,NP     
      
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
!
!     History variables
        ENGKG=ZERO
        ENGEG=ZERO
        ENGDG=ZERO
        NALL=NSTVTO+NSTVTT
!     ==================================================================
!     Material parameters
!     ==================================================================
       EMOD =PROPS(1)
       ENU  =PROPS(2)
       FTPAR=PROPS(3)
       CLPAR=PROPS(4)
       GCPAR=PROPS(5)
       PARK =PROPS(6)
       PARMS=PROPS(7)	   
       THCK =ONE	   
       C0=96.D0/5.D0
       A1=15.D0*EMOD*GCPAR/(16.D0*CLPAR*FTPAR**TWO)
       MK=SQRT(TWO*GCPAR*CLPAR**THREE/C0)	 !ONE 
       ELAMEG=EMOD/(TWO*(ONE+ENU))
       ELAMEL=ELAMEG*TWO*ENU/(ONE-TWO*ENU) 	   
!     ==================================================================
!     Initial preparetions
!     ==================================================================
       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO
!     ==================================================================
!     Local coordinates and weights
!     ==================================================================   	   
       IF (JTYPE.EQ.ONE) THEN
! U1: 3-node triangular elements    
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = 1.D0
        AINTW(1) = HALF
!
       ELSEIF (JTYPE.EQ.TWO) THEN
! U2: 4-node square elements	   
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = 4.D0
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
!		
       ELSEIF (JTYPE.EQ.THREE) THEN
! U3: 4-node linear tetrahedron elements 	
        XII(1,1) = ONE/FOUR
        XII(1,2) = ONE/FOUR
        XII(1,3) = ONE/FOUR
        INNODE=1.D0
        AINTW(1) = ONE/SIX
!
       ELSEIF (JTYPE.EQ.FOUR) THEN
! U4: 8-node linear brick elements 
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(1,3) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(2,3) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(3,3) = -ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        XII(4,3) = -ONE/THREE**HALF
        XII(5,1) = -ONE/THREE**HALF
        XII(5,2) = -ONE/THREE**HALF
        XII(5,3) = ONE/THREE**HALF
        XII(6,1) = ONE/THREE**HALF
        XII(6,2) = -ONE/THREE**HALF
        XII(6,3) = ONE/THREE**HALF
        XII(7,1) = ONE/THREE**HALF
        XII(7,2) = ONE/THREE**HALF
        XII(7,3) = ONE/THREE**HALF
        XII(8,1) = -ONE/THREE**HALF
        XII(8,2) = ONE/THREE**HALF
        XII(8,3) = ONE/THREE**HALF
        INNODE = 8.D0
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ENDIF	   	   
!
!     ==================================================================
!     Calculating properties at each integration point
!     ==================================================================   
       DO INPT=1,INNODE
!     Initializing solution dependent variables
        DO I=1,NSTV
           SDV(I)=ZERO
        END DO
        DO I=1,NSTV
          SDV(I)=SVARS(NSTV*(INPT-1)+I)
        END DO
!
!     Local coordinates of the integration point
        DO I=1,MCRD
         XI(I) = XII(INPT,I)
        END DO
! 
!     Shape functions and local derivatives
        IF (JTYPE.EQ.ONE) THEN
         CALL SHAPEFUNTRI(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.TWO) THEN
         CALL SHAPEFUNQUAD(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.THREE) THEN
         CALL SHAPEFUNTET(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.FOUR) THEN
         CALL SHAPEFUNBRICK(AN,dNdxi,XI)
        ENDIF
! 
!     Shape functions
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 		
         IY=ZERO
         DO I = 1,NNODE
          IX=IY+1
          IY=IX+1
          VNI(1,IX)=AN(I)
          VNI(1,IY)=ZERO
          VNI(2,IX)=ZERO
          VNI(2,IY)=AN(I)
         END DO
        ELSE
! --------- 3D Elements ------------  
         IZ=ZERO
         DO I = 1,NNODE
          IX=IZ+ONE
          IY=IX+ONE
          IZ=IY+ONE
          VNI(1,IX)=AN(I)
          VNI(2,IX)=ZERO
          VNI(3,IX)=ZERO
          VNI(1,IY)=ZERO
          VNI(2,IY)=AN(I)
          VNI(3,IY)=ZERO
          VNI(1,IZ)=ZERO
          VNI(2,IZ)=ZERO
          VNI(3,IZ)=AN(I)
         END DO
        ENDIF
!		
!     Jacobian
        DO I = 1,MCRD
         DO J = 1,MCRD
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO
!        
        DTM = ZERO
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------        
         DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
        ELSE
! --------- 3D Elements ------------        
        DTM = VJACOB(1,1)*VJACOB(2,2)*VJACOB(3,3)+VJACOB(1,2)*
     1   VJACOB(2,3)*VJACOB(3,1)+VJACOB(1,3)*VJACOB(2,1)*
     2   VJACOB(3,2)-VJACOB(3,1)*VJACOB(2,2)*VJACOB(1,3)-
     3   VJACOB(3,2)*VJACOB(2,3)*VJACOB(1,1)-VJACOB(3,3)*
     4   VJACOB(2,1)*VJACOB(1,2)        
        ENDIF 
!		
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT	
        ENDIF				

!     Inverse of Jacobian
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------        
         VJABOBINV(1,1)= VJACOB(2,2)/DTM
         VJABOBINV(1,2)=-VJACOB(1,2)/DTM
         VJABOBINV(2,1)=-VJACOB(2,1)/DTM
         VJABOBINV(2,2)= VJACOB(1,1)/DTM
        ELSE
! --------- 3D Elements ------------        
         VJABOBINV(1,1)=(VJACOB(2,2)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,2))/DTM
         VJABOBINV(1,2)=-(VJACOB(1,2)*VJACOB(3,3)-VJACOB(3,2)*
     1    VJACOB(1,3))/DTM
         VJABOBINV(1,3)=(VJACOB(1,2)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,2))/DTM
         VJABOBINV(2,1)=-(VJACOB(2,1)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,2)=(VJACOB(1,1)*VJACOB(3,3)-VJACOB(1,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,3)=-(VJACOB(1,1)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,1))/DTM
         VJABOBINV(3,1)=(VJACOB(2,1)*VJACOB(3,2)-VJACOB(2,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,2)=-(VJACOB(1,1)*VJACOB(3,2)-VJACOB(1,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,3)=(VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*
     1    VJACOB(2,1))/DTM      
        ENDIF 
!        
!     Derivatives of shape functions respect to global ccordinates
        DO K = 1,NNODE
         DO I = 1,MCRD
          dNdx(K,I) = ZERO
          DO J = 1,MCRD
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO
!
!     Calculating B matrix (B=LN)
       IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
        IY=0
        DO INODE=1,NNODE
         IX=IY+1
         IY=IX+1
         BB(1,IX)= dNdx(INODE,1)
         BB(1,IY)= ZERO
         BB(2,IX)= ZERO
         BB(2,IY)= dNdx(INODE,2)
         BB(3,IX)= dNdx(INODE,2)
         BB(3,IY)= dNdx(INODE,1)
        END DO
!
       ELSE
! --------- 3D Elements ------------ 
        IZ=ZERO
        DO INODE=1,NNODE
         IX=IZ+ONE
         IY=IX+ONE
         IZ=IY+ONE
         BB(1,IX)= dNdx(INODE,1)
         BB(2,IX)= ZERO
         BB(3,IX)= ZERO
         BB(4,IX)= dNdx(INODE,2)
         BB(5,IX)= dNdx(INODE,3)
         BB(6,IX)= ZERO
         BB(1,IY)= ZERO
         BB(2,IY)= dNdx(INODE,2)
         BB(3,IY)= ZERO
         BB(4,IY)= dNdx(INODE,1)
         BB(5,IY)= ZERO
         BB(6,IY)= dNdx(INODE,3)
         BB(1,IZ)= ZERO
         BB(2,IZ)= ZERO
         BB(3,IZ)= dNdx(INODE,3)
         BB(4,IZ)= ZERO
         BB(5,IZ)= dNdx(INODE,1)
         BB(6,IZ)= dNdx(INODE,2)
        END DO
       ENDIF	   
!		
       DO INODE=1,NNODE
        DO K1=1,MCRD
         BP(K1,INODE)=dNdx(INODE,K1)
        END DO
       END DO
!
       DO INODE=1,NNODE
        DO K1=1,MCRD
         BG(K1,INODE)=dNdx(INODE,K1)
        END DO
       END DO	   
	   
!
!     ==================================================================
!     Nodal displacements
!     ==================================================================
        DO J=1,MCRD
         ULOC(J)=ZERO
        END DO
        DO J=1,MCRD
         DO I=1,MCRD*NNODE
          ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
         END DO
        END DO  
!        DO J=1,MCRD
!         SDV(J)=ULOC(J)
!        END DO
! 
!     ==================================================================
!     Nodal phase-field
!     ==================================================================		
        PHASE=ZERO
        DO I=1,NNODE
        PHASE=PHASE+AN(I)*U(I+MCRD*NNODE)
        END DO
!         
        PHASE0=SDV(1)
        IF (PHASE.LT.SDV(1)) THEN
         PHASE=SDV(1)
        ENDIF
!
        IF (PHASE.LT.ZERO) THEN
         PHASE=ZERO
        ELSEIF (PHASE.GT.ONE) THEN
         PHASE=ONE
        ENDIF
        SDV(1)=PHASE      
!      
        DPHASE=ZERO 
        DPHASE=PHASE-PHASE0 
        IF (DPHASE.LT.ZERO) THEN
         DPHASE=ZERO
        ENDIF

!       Gradient
        DO I=1,MCRD
         DP(I)=ZERO
        END DO
        DO I=1,MCRD
         DO J=1,NNODE
          DP(I)=DP(I)+BP(I,J)*U(J+MCRD*NNODE)
         END DO
        END DO		
!
!     ==================================================================
!     Nodal potential field--gradient
!     ==================================================================		
        GRAD=ZERO
        DO I=1,NNODE
         GRAD=GRAD+AN(I)*U(I+MCRD*NNODE+NNODE)
        END DO
!
        SDV(2)=GRAD     
		
!       Gradient
        DO I=1,MCRD
         DG(I)=ZERO
        END DO
        DO I=1,MCRD
         DO J=1,NNODE
          DG(I)=DG(I)+BG(I,J)*U(J+MCRD*NNODE+NNODE)
         END DO
        END DO
!		
!     ==================================================================
!     Calculating strain
!     ================================================================== 
       IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
        NP=3		
       ELSE	   
! --------- 3D Elements ------------       
        NP=6
       ENDIF

        DO J=1,NP
         EPS(J)=ZERO
        END DO
        DO I=1,NP
         DO J=1,MCRD*NNODE
           EPS(I)=EPS(I)+BB(I,J)*U(J)    
         END DO
        END DO
!		
        DO K1=1,6
         EPSC(K1)=ZERO
        END DO		
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
         EPSC(1)=EPS(1)
         EPSC(2)=EPS(2)
         EPSC(4)=EPS(3) 		 
        ELSE	   
! --------- 3D Elements ------------       
         EPSC=EPS
        ENDIF		
!        DO J=1,BP
!         SDV(J)=EPS(J)
!        END DO
!     ==================================================================
!     Calculating principal materials stiffness matrix
!     ==================================================================
       IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
        DO I=1,3
         DO J=1,3
          CMAT(I,J)=ZERO
         END DO
        END DO
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
!		
       ELSE	   
! --------- 3D Elements ------------       
        DO I=1,6
         DO J=1,6
          CMAT(I,J)=ZERO
         END DO
        END DO
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
       ENDIF
!
!     ==================================================================
!     Calculating stresses
!     ==================================================================

        CALL GEOMETRICFUNC(ALPHA,DALPHA,DDALPHA,PHASE) ! GEOMETRIC FUNCTION
        CALL ENERGETICFUNC(PROPS,OMEGA,DOMEGA,DDOMEGA,PHASE) ! ENERGETIC FUNCTION
	   
        DO K1=1,NP
         STRESS(K1)=ZERO
        END DO
!
        DO K1=1,NP
         DO K2=1,NP
          STRESS(K2)=STRESS(K2)+CMAT(K2,K1)*EPS(K1)
         END DO
        END DO
		
!      Calculating von Mises Stress		
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
        VM_STRESS=(HALF*((STRESS(1)-STRESS(2))**TWO+STRESS(1)**TWO+ 
     1  STRESS(2)**TWO+SIX*(STRESS(3)**TWO)))**HALF	 
!		
        ELSE	   
! --------- 3D Elements ------------  
        VM_STRESS=(HALF*((STRESS(1)-STRESS(2))**TWO+
     1  (STRESS(2)-STRESS(3))**TWO+(STRESS(3)-STRESS(1))**TWO+ 
     2  SIX*(STRESS(4)**TWO+STRESS(5)**TWO+STRESS(6)**TWO)))**HALF	 
        ENDIF
		
        SDV(4)=VM_STRESS
		
!		
!        DO J=1,NP
!         SDV(J+6)=STRESS(J)*OMEGA
!        END DO
! 
!     ==================================================================
!     Calculating elastic ENERGY
!     ==================================================================

        ENG0=9.D0*GCPAR/(C0*CLPAR*A1)
!        ENG0=HALF*FTPAR**TWO/EMOD

        CALL EIGOWN(EPSC,EIGV,ALPHAE,ALPHAI,VECTI)
!       
        ENGP=(ELAMEL*(ALPHAE*(EIGV(1)+EIGV(2)+EIGV(3)))**TWO)/
     1       TWO+ELAMEG*((EIGV(1)*ALPHAI(1))**TWO+(EIGV(2)*
     2       ALPHAI(2))**TWO+(EIGV(3)*ALPHAI(3))**TWO)
			
        ENGN=(ELAMEL*((ONE-ALPHAE)*(EIGV(1)+EIGV(2)+EIGV(3)))**
     1   TWO)/TWO+ELAMEG*((EIGV(1)*(ONE-ALPHAI(1)))**TWO+(EIGV(2)*
     2   (ONE-ALPHAI(2)))**TWO+(EIGV(3)*(ONE-ALPHAI(3)))**TWO)
			
!        ENG=ZERO
!        DO K2=1,NP
!         ENG=ENG+STRESS(K2)*EPS(K2)*HALF
!        END DO
		
!        WRITE(7,*) 'EIGV',EIGV

        ENGM=max(ENGP,ENG0)
        ENGEG=ENGEG+((ENGP+ENGN)*OMEGA)*DTM*THCK		
!
!     ==================================================================
!     Calculating elastic ENERGY history
!     ==================================================================         
!
        HISTN=SDV(3)
        HIST=max(ENGM,HISTN)
        SDV(3)=HIST
		
        IF (NINT(PARMS)==0) THEN 
         HISTM=HIST                     ! monolithic scheme
        ELSEIF (NINT(PARMS)==1) THEN
         HISTM=max(ENG0,HISTN)          ! staggered scheme
        ENDIF
		
!
!     ==================================================================
!     Calculating fracture energy for history output
!     ==================================================================
!
        ENGD=ZERO
        ENGD=ALPHA/(C0*CLPAR)*DTM*GCPAR
!
        ENGD=ENGD+GRAD*GRAD*CLPAR**THREE*DTM/C0*GCPAR
        ENGDG=ENGDG+ENGD
!
!     ==================================================================
!     Calculating element stiffness matrix(K_uu)
!     ==================================================================
!
        DO K=1,MCRD*NNODE
         DO L=1,MCRD*NNODE
          DO I=1,NP
           DO J=1,NP
            AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1                  BB(J,L)*DTM*OMEGA*THCK
           END DO
          END DO
         END DO
        END DO
!       
!     ==================================================================
!     Internal forces (residual vector) (r_u)
!     ==================================================================

        DO K1=1,MCRD*NNODE
         DO K4=1,NP
           RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*STRESS(K4)*DTM*
     1               OMEGA*THCK
         END DO
        END DO
!  
!     ==================================================================
!     Calculating element stiffness matrix(K_dd)
!     ==================================================================
 
        PHASE_SOURCE=DOMEGA*HISTM+GCPAR/(C0*CLPAR)*DALPHA
        DPHASE_SOURCE=DDOMEGA*HISTM+GCPAR/(C0*CLPAR)*DDALPHA
!		
        NM=MCRD*NNODE
        DO I=1,NNODE
        NM=NM+1
        NN=MCRD*NNODE
         DO J=1,NNODE
           NN=NN+1
           AMATRX(NM,NN)=AMATRX(NM,NN)+AN(I)*AN(J)*AINTW(INPT)*
     1     DTM*(DPHASE_SOURCE)*THCK
         END DO
        END DO      
!
!     ==================================================================
!     Calculating element stiffness matrix(K_dg)
!     ==================================================================
!		
        NM=MCRD*NNODE
        DO I=1,NNODE
        NM=NM+1
        NN=MCRD*NNODE+NNODE
         DO J=1,NNODE
           NN=NN+1
           DO K=1,MCRD		  
            AMATRX(NM,NN)=AMATRX(NM,NN)-TWO*GCPAR*CLPAR**THREE/(C0*MK)*
     1      BP(K,I)*BG(K,J)*AINTW(INPT)*DTM*THCK
           END DO
         END DO
        END DO 
!		
!     ==================================================================
!     Internal forces (residual vector) (r_d)
!     ==================================================================

        NM=MCRD*NNODE
        DO I=1,NNODE
        NM=NM+1
         DO J=1,MCRD
           RHS(NM,1)=RHS(NM,1)+BP(J,I)*DG(J)*TWO*GCPAR*CLPAR**THREE/
     1     (C0*MK)*AINTW(INPT)*DTM*THCK
         END DO
         RHS(NM,1)=RHS(NM,1)-AN(I)*AINTW(INPT)*DTM*(PHASE_SOURCE)*THCK  
        END DO		  
!
!     ==================================================================
!     Calculating element stiffness matrix(K_gd)
!     ==================================================================
!		
        NM=MCRD*NNODE+NNODE
        DO I=1,NNODE
        NM=NM+1
        NN=MCRD*NNODE
         DO J=1,NNODE
           NN=NN+1
           DO K=1,MCRD		  
            AMATRX(NM,NN)=AMATRX(NM,NN)-MK*BG(K,I)*BP(K,J)*AINTW(INPT)*
     1      DTM*THCK
           END DO
         END DO
        END DO
!
!     ==================================================================
!     Calculating element stiffness matrix(K_gg)
!     ==================================================================
!		
        NM=MCRD*NNODE+NNODE
        DO I=1,NNODE
        NM=NM+1
        NN=MCRD*NNODE+NNODE
         DO J=1,NNODE
            NN=NN+1
            AMATRX(NM,NN)=AMATRX(NM,NN)-AN(I)*AN(J)*AINTW(INPT)*
     1      DTM*THCK
         END DO
        END DO		
!		
!     ==================================================================
!     Internal forces (residual vector) (r_g)
!     ==================================================================

        NM=MCRD*NNODE+NNODE
        DO I=1,NNODE
        NM=NM+1
         DO J=1,MCRD
           RHS(NM,1)=RHS(NM,1)+MK*BG(J,I)*DP(J)*AINTW(INPT)*DTM*THCK
         END DO
         RHS(NM,1)=RHS(NM,1)+AN(I)*AINTW(INPT)*DTM*GRAD*THCK
        END DO		
!     
!     ==================================================================
!     Uploading solution dep. variables
!     ==================================================================    
        DO I=1,NSTV
          SVARS(NSTV*(INPT-1)+I)=SDV(I)
          USRVAR(JELEM,I,INPT)=SDV(I)
        END DO
       END DO
!
        ENERGY(1)=ENGKG	 
        ENERGY(2)=ENGEG	  
        ENERGY(7)=ENGDG

      RETURN
      END
!
!----------------------------------------------------------------------
!
      SUBROUTINE ENERGETICFUNC(PROPS,OMEGA,DOMEGA,DDOMEGA,PHASE)
!
! ----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
!      
      DIMENSION PROPS(*)
      REAL*8 EMOD,ENU,FTPAR,CLPAR,GCPAR,PARK
      REAL*8 OMEGA,DOMEGA,DDOMEGA,PHASE
      REAL*8 FAC1,DFAC1,DDFAC1,FAC2,DFAC2,DDFAC2
      REAL*8 P,A1,A2,A3
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 SIX=6.D0,PI=3.1415926535897932384626433832D0)
!	  
!      Material parameters
       EMOD =PROPS(1)
       ENU  =PROPS(2)
       FTPAR=PROPS(3)
       CLPAR=PROPS(4)
       GCPAR=PROPS(5)
       PARK =PROPS(6)
      
!      Softening  parameters
!     
       A1=15.D0*EMOD*GCPAR/(16.D0*CLPAR*FTPAR**TWO)	

       IF (NINT(PARK)==1) THEN ! Linear softening
         P = 2.D0
         A2= -0.4709D0
         A3= 0.5D0
       ELSEIF (NINT(PARK)==2) THEN ! Concrete softening
         P = 2.D0
         A2= 1.3941D0
         A3= 0.D0			 
       ELSE
         WRITE (*,*) '**ERROR: MODEL LAW NO. ', PARK, 
     1                 'DOES NOT EXIST!'
         CALL XIT
       ENDIF 		  
	          
       FAC1    = (ONE-PHASE)**P
       DFAC1   = -P*(ONE-PHASE)**(P-ONE)
       DDFAC1  = P*(P-ONE)*(ONE-PHASE)**(P-TWO)
	  
       FAC2    = FAC1+ A1*PHASE+ A1*A2*PHASE**TWO+ A1*A2*A3*PHASE**THREE
       DFAC2   = DFAC1+ A1+ TWO*A1*A2*PHASE+ THREE*A1*A2*A3*PHASE**TWO
       DDFAC2  = DDFAC1+ TWO*A1*A2+ SIX*A1*A2*A3*PHASE
 
      
       OMEGA   = FAC1/FAC2+ 1.0D-6     
       DOMEGA  = (DFAC1*FAC2- FAC1*DFAC2)/(FAC2**TWO)
       DDOMEGA = ((DDFAC1*FAC2- FAC1*DDFAC2)*FAC2- TWO*
     1           (DFAC1*FAC2- FAC1*DFAC2)*DFAC2)/(FAC2**THREE)
      
      RETURN
      END
      
! ----------------------------------------------------------------------
!
      SUBROUTINE GEOMETRICFUNC(ALPHA,DALPHA,DDALPHA,PHASE)
!
! ----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
!     
      REAL*8 ALPHA,DALPHA,DDALPHA,PHASE
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,
     1 RP75=0.375D0,HALF=0.5D0)	  
!
      ALPHA  = 9.D0*PHASE
      DALPHA = 9.D0
      DDALPHA= 0.D0       

      RETURN 
      END 	  
!
! ==============================================================
! Eigenstrains from Voigt notation
! ==============================================================
!
      SUBROUTINE EIGOWN(EPS,EIGV,ALPHAE,ALPHAI,VECTI)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,
     1  TS=27.D0,THREE=3.D0,HALF=0.5D0,TOLER=1.0D-12,FOUR=4.D0,
     2  CNTN=100,TOLERE=1.0D-12)
      INTEGER I, J, K
      REAL*8 EPS(6), EIGV(3), ALPHAI(3), VECTI(3)
      REAL*8 PC, QC, ALPHAE, DISC, PI, CNT
!
       PI=FOUR*ATAN(ONE)
!   Scaling the strain vector
       VMAXE=MAXVAL(ABS(EPS))
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)/VMAXE
        END DO
       ENDIF
!    
!   Calculating eigenvalues
       VECTI(1)=EPS(1)+EPS(2)+EPS(3)
       VECTI(2)=EPS(1)*EPS(2)+EPS(2)*EPS(3)+EPS(1)*EPS(3)-
     1 EPS(4)**TWO/FOUR-EPS(5)**TWO/FOUR-EPS(6)**TWO/FOUR
       VECTI(3)=EPS(1)*EPS(2)*EPS(3)+EPS(4)*EPS(5)*EPS(6)/
     1 FOUR-EPS(1)*EPS(6)**TWO/FOUR-EPS(2)*EPS(5)**TWO/FOUR-
     2 EPS(3)*EPS(4)**TWO/FOUR
!
!   Depressed coefficients    
       PC=VECTI(2)-VECTI(1)**TWO/THREE
       QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
       DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)
!
       DO I=1,3
        EIGV(I)=ZERO
       END DO
       CNT=ZERO
       IF (ABS(DISC).LT.TOLER) THEN
        IF ((ABS(QC).LT.TOLER).AND.(ABS(PC).LT.TOLER)) THEN
         EIGV(1)=VECTI(1)/THREE
         EIGV(2)=VECTI(1)/THREE
         EIGV(3)=VECTI(1)/THREE
        ELSE
         EIGV(1)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(2)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(3)=THREE*QC/PC+VECTI(1)/THREE
         IF (EIGV(1).GT.EIGV(3)) THEN
          EONE=EIGV(1)
          EIGV(1)=EIGV(3)
          EIGV(3)=EONE
         ENDIF
        ENDIF
       ELSE
        DO I=1,3
         EIGV(I)=VECTI(1)/THREE+TWO*(MONE*PC/THREE)**HALF*
     1   COS(ONE/THREE*ACOS(MONE*QC/TWO*(TS/(MONE*PC**THREE))**
     2   HALF)+TWO*I*PI/THREE)
        END DO
       ENDIF
!       
       ALPHAE=ZERO
       IF ((EIGV(1)+EIGV(2)+EIGV(3)).GT.TOLER) THEN
        ALPHAE=ONE
       ENDIF
       DO K1=1,3
        ALPHAI(K1)=ZERO
        IF (EIGV(K1).GT.TOLER) THEN
         ALPHAI(K1)=ONE
        ENDIF
       END DO
!
!    Rescaling eigenvalues       
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)*VMAXE
        END DO
        DO K1=1,3
         EIGV(K1)=EIGV(K1)*VMAXE
        END DO
         VECTI(1)=EIGV(1)+EIGV(2)+EIGV(3)
         VECTI(2)=EIGV(1)*EIGV(2)+EIGV(1)*EIGV(3)+EIGV(3)*EIGV(2)
         VECTI(3)=EIGV(1)*EIGV(2)*EIGV(3)
       ENDIF
!   
       RETURN
       END
!
! ======================================================================
! U1: 3-node linear elements
! U2: 4-node bilinear elements
! ======================================================================
!	 	          
! ----------------------------------------------------------------------      
! Shape functions for U1: 3-node linear elements
! ----------------------------------------------------------------------      
      SUBROUTINE SHAPEFUNTRI(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(3,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = XI(1)
      AN(2) = XI(2)
      AN(3) = ONE-XI(1)-XI(2)
!
!     Derivatives of shape functions respect to local ccordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  ONE
      dNdxi(1,2) =  ZERO
      dNdxi(2,1) =  ZERO
      dNdxi(2,2) =  ONE
      dNdxi(3,1) =  MONE
      dNdxi(3,2) =  MONE
      RETURN
      END
!
! -----------------------------------------------------------------------      
! Shape functions for U2: 4-node bilinear elements
! -----------------------------------------------------------------------      
      SUBROUTINE SHAPEFUNQUAD(AN,dNdxi,xi)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)
!
!     Values of shape functions as a function of local coord.
      AN(1) = ONE/FOUR*(ONE-XI(1))*(ONE-XI(2))
      AN(2) = ONE/FOUR*(ONE+XI(1))*(ONE-XI(2))
      AN(3) = ONE/FOUR*(ONE+XI(1))*(ONE+XI(2))
      AN(4) = ONE/FOUR*(ONE-XI(1))*(ONE+XI(2))
!
!     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/FOUR*(ONE-XI(2))
      dNdxi(1,2) =  MONE/FOUR*(ONE-XI(1))
      dNdxi(2,1) =  ONE/FOUR*(ONE-XI(2))
      dNdxi(2,2) =  MONE/FOUR*(ONE+XI(1))
      dNdxi(3,1) =  ONE/FOUR*(ONE+XI(2))
      dNdxi(3,2) =  ONE/FOUR*(ONE+XI(1))
      dNdxi(4,1) =  MONE/FOUR*(ONE+XI(2))
      dNdxi(4,2) =  ONE/FOUR*(ONE-XI(1))
      RETURN
      END
!
! ======================================================================
! U3: 4-node linear tetrahedron elements
! U4: 8-node linear brick elements
! ======================================================================
! 
! ----------------------------------------------------------------------      
! Shape functions for U3: 4-node linear tetrahedron elements
! ----------------------------------------------------------------------  
      SUBROUTINE SHAPEFUNTET(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,3)
      Real*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = ONE-XI(1)-XI(2)-XI(3)
      AN(2) = XI(1)
      AN(3) = XI(2)
      AN(4) = XI(3)
!
!     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE
      dNdxi(1,2) =  MONE
      dNdxi(1,3) =  MONE
!
      dNdxi(2,1) =  ONE
      dNdxi(2,2) =  ZERO
      dNdxi(2,3) =  ZERO
!
      dNdxi(3,1) =  ZERO
      dNdxi(3,2) =  ONE
      dNdxi(3,3) =  ZERO
!
      dNdxi(4,1) =  ZERO
      dNdxi(4,2) =  ZERO
      dNdxi(4,3) =  ONE
!
      RETURN
      END
!
! --------------------------------------------------------------------      
! Shape functions for U4: 8-node linear brick elements
! --------------------------------------------------------------------     
      SUBROUTINE SHAPEFUNBRICK(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      REAL*8 AN(8),dNdxi(8,3)
      REAL*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0,EIGHT=8.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(2) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(3) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(4) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(5) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(6) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(7) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE+XI(3))
      AN(8) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE+XI(3))
      
!     Derivatives of shape functions respect to local coordinates
      DO I=1,8
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(1,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(1,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(2,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(2,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(2,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(3,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(3,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(3,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(4,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(4,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(4,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
      dNdxi(5,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(5,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(5,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(6,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(6,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(6,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(7,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(7,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(7,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(8,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(8,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(8,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
!      
      RETURN
      END
! 
! ======================================================================
! !!! NOTE: N_ELEM has to be changed according to the UEL !!!!!
! ======================================================================
!
      
       SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)

       INCLUDE 'ABA_PARAM.INC'

       CHARACTER*80 CMNAME

       DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3,3),
     3 DFGRD0(3, 3),DFGRD1(3,3)

       PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,SIX=6.D0,HALF=0.5D0,
     1 N_ELEM=5527,NSTV=4) 
       DATA NEWTON,TOLER/40,1.D-6/ 
!       
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
! 
! ---------------------------------------------------------------------- 
!
!  	   Stiffness tensor
       DDSDDE=0.D0
!	   
       NELEMAN=NOEL-N_ELEM
       IF (NPT.EQ.3) THEN
        NPT=4
       ELSEIF (NPT.EQ.4) THEN
        NPT=3
       ENDIF
!	   
       IF (NPT.EQ.7) THEN
        NPT=8
       ELSEIF (NPT.EQ.8) THEN
        NPT=7
       ENDIF
!      
       DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
       END DO
!       
       RETURN
       END 


	   