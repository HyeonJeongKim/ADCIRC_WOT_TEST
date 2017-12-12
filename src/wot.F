#ifdef CWOT
      SUBROUTINE WaveOverTopping(BARTYPE,EQTYPE,BSLOPE,BARHT,
     /                           BARCF1,BAR_DIST,BAR_ANGLE1,
     /                           ETA2,DP,BARINHT,
     /                           WAVE_H3,WAVE_A3,WAVE_T3, 
     /                           W_Q)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : G,PI,BARMIN
      
      IMPLICIT NONE
      REAL(SZ),intent(in)  ::BARTYPE,EQTYPE,BSLOPE,BARHT,BARCF1,
     /                       BAR_DIST,BAR_ANGLE1,ETA2,DP,BARINHT,
     /                       WAVE_H3,WAVE_A3,WAVE_T3
      REAL(SZ),intent(out) ::W_Q 
       
      REAL(SZ) ::W_H !impulsiveness parameter
      REAL(SZ) ::BAR_ANGLE
      REAL(SZ) ::VALID_WO,R_C 
      REAL(SZ) ::DEP_CHK
      REAL(SZ) ::WBK,S0
!     WBK: WAVE BREAKER  PARMETER, BARCF: BARRIER ROUGHNESS CF
!     S0=wave steepness 2*pi*wave height/g*wave period**2, BSLPOE: BARRIER SLOPE 
      REAL(SZ) ::DANG_A,DANG_B
      REAL(SZ) ::IN_WAVE_A2 ! FINAL INCIDENT WAVE ANGLE
      REAL(SZ) ::WQR,W_COEFF,D_H,BARCF
      REAL(SZ) ::W_Q_MAX,
     /           W_Q_5,W_Q_MAX_5,W_Q_DIFF,
     /           W_Q_7,WBK_DIFF,GV,WBK_5,WBK_7       
      REAL(SZ) ::Q0,QB     
      REAL(SZ) ::GB,RB,Z2P,BARX,BARDP           
c      write(*,*)'sub-----------------------------'  	
C....  hyl and sws 141215 ADDED EurOtop(2007) EQNS >>
C....  modyfied:  2015. 3. 20  
C....  reduction factor consideration due to incident wave angle (vertical dike)
           W_Q=0.D0               ! mean overtopped volume m3/s/m
           WBK=0.D0	              ! wave breaking parameter(IRIBARREN NUMBER) 
           S0=0.D0                ! barrier slope (e.g. slope=1:2 then input 2)              
           W_H=0.D0               ! impulsivensess parameter h*
           VALID_WO=0.D0          ! validitiy index of EurOtop Eqns application condition  
           R_C=0.D0               ! freeboard distance (crest level - swl) stilling water level                                                       
           IN_WAVE_A2=0.D0        ! initial wave incident angle  
           DANG_A=0.D0
           DANG_B=0.D0           
           DEP_CHK=0.D0              ! depth check for divergence problem in wet-dry area    DEP_CHK=DP+ETA2
           BAR_ANGLE=BAR_ANGLE1
C .... WOTR is defined in fort.15 (0: no computation,  1: wave overtopping)
C....  needed next variables in fort.14  
C....  BARTYPE: barrier (structure) type: VERTICAL WALL=1, SIMPLE SLOPED BARRIER=3 Incompleteness:(COMPOSITE VERTICAL WALL=2, COMPOSTIE SLOPE BARRIER=4) 
C .... EQTYPE: Type of computation; Deterministic EQ = 1 OR Probabilistic EQ = 2, default = 1 
C....  SLOPE: Slope of the sloped dike (only valid for BARTYPE = 3)
C....  BARCF: Friction coeff of barrier such as reduction factor of tetra-pod  (only valid for sloped dike)
C....  READ_INPUT.F  modified
C....  WAVE_H3   significant wave height (HS)
C....  WAVE_A3   mean wave DIReciton (DIR); NORTH=0 degree clockwise ==> only adcirc coipling
C....  WAVE_T3   mean wave periods(TM01)
C....  NNBB2: seaward barrier node number
C....  ETA2: SURFACE ELEVEATION
C....  DP:  water depth(fort.14)
C....  QN2:  normal flux through barrier  
 

! << Wave overtopping begins  
! depth corrected by introducing tide elevation 
         DEP_CHK=DP+ETA2
        IF(WAVE_H3.GE.0.1D0 .AND. DEP_CHK.GT.0.0)THEN  !SWITCH WAVE OVERTOPPING RATE       
! WBK < 2: BREKING WAVE, WBK >= 2: NON-BREAKING WAVE

! calculating Barrier distance in x & y-dir.
! HJK, Barrier의 마지막 node에 대한 정의가 필요함!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       barrier와 파방향에 대해 계산하지만, 파 방향이 barrier와 반대일경우에 대한 정의 필요함

C khj & sws{          
         IF(BAR_ANGLE.LT.0.D0)BAR_ANGLE=BAR_ANGLE+360.D0   !Barrier - coordinate counterclockwise, east=0
         
         IN_WAVE_A2=ABS(450.D0-WAVE_A3)
         IF(IN_WAVE_A2.GE.360.D0) 
     /      IN_WAVE_A2=IN_WAVE_A2-360.D0       !SWAN WAVE DIR - coordinate counterclockwise, east=0
           
         IF(BAR_ANGLE.GE.0.D0.AND.BAR_ANGLE.LT.180.D0)THEN
           DANG_A=BAR_ANGLE
           DANG_B=BAR_ANGLE+180.D0
         ENDIF
         IF(BAR_ANGLE.GE.180.D0.AND.BAR_ANGLE.LT.360.D0)THEN
           DANG_A=BAR_ANGLE
           DANG_B=BAR_ANGLE+180.D0-360.D0
         ENDIF      
         IN_WAVE_A2=ABS(IN_WAVE_A2-MAX(DANG_A,DANG_B))
         IF(IN_WAVE_A2.GT.90.D0) 
     /      IN_WAVE_A2=180.D0-IN_WAVE_A2
            IN_WAVE_A2=90.D0-IN_WAVE_A2  !최종 wave dir 설정, 0degree 일때 = barrier와 직각, 90degree일때 = 평행
C khj & sws} 
          
! SLOPE TYPE BARRIER 
          S0= (2.D0*PI*WAVE_H3)/(G*WAVE_T3**2) ! WAVE STEEPNESS
! breaker parameter(IRIBARREN NUMBER) 
 
	      WBK=BSLOPE/SQRT(S0)
          R_C=BARINHT-ETA2   ! R_C is crest freeboard of structure

      IF(R_C .GT. 0.0)THEN !CASE1
! >>>>>>>>>>>>WAVE OVERTOPPING CASE1: POSITIVE FREE BOARD, ONLY WAVE OVERETOPPING
! REF. EUROTOP 2007, P130~ : VERTICAL WALL
! WAVE REDUCTION FACTOR non-impulsive conditions
!        IF(WAVE_H3 .GE. BARMIN)THEN  ! start wave overtopping calculation
        VALID_WO=R_C/WAVE_H3
        W_H=1.35D0*((DP+ETA2)/WAVE_H3)*  ! impulsivensess parameter h*               
     &            ((2*PI*(DP+ETA2))/(G*WAVE_T3**2)) !EurOtop EQ. 7.1

C.... BARTYPE: VERTICAL WALL:1,  SIMPLE SLOPED BARRIER:3  
C.... INCOMPLETENESS:(COMPOSITE VERTICAL WALL:2, COMPOSTIE SLOPE BARRIER:4) 
C.... EQTYPE: Deterministic EQ=1 OR Probabilistic EQ=2, DEFAULT=1 
C.... SLOPE: slope of barrier (only valid for BARTYPE=3)
C.... BARCF: barrier friction coeff  

c        write(*,*)'sub-----------'
C***********************************************************************
c***********************************************************************
c***********************************************************************
        IF(BARTYPE .EQ. 1.0D0) THEN
c        write(*,*)'bartype=1'
            WQR=1.D0   ! in wave_a2 = 0 case, i.e. zero oblique incident wave condition 
! ONLY Probabilistic , non-impulsive conditions
!         IF( IN_WAVE_A2 .GT. 0.D0 .AND.  IN_WAVE_A2 .LT. 45.D0)  ! eurotop 2007
! improvement of OBLIQUE WAVE REDUCTION FACTOR  45-->80
! Kim et al., 2010, KSCOE 22(3) pp149-155
!         IF( IN_WAVE_A2 .GT. 45.D0)WQR=0.72D0 
         IF( IN_WAVE_A2 .GT. 0.0 .AND. IN_WAVE_A2 .LT. 80.0)THEN      !EurOtop EQ. 7.16
           WQR=1.D0-0.0062D0*IN_WAVE_A2                               ! WQR IS OBLIQUE WAVE REDUCTION FACTOR
         ENDIF
        ! improvement of OBLIQUE WAVE REDUCTION FACTOR
        ! Kim et al., 2010, KSCOE 22(3) pp149-155
        !       IF( IN_WAVE_A2 .GT. 0.D0 .AND.  IN_WAVE_A2 .LT. 80.D0)    
        !    & WQR=1.D0-0.0062D0* IN_WAVE_A2
       IF(W_H .NE. 0.D0)THEN                                          ! Don't allow  NaN,  prevent NEGATIVE DEPTHS in tidal flat  
        IF(VALID_WO .GT. 0.1 .AND. VALID_WO .LT. 3.5)THEN             ! Don't allow  NaN, if VALID_WO = 0 
         IF(W_H .GT. 0.3D0) THEN                                      ! non-impulsive conditions
          IF(EQTYPE .EQ. 2.0)W_COEFF=-2.6D0/WQR                       ! Probabilistic ,non-impulsive conditions,  CONSIDERING OBLIQUE  WAVE 
          IF(EQTYPE .EQ. 1.0)W_COEFF=-1.8D0                           ! Deterministic, non-impulsive conditions 
           W_Q=(0.04D0*exp(W_COEFF*VALID_WO))*SQRT(G*WAVE_H3**3) !EurOtop  EQ 7.3~4                            
         ENDIF
        ENDIF          
       IF(W_H .LE. 0.20 .AND. IN_WAVE_A2 .GT. 0.0 .AND.
     &    IN_WAVE_A2.LT.15.0 )THEN                                    !impulsive conditions,  LT 15 DEG
         IF(W_H*VALID_WO .GT. 0.030 .and. W_H*VALID_WO .LT. 1.0)THEN  ! Don't allow  NaN, if VALID_WO = 0                                 
           IF(EQTYPE .EQ. 2.0)W_COEFF=1.5D0                           ! Probabilistic , impulsive conditions
           IF(EQTYPE .EQ. 1.0)W_COEFF=2.8D0                           ! Deterministic, impulsive conditions 
             W_Q=(W_COEFF*0.0001D0)*((W_H*VALID_WO)**-3.1)*
     &           ((W_H**2)*SQRT(G*(DP+ETA2)**3.D0))              !EurOtop EQ. 7.6~7          
         ENDIF
         IF(W_H*VALID_WO .LT. 0.02)THEN                               !broken wave conditions
           IF(EQTYPE .EQ. 2.0)W_COEFF=2.7D0                           ! Probabilistic , non-impulsive conditions
           IF(EQTYPE .EQ. 1.0)W_COEFF=3.8D0                           ! Deterministic, non-impulsive conditions 
             W_Q=(W_COEFF*0.0001D0)*((W_H*VALID_WO)**-2.7)*
     &           ((W_H**2)*SQRT(G*(DP+ETA2)**3))                 !EurOtop EQ. 7.8~9          
         ENDIF
       ENDIF       
      ! PROBABILISTIC,  IMPULSIVE CONDITION,  OBLIQUE WAVE EFFECT 
      ! EUROTOP P144
      ! Incident angles only for (15, 30, 60 deg)     
       IF(EQTYPE .EQ. 2.0) THEN           
        IF(IN_WAVE_A2.EQ.15.0 .AND. W_H*VALID_WO.GE.0.2)THEN 
          W_Q=0.000058D0*((W_H*VALID_WO)**-3.7)*
     &          ((W_H**2)*SQRT(G*(DP+ETA2)**3))                   !EurOtop  EQ. 7.17-1
        ENDIF    
        IF(IN_WAVE_A2.EQ.15.0 .AND. W_H*VALID_WO.LT.0.2)THEN  
          W_Q=(1.5*0.0001D0)*((W_H*VALID_WO)**-3.1)*
     &        ((W_H**2)*SQRT(G*(DP+ETA2)**3.1))                   !EurOtop EQ. 7.6
        ENDIF
        IF(IN_WAVE_A2.EQ.30.0 .AND. W_H*VALID_WO.GE.0.07)THEN
          W_Q=0.000008D0*((W_H*VALID_WO)**-4.2)*
     &        ((W_H**2)*SQRT(G*(DP+ETA2)**3))                     !EurOtop EQ. 7.17-2
        ENDIF    
        IF(IN_WAVE_A2.EQ.60.0 .AND. W_H*VALID_WO.GE.0.07)THEN  
          W_Q=(0.04D0*exp(-3.61D0*(VALID_WO)))*SQRT(G*WAVE_H3**3.D0) ! EurOtop EQ 7.15 IF WAVE ANGELE 60 ~  ,  -2.6/0.72=-3.61          
        ENDIF    
       ENDIF 
      ! DETERMINISTIC,  IMPULSIVE CONDITION
       IF(EQTYPE .EQ. 1.0) THEN 
        IF(IN_WAVE_A2 .EQ. 15.0 .AND. W_H*VALID_WO.GE.0.2)THEN 
          W_Q=(2.8*0.0001D0)*((W_H*VALID_WO)**-3.1)*
     &            ((W_H**2)*SQRT(G*(DP+ETA2)**3))                   !EurOtop EQ. 7.7
        ENDIF    
        IF(IN_WAVE_A2.EQ.30.0 .AND. W_H*VALID_WO.GE.0.07)THEN  
          W_Q=(1.5*0.0001D0)*((W_H*VALID_WO)**-3.1)*
     &            ((W_H**2)*SQRT(G*(DP+ETA2)**3.1))                 !EurOtop EQ 7.17 BETA =15 DEG ==>EQ. 7.6
        ENDIF    
        IF(IN_WAVE_A2.EQ.60.0 .AND. W_H*VALID_WO.GE.0.07)THEN 
          W_Q=(0.04D0*exp(-1.8D0*(VALID_WO)))*   
     &            SQRT(G*WAVE_H3**3)                                !EurOtop  EQ. 7.4
        ENDIF 
       ENDIF        
      ENDIF !W_H .NE. 0.D0
c        write(*,*)'w_q',W_Q      
      ENDIF        
C***********************************************************************
c***********************************************************************
c***********************************************************************
        IF(BARTYPE .EQ. 2)THEN        
 	      D_H=1.35D0*((DP+ETA2-BARHT)/WAVE_H3)*
     &        ((2.D0*PI*(DP+ETA2))/(G*WAVE_T3**2))                      !EurOtop  EQ 7.2   
            IF(D_H.NE.0.0 .AND. W_H.NE.0.0 )THEN                        ! Don't allow  NaN  
             IF(D_H .LE. 0.20)THEN                                      ! impulsive conditions   
	          IF(D_H*VALID_WO.GT.0.05D0 .AND. D_H*VALID_WO.LT.1.D0 .AND.  
     &           W_H.GE.0.3D0) THEN                                     ! Don't allow  NaN           
  	           IF(EQTYPE .EQ. 2.0) W_Q=(0.00041D0*
     &                                 (D_H*VALID_WO)**-2.9)*
     &                                 (D_H**2)*sqrt(G*(DP+ETA2)**3) !EurOtop EQ 7.13   ! Probabilistic, non-impulsive conditions 
               IF(EQTYPE .EQ. 1.0) W_Q=(0.00078D0*
     &                                 (D_H*VALID_WO)**-2.6)*
     &                                 (D_H**2)*sqrt(G*(DP+ETA2)**3) ! EurOtop EQ 7.14   !  Deterministic , non-impulsive conditions 
	          ENDIF
             ENDIF

             IF(D_H .GT. 0.30)THEN   ! non-impulsive conditions
              IF(VALID_WO.GT.0.10 .and. VALID_WO.LT.3.50)THEN 
               IF(EQTYPE .EQ. 2.0)W_COEFF=-2.6D0                  ! Probabilistic , non-impulsive conditions,  CONSIDERING OBLIQUE  WAVE
               IF(EQTYPE .EQ. 1.0)W_COEFF=-1.8D0                  ! Deterministic, non-impulsive conditions 
                W_Q=(0.04D0*exp(W_COEFF*(VALID_WO)))*
     &               SQRT(G*WAVE_H3**3)                         !EurOtop  EQ 7.3~4
              ENDIF
             ENDIF 
            ENDIF   
        ENDIF
C***********************************************************************
c***********************************************************************
c***********************************************************************        
        IF(BARTYPE .EQ. 3)THEN        
    
       WBK_7=7                   ! wave breaking index 7      
       WBK_5=5                   ! wave breaking index 5
 ! BARTYPE 3 --> SIMPLIE SLOPE TYPE 
       BARCF=BARCF1
       WQR=1.D0  ! DEFAULT WAVE REDUCTION FACTOR(NO EFFECT)  
       IF(BARCF.EQ.0.0)BARCF=1.d0  ! default 1
       IF(BARCF.GT.0.6)THEN    
        IF( IN_WAVE_A2.GE.0.0 .AND. IN_WAVE_A2.LE.80.0)   !EurOtop EQ. 5.24
     &    WQR=1.D0-0.0033*IN_WAVE_A2                      ! WQR IS WAVE REDUCTION FACTOR:    EFFECT OF OBLIQUE WAVE
        IF(IN_WAVE_A2.GT.80.0)
     &    WQR=0.736D0 
       ENDIF

! classified according to BARCF 
       IF(BARCF.LE.0.60)THEN  
        IF(IN_WAVE_A2.GE.0.0 .AND. IN_WAVE_A2.LE.80.0)   !EQ. 6.8,   BARCF=1 IS NOT APPLIED
     &    WQR=1.D0-0.0063*IN_WAVE_A2                       ! WQR IS WAVE REDUCTION FACTOR
        IF(IN_WAVE_A2.GT.80.0)
     &    WQR=0.496D0                                     !1-0.0063*80
       ENDIF	 	

       IF(WBK.LE.5.0)THEN                                ! BRAEKING WAVE CONDITION
        IF(EQTYPE.EQ.2.0)W_COEFF=-4.75D0                 ! Probabilistic 
        IF(EQTYPE.EQ.1.0)W_COEFF=-4.3D0                  ! Deterministic     
                                                         ! -0.92 IS NORMALLY DISTRIBUTED STOCHASTIC FUNCTION ,10**-0.92=0.12      
	     GV=1.D0                                         ! NO EFFECT VERTICAL WALL SLOPE
	     W_Q=(0.067D0/SQRT(BSLOPE))*WBK*                 ! TAN_A , SIMPLE SLOPE
     &        EXP((W_COEFF*VALID_WO)/(WBK*BARCF*WQR*GV))*
     &        SQRT(G*WAVE_H3**3)                         !EurOtop EQ. 5.8~9
                                         ! MAXIMUM OVER TOPPING RATE
        IF(EQTYPE.EQ.2.0)W_COEFF=-2.6D0                  ! Probabilistic  
        IF(EQTYPE.EQ.1.0)W_COEFF=-2.3D0                  ! Deterministic       
         W_Q_MAX=0.2D0*EXP((W_COEFF*R_C)/(WAVE_H3* 
     &       BARCF*WQR))*SQRT(G*WAVE_H3**3)              ! EurOtop EQ.5.8~9 MAXIMUM OVER TOPPING RATE
                                         !  w_q cannot exceed w_q_max      
        IF(W_Q.GE.W_Q_MAX) W_Q=W_Q_MAX 
       ENDIF  ! END BRAEKING WAVE

! NON-BRAEKING WAVE 

!       WBK BETWEEN 5 AND 7 LINEAR INTERPOLATION PART
!        EUROTOP P82
       IF(WBK.GT.5.0 .AND. WBK.LT.7.0)THEN         
         WBK_DIFF=WBK-5.0
                                                         ! WBK= 5 PART
        IF(EQTYPE.EQ.2.0)W_COEFF=-4.75D0                      ! Probabilistic 
        IF(EQTYPE.EQ.1.0)W_COEFF=-4.3D0                       ! Deterministic     
                                                              ! -0.92 IS NORMALLY DISTRIBUTED STOCHASTIC FUNCTION ,10**-0.92=0.12      
	     GV=1.D0                                              ! NO EFFECT VERTICAL WALL SLOPE
	     W_Q_5=(0.067D0/SQRT(BSLOPE))*WBK_5*                  ! ASSUME TAN_A , SIMPLE SLOPE
     &        EXP((W_COEFF*VALID_WO)/(WBK_5*BARCF*WQR*GV))*
     &        SQRT(G*WAVE_H3**3)                              !EurOtop EQ. 5.8~9
                                                              ! MAXIMUM OVER TOPPING RATE  !  
        IF(EQTYPE.EQ.2.0)W_COEFF=-2.6D0                       ! Probabilistic  
        IF(EQTYPE.EQ.1.0)W_COEFF=-2.3D0                       ! Deterministic       
	     W_Q_MAX_5= 0.2D0*EXP((W_COEFF*R_C)/
     &             (WAVE_H3*BARCF*WQR))*SQRT(G*WAVE_H3**3)    ! EurOtop EQ.5.8~9 MAXIMUM OVER TOPPING RATE
        IF(W_Q_5 .GT. W_Q_MAX_5)W_Q_5=W_Q_MAX_5 
 
                                                          ! WBK= 7 PART
                                                               ! ASSUME WBK GE. 6.D0   ORG PARAMETER = WBK .GT. 7
                                                               !	  IF(WBK .GT. 5.D0)THEN  ! VERY SHALLOW FORESHORE
        IF(EQTYPE.EQ.2.0)W_COEFF=10.0**-0.92D0                 ! Probabilistic 
        IF(EQTYPE.EQ.1.0)W_COEFF= 0.21D0                       ! Deterministic 
                                                               ! -0.92 IS NORMALLY DISTRIBUTED STOCHASTIC FUNCTION ,10**-0.92=0.12      
	     W_Q_7= W_COEFF*EXP(-VALID_WO/(BARCF*WQR*
     &          (0.33D0+0.022D0*WBK_7)))*SQRT(G*WAVE_H3**3)    !EurOtop  EQ.5.10~11 MAXIMUM OVER TOPPING RATE
         W_Q_DIFF=W_Q_7 - W_Q_5
         W_Q=W_Q_5+(W_Q_DIFF*(WBK_DIFF/2.0))
       ENDIF  ! ! WBK BETWEEN 5 AND 7 LINEAR INTERPOLATION

                                                                ! ASSUME WBK GE. 6.D0   ORG PARAMETER = WBK .GT. 7
       IF(WBK.GE.7.0)THEN                                  ! VERY SHALLOW FORESHORE
        IF(EQTYPE .EQ. 2.0)W_COEFF=10.D0**-0.92D0                  ! Probabilistic 
        IF(EQTYPE .EQ. 1.0)W_COEFF= 0.21D0                         ! Deterministic 
                                                                   ! -0.92 IS NORMALLY DISTRIBUTED STOCHASTIC FUNCTION ,10**-0.92=0.12      
         W_Q=W_COEFF*EXP(-VALID_WO/(BARCF*WQR*
     &        (0.33D0+0.022D0*WBK)))*SQRT(G*WAVE_H3**3)      !EurOtop  EQ.5.10~11 MAXIMUM OVER TOPPING RATE
       ENDIF  ! NON-BRAEKING WAVE 
        ENDIF 
C***********************************************************************
c***********************************************************************
c***********************************************************************        
       IF(BARTYPE .EQ. 30)THEN         
 ! BARTYPE 30 --> SIMPLIE SLOPE TYPE  ! ONLY OWEN(1980) EQ. EUROTOP P82

           IF(BSLOPE .EQ. 1.0) QB= 20.1D0; Q0=0.00794D0
           IF(BSLOPE .EQ. 1.5) QB= 19.9D0; Q0=0.00884D0
           IF(BSLOPE .EQ. 2.0) QB= 21.6D0; Q0=0.00939D0
           IF(BSLOPE .EQ. 2.5) QB= 24.5D0; Q0=0.0103D0
           IF(BSLOPE .EQ. 3.0) QB= 28.7D0; Q0=0.0109D0
           IF(BSLOPE .EQ. 3.5) QB= 34.1D0; Q0=0.0112D0
           IF(BSLOPE .EQ. 4.0) QB= 41.0D0; Q0=0.0116D0
           IF(BSLOPE .EQ. 4.5) QB= 47.7D0; Q0=0.0120D0
           IF(BSLOPE .EQ. 5.0) QB= 55.6D0; Q0=0.0131D0

C    BRITISH GUIDELINES OWEN(1980) EQ. ONLY STRAIGHT SLOPE
!  W_Q(I)=Q0*EXP(-BSLOPE(I)*(R_C/(WAVE_T3*SQRT(G*WAVE_H3)))
!     & *(WAVE_T3*G*WAVE_H3)  !EurOtop EQ. 5.12
          W_Q=(Q0/SQRT(S0))*EXP(-BSLOPE*(R_C/
     &    (WAVE_H3*SQRT(S0))))*SQRT(G*WAVE_H3**3)  !EurOtop EQ. 5.13 ,  CONSIDERING WAVE BREAKING 
       ENDIF         
C***********************************************************************
c***********************************************************************
c***********************************************************************        
       IF(BARTYPE .EQ. 31)THEN         
! #################BARTYPE=31 not used ########################
! BARTYPE 31--> INCLUDE BERM     
! ASSUME SIMPLE SLOPE(AVERAGE SLOPE)
	   	      ! TAN (ALPHA)
! WBK=BSLOPE(I)/SQRT(S0) ! breaker parameter(IRIBARREN NUMBER) 
! WBK < 2 :BREKING WAVE, WBK >= 2 :NON-BREAKING WAVE
        IF(BARCF .GT. 0.60)THEN    
         IF(IN_WAVE_A2.GE.0.0 .AND. IN_WAVE_A2.LE.80.0)   !EurOtop EQ. 5.24
     &     WQR=1.D0-0.0033*IN_WAVE_A2                     ! WQR IS WAVE REDUCTION FACTOR:    EFFECT OF OBLIQUE WAVE
         IF(IN_WAVE_A2.GT.80.0)
     &     WQR=0.736D0 
        ENDIF
!BARCF   
        IF(BARCF .LE. 0.60)THEN   
         IF(IN_WAVE_A2.GE.0.0 .AND. IN_WAVE_A2.LE.80.0)   !EurOtop EQ. 6.8  ! BUT, BARCF=1 IS NOT APPLIED
     &     WQR=1.D0-0.0063*IN_WAVE_A2                     ! WQR IS WAVE REDUCTION FACTOR
         IF(IN_WAVE_A2.GT.80.0)
     &     WQR=0.496D0  !1-0.0063*80
        ENDIF	 	
                                
        IF(WBK .LE. 5.0)THEN      ! [TAW, p25 breaking condition is 2] however, EurOtop p79, it is 5  
! BRAEKING WAVE 
         IF(EQTYPE.EQ.2.0)W_COEFF=-4.75D0         ! Probabilistic 
         IF(EQTYPE.EQ.1.0)W_COEFF=-4.3D0           ! Deterministic   
          GV=1.D0                                 ! NO EFFECT VERTICAL WALL SLOPE
          RB=0.4D0                                ! TAW EQ. 14, P18B=0.4L_BERM
          Z2P=8.D0*WAVE_H3*BSLOPE                 ! SIMPLE EQ.TAW MANUAL EQ .7
         BARDP=(DP+ETA2)-BARHT
         IF(BARDP.GT.0.0  .AND. BARDP.LT.Z2P)BARX=Z2P
         IF(BARDP.LT.0.D0 .AND. BARDP.LT.-2.D0*WAVE_H3)BARX=2.D0*WAVE_H3
         IF(BARDP.GT.Z2P.OR.BARX.LT.-2.D0*WAVE_H3)BARX=1.0
         GB=1.D0-RB*(0.5D0+0.5D0*COS(PI*(BARDP/BARX)))
         W_Q=(0.067D0/SQRT(BSLOPE))*GB*WBK*    ! ASSUME TAN_A , SIMPLE SLOPE
     &        EXP((W_COEFF*VALID_WO)/(WBK*GB*BARCF*WQR*GV))*
     &        SQRT(G*WAVE_H3**3)   !EurOtop EQ. 5.8~9  INCLOUDE BERM 
         IF(EQTYPE.EQ.2.0)W_COEFF=-2.6D0   ! Probabilistic  
         IF(EQTYPE.EQ.1.0)W_COEFF=-2.3D0   ! Deterministic       
	     W_Q_MAX= 0.2D0*EXP((W_COEFF*R_C)/(WAVE_H3*
     &            BARCF*WQR))*SQRT(G*WAVE_H3**3)    ! EQ.5.8~9 MAXIMUM OVER TOPPING RATE 
! w_q cannot exceed w_q_max  
         IF(W_Q.GE.W_Q_MAX)W_Q=W_Q_MAX 
        ENDIF

! ASSUME WBK GE. 6.D0   ORG PARAMETER = WBK .GT. 7
!	  IF(WBK .GT. 7.D0)THEN  ! VERY SHALLOW FORESHORE
        IF(WBK.GE.7.0)THEN  ! VERY SHALLOW FORESHORE
         IF(EQTYPE.EQ.2.0)W_COEFF=10.D0**-0.92D0   ! Probabilistic 
         IF(EQTYPE.EQ.1.0)W_COEFF= 0.21D0  ! Deterministic 
       ! -0.92 IS NORMALLY DISTRIBUTED STOCHASTIC FUNCTION ,10**-0.92=0.12      
	     W_Q=W_COEFF*EXP(-VALID_WO/(BARCF*WQR*
     &          (0.33D0+0.022D0*WBK)))*SQRT(G*WAVE_H3**3)      ! EurOtop EQ.5.10~11 MAXIMUM OVER TOPPING RATE
        ENDIF  ! NON-BRAEKING WAVE 
       ENDIF  
      ENDIF !CASE1
C***********************************************************************
c***********************************************************************
c***********************************************************************         
C>>>>>>>>>>>>>WAVE OVERTOPPING CASE2: ZERO FREEBOARD, ONLY WAVE OVERETOPPING
      IF(R_C .LE. 0.0 .AND. R_C .GT. -BARMIN)THEN
           IF(BARTYPE.EQ.1 .OR. BARTYPE.EQ.2)THEN
	         W_Q=0.062D0*SQRT(G*WAVE_H3**3)   !EurOtop  EQ. 7.5  ONLY PROBABILISTIC, NON-IMPULSIVE CONDITION    
           ENDIF                                        ! EXCEPT +- 0.0062
         
           IF(BARTYPE.EQ.3)THEN
            IF(WBK.LE.2.0) W_Q=(0.0537D0*WBK)*(SQRT(G*WAVE_H3**3))   ! BREAKING WAVE      
	        IF(WBK.GE.2.0) W_Q=(0.136D0-(0.226D0/WBK**3))*           !NON-BREAKING WAVE
     &                          (SQRT(G*WAVE_H3**3))                 !EurOtop EQ. 5.14                    
           ENDIF !BARTYPE
      ENDIF !CASE2

C >>>>>>>>>>>>>WAVE OVERTOPPING CASE4: COMBINED FREEBOARD ,OVERFLOWS AND WAVE OVERETOPPING 
       IF(R_C .LT. 0.D0)THEN ! WAVE OVERTOPPING EFFECT VERY SMALL , IGNORE WAVE EFFECT
               W_Q=0.D0
       ENDIF !CASE4     
         
       
cc       IF(BAR_DIST .gt. 5. AND. W_Q .GT. 0.D0 ) THEN ! coupling TIME interval         
ccc        IF(MOD(TIME,RSTIMINC) .EQ. 0.0)THEN

                     
ccc        ENDIF
cc       ENDIF
      ENDIF  
          
      END SUBROUTINE WaveOverTopping 

c#####################################################################################################
      SUBROUTINE Hydraulic_gradient(BAR_DIST,W_Q,W_Q_SUM,W_Q_OUT,W_Q_SUM_OUT)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : G,PI,tdt,BAR_CFN,H0
      
      IMPLICIT NONE
      REAL(SZ),intent(in)  ::BAR_DIST,W_Q,W_Q_SUM
      REAL(SZ),intent(out) ::W_Q_OUT,W_Q_SUM_OUT    

! added HYL  HYDRAULIC GRADIENT CALCULATION
      REAL(SZ) :: BAR_H1,BAR_H2,BAR_WQ
      REAL(SZ) :: BAR_H1F, BAR_H2F     
      REAL(SZ) :: BAR_EH1,BAR_HL 
      REAL(SZ) :: BAR_EH2
      REAL(SZ) :: BAR_V1
      REAL(SZ) :: BAR_V2
      REAL(SZ) :: BAR_R1
      REAL(SZ) :: BAR_R2
      REAL(SZ) :: BAR_SF1
      REAL(SZ) :: BAR_SF2
      REAL(SZ) :: BAR_WP1
      REAL(SZ) :: BAR_WP2
      REAL(SZ) :: BAR_AREA1 
      REAL(SZ) :: BAR_AREA2
      REAL(SZ) :: BAR_FY,BAR_FYD, BAR_ALPHA 
      REAL(SZ) :: BAR_B1,BAR_B2
      REAL(SZ) :: BAR_H2D    
      REAL(SZ) :: BAR_FYD_1,BAR_FYD_2,BAR_FYD_3,BAR_FYD_4 
      REAL(SZ) :: BAR_I,DPDY          
          
      INTEGER BAR_IT	
      
! added HYL HYDRAULIC GRADIENT CALCULATION
! hyl added>      
C***********************************************************************
c***********************************************************************
c*********************************************************************** 
c**************HYDRAULIC GRADIENT CALCULATION Assumption ***************
C        STANDARD STEP METHOD, calculate downstream depth with known upstream depth  
C        USING NEWTON-RAPHSON METHOD
C        BAR_CFN: BARRIER ROUGHNESS (MANNING'S N) ==>0.015
C        BAR_H: BARRIER INITIAL HEIGHT (NO HEAD LOSS)m
C        BAR_VEL :  Velocity on barreir  
C        HEADLOSS: HEAD LOSS between two pairing points
C        BAR_HL: BAR_H-HEAD LOSS 
C        BAR_Q: Q considered HEAD LOSS         
C        added hyl 2016.4.25           
!       if barrier distnace < 1, then no headloss computation (too short distance assumption)
         ! with known BAR_H1, calculate BAR_H2 
         ! downSTREAM = BAR_H2   land-side depth 
         ! upSTRAM=BAR_H1  sea-side depth     
        BAR_WQ=0.0                       
        BAR_EH1=0.0 
        BAR_EH2=0.0 
        BAR_V1=0.0 
        BAR_V2=0.0 
        BAR_R1=0.0 
        BAR_R2=0.0 
        BAR_SF1=0.0 
        BAR_SF2=0.0 
        BAR_WP1=0.0 
        BAR_WP2=0.0 
        BAR_AREA1=0.0 
        BAR_AREA2=0.0 
        BAR_H2=0.D0
        BAR_H1=0.D0
       
        BAR_FY=0.0    
        BAR_H2D=0.0
        BAR_FYD=0.0
        BAR_B2=1.0
        BAR_B1=1.0
        BAR_IT=0
        BAR_ALPHA=1.0  ! energy correction factor 
         
        BAR_H1=((W_Q+W_Q_SUM)/SQRT(2.D0*G))**tdt  
         
        IF(BAR_H1 .LT. H0)THEN
         W_Q_SUM_OUT=W_Q_SUM+W_Q
         W_Q_OUT=0.D0
        
         ELSE
           
	   BAR_AREA1=BAR_H1*BAR_B1   ! y*b B=UNIT AREA                
         BAR_WP1=2.D0*BAR_H1+(BAR_AREA1/BAR_H1)  ! P1=2*Y+B
         BAR_R1=BAR_AREA1 / BAR_WP1
         BAR_V1=(W_Q+W_Q_SUM) / BAR_AREA1 ! mean vel          
         BAR_SF1=(BAR_V1*BAR_CFN)**2 / BAR_R1**1.3333 ! hydraulic radius
 
	   
1        DO  4 BAR_I =1,10	 ! 10 is divide dist,  dx*10 =dist 
         
         BAR_EH1=BAR_H1+(BAR_ALPHA*(BAR_V1**2)/(2.*G)) !
	   BAR_H2 = BAR_H1
	                       
         BAR_IT=0.
2        BAR_IT=BAR_IT+1   
          
         BAR_AREA2=BAR_H2*BAR_B2   ! y*b B=UNIT AREA   
         BAR_WP2=2.D0*BAR_H2+(BAR_AREA2/BAR_H2)   
         BAR_R2=BAR_AREA2 / BAR_WP2 
         BAR_V2=(W_Q+W_Q_SUM) / BAR_AREA2 ! mean vel 
         BAR_EH2=BAR_ALPHA *(BAR_V2**2) / (2.*G)  ! vel head       
         BAR_SF2=(BAR_V2*BAR_CFN)**2 / BAR_R2**1.3333 ! hydraulic radius
         BAR_AREA2=BAR_H2*BAR_B2   ! y*b B=UNIT AREA		  	 	   
			 
		 BAR_FYD_1=(BAR_ALPHA*BAR_B2*(W_Q+W_Q_SUM)**2)/(G*BAR_AREA2**3)
		 BAR_FYD_2=BAR_SF2*(BAR_B2/BAR_AREA2)
         BAR_FYD_3=(BAR_SF2/BAR_R2)
		 BAR_FYD_4=(BAR_B2/BAR_WP2)-(BAR_AREA2/BAR_WP2**2)
		 DPDY=2.D0  ! 2*SQRT(1+S**2)	 
				 
         BAR_FY=BAR_H2+BAR_EH2+0.5*BAR_SF2
     &           *(BAR_DIST*0.1D0)
     &           -BAR_EH1+0.5*BAR_SF1
     &           *(BAR_DIST*0.1D0)  ! Z1+Z2
                           
         BAR_FYD=1.D0
     &		    -BAR_FYD_1-(BAR_DIST*0.1D0)
     &          *(BAR_FYD_2+tdt*(BAR_FYD_3*BAR_FYD_4)*DPDY)  

	   BAR_H2D=BAR_H2-(BAR_FY/BAR_FYD)                      
          				
		  IF(ABS(BAR_H2D-BAR_H2).LE. 0.0001 .and. BAR_IT.le. 20)THEN	! error 0.001 m
		    IF(BAR_I .EQ. 10 .AND. BAR_H2 .GT. H0)THEN		 
		     W_Q_SUM_OUT=0.d0		
!		   write(99999,*)bar_h1, bar_h2, bar_it  
	       GOTO 10  		      
            ENDIF                
            GOTO 3 
	    ENDIF
	         
       	  IF(BAR_H2 .LT. 0. .OR. BAR_H1.LT.0.)THEN  ! NAN       	  
             BAR_H2=0.D0
		     GOTO 10       		
		  ENDIF       		
       	
	     BAR_H2=BAR_H2D 	            		       			      
       	  IF(BAR_IT.GT. 20)THEN     ! iteration max. 20						
             BAR_H2=0.D0
		     GOTO 10
		  ENDIF
				
	    GOTO 2		
3			 BAR_H1=BAR_H2 
			 BAR_SF1=BAR_SF2 	
			 
		                  
4        CONTINUE   
           
10       BAR_WQ=BAR_H2*SQRT(2.D0*G*BAR_H2)                        
         W_Q_OUT =BAR_WQ    
       ENDIF
       
       
C******************************HYDRAULIC GRADIENT CALCULATION*****************************************        
      END SUBROUTINE Hydraulic_gradient         
      
#endif