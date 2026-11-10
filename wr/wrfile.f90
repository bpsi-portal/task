! wrfile.f90

MODULE wrfile
  
  USE bpsd_kinds,     ONLY : dp
  USE bpsd_constants, ONLY : PI, AEE, AME, RKBL, VC, EPS0
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC wr_save
  PUBLIC wr_load
  PUBLIC wr_write
  PUBLIC find_root_bisection
  PUBLIC eccd_write
  PUBLIC fg_eval

CONTAINS

  FUNCTION fg_eval(gamma, y_res, n_para) RESULT(fg)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: gamma, y_res, n_para
    REAL(dp) :: fg
    
    fg = gamma**2 - 1.0_dp - ((gamma - y_res)/n_para)**2
  END FUNCTION fg_eval
  
  FUNCTION BESSEL_K2(RMD_X) RESULT(k2)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: RMD_X
    REAL(dp) :: k2, RMD_Y, RMD_Y2, RMD_Y4, POLY
    
    IF (RMD_X < 2.0_dp) THEN
      RMD_Y = RMD_X
      RMD_Y2 = RMD_Y * RMD_Y
      RMD_Y4 = RMD_Y2 * RMD_Y2
      POLY = 2.0_dp/RMD_Y2      &
           - 0.5_dp         &
           + RMD_Y2/8.0_dp      &
           - RMD_Y4/96.0_dp     &
           + RMD_Y4*RMD_Y2/3072.0_dp
      k2 = POLY
           
    ELSE
      RMD_Y = 1.0_dp / RMD_X
      POLY = 1.0_dp      &
           + 0.5_dp         &
           + RMD_Y2/8.0_dp      &
           - RMD_Y4/96.0_dp     &
           + RMD_Y4*RMD_Y2/3072.0_dp
      k2 = SQRT(ACOS(-1.0_dp)/(2.0_dp*RMD_X)) * EXP(-RMD_X) * POLY 

    END IF   
  END FUNCTION BESSEL_K2
    
!***********************************************************************
!     save ray data
!***********************************************************************

  SUBROUTINE wr_save

    USE wrcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: NFL,IERR,NRAY,I,NSTP

      NFL=21
      CALL FWOPEN(NFL,KNAMWR,0,MODEFW,'WR',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WRSAVE: FWOPEN ERROR: IERR=',IERR
         RETURN
      END IF

      WRITE(NFL,ERR=9) NRAYMAX
      DO NRAY=1,NRAYMAX
         WRITE(NFL,ERR=9) NSTPMAX_NRAY(NRAY)
      END DO
      WRITE(6,*,ERR=9) 'NRAYMAX=',NRAYMAX
      DO NRAY=1,NRAYMAX
         WRITE(6,*,ERR=9) 'NSTPMAX=',NSTPMAX_NRAY(NRAY)
      END DO
      DO NRAY=1,NRAYMAX
         WRITE(NFL,ERR=9) (RAYIN(I,NRAY),I=1,8)
         WRITE(NFL,ERR=9) (CEXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (CEYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (CEZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RKXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RKYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RKZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RAYRB1(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RAYRB2(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         DO I=0,8
            WRITE(NFL,ERR=9) (RAYS(I,NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         ENDDO
         WRITE(NFL,ERR=9) (BNXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (BNYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (BNZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (BABSS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
!         DO NSTP=NSTPMAX(NRAY)-10,NSTPMAX_NRAY(NRAY)
!            WRITE(6,'(A,I5,1PE12.4)') 'NSTP,BABSS=',NSTP,BABSS(NSTP,NRAY)
!         END DO
      ENDDO
      CLOSE(NFL)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE: ',TRIM(KNAMWR)
      RETURN

    9 WRITE(6,*) 'XX WRLOAD: File IO error detected: KNAMFR= ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE wr_save

!***********************************************************************
!     load ray data
!***********************************************************************

  SUBROUTINE wr_load(NSTAT)

    USE wrcomm
    USE wrcalpwr,ONLY: wr_calc_pwr
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: NSTAT
    INTEGER:: NFL,IERR,I,NRAY,NSTP,NSTPMAX_temp
    REAL(rkind):: RF,RP,ZP,PHI,RNK,ANGP,ANGT,UU
    INTEGER,ALLOCATABLE:: NTEMP(:)

      NSTAT=0

      NFL=21
      CALL FROPEN(NFL,KNAMWR,0,MODEFR,'WR',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WRLOAX: FROPEN ERROR: IERR=',IERR
         RETURN
      END IF

      READ(NFL,END=8,ERR=9) NRAYMAX
      ALLOCATE(NTEMP(NRAYMAX))
      NSTPMAX_temp=0
      DO NRAY=1,NRAYMAX
         READ(NFL,END=8,ERR=9) NTEMP(NRAY)
         IF(NTEMP(NRAY).GT.NSTPMAX_temp) NSTPMAX_temp=NTEMP(NRAY)
      END DO
      WRITE(6,*) '## NRAYMAX,NSTPMAX=',NRAYMAX,NSTPMAX_temp
      IF(NSTPMAX.LT.NSTPMAX_temp) NSTPMAX=NSTPMAX_temp
      CALL wr_allocate
      DO NRAY=1,NRAYMAX
         NSTPMAX_NRAY(NRAY)=NTEMP(NRAY)
      END DO
      DEALLOCATE(NTEMP)
      DO NRAY=1,NRAYMAX
         READ(NFL,END=8,ERR=9) (RAYIN(I,NRAY),I=1,8)
         READ(NFL,END=8,ERR=9) (CEXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (CEYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (CEZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RKXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RKYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RKZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RAYRB1(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RAYRB2(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         DO I=0,8
            READ(NFL,END=8,ERR=9) (RAYS(I,NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         ENDDO
         READ(NFL,END=8,ERR=9) (BNXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (BNYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (BNZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (BABSS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
      ENDDO
      CLOSE(NFL)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE: ',TRIM(KNAMWR)
      IF(RAYRB1(1,1).EQ.0.D0.AND.RAYRB2(1,1).EQ.0.D0) THEN
         NSTAT=1
      ELSE
         NSTAT=2
      ENDIF
      RF=RAYIN(1,1)
      RP=RAYIN(2,1)
      ZP=RAYIN(3,1)
      PHI=RAYIN(4,1)
      RNK=RAYIN(5,1)
      ANGP=RAYIN(6,NRAY)
      ANGT=RAYIN(7,NRAY)
      UU=RAYIN(8,NRAY)

      CALL wr_calc_pwr

      RETURN

    8 WRITE(6,*) 'XX WRLOAD: End of file detected: KNAMFR= ',TRIM(KNAMWR)
      RETURN
    9 WRITE(6,*) 'XX WRLOAD: File IO error detected: KNAMFR= ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE wr_load

!***********************************************************************
!     write ray data as ascii format
!***********************************************************************

  SUBROUTINE wr_write

    USE wrcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: NFL,IERR,NRAY,I,NSTP

    NFL=22
    CALL FWOPEN(NFL,KNAMWR,1,MODEFW,'WR',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX wr_write: FWOPEN ERROR: IERR=',IERR
       RETURN
    END IF

    WRITE(NFL,'(I8)') NRAYMAX
    DO NRAY=1,NRAYMAX
       WRITE(NFL,'(8ES16.8)') (RAYIN(I,NRAY),I=1,8)
       WRITE(NFL,'(I8)') NSTPMAX_NRAY(NRAY)
       DO NSTP=1,NSTPMAX_NRAY(NRAY)
          WRITE(NFL,'(27ES16.8)') &
               (RAYS(I,NSTP,NRAY),I=1,8), &
               CEXS(NSTP,NRAY),CEYS(NSTP,NRAY),CEZS(NSTP,NRAY), &
               RKXS(NSTP,NRAY),RKYS(NSTP,NRAY),RKZS(NSTP,NRAY), &
               RXS(NSTP,NRAY),RYS(NSTP,NRAY),RZS(NSTP,NRAY), &
               RAYRB1(NSTP,NRAY),RAYRB2(NSTP,NRAY), &
               BNXS(NSTP,NRAY),BNYS(NSTP,NRAY),BNZS(NSTP,NRAY), &
               BABSS(NSTP,NRAY),RAYS(0,NSTP,NRAY)
       END DO
    END DO
    CLOSE(NFL)

    WRITE(6,*) '# DATA WAS SUCCESSFULLY WRITTEN TO THE FILE: ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE wr_write
  
  
  SUBROUTINE find_root_bisection(a, b, y_res, n_para, tol, root)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: a, b, y_res, n_para, tol
    REAL(dp), INTENT(OUT) :: root
    REAL(dp) :: fa, fb, fm, m, aloc, bloc
    INTEGER :: iter, max_iter
    
    aloc = a
    bloc = b
      
    fa = fg_eval(aloc, y_res, n_para)
    fb = fg_eval(bloc, y_res, n_para)
    IF (fa*fb > 0.0_dp) THEN 
      root = 0.5_dp*(a + b)
      RETURN
    END IF
      
    max_iter = 100
    DO iter = 1, max_iter
      m = 0.5_dp*(aloc + bloc)
      fm = fg_eval(m, y_res, n_para)
      IF (ABS(fm) < tol .OR. (bloc - aloc)*0.5_dp < tol) EXIT
      IF (fa*fm <= 0.0_dp) THEN
        bloc = m
        fb = fm
      ELSE
        aloc = m
        fa = fm
      END IF
    END DO
      
    root = 0.5_dp*(aloc + bloc)
  END SUBROUTINE find_root_bisection
  
   
  SUBROUTINE BESK0(X, BK0, IER)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: X
    REAL(dp), INTENT(OUT) :: BK0
    INTEGER, INTENT(OUT) :: IER
 
    INTEGER :: i, j
    REAL(dp) :: A(6), B(6), C(6)
    REAL(dp) :: W, WX, WX2, WX375, BI0

    A = [ -0.57721566_dp, 0.4227842_dp, 0.23069756_dp, &
          0.0348859_dp, 0.00262698_dp, 0.0001075_dp ]

    B = [ 1.0_dp, 3.5156229_dp, 3.0899424_dp, &
          1.2067492_dp, 0.2659732_dp, 0.360768_dp ]

    C = [ 1.25331414_dp, -0.07832358_dp, 0.02189568_dp, &
         -0.01062446_dp, 0.00587872_dp, -0.0025154_dp ]

    IF (X <= 0.0_dp) THEN
      IER = 999
      BK0 = 0.0_dp
      RETURN
    END IF
    
    IER = 0
    WX = X
    IF (WX <= 2.0_dp) THEN
      WX2 = (WX / 2.0_dp)**2
      WX375 = (WX / 3.75_dp)**2
      W = 0.0000074_dp
      BI0 = 0.0045813_dp
      DO i = 1, 6
        j = 7 - i
        W   = W * WX2 + A(j)
        BI0 = BI0 * WX375 + B(j)
      END DO
      BK0 = W - LOG(WX / 2.0_dp) * BI0
    ELSE
      WX2 = 2.0_dp / WX
      W = 0.00053208_dp
      DO i = 1, 6
        j = 7 - i
        W = W * WX2 + C(j)
      END DO
      BK0 = W / SQRT(WX) / EXP(WX)
    END IF

  END SUBROUTINE BESK0
  
  
  SUBROUTINE BESK1(X, BK1, IER)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: X
    REAL(dp), INTENT(OUT) :: BK1
    INTEGER, INTENT(OUT) :: IER
 
    INTEGER :: i, j
    REAL(dp) :: A(6), B(6), C(6)
    REAL(dp) :: W, WX, WX2, WX375, BI1

    A = [ 1.0_dp, 0.15443144_dp, -0.67278579_dp, &
         -0.18156897_dp, -0.01919402_dp, -0.00110404_dp ]

    B = [ 0.5_dp, 0.87890594_dp, 0.51498869_dp, &
          0.15084934_dp, 0.02658733_dp, 0.00301532_dp ]

    C = [ 1.25331414_dp, 0.23498619_dp, -0.0365562_dp, &
          0.01504268_dp, -0.00780353_dp, 0.00325614_dp ]

    IF (X <= 0.0_dp) THEN
      IER = 999
      BK1 = 0.0_dp
      RETURN
    END IF
    
    IER = 0
    WX = X
    IF (WX <= 2.0_dp) THEN
      WX2 = (WX / 2.0_dp)**2
      WX375 = (WX / 3.75_dp)**2
      W = 0.00004686_dp
      BI1 = 0.00032411_dp
      DO i = 1, 6
        j = 7 - i
        W   = W * WX2 + A(j)
        BI1 = BI1 * WX375 + B(j)
      END DO
      BK1 = W / WX + LOG(WX / 2.0_dp) * BI1 * WX
    ELSE
      WX2 = 2.0_dp / WX
      W = -0.00068245_dp
      DO i = 1, 6
        j = 7 - i
        W = W * WX2 + C(j)
      END DO
      BK1 = W / SQRT(WX) / EXP(WX)
    END IF

  END SUBROUTINE BESK1    
      
  
!駆動電流計算  
  SUBROUTINE eccd_write
    USE wrcomm
    USE libfio
    USE plcomm
    USE wrcalpwr
    USE plprof
    USE pllocal
    !USE bpsd_kinds,  ONLY : dp
    !USE bpsd_constants,  ONLY : PI, AEE, AME, RKBL, VC, EPS0
    IMPLICIT NONE
    
    INTEGER :: NFL, IERR, IER
    INTEGER :: NRAY, NSTP, i, num, n
    INTEGER :: num_roots
    
    TYPE(pl_prf_type), DIMENSION(NSMAX) :: plf
    TYPE(pl_mag_type) :: mag

    REAL(dp) :: x, y, z, rhon
    REAL(dp) :: RK, RT, BABSS_val, RN_val, RFIN_val
    REAL(dp) :: TkeV, T
    REAL(dp) :: omg_GHz, omg, omg_c, omg_pe
    REAL(dp) :: k_para, B_RESONANCE, ne, ne_center, ne_Vave, n_para
    REAL(dp) :: v_th, v, u, w, v_para, v_perp, v_perp_norm
    REAL(dp) :: cou_log, fg, fg_box
    REAL(dp) :: u_up, u_down, u_est
    REAL(dp) :: w_up, w_down, w_est, fmj
    REAL(dp) :: v_perp_up, v_perp_down, v_perp_est
    REAL(dp) :: gamma, dgamma, gamma_min, gamma_max
    REAL(dp) :: gamma_min_local, gamma_max_local
    REAL(dp) :: JPd_nor, nu0, JPd, eta, coef, Z_eff, fv, y_res, tol
    REAL(dp) :: PWR_val, EC_CURRENT, RMD_FACTOR, fv_sum
    REAL(dp) :: BK0, BK1, BK2, momentum, dmom, momentum_min, momentum_max
    REAL(dp) :: fgamma, BMIN, BMAX, alpha, B_BOUNCE, mu_t, Xi
    REAL(dp) :: coef_capture, JPd_capture, eta_capture, EC_CURRENT_CAPTURE
    REAL(dp) :: major_radius, delta_major_radius, EC_CURRENT_per_major_radius
    REAL(dp) :: EC_CURRENT_per_major_radius_local, TI, TOTAL_CURRENT
    REAL(dp) :: TOTAL_ETA, TOTAL_JPd
    REAL(dp) :: EC_CURRENT_CAP, EC_CURRENT_CAP_per_major_radius
    REAL(dp) :: EC_CURRENT_CAP_per_major_radius_local, TI_CAP, TOTAL_CURRENT_CAP
    REAL(dp) :: TOTAL_ETA_CAP, TOTAL_JPd_CAP
    REAL(dp), DIMENSION(2) :: a_root, b_root

    NFL=26
    CALL FWOPEN(NFL,KNAMWR,1,MODEFW,'WR',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX wr_write: FWOPEN ERROR: IERR=',IERR
       RETURN
    END IF
          

    OPEN(UNIT=4, FILE='eta.out', STATUS='UNKNOWN')
    
    num = 19000
    dgamma = 0.001_dp
    fv_sum = 0.0_dp
    DO NRAY=1,NRAYMAX
       TOTAL_CURRENT = 0.0_dp
       EC_CURRENT_per_major_radius_local = 0.0_dp
       DO NSTP=1,NSTPMAX_NRAY(NRAY)
          !WRITE(6, *) '@@@point000', NSTP, NSTPMAX_NRAY(NRAY)
          !温度と密度用
          x = RXS(NSTP, NRAY) !
          y = RYS(NSTP, NRAY) !
          z = RZS(NSTP, NRAY) !
          major_radius = SQRT(x**2 + y**2)
          CALL pl_mag(X,Y,Z,mag)
          rhon=mag%rhon
          CALL pl_prof(rhon,plf)  
          CALL wr_calc_pwr    
          CALL pl_bminmax(RHON, BMIN, BMAX)
          
          RK = RKXS(NSTP,NRAY) * BNXS(NSTP,NRAY) &
              +RKYS(NSTP,NRAY) * BNYS(NSTP,NRAY) &
              +RKZS(NSTP,NRAY) * BNZS(NSTP,NRAY)         
          RT = SQRT((plf(1)%RTPR)**2  + (plf(1)%RTPP)**2)
          BABSS_val = BABSS(NSTP,NRAY)
          RN_val = plf(1)%RN
          RFIN_val = RFIN(NRAY) ![MHz]
          PWR_val = pwr_nsa_nstp_nray(1,NSTP,NRAY) * 1.0D6
       
          !WRITE(NFL,'(I8,3ES16.8)') &
               !NSTP, RK, RFIN_val
               
          
          !CALL eccd_calc(RT, RK, BABSS_val, RN_val, RFIN_val)
          
          !pi = ATAN(1.0_dp) * 4.0_dp
          !q = 1.602176634_dp * 1.0D-19
          !kB = 1.380649_dp * 1.0D-23
          TkeV = plf(1)%RTPR
          T = TkeV * 1.0D3 * AEE / RKBL
          !me = 9.1093837139_dp * 1.0D-31

          omg_GHz = RFIN_val * 1.0D-3
          k_para = RK
          omg = omg_GHz * 1.0D9 * 2.0_dp * pi
          n = 1
          B_RESONANCE = BABSS_val
          ne = RN_val * 1.0D20
          ne_center = 0.9D20
          ne_Vave = ne_center / 1.15_dp
          !eps0 = 8.8541878176_dp * 1.0D-12
          !c = 2.99792458_dp * 1.0D8
          omg_c = AEE * B_RESONANCE / AME
          y_res = n * omg_c / omg
          omg_pe = SQRT(ne * AEE**2 / AME / EPS0)
          n_para = k_para * VC / omg
          IF (ne == 0.0_dp) THEN
            cou_log = 0.0_dp
          ELSE
            cou_log = 8.3_dp + 2.3_dp * LOG10(TkeV * 1.0D3 / SQRT(ne / 1.0D20))
          END IF
          v_th = VC * SQRT(1.0_dp - 1.0_dp / (RKBL * T / (AME * VC**2) + 1.0_dp)**2)
          
          IF (NSTP==1) THEN
            delta_major_radius = 1.0_dp
          ELSE
            delta_major_radius = ABS(SQRT(RXS(NSTP, NRAY)**2 + RYS(NSTP, NRAY)**2) &
                                - SQRT(RXS(NSTP-1, NRAY)**2 + RYS(NSTP-1, NRAY)**2))
          END IF
          
          !γの範囲
          fg_box = 0.0_dp
          num_roots = 0
          tol = 1.0e-6_dp
          
          DO i = 0, num
            gamma = 1.0_dp + dgamma * REAL(i, dp)
            fg = fg_eval(gamma, y_res, n_para)            
            IF (fg_box * fg < 0.0_dp) THEN
              num_roots = num_roots + 1
              IF (num_roots <= 2) THEN
                a_root(num_roots) = gamma - dgamma
                b_root(num_roots) = gamma
              END IF
            END IF
            fg_box = fg
          END DO
          
          IF (num_roots < 2) THEN
            gamma_min = 0.0_dp
            gamma_max = 0.0_dp
          ELSE 
            CALL find_root_bisection(a_root(1), b_root(1), y_res, n_para, tol, gamma_min)
            CALL find_root_bisection(a_root(2), b_root(2), y_res, n_para, tol, gamma_max)
          END IF
          
          IF (gamma_max == 0.0_dp) gamma_max = gamma_min
          
          gamma_min_local = gamma_min
          gamma_max_local = gamma_max
          
          gamma = gamma_min_local
          
          !共鳴しない場合
          IF (gamma == 0.0_dp) THEN
            JPd = 0.0_dp
            eta = 0.0_dp
            EC_CURRENT = 0.0_dp
            EC_CURRENT_per_major_radius = 0.0_dp
            JPd_capture = 0.0_dp
            eta_capture = 0.0_dp
            EC_CURRENT_CAPTURE = 0.0_dp
            w = 0.0_dp
            u = 0.0_dp
            WRITE(4,'(I8,1X,3ES16.8)') &
                 !NSTP, x, eta / 1.0D20, ABS(w), u
                 !NSTP, major_radius, eta/1.0D20, ABS(k_para)
                 !NSTP, x, y, z
                 NSTP, BABSS_val, BMAX, TkeV
            CYCLE
          END IF
          
          !<w>
          w_up = 0.0_dp
          w_down = 0.0_dp
          momentum_min = AME * VC * SQRT(gamma_min_local**2 - 1.0_dp)
          momentum_max = AME * VC * SQRT(gamma_max_local**2 - 1.0_dp)
          dmom = (momentum_max - momentum_min) / 1000
          momentum = momentum_min
          RMD_FACTOR = AME * VC**2 / (RKBL * T)
          !RMD_FACTOR = (VC / v_th)**2
          CALL BESK0(RMD_FACTOR, BK0, IER)
          CALL BESK1(RMD_FACTOR, BK1, IER)
          BK2 = (2.0_dp / RMD_FACTOR) * BK1 + BK0
          
          DO WHILE (momentum <= momentum_max)
            gamma = SQRT(1 + momentum**2/(AME * VC)**2)
            v_para = VC * ((gamma - y_res) / n_para) / gamma
            w = v_para / v_th
            !fv = 4.0_dp * pi * v**2 * (AME / (2.0_dp * pi * RKBL * T))**(1.5_dp) &
                 !* EXP(-AME * v**2 / (2.0_dp * RKBL * T))
            fv = RMD_FACTOR * EXP(-RMD_FACTOR * gamma) / &
                 (4.0_dp * PI * BK2)
            w_up = w_up + fv * w
            w_down = w_down + fv
           
            momentum = momentum + dmom
            
          END DO
          
          !<u>, <v_perp>
          gamma = gamma_min_local
          u_up = 0.0_dp
          u_down = 0.0_dp  
          v_perp_up = 0.0_dp
          v_perp_down = 0.0_dp   
          DO WHILE (gamma <= gamma_max_local)
            v = VC * SQRT(gamma**2 - 1.0_dp)
            !v_perp = VC * SQRT(gamma**2 - 1.0_dp - ((gamma - y_res)/n_para)**2)
            u = v / v_th
            !v_perp_norm = v_perp / v_th
            !WRITE(6, *) '@@@point2', NSTP
            fmj = (AME*VC)**3 * gamma * SQRT(gamma**2 - 1.0_dp) * RMD_FACTOR &
                     * EXP(-RMD_FACTOR * gamma) / BK2
            u_up = u_up + fmj * u
            u_down = u_down + fmj
            !v_perp_up = v_perp_up + fmj * v_perp_norm
            !v_perp_down = v_perp_down + fmj
           
            gamma = gamma + dgamma
          END DO
          
          !電流駆動効率計算
          u_est = u_up / u_down    
          u = u_est
          w_est = w_up / w_down    
          w = w_est
          !v_perp_est = v_perp_up / v_perp_down
          !v_perp_norm = v_perp_est
          Z_eff = 1.0_dp
          JPd_nor = u * w * 3.0_dp / (5.0_dp + Z_eff)          
          
          IF (cou_log == 0.0_dp) THEN
            JPd = 0.0_dp
          ELSE
            coef = -4.0_dp * EPS0**2 * 1.0D3 / (AEE**2 * 1.0D20)
            JPd = coef * TkeV / (RR * (ne / 1.0D20) * cou_log) * JPd_nor
          END IF 
          eta = RR * ne_Vave * JPd
          EC_CURRENT = JPd * PWR_val
          EC_CURRENT_per_major_radius = EC_CURRENT / delta_major_radius   
          TI = delta_major_radius / 2.0_dp &
               * (EC_CURRENT_per_major_radius+EC_CURRENT_per_major_radius_local)
          TOTAL_CURRENT = TOTAL_CURRENT + TI
          EC_CURRENT_per_major_radius_local = EC_CURRENT_per_major_radius
          
          !電子捕捉効果
          !B_BOUNCE = B_RESONANCE * (1.0_dp + (w / v_perp_norm)**2)
          B_BOUNCE = B_RESONANCE
          IF (B_BOUNCE / BMAX > 1.0_dp) THEN
            mu_t = 0.0_dp
          ELSE
            mu_t = (1.0_dp - B_BOUNCE / BMAX)**(0.5_dp)
          END IF
          !WRITE(6, *) '@@@point2', NSTP, B_BOUNCE, BMAX 
          Xi = u * mu_t / ABS(w)
          !coef_capture = - RKBL * T / (2.0_dp * PI * AEE**3 * ne * EXP(cou_log))
          !coef_capture = - AME * v_th**2 / (2.0_dp * PI * AEE**3 * ne * EXP(cou_log)) 
          alpha = (5.0_dp + Z_eff) / (1.0_dp + Z_eff)
          coef_capture = 1.0_dp - (Xi**(alpha)) * (1.0_dp + alpha/3.0_dp)
          !JPd_capture = (1 - Xi**alpha) * JPd &
                       !- coef_capture * Xi**alpha / &
                         !(2.0_dp * (1.0_dp + Z_eff)) * w * u
          JPd_capture = coef_capture * JPd
          eta_capture = RR * ne * JPd_capture
          EC_CURRENT_CAP = JPd_capture * PWR_val
          EC_CURRENT_CAP_per_major_radius = EC_CURRENT_CAP / delta_major_radius   
          TI_CAP = delta_major_radius / 2.0_dp &
               * (EC_CURRENT_CAP_per_major_radius+EC_CURRENT_CAP_per_major_radius_local)
          TOTAL_CURRENT_CAP = TOTAL_CURRENT_CAP + TI_CAP
          EC_CURRENT_CAP_per_major_radius_local = EC_CURRENT_CAP_per_major_radius
          WRITE(4,'(I8,1X,3ES16.8)') &
                !NSTP, major_radius, eta / 1.0D20, ABS(w), u
                !NSTP, major_radius, eta / 1.0D20, ABS(k_para)
                !NSTP, x, y, z, major_radius
                !NSTP, gamma_min, gamma_max
                NSTP, BABSS_val, BMAX, TkeV
                
                !NSTP, TkeV, omg_c, omg_pe
          !WRITE(NFL,'(I8,3ES16.8)') &
               !NSTP, u, w, v
       END DO
       TOTAL_JPd = TOTAL_CURRENT / 1.0D6
       TOTAL_ETA = TOTAL_CURRENT / 1.0D6 * RR * ne
       TOTAL_JPd_CAP = TOTAL_CURRENT_CAP / 1.0D6
       TOTAL_ETA_CAP = TOTAL_CURRENT_CAP / 1.0D6 * RR * ne
       WRITE(6, *) 'NRAY =', NRAY, 'ETA_cd =', TOTAL_JPd * 1.0D3, '[kA/MW]'
    END DO
    
    !分布確認
    !f0(p)
    OPEN(UNIT=7, FILE='fv.out', STATUS='UNKNOWN') 
    momentum_min = AME * VC * SQRT(gamma_min_local**2 - 1.0_dp)
    momentum_max = AME * VC * SQRT(gamma_max_local**2 - 1.0_dp)
    fv_sum = 0.0_dp
    momentum =  -momentum_max
    
    T = 25.0_dp * 1.0D3 * AEE / RKBL
    RMD_FACTOR = AME * VC**2 / (RKBL * T)
    !RMD_FACTOR = (VC / v_th)**2
    CALL BESK0(RMD_FACTOR, BK0, IER)
    CALL BESK1(RMD_FACTOR, BK1, IER)
    BK2 = (2.0_dp / RMD_FACTOR) * BK1 + BK0
    DO WHILE (momentum <= momentum_max)
      gamma = SQRT(1 + momentum**2/(AME * VC)**2)
      fv = RMD_FACTOR * EXP(-RMD_FACTOR * gamma) / (4.0_dp * PI * BK2)
      !fv = RMD_FACTOR * gamma**2 * SQRT(1 - 1/gamma**2) * EXP(-gamma*RMD_FACTOR) / BK2
      IF (i==0 .OR. i==num) THEN
        fv_sum = fv_sum + 0.5_dp * fv
      ELSE
        fv_sum = fv_sum + fv
      END IF 
      !WRITE(NFL,'(2ES16.8)') gamma, fv_sum
      WRITE(7,'(2ES16.8)') momentum*1.0D22, fv
      
      momentum = momentum + dmom
    END DO
    
    !f0(gamma)
    OPEN(UNIT=8, FILE='f_gamma.out', STATUS='UNKNOWN') 
    fv_sum = 0.0_dp
    gamma =1.0_dp
    
    T = 25.0_dp * 1.0D3 * AEE / RKBL
    RMD_FACTOR = AME * VC**2 / (RKBL * T)
    !RMD_FACTOR = (VC / v_th)**2
    CALL BESK0(RMD_FACTOR, BK0, IER)
    CALL BESK1(RMD_FACTOR, BK1, IER)
    BK2 = (2.0_dp / RMD_FACTOR) * BK1 + BK0
    DO WHILE (gamma <= 2.0_dp)
      fgamma = gamma * SQRT(gamma**2 - 1) * RMD_FACTOR &
           * EXP(-RMD_FACTOR * gamma) / BK2
      !fv = RMD_FACTOR * gamma**2 * SQRT(1 - 1/gamma**2) * EXP(-gamma*RMD_FACTOR) / BK2
      IF (i==0 .OR. i==num) THEN
        fv_sum = fv_sum + 0.5_dp * fv
      ELSE
        fv_sum = fv_sum + fv
      END IF 
      !WRITE(NFL,'(2ES16.8)') gamma, fv_sum
      WRITE(8,'(2ES16.8)') gamma, fgamma
      
      gamma = gamma + dgamma
    END DO
    
    
    CALL BESK1(10.0_dp, BK1, IER)
    fv_sum = fv_sum * dgamma
    WRITE(NFL,'(2ES16.8)') fv_sum, BK1 
    CLOSE(7)
    CLOSE(8)
 
    CLOSE(NFL)
    CLOSE(4)
    WRITE(6,*) '# DATA WAS SUCCESSFULLY WRITTEN TO THE FILE: ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE eccd_write
  

END MODULE wrfile
