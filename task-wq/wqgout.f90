! wqgout.f90

MODULE wqgout

  PRIVATE
  PUBLIC wq_gout

CONTAINS

  ! *** graphic menu ***
  
  SUBROUTINE wq_gout
    USE wqcomm
    USE wqparm,ONLY: wq_parm,wq_broadcast
    USE wqview,ONLY: wq_view
    USE libchar
    USE libkio
    USE libmpi
    IMPLICIT NONE
    CHARACTER(LEN=80):: line
    CHARACTER(LEN=1):: kid
    INTEGER:: mode

1   CONTINUE
    
    IF(nrank.EQ.0) THEN
       WRITE(6,'(A)') &
            '## wqgout menu: g:global p:profile t:time  parm  x/end'
       CALL TASK_KLIN(line,kid,mode,wq_parm)
       CALL ToUpper(kid)
    END IF
    CALL mtx_broadcast1_character(kid)
    CALL mtx_broadcast1_integer(mode)

    IF(mode.EQ.2) CALL wq_broadcast
    IF(mode.NE.1) GO TO 1

    SELECT CASE(kid)
    CASE('G')
       CALL wq_gsub_global
    CASE('P')
       CALL wq_gsub_profile
    CASE('T')
       CALL wq_gsub_time
    CASE('X')
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wq_gout

  ! *** time history of global quantities ***
  
  SUBROUTINE wq_gsub_global

    USE wqcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER:: ngt
    REAL(rkind):: tn(ngt_m)

    DO ngt=1,ngt_max
       tn(ngt)=t_ngt(ngt)/period
    END DO

    CALL PAGES
    WRITE(6,'(A,I8,4ES12.4)') 'ngt_max,t,E:', &
         ngt_max,t_ngt(ngt_max), &
         EX_max_ngt(ngt_max),EY_max_ngt(ngt_max),EZ_max_ngt(ngt_max)
    CALL grd1d(1,tn,EX_max_ngt,ngt_m,ngt_max,1,'@Ex_max vs tn@',FMIN=0.D0)
    CALL grd1d(2,tn,EY_max_ngt,ngt_m,ngt_max,1,'@Ey_max vs tn@',FMIN=0.D0)
    CALL grd1d(3,tn,EZ_max_ngt,ngt_m,ngt_max,1,'@Ez_max vs tn@',FMIN=0.D0)
    CALL grd1d(4,tn,pabs_tot_ngt,ngt_m,ngt_max,1,'@Pabs_tot vs tn@',FMIN=0.D0)

    CALL PAGEE
    RETURN
  END SUBROUTINE wq_gsub_global

  ! *** profile of elecric field and pabs ***
  
  SUBROUTINE wq_gsub_profile
  
    USE wqcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER:: nx,ny
    REAL(rkind):: fr(nxmax,nymax),fi(nxmax,nymax),fa(nxmax,nymax)
    REAL(rkind):: line_thickness(1)
    
    line_thickness(1)=0.035D0

    CALL PAGES

    DO ny=1,nymax
       DO nx=1,nxmax
          fr(nx,ny)=REAL(EX(nx,ny))
          fi(nx,ny)=AIMAG(EX(nx,ny))
          fa(nx,ny)=ABS(EX(nx,ny))
       END DO
    END DO
    CALL grd2d_wq(5,xn_nx,yn_ny,fr,nxmax,nxmax,nymax,'@Real Ex@',ASPECT=0.D0)
    CALL grd2d_wq(6,xn_nx,yn_ny,fi,nxmax,nxmax,nymax,'@Imag Ex@',ASPECT=0.D0)
    CALL grd2d_wq(7,xn_nx,yn_ny,fa,nxmax,nxmax,nymax,'@Abs  Ex@',ASPECT=0.D0)
    
    DO ny=1,nymax
       DO nx=1,nxmax
          fr(nx,ny)=REAL(EY(nx,ny))
          fi(nx,ny)=AIMAG(EY(nx,ny))
          fa(nx,ny)=ABS(EY(nx,ny))
       END DO
    END DO
    CALL grd2d_wq( 8,xn_nx,yn_ny,fr,nxmax,nxmax,nymax,'@Real Ey@',ASPECT=0.D0)
    CALL grd2d_wq( 9,xn_nx,yn_ny,fi,nxmax,nxmax,nymax,'@Imag Ey@',ASPECT=0.D0)
    CALL grd2d_wq(10,xn_nx,yn_ny,fa,nxmax,nxmax,nymax,'@Abs  Ey@',ASPECT=0.D0)
    
    DO ny=1,nymax
       DO nx=1,nxmax
          fr(nx,ny)=REAL(EZ(nx,ny))
          fi(nx,ny)=AIMAG(EZ(nx,ny))
          fa(nx,ny)=ABS(EZ(nx,ny))
       END DO
    END DO
    CALL grd2d_wq(11,xn_nx,yn_ny,fr,nxmax,nxmax,nymax,'@Real Ez@',ASPECT=0.D0)
    CALL grd2d_wq(12,xn_nx,yn_ny,fi,nxmax,nxmax,nymax,'@Imag Ez@',ASPECT=0.D0)
    CALL grd2d_wq(13,xn_nx,yn_ny,fa,nxmax,nxmax,nymax,'@Abs  Ez@',ASPECT=0.D0)

    CALL PAGEE

    IF(ngt_max.EQ.0) RETURN
    IF(pabs_tot_ngt(ngt_max).LE.0.D0) RETURN

    CALL PAGES
    
    DO ny=1,nymax
       DO nx=1,nxmax
          fr(nx,ny)=REAL(pabs(nx,ny))
       END DO
    END DO
    CALL grd2d_wq(5,xn_nx,yn_ny,fr,nxmax,nxmax,nymax,'@Pabs@',ASPECT=0.D0)

    CALL PAGEE
    RETURN
  END SUBROUTINE wq_gsub_profile

  ! *** time history of elecric field and pabs ***
  
  SUBROUTINE wq_gsub_time
  
    USE wqcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER:: id,ngr,nx,ny,npage,npage_max,ng_id,ng
    CHARACTER(LEN=9):: title
    CHARACTER(LEN=16):: title_t
    REAL(rkind):: f(nxmax,nymax)
    REAL(rkind):: line_thickness(1)
    
    line_thickness(1)=0.035D0

1   CONTINUE
    WRITE(6,'(A)') '## Input : 1:Exr 2:Exi 3:Exa 4:Eyr 5:Eyi 6:Eya'
    WRITE(6,'(A)') '           7:Ezr 8:Ezi 9:Eza 10:Pabs 11:Ea 0:end'
    READ(5,*,END=9000,ERR=1) id
    IF(id.EQ.0) GO TO 9000

    SELECT CASE(id)
    CASE(1)
       title='@Real Ex:'
    CASE(2)
       title='@Imag Ex:'
    CASE(3)
       title='@Abs  Ex:'
    CASE(4)
       title='@Real Ey:'
    CASE(5)
       title='@Imag Ey:'
    CASE(6)
       title='@Abs  Ey:'
    CASE(7)
       title='@Real Ez:'
    CASE(8)
       title='@Imag Ez:'
    CASE(9)
       title='@Abs  Ez:'
    CASE(10)
       title='@Pabs   :'
    CASE(11)
       title='@Eabs   :'
    END SELECT

    npage_max=(ngr_max-1)/16+1
    IF(npage_max.EQ.1) THEN
       SELECT CASE(ngr_max)
       CASE(1)
          ng_id=0
       CASE(2:4)
          ng_id=1
       CASE(5:9)
          ng_id=5
       CASE(10:16)
          ng_id=14
       END SELECT
    ELSE
       ng_id=14
    END IF

    DO npage=1,npage_max
       CALL PAGES
       DO ng=1,MIN(ngr_max-16*(npage-1),16)
          ngr=16*(npage-1)+ng
          WRITE(title_t,'(A,ES12.4,A)') ' t=',t_ngr(ngr),'@'
          DO ny=1,nymax
             DO nx=1,nxmax
                SELECT CASE(id)
                CASE(1)
                   f(nx,ny)=REAL(EX_save(nx,ny,ngr))
                CASE(2)
                   f(nx,ny)=AIMAG(EX_save(nx,ny,ngr))
                CASE(3)
                   f(nx,ny)=ABS(EX_save(nx,ny,ngr))
                CASE(4)
                   f(nx,ny)=REAL(EY_save(nx,ny,ngr))
                CASE(5)
                   f(nx,ny)=AIMAG(EY_save(nx,ny,ngr))
                CASE(6)
                   f(nx,ny)=ABS(EY_save(nx,ny,ngr))
                CASE(7)
                   f(nx,ny)=REAL(EZ_save(nx,ny,ngr))
                CASE(8)
                   f(nx,ny)=AIMAG(EZ_save(nx,ny,ngr))
                CASE(9)
                   f(nx,ny)=ABS(EZ_save(nx,ny,ngr))
                CASE(10)
                   f(nx,ny)=pabs_save(nx,ny,ngr)
                CASE(11)
                   f(nx,ny)=SQRT(ABS(EX_save(nx,ny,ngr))**2 &
                                +ABS(EY_save(nx,ny,ngr))**2 &
                                +ABS(EZ_save(nx,ny,ngr))**2)
                END SELECT
             END DO ! nx
          END DO !ny

          CALL grd2d_wq(ng-1+ng_id,xn_nx,yn_ny,f,nxmax,nxmax,nymax, &
               title//title_t,ASPECT=0.D0,LINE_THICKNESS=line_thickness)
       END DO
       CALL PAGEE
    END DO
    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wq_gsub_time

  SUBROUTINE grd2d_wq(ngp,fx,fy,ff,nxm,nxmax,nymax,title,aspect,line_thickness)

    USE wqcomm_parm
    USE libgrf
    IMPLICIT NONE

    INTEGER,INTENT(IN):: NGP                       ! Graph size and postion ID
    INTEGER,INTENT(IN):: NXM,NXMAX,NYMAX           ! Number of data  
    REAL(dp),INTENT(IN):: FX(NXMAX),FY(NYMAX)       ! 1D data
    REAL(dp),INTENT(IN):: FF(NXM,NYMAX)             ! 1D data
    REAL(dp),INTENT(IN),OPTIONAL:: ASPECT
    REAL(dp),INTENT(IN),OPTIONAL:: LINE_THICKNESS(1)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: TITLE

    TYPE(grf_attr_type):: A

    INTEGER:: NX,NY
    REAL,DIMENSION(:),ALLOCATABLE:: GX,GY
    REAL,DIMENSION(:,:),ALLOCATABLE:: GF
    INTEGER:: medium
    REAL:: x1,x2,y1,y2

    INTERFACE
       SUBROUTINE SETLNW(VALUE)
         REAL,INTENT(IN):: VALUE
       END SUBROUTINE SETLNW
       SUBROUTINE SETRGB(R,G,B)
         REAL,INTENT(IN):: R,G,B
       END SUBROUTINE SETRGB
       FUNCTION GUCLIP(VALUE)
         USE task_kinds,ONLY: dp
         REAL(dp),INTENT(IN):: VALUE
         REAL:: GUCLIP
       END FUNCTION GUCLIP
    END INTERFACE
 
    ALLOCATE(GX(NXMAX),GY(NYMAX),GF(NXM,NYMAX))
    A%MODE_2D=1

    DO NX=1,NXMAX
       GX(NX)=GUCLIP(FX(NX))
    ENDDO
    DO NY=1,NYMAX
       GY(NY)=GUCLIP(FY(NY))
    END DO
    DO NY=1,NYMAX
       DO NX=1,NXMAX
          GF(NX,NY)=GUCLIP(FF(NX,NY))
       ENDDO
    ENDDO

    CALL GRD_CONVERT(NGP,GX,GY,GF,NXM,NXMAX,NYMAX, &
         TITLE=TITLE,ASPECT=ASPECT,LINE_THICkNESS=LINE_THICKNESS, &
         A=A)

    CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                A%XMIN,A%XMAX,A%YMIN,A%YMAX)
    CALL SETCLP(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX)

    CALL SETRGB(0.0,1.0,0.0)
    CALL SETLNW(A%FRAME_THICKNESS)
    
    DO medium=1,medium_max
       x1=xnmin_medium(medium)
       x2=xnmax_medium(medium)
       y1=ynmin_medium(medium)
       y2=ynmax_medium(medium)
       CALL MOVE2D(x1,y1)
       CALL DRAW2D(x2,y1)
       CALL DRAW2D(x2,y2)
       CALL DRAW2D(x1,y2)
       CALL DRAW2D(x1,y1)
    END DO
    
    CALL grd2d(ngp,fx,fy,ff,nxm,nxmax,nymax, &
               TITLE=title,ASPECT=aspect,LINE_THICKNESS=line_thickness)

  END SUBROUTINE grd2d_wq
END MODULE wqgout
