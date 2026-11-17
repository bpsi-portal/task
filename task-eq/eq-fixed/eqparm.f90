! eqparm.f90

MODULE eqparm

  PRIVATE
  PUBLIC eq_parm
  
CONTAINS
  
!     ****** INPUT PARAMETERS ******

  SUBROUTINE eq_parm(MODE,KIN,IERR)

!     MODE=0 : standard namelist input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

    USE libkio
    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    CHARACTER(LEN=*),INTENT(IN)::  kin
    INTEGER,INTENT(OUT):: ierr
    
1   CONTINUE
    CALL TASK_PARM(mode,'EQ',kin,eq_nlin,eq_plst,ierr)
    IF(IERR.NE.0) RETURN

    CALL eq_check(ierr)
    IF(mode.EQ.0.AND.ierr.NE.0) GOTO 1
    IF(ierr.NE.0) ierr=ierr+100

    RETURN
  END SUBROUTINE eq_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE eq_nlin(nid,ist,ierr)

    USE eqcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nid
    INTEGER,INTENT(OUT):: ist,ierr
    
    NAMELIST /EQ/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,
                  RHOMIN,QMIN,MODELG,MODELQ,RHOITB, &
                  IDEBUG,KNAMEQ,KNAMEQ2, &
                  MODEFR,MODEFW, &
                  PP0,PP1,PP2,PROFP0,PROFP1,PROFP2, &
                  FF0,FF1,FF2,PROFF0,PROFF1,PROFF2, &
                  PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2, &
                  PT0,PT1,PT2,PROFTP0,PROFTP1,PROFTP2,PTSEQ, &
                  PV0,PV1,PV2,PROFV0,PROFV1,PROFV2,PN0EQ, &
                  PROFR0,PROFR1,PROFR2, &
                  RGMIN,RGMAX,ZGMIN,ZGMAX, &
                  ZLIMP,ZLIMM, &
                  PSIB,RIPFC,RPFC,ZPFC,WPFC,FRBIN, &
                  EPSEQ,EPSNW,DELNW, &
                  NSGMAX,NTGMAX,NUGMAX, &
                  NRGMAX,NZGMAX, &
                  NPSMAX,NRVMAX,NTVMAX, &
                  NRMAX,NTHMAX,NSUMAX, &
                  MDLEQF,MDLEQC,MDLEQA,MDLEQX,MDLEQV,NPRINT, &
                  NLPMAX,NPFCMAX,NLPNW&
                  
    READ(nid,EQ,IOSTAT=IST,ERR=9800,END=9900)
                  
    RBRA=RB/RA
    ierr=0
    RETURN

9800 IERR=8
    WRITE(6,*) 'XX eq_parm READ ERROR: IST=',IST
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE eq_nlin

  !     ***** INPUT PARAMETER LIST *****

  SUBROUTINE eq_plst

    IMPLICIT NONE

    WRITE(6,*) '&EQ : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP'
    WRITE(6,*) 'RHOMIN,QMIN,MODELG,MODELQ,RHOITB,'
    WRITE(6,*) 'IDEBUG,KNAMEQ,KNAMEQ2,'
    WRITE(6,*) 'PP0,PP1,PP2,PROFP0,PROFP1,PROFP2'
    WRITE(6,*) 'FF0,FF1,FF2,PROFF0,PROFF1,PROFF2'
    WRITE(6,*) 'PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2'
    WRITE(6,*) 'PT0,PT1,PT2,PROFTP0,PROFTP1,PROFTP2,PTSEQ'
    WRITE(6,*) 'PV0,PV1,PV2,PROFV0,PROFV1,PROFV2,PN0EQ,HM'
    WRITE(6,*) 'PROFR0,PROFR1,PROFR2'
    WRITE(6,*) 'NSGMAX,NTGMAX,NUGMAX,NRGMAX,NZGMAX,NPSMAX'
    WRITE(6,*) 'NRMAX,NTHMAX,NSUMAX,NRVMAX,NTVMAX'
    WRITE(6,*) 'MDLEQF,MDLEQC,MDLEQA,MDLEQX,MDLEQV,NPRINT'
    WRITE(6,*) 'EPSEQ,NLPMAX,EPSNW,DELNW,NLPNW'
    WRITE(6,*) 'RGMIN,RGMAX,RZMIN,RZMAX,ZLIMP,ZLIMM'
    WRITE(6,*) 'MODEFR,MODEFW,'
    WRITE(6,*) 'PSIB,NPFCMAX,RIPFC,RPFC,ZPFC,WPFC,FRBIN'
    RETURN
  END SUBROUTINE eq_plst

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE eq_check(ierr)

    USE eqcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    ierr=0

    IF(RB.LT.RA) THEN
       WRITE(6,'(A,1P2E12.4)') '!! RB.LT.RA: set RB=RA: RA,RB=',RA,RB
       RB=RA
    ENDIF

    RETURN
  END SUBROUTINE eq_check

END MODULE eqparm
