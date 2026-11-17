! eqview.f90

MODULE eqview

CONTAINS
  
!     ****** SHOW PARAMETERS ******

  SUBROUTINE eq_view

    USE eqcomm_parm
    IMPLICIT NONE

    WRITE(6,601) 'RR    ',RR,     'RA    ',RA, &
                 'RKAP  ',RKAP,   'RDLT  ',RDLT
    WRITE(6,601) 'BB    ',BB,     'RIP   ',RIP, &
                 'Q0    ',Q0,     'QA    ',QA
    WRITE(6,601) 'RHOMIN',RHOMIN, 'QMIN  ',QMIN, &
                 'RB    ',RB
    WRITE(6,601) 'RGMIN ',RGMIN,  'RGMAX ',RGMAX, &
                 'ZGMIN ',ZGMIN,  'ZGMAX ',ZGMAX
    WRITE(6,601) 'ZLIMP ',ZLIMP,  'ZLIMM ',ZLIMM, &
                 'FRBIN ',FRBIN
    WRITE(6,601) 'PP0   ',PP0,    'PROFP0',PROFP0, &
                 'PJ0   ',PJ0,    'PROFJ0',PROFJ0
    WRITE(6,601) 'PP1   ',PP1,    'PROFP1',PROFP1, &
                 'PJ1   ',PJ1,    'PROFJ1',PROFJ1
    WRITE(6,601) 'FF0   ',FF0,    'PROFF0',PROFF0
    WRITE(6,601) 'FF1   ',FF1,    'PROFF1',PROFF1, &
                 'FF2   ',FF2,    'PROFF2',PROFF2
    WRITE(6,601) 'PT0   ',PT0,    'PROFT0',PROFTP0, &
                 'PV0   ',PV0,    'PROFV0',PROFV0
    WRITE(6,601) 'PT1   ',PT1,    'PROFT1',PROFTP1, &
                 'PV1   ',PV1,    'PROFV1',PROFV1
    WRITE(6,601) 'PT2   ',PT2,    'PROFT2',PROFTP2, &
                 'PV2   ',PV2,    'PROFV2',PROFV2
    WRITE(6,601) 'PTSEQ ',PTSEQ,  'PN0EQ ',PN0EQ
    WRITE(6,601) 'PROFR0',PROFR0, 'PROFR1',PROFR1, &
                 'PROFR2',PROFR2
    WRITE(6,601) 'EPSEQ ',EPSEQ,  'EPSNW ',EPSNW, &
                 'DELNW ',DELNW
    IF(MDLEQF.GE.10.AND.MDLEQF.LT.20) THEN
       WRITE(6,601) 'PSIB:0',PSIB(0),  'PSIB:1',PSIB(1), &
                    'PSIB:2',PSIB(2),  'PSIB:3',PSIB(3)
       WRITE(6,604) &
            (NPFC,RIPFC(NPFC),RPFC(NPFC),ZPFC(NPFC),WPFC(NPFC), &
             NPFC=1,NPFCMAX)
    ENDIF
    WRITE(6,602) 'NSGMAX',NSGMAX, 'NTGMAX',NTGMAX, &
                 'NUGMAX',NUGMAX, 'MODELG',MODELG
    WRITE(6,602) 'NRGMAX',NRGMAX, 'NZGMAX',NZGMAX, &
                 'NRVMAX',NRVMAX, 'NTVMAX',NTVMAX
    WRITE(6,602) 'NRMAX ',NRMAX,  'NTHMAX',NTHMAX, &
                 'NPSMAX',NPSMAX, 'NSUMAX',NSUMAX
    WRITE(6,602) 'MDLEQF',MDLEQF, 'MDLEQC',MDLEQC, &
                 'MDLEQA',MDLEQA, 'MODELQ',MODELQ
    WRITE(6,602) 'MDLEQX',MDLEQX, 'MDLEQV',MDLEQV, &
                 'MODEFR',MODEFR, 'MDDEFW',MODEFW
    WRITE(6,602) 'NPRINT',NPRINT, 'NLPMAX',NLPMAX, &
                 'NLPNW ',NLPNW
    RETURN
601 FORMAT(4(A6,'=',1PE11.2:2X))
602 FORMAT(4(A6,'=',I7:6X))
604 FORMAT(' NPFC ','RIPFC ',6X,'RPFC  ',6X,'ZPFC  ',6X,'WPFC'/ &
         (I6,1P4E12.4))
  END SUBROUTINE eq_VIEW
END MODULE eqview
