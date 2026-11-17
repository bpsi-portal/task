! eqinit.f90

MODULE eqinit

  PRiVATE
  PUBLIC eq_init

CONTAINS

!  ****** INITIALIZE CONSTATNS AND DEFAULT VALUES ******

  SUBROUTINE eq_init

    USE eqcomm
    IMPLICIT NONE

!        RIP   : Plasma current                                 (MA)
!        RBRA  : default ratio of RB/RA
!        FRBIN : (RB_inside-RA)/(RB_outside-RA)

      RIP   = 3.D0
      RBRA  = RB/RA
      FRBIN = 1.D0

!     ======( MODEL PARAMETERS )======

!        MODELG: Control plasma geometry model
!                   0: Slab geometry
!                   1: Cylindrical geometry
!                   2: Toroidal geometry
!                   3: TASK/EQ output geometry
!                   4: VMEC output geometry
!                   5: EQDSK output geometry
!                   6: Boozer output geometry
!        MODELN: Control plasma profile
!                   0: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; PN=0 in SOL
!                   1: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; PNS in SOL
!                   7: Read from file by means of WMDPRF routine (DIII-D)
!                   8: Read from file by means of WMXPRF routine (JT-60)
!                   9: Read from file KNAMTR (TASK/TR)
!        MODELQ: Control safety factor profile (for MODELG=0,1,2)
!                   0: Parabolic q profile (Q0,QA,RHOMIN,RHOITB)
!                   1: Given current profile (RIP,PROFJ0,PROFJ1,PROFJ2)

      MODELG= 2  ! @plcomm
      MODELN= 0  ! model_prof @plcomm
      MODELQ= 0  ! @plcomm

!     ======( IO file names )======

!        KNAMEQ2:Filename of addisional equilibrium data

      KNAMEQ2= 'eqdata2'

      NRMAXPL= 100
      NSMAXPL= NSMAX

      IDEBUG = 0

!     *** PROFILE PARAMETERS ***

!        PP0   : Plasma pressure (main component)              (MPa)
!        PP1   : Plasma pressure (sub component)               (MPa)
!        PP2   : Plasma pressure (increment within ITB)        (MPa)
!        PROFP0: Pressure profile parameter
!        PROFP1: Pressure profile parameter
!        PROFP2: Pressure profile parameter

!        PPSI=PP0*(1.D0-PSIN**PROFR0)**PROFP0
!    &       +PP1*(1.D0-PSIN**PROFR1)**PROFP1
!    &       +PP2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFP2

!        The third term exists for RHO < RHOITB

      PP0    = 0.001D0
      PP1    = 0.0D0
      PP2    = 0.0D0
      PROFP0 = 1.5D0
      PROFP1 = 1.5D0
      PROFP2 = 2.0D0

!        PJ0   : Current density at R=RR (main component) : Fixed to 1
!        PJ1   : Current density at R=RR (sub component)       (arb)
!        PJ2   : Current density at R=RR (sub component)       (arb)
!        PROFJ0: Current density profile parameter
!        PROFJ1: Current density profile parameter
!        PROFJ2: Current density profile parameter

!      HJPSI=-PJ0*(1.D0-PSIN**PROFR0)**PROFJ0
!     &                *PSIN**(PROFR0-1.D0)
!     &      -PJ1*(1.D0-PSIN**PROFR1)**PROFJ1
!     &                *PSIN**(PROFR1-1.D0)
!     &      -PJ2*(1.D0-PSIN**PROFR2)**PROFJ2
!     &                *PSIN**(PROFR2-1.D0)

!        The third term exists for RHO < RHOITB

      PJ0    = 1.00D0
      PJ1    = 0.0D0
      PJ2    = 0.0D0
      PROFJ0 = 1.5D0
      PROFJ1 = 1.5D0
      PROFJ2 = 1.5D0

!        FF0   : Current density at R=RR (main component) : Fixed to 1
!        FF1   : Current density at R=RR (sub component)       (arb)
!        FF2   : Current density at R=RR (sub component)       (arb)
!        PROFF0: Current density profile parameter
!        PROFF1: Current density profile parameter
!        PROFF2: Current density profile parameter

!      FPSI=BB*RR
!     &      +FF0*(1.D0-PSIN**PROFR0)**PROFF0
!     &      +FF1*(1.D0-PSIN**PROFR1)**PROFF1
!     &      +FF2*(1.D0-PSIN**PROFR2)**PROFF2

!        The third term exists for RHO < RHOITB

      FF0    = 1.0D0
      FF1    = 0.0D0
      FF2    = 0.0D0
      PROFF0 = 1.5D0
      PROFF1 = 1.5D0
      PROFF2 = 1.5D0

!        PT0   : Plasma temperature (main component)           (keV)
!        PT1   : Plasma temperature (sub component)            (keV)
!        PT2   : Plasma temperature (increment within ITB)     (keV)
!        PTSEQ : Plasma temperature (at surface)               (keV)
!        PROFTP0: Temperature profile parameter
!        PROFTP1: Temperature profile parameter
!        PROFTP2: Temperature profile parameter

!        TPSI=PTSEQ+(PT0-PTSEQ)*(1.D0-PSIN**PROFR0)**PROFTP0
!    &       +PT1*(1.D0-PSIN**PROFR1)**PROFTP1
!    &       +PT2*(1.D0-PSIN/PSIITB)**PROFR2)**PROFTP2
!    &       +PTSEQ

!        The third term exits for RHO < RHOITB

      PT0    = 1.0D0
      PT1    = 0.0D0
      PT2    = 0.0D0
      PTSEQ  = 0.05D0

      PROFTP0 = 1.5D0
      PROFTP1 = 1.5D0
      PROFTP2 = 2.0D0
!----
!        PV0   : Toroidal rotation (main component)              (m/s)
!        PV1   : Toroidal rotation (sub component)               (m/s)
!        PV2   : Toroidal rotation (increment within ITB)        (m/s)
!        PROFV0: Velocity profile parameter
!        PROFV1: Velocity profile parameter
!        PROFV2: Velocity profile parameter

!        PVSI=PV0*(1.D0-PSIN**PROFR0)**PROFV0
!    &       +PV1*(1.D0-PSIN**PROFR1)**PROFV1
!    &       +PV2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFV2

!        The third term exits for RHO < RHOITB

      PV0    = 0.0D0
      PV1    = 0.0D0
      PV2    = 0.0D0
      PROFV0 = 1.5D0
      PROFV1 = 1.5D0
      PROFV2 = 2.0D0

!        PN0EQ : Plasma number density(constant)

      PN0EQ  = 1.D20

!        PROFR0: Profile parameter
!        PROFR1: Profile parameter
!        PROFR2: Profile parameter
!        RHOITB: Normalized radius SQRT(PSI/PSIA) at ITB

      PROFR0 = 1.D0
      PROFR1 = 2.D0
      PROFR2 = 2.D0

!        OTC   : Constant OMEGA**2/TPSI
!        HM    : Constant                                       (Am)

      OTC = 0.15D0
      HM  = 1.D6

!     *** MESH PARAMETERS ***

!        NSGMAX: Number of radial mesh points for Grad-Shafranov eq.
!        NTGMAX: Number of poloidal mesh points for Grad-Shafranov eq.
!        NUGMAX: Number of radial mesh points for flux-average quantities
!        NRGMAX: Number of horizontal mesh points in R-Z plane
!        NZGMAX: Number of vertical mesh points in R-Z plane
!        NPSMAX: Number of flux surfaces
!        NRMAX : Number of radial mesh points for flux coordinates
!        NTHMAX: Number of poloidal mesh points for flux coordinates
!        NSUMAX: Number of boundary points
!        NRVMAX: Number of radial mesh of surface average
!        NTVMAX: Number of poloidal mesh for surface average

      NSGMAX = 32
      NTGMAX = 32
      NUGMAX = 32

      NRGMAX = 33
      NZGMAX = 33
      NPSMAX = 21

      NRMAX  = 50
      NTHMAX = 64
      NSUMAX = 65

      NRVMAX = 50
      NTVMAX = 400

!     *** CONTROL PARAMETERS ***

!        EPSEQ  : Convergence criterion for equilibrium
!        NLPMAX : Maximum iteration number of EQ
!        EPSNW  : Convergence criterion for newton method
!        DELNW  : Increment for derivative in newton method
!        NLPNW  : Maximum iteration number in newton method

      EPSEQ  = 1.D-6
      NLPMAX = 100
      EPSNW  = 1.D-2
      DELNW  = 1.D-2
      NLPNW  = 20

!        MDLEQF : Profile parameter
!            0: given analytic profile  P,J_tor,T,Vph + Ip
!            1: given analytic profile  P,F           + Ip
!            2: given analytic profile  P,J_para      + Ip
!            3: given analytic profile  P,J_para
!            4: given analytic profile  P,q
!            5: given spline profile    P,J_tor,T,Vph + Ip
!            6: given spline profile    P,F           + Ip
!            7: given spline profile    P,J_para      + Ip
!            8: given spline profile    P,J_para
!            9: given spline profile    P,q

      MDLEQF = 0

!        MDLEQA : Rho in P(rho), F(rho), q(rho),...
!            0: SQRT(PSIP/PSIPA)
!            1: SQRT(PSIT/PSITA)

      MDLEQA = 0

!        MDLEQC : Poloidal coordinate parameter
!            0: Poloidal length coordinate
!            1: Boozer coordinate

      MDLEQC = 0

!        MDLEQX : Free boundary calculation
!            0: Given PSIB and RIPFC
!            1: PSIB adjusted after loop for given RR,RA,RKAP,RDLT
!            2: PSIB adjusted eqch loop for given RR,RA,RKAP,RDLT

      MDLEQX = 0

!        MDLEQV : Order of extrapolation of psi in the vacuum region
!                 if positive, wall is linearly extended from plasma surface
!                 if negative, wall is extraporated by polynomials

      MDLEQV = 3

!        NPRINT: Level print out
!            0: no print
!            1: print first and last loop
!            2: print all loop

      NPRINT= 0

!        RGMIN: Minimum R of computation region [m]
!        RGMAX: Maxmum  R of computation region [m]
!        ZGMIN: Minimum Z of computation region [m]
!        ZGMAX: Maxmum  Z of computation region [m]
!        ZLIMP: Position of upper X points
!        ZLIMM: Position of lower X points

      RGMIN = 1.5D0
      RGMAX = 4.5D0
      ZGMIN =-2.0D0
      ZGMAX = 2.0D0
      ZLIMM =-2.5D0
      ZLIMP = 2.5D0

!        PSIB(0:5): Multipole moments of poloidal flux PSIRZ on boundary

      PSIB(0) =  2.0D0
      PSIB(1) =  0.5D0
      PSIB(2) =  0.D0
      PSIB(3) =  0.D0
      PSIB(4) =  0.D0
      PSIB(5) =  0.D0

!        NPFCMAX : Number of poloidal field coils (PFXs)
!        RIPFC(NPFC) : PFC coil current    [MA]
!        RPFC(NPFC)  : PFC coil position R [m]
!        ZPFC(NPFC)  : PFC coil position Z [m]
!        WPFC(NPFC)  : PFC coil width      [m]

      NPFCMAX = 0
      DO NPFC=1,NPFCM
         RIPFC(NPFC) = 0.D0
         RPFC(NPFC)  = 3.D0
         ZPFC(NPFC)  =-1.75D0
         WPFC(NPFC)  = 0.75D0
      ENDDO

      MODEFW=0 ! dangerous setting
      MODEFR=0 ! dangerous setting

   RETURN
 END SUBROUTINE eq_init
END MODULE eqinit
