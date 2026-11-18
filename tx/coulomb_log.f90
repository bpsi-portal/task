module mod_coulomb
  use, intrinsic :: iso_fortran_env, only : int32, dp => real64
  implicit none

contains

!*****************************************************************************
!   Coulomb logarithm Formulae given by Mitsuru HONDA (2011/06/13)
!     References: M. Honda, Jpn. J. Appl. Phys. 52 (2013) 108002.
!                 D.V. Sivukhin, Rev. Plasma Phys. Vol. 1 (1966)
!
! Note: Beam temperature defined by Tb = (2/3)Eb, where Eb: injection energy
!*****************************************************************************

!     Note: Coulog assumes that the ions share a common temperature, ti.
!           Set "2" to beam ions if "1" is the ion species.
  function coulog( zeff, ne, te, ti, A1, Z1, A2, Z2, tb ) result( lambda )
    ! arguments
    
    real(dp), intent(in) :: zeff, ne, te, ti, A1, Z1, A2, Z2
    real(dp), intent(in), optional :: tb
    ! zeff   : effective charge
    ! ne     : electron density in m^{-3} or 10^20 m^{-3}
    ! te, ti : electron and ion temperatures in eV or keV
    ! A1, A2 : mass number
    ! Z1, Z2 : ABSOLUTE charge number, e.g. Ze = 1
    ! tb     : mean temperature of beam ions defined by Tb = (2/3)Eb in eV or keV
    real(dp), parameter :: Ae = 5.446170219e-4_dp ! from CODATA 2010
    real(dp), parameter :: eps_max = 1.e3_dp ! margin
    real(dp) :: lambda, CD, t1, t2, fac_n, fac_t, eps
    real(dp) :: coef, coef1, coef2, thres_l, thres_r, Z1L, Z2L

    Z1L = abs(Z1)
    Z2L = abs(Z2)

    eps = eps_max * epsilon(1.0_dp)

    if (ne > 1.e10_dp) then
       ! density in m^{-3}, temperature in eV
       fac_n = 1.0_dp  ; fac_t = 1.0_dp
    else
       ! density in 10^20 m^{-3}, temperature in keV
       fac_n = 1.e20_dp ; fac_t = 1.e3_dp
    end if

    CD = ne * ( 1.0_dp / te + zeff / ti ) * fac_n / fac_t

    if ( present( tb ) ) then
       if ( abs( A1 - Ae ) < eps ) then
          t1 = te ; t2 = tb
          if ( abs( A2 - Ae ) < eps ) stop "Error!"
       else if ( abs( A2 - Ae ) < eps ) then
          t1 = tb ; t2 = te
          if ( abs( A1 - Ae ) < eps ) stop "Error!"
       else
          t1 = ti ; t2 = tb
       end if
    else
       t1 = te ; t2 = te
       if ( abs( A1 - Ae ) > eps ) t1 = ti
       if ( abs( A2 - Ae ) > eps ) t2 = ti
    end if

    coef1 = ( A1 + A2 ) / ( A2 * t1 + A1 * t2 ) / fac_t
    coef2 = A1 * A2 / ( A1 + A2 )
    coef  = log( coef1 * sqrt(CD) )

    thres_l = 1.0_dp / coef1
    thres_r = 2.45e4_dp * ( Z1L * Z2L )**2 * coef2
    if ( thres_l <= thres_r ) then
       !     Classical formula
       lambda = 30.4_dp - log( Z1L * Z2L ) - coef
    else
       !     Quantum-mechanical formula
       lambda = 35.4_dp - coef + 0.5_dp * log( coef1 * coef2 )
    end if

  end function coulog

!

  pure function coulog_gen( ne, te, CDi, A1, Z1, t1, A2, Z2, t2 ) result( lambda )
    ! arguments
    real(dp), intent(in) :: ne, te, CDi, A1, Z1, t1, A2, Z2, t2
    ! ne     : electron density in m^{-3} or 10^20 m^{-3}
    ! te     : electron temperature in eV or keV
    ! CDi    : sum_j (Z_j^2 n_j/T_j), where j denotes the ion species
    ! A1, A2 : mass number
    ! Z1, Z2 : ABSOLUTE charge number, e.g. Ze = 1
    ! t1, t2 : temperature in eV or keV
    real(dp) :: lambda, CD, fac_n, fac_t
    real(dp) :: coef, coef1, coef2, thres_l, thres_r, Z1L, Z2L

    Z1L = abs(Z1)
    Z2L = abs(Z2)

    if (ne > 1.e10_dp) then
       ! density in m^{-3}, temperature in eV
       fac_n = 1.0_dp  ; fac_t = 1.0_dp
    else
       ! density in 10^20 m^{-3}, temperature in keV
       fac_n = 1.e20_dp ; fac_t = 1.e3_dp
    end if

    CD = ( ne / te + CDi ) * fac_n / fac_t

    coef1 = ( A1 + A2 ) / ( A2 * t1 + A1 * t2 ) / fac_t
    coef2 = A1 * A2 / ( A1 + A2 )
    coef  = log( coef1 * sqrt(CD) )

    thres_l = 1.0_dp / coef1
    thres_r = 2.45e4_dp * ( Z1L * Z2L )**2 * coef2
    if ( thres_l <= thres_r ) then
       !     Classical formula
       lambda = 30.4_dp - log( Z1L * Z2L ) - coef
    else
       !     Quantum-mechanical formula
       lambda = 35.4_dp - coef + 0.5_dp * log( coef1 * coef2 )
    end if

  end function coulog_gen

!***************************************************************
!
!   Coulomb logarithm
!     Reference: NRL Formulary (2009)
!
!     imodel = 1 : electron - electron
!              2 : electron - ion
!              3 : ion - ion
!
!***************************************************************

  pure real(dp) function coulog_NRL(imodel, Ne, Te, Ni, Ti, PA, PZ) result(f)
    integer(int32), intent(in) :: imodel
    real(dp), intent(in) :: Ne, Te
    real(dp), intent(in), optional :: Ni, Ti, PA, PZ
    real(dp), parameter :: memp = 5.446170219e-4_dp ! me/mp
    real(dp) :: Ne_m3, Ni_m3, Te_eV, Ti_eV, rat_mass

    !     (NRL Plasma Formulary p34,35 (2007))
    Ne_m3 = Ne * 1.e20_dp ; Te_eV = Te * 1.e3_dp

    if (imodel == 1) then
       f = 30.4_dp - log(sqrt(Ne_m3)/(Te_eV**1.25_dp)) &
            &      - sqrt(1.e-5_dp+(log(Te_eV)-2.0_dp)**2/16.0_dp)
    else
       if(present(Ni) .and. present(Ti) .and. present(PA) .and. present(PZ)) then
          Ni_m3 = Ni * 1.e20_dp ; Ti_eV = Ti * 1.e3_dp
          if (imodel == 2) then
             rat_mass = memp / PA
             if(10.0_dp*PZ**2 > Ti_eV*rat_mass .and. 10.0_dp*PZ**2 < Te_eV) then ! usual case
                f = 30.9_dp - log(sqrt(Ne_m3)/Te_eV)
             else if(Te_eV > Ti_eV*rat_mass .and. Te_eV < 10.0_dp*PZ**2) then ! rare case
                f = 29.9_dp - log(sqrt(Ne_m3)*PZ/Te_eV**1.5_dp)
             else if(Te_eV < Ti_eV*PZ*rat_mass) then ! very rare case
                f = 36.9_dp - log(sqrt(Ni_m3)/Ti_eV**1.5_dp*PZ**2/PA)
             else ! assumption
                f = 30.9_dp - log(sqrt(Ne_m3)/Te_eV)
             end if
             !tokamaks       f = 37.8_dp - LOG(SQRT(Ne_m3)/(Te))
          else ! imodel = 3
             f = 29.9_dp - log(PZ**2/Ti_eV*sqrt(2.0_dp*Ni_m3*PZ**2/Ti_eV))
             !tokamaks       f = 40.3_dp - LOG(PZ**2/Ti*SQRT(2.0_dp*Ni_m3*PZ**2/Ti))
          end if
       else
        
       end if
    end if

  end function coulog_NRL

end module mod_coulomb
