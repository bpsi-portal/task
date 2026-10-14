module tx_interface

  !****************!
  !   txmmm.f90    !
  !****************!

  interface
     subroutine txmmm95(dNsdr,dTsdr,dQdr,cexb,gamma)
       real(8), intent(in) :: cexb
       real(8), dimension(:,:), intent(in) :: dNsdr,dTsdr
       real(8), dimension(:),   intent(in) :: dQdr
       real(8), dimension(:),   intent(out), optional :: gamma
     end subroutine txmmm95
  end interface

  !****************!
  !   txlib.f90    !
  !****************!

  interface APTOS
     subroutine APITOS(STR, NSTR, I)
       character(len=*), intent(INOUT) :: STR
       integer(4),       intent(INOUT) :: NSTR
       integer(4),       intent(IN)    :: I
     end subroutine APITOS

     subroutine APSTOS(STR, NSTR, INSTR, NINSTR)
       character(len=*), intent(INOUT) :: STR
       integer(4),       intent(INOUT) :: NSTR
       character(len=*), intent(IN)    :: INSTR
       integer(4),       intent(IN)    :: NINSTR
     end subroutine APSTOS

     subroutine APDTOS(STR, NSTR, D, FORM)
       character(len=*), intent(INOUT) :: STR
       integer(4),       intent(INOUT) :: NSTR
       real(8),          intent(IN)    :: D
       character(len=*), intent(IN)    :: FORM
     end subroutine APDTOS

     subroutine APRTOS(STR, NSTR, GR, FORM)
       character(len=*), intent(INOUT) :: STR
       integer(4),       intent(INOUT) :: NSTR
       real(4),          intent(IN)    :: GR
       character(len=*), intent(IN)    :: FORM
     end subroutine APRTOS
  end interface

  interface
     subroutine TOUPPER(KTEXT)
       character(len=*), intent(INOUT) ::  KTEXT
     end subroutine TOUPPER
  end interface

  interface
     subroutine KSPLIT_TX(KKLINE,KID,KKLINE1,KKLINE2)
       character(LEN=*),  intent(IN)  :: KKLINE
       character(LEN=1),  intent(IN)  :: KID
       character(LEN=*), intent(OUT) :: KKLINE1, KKLINE2
     end subroutine KSPLIT_TX
  end interface

  interface
     pure real(8) function DERIVF(NR,R,F,NRMAX)
       integer(4), intent(in) :: NR, NRMAX
       real(8), dimension(0:NRMAX), intent(in)  :: R
       real(8), dimension(0:NRMAX), intent(in)  :: F
     end function DERIVF
  end interface

  interface
     function dfdx(x,f,nmax,mode,daxs,dbnd)
       integer(4), intent(in) :: nmax, mode
       real(8), dimension(0:nmax), intent(in) :: x, f
       real(8), dimension(0:nmax) :: dfdx
       real(8), optional :: daxs, dbnd
     end function dfdx
  end interface

  interface
     subroutine INTDERIV3(X,R,intX,FVAL,NRMAX,ID)
       integer(4), intent(in) :: NRMAX, ID
       real(8), intent(in), dimension(0:NRMAX) :: X, R
       real(8), intent(in) :: FVAL
       real(8), intent(out), dimension(0:NRMAX) :: intX
     end subroutine INTDERIV3
  end interface

  interface
     pure real(8) function LORENTZ(R,C1,C2,W1,W2,RC1,RC2,AMP)
       real(8), intent(in) :: r, c1, c2, w1, w2, rc1, rc2
       real(8), intent(in), optional :: AMP
     end function LORENTZ
  end interface

  interface
     pure real(8) function LORENTZ_PART(R,W1,W2,RC1,RC2,ID)
       real(8), intent(in) :: r, w1, w2, rc1, rc2
       integer(4), intent(in) :: ID
     end function LORENTZ_PART
  end interface

  interface
     subroutine BISECTION(f,cl1,cl2,w1,w2,rc1,rc2,amp,s,valmax,val,valmin)
       real(8), external :: f
       real(8), intent(in) :: cl1, cl2, w1, w2, rc1, rc2, amp, s, valmax
       real(8), intent(in), optional :: valmin
       real(8), intent(out) :: val
     end subroutine BISECTION
  end interface

  interface
     subroutine inexpolate(nmax_in,r_in,dat_in,nmax_std,r_std,iedge,dat_out,ideriv,nrbound,idx)
       integer(4), intent(in) :: nmax_in, nmax_std, iedge
       integer(4), intent(in),  optional :: ideriv, idx
       integer(4), intent(out), optional :: nrbound
       real(8), dimension(1:nmax_in), intent(in) :: r_in, dat_in
       real(8), dimension(0:nmax_std), intent(in) :: r_std
       real(8), dimension(0:nmax_std), intent(out) :: dat_out
     end subroutine inexpolate
  end interface

  interface
     pure real(8) function fgaussian(x,mu,sigma,norm)
       real(8), intent(in) :: x, mu, sigma
       integer(4), intent(in), optional :: norm
     end function fgaussian
  end interface

  interface
     pure real(8) function moving_average(i,f,imax,iend)
       integer(4), intent(in) :: i, imax
       integer(4), intent(in), optional :: iend
       real(8), dimension(0:imax), intent(in) :: f
     end function moving_average
  end interface

  interface
     subroutine replace_interpolate_value(val,index,xarray,varray)
       use mod_spln
       integer(4), intent(in) :: index
       real(8), dimension(:), intent(in) :: xarray, varray ! 0:NRMAX
       real(8), intent(out) :: val
     end subroutine replace_interpolate_value
  end interface

  !****************!
  !   txfile.f90   !
  !****************!

  interface
     subroutine TXLOAD(IST)
       integer(4), intent(out) :: IST
     end subroutine TXLOAD
  end interface

  interface
     subroutine TXEXTCHK(mode)
       integer(4), intent(in) :: mode
     end subroutine TXEXTCHK
  end interface

  interface
     subroutine TXGLOD(IST)
       integer(4), intent(out) :: IST
     end subroutine TXGLOD
  end interface

  interface
     real(8) function rLINEAVE(Rho)
       real(8), intent(IN) :: Rho
     end function rLINEAVE
  end interface

  interface
     integer(4) function detect_datatype(kchar)
       character(len=*), intent(in) :: kchar
     end function detect_datatype
  end interface

  interface
     subroutine initprof_input(nr, idx, out)
       integer(4), optional :: nr, idx
       real(8), optional :: out
     end subroutine initprof_input
  end interface

#ifndef nonGSAF
  !***************************!
  !   txg3d.f90, txg2d.f90    !
  !***************************!

  interface
     subroutine TXGRUR(GX,GTX,GYL,NRMAX,NGT,NGTM)!,STR,KV,INQ)
       integer(4),                  intent(in) :: NRMAX, NGT, NGTM!, INQ
       real(4), dimension(0:NRMAX), intent(in) :: GX
       real(4), dimension(0:NGT),   intent(in) :: GTX
       real(4), dimension(0:NRMAX,0:NGTM), intent(in) :: GYL
!       character(LEN=80),intent(IN):: STR, KV
     end subroutine TXGRUR
  end interface

  interface
     subroutine TXGRURA(GX,GTX,GYL,NRMAX,NGT,NGTM)!,STR,KV,INQ)
       integer(4),                  intent(in) :: NRMAX, NGT, NGTM!, INQ
       real(4), dimension(0:NRMAX), intent(in) :: GX
       real(4), dimension(0:NGT),   intent(in) :: GTX
       real(4), dimension(0:NRMAX,0:NGTM), intent(in) :: GYL
!       character(LEN=80),intent(IN):: STR, KV
     end subroutine TXGRURA
  end interface

  interface
     subroutine TXGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NXM,NXMAX,NYMAX,STR,KV,MODE)
       real(4),    intent(IN) :: GX1, GX2, GY1, GY2
       integer(4), intent(IN) :: NXM, NXMAX, NYMAX, MODE
       real(4), dimension(NXMAX),     intent(IN) :: GX
       real(4), dimension(NYMAX),     intent(IN) :: GY
       real(4), dimension(NXM,NYMAX), intent(IN) :: GZ
       character(LEN=80) :: STR, KV
     end subroutine TXGR3D
  end interface
#endif

end module tx_interface
