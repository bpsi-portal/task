module mod_savgol
  use, intrinsic :: iso_fortran_env, only : int32, dp => real64
  implicit none
  private
  public :: savgol_filter

contains

  subroutine savgol_filter(nl,nr,ld,m,n1,y,flag)
    !-----------------------------------------------------------------------------------
    ! This routine is used to perform the Savitzky-Golay algorithm.
    !-----------------------------------------------------------------------------------
    !    nl:: input, integer, the number of leftward data points used.
    !    nr:: input, integer, the number of rightward data points used.
    !    ld:: input, integer, the order of the derivative desired.
    !     m:: input, integer, the order of the smoothing polynomial.
    !    n1:: input, integer, the number of data points.
    ! y(n1):: input/output, real values, the data to be smoothed.
    !  flag:: output, integer, error message, 0=success, 1=failure.
    !-----------------------------------------------------------------------------------
    ! The points are taken as equally spaced: on a non-uniform mesh the polynomial is
    ! fitted in the index, not in the coordinate.  For ld = 0 every point is smoothed,
    ! the first nl and the last nr points with the window shifted inside the data.
    !-----------------------------------------------------------------------------------
    ! Author: Peng Jun, 2019.03.26.
    ! Modified by HONDA Mitsuru
    !-----------------------------------------------------------------------------------
    ! Dependence:: subroutine savgol.
    ! -----------------------------------------------------------------------------------

    implicit none
    integer(int32), intent(in)    :: nl, nr, ld, m, n1
    real   (dp),    intent(inout) :: y(n1)
    integer(int32), intent(out)   :: flag
    ! Local variables.
    integer(int32) :: i, j, xl(nl+nr+1)
    real   (dp)    :: y0(n1), coef(nl+nr+1)

    xl(1) = 0
    y0 = y

    do i=1, nl
       xl(i+1) = -i
    end do

    do i=1, nr
       xl(1+nl+i) = nr-i+1
    end do

    call savgol(nl,nr,ld,m,coef,flag)

    if (flag/=0) return

    do i=1, n1-nr
       y(i) = 0.0_dp
       do j=1, nl+nr+1
          if (i+xl(j) > 0) then
!             write(6,*) j, y(i),coef(j)*y0(i+xl(j))
             y(i) = y(i) + coef(j)*y0(i+xl(j))
          end if
       end do
    end do

!    if (ld==0) then
!       y(1:nl) = y0(1:nl)
!       y(n1-nr+1:n1) = y0(n1-nr+1:n1)
!    else 
    if( ld /= 0 ) then
       y(1:nl) = y(nl+1)
       y(n1-nr+1:n1) = y(n1-nr)
       return
    end if

    ! The first nl and the last nr points: the loop above truncated the window at the
    ! left end without renormalizing it (a constant came out as 0.67 of itself at the
    ! first point) and did not reach the last nr points.  Use there the same number of
    ! points and the same polynomial order with the window shifted to lie inside the
    ! data: the least-squares polynomial of the nl+nr+1 points next to the end, evaluated
    ! at the point.  It reproduces polynomials up to order m, constants in particular;
    ! the interior points are those of the loop above.  With fewer than nl+nr+1 points
    ! there is no such window and the data are returned unchanged.
    if (n1 < nl+nr+1) then
       y(1:min(nl,n1)) = y0(1:min(nl,n1))
       return
    end if
    do i=1, nl
       call smooth_at(i, i-1, nl+nr-(i-1))
    end do
    do i=n1-nr+1, n1
       call smooth_at(i, nl+nr-(n1-i), n1-i)
    end do

  contains

    subroutine smooth_at(ip, nlb, nrb)
      integer(int32), intent(in) :: ip, nlb, nrb
      integer(int32) :: jj, xlb(nlb+nrb+1)
      real   (dp)    :: cb(nlb+nrb+1)
      xlb(1) = 0
      do jj=1, nlb
         xlb(jj+1) = -jj
      end do
      do jj=1, nrb
         xlb(1+nlb+jj) = nrb-jj+1
      end do
      call savgol(nlb,nrb,0,m,cb,flag)
      y(ip) = 0.0_dp
      do jj=1, nlb+nrb+1
         y(ip) = y(ip) + cb(jj)*y0(ip+xlb(jj))
      end do
    end subroutine smooth_at

  end subroutine savgol_filter

  subroutine savgol(nl,nr,ld,m,coef,flag)
    !-----------------------------------------------------------------------------------
    ! This routine is used to calculate a set of Savitzky-Golay filter coefficients.
    !-----------------------------------------------------------------------------------
    !            nl:: input, integer, the number of leftward data points used.
    !            nr:: input, integer, the number of rightward data points used.
    !            ld:: input, integer, the order of the derivative desired.
    !             m:: input, integer, the order of the smoothing polynomial.
    ! coef(nl+nr+1):: output, real values, calculated coefficents in wrap-around order.
    !          flag:: output, integer, error message, 0=success, 1=failure.
    !-----------------------------------------------------------------------------------
    ! Author: Peng Jun, 2019.03.20.
    ! Modified by HONDA Mitsuru
    !-----------------------------------------------------------------------------------
    ! Dependence:: subroutine ludcmp;
    !              subroutine lubksb.
    !-----------------------------------------------------------------------------------
    ! Reference: Press et al, 1986. Numerical recipes in Fortran 77, 
    !            the Art of Scientific Computing, second edition. 
    ! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.646 IN Press et al.
    ! -----------------------------------------------------------------------------------
    implicit none
    integer(int32), intent(in)    :: nl, nr, ld, m
    real   (dp),    intent(inout) :: coef(nl+nr+1)
    integer(int32), intent(out)   :: flag
    ! Local variables.
    integer(int32) :: imj, ipj, k, kk, mm, indx(m+1)
    real   (dp)    :: d, fac, summ, a(m+1,m+1), b(m+1)

    flag = 0

    if (nl < 0 .or. nr < 0 .or. ld > m .or. nl+nr < m) then
       flag = 1
       return
    end if

    do ipj=0, 2*m
       summ = 0.0_dp
       if (ipj .eq. 0) summ = 1.0_dp

       do k=1, nr
          summ = summ + (real(k,dp))**ipj
       end do

       do k=1, nl
          summ = summ + (real(-k,dp))**ipj
       end do

       mm = min(ipj, 2*m-ipj)
       do imj=-mm, mm, 2
          a(1+(ipj+imj)/2,1+(ipj-imj)/2) = summ
       end do
    end do

    call ludcmp(a,m+1,indx,d,flag)

    if (flag .ne. 0) return
    b = 0.0_dp
    b(ld+1) = 1.0_dp

    call lubksb(a,m+1,indx,b)

    coef = 0.0_dp
    do k=-nl, nr
       summ = b(1)
       fac = 1.0_dp
       do mm=1, m
          fac = fac * k
          summ = summ + b(mm+1) * fac
       end do
       kk = mod(nl+nr+1-k, nl+nr+1) + 1
       coef(kk) = summ
    end do

  end subroutine savgol

  subroutine ludcmp(a,n,indx,d,flag)
    !-------------------------------------------------------------------------
    !This routine is used in combination with lubksb to solve 
    !linear equations or invert a matrix.
    !-------------------------------------------------------------------------
    !  a(n,n):: input, real values, a matrix to be decomposed.
    !       n:: input, integer, the dimension of the matrix.
    ! indx(n):: output, integer values, vector that records the row 
    !           permutation effected by the partial pivoting.
    !       d:: output, integer, output as 1 or -1 depending on whether 
    !           the number of row interchanges was even or odd.
    !    flag:: output, integer, error message, 0=success, 1=singular matrix.
    !-------------------------------------------------------------------------
    ! Author: Peng Jun, 2019.03.20.
    ! Modified by HONDA Mitsuru
    !-------------------------------------------------------------------------
    ! Dependence:: No.--------------------------------------------------------
    !-------------------------------------------------------------------------
    ! Reference: Press et al, 1986. Numerical recipes in Fortran 77, 
    !            the Art of Scientific Computing, second edition. 
    ! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.38 IN Press et al.
    ! ------------------------------------------------------------------------
    implicit none
    integer(int32), intent(in)    :: n
    integer(int32), intent(out)   :: indx(n), flag
    real   (dp),    intent(inout) :: a(n,n)
    real   (dp),    intent(out)   :: d
    ! Local variables.
    integer(int32) :: i, j, k, imax
    real   (dp)    :: aamax, dum, summ, vv(n)

    indx = 0
    flag = 0
    d = 1.0_dp

    do i=1, n
       aamax = 0.0_dp
       do j=1, n
          aamax = max(abs(a(i,j)), aamax)
       end do

       if (aamax == 0.0_dp) then 
          flag = 1
          return
       end if

       vv(i) = 1.0_dp/aamax
    end do

    do j=1, n
       do i=1, j-1
          summ = a(i,j)
          do k=1, i-1
             summ = summ - a(i,k) * a(k,j)
          end do
          a(i,j) = summ
       end do

       aamax = 0.0_dp

       do i=j, n
          summ = a(i,j)
          do k=1, j-1
             summ = summ - a(i,k) * a(k,j)
          end do

          a(i,j) = summ

          dum = vv(i) * abs(summ)

          if (dum >= aamax) then
             imax = i
             aamax = dum
          end if
       end do

       if (j /= imax) then
          do k=1, n
             dum = a(imax,k)
             a(imax,k) = a(j,k)
             a(j,k) = dum
          end do

          d = -d
          vv(imax) = vv(j)
       end if

       indx(j) = imax

       if (a(j,j) == 0.0_dp) a(j,j) = tiny(0.0_dp)

       if (j /= n) then
          dum = 1.0_dp / a(j,j)
          do i=j+1, n
             a(i,j) = a(i,j) * dum
          end do
       end if
    end do

  end subroutine ludcmp

  subroutine lubksb(a,n,indx,b)
    !-------------------------------------------------------------------------
    !  a(n,n):: input, real values, the LU decomposition of a matrix.
    !       n:: input, integer, the dimenstion of the matrix.
    ! indx(n):: input,  integer values, vector that records the row 
    !           permutation effected by the partial pivoting.
    !    b(n):: output, real values, the solution vector X for 
    !                   linear equations A*X=B.
    !-------------------------------------------------------------------------
    ! Author: Peng Jun, 2019.03.18.
    ! Modified by HONDA Mitsuru
    !-------------------------------------------------------------------------
    ! Dependence:: No.--------------------------------------------------------
    !-------------------------------------------------------------------------
    ! Reference: Press et al, 1986. Numerical recipes in Fortran 77, 
    !            the Art of Scientific Computing, second edition. 
    ! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.39 IN Press et al.
    ! -------------------------------------------------------------------------
    implicit none
    integer(int32), intent(in)    :: n, indx(n)
    real   (dp),    intent(in)    :: a(n,n)
    real   (dp),    intent(inout) :: b(n)
    ! Local variables.
    integer(int32) :: i, ii, j, ll
    real   (dp)    :: summ

    ii = 0

    do i=1, n
       ll = indx(i)
       summ = b(ll)
       b(ll) = b(i)

       if (ii /= 0) then
          do j=ii, i-1
             summ = summ - a(i,j) * b(j)
          end do
       else if (summ /= 0.0_dp) then
          ii = i
       end if

       b(i) = summ
    end do

    do i=n, 1, -1
       summ = b(i)
       do j=i+1, n
          summ = summ - a(i,j) * b(j)
       end do
       b(i) = summ / a(i,i)
    end do

  end subroutine lubksb

end module mod_savgol
