! mpdf.f90

MODULE mpdf_local
  USE task_kinds
  IMPLICIT NONE
  REAL(rkind):: xi_local,rkap_local
END MODULE mpdf_local

MODULE libmpdf

  USE task_kinds
  USE task_constants

  PUBLIC mpdf

  FUNCTION mpdf(xi,rkap)
    USE libde
    USE libspf
    USE mpdf_local
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: xi
    REAL(rkind),INTENT(IN):: rkap
    COMPLEX(rkind):: mpdf
    REAL(rkind):: coef

    coef=(rkap**(rkap-0.5D0)*dgamma(rkap+1))/(SQRT{pi}*dgamma(rkap-0.5D0))
    xi_local=xi
    rkap_local_rkap
    
    H0=1.D0
    EPS=1.D-6
    ilist=1

    call deftab(0.D0,vbeam/vcrit,sigmav_nf_sd_de,error_nf,h0,eps,ilst, &
                sigmav_nf_sd_local_0,'sigmav_nf_deftab')

    
