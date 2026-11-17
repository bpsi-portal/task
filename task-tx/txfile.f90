!***************************************************************
!
!   Extension of the save file
!
!***************************************************************
! A save file is the base layout of TXSAVE (unchanged since 2025-12: 33 records, the
! last one thrp) followed by one extension, which carries every input and every state
! added after that layout, item by item and by name:
!   header     : EXTSENT (-1), EXTTAG ('TXEXTREC'), EXTVER, the number of items
!   each item  : a record (name, kind, length) and a record of the values
!                (kind 1: integer(4), 2: real(8), 3: character; a logical is 0 or 1)
!   end        : a record ('END', 0, 0)
! The items are listed once, in savext_items, which TXSAVE runs to write them and TXLOAD
! to take them from the items it has read.  An item that a file does not have keeps its
! current value and is reported; an item of the file that the build does not know is
! reported and skipped.  So an item is added by one line in savext_items, and EXTVER
! changes only when an item changes its meaning.
! A TXLOAD of the builds with the development extension records (1-7, BR-43 to BR-51)
! reads the leading -1 as the version of those records and sees no extension.  This
! TXLOAD does not read those records: a file that has them loads as one of the base
! layout.
module tx_savext
  implicit none
  private
  public :: EXTVER, EXTSENT, EXTTAG, XM_COUNT, XM_WRITE, XM_READ, xmode, xio, xcount, &
       &    xreg, xreg_ia, xreg_ra, xclear, xread_items, xhas, xok, xreport
  integer(4), parameter :: EXTVER = 1, EXTSENT = -1
  character(len=8), parameter :: EXTTAG = 'TXEXTREC'
  integer(4), parameter :: XM_COUNT = 1, XM_WRITE = 2, XM_READ = 3
  integer(4), parameter :: KINT = 1, KREAL = 2, KCHAR = 3
  integer(4) :: xmode = 0, xio = 0, xcount = 0
  type xitem
     character(len=16) :: name
     integer(4) :: kind, n
     logical :: used, bad
     integer(4), allocatable :: iv(:)
     real(8), allocatable :: rv(:)
     character(len=:), allocatable :: cv
  end type xitem
  type(xitem), allocatable :: xtab(:)
  integer(4) :: nx = 0, nmiss = 0, nbad = 0
  character(len=3000) :: xmiss = ' '
  character(len=1000) :: xbad = ' '
  interface xreg
     module procedure xreg_i, xreg_ia, xreg_r, xreg_ra, xreg_l, xreg_c
  end interface xreg

contains

  subroutine xclear
    if(allocated(xtab)) deallocate(xtab)
    nx = 0 ; nmiss = 0 ; nbad = 0 ; xmiss = ' ' ; xbad = ' '
  end subroutine xclear

  ! The nitem items after the header and the end record.  ist /= 0 on a read error or
  ! a missing end record; the items read before it are kept.
  subroutine xread_items(io, nitem, ist)
    integer(4), intent(in) :: io, nitem
    integer(4), intent(out) :: ist
    character(len=16) :: name
    integer(4) :: kind, n, i
    call xclear
    allocate(xtab(max(nitem,0)))
    ist = 0
    do i = 1, nitem
       read(io, iostat=ist) name, kind, n
       if(ist /= 0) return
       xtab(i)%name = name ; xtab(i)%kind = kind ; xtab(i)%n = n ; xtab(i)%used = .false. ; xtab(i)%bad = .false.
       select case(kind)
       case(KINT)
          allocate(xtab(i)%iv(n)) ; read(io, iostat=ist) xtab(i)%iv
       case(KREAL)
          allocate(xtab(i)%rv(n)) ; read(io, iostat=ist) xtab(i)%rv
       case(KCHAR)
          allocate(character(len=n) :: xtab(i)%cv) ; read(io, iostat=ist) xtab(i)%cv
       case default
          read(io, iostat=ist)
       end select
       if(ist /= 0) return
       nx = i
    end do
    read(io, iostat=ist) name, kind, n
    if(ist == 0 .and. name /= 'END') ist = 1
  end subroutine xread_items

  integer(4) function xfind(name)
    character(len=*), intent(in) :: name
    integer(4) :: i
    xfind = 0
    do i = 1, nx
       if(xtab(i)%name == name) then
          xfind = i
          return
       end if
    end do
  end function xfind

  ! whether the file read last has the item
  logical function xhas(name)
    character(len=*), intent(in) :: name
    xhas = xfind(name) > 0
  end function xhas

  ! whether the item was taken from the file read last (there, and of the kind and length of
  ! the build); valid after the pass of savext_items that read the file
  logical function xok(name)
    character(len=*), intent(in) :: name
    integer(4) :: k
    k = xfind(name)
    xok = .false.
    if(k > 0) xok = .not. xtab(k)%bad
  end function xok

  ! the item k of the file for the item name of kind kind and length n, or 0
  integer(4) function xget(name, kind, n)
    character(len=*), intent(in) :: name
    integer(4), intent(in) :: kind, n
    xget = xfind(name)
    if(xget == 0) then
       nmiss = nmiss + 1
       if(len_trim(xmiss) + len_trim(name) + 2 < len(xmiss)) xmiss = trim(xmiss)//' '//trim(name)
    else
       xtab(xget)%used = .true.
       if(xtab(xget)%kind /= kind .or. (kind /= KCHAR .and. xtab(xget)%n /= n)) then
          nbad = nbad + 1
          if(len_trim(xbad) + len_trim(name) + 2 < len(xbad)) xbad = trim(xbad)//' '//trim(name)
          xtab(xget)%bad = .true.
          xget = 0
       end if
    end if
  end function xget

  subroutine xhead(name, kind, n)
    character(len=*), intent(in) :: name
    integer(4), intent(in) :: kind, n
    character(len=16) :: name16
    name16 = name
    write(xio) name16, kind, n
  end subroutine xhead

  subroutine xreg_ia(name, v, n)
    character(len=*), intent(in) :: name
    integer(4), intent(in) :: n
    integer(4), intent(inout) :: v(n)
    integer(4) :: k
    select case(xmode)
    case(XM_COUNT)
       xcount = xcount + 1
    case(XM_WRITE)
       call xhead(name, KINT, n) ; write(xio) v
    case(XM_READ)
       k = xget(name, KINT, n)
       if(k > 0) v(1:n) = xtab(k)%iv(1:n)
    end select
  end subroutine xreg_ia

  subroutine xreg_ra(name, v, n)
    character(len=*), intent(in) :: name
    integer(4), intent(in) :: n
    real(8), intent(inout) :: v(n)
    integer(4) :: k
    select case(xmode)
    case(XM_COUNT)
       xcount = xcount + 1
    case(XM_WRITE)
       call xhead(name, KREAL, n) ; write(xio) v
    case(XM_READ)
       k = xget(name, KREAL, n)
       if(k > 0) v(1:n) = xtab(k)%rv(1:n)
    end select
  end subroutine xreg_ra

  subroutine xreg_i(name, v)
    character(len=*), intent(in) :: name
    integer(4), intent(inout) :: v
    integer(4) :: w(1)
    w(1) = v ; call xreg_ia(name, w, 1) ; v = w(1)
  end subroutine xreg_i

  subroutine xreg_r(name, v)
    character(len=*), intent(in) :: name
    real(8), intent(inout) :: v
    real(8) :: w(1)
    w(1) = v ; call xreg_ra(name, w, 1) ; v = w(1)
  end subroutine xreg_r

  subroutine xreg_l(name, v)
    character(len=*), intent(in) :: name
    logical, intent(inout) :: v
    integer(4) :: w(1)
    w(1) = merge(1, 0, v) ; call xreg_ia(name, w, 1) ; v = w(1) /= 0
  end subroutine xreg_l

  subroutine xreg_c(name, v)
    character(len=*), intent(in) :: name
    character(len=*), intent(inout) :: v
    integer(4) :: k
    select case(xmode)
    case(XM_COUNT)
       xcount = xcount + 1
    case(XM_WRITE)
       call xhead(name, KCHAR, len(v)) ; write(xio) v
    case(XM_READ)
       k = xget(name, KCHAR, len(v))
       if(k > 0) v = xtab(k)%cv
    end select
  end subroutine xreg_c

  ! after a LOAD: the items the file did not have (kept at their current values), those
  ! of another kind or length (kept) and those this build does not know (skipped)
  subroutine xreport
    integer(4) :: i, nunk
    character(len=1000) :: unk
    if(nmiss > 0) then
       write(6,'(A,I3,A)') ' ##', nmiss, ' items not in the file, kept at their current values:'
       call xwrap(xmiss)
    end if
    if(nbad > 0) then
       write(6,'(A,I3,A)') ' ## WARNING:', nbad, ' items of the file of another kind or length, not read:'
       call xwrap(xbad)
    end if
    nunk = 0 ; unk = ' '
    do i = 1, nx
       if(.not. xtab(i)%used) then
          nunk = nunk + 1
          if(len_trim(unk) + len_trim(xtab(i)%name) + 2 < len(unk)) unk = trim(unk)//' '//trim(xtab(i)%name)
       end if
    end do
    if(nunk > 0) then
       write(6,'(A,I3,A)') ' ##', nunk, ' items of the file not known to this build, skipped:'
       call xwrap(unk)
    end if
  end subroutine xreport

  subroutine xwrap(list)
    character(len=*), intent(in) :: list
    integer(4) :: i, j, n
    n = len_trim(list) ; i = 1
    do while(i <= n)
       j = min(n, i + 99)
       if(j < n) then
          do while(j > i .and. list(j+1:j+1) /= ' ')
             j = j - 1
          end do
       end if
       write(6,'(4X,A)') trim(adjustl(list(i:j)))
       i = j + 1
    end do
  end subroutine xwrap

end module tx_savext

!***************************************************************
!
!   Identification of the external data of a process
!
!***************************************************************
! The NBI data read by F (ascii_input: 1 OFMC OrbitEffectDist.dat, 2 OFMC Torque.txt,
! 3 an ASCII file) and the ripple data read by F R stay in the process: a save file does
! not carry them, and a LOAD does not change them.  The extension of the save file
! records what the run that wrote it had (the kind, the file, the number of data and of
! lines, the total power and a hash of the values read), and a LOAD compares it with
! what the process has, again at the first CONT after the LOAD (F may come in between).
! A difference is reported only; nothing is added, removed or stopped.
module tx_extid
  use tx_commons, only : iflag_file, infiles, n_infiles, infile_name, IRPIN, nmax_file
  implicit none
  private
  public :: extid, xid, extid_process, extid_loaded, extid_check, extid_clear
  type extid
     integer(4) :: kind = 0, ndata = 0, nol = 0, irpin = 0
     real(8) :: totp = 0.d0
     character(len=100) :: fname = ' '
     character(len=16) :: hash = ' '
  end type extid
  type(extid) :: xid                 ! the items of the save file (savext_items)
  type(extid) :: id_file, id_load    ! of the file loaded last, of the process at that LOAD
  logical :: have_file = .false.     ! the file had the record
  logical :: pending = .false.       ! a LOAD not yet followed by CONT or RUN
  logical :: differed = .false.      ! the check at that LOAD reported a difference
  real(8) :: pnbhex_file = 0.d0      ! PNBH - (PNBHP + PNBHT1 + PNBHT2) of the file [MW]

contains

  ! the external data of the process; the hash is over the values that ascii_input
  ! sets for the kind (the others are not defined)
  type(extid) function extid_process() result(id)
    integer(4) :: i, j, n
    integer(8) :: h
    id%irpin = IRPIN
    if(iflag_file < 1 .or. iflag_file > 3 .or. .not. allocated(infiles)) return
    id%kind = iflag_file ; id%ndata = n_infiles ; id%fname = infile_name
    id%nol = infiles(1)%nol
    h = int(iflag_file, 8)
    do i = 1, min(n_infiles, size(infiles))
       n = min(max(infiles(i)%nol, 0), nmax_file)
       h = mix(h, int(n, 8))
       do j = 1, n
          h = mix(h, transfer(infiles(i)%r(j), 0_8))
          h = mix(h, transfer(infiles(i)%data(j), 0_8))
       end do
       select case(iflag_file)
       case(1)
          h = mix(h, transfer(infiles(i)%totS, 0_8))
          h = mix(h, transfer(infiles(i)%totP, 0_8))
       case(2)
          h = mix(h, transfer(infiles(i)%totS, 0_8))
       case(3)
          h = mix(h, int(infiles(i)%ncol_mesh, 8))
          h = mix(h, int(infiles(i)%ncol_data, 8))
          do j = 1, len(infiles(i)%name)
             h = mix(h, int(ichar(infiles(i)%name(j:j)), 8))
          end do
       end select
    end do
    if(iflag_file == 1 .and. n_infiles >= 2) id%totp = infiles(2)%totP * 1.d-6   ! PNBHex (txauxs)
    write(id%hash,'(Z16.16)') h
  end function extid_process

  integer(8) function mix(h, v)
    integer(8), intent(in) :: h, v
    mix = ieor(ishftc(h, 7), v)
  end function mix

  ! at a LOAD, from the items of the file (xid) and its NBI powers
  subroutine extid_loaded(have, pnbhex)
    logical, intent(in) :: have
    real(8), intent(in) :: pnbhex
    have_file = have ; id_file = xid ; pnbhex_file = pnbhex
    id_load = extid_process()
    pending = .true. ; differed = .false.
  end subroutine extid_loaded

  ! after a RUN: no loaded state to check any more
  subroutine extid_clear
    pending = .false.
  end subroutine extid_clear

  ! when = 1 at a LOAD, 2 at the first CONT after it (checked again only if the LOAD
  ! check reported a difference or the external data changed since, by F)
  subroutine extid_check(when)
    integer(4), intent(in) :: when
    type(extid) :: p
    character(len=8) :: at
    if(.not. pending) return
    p = extid_process()
    if(when == 2) then
       pending = .false.
       if(.not. differed .and. same(p, id_load)) return
    end if
    at = merge('LOAD    ', 'CONT    ', when == 1)
    differed = .false.
    if(have_file) then
       if(id_file%kind /= 0 .and. p%kind == 0) then
          call warn(trim(at)//': the run that wrote the file had external NBI data ('//trim(desc(id_file)) &
               & //'); this process has none, and the run goes on without them unless F reads them')
       else if(id_file%kind == 0 .and. p%kind /= 0) then
          call warn(trim(at)//': the run that wrote the file had no external NBI data; those of this process (' &
               & //trim(desc(p))//') are added to the loaded run')
       else if(p%kind /= 0 .and. .not. same_nbi(p, id_file)) then
          call warn(trim(at)//': the external NBI data of this process ('//trim(desc(p)) &
               & //') are not those of the run that wrote the file ('//trim(desc(id_file))//')')
       else if(p%kind /= 0) then
          write(6,'(A)') ' ## '//trim(at)//': the external NBI data of this process are those of the file (' &
               & //trim(desc(p))//')'
       end if
       if(id_file%irpin /= p%irpin) then
          if(id_file%irpin /= 0) then
             call warn(trim(at)//': the run that wrote the file had ripple data (F R); this process has none')
          else
             call warn(trim(at)//': the run that wrote the file had no ripple data; those of this process (F R) are used')
          end if
       end if
    else
       ! a file without the record: the NBI power of the file beyond its inputs tells
       ! whether its run had OFMC data (kind 1, the only kind with a power of its own)
       if(pnbhex_file > 1.d-6 .and. p%kind == 0) then
          call warn(trim(at)//': the file has an NBI power beyond PNBHP + PNBHT1 + PNBHT2 (' &
               & //trim(fmt(pnbhex_file))//' MW), external NBI data (F) of its run, which this process does not have')
       else if(pnbhex_file <= 1.d-6 .and. p%kind == 1) then
          call warn(trim(at)//': the file has no NBI power beyond its inputs; the external NBI data of this process (' &
               & //trim(desc(p))//') are added to the loaded run')
       else if(p%kind == 1 .and. abs(p%totp - pnbhex_file) > 1.d-6 * max(1.d0, abs(pnbhex_file))) then
          call warn(trim(at)//': the NBI power of the file beyond its inputs ('//trim(fmt(pnbhex_file)) &
               & //' MW) is not that of the external NBI data of this process ('//trim(desc(p))//')')
       else if(p%kind /= 0) then
          write(6,'(A)') ' ## '//trim(at)//': external NBI data of this process ('//trim(desc(p)) &
               & //'); the file does not record those of its run (its NBI power beyond the inputs agrees)'
       end if
    end if
  end subroutine extid_check

  subroutine warn(text)
    character(len=*), intent(in) :: text
    write(6,'(A)') ' ## WARNING: '//text
    differed = .true.
  end subroutine warn

  logical function same(a, b)
    type(extid), intent(in) :: a, b
    same = same_nbi(a, b) .and. a%irpin == b%irpin
  end function same

  logical function same_nbi(a, b)
    type(extid), intent(in) :: a, b
    same_nbi = a%kind == b%kind .and. a%ndata == b%ndata .and. a%nol == b%nol .and. a%hash == b%hash &
         &     .and. a%fname == b%fname
  end function same_nbi

  character(len=200) function desc(a)
    type(extid), intent(in) :: a
    character(len=*), parameter :: KINDS(3) = [character(len=24) :: 'OFMC OrbitEffectDist.dat', &
         & 'OFMC Torque.txt', 'ASCII']
    write(desc,'(A,I1,3A,I2,A,I3,A)') 'kind ', a%kind, ' ', trim(KINDS(max(1,min(3,a%kind)))), &
         & ', file '//trim(a%fname)//',', a%ndata, ' data,', a%nol, ' lines'
    if(a%kind == 1) desc = trim(desc)//', '//trim(fmt(a%totp))//' MW'
    desc = trim(desc)//', hash '//a%hash
  end function desc

  character(len=12) function fmt(x)
    real(8), intent(in) :: x
    write(fmt,'(F12.4)') x
    fmt = adjustl(fmt)
  end function fmt

end module tx_extid

!***************************************************************
!
!   Check of the external data after a LOAD (TXMENU: CONT, RUN)
!
!***************************************************************
! mode 2: at CONT (tx_extid, extid_check); 0: at RUN (the loaded state is gone)

subroutine TXEXTCHK(mode)
  use tx_extid, only : extid_check, extid_clear
  implicit none
  integer(4), intent(in) :: mode
  if(mode == 0) then
     call extid_clear
  else
     call extid_check(mode)
  end if
end subroutine TXEXTCHK


!***************************************************************
!
!   The items of the extension of the save file (tx_savext)
!
!***************************************************************
! TXSAVE writes them in this order and TXLOAD takes them by name (xreg).

subroutine savext_items
  use tx_commons, only : NSM, MDLNSUP, RWALL, RWALLZ, RWALLIZ, FSTPTM, PROFD2C, PROFD2W, &
       & MDLSOLT, DSOLT, CHISOLT, RMUSOLT, WSOLT, MDLRCY, FcoefETB, ETBheld, &
       & Zeffin, FSCBAL, FSCBEL, RhoETB, PCHRLN, PCHSLP, PCHWID, RhoPch1, RhoPch2, Tqt0, Tqp0, &
       & RDltRPn, ICCALC, ieqread, ipbtdir, DMAG0, RMAGMN, RMAGMX, oldmixa, iaxproj, rhoaxproj, &
       & Dnfac, MDBSETA, MDLNZD, MDLN0FL, FLlim, RLXN0, MDLSOL, MDLWBC, NANDER, NRETRY, MDLNCA, &
       & EpsHM, HPN, irktrc, DelRho, DelN, Q0, QA, MDITSN, MDITST, MDINTN, MDINTT, MDINTC, midbg, &
       & rMUb1, CPsi_old, PoyntI, WPDOT, WPPRE, TPRE, XOLD, bdf_history_ready, bdf_history_dt
  use tx_glob, only : T_TX_old, zIp_old, zIpol_old, VLOOP_old, VLOOPpol_old, PS2_old, &
       & WFEengy_old, WFMengy_old, WJengy_old, WFMpengy_old, WFMtengy_old
  use tx_savext, only : xreg, xreg_ia, xreg_ra
  use tx_extid, only : xid
  implicit none

  ! --- Inputs added after the base layout (the extension records 1-7 of BR-43 to BR-51) ---
  call xreg('MDLNSUP', MDLNSUP)
  call xreg('RWALL', RWALL)
  call xreg('RWALLZ', RWALLZ)
  call xreg('RWALLIZ', RWALLIZ)
  call xreg('FSTPTM', FSTPTM, NSM)
  call xreg('PROFD2C', PROFD2C)
  call xreg('PROFD2W', PROFD2W)
  call xreg('MDLSOLT', MDLSOLT)
  call xreg('DSOLT', DSOLT)
  call xreg('CHISOLT', CHISOLT, NSM)
  call xreg('RMUSOLT', RMUSOLT, NSM)
  call xreg('WSOLT', WSOLT)
  call xreg('MDLRCY', MDLRCY)

  ! --- State of the run: the ETB factor and whether it is held (MDLETB) ---
  call xreg('FCOEFETB', FcoefETB)
  call xreg('ETBHELD', ETBheld)
  ! --- State of the run: the form of psi, X(:,LQm4) = psi / rMUb2 (rMUb1 = 1, rMUb2 = rMU0
  !     once a beam was given in the process, else rMU0 and 1; TXMENU) ---
  call xreg('RMUB1', rMUb1)
  ! --- State of the run: the values of TXGLOB at the saved time that its next call takes
  !     differences and time integrals against (the Poynting flux PoyntI, the flux
  !     consumption CPsi, dW/dt); CPsi and VPoynt are in the base layout (TXLOAD: all
  !     of HISTITEMS or none) ---
  call xreg('CPSI_OLD', CPsi_old, 4)
  call xreg('POYNTI', PoyntI, 4)
  call xreg('WPDOT', WPDOT)
  call xreg('WPPRE', WPPRE)
  call xreg('TPRE', TPRE)
  call xreg('T_TX_OLD', T_TX_old)
  call xreg('ZIP_OLD', zIp_old)
  call xreg('ZIPOL_OLD', zIpol_old)
  call xreg('VLOOP_OLD', VLOOP_old)
  call xreg('VLOOPPOL_OLD', VLOOPpol_old)
  call xreg('PS2_OLD', PS2_old)
  call xreg('WFEENGY_OLD', WFEengy_old)
  call xreg('WFMENGY_OLD', WFMengy_old)
  call xreg('WJENGY_OLD', WJengy_old)
  call xreg('WFMPENGY_OLD', WFMpengy_old)
  call xreg('WFMTENGY_OLD', WFMtengy_old)
  ! --- State of the run: the history of the two-step BDF (IGBDF /= 0; TXLOOP), X one step
  !     before the saved one, the step it was taken with and whether it is there, so that the
  !     first step after LOAD is the one the run would have taken (TXLOAD: all three or none;
  !     XOLD has the length of the grid of the file) ---
  call xreg('BDF_READY', bdf_history_ready)
  call xreg('BDF_DT', bdf_history_dt)
  call xreg_ra('XOLD', XOLD, size(XOLD))
  ! --- External data of the process at the saved time, not carried by the file (tx_extid):
  !     NBI data read by F (kind: 0 none, 1-3 as ascii_input) and ripple data (F R) ---
  call xreg('EXT_NBI_KIND', xid%kind)
  call xreg('EXT_NBI_FILE', xid%fname)
  call xreg('EXT_NBI_NDATA', xid%ndata)
  call xreg('EXT_NBI_NOL', xid%nol)
  call xreg('EXT_NBI_PNBHEX', xid%totp)
  call xreg('EXT_NBI_HASH', xid%hash)
  call xreg('EXT_RIPPLE', xid%irpin)

  ! --- Inputs of the namelist that the base layout missed (BR-52) ---
  ! used in the time evolution (Zeffin also by TXPARM_CHECK at a LOAD: PNsSUP(2:3) of MDLNSUP = 1)
  call xreg('ZEFFIN', Zeffin)
  call xreg('FSCBAL', FSCBAL)
  call xreg('FSCBEL', FSCBEL)
  call xreg('RHOETB', RhoETB, size(RhoETB))
  call xreg('PCHRLN', PCHRLN)
  call xreg('PCHSLP', PCHSLP)
  call xreg('PCHWID', PCHWID)
  call xreg('RHOPCH1', RhoPch1)
  call xreg('RHOPCH2', RhoPch2)
  call xreg('TQT0', Tqt0)
  call xreg('TQP0', Tqp0)
  call xreg('RDLTRPN', RDltRPn)
  call xreg('ICCALC', ICCALC)
  call xreg('IEQREAD', ieqread)
  call xreg('IPBTDIR', ipbtdir)
  call xreg('DMAG0', DMAG0)
  call xreg('RMAGMN', RMAGMN)
  call xreg('RMAGMX', RMAGMX)
  call xreg('OLDMIXA', oldmixa)
  call xreg('IAXPROJ', iaxproj)
  call xreg('RHOAXPROJ', rhoaxproj)
  call xreg('DNFAC', Dnfac)
  call xreg('MDBSETA', MDBSETA)
  call xreg('MDLNZD', MDLNZD)
  call xreg('MDLN0FL', MDLN0FL)
  call xreg('FLLIM', FLlim)
  call xreg('RLXN0', RLXN0)
  call xreg('MDLSOL', MDLSOL)
  call xreg('MDLWBC', MDLWBC)
  call xreg('NANDER', NANDER)
  call xreg('NRETRY', NRETRY)
  call xreg('MDLNCA', MDLNCA)
  call xreg_ra('EPSHM', EpsHM, size(EpsHM))   ! rank 2, as a sequence
  call xreg_ia('HPN', HPN, size(HPN))
  ! used only at a RUN (initial profiles, equilibrium) or by a menu command (N: DelRho, DelN)
  call xreg('IRKTRC', irktrc)
  call xreg('DELRHO', DelRho)
  call xreg('DELN', DelN)
  call xreg('Q0', Q0)
  call xreg('QA', QA)
  call xreg('MDITSN', MDITSN)
  call xreg('MDITST', MDITST)
  call xreg('MDINTN', MDINTN)
  call xreg('MDINTT', MDINTT)
  call xreg('MDINTC', MDINTC)
  call xreg('MIDBG', midbg, size(midbg))

end subroutine savext_items

!***************************************************************
!
!   Write Data
!
!***************************************************************

subroutine TXWDAT
  use tx_commons, only : PI, RB, PN01V, NRMAX, WPT, Var
  use tx_interface, only : rLINEAVE
  use tx_core_module, only : intg_area

  implicit none
  real(8) :: rNbar

  !     ***** Volume-averaged density *****

  rNbar = intg_area(Var(:,1)%n) * 1.D20
  rNbar = rNbar / (PI * RB**2)

  write(6,'((1X,A," =",ES9.2,3(2X,A,"=",ES9.2)))') &
       &     'Ne(0)',    Var(0,1)%n,  &
       &     'UePhi(0)', Var(0,1)%Uph / 1.D3,  &
       &     'UiPhi(0)', Var(0,2)%Uph / 1.D3,  &
       &     'N0(RB)',   PN01V(NRMAX) * 1.D20, &
       &     'NB(0)',    rLINEAVE(0.D0)   / 1.D20,  &
       &     'NB(0.24)', rLINEAVE(0.24D0) / 1.D20,  &
       &     'NB(0.60)', rLINEAVE(0.6D0)  / 1.D20,  &
       &     'PF    ',   Var(0,1)%n * 1.D20 / rNbar
  write(6,'(1X,A," =",ES9.2,2(2X,A,"=",ES9.2))') &
       &     'Te(0)',    Var(0,1)%T,  &
       &     'Ti(0)   ', Var(0,2)%T,  &
       &     'Wst     ', WPT

end subroutine TXWDAT

!***************************************************************
!
!   Write Data 2
!
!***************************************************************

subroutine TXWDAT2

  use tx_graphic, only : GTY, NGT
  implicit none
  real(8) :: PNeMIN, PNeMAX, NB0MIN, NB0MAX, UiphMIN, UiphMAX

  PNeMIN  = minval(real(GTY(0:NGT+1, 1),8))
  PNeMAX  = maxval(real(GTY(0:NGT+1, 1),8))
  NB0MIN  = minval(real(GTY(0:NGT+1, 2),8))
  NB0MAX  = maxval(real(GTY(0:NGT+1, 2),8))
  UiphMIN = minval(real(GTY(0:NGT+1,11),8))
  UiphMAX = maxval(real(GTY(0:NGT+1,11),8))

  write(6,'((1X,A,"=",ES9.2,2(2X,A,"=",ES9.2)))') &
       &     'MAX(Ne(0))',     PNeMAX / 1.E20,  &
       &     'MAX(NB(0))',     NB0MAX / 1.E20,  &
       &     'MAX(UiPhi(0))', UiphMAX / 1.E3,   &
       &     'MIN(Ne(0))',     PNeMIN / 1.E20,  &
       &     'MIN(NB(0))',     NB0MIN / 1.E20,  &
       &     'MIN(UiPhi(0))', UiphMIN / 1.E3

end subroutine TXWDAT2

!***************************************************************
!
!   Write Statistic Data
!
!***************************************************************

subroutine TXSTAT
  use tx_commons, only : VOLAVN, ALI, VLOOP, TAUE1, TAUE2, TAUEP, TAUEH, BETAA, &
       &                 BETAPA, BETAN, Q, ANSAV, rIp, PI, RA, NRA, NRMAX, &
       &                 rMus, Chis, TSAV, Gamma_a, TAUPA, achg, &
       &                 rKeV, amas, amp, RR, rNuei, rho, Var, vlt
  use tx_glob, only : cal_flux, particle_source_balance
  use libitp, only: aitken2p
  implicit none
  integer(4) :: NR, NRL1, NRL2
  real(8) :: rhol1, rhol2, rmuil, chiil, uiphl, PTeVL, WDe, rNueiL
  character(len=21) :: cfmt = '(1X,2(A27,ES10.3,3X))'

  rhol1 = 0.3d0 ; rhol2 = 0.5d0
  do NR = 0, NRMAX-1
     if(rho(NR) <= rhol1 .and. rho(NR+1) >= rhol1) then
        NRL1 = NR
     else if(rho(NR) <= rhol2 .and. rho(NR+1) >= rhol2) then
        NRL2 = NR
        exit
     end if
  end do

  call cal_flux

  rmuil = aitken2p(rhol1,rmus(NRL1,2),rmus(NRL1+1,2),rmus(NRL1+2,2),rho(NRL1),rho(NRL1+1),rho(NRL1+2))
  chiil = aitken2p(rhol1,chis(NRL1,2),chis(NRL1+1,2),chis(NRL1+2,2),rho(NRL1),rho(NRL1+1),rho(NRL1+2))
  uiphl = aitken2p(rhol1,Var(NRL1,2)%Uph,Var(NRL1+1,2)%Uph,Var(NRL1+2,2)%Uph,rho(NRL1),rho(NRL1+1),rho(NRL1+2))

  ! For effective collision frequency at rho = 0.5

  PTeVL = aitken2p(rhol2,Var(NRL2,1)%T,Var(NRL2+1,1)%T,Var(NRL2+2,1)%T,rho(NRL2),rho(NRL2+1),rho(NRL2+2))
  ! WDe : curvature drift frequency
  WDe   = 2.d0 * sqrt(0.1d0) * sqrt(PTeVL * rKeV / (amas(2)*amp)) / RR
  rNueiL = aitken2p(rhol2,rNuei(NRL2),rNuei(NRL2+1),rNuei(NRL2+2),rho(NRL2),rho(NRL2+1),rho(NRL2+2))

  write(6,cfmt) "Vol. ave. of neutrality  = ", VOLAVN
  write(6,cfmt) "Inductance               = ", ALI, &
       &        "Loop voltage             = ", VLOOP
  write(6,cfmt) "Confinement time 1       = ", TAUE1, &
       &        "Confinement time 2       = ", TAUE2
  write(6,cfmt) "L-mode scaling time      = ", TAUEP, &
       &        "IPB98(y,2) scaling time  = ", TAUEH
  write(6,cfmt) "Beta                     = ", BETAA, &
       &        "Poloidal beta            = ", BETAPA
  write(6,cfmt) "Normalized beta          = ", BETAN, &
       &        "tau_p inside sep.        = ", TAUPA
  write(6,cfmt) "Ion flux thru sep.       = ", achg(2)*Var(NRA,2)%n*Var(NRA,2)%Ur*1.D20, &
       &        "Ion flux thru sep. (est) = ", Gamma_a
  write(6,cfmt) "Vol. averaged e density  = ", ANSAV(1), &
       &        "Greenwald density        = ", rIp / (PI * RA**2)
  write(6,cfmt) "Vol. averaged e temp.    = ", TSAV(1), &
       &        "Vol. averaged i temp.    = ", TSAV(2)
  write(6,cfmt) "Safety factor on axis    = ", Q(0), &
       &        "Safety factor at sep.    = ", Q(NRA)
  write(6,cfmt) "Ion Prandtl num. at 0.3  = ", rmuil/chiil, &
       &        "Ion tor. velocity at 0.3 = ", uiphl
  write(6,cfmt) "Eff .col. freq. at 0.5   = ", rNueiL/WDe, &
       &        "Plasma volume inside sep.= ", vlt(NRA)
  call particle_source_balance(6)

end subroutine TXSTAT

!***************************************************************
!
!   Steady state check
!
!***************************************************************

subroutine steady_check
  use tx_commons, only : NRMAX, Var, T_TX, rho
  use libitp, only: aitken2p

  implicit none
  integer(4) :: nr
  integer(4), save :: nrl = 0
  real(8) :: rhol, uiphl, pnevl, pevl
  real(8), save :: uiphl_old = 0.d0, pnevl_old = 0.d0, pevl_old = 0.d0

  ! Seel a grid number "nrl" nearest rho=0.3
  if(nrl == 0) then
     rhol = 0.3d0
     do NR = 0, NRMAX-1
        if(rho(NR) <= rhol .and. rho(NR+1) >= rhol) then
           NRL = NR
           exit
        end if
     end do
  end if
  
  ! Ion toroidal velocity
  uiphl = aitken2p(rhol,Var(nrl,2)%Uph,Var(nrl+1,2)%Uph,Var(nrl+2,2)%Uph,rho(nrl),rho(nrl+1),rho(nrl+2))

  ! Electron density
  PNeVl = aitken2p(rhol,Var(nrl,1)%n,  Var(nrl+1,1)%n,  Var(nrl+2,1)%n,  rho(nrl),rho(nrl+1),rho(nrl+2))

  ! Electron pressure
  PeVl  = aitken2p(rhol,Var(nrl,1)%p,  Var(nrl+1,1)%p,  Var(nrl+2,1)%p,  rho(nrl),rho(nrl+1),rho(nrl+2))

  ! Dispaly
  if(uiphl /= 0.d0) then
     write(6,*) real(t_tx),real(abs(uiphl - uiphl_old)/uiphl), &
          &                real(abs(pnevl - pnevl_old)/pnevl), &
          &                real(abs(pevl  - pevl_old )/pevl ), real(uiphl)
  end if
  uiphl_old = uiphl
  pnevl_old = pnevl
  pevl_old  = pevl

end subroutine steady_check

!***********************************************************
!
!  LINE AVERAGE OF rN
!
!***********************************************************

REAL(8) FUNCTION rLINEAVE(Rho)

  use tx_commons, only : RA, NRMAX, Var
  implicit none
  real(8), intent(IN) :: Rho
  integer(4) :: I, IR, NY = 100
  real(8) :: D, DY, Y, RL, SUML

  SUML = 0.D0
  D = Rho * RA
  DY = sqrt(RA*RA - D*D) / NY
  do I = 0, NY
     Y = DY * I
     RL = sqrt(Y*Y + D*D)
     IR = nint(RL * NRMAX / RA)
     SUML = SUML + Var(IR,1)%n * 1.D20 * DY
  end do
  rLINEAVE = SUML / sqrt(RA**2 - D**2)

end FUNCTION rLINEAVE

!***************************************************************
!
!   Save transport data
!
!***************************************************************

subroutine TXSAVE
  use tx_commons, only : &
       & SLID,stdin_ended,RA,rhob,rhoaccum,RR,BB,rbvt,ravl,rbvl,amas,achg,amb,achgb,rIPs,rIPe, &
       & PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza, &
       & PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0,PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB, &
       & Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC,FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1, &
       & FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV,FSCX,FSLC,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSNCOL, &
       & FSD01,FSD02,FSD03,FSD0z,FSG0iz,FSRP,FSNF, &
       & FSADV,FSADVB,FSUG,FSHL,MDLC,rLn,rLT,Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,&
       & PNBHT1,PNBHT2,rNRFe,RRFew,RRFe0,PRFHe,rNRFi,RRFiw,RRFi0,PRFHi,PNBCD, &
       & PN0s,V0,rGamm0_in,rGASPF,PN0zs,V0z,rGamm0zz_in,rGamm0iz_in,rGASPFz,PNsSUP,PTsSUP, &
       & DT,EPS,ADV,tiny_cap,CMESH0,WMESH0,CMESH,WMESH, &
       & ICMAX,NRMAX,NTMAX,NTSTEP,T_TX,TMAX,NT,NTCUM,NQMAX,IERR,X, &
       & NLCMAX,NCM,NTCOIL,DltRPn,m_pol,n_tor, &
       & MODEAV,IDIAG,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab, &
       & IGBDF,ISMTHD,iprestab,MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDANOM,iPoyntpol,iPoynttor, &
       & MDLETB,MDLNBD,MDBMCX,imodel_neo,PNBMPD,PNBPTC,thrp,kappa,CPsi,VPoynt,PoyntI,imodel_neo
  use tx_graphic, only : NGYTM,NGYVM,MODEG,MODEGL,NGT,NGVV,NGRSTP,NGTSTP,NGVSTP,GTY,GVY,GQY,GTX,GVX
  use tx_interface, only : TOUPPER
  use tx_savext, only : EXTSENT, EXTTAG, EXTVER, XM_COUNT, XM_WRITE, xmode, xio, xcount
  use tx_extid, only : xid, extid_process

  implicit none
  integer(4) :: IST, NQ, NR, NC, I, IGYT, IGYV, io
  character(len=100) :: TXFNAM, RCSId
  character(len=1) :: STR
  logical :: LEX

  RCSId = ' '

  do 
     write(6,*) '# INPUT : SAVE FILE NAME'
     flush(6)
     if (stdin_ended) return
     read(*,'(A100)',iostat=IST) TXFNAM
     if (IST > 0) then
        cycle
     else if (IST < 0) then
        stdin_ended = .true.
        return
     end if
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
        write(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',  &
             &              'ARE YOU SURE {Y/N} ?'
        flush(6)
        read(*,'(A1)') STR
        call TOUPPER(STR)
        if (STR == 'Y') then
           open(newunit=io,file=TXFNAM,iostat=IST,status='old',form='unformatted')
           if (IST == 0) then
              write(6,*) '# OLD FILE ( ', TXFNAM(1:len_trim(TXFNAM)), &
                   &     ' ) IS ASSIGNED FOR OUTPUT.'
              exit
           else if (IST > 0) then
              write(6,*) 'XX  OLD FILE open ERROR !, iostat = ', IST
           end if
        end if
     else
        open(newunit=io,file=TXFNAM,iostat=IST,status='new',form='unformatted')
        if (IST == 0) then
           write(6,*) '# NEW FILE ( ', TXFNAM(1:len_trim(TXFNAM)), &
                &     ' ) IS CREATED FOR OUTPUT.'
           exit
        else if (IST > 0) then
           write(6,*) 'XX  NEW FILE open ERROR !, iostat = ', IST
        end if
     end if
  end do

  ! The base layout below is followed by the extension (tx_savext, savext_items), which
  ! carries every input of the namelist that the base layout does not and the state of
  ! the run added since.

  write(io) SLID
  write(io) RCSId

  write(io) RA,rhob,rhoaccum,RR,BB,rbvt,ravl,rbvl
  write(io) amas,achg,amb,achgb,rIPs,rIPe
  write(io) PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza,PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0
  write(io) PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB
  write(io) Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC
  write(io) FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1,FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV
  write(io) FSCX,FSLC,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSNCOL,FSD01,FSD02,FSD03,FSD0z,FSG0iz,FSRP,FSNF
  write(io) FSADV,FSADVB,FSUG,FSHL
  write(io) rLn,rLT
  write(io) Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2
  write(io) rNRFe,RRFew,RRFe0,PRFHe,rNRFi,RRFiw,RRFi0,PRFHi,PNBCD,PNBMPD,PNBPTC
  write(io) PN0s,V0,rGamm0_in,rGASPF,PN0zs,V0z,rGamm0zz_in,rGamm0iz_in,rGASPFz,PNsSUP,PTsSUP
  write(io) DltRPn,kappa
  write(io) DT,EPS,ADV,tiny_cap,CMESH0,WMESH0,CMESH,WMESH
  write(io) ICMAX,NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  write(io) MODEG,MODEAV,MODEGL,IDIAG,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab,IGBDF
  write(io) ISMTHD,iprestab,MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDANOM,iPoyntpol,iPoynttor
  write(io) MDLETB,MDLNBD,MDBMCX,imodel_neo,MDLC,NTCOIL,m_pol,n_tor

  write(io) T_TX,TMAX,NT,NTCUM,NQMAX,IERR
  write(io) ((X(NR,NQ), NR=0, NRMAX), NQ=1, NQMAX)

  write(io) NGT
  write(io) (GTX(I), I=0, NGT)
  write(io) (GVX(I), I=0, NGVV)
  write(io) ((GTY(I,IGYT), I=0, NGT),  IGYT =1, NGYTM)
  write(io) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  write(io) (NLCMAX(NQ), NQ=1,NQMAX)
  write(io) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQMAX)
  write(io) (PoyntI(I), I=1, 2)
  write(io) (CPsi(I), I=0, 3)
  write(io) (VPoynt(I), I=0, 3)
  write(io) (thrp(I), I=1, 2*NRMAX)
  ! The extension: header, the items of savext_items, end record (tx_savext)
  xid = extid_process()
  xmode = XM_COUNT ; xcount = 0
  call savext_items
  write(io) EXTSENT, EXTTAG, EXTVER, xcount
  xmode = XM_WRITE ; xio = io
  call savext_items
  write(io) 'END             ', 0, 0
  xmode = 0
  close(io)
  write(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'

end SUBROUTINE TXSAVE

!***************************************************************
!
!   Load transport data
!
!***************************************************************

subroutine TXLOAD(IST)
  use tx_commons, only : &
       & allocate_txcomm, deallocate_txcomm, &
       & RA,RB,rhob,rhoaccum,RR,BB,rbvt,ravl,rbvl,amas,achg,amb,achgb,rIPs,rIPe, &
       & PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza, &
       & PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0,PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB, &
       & Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC,FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1, &
       & FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV,FSCX,FSLC,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSNCOL, &
       & FSD01,FSD02,FSD03,FSD0z,FSG0iz,FSRP,FSNF, &
       & FSADV,FSADVB,FSUG,FSHL,MDLC,rLn,rLT,Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP, &
       & PNBHT1,PNBHT2,rNRFe,RRFew,RRFe0,PRFHe,rNRFi,RRFiw,RRFi0,PRFHi,PNBCD, &
       & PN0s,V0,rGamm0_in,rGASPF,PN0zs,V0z,rGamm0zz_in,rGamm0iz_in,rGASPFz,PNsSUP,PTsSUP, &
       & DT,EPS,ADV,tiny_cap,CMESH0,WMESH0,CMESH,WMESH, &
       & ICMAX,NRMAX,NTMAX,NTSTEP,T_TX,TMAX,NT,NTCUM,NQMAX,IERR,X, &
       & NLCMAX,NCM,NTCOIL,DltRPn,m_pol,n_tor, &
       & MODEAV,IDIAG,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab, &
       & IGBDF,ISMTHD,iprestab,MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDANOM,iPoyntpol,iPoynttor, &
       & MDLETB,MDLNBD,MDBMCX,imodel_neo,PNBMPD,PNBPTC,rIP,thrp,kappa,CPsi,VPoynt,PoyntI, &
       & ErV,PTsV_FIX,PNsV_FIX,ErV_FIX, &
       & rMU0,rMUb1,rMUb2,NEMAX,ICONT,TAUE2,LQb1,Var,irestart, NSM,bdf_history_ready,bdf_history_dt,stdin_ended, &
       & PNsSUPlast, FcoefETB, FrdcETB, ETBheld, &
       & LQm4, ckt, vv, Pi
  use tx_init, only : TXCALM
  use tx_graphic, only : allocate_txgraf,deallocate_txgraf, &
       &                 MODEG,MODEGL,NGYTM,NGYVM,NGR,NGT,NGVV,NGRSTP,NGTSTP,NGVSTP, &
       &                 GTY,GVY,GQY,GTX,GVX
  use tx_variables
  use tx_coefficients, only : TXCALA
  use tx_parameter_control, only : TXPARM_CHECK
  use tx_glob, only : TXGLOB, glob_hist_loaded
  use mod_eqneo, only : wrap_eqneo
  use mod_cross_section, only : spline_table_carbon_rate_coef_adas, spline_table_beam_rate_coef
  use tx_savext, only : EXTSENT, EXTTAG, EXTVER, XM_READ, xmode, xclear, xread_items, xhas, xok, xreport
  use tx_extid, only : extid, xid, extid_loaded, extid_check
  implicit none
  integer(4), intent(out) :: IST
  integer(4) :: NQ, NR, NC, I, IGYT, IGYV, io
  character(len=100) ::  TXFNAM, RCSId
  character(len=8) :: LOADSLID, exttag_file
  integer(4) :: istx, iextver, ifirst, nitem, ipot_cur, ierr_file, ier_alloc
  real(8) :: rmub2inv, rmub1_proc
  logical :: LEX, etb_loaded, beam_form, beam_form_est
  ! the history of TXGLOB in the extension (savext_items): used only when all are there
  character(len=12), parameter :: HISTITEMS(16) = [character(len=12) :: 'CPSI_OLD', 'POYNTI', &
       & 'WPDOT', 'WPPRE', 'TPRE', 'T_TX_OLD', 'ZIP_OLD', 'ZIPOL_OLD', 'VLOOP_OLD', 'VLOOPPOL_OLD', &
       & 'PS2_OLD', 'WFEENGY_OLD', 'WFMENGY_OLD', 'WJENGY_OLD', 'WFMPENGY_OLD', 'WFMTENGY_OLD']

  do 
     write(6,*) '# INPUT : LOAD FILE NAME'
     flush(6)
     if (stdin_ended) return
     read(*,'(A100)',iostat=IST) TXFNAM
     if (IST > 0) then
        cycle
     else if (IST < 0) then
        stdin_ended = .true.
        return
     end if
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
        open(newunit=io,file=TXFNAM,iostat=IST,status='old',form='unformatted')
        if (IST == 0) then
           write(6,*) '# OLD FILE ( ', TXFNAM(1:len_trim(TXFNAM)),  &
                &     ' ) IS ASSIGNED FOR INPUT.'
           exit
        else if (IST > 0) then
           write(6,*) 'XX  OLD FILE open ERROR !, iostat = ', IST
        end if
     else
        write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
     end if
  end do

  read(io,iostat=IST) LOADSLID
  if (IST > 0) then
     write(6,*) 'XX read ERROR in TXLOAD !'
     close(io)
     return
  end if
  read(io) RCSId

  read(io) RA,rhob,rhoaccum,RR,BB,rbvt,ravl,rbvl
  read(io) amas,achg,amb,achgb,rIPs,rIPe
  read(io) PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza,PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0
  read(io) PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB
  read(io) Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC
  read(io) FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1,FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV
  read(io) FSCX,FSLC,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSNCOL,FSD01,FSD02,FSD03,FSD0z,FSG0iz,FSRP,FSNF
  read(io) FSADV,FSADVB,FSUG,FSHL
  read(io) rLn,rLT
  read(io) Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2
  read(io) rNRFe,RRFew,RRFe0,PRFHe,rNRFi,RRFiw,RRFi0,PRFHi,PNBCD,PNBMPD,PNBPTC
  read(io) PN0s,V0,rGamm0_in,rGASPF,PN0zs,V0z,rGamm0zz_in,rGamm0iz_in,rGASPFz,PNsSUP,PTsSUP
  read(io) DltRPn,kappa
  read(io) DT,EPS,ADV,tiny_cap,CMESH0,WMESH0,CMESH,WMESH
  read(io) ICMAX,NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  read(io) MODEG,MODEAV,MODEGL,IDIAG,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab,IGBDF
  read(io) ISMTHD,iprestab,MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDANOM,iPoyntpol,iPoynttor
  ! imodel_neo (with the potato closure imodel_neo(4)) travels in this record: a file keeps the
  ! closure it was run with (0, the matrix inversion without it, in files written before the closure).
  ipot_cur = imodel_neo(4)
  read(io) MDLETB,MDLNBD,MDBMCX,imodel_neo,MDLC,NTCOIL,m_pol,n_tor
  if( imodel_neo(4) /= ipot_cur ) write(6,'(A,I2,A,I2,A)') &
       & ' ## Potato closure from the file: imodel_neo(4) =', imodel_neo(4), ' (was', ipot_cur, ')'

  ! The error flag of the run that wrote the file: a file written after a failed step
  ! holds the last accepted step, and a LOAD clears the flag (TXMENU), as it always
  ! did in a new process.  It was read into IERR, which then also took the status of
  ! the allocation; with the arrays already allocated that status was not set, so
  ! such a file made a LOAD in the same process free every array and stop.
  read(io) T_TX,TMAX,NT,NTCUM,NQMAX,ierr_file
  if(ierr_file /= 0) write(6,'(A,I3,A)') ' ## The file was written after a failed step (IERR =', &
       & ierr_file, '); it holds the last accepted step'

  call allocate_txcomm(ier_alloc)
  if(ier_alloc == 0) call allocate_txgraf(ier_alloc)
  if(ier_alloc /= 0) then
     call deallocate_txcomm
     call deallocate_txgraf
     write(6,*) "XX Allocation error : TXLOAD"
     close(io)
     IST = 1
     return
  end if

  read(io) ((X(NR,NQ), NR=0, NRMAX), NQ=1, NQMAX)
  ! The two-step BDF history is taken from the extension below if the file has it; else the
  ! BDF restarts from a backward-Euler step (TXCALB).
  bdf_history_ready = .false.

  read(io) NGT
  read(io) (GTX(I), I=0, NGT)
  ! The GV history is not restored (NGVV = -1 below) and its length is not in the
  ! file, so its two records are skipped.  Read with the NGVV of a run still in
  ! memory, they would ask for more data than the file has after that run went
  ! beyond the saved time.
  NGVV = -1
  read(io) (GVX(I), I=0, NGVV)
  read(io) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
  read(io) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  read(io) (NLCMAX(NQ), NQ=1,NQMAX)
  read(io) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQMAX)
  read(io) (PoyntI(I), I=1, 2)
  read(io) (CPsi(I), I=0, 3)
  read(io) (VPoynt(I), I=0, 3)
  read(io) (thrp(I), I=1, 2*NRMAX)
  ! The extension (tx_savext): every input and state added after the base layout, by name.
  ! A file of the base layout only (before the extension records) ends here; a file of a
  ! development build of BR-43 to BR-51 has its extension records 1-7 here, which are not
  ! read.  Items not in the file keep their current values (reported by xreport).
  call xclear
  read(io,iostat=istx) ifirst
  if( istx < 0 ) then
     write(6,'(A)') ' ## The file has the base layout only (written before the extension records)'
  else if( istx == 0 .and. ifirst == EXTSENT ) then
     backspace(io)
     read(io,iostat=istx) ifirst, exttag_file, iextver, nitem
     if( istx == 0 .and. exttag_file == EXTTAG ) then
        call xread_items(io, nitem, istx)
        if( istx /= 0 ) write(6,'(A)') &
             & ' ## WARNING: the extension of the file ends early; the items read before are used'
        if( iextver > EXTVER ) write(6,'(A,I3,A,I3,A)') ' ## The extension of the file is of version', &
             & iextver, ' (this build:', EXTVER, '); the items are taken by name'
     else
        write(6,'(A)') ' ## WARNING: the extension of the file is not readable; it is not read'
     end if
  else if( istx == 0 .and. ifirst >= 1 .and. ifirst <= 7 ) then
     write(6,'(A)') ' ## The file has the extension records 1-7 of a development build (BR-43 to BR-51);' &
          & //' they are not read'
  else
     write(6,'(A)') ' ## WARNING: the records after the base layout of the file are not recognized; not read'
  end if
  rmub1_proc = rMUb1
  xid = extid()
  xmode = XM_READ
  call savext_items
  xmode = 0
  call xreport
  call extid_loaded(xhas('EXT_NBI_KIND'), PNBH - (PNBHP + PNBHT1 + PNBHT2))
  etb_loaded = xhas('FCOEFETB') .and. xhas('ETBHELD')
  glob_hist_loaded = all([(xok(trim(HISTITEMS(i))), i = 1, size(HISTITEMS))])
  ! The two-step BDF history: all of it or none.  A file without it (written before
  ! 2026-10-14) restarts from a backward-Euler step, as before.
  if( .not. (xok('XOLD') .and. xok('BDF_READY') .and. xok('BDF_DT')) ) bdf_history_ready = .false.
  if( IGBDF /= 0 ) then
     if( bdf_history_ready ) then
        write(6,'(A,ES11.4,A)') ' ## Two-step BDF history from the file (its step', bdf_history_dt, &
             & '): the next step continues it'
     else
        write(6,'(A)') ' ## No two-step BDF history in the file: the first step after LOAD is backward Euler'
     end if
  end if
  if( .not. etb_loaded ) then
     ! A file without the ETB factor: from the step count of the run that wrote it, as
     ! TXCALV forms it (held if that run reached its last step).  A
     ! factor held by an earlier run of a run that then stopped early (NT < NTMAX)
     ! cannot be told from the file and comes out unheld.
     ETBheld = NTMAX > 0 .and. NT >= NTMAX
     if( ETBheld ) then
        FcoefETB = FrdcETB
     else if( NTMAX > 0 ) then
        FcoefETB = 1.d0 - NT * ((1.d0 - FrdcETB) / NTMAX)
     else
        FcoefETB = 1.d0
     end if
     if( MDLETB /= 0 ) write(6,'(A,L2,F9.5)') &
          & ' ## No record of the ETB factor in the file: from NT and NTMAX, held, factor =', ETBheld, FcoefETB
  end if
  close(io)
  ! The loaded PNsSUP are the values last set, so TXPARM_CHECK does not take PNsSUP(2:3)
  ! of the file for explicit input; with MDLNSUP = 1 it derives them again from PNsSUP(1).
  PNsSUPlast(:) = PNsSUP(:)
  write(6,'(2A)') '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE. ID = ',LOADSLID

  rb = rhob * ra

  rIP=rIPs
  ICONT = 1
  irestart = 1

  call TXPARM_CHECK

  !  Read spline table for rate coefficients
  call spline_table_carbon_rate_coef_adas
  call spline_table_beam_rate_coef

  NEMAX = NRMAX
  ! Create the mesh and read an equilibrium
  call TXCALM

  ! X(:,LQm4) is psi / rMUb2 in the form of the run that wrote the file: rMUb1 = 1,
  ! rMUb2 = rMU0 once a beam was given in that process (TXMENU), rMU0 and 1 otherwise.
  ! The extension of the file records it (RMUB1).  For a file without the record it is
  ! read from psi, as before: the boundary condition of Ampere's law,
  ! dpsi/dV = 2 pi rMU0 Ip / ckt at the wall, makes the edge gradient of X(:,LQm4) over
  ! that value 1 / rMUb2, about 1 or 1 / rMU0 (the beam inputs of the file do not tell
  ! the form: a beam switched off later leaves it).  With the record this estimate is a
  ! check, and a disagreement is reported.
  if(rIP /= 0.d0) then
     rmub2inv = ckt(NRMAX) * (X(NRMAX,LQm4) - X(NRMAX-1,LQm4)) / (vv(NRMAX) - vv(NRMAX-1)) &
          &   / (2.d0 * Pi * rMU0 * rIP * 1.d6)
     beam_form_est = abs(rmub2inv) > 1.d0 / sqrt(rMU0)
  else
     beam_form_est = rmub1_proc == 1.d0 .or. maxval(X(:,LQb1)) > epsilon(1.d0)
  end if
  if(xhas('RMUB1')) then
     beam_form = rMUb1 == 1.d0
     if(beam_form .neqv. beam_form_est) write(6,'(A,L2,A)') &
          & ' ## WARNING: the form of psi recorded in the file (with a beam:', beam_form, &
          & ') is not the one its edge gradient gives; the record is used'
  else
     beam_form = beam_form_est
  end if
  if(beam_form) then
     if(rmub1_proc /= 1.d0) write(6,'(A)') ' ## psi of the file in the form with a beam: rMUb1 = 1, rMUb2 = rMU0'
     rMUb1 = 1.D0
     rMUb2 = rMU0
  else
     if(rmub1_proc /= rMU0) write(6,'(A)') ' ## psi of the file in the form without a beam: rMUb1 = rMU0, rMUb2 = 1'
     rMUb1 = rMU0
     rMUb2 = 1.D0
  end if

  call TXCALV(X,0)

  do i = 1, NSM
     PNsV_FIX(:,i) = Var(:,i)%n
     PTsV_FIX(:,i) = Var(:,i)%T
  end do
  ErV_FIX (:) = ErV (:)

  call wrap_eqneo

  call TXCALC(0)
  call TXCALA

  ! TAUE2, PoyntI and VPoynt of the LOAD point follow the first-point rule of TXGLOB
  ! (irestart = 1): no previous step; with the history of the extension (glob_hist_loaded)
  ! they, CPsi and the values the next step takes differences against are those of the
  ! run at the saved time.  (TAUE2 was set here from GTY(NGT,34), which is WST(2), and
  ! TXGLOB overwrote it.)

  ! Reset start point of graphics
  NGT=-1
  NGR=-1
  NGVV=-1

  call TXGLOB
  glob_hist_loaded = .false.
  call TXWDAT
  call TXWDAT2

  irestart = 0

  ! The external data of the process against those of the run that wrote the file
  call extid_check(1)

end subroutine TXLOAD

!***************************************************************
!
!   Save graphic data
!
!***************************************************************

subroutine TXGSAV

  use tx_commons, only : &
       & SLID,RA,rhob,RR,BB,rbvt,ravl,rbvl,amas,achg,amb,achgb,PTe0,PTea,PTi0,PTia,PTz0,PTza, &
       & Dfs0,rMus0,WPM0,Chis0,FSDFIX,FSANOM,FSCBKP,FSCBSH, &
       & FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV,PROFD,PROFC, &
       & FSCX,FSLC,FSRP,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSNCOL, &
       & FSD01,FSD02,FSD03,FSD0z,FSG0iz,Ebmax,esps,PNBH,PNBHP,PNBHT1,PNBHT2,PRFHe,PRFHi,PNBCD,PNBMPD,PNBPTC, &
       & V0,rGamm0,rGASPF,V0z,rGamm0zz,rGamm0iz,rGASPFz,PNsSUP,PTsSUP, &
       & DT,NRMAX,NTMAX,NTSTEP,rG1,T_TX,TMAX,NT,NTCUM,NQMAX,IERR,X, &
       & NLCMAX,NCM,DltRPn,thrp,kappa,CPsi,VPoynt,PoyntI
  use tx_graphic, only : MODEG,MODEGL,NGYRM,NGYTM,NGYVM,NGR,NGT,NGVV,NGRSTP,NGTSTP,NGVSTP, &
       &                 GTY,GVY,GQY,GY,GYT,GTX,GVX,GT
  implicit none
  integer(4) :: IST, NQ, NR, NC, IGR, I, IGYR, IGYT, IGYV, io
  character(len=100) :: TXFNAM, RCSId
  character(len=1) :: STR
  logical :: LEX

  RCSId = ' '

  do
     write(6,*) '# INPUT : SAVE FILE NAME'
     flush(6)
     read(*,'(A100)',iostat=IST) TXFNAM
     if (IST > 0) then
        cycle
     else if(IST < 0) then
        return
     end if
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
        write(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',  &
             &              'ARE YOU SURE {Y/N} ?'
        flush(6)
        read(*,'(A1)') STR
        if (STR /= 'Y' .and. STR /= 'y') then
        else
           open(newunit=io,file=TXFNAM,iostat=IST,status='old',form='unformatted')
           if (IST == 0) then
              write(6,*) '# OLD FILE ( ', TXFNAM(1:len_trim(TXFNAM)), &
                   &     ' ) IS ASSIGNED FOR OUTPUT.'
              exit
           else if (IST > 0) then
              write(6,*) 'XX  OLD FILE open ERROR !, iostat = ', IST
           end if
        end if
     else
        open(newunit=io,file=TXFNAM,iostat=IST,status='new',form='unformatted')
        if (IST == 0) then
           write(6,*) '# NEW FILE ( ', TXFNAM(1:len_trim(TXFNAM)), &
                &     ' ) IS CREATED FOR OUTPUT.'
           exit
        else if (IST > 0) then
           write(6,*) 'XX  NEW FILE open ERROR !, iostat = ', IST
        end if
     end if
  end do

    write(io) SLID
!!$    write(io) RCSId

  write(io) RA,rhob,RR,BB,rbvt,ravl,rbvl
  write(io) amas,achg,amb,achgb
  write(io) PTe0,PTea,PTi0,PTia,PTz0,PTza
  write(io) Dfs0,rMus0,WPM0,Chis0
  write(io) FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1,FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV
  write(io) PROFD,PROFC
  write(io) FSCX,FSLC,FSRP,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSNCOL,FSD01,FSD02,FSD03,FSD0z,FSG0iz
  write(io) Ebmax,esps,PNBH,PNBHP,PNBHT1,PNBHT2
  write(io) PRFHe,PRFHi,PNBCD,PNBMPD,PNBPTC
  write(io) V0,rGamm0,rGASPF,V0z,rGamm0zz,rGamm0iz,rGASPFz,PNsSUP,PTsSUP
  write(io) DltRPn,kappa
  write(io) DT
  write(io) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  write(io) MODEG,MODEGL

  write(io) T_TX,TMAX,NT,NTCUM,NQMAX,IERR
  write(io) ((X(NR,NQ), NR=0, NRMAX), NQ=1, NQMAX)

  write(io) NGR,NGT,NGVV
  write(io) (GT(IGR), IGR=0, NGR)
  write(io)(((GY%v(I,IGR,IGYR), I=0,NRMAX), IGR=0,NGR), IGYR=1,NGYRM)
  write(io) (GY%gnrm(IGYR), IGYR=1,NGYRM)
  write(io) (GTX(I), I=0, NGT)
  write(io) (GVX(I), I=0, NGVV)
  write(io) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
  write(io) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  write(io) (NLCMAX(NQ), NQ=1,NQMAX)
  write(io) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQMAX)
  write(io) (((GYT%v(NR,I,IGYR), NR=0,NRMAX), I=0,NGT), IGYR=1,NGYRM)
  write(io) (PoyntI(I), I=1, 2)
  write(io) (CPsi(I), I=0, 3)
  write(io) (VPoynt(I), I=0, 3)
  write(io) (thrp(I), I=1, 2*NRMAX)
  close(io)
  write(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'

end subroutine TXGSAV

!***************************************************************
!
!   Load graphic data
!
!***************************************************************

subroutine TXGLOD(IST)

  use tx_commons, only : &
       & allocate_txcomm, deallocate_txcomm, &
       & RA,RB,rhob,RR,BB,rbvt,ravl,rbvl,amas,achg,amb,achgb,PTe0,PTea,PTi0,PTia,PTz0,PTza, &
       & Dfs0,rMus0,WPM0,Chis0,FSDFIX,FSANOM,FSCBKP,FSCBSH, &
       & FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV,PROFD,PROFC, &
       & FSCX,FSLC,FSRP,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSNCOL, &
       & FSD01,FSD02,FSD03,FSD0z,FSG0iz,Ebmax,esps,PNBH,PNBHP,PNBHT1,PNBHT2,PRFHe,PRFHi,PNBCD,PNBMPD,PNBPTC, &
       & V0,rGamm0,rGASPF,V0z,rGamm0zz,rGamm0iz,rGASPFz,PNsSUP,PTsSUP, &
       & DT,NRMAX,NTMAX,NTSTEP,rG1,T_TX,TMAX,NT,NTCUM,NQMAX,IERR,X, &
       & NLCMAX,NCM,DltRPn,thrp,kappa,CPsi,VPoynt,PoyntI!,rho
  use tx_init, only : TXCALM
  use tx_graphic, only : allocate_txgraf,deallocate_txgraf, &
       &                 MODEG,MODEGL,NGYRM,NGYTM,NGYVM,NGR,NGT,NGVV,NGRSTP,NGTSTP,NGVSTP, &
       &                 GTY,GVY,GQY,GY,GYT,GTX,GVX,GT
  use tx_variables

  implicit none
  integer(4), intent(out) :: IST
  integer(4) :: NQ, NR, NC, IGR, I, IGYR, IGYT, IGYV, io
  character(len=100) :: TXFNAM
!!  character(len=100) :: RCSId
  character(len=8) :: LOADSLID
  logical :: LEX

  do 
     write(6,*) '# INPUT : LOAD FILE NAME'
     flush(6)
     read(*,'(A100)',iostat=IST) TXFNAM
     if (IST > 0) then
        cycle
     else if(IST < 0) then
        return
     end if
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
        open(newunit=io,file=TXFNAM,iostat=IST,status='old',form='unformatted')
        if (IST == 0) then
           write(6,*) '# OLD FILE ( ', TXFNAM(1:len_trim(TXFNAM)),&
                &     ' ) IS ASSIGNED FOR INPUT.'
           exit
        else if (IST > 0) then
           write(6,*) 'XX  OLD FILE open ERROR !, iostat = ', IST
        end if
     else
        write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
     end if
  end do

  read(io,iostat=IST) LOADSLID
  if (IST > 0) then
     write(6,*) 'XX read ERROR in TXGLOD !'
     close(io)
     return
  end if
!!$    !  if(LOADSLID(1:5) == 'tx459') then
!!$    read(io) RCSId

  read(io) RA,rhob,RR,BB,rbvt,ravl,rbvl
  read(io) amas,achg,amb,achgb
  read(io) PTe0,PTea,PTi0,PTia,PTz0,PTza
  read(io) Dfs0,rMus0,WPM0,Chis0
  read(io) FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1,FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV
  read(io) PROFD,PROFC
  read(io) FSCX,FSLC,FSRP,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSNCOL,FSD01,FSD02,FSD03,FSD0z,FSG0iz
  read(io) Ebmax,esps,PNBH,PNBHP,PNBHT1,PNBHT2
  read(io) PRFHe,PRFHi,PNBCD,PNBMPD,PNBPTC
  read(io) V0,rGamm0,rGASPF,V0z,rGamm0zz,rGamm0iz,rGASPFz,PNsSUP,PTsSUP
  read(io) DltRPn,kappa
  read(io) DT
  read(io) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  read(io) MODEG,MODEGL

  read(io) T_TX,TMAX,NT,NTCUM,NQMAX,IERR

  call allocate_txcomm(ierr)
  call allocate_txgraf(ierr)
  if(ierr /= 0) then
     call deallocate_txcomm
     call deallocate_txgraf
     write(6,*) "XX Allocation error : TXGLOD"
  end if

  read(io) ((X(NR,NQ), NR=0, NRMAX), NQ=1, NQMAX)

  read(io) NGR,NGT,NGVV
  read(io) (GT(IGR), IGR=0, NGR)
  read(io) (((GY%v(I,IGR,IGYR), I=0, NRMAX), IGR=0, NGR), IGYR=1, NGYRM)
  read(io) (GY%gnrm(IGYR), IGYR=1,NGYRM)
  read(io) (GTX(I), I=0, NGT)
  read(io) (GVX(I), I=0, NGVV)
  read(io) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
  read(io) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  read(io) (NLCMAX(NQ), NQ=1,NQMAX)
  read(io) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQMAX)
  read(io) (((GYT%v(NR,I,IGYR), NR=0,NRMAX), I=0,NGT), IGYR=1,NGYRM)
  read(io) (PoyntI(I), I=1, 2)
  read(io) (CPsi(I), I=0, 3)
  read(io) (VPoynt(I), I=0, 3)
  read(io) (thrp(I), I=1, 2*NRMAX)
  !  end if
  close(io)
  write(6,'(2A)') '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE. ID = ',LOADSLID

  GYT%gnrm(:) = GY%gnrm(:)

  rb = rhob * ra

  call TXCALM
!!$    call TXCALV(X)
!!$    call TXCALC(0)
!!$    call TXGLOB
!!$    call TXWDAT
!!$    call TXWDAT2

!!$  do nr = 0, nrmax
!!$     write(6,*) real(RHO(NR)), GY%v(NR,NGR,1), GY%v(NR,NGR,37), GY%v(NR,NGR,14)*1.e3, GY%v(NR,NGR,15)*1.e3, GY%v(NR,NGR,35), GY%v(NR,NGR,36)
!!$  end do

  RETURN
end SUBROUTINE TXGLOD

!***************************************************************
!
!   Read ASCII data file to substitute it into array
!
!***************************************************************

subroutine ascii_input

  use tx_commons, only : infiles, nmax_file, n_infiles, iflag_file, datatype, infile_name
  use tx_interface, only : KSPLIT_TX
  implicit none
  integer(4) :: IST, i, j, k, nrho, i_start, jshift, io
  character(len=100) :: TXFNAM
  character(len=180) :: kline, kline1, kline2
  character(len=20)  :: kmesh, kdata
  integer(4) :: nol, itype, kfile
  LOGICAL :: LEX

  ! The kind of the data (1-3) is read into kfile.  iflag_file, the kind of the external
  ! data the process holds, is set only when the data have been read: an answer whose file
  ! does not exist, or an end of input before the data, leaves the external data of the
  ! process as they were (kind and data); a failure after the data held before were freed
  ! leaves the process without external data (iflag_file = 0).
  ! Check a mode of input file
  do 
     write(6,*) '# NBI input from OFMC (1)(2) or arbitrary input (3) ?'
     flush(6)
     read(*,'(I1)',iostat=ist) kfile
     if(ist > 0) then
        cycle
     else if(ist < 0) then
        return
     end if
     if(kfile >= 1 .and. kfile <= 3) then
        exit
     else
        cycle
     end if
  end do

  ! *** Pre-defined input (OFMC: OrbitEffectDist.dat) ***************
  if(kfile == 1) then
     TXFNAM = 'OrbitEffectDist.dat'
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
     else
        write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
        write(6,'(A)') ' ## Nothing read: the external data of the process are as they were'
        return
     end if

     n_infiles = 8

  ! *** Pre-defined input (OFMC: Torque.txt) ************************
  else if(kfile == 2) then
     TXFNAM = 'Torque.txt'
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
     else
        write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
        write(6,'(A)') ' ## Nothing read: the external data of the process are as they were'
        return
     end if

     n_infiles = 1 ! only the TOTAL torque

  ! *** Arbitrary input *********************************************
  else if(kfile == 3) then
     ! Check ASCII file
     do
        write(6,*) '# INPUT : LOAD ASCII FILE NAME'
        flush(6)
        read(*,'(A100)',iostat=IST) TXFNAM
        if (IST > 0) then
           cycle
        else if (IST < 0) then
           return
        end if
        inquire(file=TXFNAM,exist=LEX)
        if (LEX) then
           exit
        else
           write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
        end if
     end do

     ! Number of data
     do
        write(6,*) '# Number of data which you would like to use as inputs ?'
        flush(6)
        read(*,'(I2)',iostat=ist) n_infiles
        if(ist > 0) then
           cycle
        else if(ist < 0) then
           return
        else
           exit   
        end if
     end do
  end if

  ! Deallocate if already allocated
  if(allocated(infiles)) deallocate(infiles)

  ! Allocate derived type
  allocate(infiles(1:n_infiles))

  ! *** Pre-defined input (OFMC: OrbitEffectDist.dat) ***************
  if(kfile == 1) then
     ! Read data
     open(newunit=io,file=TXFNAM,iostat=IST,status='old',form='formatted')
     if (IST == 0) then
        write(6,*) '# ASCII FILE ( ', TXFNAM(1:len_trim(TXFNAM)),  &
             &     ' ) IS ASSIGNED FOR INPUT.'
     else if (IST > 0) then
        write(6,*) 'XX  ASCII FILE open ERROR !, iostat = ', IST
     else if (IST < 0) then
        deallocate(infiles) ; iflag_file = 0 ; infile_name = ' '
        return
     end if

     ! Name
     !   In case of "NBI input from OFMC (1)", sequence of data is already
     !   defined as follows:
     !     (1) S_birth_ele,  (2) S_birth_tot, (3) S_birth_trap, (4) S_birth_pass,
     !     (5) S_birth_loss, (6) S_orbit_tot, (7) S_orbit_trap, (8) S_orbit_pass
     !   Therefore name of the data is not necessary; hence null is substituted.
     !   Note: The following relations are satisfied:
     !     (1) = (3) + (4) + (5), (2) = (3) + (4), (6) = (7) + (8)
     !     (3) /= (7), (4) /= (8), (2) = (6)
     infiles(1:n_infiles)%name = ' '

     ! Read data from the file
     i    = 0  ! outermost loop count
     nol  = 0  ! number of lines
     nrho = 0  ! number of times that "RHO" appears in kline, i.e. separator
     k    = 0  ! data is taken during k/=0
     do 
        read(io,'(A)',iostat=IST) kline
        if(ist > 0) then
           stop 'Read file error'
        else if(ist < 0) then
           exit ! detect the end of the file
        end if
        i = i + 1
        if(k == 0) then
           if(index(kline,"RHO") /= 0) then ! detect the start position of the data chunk
              nrho = nrho + 1
              k = 1
              i_start = i + 1

              jshift = 2 ! denotes the index "I" and the radial coordinates "RHO"
              if(index(kline,"PSI") /= 0) then ! detect whether the file is new or old
                 !== jshift ================================================================
                 !  Old type "OrbitEffctDist.dat" includes only the "RHO" coordinate.
                 !  New type "OrbitEffctDist.dat" includes the "PSI" and "RHO" coordinates.
                 !==========================================================================
                 jshift = 3 ! denotes the radial coordinates "PSI", which will be discarded
                            !   as well as "I" and "RHO"
              end if
           end if
           cycle
        end if

        if(nrho == 1) then ! === Get "TOTAL" and "TOTAL POWER" data ===
           if(index(kline,"sec") /= 0) then
              ! Total number of ions per second
              j = 0
              do
                 kline = trim(adjustl(kline))
                 call ksplit_tx(kline,' ',kline1,kline2)
                 if(index(kline1,"E+") == 0) then
                    kline = kline2
                    cycle
                 else
                    j = j + 1
                    kdata = trim(kline1)
                    read(kdata,'(E15.7)') infiles(j)%totS
                    kline = kline2
                 end if
                 if(len_trim(kline2) == 0) exit ! all the data has been already taken.
              end do

           else if(index(kline,"W") /= 0) then
              ! Total power of ions
              j = 0
              do
                 kline = trim(adjustl(kline))
                 call ksplit_tx(kline,' ',kline1,kline2)
                 if(index(kline1,"E+") == 0) then
                    kline = kline2
                    cycle
                 else
                    j = j + 1
                    kdata = trim(kline1)
                    read(kdata,'(E15.7)') infiles(j)%totP
                    kline = kline2
                 end if
                 if(len_trim(kline2) == 0) then
                    exit ! all the data has been already taken.
                 end if
              end do
              
              k = 0 ! flag for the next data chunk
           end if

       else if(nrho == 2) then ! === Get "SMOOTHING" data ===
           if(len_trim(kline) == 0) then ! detect the end of the data chunk
              k = 0
           end if
           nol = i - i_start + 1

           do j = 1, n_infiles + jshift
              ! jshift denotes the index number and the mesh data columns
              kline = trim(adjustl(kline))
              call ksplit_tx(kline,' ',kline1,kline2)
              if(j == 1) then
                 ! First line (index I) is discarded.
                 kline = kline2
                 cycle
              else if(j == jshift) then
                 ! Second (jshift == 0) or third (jshift == 1) line is defined as a mesh data, rho.
                 kmesh = trim(kline1)
                 kline = kline2
              else if(j > jshift) then
                 ! Third to eighth lines are defined as data, in case of jshift == 0
                 ! Fourth to ninth lines are defined as data, in case of jshift == 1
                 ! Please see "Name" shown above
                 kdata = trim(kline1)
                 read(kmesh,'(E15.7)') infiles(j-jshift)%r(nol)
                 read(kdata,'(E15.7)') infiles(j-jshift)%data(nol)
                 kline = kline2
              else
                 kline = kline2
              end if
           end do

        else if(nrho == 4) then ! === Get "V//" data ===
           if(len_trim(kline) == 0) then ! detect the end of the data chunk
              k = 0
           else
              nol = i - i_start + 1
              j = 5 ! Vb data exist for (6)(7)(8) only
              infiles(1:j)%vb(nol) = 0.d0
              do
                 kline = trim(adjustl(kline))
                 call ksplit_tx(kline,' ',kline1,kline2)
                 if(index(kline1,"E+") == 0) then
                    kline = kline2
                    cycle
                 else
                    j = j + 1
                    ! Parallel velocity
                    kdata = trim(kline1)
                    read(kdata,'(E15.7)') infiles(j)%vb(nol)
                    kline = kline2
                 end if
                 if(len_trim(kline2) == 0) exit ! all the data has been already taken.
              end do
           end if

        else
           ! detect the end of the data chunk
           if(len_trim(kline) == 0) k = 0

        end if

     end do

     if(nol > nmax_file) then
        write(6,*) 'XX The number of lines in the file exceeds ',nmax_file
        write(6,*) '   Program is terminated.'
        stop
     else
        ! All the data should have same mesh data.
        infiles(:)%nol = nol
     end if

     close(io)

!!$     write(6,*) "===== data ====="
!!$     do i = 1, nol
!!$        write(6,'(I3,F6.3,8(ES12.4))') i,infiles(1)%r(i),(infiles(j)%data(i),j=1,8)
!!$     end do
!!$     write(6,*) "===== V// ====="
!!$     do i = 1, nol
!!$        write(6,'(I3,F6.3,8(ES12.4))') i,infiles(1)%r(i),(infiles(j)%vb(i),j=1,8)
!!$     end do
!!$     write(6,*) "===== total S ====="
!!$     write(6,'(9X,8(ES12.4))') (infiles(i)%totS, i = 1, 8)
!!$     write(6,*) "===== total W ====="
!!$     write(6,'(9X,8(ES12.4))') (infiles(i)%totP, i = 1, 8)
!!$     stop

  ! *** Pre-defined input (OFMC: Torque.txt) ************************
  else if(kfile == 2) then
     ! Read data
     open(newunit=io,file=TXFNAM,iostat=IST,status='old',form='formatted',position='rewind')
     if (IST == 0) then
        write(6,*) '# ASCII FILE ( ', TXFNAM(1:len_trim(TXFNAM)),  &
             &     ' ) IS ASSIGNED FOR INPUT.'
     else if (IST > 0) then
        write(6,*) 'XX  ASCII FILE open ERROR !, iostat = ', IST
     else if (IST < 0) then
        deallocate(infiles) ; iflag_file = 0 ; infile_name = ' '
        return
     end if

     infiles(1:n_infiles)%name = 'TOTAL'

     k  = 0  ! data is taken during k/=0
     do 
        read(io,'(A)',iostat=IST) kline
        if(ist > 0) then
           stop 'Read file error'
        else if(ist < 0) then
           exit ! detect the end of the file
        end if
        
        if(k == 0) then
           if(index(kline,"I,") /= 0) then ! detect the start position of the data chunk
              k = 1
              nol = 0
           end if
           cycle
        end if

        j = 0
        do
           kline = trim(adjustl(kline))
           if(kline == ' ') then
              k = 0
              exit
           end if
           call ksplit_tx(kline,',',kline1,kline2)
           if(index(kline1,"E") == 0) then
              kdata = trim(kline1)
              if(kdata == ' ' .or. kdata == "TOTAL") then
                 nol = 0
              else
                 nol = nol + 1 ! Not use "I" index
!                 read(kdata,'(I3)') nol ! Use "I" index
              end if
              kline = kline2
              cycle
           else
              j = j + 1
              kdata = trim(kline1)
!              write(6,*) "nol",nol,kdata,"@",kline2,"@"
              if(j == 1 .and. nol /= 0) then
                 infiles(1)%nol = nol
!                 write(6,*) nol,kdata
                 if(nol /= 0) read(kdata,'(E15.7)') infiles(1)%r(nol)
!                 write(6,*) "rho",nol,infiles(1)%r(nol)
                 kline = kline2
              else
                 if(index(kline2,"E") /= 0) then
                    kdata = trim(kline2)
                    if(nol /= 0) then
                       read(kdata,'(E15.7)') infiles(1)%data(nol)
!                       write(6,*) "data",nol,infiles(1)%data(nol)
                    else
                       read(kdata,'(E15.7)') infiles(1)%totS
!                       write(6,*) "total",nol,infiles(1)%totS
                    end if
                    kline = kline2
                 else
                    exit
                 end if
              end if
!              pause
              cycle
           end if
        end do
     end do

!!$     write(6,*) infiles(1)%nol
!!$     do nol = 1, infiles(1)%nol
!!$        write(6,*) nol,infiles(1)%r(nol),infiles(1)%data(nol)
!!$     end do
!!$     write(6,*) infiles(1)%totS

     close(io)

  ! *** Arbitrary input *********************************************
  else if(kfile == 3) then

     ! Read data
     do i = 1, n_infiles
        open(newunit=io,file=TXFNAM,iostat=IST,status='old',form='formatted')
        if (IST == 0) then
           write(6,*) '# ASCII FILE ( ', TXFNAM(1:len_trim(TXFNAM)),  &
                &     ' ) IS ASSIGNED FOR INPUT.'
        else if (IST > 0) then
           write(6,*) 'XX  ASCII FILE open ERROR !, iostat = ', IST
        else if (IST < 0) then
           deallocate(infiles) ; iflag_file = 0 ; infile_name = ' '
           return
        end if

        ! Select a type of input data from ascii files
        do 
           write(6,'(A,I2,A)') '# Select a type of input data from ascii file for ',i,' data.'
           write(6,*) '# 1: PNBP, 2: PNBT1, 3: PNBT2, 4: PRF, 5: LQe4'
           flush(6)
           read(*,'(I2)',iostat=ist) itype
           if(ist > 0) then
              cycle
           else if(ist < 0) then
              deallocate(infiles) ; iflag_file = 0 ; infile_name = ' '
              return
           else
              if(itype <= 0 .or. itype > 5) cycle
              exit
           end if
        end do
        infiles(i)%name = datatype(itype)

        ! Which column indicates the mesh data
        do
           write(6,*) '# Which column indicates the mesh data in ', &
                &     TXFNAM(1:len_trim(TXFNAM)), ' ?'
           read(*,'(I2)', iostat=ist) infiles(i)%ncol_mesh
           if(infiles(i)%ncol_mesh == 0) cycle
           if (ist > 0) then
              cycle
           else if(ist < 0) then
              deallocate(infiles) ; iflag_file = 0 ; infile_name = ' '
              return
           else
              exit
           end if
        end do

        ! Which column indicates the data which you would like to use
        do
           write(6,*) '# Which column indicates the data which you would like to use in ', &
                &     TXFNAM(1:len_trim(TXFNAM)), ' ?'
           read(*,'(I2)', iostat=ist) infiles(i)%ncol_data
           if(infiles(i)%ncol_data == 0) cycle
           if (ist > 0) then
              cycle
           else if(ist < 0) then
              deallocate(infiles) ; iflag_file = 0 ; infile_name = ' '
              return
           else
              exit
           end if
        end do

        ! Read data from the file
        nol = 0 ! number of lines
        do 
           read(io,'(A)',iostat=IST) kline
           if(ist < 0) exit ! detect the end of the file
           nol = nol + 1

           j = 0
           do
              j = j + 1
              kline = trim(adjustl(kline))
              call ksplit_tx(kline,' ',kline1,kline2)
              if(j == infiles(i)%ncol_mesh) then
                 kmesh = trim(kline1)
              else if(j == infiles(i)%ncol_data) then
                 kdata = trim(kline1)
              end if
              kline = kline2
              if(j >= infiles(i)%ncol_mesh .and. j >= infiles(i)%ncol_data) exit
           end do

           read(kmesh,'(E15.7)') infiles(i)%r(nol)
           read(kdata,'(E15.7)') infiles(i)%data(nol)

!           write(6,*) infiles(i)%r(nol),infiles(i)%data(nol)
        end do
        if(nol > nmax_file) then
           write(6,*) 'XX The number of lines in the file exceeds ',nmax_file
           write(6,*) '   Program is terminated.'
           stop
        else
           infiles(i)%nol = nol
        end if

        close(io)

     end do

  end if
  iflag_file = kfile ; infile_name = TXFNAM

end subroutine ascii_input

!***************************************************************
!
!   Return the number of the datatype in "infiles"
!
!***************************************************************

integer(4) function detect_datatype(kchar)

  use tx_commons, only : infiles, n_infiles
  use libchar, only : kmatch
  implicit none
  character(len=*), intent(in) :: kchar
  integer(4) :: i

  do i = 1, n_infiles
     if(kmatch(infiles(i)%name,kchar)) then
        detect_datatype = i
        return
     end if
  end do

  detect_datatype = 0

end function detect_datatype

!***************************************************************
!
!   Control routine for writing ASCII data file
!
!***************************************************************

subroutine outfile
  use tx_commons, only : stdin_ended
  implicit none
  integer(4) :: n, ist

  do
     write(6,'(1X,A,1X,I1)') '## outfile: input ?'
     write(6,'(1X,A)') '##          1: profiles, 2: 0d quantities, 3: 0d variables, ' &
          & // '4: for OFMC, 5: for TOPICS/NTMAIN, 6: for CHARROT, 7: particle sources, 9: exit'
     if (stdin_ended) return
     read(*,'(I1)',iostat=ist) n
     if (ist > 0) then
        cycle
     else if (ist < 0) then
        stdin_ended = .true.
        return
     end if

     select case(n)
     case(1)
        call write_profs
     case(2)
        call write_quantities
     case(3)
        call write_vars
     case(4)
        call for_ofmc
     case(5)
        call for_ntmain
     case(6)
        call for_charrot
     case(7)
        call write_sources
     case default
        exit
     end select
  end do

end subroutine outfile

!***************************************************************
!
!   Append the particle sources through the SOL and the wall
!
!***************************************************************

subroutine write_sources
  use tx_glob, only : particle_source_balance
  implicit none
  integer(4) :: iu

  open(newunit=iu, file='tx_sources.dat', position='append', action='write')
  call particle_source_balance(iu)
  close(iu)
  write(6,*) '# PARTICLE SOURCES APPENDED TO tx_sources.dat.'

end subroutine write_sources

!***************************************************************
!
!   Write ASCII data file in terms of profiles
!
!***************************************************************

subroutine write_profs
  use tx_commons, only : NRMAX
  use tx_graphic, only : GT, GX, GY, NGR, NGYRM
  use libfio, only : FWopen
  implicit none
  integer(4) :: NG, NR, NGYR, IERR
  character(len=80) :: fmt = ''

  call FWopen(21,'tx_profs.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX write_profs: FWopen: IERR=', IERR
     return
  end if

  write(fmt, *) '(',NGYRM+2,'ES15.7)'

  do NG = 0, NGR
     do NR = 0, NRMAX
        write(21,fmt)  GT(NG), GX(NR), (GY%v(NR,NG,NGYR), NGYR=1,NGYRM)
     end do
     write(21,*) 
  end do

  close(21)

end subroutine write_profs

!***************************************************************
!
!   Write ASCII data file in terms of 0d quantities
!
!***************************************************************

subroutine write_quantities
  use tx_graphic, only : GTX, GTY, NGT, NGYTM
  use libfio, only : FWopen
  implicit none
  integer(4) :: NG, NGYT, IERR
  character(len=80) :: fmt = ''

  call FWopen(21,'tx_quantities.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX write_quantities: FWopen: IERR=', IERR
     return
  end if

  write(fmt, *) '(',NGYTM+1,'ES15.7)'

  do NG = 0, NGT
     write(21,fmt)  GTX(NG), (GTY(NG,NGYT), NGYT=1,NGYTM)
  end do

  close(21)

end subroutine write_quantities

!***************************************************************
!
!   Write ASCII data file in terms of 0d variables
!
!***************************************************************

subroutine write_vars
  use tx_graphic, only : GVX, GVY, NGVV, NGYVM
  use libfio, only : FWopen
  implicit none
  integer(4) :: NG, NGYV, IERR
  character(len=80) :: fmt = ''

  call FWopen(21,'tx_vars.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX write_vars: FWopen: IERR=', IERR
     return
  end if

  write(fmt, *) '(',NGYVM+1,'ES15.7)'

  do NG = 0, NGVV
     write(21,fmt)  GVX(NG), (GVY(NG,NGYV), NGYV=1,NGYVM)
  end do

  close(21)

end subroutine write_vars

!***************************************************************
!
!   Write ASCII data file for TOPICS/NTMAIN
!
!***************************************************************

subroutine for_ntmain
  use tx_commons, only : NRMAX, Rho, PN01V, PN02V, Var
  use libfio, only : FWopen
  implicit none
  integer(4) :: IERR, NR

  call FWopen(21,'tx_ntmain.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX for_ntmain: FWopen: IERR=', IERR
     return
  end if

  write(21,'(4X,I3)') NRMAX+1
  do nr = 0, nrmax
     write(21,'(7ES17.9)') Rho(NR), Var(NR,1)%n*1.D20, Var(NR,2)%n*1.D20, &
          & Var(NR,1)%T*1.D3, Var(NR,2)%T*1.D3, PN01V(NR)*1.D20, PN02V(NR)*1.D20
  end do

  close(21)

end subroutine for_ntmain

!***************************************************************
!
!   Write ASCII data file for OFMC
!
!      Prepare three profiles of electron density and
!         electron and ion temperatures from TASK/TX
!
!      For OFMC code, the ele. density should be in m^{-3} and
!                     the temperatures in eV.
!
!      This routine was originally built as task/data_spline.
!
!***************************************************************

subroutine for_ofmc
  use tx_commons, only : NRMAX, NRA, RR, RA, BB, rIp, psiV, Var, array_init_NR
  use libfio, only : FWopen
  use libspl1d, only : spl1d, spl1df
  implicit none

  integer(4) :: nmax, ist, i, j, NR, ierr
  real(8) :: dpsi, psirho_a
  real(8), dimension(:), allocatable :: psirho, psi_out, data, data_out, deriv
  real(8), dimension(:,:), allocatable :: u

  do
     write(6,'(1X,A)') '## for_ofmc: number of grids ? [Default=30]'
     read(*,'(I3)',iostat=ist) nmax
     if (ist > 0) then
        cycle
     else if (ist < 0) then
        return
     end if
     if(nmax == 0) then
        nmax = 30
     else if(nmax < 0) then
        cycle
     end if
     exit
  end do

  ! Output file open
  call FWopen(21,'tx_ofmc.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX for_ofmc: FWopen: IERR=', IERR
     return
  end if

  write(21,'(A)') 'PARM'
  write(21,'(1X,A5,1X,I5,2(1X,A5))') '*****',-1,('*****', j = 1, 2)
  write(21,'(1X,I5,11(1X,A5))') nmax,('*****', j = 1, 11)
  write(21,'(4(ES12.5))') RR,RA,BB,rIp
  write(21,'(3(1X,A11),2(ES12.5))') ('************', j = 1, 3), 1.d0, 0.d0
  write(21,'(5(1X,A11))') ('************', j = 1, 5)

  allocate(psirho, data, deriv, source=array_init_NR)
  allocate(u(1:4,0:NRMAX))
  allocate(psi_out(1:nmax),data_out(1:nmax+1))

  ! normalized psi with respect to rho
  do NR = 1, NRMAX
     psirho(NR) = psiV(NR) - psiV(0) !- RR * (AphV(NR) - AphV(0))
  end do
  psirho_a = psirho(NRA)
  psirho(:) = psirho(:) / psirho_a

  ! Equally-spaced psi on a half mesh
  dpsi = 1.d0 / real(nmax,8)
  do i = 1, nmax
     psi_out(i) = i * dpsi - 0.5d0 * dpsi
  end do

  do i = 1, 3
     if     (i == 1) then ! Ne
        data(:) = Var(:,1)%n * 1.D20 ! in m^{-3}
     else if(i == 2) then ! Te
        data(:) = Var(:,1)%T * 1.D3  ! in eV
     else if(i == 3) then ! Ti
        data(:) = Var(:,2)%T * 1.D3  ! in eV
     end if

     ! Output data on the magnetic axis
     write(21,'(2(ES12.5))') data(0), -1.d0

     ! Spline data for equally-spaced psi
     call spl1d(psirho,data,deriv,u,NRMAX+1,0,ierr)
     if(ierr /= 0) stop 'Error at spl1d'
     do j = 1, nmax
        call spl1df(psi_out(j),data_out(j),psirho,u,NRMAX+1,ierr)
        if(ierr /= 0) stop 'Error at spl1df'
     end do

     ! Add data on the separatrix to output data array
     data_out(nmax+1) = data(NRA)

     ! Output data
     write(21,'(4(ES18.10))') (data_out(j), j = 1, nmax+1)
  end do

  deallocate(psirho,data,deriv,u,psi_out,data_out)

  close(21)

end subroutine for_ofmc

!***************************************************************
!
!   Write ASCII data file for CHARROT
!
!     Output profile file for CHARROT code. 
!     The file format is consistent with output from JT-60U analysis tools.
!
!***************************************************************

subroutine for_charrot
  use tx_commons, only : NRA
  use tx_graphic, only : GX, GY, NGR
  use libfio, only : FWopen
  implicit none
  integer(4) :: NR, IERR
  character(len=80) :: fmt = ''

  call FWopen(21,'tx_charrot.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX for_charrot: FWopen: IERR=', IERR
     return
  end if

  write(21,'(A)') 'RO/A,Ti[keV],Vt[km/s],Te[keV],ne[1e+19m^-3],Zeff,BEpara[TV/m]'

  write(fmt, *) "(ES12.6,6(', ',ES13.6))"

  do NR = 0, NRA
     write(21,fmt)  GX(NR) &
          & , GY%v(NR,NGR, 66)        & ! Ti [keV]
          & , GY%v(NR,NGR, 89)*1.e-3  & ! Vt_Z [km/s]
          & , GY%v(NR,NGR, 53)        & ! Te [keV]
          & , GY%v(NR,NGR, 52)*1.e-19 & ! ne [10^{19}m^{-3}]
          & , GY%v(NR,NGR,202)        & ! Zeff
          & , GY%v(NR,NGR, 96)          ! BEpara [TV/m]
  end do

  close(21)

end subroutine for_charrot

!***************************************************************
!
!   Read ASCII data file for initial profiles
!
!      Prepare profiles of 
!         PNeV          : electron density [/m^3], 
!         PTeV, PTiV    : electron and ion temperatures [eV], and
!         AEE*PNeV*UephV: toroidal current density [A/m^2]
!         for TASK/TX
!
!      List structure has been used.
!
!***************************************************************

subroutine initprof_input(nr, idx, out)

  use tx_commons, only : Rho, AEE
  use libitp, only : aitken2p, fctr4pt
  use libfio, only : fropen
  use libspl1d, only : spl1d, spl1df

  integer(4), optional, intent(in) :: nr, idx
  real(8), optional, intent(out) :: out
  integer(4) :: nintin, ier, ist, k, iflag, j
  integer(4), save :: nrinmax
  real(8), dimension(:), allocatable, save :: rho_in, deriv
  real(8), dimension(:,:), allocatable, save :: prof_in, u1, u2, u3, u4

  type unit
     integer(4)          :: l_p
     real(8)             :: rho_p, prof1_p, prof2_p, prof3_p, prof4_p
     type(unit), pointer :: next
  end type unit
  type(unit), pointer :: ent, new, p
  integer(4) :: l
  real(8)    :: rho_r, prof1_r, prof2_r, prof3_r, prof4_r

  if((present(nr) .eqv. .false.) .and. (present(idx) .eqv. .false.)) then
     allocate(ent); nullify(ent%next)

     nintin = 24
     call FRopen(nintin,'initprof.dat',1,0,'INIT PROF',ier)
     if(ier /= 0) return

     do
!!!        read(nintin,'(1X,I3,5(E9.3))',iostat=ist) l, rho_r, prof1_r, prof2_r, prof3_r, prof4_r
        read(nintin,'(1X,I3,5(E10.3))',iostat=ist) l, rho_r, prof1_r, prof2_r, prof3_r, prof4_r
        if(ist > 0) then
           write(6,*) 'XX file format is invalid. Aborting.'
           close(nintin)
           stop
        else if (ist < 0) then
           exit
        end if
        call rearrange
     end do

     iflag = 0
     p => ent%next
     k = 0
     do
        if(associated(p)) then 
!           write(6,*) p%l_p, p%rho_p, p%prof1_p, p%prof2_p, p%prof3_p, p%prof4_p
           if(k == 0 .and. p%rho_p /= 1.d0) iflag = iflag + 1
           if(p%l_p == 1 .and. p%rho_p /= 0.d0) iflag = iflag + 2
           k = k + 1
           p => p%next
        else
           exit
        end if
     end do
     if(iflag == 0) then
        nrinmax = k
     else if(iflag == 1 .or. iflag == 2) then
        nrinmax = k + 1
     else
        nrinmax = k + 2
     end if

     allocate(rho_in(1:nrinmax),prof_in(1:nrinmax,1:4))
     p => ent%next
     do
        if(associated(p)) then
           rho_in (p%l_p)   = p%rho_p
           prof_in(p%l_p,1) = p%prof1_p
           prof_in(p%l_p,2) = p%prof2_p
           prof_in(p%l_p,3) = p%prof3_p
           prof_in(p%l_p,4) = p%prof4_p
           p => p%next
        else
           exit
        end if
     end do

     if(iflag == 1) then ! Extrapolate values at r/a = 1.
        rho_in(k+1) = 1.d0
        do j = 1, 4
           prof_in(nrinmax,j) = AITKEN2P(rho_in(nrinmax),prof_in(nrinmax-1,j), &
                & prof_in(nrinmax-2,j),prof_in(nrinmax-3,j),rho_in(nrinmax-1), &
                & rho_in(nrinmax-2),rho_in(nrinmax-3))
        end do
     else if(iflag == 2) then ! Extrapolate values at r/a = 0.
        rho_in(2:k+1) = rho_in(1:k)
        rho_in(1) = 0.d0
        prof_in(2:k+1,1:4) = prof_in(1:k,1:4)
        do j = 1, 4
           prof_in(1,j) = FCTR4pt(rho_in(2),rho_in(3),rho_in(4), &
                &                 prof_in(2,j),prof_in(3,j),prof_in(4,j))
        end do
     else if(iflag == 3) then ! Extrapolate values at r/a = 0 and 1.
        rho_in(2:k+1) = rho_in(1:k)
        rho_in(1) = 0.d0
        rho_in(k+2) = 1.d0
        prof_in(2:k+1,1:4) = prof_in(1:k,1:4)
        do j = 1, 4
           prof_in(1,j) = FCTR4pt(rho_in(2),rho_in(3),rho_in(4), &
                &                 prof_in(2,j),prof_in(3,j),prof_in(4,j))
           prof_in(nrinmax,j) = AITKEN2P(rho_in(nrinmax),prof_in(nrinmax-1,j), &
                & prof_in(nrinmax-2,j),prof_in(nrinmax-3,j),rho_in(nrinmax-1), &
                & rho_in(nrinmax-2),rho_in(nrinmax-3))
        end do
     end if

!!$     do k = 1, nrinmax
!!$        write(6,'(1X,I3,5ES10.3)') k, rho_in(k),(prof_in(k,j),j=1,4)
!!$     end do
!!$
     close(nintin)
     if(associated(ent)) deallocate(ent)
     if(associated(new)) deallocate(new)

     allocate(deriv(1:nrinmax),u1(1:4,1:nrinmax),u2(1:4,1:nrinmax),u3(1:4,1:nrinmax),u4(1:4,1:nrinmax))
     call spl1d(rho_in,prof_in(1:nrinmax,1),deriv,u1,nrinmax,0,ier)
     if(ier /= 0) stop 'Error at spl1d in initprof_input: prof_in(1)'
     call spl1d(rho_in,prof_in(1:nrinmax,2),deriv,u2,nrinmax,0,ier)
     if(ier /= 0) stop 'Error at spl1d in initprof_input: prof_in(2)'
     call spl1d(rho_in,prof_in(1:nrinmax,3),deriv,u3,nrinmax,0,ier)
     if(ier /= 0) stop 'Error at spl1d in initprof_input: prof_in(3)'
     call spl1d(rho_in,prof_in(1:nrinmax,4),deriv,u4,nrinmax,0,ier)
     if(ier /= 0) stop 'Error at spl1d in initprof_input: prof_in(4)'

  else if((present(nr) .eqv. .true.) .and. (present(idx) .eqv. .true.)) then
     if(idx == 1) then ! Electron density
        call spl1df(Rho(nr),out,rho_in,u1,nrinmax,ier)
        if(ier /= 0) stop 'Error at spl1df in initprof_input: idx = 1'
        out = out * 1.d-20
     else if(idx == 2) then ! Electron temperature
        call spl1df(Rho(nr),out,rho_in,u2,nrinmax,ier)
        if(ier /= 0) stop 'Error at spl1df in initprof_input: idx = 2'
        out = out * 1.d-3
     else if(idx == 3) then ! Ion temperature
        call spl1df(Rho(nr),out,rho_in,u3,nrinmax,ier)
        if(ier /= 0) stop 'Error at spl1df in initprof_input: idx = 3'
        out = out * 1.d-3
     else if(idx == 4) then ! Toroidal current
        call spl1df(Rho(nr),out,rho_in,u4,nrinmax,ier)
        if(ier /= 0) stop 'Error at spl1df in initprof_input: idx = 4'
        out = out
     else 
        stop 'initprof_input: wrong input!'
     end if

  else if((present(nr) .eqv. .false.) .and. (present(idx) .eqv. .true.)) then
     if(idx == 0) then
        deallocate(rho_in,prof_in,deriv,u1,u2,u3,u4)
     else
        stop 'initprof_input: wrong input!'
     end if
  else
     stop 'initprof_input: wrong input!'
  end if

contains
  subroutine rearrange
    allocate(new)
    new = unit(l, rho_r, prof1_r, prof2_r, prof3_r, prof4_r, ent%next)
    ent%next => new
  end subroutine rearrange

end subroutine initprof_input
