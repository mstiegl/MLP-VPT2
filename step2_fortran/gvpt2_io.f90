!=====================================================================
! gvpt2_io
!
! File readers for the Fortran translation of gvpt2.py:
!   read_xyz              -> geometry file  {mol}_opt.out   (xyz, Angstrom)
!   read_freqs_and_vecs   -> {mol}_nma.out  (harmonic freqs + normal modes)
!   parse_cubic_tensor    -> {mol}_cubics.out
!   parse_quartic_tensor  -> {mol}_quartics.out
!
! The parsers follow the line-oriented logic of the Python originals
! (regular expressions replaced by explicit tokenisation) so that the
! same files produced by step1_python or step1_fortran are accepted.
!=====================================================================
module gvpt2_io
  use gvpt2_constants, only: dp, to_upper
  use iso_fortran_env, only: error_unit
  implicit none
  private

  public :: read_xyz, read_freqs_and_vecs, parse_cubic_tensor, parse_quartic_tensor
  public :: split_tokens, is_int_token, is_real_token, read_lines
  public :: pyf, fmt_signed, fix_f0, py_repr, np_array_str

  integer, parameter, public :: LINE_LEN = 4096
  integer, parameter, public :: TOK_LEN  = 64
  integer, parameter :: MAX_TOK = 64

contains

  !-------------------------------------------------------------------
  ! Read an entire text file into an array of lines.
  !-------------------------------------------------------------------
  subroutine read_lines(path, lines, nlines)
    character(len=*), intent(in) :: path
    character(len=LINE_LEN), allocatable, intent(out) :: lines(:)
    integer, intent(out) :: nlines
    integer :: u, ios, i
    character(len=LINE_LEN) :: buf

    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,'(A)') 'ERROR: cannot open file "' // trim(path) // '"'
      stop 1
    end if

    nlines = 0
    do
      read(u, '(A)', iostat=ios) buf
      if (ios /= 0) exit
      nlines = nlines + 1
    end do

    allocate(lines(max(nlines, 1)))
    rewind(u)
    do i = 1, nlines
      read(u, '(A)') lines(i)
    end do
    close(u)
  end subroutine read_lines

  !-------------------------------------------------------------------
  ! Split a line on blanks/tabs into tokens (like str.split()).
  !-------------------------------------------------------------------
  subroutine split_tokens(line, tokens, ntok)
    character(len=*), intent(in) :: line
    character(len=TOK_LEN), intent(out) :: tokens(MAX_TOK)
    integer, intent(out) :: ntok
    character(len=len(line)) :: s
    integer :: i, n, start
    logical :: intok

    s = line
    ! treat tabs and carriage returns as blanks
    do i = 1, len(s)
      if (s(i:i) == achar(9) .or. s(i:i) == achar(13)) s(i:i) = ' '
    end do

    n = len_trim(s)
    ntok = 0
    tokens = ''
    intok = .false.
    start = 1
    do i = 1, n + 1
      if (i <= n .and. s(i:i) /= ' ') then
        if (.not. intok) then
          intok = .true.
          start = i
        end if
      else
        if (intok) then
          ntok = ntok + 1
          if (ntok > MAX_TOK) then
            write(error_unit,'(A)') 'ERROR: too many tokens on a line'
            stop 1
          end if
          tokens(ntok) = s(start:i-1)
          intok = .false.
        end if
      end if
    end do
  end subroutine split_tokens

  !-------------------------------------------------------------------
  ! True if the token is an unsigned decimal integer (regex \d+).
  !-------------------------------------------------------------------
  pure function is_int_token(tok) result(ok)
    character(len=*), intent(in) :: tok
    logical :: ok
    integer :: i, n
    n = len_trim(tok)
    ok = n > 0
    do i = 1, n
      if (tok(i:i) < '0' .or. tok(i:i) > '9') then
        ok = .false.
        return
      end if
    end do
  end function is_int_token

  !-------------------------------------------------------------------
  ! True if the token looks like a floating-point number:
  !   [-+]? (digits [. digits?] | . digits) ([eEdD] [-+]? digits)?
  ! This accepts everything the Python regexes accept (plus Fortran
  ! style D exponents and plain integers).
  !-------------------------------------------------------------------
  pure function is_real_token(tok) result(ok)
    character(len=*), intent(in) :: tok
    logical :: ok
    integer :: n, pos, nd1, nd2, ne

    n = len_trim(tok)
    ok = .false.
    if (n == 0) return
    pos = 1
    if (tok(pos:pos) == '+' .or. tok(pos:pos) == '-') pos = pos + 1

    nd1 = 0
    do while (pos <= n)
      if (tok(pos:pos) < '0' .or. tok(pos:pos) > '9') exit
      nd1 = nd1 + 1
      pos = pos + 1
    end do

    nd2 = 0
    if (pos <= n) then
      if (tok(pos:pos) == '.') then
        pos = pos + 1
        do while (pos <= n)
          if (tok(pos:pos) < '0' .or. tok(pos:pos) > '9') exit
          nd2 = nd2 + 1
          pos = pos + 1
        end do
      end if
    end if
    if (nd1 + nd2 == 0) return

    if (pos <= n) then
      if (index('eEdD', tok(pos:pos)) > 0) then
        pos = pos + 1
        if (pos <= n) then
          if (tok(pos:pos) == '+' .or. tok(pos:pos) == '-') pos = pos + 1
        end if
        ne = 0
        do while (pos <= n)
          if (tok(pos:pos) < '0' .or. tok(pos:pos) > '9') exit
          ne = ne + 1
          pos = pos + 1
        end do
        if (ne == 0) return
      end if
    end if

    ok = (pos == n + 1)
  end function is_real_token

  !-------------------------------------------------------------------
  ! Convert a validated real token to a double (D and E exponents).
  !-------------------------------------------------------------------
  function token_to_real(tok) result(val)
    character(len=*), intent(in) :: tok
    real(dp) :: val
    integer :: ios
    read(tok, *, iostat=ios) val
    if (ios /= 0) then
      write(error_unit,'(A)') 'ERROR: cannot convert "' // trim(tok) // '" to a real number'
      stop 1
    end if
  end function token_to_real

  !-------------------------------------------------------------------
  ! Read an xyz geometry file (natoms / comment / symbol x y z ...).
  ! Positions are returned as given in the file (Angstrom).
  !-------------------------------------------------------------------
  subroutine read_xyz(path, natoms, symbols, pos)
    character(len=*), intent(in) :: path
    integer, intent(out) :: natoms
    character(len=8), allocatable, intent(out) :: symbols(:)
    real(dp), allocatable, intent(out) :: pos(:,:)     ! (natoms, 3)
    character(len=LINE_LEN), allocatable :: lines(:)
    character(len=TOK_LEN) :: tok(MAX_TOK)
    integer :: nlines, ntok, i, ios, k

    call read_lines(path, lines, nlines)
    if (nlines < 1) then
      write(error_unit,'(A)') 'ERROR: empty geometry file "' // trim(path) // '"'
      stop 1
    end if

    read(lines(1), *, iostat=ios) natoms
    if (ios /= 0 .or. natoms < 1) then
      write(error_unit,'(A)') 'ERROR: first line of "' // trim(path) // '" must be the number of atoms'
      stop 1
    end if
    if (nlines < natoms + 2) then
      write(error_unit,'(A)') 'ERROR: geometry file "' // trim(path) // '" has too few lines'
      stop 1
    end if

    allocate(symbols(natoms), pos(natoms, 3))
    do i = 1, natoms
      call split_tokens(lines(i + 2), tok, ntok)
      if (ntok < 4) then
        write(error_unit,'(A,I0,A)') 'ERROR: atom line ', i, ' in "' // trim(path) // '" needs symbol x y z'
        stop 1
      end if
      symbols(i) = tok(1)(1:8)
      do k = 1, 3
        if (.not. is_real_token(tok(k + 1))) then
          write(error_unit,'(A,I0,A)') 'ERROR: bad coordinate on atom line ', i, ' in "' // trim(path) // '"'
          stop 1
        end if
        pos(i, k) = token_to_real(tok(k + 1))
      end do
    end do
  end subroutine read_xyz

  !-------------------------------------------------------------------
  ! True if the line is a "Frequencies (cm^-1):" header
  ! (regex ^\s*Frequencies\s*\(cm\^-1\)\s*:\s*$).
  !-------------------------------------------------------------------
  pure function is_freq_header(line) result(ok)
    character(len=*), intent(in) :: line
    logical :: ok
    character(len=len(line)) :: s
    integer :: i, n
    ! remove all blanks/tabs and compare
    s = ''
    n = 0
    do i = 1, len_trim(line)
      if (line(i:i) /= ' ' .and. line(i:i) /= achar(9) .and. line(i:i) /= achar(13)) then
        n = n + 1
        s(n:n) = line(i:i)
      end if
    end do
    ok = (s(1:max(n,1)) == 'Frequencies(cm^-1):' .and. n == len('Frequencies(cm^-1):'))
  end function is_freq_header

  !-------------------------------------------------------------------
  ! If the line is a "Mode N" header (regex ^\s*mode\s+(\d+)\s*$,
  ! case-insensitive) return N, otherwise 0.
  !-------------------------------------------------------------------
  function mode_header_index(line) result(idx)
    character(len=*), intent(in) :: line
    integer :: idx
    character(len=TOK_LEN) :: tok(MAX_TOK)
    integer :: ntok, ios
    idx = 0
    call split_tokens(line, tok, ntok)
    if (ntok /= 2) return
    if (to_upper(trim(tok(1))) /= 'MODE') return
    if (.not. is_int_token(tok(2))) return
    read(tok(2), *, iostat=ios) idx
    if (ios /= 0) idx = 0
  end function mode_header_index

  !-------------------------------------------------------------------
  ! Parse the normal-mode file. Harmonic frequencies are taken from the
  ! LAST "Frequencies (cm^-1):" block; the mass-weighted normal-mode
  ! vectors from the "Mode N" blocks. As in gvpt2.py, the last n_modes
  ! entries are returned (external modes, if present, are dropped).
  !
  !   omega_cm(n_modes)          harmonic frequencies, cm^-1
  !   modes(n_modes, natom, 3)   mass-weighted normal-mode vectors
  !-------------------------------------------------------------------
  subroutine read_freqs_and_vecs(path, n_modes, omega_cm, modes, natom)
    character(len=*), intent(in) :: path
    integer, intent(in) :: n_modes
    real(dp), allocatable, intent(out) :: omega_cm(:)
    real(dp), allocatable, intent(out) :: modes(:,:,:)
    integer, intent(out) :: natom

    character(len=LINE_LEN), allocatable :: lines(:)
    character(len=TOK_LEN) :: tok(MAX_TOK)
    integer :: nlines, ntok, i, j, hdr, nfreq, idx, nhdr
    integer :: first, last, nrow, m, k
    real(dp), allocatable :: freqs(:), blocks(:,:,:)
    integer, allocatable :: hdr_line(:), hdr_mode(:)
    logical, allocatable :: have(:)

    call read_lines(path, lines, nlines)

    ! ---- locate the last "Frequencies (cm^-1):" header ----
    hdr = 0
    do i = 1, nlines
      if (is_freq_header(lines(i))) hdr = i
    end do
    if (hdr == 0) then
      write(error_unit,'(A)') 'ERROR: Could not find "Frequencies (cm^-1):" block in ' // trim(path)
      stop 1
    end if

    ! ---- parse "index  value" lines that follow ----
    allocate(freqs(max(nlines - hdr, 1)))
    nfreq = 0
    do i = hdr + 1, nlines
      if (len_trim(lines(i)) == 0) then
        if (nfreq > 0) exit
        cycle
      end if
      call split_tokens(lines(i), tok, ntok)
      if (ntok == 2) then
        if (is_int_token(tok(1)) .and. is_real_token(tok(2))) then
          nfreq = nfreq + 1
          freqs(nfreq) = token_to_real(tok(2))
          cycle
        end if
      end if
      if (nfreq > 0) exit
    end do

    if (nfreq < n_modes) then
      write(error_unit,'(A,I0,A,I0,A)') 'ERROR: Found only ', nfreq, ' freqs in last block, need ', n_modes, '.'
      stop 1
    end if

    ! ---- locate all "Mode N" headers ----
    allocate(hdr_line(nlines), hdr_mode(nlines))
    nhdr = 0
    do i = 1, nlines
      idx = mode_header_index(lines(i))
      if (idx > 0) then
        nhdr = nhdr + 1
        hdr_line(nhdr) = i
        hdr_mode(nhdr) = idx
      end if
    end do
    if (nhdr == 0) then
      write(error_unit,'(A)') 'ERROR: no "Mode N" blocks found in ' // trim(path)
      stop 1
    end if

    ! ---- number of atoms from the first block ----
    natom = 0
    first = hdr_line(1) + 1
    if (nhdr > 1) then
      last = hdr_line(2) - 1
    else
      last = nlines
    end if
    do i = first, last
      if (len_trim(lines(i)) > 0) natom = natom + 1
    end do
    if (natom == 0) then
      write(error_unit,'(A)') 'ERROR: empty normal-mode block in ' // trim(path)
      stop 1
    end if

    allocate(blocks(nfreq, natom, 3), have(nfreq))
    blocks = 0.0_dp
    have = .false.

    ! ---- parse each block ----
    do m = 1, nhdr
      idx = hdr_mode(m)
      if (idx > nfreq) then
        write(error_unit,'(A,I0,A,I0,A)') 'ERROR: "Mode ', idx, '" header exceeds the number of frequencies (', nfreq, ')'
        stop 1
      end if
      first = hdr_line(m) + 1
      if (m < nhdr) then
        last = hdr_line(m + 1) - 1
      else
        last = nlines
      end if
      nrow = 0
      do i = first, last
        if (len_trim(lines(i)) == 0) cycle
        call split_tokens(lines(i), tok, ntok)
        if (ntok /= 3) then
          write(error_unit,'(A,I0,A,I0)') 'ERROR: Expected 3 numbers in mode ', idx, ', got ', ntok
          stop 1
        end if
        nrow = nrow + 1
        if (nrow > natom) then
          write(error_unit,'(A,I0,A)') 'ERROR: mode ', idx, ' has more rows than the first mode block'
          stop 1
        end if
        do k = 1, 3
          if (.not. is_real_token(tok(k))) then
            write(error_unit,'(A,I0)') 'ERROR: non-numeric entry in mode ', idx
            stop 1
          end if
          blocks(idx, nrow, k) = token_to_real(tok(k))
        end do
      end do
      if (nrow /= natom) then
        write(error_unit,'(A,I0,A)') 'ERROR: mode ', idx, ' has a different number of rows than the first mode block'
        stop 1
      end if
      have(idx) = .true.
    end do

    ! ---- keep the last n_modes entries (as freqs[-n_modes:]) ----
    allocate(omega_cm(n_modes), modes(n_modes, natom, 3))
    do j = 1, n_modes
      idx = nfreq - n_modes + j
      omega_cm(j) = freqs(idx)
      if (.not. have(idx)) then
        write(error_unit,'(A,I0,A)') 'ERROR: normal-mode vector for mode ', idx, ' not found in ' // trim(path)
        stop 1
      end if
      modes(j, :, :) = blocks(idx, :, :)
    end do

  end subroutine read_freqs_and_vecs

  !-------------------------------------------------------------------
  ! Cubic force constants in mass-weighted normal coordinates Q.
  ! Lines of the form "i j k value" are accepted (anything else, e.g.
  ! "#" comments, is skipped); the value is stored for all permutations
  ! of (i,j,k). Entries not present in the file are zero.
  !-------------------------------------------------------------------
  subroutine parse_cubic_tensor(path, n, phi3)
    character(len=*), intent(in) :: path
    integer, intent(in) :: n
    real(dp), intent(out) :: phi3(n, n, n)
    character(len=LINE_LEN), allocatable :: lines(:)
    character(len=TOK_LEN) :: tok(MAX_TOK)
    integer :: nlines, ntok, l, i, j, k
    real(dp) :: val

    phi3 = 0.0_dp
    call read_lines(path, lines, nlines)
    do l = 1, nlines
      call split_tokens(lines(l), tok, ntok)
      if (ntok /= 4) cycle
      if (.not. (is_int_token(tok(1)) .and. is_int_token(tok(2)) .and. &
                 is_int_token(tok(3)) .and. is_real_token(tok(4)))) cycle
      read(tok(1), *) i
      read(tok(2), *) j
      read(tok(3), *) k
      val = token_to_real(tok(4))
      ! indices outside 1..n are never used by gvpt2.py; ignore them
      if (min(i, j, k) < 1 .or. max(i, j, k) > n) cycle
      phi3(i, j, k) = val
      phi3(i, k, j) = val
      phi3(j, i, k) = val
      phi3(j, k, i) = val
      phi3(k, i, j) = val
      phi3(k, j, i) = val
    end do
  end subroutine parse_cubic_tensor

  !-------------------------------------------------------------------
  ! Semi-diagonal quartic force constants Phi_iikk. Lines of the form
  ! "i j k l value" are accepted; as in gvpt2.py only (i,k) are used
  ! (the file is assumed to contain i i k k lines) and the value is
  ! stored symmetrically in phi4(i,k) = phi4(k,i).
  !-------------------------------------------------------------------
  subroutine parse_quartic_tensor(path, n, phi4)
    character(len=*), intent(in) :: path
    integer, intent(in) :: n
    real(dp), intent(out) :: phi4(n, n)
    character(len=LINE_LEN), allocatable :: lines(:)
    character(len=TOK_LEN) :: tok(MAX_TOK)
    integer :: nlines, ntok, l, i, j, k, m
    real(dp) :: val

    phi4 = 0.0_dp
    call read_lines(path, lines, nlines)
    do l = 1, nlines
      call split_tokens(lines(l), tok, ntok)
      if (ntok /= 5) cycle
      if (.not. (is_int_token(tok(1)) .and. is_int_token(tok(2)) .and. &
                 is_int_token(tok(3)) .and. is_int_token(tok(4)) .and. &
                 is_real_token(tok(5)))) cycle
      read(tok(1), *) i
      read(tok(2), *) j
      read(tok(3), *) k
      read(tok(4), *) m
      val = token_to_real(tok(5))
      if (min(i, k) < 1 .or. max(i, k) > n) cycle
      phi4(i, k) = val
      phi4(k, i) = val
    end do
  end subroutine parse_quartic_tensor

  !===================================================================
  ! Python / numpy style number formatting used for the printout
  !===================================================================
  !-------------------------------------------------------------------
  ! Python-style "{x:W.Df}": fixed D decimals, right-justified in a
  ! field of at least W characters. Unlike a Fortran Fw.d descriptor the
  ! field grows instead of printing asterisks when the number is wider.
  !-------------------------------------------------------------------
  function pyf(x, w, d) result(s)
    real(dp), intent(in) :: x
    integer, intent(in) :: w, d
    character(len=:), allocatable :: s
    character(len=64) :: buf, fmt
    write(fmt, '(A,I0,A)') '(F0.', d, ')'
    write(buf, fmt) x
    buf = fix_f0(buf)
    if (len_trim(buf) >= w) then
      s = trim(buf)
    else
      s = repeat(' ', w - len_trim(buf)) // trim(buf)
    end if
  end function pyf

  !-------------------------------------------------------------------
  ! Python-style "{x: .Nf}" : fixed N decimals, leading blank for
  ! non-negative numbers, no field padding.
  !-------------------------------------------------------------------
  function fmt_signed(x, nd) result(str)
    real(dp), intent(in) :: x
    integer, intent(in) :: nd
    character(len=48) :: str, buf, fmt
    write(fmt, '(A,I0,A)') '(F0.', nd, ')'
    write(buf, fmt) x
    buf = fix_f0(buf)
    if (buf(1:1) == '-') then
      str = buf
    else
      str = ' ' // trim(buf)
    end if
  end function fmt_signed

  ! ensure a leading zero before the decimal point (".5" -> "0.5")
  function fix_f0(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    character(len=len(s)) :: u
    u = adjustl(s)
    if (u(1:1) == '.') then
      t = '0' // trim(u)
    else if (u(1:2) == '-.') then
      t = '-0' // trim(u(2:))
    else
      t = u
    end if
  end function fix_f0

  !-------------------------------------------------------------------
  ! Python repr() of a double: shortest decimal string that round-trips,
  ! fixed notation for 1e-4 <= |x| < 1e16, exponential otherwise
  ! (repr switches to exponent form when the decimal exponent is
  ! below -4 or at least 16).
  !-------------------------------------------------------------------
  function py_repr(x) result(str)
    real(dp), intent(in) :: x
    character(len=40) :: str
    character(len=48) :: buf, fmt
    character(len=24) :: digits, expstr
    real(dp) :: y
    integer :: p, ios, e, nd, epos, i
    logical :: neg

    if (x == 0.0_dp) then
      str = '0.0'
      if (sign(1.0_dp, x) < 0.0_dp) str = '-0.0'
      return
    end if
    if (x /= x) then
      str = 'nan'
      return
    end if

    do p = 1, 17
      write(fmt, '(A,I0,A)') '(ES40.', p - 1, 'E3)'
      write(buf, fmt) x
      read(buf, *, iostat=ios) y
      if (ios == 0 .and. y == x) exit
    end do
    if (p > 17) p = 17

    buf = adjustl(buf)
    neg = buf(1:1) == '-'
    if (neg) buf = buf(2:)
    epos = index(buf, 'E')
    if (epos == 0) epos = index(buf, 'e')
    ! mantissa digits without the decimal point
    digits = ''
    nd = 0
    do i = 1, epos - 1
      if (buf(i:i) >= '0' .and. buf(i:i) <= '9') then
        nd = nd + 1
        digits(nd:nd) = buf(i:i)
      end if
    end do
    read(buf(epos + 1:), *) e
    ! strip trailing zeros (keep at least one digit)
    do while (nd > 1 .and. digits(nd:nd) == '0')
      nd = nd - 1
    end do

    if (e >= -4 .and. e < 16) then
      if (e >= 0) then
        if (nd >= e + 2) then
          str = digits(1:e + 1) // '.' // digits(e + 2:nd)
        else
          str = digits(1:nd) // repeat('0', e + 1 - nd) // '.0'
        end if
      else
        str = '0.' // repeat('0', -e - 1) // digits(1:nd)
      end if
    else
      if (nd > 1) then
        str = digits(1:1) // '.' // digits(2:nd)
      else
        str = digits(1:1)
      end if
      if (e < 0) then
        write(expstr, '(A,I2.2)') 'e-', -e
      else
        write(expstr, '(A,I2.2)') 'e+', e
      end if
      str = trim(str) // trim(expstr)
    end if
    if (neg) str = '-' // trim(str)
  end function py_repr

  !-------------------------------------------------------------------
  ! numpy-style print of a small 1-D float array, e.g.
  !   [0.90107718 0.0343634  0.03351339]
  ! (precision 8, trailing zeros trimmed, columns aligned)
  !-------------------------------------------------------------------
  function np_array_str(v) result(str)
    real(dp), intent(in) :: v(:)
    character(len=256) :: str
    character(len=32) :: parts(size(v)), buf
    integer :: a, dot, ndec, maxdec, maxint, lint, pos
    real(dp) :: amax, amin
    logical :: sci

    amax = maxval(abs(v))
    amin = minval(abs(v), mask=(v /= 0.0_dp))
    sci = .false.
    if (amax >= 1.0e8_dp) sci = .true.
    if (any(v /= 0.0_dp)) then
      if (amin < 1.0e-4_dp .or. amax / amin > 1000.0_dp) sci = .true.
    end if

    if (sci) then
      do a = 1, size(v)
        write(parts(a), '(ES15.8E2)') v(a)
        parts(a) = adjustl(parts(a))
      end do
    else
      maxdec = 0
      maxint = 0
      do a = 1, size(v)
        write(buf, '(F0.8)') v(a)
        buf = fix_f0(buf)
        dot = index(buf, '.')
        ndec = len_trim(buf) - dot
        do while (ndec > 1 .and. buf(dot + ndec:dot + ndec) == '0')
          ndec = ndec - 1
        end do
        parts(a) = buf(1:dot + ndec)
        maxdec = max(maxdec, ndec)
        maxint = max(maxint, dot - 1)
      end do
      do a = 1, size(v)
        dot = index(parts(a), '.')
        lint = dot - 1
        ndec = len_trim(parts(a)) - dot
        buf = repeat(' ', maxint - lint) // trim(parts(a)) // repeat(' ', maxdec - ndec)
        parts(a) = buf
      end do
    end if

    str = '['
    pos = 2
    do a = 1, size(v)
      if (a > 1) then
        str(pos:pos) = ' '
        pos = pos + 1
      end if
      if (sci) then
        buf = trim(parts(a))
      else
        buf = parts(a)(1:maxint + 1 + maxdec)
      end if
      lint = len_trim(buf)
      if (.not. sci) lint = maxint + 1 + maxdec
      str(pos:pos + lint - 1) = buf(1:lint)
      pos = pos + lint
    end do
    str(pos:pos) = ']'
  end function np_array_str

end module gvpt2_io
