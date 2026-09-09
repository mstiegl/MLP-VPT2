!=====================================================================
! gvpt2_polyad
!
! GVPT2 resonance polyads: Fortran translation of the State,
! Interaction, Polyad and fermi_solver code in gvpt2.py (which was
! adapted from PyVPT2, Copyright 2021-2024 Philip Nelson, BSD 3-Clause,
! see LICENSES/PyVPT2_LICENSE.md).
!
! A vibrational state is either a fundamental |k> (a=k, b=0) or a
! two-quantum state |ii> / |ij> (a=i, b=j with i<=j).  Each Fermi
! resonance is an interaction between a fundamental (left) and an
! overtone/combination (right) with cubic coupling phi and type 1 or 2.
! Interactions sharing a state are grouped into polyads whose small
! effective Hamiltonians are diagonalised (LAPACK dsyev); every basis
! state is then assigned to an eigenvalue by maximising the total
! squared overlap, i.e. a linear sum assignment problem solved here
! with the Hungarian algorithm (scipy.optimize.linear_sum_assignment
! in the Python original).
!=====================================================================
module gvpt2_polyad
  use gvpt2_constants, only: dp
  use gvpt2_io, only: pyf
  use iso_fortran_env, only: error_unit
  implicit none
  private

  public :: vstate, interaction_t, polyad_t, state_energy_t
  public :: state_equal, state_less, state_to_string
  public :: fermi_solver, hungarian

  type :: vstate
    integer :: a = 0
    integer :: b = 0        ! 0 for a fundamental (a,), else (a, b)
  end type vstate

  type :: interaction_t
    type(vstate) :: left            ! fundamental |k>
    type(vstate) :: right           ! overtone |ii> or combination |ij>
    real(dp) :: nu_left = 0.0_dp    ! deperturbed energies (cm^-1)
    real(dp) :: nu_right = 0.0_dp
    real(dp) :: phi = 0.0_dp        ! cubic constant phi_kij (cm^-1)
    integer :: ftype = 0            ! 1: 2 w_i ~ w_k,  2: w_i + w_j ~ w_k
  end type interaction_t

  type :: polyad_t
    integer :: nstates = 0
    type(vstate), allocatable :: states(:)      ! basis states (insertion order)
    real(dp), allocatable :: nu(:)              ! diagonal energies
    integer :: nint = 0
    integer, allocatable :: ileft(:), iright(:) ! coupling end points (indices into states)
    real(dp), allocatable :: phi(:)
    integer, allocatable :: ftype(:)
  end type polyad_t

  type :: state_energy_t
    type(vstate) :: state
    real(dp) :: energy = 0.0_dp
  end type state_energy_t

contains

  !===================================================================
  ! state helpers
  !===================================================================
  pure function state_equal(s, t) result(eq)
    type(vstate), intent(in) :: s, t
    logical :: eq
    eq = (s%a == t%a) .and. (s%b == t%b)
  end function state_equal

  ! Ordering used when gvpt2.py prints the solved state map:
  ! key = (len(state), state) -> fundamentals first, then by tuple.
  pure function state_less(s, t) result(lt)
    type(vstate), intent(in) :: s, t
    logical :: lt
    integer :: ls, ltt
    ls = 1
    if (s%b /= 0) ls = 2
    ltt = 1
    if (t%b /= 0) ltt = 2
    if (ls /= ltt) then
      lt = ls < ltt
    else if (s%a /= t%a) then
      lt = s%a < t%a
    else
      lt = s%b < t%b
    end if
  end function state_less

  ! Python str() of the state tuple: "(8,)" or "(5, 5)".
  function state_to_string(s) result(str)
    type(vstate), intent(in) :: s
    character(len=32) :: str
    character(len=16) :: sa, sb
    write(sa, '(I0)') s%a
    if (s%b == 0) then
      str = '(' // trim(sa) // ',)'
    else
      write(sb, '(I0)') s%b
      str = '(' // trim(sa) // ', ' // trim(sb) // ')'
    end if
  end function state_to_string

  !===================================================================
  ! polyad construction
  !===================================================================
  pure function polyad_find(p, s) result(idx)
    type(polyad_t), intent(in) :: p
    type(vstate), intent(in) :: s
    integer :: idx, i
    idx = 0
    do i = 1, p%nstates
      if (state_equal(p%states(i), s)) then
        idx = i
        return
      end if
    end do
  end function polyad_find

  ! Add a state if new; record (or overwrite) its energy. Returns index.
  subroutine polyad_put_state(p, s, nu, idx)
    type(polyad_t), intent(inout) :: p
    type(vstate), intent(in) :: s
    real(dp), intent(in) :: nu
    integer, intent(out) :: idx
    type(vstate), allocatable :: tmp_s(:)
    real(dp), allocatable :: tmp_nu(:)

    idx = polyad_find(p, s)
    if (idx == 0) then
      allocate(tmp_s(p%nstates + 1), tmp_nu(p%nstates + 1))
      if (p%nstates > 0) then
        tmp_s(1:p%nstates) = p%states(1:p%nstates)
        tmp_nu(1:p%nstates) = p%nu(1:p%nstates)
      end if
      p%nstates = p%nstates + 1
      tmp_s(p%nstates) = s
      tmp_nu(p%nstates) = nu
      call move_alloc(tmp_s, p%states)
      call move_alloc(tmp_nu, p%nu)
      idx = p%nstates
    else
      p%nu(idx) = nu        ! Python: nu_list.update({state: nu})
    end if
  end subroutine polyad_put_state

  ! Polyad.add(): register both states and the coupling. A coupling for
  ! an existing (left, right) pair is overwritten (Python dict update).
  subroutine polyad_add(p, inter)
    type(polyad_t), intent(inout) :: p
    type(interaction_t), intent(in) :: inter
    integer :: il, ir, k
    integer, allocatable :: tl(:), tr(:), tf(:)
    real(dp), allocatable :: tp(:)

    call polyad_put_state(p, inter%left, inter%nu_left, il)
    call polyad_put_state(p, inter%right, inter%nu_right, ir)

    do k = 1, p%nint
      if (p%ileft(k) == il .and. p%iright(k) == ir) then
        p%phi(k) = inter%phi
        p%ftype(k) = inter%ftype
        return
      end if
    end do

    allocate(tl(p%nint + 1), tr(p%nint + 1), tf(p%nint + 1), tp(p%nint + 1))
    if (p%nint > 0) then
      tl(1:p%nint) = p%ileft(1:p%nint)
      tr(1:p%nint) = p%iright(1:p%nint)
      tf(1:p%nint) = p%ftype(1:p%nint)
      tp(1:p%nint) = p%phi(1:p%nint)
    end if
    p%nint = p%nint + 1
    tl(p%nint) = il
    tr(p%nint) = ir
    tf(p%nint) = inter%ftype
    tp(p%nint) = inter%phi
    call move_alloc(tl, p%ileft)
    call move_alloc(tr, p%iright)
    call move_alloc(tf, p%ftype)
    call move_alloc(tp, p%phi)
  end subroutine polyad_add

  !===================================================================
  ! Polyad.solve(): build and diagonalise the effective Hamiltonian,
  ! print eigenvalues/CI coefficients, and assign each basis state to
  ! an eigenvalue by maximum total squared overlap.
  !===================================================================
  subroutine polyad_solve(p, energies)
    type(polyad_t), intent(in) :: p
    real(dp), intent(out) :: energies(p%nstates)   ! energy assigned to states(i)
    integer :: dim, i, j, k, info, lwork, root
    real(dp), allocatable :: H(:,:), evals(:), work(:), cost(:,:)
    integer, allocatable :: assign(:)
    real(dp) :: wq(1), coef
    character(len=32) :: sname

    dim = p%nstates
    allocate(H(dim, dim), evals(dim), cost(dim, dim), assign(dim))
    H = 0.0_dp
    do i = 1, dim
      H(i, i) = p%nu(i)
    end do
    do k = 1, p%nint
      i = p%ileft(k)
      j = p%iright(k)
      if (p%ftype(k) == 1) then
        H(i, j) = 1.0_dp / 4.0_dp * p%phi(k)
        H(j, i) = H(i, j)
      else if (p%ftype(k) == 2) then
        H(i, j) = 1.0_dp / (sqrt(2.0_dp) * 2.0_dp) * p%phi(k)
        H(j, i) = H(i, j)
      end if
    end do

    ! symmetric eigenproblem; eigenvalues ascending, eigenvectors in columns
    lwork = -1
    call dsyev('V', 'U', dim, H, dim, evals, wq, lwork, info)
    lwork = max(1, int(wq(1)))
    allocate(work(lwork))
    call dsyev('V', 'U', dim, H, dim, evals, work, lwork, info)
    if (info /= 0) then
      write(error_unit,'(A,I0)') 'ERROR: dsyev failed in polyad diagonalisation, info = ', info
      stop 1
    end if

    ! ---- printout (same layout as gvpt2.py) ----
    write(*,'(A)') ''
    write(*,'(A)') '==================================================='
    write(*,'(A)') 'GVPT2 Polyad Prints Eigen values and vecs'
    write(*,'(A)') '==================================================='
    write(*,'(A)') 'Basis states:'
    do i = 1, dim
      sname = state_to_string(p%states(i))
      write(*,'(A,I0,A)') '  ', i - 1, ': ' // trim(sname)
    end do
    write(*,'(A)') ''
    write(*,'(A)') 'Eigenvalues and CI coefficients:'
    do root = 1, dim
      write(*,'(A)') ''
      write(*,'(A,I0,A,A,A)') 'Root ', root, ': ', pyf(evals(root), 12, 6), ' cm^-1'
      write(*,'(A)') '                    CI coef      coef**2'
      do i = 1, dim
        coef = H(i, root)
        sname = state_to_string(p%states(i))
        write(*,'(A,A,A,A,A,A)') '   ', sname(1:max(12, len_trim(sname))), '  ', pyf(coef, 10, 6), &
             '   ', pyf(coef**2, 10, 6)
      end do
    end do
    write(*,'(A)') '==================================================='

    ! ---- assignment: maximise sum of squared overlaps ----
    do i = 1, dim
      do j = 1, dim
        cost(i, j) = -(H(i, j)**2)      ! rows: basis states, cols: eigenstates
      end do
    end do
    call hungarian(cost, dim, assign)
    do i = 1, dim
      energies(i) = evals(assign(i))
    end do
  end subroutine polyad_solve

  !===================================================================
  ! fermi_solver(): group interactions into polyads and solve them.
  !
  ! Mirrors the Python loop exactly: an interaction is added to every
  ! existing polyad that already contains its left or right state; a
  ! new polyad is created only when none does. Polyads are then solved
  ! in order and their results merged into one state -> energy map,
  ! later polyads overwriting earlier entries for a shared state.
  !===================================================================
  subroutine fermi_solver(inters, ninter, smap, nmap)
    type(interaction_t), intent(in) :: inters(:)
    integer, intent(in) :: ninter
    type(state_energy_t), allocatable, intent(out) :: smap(:)
    integer, intent(out) :: nmap

    type(polyad_t), allocatable :: polyads(:)
    integer :: npoly, ip, k, i, j
    logical :: flag
    real(dp), allocatable :: energies(:)

    allocate(polyads(max(ninter, 1)))
    npoly = 0
    do k = 1, ninter
      flag = .false.
      do ip = 1, npoly
        if (polyad_find(polyads(ip), inters(k)%left) /= 0) then
          call polyad_add(polyads(ip), inters(k))
          flag = .true.
        else if (polyad_find(polyads(ip), inters(k)%right) /= 0) then
          call polyad_add(polyads(ip), inters(k))
          flag = .true.
        end if
      end do
      if (.not. flag) then
        npoly = npoly + 1
        call polyad_add(polyads(npoly), inters(k))
      end if
    end do

    ! upper bound on the number of distinct states
    nmap = 0
    do ip = 1, npoly
      nmap = nmap + polyads(ip)%nstates
    end do
    allocate(smap(max(nmap, 1)))
    nmap = 0

    do ip = 1, npoly
      allocate(energies(polyads(ip)%nstates))
      call polyad_solve(polyads(ip), energies)
      do i = 1, polyads(ip)%nstates
        ! dict.update semantics: overwrite if present, else append
        do j = 1, nmap
          if (state_equal(smap(j)%state, polyads(ip)%states(i))) exit
        end do
        if (j > nmap) then
          nmap = nmap + 1
          smap(nmap)%state = polyads(ip)%states(i)
        end if
        smap(j)%energy = energies(i)
      end do
      deallocate(energies)
    end do
  end subroutine fermi_solver

  !===================================================================
  ! Hungarian (Kuhn-Munkres) algorithm for the square linear sum
  ! assignment problem: minimise sum_i cost(i, assign(i)) over
  ! permutations. O(n^3), potentials formulation.
  !===================================================================
  subroutine hungarian(cost, n, assign)
    integer, intent(in) :: n
    real(dp), intent(in) :: cost(n, n)
    integer, intent(out) :: assign(n)      ! row i -> column assign(i)
    real(dp) :: u(0:n), v(0:n), minv(0:n), delta, cur
    integer :: p(0:n), way(0:n)
    logical :: used(0:n)
    integer :: i, j, i0, j0, j1

    u = 0.0_dp
    v = 0.0_dp
    p = 0
    way = 0
    do i = 1, n
      p(0) = i
      j0 = 0
      minv = huge(1.0_dp)
      used = .false.
      do
        used(j0) = .true.
        i0 = p(j0)
        delta = huge(1.0_dp)
        j1 = 0
        do j = 1, n
          if (.not. used(j)) then
            cur = cost(i0, j) - u(i0) - v(j)
            if (cur < minv(j)) then
              minv(j) = cur
              way(j) = j0
            end if
            if (minv(j) < delta) then
              delta = minv(j)
              j1 = j
            end if
          end if
        end do
        do j = 0, n
          if (used(j)) then
            u(p(j)) = u(p(j)) + delta
            v(j) = v(j) - delta
          else
            minv(j) = minv(j) - delta
          end if
        end do
        j0 = j1
        if (p(j0) == 0) exit
      end do
      do
        j1 = way(j0)
        p(j0) = p(j1)
        j0 = j1
        if (j0 == 0) exit
      end do
    end do

    assign = 0
    do j = 1, n
      if (p(j) /= 0) assign(p(j)) = j
    end do
  end subroutine hungarian

end module gvpt2_polyad
