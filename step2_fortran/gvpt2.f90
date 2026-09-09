!=====================================================================
! gvpt2.x  --  Fortran translation of MLP-VPT2 step2/gvpt2.py
!
! Usage:  ./gvpt2.x --mol_name=MOLECULE [--num_modes=N] [--enable_intensity]
!         ./gvpt2.x MOLECULE [N]
!
! MOLECULE : prefix of the step-1 output files
!            MOLECULE_opt.out, MOLECULE_nma.out,
!            MOLECULE_cubics.out, MOLECULE_quartics.out
! N        : number of vibrational modes kept (default 3*N_atoms - 6)
!
! Computes VPT2 fundamentals, deperturbed (DVPT2) fundamentals and
! GVPT2 fundamentals with Fermi resonances treated variationally in
! small effective Hamiltonians, plus chi0 and the anharmonic ZPVE.
! Output goes to standard output in the same layout as gvpt2.py.
!
! The optional double-harmonic IR/Raman intensities of gvpt2.py require
! the MACE-MDP PyTorch model and are NOT available in this Fortran
! version; --enable_intensity is accepted but only prints a notice.
!
! Workflow (numbering follows gvpt2.py):
!   0. read geometry, harmonic frequencies and normal modes;
!      rotational constants and Coriolis zeta
!   2.1 omega -> Eh for the Q -> q scaling
!   2.2 read Phi(Q) cubic and semi-diagonal quartic constants
!   2.3 Q -> dimensionless q, Eh -> cm^-1
!   2.4 Fermi-1 / Fermi-2 resonance detection
!   2.5 chi (VPT2, no deperturbation) and chi (DVPT2)
!   2.6 deperturbed fundamentals
!   2.7 chi0 and ZPVE
!   2.8 GVPT2 polyads
!=====================================================================
program gvpt2
  use gvpt2_constants
  use gvpt2_io
  use gvpt2_polyad
  use gvpt2_core
  use iso_fortran_env, only: error_unit
  implicit none

  character(len=256) :: mol_name, arg, val
  character(len=512) :: xyz_file, nma_file, cubic_file, quartic_file
  integer :: nargs, ia, n_modes, ios, natoms, natom_nma, n, i, j, k, nmap
  logical :: enable_intensity, have_mol, next_is_mol, next_is_nm

  character(len=8), allocatable :: symbols(:)
  real(dp), allocatable :: pos_ang(:,:), x_bohr(:,:), masses_me(:)
  real(dp), allocatable :: omega_cm(:), omega_eh(:), modes(:,:,:), zeta(:,:,:)
  real(dp), allocatable :: phi3_Qc(:,:,:), phi4_Qc(:,:), phi3_cm(:,:,:), phi4_cm(:,:)
  real(dp), allocatable :: chi_vpt2(:,:), chi_cm(:,:)
  real(dp), allocatable :: nu_vpt2(:), nu_depert(:), nu_gvpt2(:)
  logical, allocatable :: fermi(:,:,:), no_fermi(:,:,:)
  type(state_energy_t), allocatable :: smap(:)
  real(dp) :: rotc(3), chi0_cm, harmonic_zpve, anharmonic_zpve, zpve_corr
  integer :: nfermi
  character(len=32) :: sname
  character(len=64) :: s1

  ! ------------------------------------------------------------------
  ! command line
  ! ------------------------------------------------------------------
  mol_name = ''
  have_mol = .false.
  n_modes = -1
  enable_intensity = .false.
  next_is_mol = .false.
  next_is_nm = .false.
  nargs = command_argument_count()
  do ia = 1, nargs
    call get_command_argument(ia, arg)
    if (next_is_mol) then
      mol_name = arg
      have_mol = .true.
      next_is_mol = .false.
    else if (next_is_nm) then
      read(arg, *, iostat=ios) n_modes
      if (ios /= 0) call usage_error('bad value for --num_modes: ' // trim(arg))
      next_is_nm = .false.
    else if (arg(1:11) == '--mol_name=') then
      mol_name = arg(12:)
      have_mol = .true.
    else if (trim(arg) == '--mol_name') then
      next_is_mol = .true.
    else if (arg(1:12) == '--num_modes=') then
      val = arg(13:)
      read(val, *, iostat=ios) n_modes
      if (ios /= 0) call usage_error('bad value for --num_modes: ' // trim(val))
    else if (trim(arg) == '--num_modes') then
      next_is_nm = .true.
    else if (trim(arg) == '--enable_intensity') then
      enable_intensity = .true.
    else if (trim(arg) == '-h' .or. trim(arg) == '--help') then
      call print_usage()
      stop
    else if (arg(1:1) == '-') then
      call usage_error('unknown option ' // trim(arg))
    else if (.not. have_mol) then
      mol_name = arg
      have_mol = .true.
    else if (n_modes < 0) then
      read(arg, *, iostat=ios) n_modes
      if (ios /= 0) call usage_error('bad value for num_modes: ' // trim(arg))
    else
      call usage_error('unexpected argument ' // trim(arg))
    end if
  end do
  if (next_is_mol .or. next_is_nm) call usage_error('missing value after last option')
  if (.not. have_mol) call usage_error('--mol_name is required')

  xyz_file     = trim(mol_name) // '_opt.out'
  nma_file     = trim(mol_name) // '_nma.out'
  cubic_file   = trim(mol_name) // '_cubics.out'
  quartic_file = trim(mol_name) // '_quartics.out'

  ! ------------------------------------------------------------------
  ! Step 0: molecule set-up
  ! 0.1) Cartesian coordinates, harmonic frequencies and normal modes
  ! ------------------------------------------------------------------
  call read_xyz(xyz_file, natoms, symbols, pos_ang)
  if (n_modes < 0) n_modes = 3 * natoms - 6
  if (n_modes < 1) call usage_error('number of modes must be positive')
  n = n_modes

  call read_freqs_and_vecs(nma_file, n_modes, omega_cm, modes, natom_nma)
  if (natom_nma /= natoms) then
    write(error_unit,'(A,I0,A,I0,A)') 'ERROR: normal-mode file has ', natom_nma, &
         ' atoms per mode but the geometry file has ', natoms, ' atoms'
    stop 1
  end if

  ! 0.2) rotational constants and Coriolis coupling
  allocate(x_bohr(natoms, 3), masses_me(natoms))
  x_bohr = pos_ang / BOHR_TO_ANG
  do i = 1, natoms
    masses_me(i) = atomic_mass_amu(symbols(i)) * AMU_TO_ME
  end do
  rotc = rotational_constants(x_bohr, masses_me) * CM_PER_EH
  write(*,'(A)') 'Rotational constants (cm-1):'
  write(*,'(A)') trim(np_array_str(rotc))

  allocate(zeta(3, n, n))
  call coriolis_zeta(modes, zeta)

  ! Step 1 (optional) in gvpt2.py: MACE-MDP intensities -- not available here
  if (enable_intensity) then
    write(*,'(A)') ''
    write(*,'(A)') 'NOTE: --enable_intensity is not supported by the Fortran version;'
    write(*,'(A)') '      IR/Raman intensities need the MACE-MDP model (use gvpt2.py).'
    write(*,'(A)') '      Continuing with the energy calculation only.'
  end if

  ! ------------------------------------------------------------------
  ! Step 2: anharmonic frequencies (DVPT2/GVPT2)
  ! ------------------------------------------------------------------
  ! 2.1) omega in Eh for the Q -> q scaling only
  allocate(omega_eh(n))
  omega_eh = omega_cm * EH_PER_CM

  ! 2.2) Phi(Q) in Eh/Q^3, Eh/Q^4
  allocate(phi3_Qc(n, n, n), phi4_Qc(n, n))
  call parse_cubic_tensor(cubic_file, n, phi3_Qc)
  call parse_quartic_tensor(quartic_file, n, phi4_Qc)

  write(*,'(A)') ''
  write(*,'(A)') 'RAW TENSOR MAGNITUDES (before any scaling)'
  call stats3('Phi3_Q', phi3_Qc)
  call stats2('Phi4_Q', phi4_Qc)
  write(*,'(A)') 'omega_cm min/max: ' // trim(py_repr(minval(omega_cm))) // ' ' // trim(py_repr(maxval(omega_cm)))

  ! 2.3) Q -> dimensionless q, in cm^-1
  allocate(phi3_cm(n, n, n), phi4_cm(n, n))
  call convert_Q_to_q(omega_eh, phi3_Qc, phi4_Qc, phi3_cm, phi4_cm)

  write(*,'(A)') ''
  write(*,'(A)') 'RAW TENSOR MAGNITUDES (after Q->q scaling, still Eh)'
  call stats3('Phi3_q_cm', phi3_cm)
  call stats2('Phi4_q_cm', phi4_cm)

  ! 2.4) Fermi resonances (for deperturbation + GVPT2)
  allocate(fermi(n, n, n), no_fermi(n, n, n))
  call identify_fermi_resonances(omega_cm, phi3_cm, FERMI_OMEGA_THRESH, FERMI_K_THRESH, fermi)
  no_fermi = .false.

  nfermi = count(fermi)
  write(*,'(A)') ''
  write(*,'(A,I0,A)') 'Detected ', nfermi, ' Fermi candidates:'
  do k = 1, n
    do i = 1, n
      do j = 1, n
        if (.not. fermi(k, i, j)) cycle
        if (i == j) then
          write(s1, '(F0.3)') abs(2.0_dp * omega_cm(i) - omega_cm(k))
          write(*,'(A,I0,A,I0,A,I0,A)') '  (k=', k, ', (i,i)=(', i, ',', i, '))  |2ω_i-ω_k|=' // &
               trim(fix_f0(s1)) // ' cm^-1'
        else
          write(s1, '(F0.3)') abs(omega_cm(i) + omega_cm(j) - omega_cm(k))
          write(*,'(A,I0,A,I0,A,I0,A)') '  (k=', k, ', (i,j)=(', i, ',', j, '))  |ω_i+ω_j-ω_k|=' // &
               trim(fix_f0(s1)) // ' cm^-1'
        end if
      end do
    end do
  end do

  ! 2.5a) ordinary VPT2 chi: no deperturbation
  allocate(chi_vpt2(n, n), chi_cm(n, n), nu_vpt2(n), nu_depert(n), nu_gvpt2(n))
  call compute_chi_deperturbed(omega_cm, phi3_cm, phi4_cm, rotc, zeta, no_fermi, chi_vpt2)
  call fundamentals_from_chi(omega_cm, chi_vpt2, nu_vpt2)

  ! 2.5b) DVPT2 chi (cm^-1)
  call compute_chi_deperturbed(omega_cm, phi3_cm, phi4_cm, rotc, zeta, fermi, chi_cm)

  ! 2.6) deperturbed fundamentals
  call fundamentals_from_chi(omega_cm, chi_cm, nu_depert)

  ! 2.7) zero-point vibrational energy with chi0
  chi0_cm = compute_chi0(omega_cm, phi3_cm, phi4_cm, rotc, zeta, fermi)
  call compute_zpve(omega_cm, chi_cm, chi0_cm, harmonic_zpve, anharmonic_zpve, zpve_corr)

  ! 2.8) GVPT2: variational mixing in resonant polyads
  call apply_gvpt2_fermi(omega_cm, nu_depert, chi_cm, phi3_cm, fermi, nu_gvpt2, smap, nmap)

  ! ------------------------------------------------------------------
  ! summary
  ! ------------------------------------------------------------------
  write(*,'(A)') ''
  write(*,'(A)') 'Harmonic ω (cm^-1):'
  do i = 1, n
    write(*,'(A,I2,A,A)') '  ω_', i, ' = ', pyf(omega_cm(i), 12, 6)
  end do

  write(*,'(A)') ''
  write(*,'(A)') 'Deperturbed χ_ii (cm^-1):'
  do i = 1, n
    write(*,'(A,I2.2,I2.2,A,A)') '  χ_', i, i, ' = ', trim(fmt_signed(chi_cm(i, i), 8))
  end do

  write(*,'(A)') ''
  write(*,'(A)') 'Fundamentals (cm^-1):'
  write(*,'(A)') '  mode     Harmonic     VPT2      DVPT2      GVPT2       GVPT2-HOshift   GVPT2-VPT2shift'
  do i = 1, n
    write(*,'(A,I4,6(A,A))') '  ', i, '  ', pyf(omega_cm(i), 10, 4), '  ', pyf(nu_vpt2(i), 10, 4), &
         '  ', pyf(nu_depert(i), 10, 4), '  ', pyf(nu_gvpt2(i), 10, 4), &
         '  ', pyf(nu_gvpt2(i) - omega_cm(i), 10, 4), '  ', pyf(nu_gvpt2(i) - nu_vpt2(i), 10, 4)
  end do

  if (nmap > 0) then
    call sort_state_map(smap, nmap)
    write(*,'(A)') ''
    write(*,'(A)') 'GVPT2 solved state energies (cm^-1) in resonant polyads (includes combo/overtones):'
    do i = 1, nmap
      sname = state_to_string(smap(i)%state)
      write(*,'(A,A,A,A)') '  ', sname(1:max(12, len_trim(sname))), '  ', pyf(smap(i)%energy, 12, 6)
    end do
  end if

  write(*,'(A)') ''
  write(*,'(A)') 'Zero-Point Vibrational Energy:'
  write(*,'(A)') '  chi0 = ' // trim(fmt_signed(chi0_cm, 8)) // ' cm^-1'
  write(*,'(A)') '  unit          harmonic ZPVE      correction       anharmonic ZPVE'
  call print_zpve_line('cm^-1   ', 1.0_dp)
  call print_zpve_line('kcal/mol', CM_TO_KCAL)
  call print_zpve_line('kJ/mol  ', CM_TO_KJ)

  write(*,'(A)') ''
  write(*,'(A)') 'Done.'
  write(*,'(A)') ''
  write(*,'(A)') 'GVPT2 code run is now complete'

contains

  !-------------------------------------------------------------------
  subroutine print_usage()
    write(*,'(A)') 'Usage: gvpt2.x --mol_name=MOLECULE [--num_modes=N] [--enable_intensity]'
    write(*,'(A)') '       gvpt2.x MOLECULE [N]'
    write(*,'(A)') ''
    write(*,'(A)') '  MOLECULE  prefix of MOLECULE_opt.out, _nma.out, _cubics.out, _quartics.out'
    write(*,'(A)') '  N         number of vibrational modes (default 3*N_atoms-6)'
    write(*,'(A)') '  --enable_intensity  accepted for compatibility; intensities need gvpt2.py'
  end subroutine print_usage

  subroutine usage_error(msg)
    character(len=*), intent(in) :: msg
    write(error_unit,'(A)') 'ERROR: ' // trim(msg)
    call print_usage()
    stop 1
  end subroutine usage_error

  !-------------------------------------------------------------------
  subroutine print_zpve_line(unit_name, factor)
    character(len=*), intent(in) :: unit_name
    real(dp), intent(in) :: factor
    write(*,'(A,A,A,A,A,A,A,A)') '  ', unit_name, '  ', pyf(factor * harmonic_zpve, 14, 6), &
         '  ', pyf(factor * zpve_corr, 14, 6), '  ', pyf(factor * anharmonic_zpve, 14, 6)
  end subroutine print_zpve_line

  !-------------------------------------------------------------------
  ! sort the solved state map like sorted(state_map.items(),
  ! key=lambda x: (len(x[0]), x[0]))
  !-------------------------------------------------------------------
  subroutine sort_state_map(m, nm)
    type(state_energy_t), intent(inout) :: m(:)
    integer, intent(in) :: nm
    type(state_energy_t) :: t
    integer :: a, b
    do a = 2, nm
      t = m(a)
      b = a - 1
      do while (b >= 1)
        if (.not. state_less(t%state, m(b)%state)) exit
        m(b + 1) = m(b)
        b = b - 1
      end do
      m(b + 1) = t
    end do
  end subroutine sort_state_map

  !-------------------------------------------------------------------
  ! debug statistics of a tensor, like stats() in gvpt2.py:
  ! count / min / median / max of |v| over the nonzero entries
  !-------------------------------------------------------------------
  subroutine stats3(name, t)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: t(:,:,:)
    real(dp), allocatable :: v(:)
    integer :: cnt, a, b, c
    allocate(v(size(t)))
    cnt = 0
    do c = 1, size(t, 3)
      do b = 1, size(t, 2)
        do a = 1, size(t, 1)
          if (t(a, b, c) /= 0.0_dp) then
            cnt = cnt + 1
            v(cnt) = abs(t(a, b, c))
          end if
        end do
      end do
    end do
    call print_stats(name, v, cnt)
  end subroutine stats3

  subroutine stats2(name, t)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: t(:,:)
    real(dp), allocatable :: v(:)
    integer :: cnt, a, b
    allocate(v(size(t)))
    cnt = 0
    do b = 1, size(t, 2)
      do a = 1, size(t, 1)
        if (t(a, b) /= 0.0_dp) then
          cnt = cnt + 1
          v(cnt) = abs(t(a, b))
        end if
      end do
    end do
    call print_stats(name, v, cnt)
  end subroutine stats2

  subroutine print_stats(name, v, cnt)
    character(len=*), intent(in) :: name
    real(dp), intent(inout) :: v(:)
    integer, intent(in) :: cnt
    character(len=16) :: sc
    if (cnt == 0) then
      write(*,'(A)') trim(name) // ' no nonzero values'
      return
    end if
    call sort_reals(v, cnt)
    write(sc, '(I0)') cnt
    write(*,'(A)') trim(name) // ' count = ' // trim(sc) // ' min = ' // trim(py_repr(v(1))) // &
         ' median = ' // trim(py_repr(v(cnt / 2 + 1))) // ' max = ' // trim(py_repr(v(cnt)))
  end subroutine print_stats

  !-------------------------------------------------------------------
  ! in-place heapsort (ascending)
  !-------------------------------------------------------------------
  subroutine sort_reals(v, cnt)
    real(dp), intent(inout) :: v(:)
    integer, intent(in) :: cnt
    integer :: a, b, c
    real(dp) :: t
    do a = cnt / 2, 1, -1
      call sift_down(v, a, cnt)
    end do
    do b = cnt, 2, -1
      t = v(1)
      v(1) = v(b)
      v(b) = t
      c = b - 1
      call sift_down(v, 1, c)
    end do
  end subroutine sort_reals

  subroutine sift_down(v, start, last)
    real(dp), intent(inout) :: v(:)
    integer, intent(in) :: start, last
    integer :: root, child
    real(dp) :: t
    root = start
    do while (2 * root <= last)
      child = 2 * root
      if (child < last) then
        if (v(child) < v(child + 1)) child = child + 1
      end if
      if (v(root) < v(child)) then
        t = v(root)
        v(root) = v(child)
        v(child) = t
        root = child
      else
        exit
      end if
    end do
  end subroutine sift_down





end program gvpt2
