!=====================================================================
! gvpt2_core
!
! Physics of the Fortran translation of gvpt2.py:
!   rotational_constants      A, B, C from principal-axis coordinates
!   coriolis_zeta             Coriolis coupling constants zeta^alpha_kl
!   convert_Q_to_q            Phi(Q) [Eh] -> phi(q) [cm^-1], dimensionless q
!   identify_fermi_resonances Fermi-1 (2 w_i ~ w_k) and Fermi-2 (w_i + w_j ~ w_k)
!   compute_chi_deperturbed   anharmonic constants chi_ij (DVPT2 when
!                             resonant terms are removed)
!   fundamentals_from_chi     nu_i = w_i + 2 chi_ii + 1/2 sum_j chi_ij
!   compute_chi0              constant term chi0 for the ZPVE
!   compute_zpve              harmonic / anharmonic ZPVE
!   build_explicit_state_energies  |i>, |ii>, |ij> deperturbed energies
!   apply_gvpt2_fermi         build polyads from the Fermi set and solve
!
! Mode indices are 1-based throughout. The Fermi set is stored as a
! logical array fermi(k, i, j) == .true.  <=>  (k, (i, j)) in fermi_set
! of the Python code; Fermi-2 entries are present for both (i,j) and
! (j,i), exactly as produced by itertools.permutations in gvpt2.py.
!=====================================================================
module gvpt2_core
  use gvpt2_constants, only: dp, CM_PER_EH, EPS_DIV
  use gvpt2_polyad
  use iso_fortran_env, only: error_unit
  implicit none
  private

  public :: rotational_constants, coriolis_zeta, convert_Q_to_q
  public :: identify_fermi_resonances, compute_chi_deperturbed
  public :: fundamentals_from_chi, compute_chi0, compute_zpve
  public :: build_explicit_state_energies, apply_gvpt2_fermi

contains

  !-------------------------------------------------------------------
  ! Rotational constants (Hartree, hbar = 1) from coordinates in the
  ! principal-axis frame (bohr) and masses (electron masses). Only the
  ! diagonal moments of inertia are used, as in gvpt2.py.
  !-------------------------------------------------------------------
  function rotational_constants(x, masses) result(rotc)
    real(dp), intent(in) :: x(:,:)        ! (natom, 3)
    real(dp), intent(in) :: masses(:)     ! (natom)
    real(dp) :: rotc(3)
    real(dp) :: I_x, I_y, I_z
    I_x = sum(masses * (x(:,2)**2 + x(:,3)**2))
    I_y = sum(masses * (x(:,1)**2 + x(:,3)**2))
    I_z = sum(masses * (x(:,1)**2 + x(:,2)**2))
    rotc(1) = 1.0_dp / (2.0_dp * I_x)
    rotc(2) = 1.0_dp / (2.0_dp * I_y)
    rotc(3) = 1.0_dp / (2.0_dp * I_z)
  end function rotational_constants

  !-------------------------------------------------------------------
  ! Coriolis coupling constants zeta(alpha, k, l) from mass-weighted
  ! normal modes vec(k, atom, xyz).
  !-------------------------------------------------------------------
  subroutine coriolis_zeta(vec, zeta)
    real(dp), intent(in) :: vec(:,:,:)              ! (nmode, natom, 3)
    real(dp), intent(out) :: zeta(3, size(vec,1), size(vec,1))
    integer :: nmode, k, l
    nmode = size(vec, 1)
    zeta = 0.0_dp
    do k = 1, nmode
      do l = 1, nmode
        zeta(1, k, l) = sum(vec(k,:,2) * vec(l,:,3) - vec(k,:,3) * vec(l,:,2))
        zeta(2, k, l) = sum(vec(k,:,3) * vec(l,:,1) - vec(k,:,1) * vec(l,:,3))
        zeta(3, k, l) = sum(vec(k,:,1) * vec(l,:,2) - vec(k,:,2) * vec(l,:,1))
      end do
    end do
  end subroutine coriolis_zeta

  !-------------------------------------------------------------------
  ! Phi(Q) in Eh/Q^n  ->  phi(q) in cm^-1 with dimensionless q:
  !   phi3(q) = Phi3(Q) / sqrt(w_i w_j w_k) * CM_PER_EH
  !   phi4(q) = Phi4(Q) / (w_i w_k)         * CM_PER_EH
  ! with omega_eh the harmonic frequencies in Hartree.
  !-------------------------------------------------------------------
  subroutine convert_Q_to_q(omega_eh, phi3_in, phi4_in, phi3_q, phi4_q)
    real(dp), intent(in) :: omega_eh(:)
    real(dp), intent(in) :: phi3_in(:,:,:), phi4_in(:,:)   ! Phi(Q), Eh
    real(dp), intent(out) :: phi3_q(size(omega_eh), size(omega_eh), size(omega_eh))
    real(dp), intent(out) :: phi4_q(size(omega_eh), size(omega_eh))
    real(dp) :: sqrtw(size(omega_eh))
    integer :: n, i, j, k

    n = size(omega_eh)
    do i = 1, n
      if (omega_eh(i) < 0.0_dp) then
        write(error_unit,'(A,I0,A)') 'ERROR: negative (imaginary) harmonic frequency for mode ', i, &
             '; cannot convert to dimensionless coordinates'
        stop 1
      end if
      sqrtw(i) = sqrt(omega_eh(i))
    end do

    do k = 1, n
      do j = 1, n
        do i = 1, n
          phi3_q(i, j, k) = phi3_in(i, j, k) / (sqrtw(i) * sqrtw(j) * sqrtw(k)) * CM_PER_EH
        end do
      end do
    end do

    do k = 1, n
      do i = 1, n
        phi4_q(i, k) = phi4_in(i, k) / (sqrtw(i)**2 * sqrtw(k)**2) * CM_PER_EH
      end do
    end do
  end subroutine convert_Q_to_q

  !-------------------------------------------------------------------
  ! Heuristic Fermi resonance detection.
  !   Fermi-1: |2 w_i - w_k| <= omega_thresh and
  !            phi_iik^4 / (256 d^3) >= K_thresh      -> fermi(k, i, i)
  !   Fermi-2: |w_i + w_j - w_k| <= omega_thresh and
  !            phi_ijk^4 / (64 d^3) >= K_thresh       -> fermi(k, i, j)
  !-------------------------------------------------------------------
  subroutine identify_fermi_resonances(omega_cm, phi3_cm, omega_thresh, K_thresh, fermi)
    real(dp), intent(in) :: omega_cm(:)
    real(dp), intent(in) :: phi3_cm(:,:,:)
    real(dp), intent(in) :: omega_thresh, K_thresh
    logical, intent(out) :: fermi(size(omega_cm), size(omega_cm), size(omega_cm))
    integer :: n, i, j, k
    real(dp) :: d_omega, phi, dK

    n = size(omega_cm)
    fermi = .false.

    ! Fermi-1: 2*i ~ k   (ordered pairs i /= k)
    do i = 1, n
      do k = 1, n
        if (k == i) cycle
        d_omega = abs(2.0_dp * omega_cm(i) - omega_cm(k))
        if (d_omega <= omega_thresh) then
          phi = phi3_cm(i, i, k)
          dK = (phi**4) / (256.0_dp * (d_omega**3 + 1.0e-30_dp))
          if (dK >= K_thresh) fermi(k, i, i) = .true.
        end if
      end do
    end do

    ! Fermi-2: i + j ~ k   (ordered triples, all distinct)
    do i = 1, n
      do j = 1, n
        if (j == i) cycle
        do k = 1, n
          if (k == i .or. k == j) cycle
          d_omega = abs(omega_cm(i) + omega_cm(j) - omega_cm(k))
          if (d_omega <= omega_thresh) then
            phi = phi3_cm(i, j, k)
            dK = (phi**4) / (64.0_dp * (d_omega**3 + 1.0e-30_dp))
            if (dK >= K_thresh) fermi(k, i, j) = .true.
          end if
        end do
      end do
    end do
  end subroutine identify_fermi_resonances

  !-------------------------------------------------------------------
  pure function delta_ijk(wi, wj, wk) result(d)
    real(dp), intent(in) :: wi, wj, wk
    real(dp) :: d
    d = (wi + wj + wk) * (wi + wj - wk) * (wi - wj + wk) * (wi - wj - wk)
  end function delta_ijk

  !-------------------------------------------------------------------
  ! Anharmonic constants chi_ij (cm^-1). With an empty Fermi set this
  ! is ordinary VPT2; with resonant terms flagged the near-singular
  ! denominators are dropped (DVPT2).
  !-------------------------------------------------------------------
  subroutine compute_chi_deperturbed(omega_cm, phi3_cm, phi4_cm, rotc, zeta, fermi, chi)
    real(dp), intent(in) :: omega_cm(:)
    real(dp), intent(in) :: phi3_cm(:,:,:), phi4_cm(:,:)
    real(dp), intent(in) :: rotc(3)
    real(dp), intent(in) :: zeta(:,:,:)
    logical, intent(in) :: fermi(:,:,:)
    real(dp), intent(out) :: chi(size(omega_cm), size(omega_cm))
    integer :: n, i, j, k
    real(dp) :: wi, wj, wk, term, s, s1, s2, s3, phi, num, den, delta_ij, D

    n = size(omega_cm)
    chi = 0.0_dp

    ! diagonal chi_ii
    do i = 1, n
      wi = omega_cm(i)
      term = phi4_cm(i, i)
      s = 0.0_dp
      do k = 1, n
        wk = omega_cm(k)
        phi = phi3_cm(i, i, k)
        if (fermi(k, i, i)) then
          ! safe form (no 4 wi^2 - wk^2)
          s = s + 0.5_dp * (phi**2) * (1.0_dp / (2.0_dp * wi + wk + EPS_DIV) + 4.0_dp / (wk + EPS_DIV))
        else
          num = (8.0_dp * wi**2 - 3.0_dp * wk**2) * (phi**2)
          den = wk * (4.0_dp * wi**2 - wk**2)
          s = s + num / (den + EPS_DIV)
        end if
      end do
      chi(i, i) = (term - s) / 16.0_dp
    end do

    ! off-diagonal chi_ij
    do i = 1, n
      wi = omega_cm(i)
      do j = 1, n
        if (i == j) cycle
        wj = omega_cm(j)
        term = phi4_cm(i, j)

        ! - sum_k phi_iik phi_jjk / wk
        s1 = 0.0_dp
        do k = 1, n
          wk = omega_cm(k)
          s1 = s1 + (phi3_cm(i, i, k) * phi3_cm(j, j, k)) / (wk + EPS_DIV)
        end do

        ! + sum_k phi_ijk^2 * delta_ij
        s2 = 0.0_dp
        do k = 1, n
          wk = omega_cm(k)
          phi = phi3_cm(i, j, k)
          if (fermi(k, i, j)) then
            ! i + j = k : drop 1/(wi+wj-wk)
            delta_ij = ( 1.0_dp / (wi + wj + wk + EPS_DIV) + &
                         1.0_dp / (-wi + wj + wk + EPS_DIV) + &
                         1.0_dp / (wi - wj + wk + EPS_DIV) ) / (-2.0_dp)
          else if (fermi(i, j, k)) then
            ! j + k = i : drop 1/(-wi+wj+wk)
            delta_ij = ( 1.0_dp / (wi + wj + wk + EPS_DIV) + &
                         1.0_dp / (-wi - wj + wk + EPS_DIV) + &
                         1.0_dp / (wi - wj + wk + EPS_DIV) ) / (-2.0_dp)
          else if (fermi(j, i, k)) then
            ! i + k = j : drop 1/(wi-wj+wk)
            delta_ij = ( 1.0_dp / (wi + wj + wk + EPS_DIV) + &
                         1.0_dp / (-wi - wj + wk + EPS_DIV) + &
                         1.0_dp / (-wi + wj + wk + EPS_DIV) ) / (-2.0_dp)
          else
            D = delta_ijk(wi, wj, wk)
            delta_ij = 2.0_dp * wk * (wi**2 + wj**2 - wk**2) / (D + EPS_DIV)
          end if
          s2 = s2 + (phi**2) * delta_ij
        end do

        ! + 4 (wi^2 + wj^2) / (wi wj) sum_alpha B_alpha zeta^alpha_ij^2
        s3 = 0.0_dp
        do k = 1, 3
          s3 = s3 + 4.0_dp * (wi**2 + wj**2) / (wi * wj) * rotc(k) * zeta(k, i, j)**2
        end do

        chi(i, j) = (term - s1 + s2 + s3) / 4.0_dp
      end do
    end do
  end subroutine compute_chi_deperturbed

  !-------------------------------------------------------------------
  ! Nondegenerate fundamentals nu_i = w_i + 2 chi_ii + 1/2 sum_{j/=i} chi_ij
  !-------------------------------------------------------------------
  subroutine fundamentals_from_chi(omega_cm, chi, nu)
    real(dp), intent(in) :: omega_cm(:), chi(:,:)
    real(dp), intent(out) :: nu(size(omega_cm))
    integer :: n, i, j
    real(dp) :: offsum
    n = size(omega_cm)
    do i = 1, n
      offsum = 0.0_dp
      do j = 1, n
        if (j /= i) offsum = offsum + chi(i, j)
      end do
      nu(i) = omega_cm(i) + 2.0_dp * chi(i, i) + 0.5_dp * offsum
    end do
  end subroutine fundamentals_from_chi

  !-------------------------------------------------------------------
  ! chi0 (cm^-1):
  !   chi0  = sum_i [ phi_iiii - (7/9) phi_iii^2 / w_i
  !                   + sum_{j/=i} 3 w_i phi_ijj^2 / (4 w_j^2 - w_i^2) ]
  !         + sum_{i<j<k} 2 phi_ijk^2 delta_0
  !         - sum_alpha 16 B_alpha (1 + 2 sum_{i<j} zeta^alpha_ij^2)
  !   chi0 /= 64
  !-------------------------------------------------------------------
  function compute_chi0(omega_cm, phi3_cm, phi4_cm, rotc, zeta, fermi) result(chi0)
    real(dp), intent(in) :: omega_cm(:)
    real(dp), intent(in) :: phi3_cm(:,:,:), phi4_cm(:,:)
    real(dp), intent(in) :: rotc(3)
    real(dp), intent(in) :: zeta(:,:,:)
    logical, intent(in) :: fermi(:,:,:)
    real(dp) :: chi0
    integer :: n, i, j, k
    real(dp) :: wi, wj, wk, den, delta_0, D, sr

    n = size(omega_cm)
    chi0 = 0.0_dp

    ! one-mode and two-mode non-rotational terms
    do i = 1, n
      wi = omega_cm(i)
      chi0 = chi0 + phi4_cm(i, i)
      chi0 = chi0 - (7.0_dp / 9.0_dp) * (phi3_cm(i, i, i)**2) / (wi + EPS_DIV)
      do j = 1, n
        if (j == i) cycle
        wj = omega_cm(j)
        den = 4.0_dp * (wj**2) - (wi**2)
        chi0 = chi0 + 3.0_dp * wi * (phi3_cm(i, j, j)**2) / (den + EPS_DIV)
      end do
    end do

    ! three-mode terms, i < j < k
    do i = 1, n
      do j = i + 1, n
        do k = j + 1, n
          wi = omega_cm(i)
          wj = omega_cm(j)
          wk = omega_cm(k)
          if (fermi(k, i, j)) then
            ! i + j = k: drop 1/(wi + wj - wk)
            delta_0 = 1.0_dp / (wi + wj + wk + EPS_DIV) &
                    - 1.0_dp / (wi - wj + wk + EPS_DIV) &
                    - 1.0_dp / (-wi + wj + wk + EPS_DIV)
          else if (fermi(i, j, k)) then
            ! j + k = i: drop 1/(-wi + wj + wk)
            delta_0 = 1.0_dp / (wi + wj + wk + EPS_DIV) &
                    - 1.0_dp / (wi + wj - wk + EPS_DIV) &
                    - 1.0_dp / (wi - wj + wk + EPS_DIV)
          else if (fermi(j, i, k)) then
            ! i + k = j: drop 1/(wi - wj + wk)
            delta_0 = 1.0_dp / (wi + wj + wk + EPS_DIV) &
                    - 1.0_dp / (wi + wj - wk + EPS_DIV) &
                    - 1.0_dp / (-wi + wj + wk + EPS_DIV)
          else
            D = delta_ijk(wi, wj, wk)
            delta_0 = -8.0_dp * wi * wj * wk / (D + EPS_DIV)
          end if
          chi0 = chi0 + 2.0_dp * (phi3_cm(i, j, k)**2) * delta_0
        end do
      end do
    end do

    ! rotational terms
    do k = 1, 3
      sr = 0.0_dp
      do i = 1, n
        do j = i + 1, n
          sr = sr + zeta(k, i, j)**2
        end do
      end do
      chi0 = chi0 - 16.0_dp * rotc(k) * (1.0_dp + 2.0_dp * sr)
    end do

    chi0 = chi0 / 64.0_dp
  end function compute_chi0

  !-------------------------------------------------------------------
  ! Harmonic and anharmonic ZPVE (cm^-1):
  !   ZPVE_harm = 1/2 sum_i w_i
  !   ZPVE_anh  = chi0 + 1/2 sum_i (w_i + 1/2 chi_ii) + 1/4 sum_{i<j} chi_ij
  !-------------------------------------------------------------------
  subroutine compute_zpve(omega_cm, chi, chi0, harmonic_zpve, anharmonic_zpve, correction)
    real(dp), intent(in) :: omega_cm(:), chi(:,:), chi0
    real(dp), intent(out) :: harmonic_zpve, anharmonic_zpve, correction
    integer :: n, i, j
    n = size(omega_cm)
    harmonic_zpve = 0.0_dp
    anharmonic_zpve = chi0
    do i = 1, n
      harmonic_zpve = harmonic_zpve + 0.5_dp * omega_cm(i)
      anharmonic_zpve = anharmonic_zpve + 0.5_dp * (omega_cm(i) + 0.5_dp * chi(i, i))
    end do
    do i = 1, n
      do j = i + 1, n
        anharmonic_zpve = anharmonic_zpve + 0.25_dp * chi(i, j)
      end do
    end do
    correction = anharmonic_zpve - harmonic_zpve
  end subroutine compute_zpve

  !-------------------------------------------------------------------
  ! Deperturbed state energies (nondegenerate, g = 0 limit):
  !   nu_i     = w_i + 2 chi_ii + 1/2 sum_{j/=i} chi_ij
  !   overtone = 2 w_i + 6 chi_ii + sum_{j/=i} chi_ij
  !   band_ij  = w_i + w_j + 2 chi_ii + 2 chi_jj + 2 chi_ij
  !              + 1/2 sum_{k/=i,j} (chi_ik + chi_jk)          (i < j)
  !-------------------------------------------------------------------
  subroutine build_explicit_state_energies(omega_cm, chi, nu, overtone, band)
    real(dp), intent(in) :: omega_cm(:), chi(:,:)
    real(dp), intent(out) :: nu(size(omega_cm))
    real(dp), intent(out) :: overtone(size(omega_cm))
    real(dp), intent(out) :: band(size(omega_cm), size(omega_cm))
    integer :: n, i, j, k
    real(dp) :: nu_i, ov_i, bij

    n = size(omega_cm)
    band = 0.0_dp
    do i = 1, n
      nu_i = omega_cm(i) + 2.0_dp * chi(i, i)
      ov_i = 2.0_dp * omega_cm(i) + 6.0_dp * chi(i, i)
      do j = 1, n
        if (j == i) cycle
        nu_i = nu_i + 0.5_dp * chi(i, j)
        ov_i = ov_i + chi(i, j)
      end do
      nu(i) = nu_i
      overtone(i) = ov_i
    end do

    do i = 1, n
      do j = i + 1, n
        bij = omega_cm(i) + omega_cm(j) + 2.0_dp * chi(i, i) + 2.0_dp * chi(j, j) + 2.0_dp * chi(i, j)
        do k = 1, n
          if (k == i .or. k == j) cycle
          bij = bij + 0.5_dp * (chi(i, k) + chi(j, k))
        end do
        band(i, j) = bij
      end do
    end do
  end subroutine build_explicit_state_energies

  !-------------------------------------------------------------------
  ! GVPT2: build the fundamental <-> overtone/combination interactions
  ! from the Fermi set (in the sorted order gvpt2.py iterates them),
  ! solve the polyads and update the fundamentals.
  !
  !   nu_gvpt2(n)      updated fundamentals (cm^-1)
  !   smap, nmap       all states solved in polyads and their energies
  !-------------------------------------------------------------------
  subroutine apply_gvpt2_fermi(omega_cm, nu_depert, chi, phi3_cm, fermi, nu_gvpt2, smap, nmap)
    real(dp), intent(in) :: omega_cm(:), nu_depert(:), chi(:,:), phi3_cm(:,:,:)
    logical, intent(in) :: fermi(:,:,:)
    real(dp), intent(out) :: nu_gvpt2(size(omega_cm))
    type(state_energy_t), allocatable, intent(out) :: smap(:)
    integer, intent(out) :: nmap

    integer :: n, i, j, k, ninter, m, ii, jj
    real(dp) :: nu_tmp(size(omega_cm)), overtone(size(omega_cm))
    real(dp) :: band(size(omega_cm), size(omega_cm))
    type(interaction_t), allocatable :: inters(:)

    n = size(omega_cm)
    call build_explicit_state_energies(omega_cm, chi, nu_tmp, overtone, band)

    ninter = count(fermi)
    allocate(inters(max(ninter, 1)))
    m = 0
    ! sorted(fermi_set): by k, then (i, j)
    do k = 1, n
      do i = 1, n
        do j = 1, n
          if (.not. fermi(k, i, j)) cycle
          m = m + 1
          inters(m)%left = vstate(k, 0)
          inters(m)%nu_left = nu_depert(k)
          if (i == j) then
            ! Fermi-1: 2*i ~ k
            inters(m)%right = vstate(i, i)
            inters(m)%nu_right = overtone(i)
            inters(m)%phi = phi3_cm(k, i, i)
            inters(m)%ftype = 1
          else
            ! Fermi-2: i + j ~ k
            ii = min(i, j)
            jj = max(i, j)
            inters(m)%right = vstate(ii, jj)
            inters(m)%nu_right = band(ii, jj)
            inters(m)%phi = phi3_cm(k, ii, jj)
            inters(m)%ftype = 2
          end if
        end do
      end do
    end do

    nu_gvpt2 = nu_depert
    if (m == 0) then
      allocate(smap(1))
      nmap = 0
      return
    end if

    call fermi_solver(inters, m, smap, nmap)

    ! update only the fundamentals
    do i = 1, nmap
      if (smap(i)%state%b == 0) nu_gvpt2(smap(i)%state%a) = smap(i)%energy
    end do
  end subroutine apply_gvpt2_fermi

end module gvpt2_core
