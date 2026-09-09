!=====================================================================
! gvpt2_constants
!
! Precision kind, physical constants, resonance thresholds and the
! atomic-mass table used by the Fortran translation of gvpt2.py
! (MLP-VPT2, step 2: DVPT2 + GVPT2 fundamentals).
!
! The numerical constants are exactly those hard-coded in gvpt2.py.
! The atomic masses reproduce ase.data.atomic_masses (IUPAC 2016
! standard atomic weights), which gvpt2.py obtains through ASE.
! Note that these differ slightly from step1_fortran/constants.f90,
! which uses most-abundant-isotope masses; gvpt2.py uses the ASE
! values for the rotational constants, so we do the same here.
!=====================================================================
module gvpt2_constants
  use iso_fortran_env, only: error_unit
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)

  ! unit conversion factors (same values as gvpt2.py)
  real(dp), parameter, public :: CM_PER_EH   = 219474.6313705_dp   ! cm^-1 per Hartree
  real(dp), parameter, public :: EH_PER_CM   = 1.0_dp / CM_PER_EH
  real(dp), parameter, public :: BOHR_TO_ANG = 0.529177210903_dp
  real(dp), parameter, public :: AMU_TO_ME   = 1822.888486209_dp

  ! Fermi resonance thresholds (gvpt2.py: FERMI_OMEGA_THRESH, FERMI_K_THRESH)
  real(dp), parameter, public :: FERMI_OMEGA_THRESH = 100.0_dp   ! cm^-1
  real(dp), parameter, public :: FERMI_K_THRESH     = 1.0_dp

  ! small number added to denominators (gvpt2.py: eps_div)
  real(dp), parameter, public :: EPS_DIV = 1.0e-30_dp

  ! ZPVE unit conversion factors printed by gvpt2.py
  real(dp), parameter, public :: CM_TO_KCAL = 0.002859144_dp
  real(dp), parameter, public :: CM_TO_KJ   = 0.011962656_dp

  ! ---- atomic mass table (ase.data.atomic_masses, index = atomic number) ----
  integer, parameter, public :: MAX_Z = 118

  character(len=2), parameter, public :: ase_symbol(0:MAX_Z) = [character(len=2) :: &
      'X ', 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', &
      'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', &
      'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', &
      'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', &
      'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', &
      'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', &
      'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', &
      'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', &
      'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', &
      'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', &
      'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', &
      'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og' ]

  real(dp), parameter, public :: ase_mass(0:MAX_Z) = [ &
      1.0_dp, 1.008_dp, 4.002602_dp, 6.94_dp, &
      9.0121831_dp, 10.81_dp, 12.011_dp, 14.007_dp, &
      15.999_dp, 18.998403163_dp, 20.1797_dp, 22.98976928_dp, &
      24.305_dp, 26.9815385_dp, 28.085_dp, 30.973761998_dp, &
      32.06_dp, 35.45_dp, 39.948_dp, 39.0983_dp, &
      40.078_dp, 44.955908_dp, 47.867_dp, 50.9415_dp, &
      51.9961_dp, 54.938044_dp, 55.845_dp, 58.933194_dp, &
      58.6934_dp, 63.546_dp, 65.38_dp, 69.723_dp, &
      72.63_dp, 74.921595_dp, 78.971_dp, 79.904_dp, &
      83.798_dp, 85.4678_dp, 87.62_dp, 88.90584_dp, &
      91.224_dp, 92.90637_dp, 95.95_dp, 97.90721_dp, &
      101.07_dp, 102.9055_dp, 106.42_dp, 107.8682_dp, &
      112.414_dp, 114.818_dp, 118.71_dp, 121.76_dp, &
      127.6_dp, 126.90447_dp, 131.293_dp, 132.90545196_dp, &
      137.327_dp, 138.90547_dp, 140.116_dp, 140.90766_dp, &
      144.242_dp, 144.91276_dp, 150.36_dp, 151.964_dp, &
      157.25_dp, 158.92535_dp, 162.5_dp, 164.93033_dp, &
      167.259_dp, 168.93422_dp, 173.054_dp, 174.9668_dp, &
      178.49_dp, 180.94788_dp, 183.84_dp, 186.207_dp, &
      190.23_dp, 192.217_dp, 195.084_dp, 196.966569_dp, &
      200.592_dp, 204.38_dp, 207.2_dp, 208.9804_dp, &
      208.98243_dp, 209.98715_dp, 222.01758_dp, 223.01974_dp, &
      226.02541_dp, 227.02775_dp, 232.0377_dp, 231.03588_dp, &
      238.02891_dp, 237.04817_dp, 244.06421_dp, 243.06138_dp, &
      247.07035_dp, 247.07031_dp, 251.07959_dp, 252.083_dp, &
      257.09511_dp, 258.09843_dp, 259.101_dp, 262.11_dp, &
      267.122_dp, 268.126_dp, 271.134_dp, 270.133_dp, &
      269.1338_dp, 278.156_dp, 281.165_dp, 281.166_dp, &
      285.177_dp, 286.182_dp, 289.19_dp, 289.194_dp, &
      293.204_dp, 293.208_dp, 294.214_dp ]

  public :: atomic_number, atomic_mass_amu, to_upper

contains

  !-------------------------------------------------------------------
  ! Uppercase copy of a string (ASCII only).
  !-------------------------------------------------------------------
  pure function to_upper(s) result(u)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: u
    integer :: i, c
    u = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('a') .and. c <= iachar('z')) u(i:i) = achar(c - 32)
    end do
  end function to_upper

  !-------------------------------------------------------------------
  ! Atomic number from an element symbol ("C", "Br", case-insensitive)
  ! or from a string holding the atomic number itself ("6").
  ! Returns 0 if the symbol is not recognised.
  !-------------------------------------------------------------------
  function atomic_number(symbol) result(z)
    character(len=*), intent(in) :: symbol
    integer :: z
    character(len=len(symbol)) :: s
    integer :: i, ios
    logical :: all_digits

    s = adjustl(symbol)
    z = 0
    if (len_trim(s) == 0) return

    all_digits = .true.
    do i = 1, len_trim(s)
      if (s(i:i) < '0' .or. s(i:i) > '9') then
        all_digits = .false.
        exit
      end if
    end do

    if (all_digits) then
      read(s, *, iostat=ios) z
      if (ios /= 0 .or. z < 0 .or. z > MAX_Z) z = 0
      return
    end if

    do i = 1, MAX_Z
      if (to_upper(trim(s)) == to_upper(trim(ase_symbol(i)))) then
        z = i
        return
      end if
    end do
  end function atomic_number

  !-------------------------------------------------------------------
  ! Atomic mass in amu for an element symbol (ASE table). Stops the
  ! program with a message if the symbol is unknown.
  !-------------------------------------------------------------------
  function atomic_mass_amu(symbol) result(m)
    character(len=*), intent(in) :: symbol
    real(dp) :: m
    integer :: z
    z = atomic_number(symbol)
    if (z == 0) then
      write(error_unit,'(A)') 'ERROR: unknown element symbol "' // trim(symbol) // '" in geometry file'
      stop 1
    end if
    m = ase_mass(z)
  end function atomic_mass_amu

end module gvpt2_constants
