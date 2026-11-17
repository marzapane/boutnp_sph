! =============================================================================
! SHT_3D_POSTPROCESS
!
! This utility performs a Spherical Harmonic Transform (SHT) on a 3D vector
! velocity field using the SHTOOLS library.
!
! Features:
! - Modular design separating I/O and SHT processing.
! - Computes the total kinetic energy spectrum, and separate m=0 / m!=0 spectra.
! - Properly handles both polar and Cartesian input velocity fields.
!
! Compilation (example):
! gfortran -O3 sht_3d_postprocess.f90 -o sht_3d_utility -lshtools -lfftw3 -lm
! gfortran -O3 src/sht_3d_postprocess.f90 -o sht_3d_utility -I/usr/local/include -m64 -fPIC -O3 -std=gnu -ffast-math -L/usr/local/lib -lSHTOOLS -L/usr/local -lfftw3 -lm -llapack -lblas
! =============================================================================

module data_io_mod
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    implicit none

contains

    subroutine read_simulation_data(file_path, n1, n2, n3, sim_time, vt, vr, vp, coord_type)
        character(len=*), intent(in) :: file_path
        integer, intent(out) :: n1, n2, n3
        real(kind=real64), intent(out) :: sim_time
        real(kind=real64), allocatable, intent(out) :: vt(:,:,:), vr(:,:,:), vp(:,:,:)
        character(len=*), intent(in) :: coord_type

        integer :: in_unit = 10, io_status, dummy_int, k
        real(kind=real64) :: epsil, re
        real(kind=real64), allocatable :: dummy_slice(:,:)

        write(*,*) "-> Reading data from: ", trim(file_path)
        open(unit=in_unit, file=trim(file_path), form='unformatted', status='old', action='read', iostat=io_status)
        if (io_status /= 0) then
            write(error_unit, '("Error: Could not open ", A)') trim(file_path)
            stop 1
        end if

        read(in_unit) dummy_int
        read(in_unit) n1, n2, n3
        read(in_unit) epsil, epsil, re, sim_time

        allocate(dummy_slice(n1, n2))
        allocate(vt(n1, n2, n3), vr(n1, n2, n3), vp(n1, n2, n3))

        ! Read fields (dummy_slice skips field_temp; pr is skipped at end of record)
        read(in_unit) (dummy_slice, k = 1, n3), vt, vr, vp
        deallocate(dummy_slice)
        close(in_unit)

        call convert_to_spherical(vt, vr, vp, n1, n2, n3, coord_type)

    end subroutine read_simulation_data

    subroutine convert_to_spherical(vt, vr, vp, n1, n2, n3, coord_type)
        real(kind=real64), intent(inout) :: vt(:,:,:), vr(:,:,:), vp(:,:,:)
        integer, intent(in) :: n1, n2, n3
        character(len=*), intent(in) :: coord_type

        real(kind=real64), parameter :: pi = 2.0_real64 * asin(1.0_real64)
        real(kind=real64) :: phi, cp, sp, theta, ct, st, vx, vy, vz
        integer :: i, j, k

        if (trim(coord_type) == '-p') then
            write(*,*) "-> Using spherical velocity components directly (vt, vr, vp)."
        else if (trim(coord_type) == '-c') then
            write(*,*) "-> Converting Cartesian velocity components to spherical..."
            do k = 1, n3
                ! Colatitude (polar angle phi) varies along k
                phi = real(k-1, real64) * pi / real(n3-1, real64)
                cp = cos(phi)
                sp = sin(phi)
                do j = 1, n2
                    do i = 1, n1
                        ! Longitude (azimuthal angle theta) varies along i
                        theta = real(i-1, real64) * 2.0_real64 * pi / real(n1-1, real64)
                        ct = cos(theta)
                        st = sin(theta)

                        vx = vt(i, j, k)
                        vy = vr(i, j, k)
                        vz = vp(i, j, k)

                        vt(i, j, k) = -vx*st    + vy*ct
                        vr(i, j, k) =  vx*sp*ct + vy*sp*st + vz*cp
                        vp(i, j, k) =  vx*cp*ct + vy*cp*st - vz*sp
                    end do
                end do
            end do
        else
            write(error_unit, '("Error: Unknown coordinate type option. Must be ''-p'' or ''-c''.")')
            stop 1
        end if
    end subroutine convert_to_spherical

end module data_io_mod


module binner_mod
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    real(kind=real64), parameter, private :: pi = 2.0_real64 * asin(1.0_real64)

    type :: spectral_binner_t
        integer :: lmax = 0, n_shells = 0, n_bins = 0
        real(kind=real64) :: r_min = 0.0_real64, r_max = 0.0_real64, dr = 0.0_real64
        real(kind=real64) :: log_K_min = 0.0_real64, log_K_max = 0.0_real64
        real(kind=real64) :: delta_log = 0.0_real64, inv_delta_log = 0.0_real64
        real(kind=real64), allocatable :: K_bins(:)
        real(kind=real64), allocatable :: val_tot(:), val_1(:), val_2(:)
        real(kind=real64), allocatable :: count_tot(:), count_1(:), count_2(:)
    contains
        procedure :: init => binner_init
        procedure :: accumulate_volume_energy => binner_accumulate_volume_energy
        procedure :: accumulate_horizontal_energy => binner_accumulate_horizontal_energy
        procedure :: normalize => binner_normalize
        procedure :: compute_flux => binner_compute_flux
        procedure :: destroy => binner_destroy
    end type spectral_binner_t

contains

    subroutine binner_init(self, lmax, n_shells, r_min, r_max, n_bins)
        class(spectral_binner_t), intent(inout) :: self
        integer, intent(in) :: lmax, n_shells, n_bins
        real(kind=real64), intent(in) :: r_min, r_max

        integer :: j

        self%lmax     = lmax
        self%n_shells = n_shells
        self%r_min    = r_min
        self%r_max    = r_max
        self%n_bins   = n_bins

        self%dr = (r_max - r_min) / real(n_shells - 1, real64)

        self%log_K_min = log(min(sqrt(2.0_real64)/r_max, 2.0_real64*pi/(r_max-r_min)))
        self%log_K_max = 0.5_real64 * log(real(lmax*(lmax+1), real64)/r_min**2 &
            + (2.0_real64*sqrt(real(n_shells*(n_shells+1), real64))*2.0_real64*pi/(r_max-r_min))**2)
        self%delta_log = (self%log_K_max - self%log_K_min) / real(n_bins, real64)
        self%inv_delta_log = real(n_bins, real64) / (self%log_K_max - self%log_K_min)

        allocate(self%K_bins(0:n_bins))
        allocate(self%val_tot(0:n_bins), self%val_1(0:n_bins), self%val_2(0:n_bins))
        allocate(self%count_tot(0:n_bins), self%count_1(0:n_bins), self%count_2(0:n_bins))

        self%K_bins(0) = 0.0_real64
        do j = 1, n_bins
            self%K_bins(j) = exp(self%log_K_min + (real(j, real64) - 0.5_real64) * self%delta_log)
        end do

        self%val_tot   = 0.0_real64
        self%val_1     = 0.0_real64
        self%val_2     = 0.0_real64
        self%count_tot = 0.0_real64
        self%count_1   = 0.0_real64
        self%count_2   = 0.0_real64

    end subroutine binner_init

    subroutine binner_accumulate_volume_energy(self, se_diff, k, comp)
        class(spectral_binner_t), intent(inout) :: self
        real(kind=real64), intent(in) :: se_diff(0:self%lmax, self%n_shells)
        integer, intent(in) :: k, comp

        logical :: is_2d_mode
        real(kind=real64) :: k_rad, k_ang, K_tot, r_j, E_contrib
        real(kind=real64) :: real_idx, w_lo, w_hi
        integer :: j, l, idx

        is_2d_mode = (k == 0 .and. comp /= 2)

        k_rad = 2.0_real64 * sqrt(real(k*(k+1), real64)) * 2.0_real64 * pi / (self%r_max - self%r_min)

        do j = 1, self%n_shells
            r_j = self%r_min + real(j - 1, real64) * self%dr
            do l = 0, self%lmax
                E_contrib = 2.0_real64 * pi * se_diff(l, j) * r_j**2 * self%dr

                k_ang = sqrt(real(l*(l+1), real64)) / r_j
                K_tot = sqrt(k_ang**2 + k_rad**2)

                ! Guard against log(0) for the DC / zero-wavenumber mode (k=0, l=0)
                if (K_tot <= 0.0_real64) then
                    self%val_tot(0)   = self%val_tot(0)   + E_contrib
                    self%count_tot(0) = self%count_tot(0) + 1.0_real64
                    if (is_2d_mode) then
                        self%val_1(0)   = self%val_1(0)   + E_contrib
                        self%count_1(0) = self%count_1(0) + 1.0_real64
                    else
                        self%val_2(0)   = self%val_2(0)   + E_contrib
                        self%count_2(0) = self%count_2(0) + 1.0_real64
                    end if
                    cycle
                end if

                ! Direct O(1) index relative to bin centres (offset by +0.5)
                real_idx = (log(K_tot) - self%log_K_min) * self%inv_delta_log + 0.5_real64
                idx      = int(real_idx)

                ! CIC: interpolate in logarithmic scale between the two nearest bin centres
                if (idx < 1) then
                    self%val_tot(1)   = self%val_tot(1)   + E_contrib
                    self%count_tot(1) = self%count_tot(1) + 1.0_real64
                    if (is_2d_mode) then
                        self%val_1(1)   = self%val_1(1)   + E_contrib
                        self%count_1(1) = self%count_1(1) + 1.0_real64
                    else
                        self%val_2(1)   = self%val_2(1)   + E_contrib
                        self%count_2(1) = self%count_2(1) + 1.0_real64
                    end if
                else if (idx >= self%n_bins) then
                    self%val_tot(self%n_bins)   = self%val_tot(self%n_bins)   + E_contrib
                    self%count_tot(self%n_bins) = self%count_tot(self%n_bins) + 1.0_real64
                    if (is_2d_mode) then
                        self%val_1(self%n_bins)   = self%val_1(self%n_bins)   + E_contrib
                        self%count_1(self%n_bins) = self%count_1(self%n_bins) + 1.0_real64
                    else
                        self%val_2(self%n_bins)   = self%val_2(self%n_bins)   + E_contrib
                        self%count_2(self%n_bins) = self%count_2(self%n_bins) + 1.0_real64
                    end if
                else
                    ! CIC Interpolation
                    w_hi = real_idx - real(idx, real64)
                    w_lo = 1.0_real64 - w_hi

                    self%val_tot(idx)     = self%val_tot(idx)     + E_contrib * w_lo
                    self%val_tot(idx+1)   = self%val_tot(idx+1)   + E_contrib * w_hi
                    self%count_tot(idx)   = self%count_tot(idx)   + w_lo
                    self%count_tot(idx+1) = self%count_tot(idx+1) + w_hi

                    if (is_2d_mode) then
                        self%val_1(idx)     = self%val_1(idx)     + E_contrib * w_lo
                        self%val_1(idx+1)   = self%val_1(idx+1)   + E_contrib * w_hi
                        self%count_1(idx)   = self%count_1(idx)   + w_lo
                        self%count_1(idx+1) = self%count_1(idx+1) + w_hi
                    else
                        self%val_2(idx)     = self%val_2(idx)     + E_contrib * w_lo
                        self%val_2(idx+1)   = self%val_2(idx+1)   + E_contrib * w_hi
                        self%count_2(idx)   = self%count_2(idx)   + w_lo
                        self%count_2(idx+1) = self%count_2(idx+1) + w_hi
                    end if
                end if
            end do
        end do

    end subroutine binner_accumulate_volume_energy

    subroutine binner_accumulate_horizontal_energy(self, spec_m0, spec_rest)
        class(spectral_binner_t), intent(inout) :: self
        real(kind=real64), intent(in) :: spec_m0(0:self%lmax, self%n_shells)
        real(kind=real64), intent(in) :: spec_rest(0:self%lmax, self%n_shells)

        real(kind=real64) :: r_j, K_tot, E_m0, E_rest, E_tot
        real(kind=real64) :: real_idx, w_lo, w_hi
        integer :: j, l, idx

        do j = 1, self%n_shells
            r_j = self%r_min + real(j - 1, real64) * self%dr
            do l = 0, self%lmax
                E_m0   = 2.0_real64 * pi * spec_m0(l, j)   * r_j**2 * self%dr
                E_rest = 2.0_real64 * pi * spec_rest(l, j) * r_j**2 * self%dr
                E_tot  = E_m0 + E_rest

                if (E_m0 == 0.0_real64 .and. E_rest == 0.0_real64) cycle

                K_tot = sqrt(real(l * (l + 1), real64)) / r_j

                ! Guard against log(0) for the DC / zero-wavenumber mode (l=0)
                if (K_tot <= 0.0_real64) then
                    self%val_tot(0)   = self%val_tot(0)   + E_tot
                    self%val_1(0)     = self%val_1(0)     + E_m0
                    self%val_2(0)     = self%val_2(0)     + E_rest
                    self%count_tot(0) = self%count_tot(0) + 1.0_real64
                    self%count_1(0)   = self%count_1(0)   + 1.0_real64
                    self%count_2(0)   = self%count_2(0)   + 1.0_real64
                    cycle
                end if

                real_idx = (log(K_tot) - self%log_K_min) * self%inv_delta_log + 0.5_real64
                idx      = int(real_idx)

                if (idx < 1) then
                    self%val_tot(1)   = self%val_tot(1)   + E_tot
                    self%val_1(1)     = self%val_1(1)     + E_m0
                    self%val_2(1)     = self%val_2(1)     + E_rest
                    self%count_tot(1) = self%count_tot(1) + 1.0_real64
                    self%count_1(1)   = self%count_1(1)   + 1.0_real64
                    self%count_2(1)   = self%count_2(1)   + 1.0_real64
                else if (idx >= self%n_bins) then
                    self%val_tot(self%n_bins)   = self%val_tot(self%n_bins)   + E_tot
                    self%val_1(self%n_bins)     = self%val_1(self%n_bins)     + E_m0
                    self%val_2(self%n_bins)     = self%val_2(self%n_bins)     + E_rest
                    self%count_tot(self%n_bins) = self%count_tot(self%n_bins) + 1.0_real64
                    self%count_1(self%n_bins)   = self%count_1(self%n_bins)   + 1.0_real64
                    self%count_2(self%n_bins)   = self%count_2(self%n_bins)   + 1.0_real64
                else
                    ! CIC Interpolation
                    w_hi         = real_idx - real(idx, real64)
                    w_lo         = 1.0_real64 - w_hi

                    self%val_tot(idx)     = self%val_tot(idx)     + E_tot  * w_lo
                    self%val_tot(idx+1)   = self%val_tot(idx+1)   + E_tot  * w_hi
                    self%val_1(idx)       = self%val_1(idx)       + E_m0   * w_lo
                    self%val_1(idx+1)     = self%val_1(idx+1)     + E_m0   * w_hi
                    self%val_2(idx)       = self%val_2(idx)       + E_rest * w_lo
                    self%val_2(idx+1)     = self%val_2(idx+1)     + E_rest * w_hi

                    self%count_tot(idx)   = self%count_tot(idx)   + w_lo
                    self%count_tot(idx+1) = self%count_tot(idx+1) + w_hi
                    self%count_1(idx)     = self%count_1(idx)     + w_lo
                    self%count_1(idx+1)   = self%count_1(idx+1)   + w_hi
                    self%count_2(idx)     = self%count_2(idx)     + w_lo
                    self%count_2(idx+1)   = self%count_2(idx+1)   + w_hi
                end if
            end do
        end do

    end subroutine binner_accumulate_horizontal_energy

    subroutine binner_normalize(self)
        class(spectral_binner_t), intent(inout) :: self

        real(kind=real64) :: w_lo
        integer :: j

        w_lo = exp(0.5_real64 * self%delta_log) - exp(-0.5_real64 * self%delta_log)
        do j = 1, self%n_bins
            self%val_tot(j)   = self%val_tot(j)   / (self%K_bins(j) * w_lo)
            self%val_1(j)     = self%val_1(j)     / (self%K_bins(j) * w_lo)
            self%val_2(j)     = self%val_2(j)     / (self%K_bins(j) * w_lo)
            self%count_tot(j) = self%count_tot(j) / (self%K_bins(j) * w_lo)
            self%count_1(j)   = self%count_1(j)   / (self%K_bins(j) * w_lo)
            self%count_2(j)   = self%count_2(j)   / (self%K_bins(j) * w_lo)
        end do

    end subroutine binner_normalize

    subroutine binner_compute_flux(self, flux_tot, flux_1, flux_2)
        class(spectral_binner_t), intent(in) :: self
        real(kind=real64), intent(out) :: flux_tot(0:self%n_bins)
        real(kind=real64), intent(out) :: flux_1(0:self%n_bins)
        real(kind=real64), intent(out) :: flux_2(0:self%n_bins)

        real(kind=real64) :: w_lo, dK, cum_tot, cum_1, cum_2
        integer :: i

        w_lo = exp(0.5_real64 * self%delta_log) - exp(-0.5_real64 * self%delta_log)

        ! DC mode (i = 0) was not divided by dK in normalize()
        cum_tot = -self%val_tot(0)
        cum_1   = -self%val_1(0)
        cum_2   = -self%val_2(0)

        flux_tot(0) = cum_tot
        flux_1(0)   = cum_1
        flux_2(0)   = cum_2

        do i = 1, self%n_bins
            dK = self%K_bins(i) * w_lo
            ! Evaluate cumulative flux at bin centre K_bins(i) (midpoint rule),
            ! matching the wavenumber where T(K_i) is defined.
            flux_tot(i) = cum_tot - 0.5_real64 * self%val_tot(i) * dK
            flux_1(i)   = cum_1   - 0.5_real64 * self%val_1(i)   * dK
            flux_2(i)   = cum_2   - 0.5_real64 * self%val_2(i)   * dK

            ! Advance cumulative sum to the upper edge of bin i
            cum_tot = cum_tot - self%val_tot(i) * dK
            cum_1   = cum_1   - self%val_1(i)   * dK
            cum_2   = cum_2   - self%val_2(i)   * dK
        end do
    end subroutine binner_compute_flux

    subroutine binner_destroy(self)
        class(spectral_binner_t), intent(inout) :: self

        if (allocated(self%K_bins))    deallocate(self%K_bins)
        if (allocated(self%val_tot))   deallocate(self%val_tot)
        if (allocated(self%val_1))     deallocate(self%val_1)
        if (allocated(self%val_2))     deallocate(self%val_2)
        if (allocated(self%count_tot)) deallocate(self%count_tot)
        if (allocated(self%count_1))   deallocate(self%count_1)
        if (allocated(self%count_2))   deallocate(self%count_2)

    end subroutine binner_destroy

end module binner_mod


module advection_mod
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    real(kind=real64), parameter, private :: pi = 2.0_real64 * asin(1.0_real64)

contains

    ! Computes one component of the nonlinear advection acceleration vector on a single
    ! radial shell j:
    !   N = -(u . grad)u in spherical coordinates (r, theta, phi)
    ! including metric curvature terms.
    !   comp = 1: N_theta (azimuthal)
    !   comp = 2: N_r     (radial)
    !   comp = 3: N_phi   (meridional / polar)
    subroutine compute_advection_shell(vt, vr, vp, n1, n2, n3, rint, dr, j, comp, N_shell)
        integer,           intent(in)  :: n1, n2, n3, j, comp
        real(kind=real64), intent(in)  :: vt(n1, n2, n3), vr(n1, n2, n3), vp(n1, n2, n3)
        real(kind=real64), intent(in)  :: rint, dr
        real(kind=real64), intent(out) :: N_shell(:,:)

        real(kind=real64) :: dtheta, dphi, inv_2dtheta, inv_2dphi, inv_2dr
        real(kind=real64) :: r_j, inv_r, phi_k, sin_phi, cos_phi, cot_phi, inv_r_sin_phi
        real(kind=real64) :: df_dr, df_dphi, df_dtheta, u_dot_grad_f
        integer :: i, k, i_prev, i_next

        dtheta = 2.0_real64 * pi / real(n1 - 1, real64)
        dphi   = pi / real(n3 - 1, real64)

        inv_2dtheta = 1.0_real64 / (2.0_real64 * dtheta)
        inv_2dphi   = 1.0_real64 / (2.0_real64 * dphi)
        inv_2dr     = 1.0_real64 / (2.0_real64 * dr)

        r_j   = rint + real(j - 1, real64) * dr
        inv_r = 1.0_real64 / r_j

        N_shell = 0.0_real64

        !$OMP PARALLEL DO PRIVATE(k, i, i_prev, i_next, phi_k, sin_phi, cos_phi, cot_phi, &
        !$OMP& inv_r_sin_phi, df_dr, df_dphi, df_dtheta, u_dot_grad_f) SCHEDULE(STATIC)
        do k = 2, n3 - 1
            phi_k = real(k - 1, real64) * dphi
            sin_phi = sin(phi_k)
            cos_phi = cos(phi_k)
            cot_phi = cos_phi / sin_phi
            inv_r_sin_phi = 1.0_real64 / (r_j * sin_phi)

            do i = 1, n1 - 1
                i_prev = merge(n1 - 1, i - 1, i == 1)
                i_next = merge(1,      i + 1, i == n1 - 1)

                select case (comp)
                case (1)
                    ! Component 1: theta / azimuthal (vt)
                    if (j == 1) then
                        df_dr = (-3.0_real64 * vt(i, 1, k) + 4.0_real64 * vt(i, 2, k) - vt(i, 3, k)) * inv_2dr
                    else if (j == n2) then
                        df_dr = (3.0_real64 * vt(i, n2, k) - 4.0_real64 * vt(i, n2 - 1, k) + vt(i, n2 - 2, k)) * inv_2dr
                    else
                        df_dr = (vt(i, j + 1, k) - vt(i, j - 1, k)) * inv_2dr
                    end if
                    df_dphi   = (vt(i, j, k + 1) - vt(i, j, k - 1)) * inv_2dphi
                    df_dtheta = (vt(i_next, j, k) - vt(i_prev, j, k)) * inv_2dtheta

                    u_dot_grad_f = vr(i, j, k) * df_dr + &
                                   (vp(i, j, k) * inv_r) * df_dphi + &
                                   (vt(i, j, k) * inv_r_sin_phi) * df_dtheta

                    N_shell(i, k) = -( u_dot_grad_f + &
                                      (vt(i, j, k) * vr(i, j, k) * inv_r) + &
                                      (vt(i, j, k) * vp(i, j, k) * cot_phi * inv_r) )

                case (2)
                    ! Component 2: radial (vr)
                    if (j == 1) then
                        df_dr = (-3.0_real64 * vr(i, 1, k) + 4.0_real64 * vr(i, 2, k) - vr(i, 3, k)) * inv_2dr
                    else if (j == n2) then
                        df_dr = (3.0_real64 * vr(i, n2, k) - 4.0_real64 * vr(i, n2 - 1, k) + vr(i, n2 - 2, k)) * inv_2dr
                    else
                        df_dr = (vr(i, j + 1, k) - vr(i, j - 1, k)) * inv_2dr
                    end if
                    df_dphi   = (vr(i, j, k + 1) - vr(i, j, k - 1)) * inv_2dphi
                    df_dtheta = (vr(i_next, j, k) - vr(i_prev, j, k)) * inv_2dtheta

                    u_dot_grad_f = vr(i, j, k) * df_dr + &
                                   (vp(i, j, k) * inv_r) * df_dphi + &
                                   (vt(i, j, k) * inv_r_sin_phi) * df_dtheta

                    N_shell(i, k) = -( u_dot_grad_f - &
                                      (vt(i, j, k)**2 + vp(i, j, k)**2) * inv_r )

                case (3)
                    ! Component 3: phi / meridional (vp)
                    if (j == 1) then
                        df_dr = (-3.0_real64 * vp(i, 1, k) + 4.0_real64 * vp(i, 2, k) - vp(i, 3, k)) * inv_2dr
                    else if (j == n2) then
                        df_dr = (3.0_real64 * vp(i, n2, k) - 4.0_real64 * vp(i, n2 - 1, k) + vp(i, n2 - 2, k)) * inv_2dr
                    else
                        df_dr = (vp(i, j + 1, k) - vp(i, j - 1, k)) * inv_2dr
                    end if
                    df_dphi   = (vp(i, j, k + 1) - vp(i, j, k - 1)) * inv_2dphi
                    df_dtheta = (vp(i_next, j, k) - vp(i_prev, j, k)) * inv_2dtheta

                    u_dot_grad_f = vr(i, j, k) * df_dr + &
                                   (vp(i, j, k) * inv_r) * df_dphi + &
                                   (vt(i, j, k) * inv_r_sin_phi) * df_dtheta

                    N_shell(i, k) = -( u_dot_grad_f + &
                                      (vp(i, j, k) * vr(i, j, k) * inv_r) - &
                                      (vt(i, j, k)**2 * cot_phi * inv_r) )
                end select
            end do
            N_shell(n1, k) = N_shell(1, k)
        end do
        !$OMP END PARALLEL DO

        ! Coordinate singularity guard at poles (k = 1, n3):
        ! A smooth horizontal vector at the pole has theta-dependent components (vt, vp),
        ! and 1/sin(phi) diverges. Finite placeholders here prevent NaN/Inf
        ! in SHExpandDH while having zero effect on results because:
        !   1. the physical volume element dOmega = sin(phi) dtheta dphi vanishes at the poles,
        !   2. the Driscoll-Healy quadrature weight at k = 1 (phi = 0) is identically zero, and
        !   3. the south pole row k = n3 lies outside the 1..n_lat transform slice.
        if (comp == 2) then
            N_shell(:, 1)  = N_shell(:, 2)
            N_shell(:, n3) = N_shell(:, n3 - 1)
        else
            N_shell(:, 1)  = 0.0_real64
            N_shell(:, n3) = 0.0_real64
        end if

    end subroutine compute_advection_shell

    ! Thin wrapper computing one full 3D component of N = -(u . grad)u across all shells.
    subroutine compute_advection_component(vt, vr, vp, n1, n2, n3, rint, dr, comp, N_comp)
        integer,           intent(in)  :: n1, n2, n3, comp
        real(kind=real64), intent(in)  :: vt(n1, n2, n3), vr(n1, n2, n3), vp(n1, n2, n3)
        real(kind=real64), intent(in)  :: rint, dr
        real(kind=real64), intent(out) :: N_comp(n1, n2, n3)

        integer :: j

        do j = 1, n2
            call compute_advection_shell(vt, vr, vp, n1, n2, n3, rint, dr, j, comp, N_comp(:, j, :))
        end do
    end subroutine compute_advection_component

    ! Computes the physical-space nonlinear transfer integral contribution from a single
    ! radial shell j: sum_{k, i<n1} u_shell(i,k) * N_shell(i,k) * r_j^2 * dr * d_theta * d_phi * sin(phi_k).
    pure real(kind=real64) function shell_transfer(u_shell, N_shell, n1, n3, r_j, dr, d_theta, d_phi) result(T_shell)
        integer,           intent(in) :: n1, n3
        real(kind=real64), intent(in) :: u_shell(:,:), N_shell(:,:)
        real(kind=real64), intent(in) :: r_j, dr, d_theta, d_phi

        real(kind=real64) :: phi, sp
        integer :: i, k

        T_shell = 0.0_real64
        do k = 1, n3
            phi = real(k - 1, real64) * d_phi
            sp  = sin(phi)
            do i = 1, n1 - 1
                T_shell = T_shell + u_shell(i, k) * N_shell(i, k) * r_j**2 * dr * d_theta * d_phi * sp
            end do
        end do
    end function shell_transfer

end module advection_mod


module sht_processing_mod
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    use SHTOOLS, only: SHExpandDH
    use binner_mod, only: spectral_binner_t
    implicit none

    real(kind=real64), parameter, private :: pi = 2.0_real64 * asin(1.0_real64)

contains

    ! 1. Computes the transform and returns the spherical harmonics coefficients
    subroutine get_sht_coeffs(grid, n_lat, lmax_calc, coeffs)
        integer, intent(in) :: n_lat, lmax_calc
        integer             :: lmax
        real(kind=real64), intent(in) :: grid(n_lat, 2*n_lat)
        real(kind=real64), intent(out) :: coeffs(2, lmax_calc+1, lmax_calc+1)

        ! Perform SHT (Grid -> Coefficients)
        ! sampling=2: grid is equally spaced (n x 2n)
        call SHExpandDH(grid, n_lat, coeffs, lmax, lmax_calc=lmax_calc, sampling=2)
    end subroutine get_sht_coeffs

    ! 2. Packs the 2D spherical harmonic coefficients for 0 <= m <= l into a contiguous 1D array
    !    where mode (l, m) is stored at index lm = (l * (l + 1)) / 2 + m + 1.
    !    Entries with m > l in the SHTOOLS square output are mathematically zero and are omitted,
    !    cutting coefficient storage by ~50%.
    subroutine pack_sht_coeffs(coeffs_shell, lmax, n_lm, coeffs_packed)
        integer,           intent(in)  :: lmax, n_lm
        real(kind=real64), intent(in)  :: coeffs_shell(2, 0:lmax, 0:lmax)
        real(kind=real64), intent(out) :: coeffs_packed(2, n_lm)

        integer :: l, m, lm

        do l = 0, lmax
            lm = (l * (l + 1)) / 2 + 1
            do m = 0, l
                coeffs_packed(1, lm + m) = coeffs_shell(1, l, m)
                coeffs_packed(2, lm + m) = coeffs_shell(2, l, m)
            end do
        end do
    end subroutine pack_sht_coeffs

    ! Wrapper for the SHT
    subroutine compute_sht_coeffs(field_slice, n_lat, n_lon, lmax, n_lm, &
                                  coeffs_packed)
        integer,           intent(in)    :: n_lat, n_lon, lmax, n_lm
        real(kind=real64), intent(in)    :: field_slice(:,:)
        real(kind=real64), intent(inout) :: coeffs_packed(2, n_lm)

        real(kind=real64), allocatable :: grid_2d(:,:)
        real(kind=real64), allocatable :: coeffs_shell(:,:,:)    
        
        allocate(grid_2d(n_lat, n_lon))
        allocate(coeffs_shell(2, 0:lmax, 0:lmax))

        grid_2d(:, :) = transpose(field_slice(1:n_lon, 1:n_lat))
        call get_sht_coeffs(grid_2d, n_lat, lmax, coeffs_shell)
        call pack_sht_coeffs(coeffs_shell, lmax, n_lm, coeffs_packed)
        
        deallocate(grid_2d, coeffs_shell)

    end subroutine compute_sht_coeffs

    ! 3. Computes the horizontal spectra and energy transfer for ONE velocity component 
    !    across all shells, accumulating on the fly into binner using 
    !    triangularly packed (l, m) coefficients.
    subroutine compute_horizontal_spectrum_comp(coeffs_u, lmax, n_lm, n_shells, binner_H, &
                                                coeffs_N, binner_TH)
        integer,                 intent(in)    :: lmax, n_lm, n_shells
        real(kind=real64),       intent(in)    :: coeffs_u(2, n_lm, n_shells)
        type(spectral_binner_t), intent(inout) :: binner_H
        real(kind=real64),       intent(in),    optional :: coeffs_N(2, n_lm, n_shells)
        type(spectral_binner_t), intent(inout), optional :: binner_TH

        real(kind=real64), allocatable :: spec_m0(:,:), spec_rest(:,:)
        real(kind=real64), allocatable :: trans_m0(:,:), trans_rest(:,:)
        integer :: j, l, m, lm_base

        allocate(spec_m0(0:lmax, n_shells), spec_rest(0:lmax, n_shells))
        spec_m0   = 0.0_real64
        spec_rest = 0.0_real64

        do j = 1, n_shells
            do l = 0, lmax
                lm_base = (l * (l + 1)) / 2 + 1
                ! The power at m=0 is the square of the cosine coefficient (1, lm_base, j)
                ! (assuming default 4-pi normalized spherical harmonics)
                spec_m0(l, j) = coeffs_u(1, lm_base, j)**2
                do m = 1, l
                    ! The rest is the total minus the m=0 component
                    spec_rest(l, j) = spec_rest(l, j) + &
                        coeffs_u(1, lm_base + m, j)**2 + coeffs_u(2, lm_base + m, j)**2
                end do
            end do
        end do

        call binner_H%accumulate_horizontal_energy(spec_m0, spec_rest)
        deallocate(spec_m0, spec_rest)

        if (present(coeffs_N) .and. present(binner_TH)) then
            allocate(trans_m0(0:lmax, n_shells), trans_rest(0:lmax, n_shells))
            trans_m0   = 0.0_real64
            trans_rest = 0.0_real64

            do j = 1, n_shells
                do l = 0, lmax
                    lm_base = (l * (l + 1)) / 2 + 1
                    ! Factor of 2.0 accounts for Plancherel theorem inner product (u . N)
                    trans_m0(l, j) = 2.0_real64 * coeffs_u(1, lm_base, j) * coeffs_N(1, lm_base, j)
                    do m = 1, l
                        trans_rest(l, j) = trans_rest(l, j) + 2.0_real64 * ( &
                            coeffs_u(1, lm_base + m, j) * coeffs_N(1, lm_base + m, j) + &
                            coeffs_u(2, lm_base + m, j) * coeffs_N(2, lm_base + m, j) )
                    end do
                end do
            end do

            call binner_TH%accumulate_horizontal_energy(trans_m0, trans_rest)
            deallocate(trans_m0, trans_rest)
        end if
    end subroutine compute_horizontal_spectrum_comp

    ! =========================================================================
    ! VOLUME SPECTRUM ROUTINES
    ! =========================================================================

    ! 4. Applies radial coarse-graining filter to packed SHT coefficients for ONE component.
    !
    !    Filter scale behavior controlled by ell:
    !      - ell < 0.0 : Infinite filter scale (k = 0), radial average over all shells.
    !                    Represents the depth-independent (radially uniform) mean flow.
    !      - ell == 0.0: Unfiltered scale (k = n_shells, k_max), raw coefficients.
    !      - ell > 0.0 : Gaussian low-pass filter with width ell along the radial coordinate.
    subroutine compute_filtered_sht_coeffs(coeffs_comp, n_lm, n_shells, ell, coeffs_filt)
        integer,           intent(in)  :: n_lm, n_shells
        real(kind=real64), intent(in)  :: coeffs_comp(2, n_lm, n_shells)
        real(kind=real64), intent(in)  :: ell
        real(kind=real64), intent(out) :: coeffs_filt(2, n_lm, n_shells)

        real(kind=real64) :: sigma, inv_2sig2, norm
        integer :: j_max, j_lo, j_hi, i, j, k_sh, lm
        real(kind=real64), allocatable :: w(:)

        if (ell < 0.0_real64) then
            ! Case 1: k = 0 infinite filter scale (radial average over all shells).
            coeffs_filt(:, :, 1) = 0.0_real64
            do j = 1, n_shells
                do lm = 1, n_lm
                    coeffs_filt(1, lm, 1) = coeffs_filt(1, lm, 1) + coeffs_comp(1, lm, j)
                    coeffs_filt(2, lm, 1) = coeffs_filt(2, lm, 1) + coeffs_comp(2, lm, j)
                end do
            end do
            coeffs_filt(:, :, 1) = coeffs_filt(:, :, 1) / real(n_shells, real64)
            do j = 2, n_shells
                coeffs_filt(:, :, j) = coeffs_filt(:, :, 1)
            end do

        else if (abs(ell) <= 1.0e-8_real64) then
            ! Case 2: k = n_shells unfiltered scale (ell = 0, highest radial wavenumber).
            coeffs_filt(:, :, :) = coeffs_comp(:, :, :)

        else
            ! Case 3: Finite Gaussian filter scale (1 <= k <= n_shells - 1).
            ! Gaussian low-pass kernel: sigma = ell / sqrt(12), truncated at 3*sigma.
            sigma = ell / sqrt(12.0_real64)
            j_max = min(3 * ceiling(sigma), n_shells - 1)
            inv_2sig2 = 1.0_real64 / (2.0_real64 * sigma * sigma)

            ! Precompute symmetric Gaussian kernel weights w(0:2*j_max):
            !   w(k) = exp( -(k - j_max)^2 / (2 sigma^2) )
            allocate(w(0:2*j_max))
            do j = 0, 2*j_max
                w(j) = exp(-real((j - j_max)**2, real64) * inv_2sig2)
            end do

            ! Convolve each 1D radial pencil across shells along contiguous packed (l, m) modes.
            ! Boundaries are handled by truncating the kernel to [j_lo, j_hi] and renormalizing by norm.
            !$OMP PARALLEL DO PRIVATE(i, j_lo, j_hi, norm, j, k_sh, lm) SCHEDULE(STATIC)
            do i = 1, n_shells
                ! Clamp kernel window to stay within valid shell range [1, n_shells]
                j_lo = max(0, j_max - (i - 1))
                j_hi = min(2*j_max, j_max + (n_shells - i))

                norm = 0.0_real64
                do j = j_lo, j_hi
                    norm = norm + w(j)
                end do

                coeffs_filt(:, :, i) = 0.0_real64
                do j = j_lo, j_hi
                    k_sh = i + (j - j_max)
                    do lm = 1, n_lm
                        coeffs_filt(1, lm, i) = coeffs_filt(1, lm, i) + w(j) * coeffs_comp(1, lm, k_sh)
                        coeffs_filt(2, lm, i) = coeffs_filt(2, lm, i) + w(j) * coeffs_comp(2, lm, k_sh)
                    end do
                end do
                coeffs_filt(:, :, i) = coeffs_filt(:, :, i) / norm
            end do
            !$OMP END PARALLEL DO

            deallocate(w)
        end if
    end subroutine compute_filtered_sht_coeffs

    ! 5. Computes filtered spectral energy density fse(0:lmax, 1:n_shells) for ONE component
    !    from packed SHT coefficients, and optionally computes the scale-to-scale energy
    !    transfer fst(0:lmax, 1:n_shells) when coeffs_N_filt and fst are provided:
    !      fse(l, i) = sum_{m=0}^l [ C_filt(l,m,i)^2 + S_filt(l,m,i)^2 ]
    !      fst(l, i) = sum_{m=0}^l 2 * [ C_filt*CN_filt + S_filt*SN_filt ]
    subroutine compute_filtered_spectral_energy(coeffs_filt, lmax, n_lm, n_shells, fse, &
                                                coeffs_N_filt)
        integer,           intent(in)  :: lmax, n_lm, n_shells
        real(kind=real64), intent(in)  :: coeffs_filt(2, n_lm, n_shells)
        real(kind=real64), intent(out) :: fse(0:lmax, n_shells)
        real(kind=real64), intent(in),  optional :: coeffs_N_filt(2, n_lm, n_shells)

        real(kind=real64) :: e_sum
        integer :: i, l, m, lm_base
        logical :: calc_flux_local

        calc_flux_local = present(coeffs_N_filt)

        !$OMP PARALLEL DO PRIVATE(i, l, lm_base, e_sum, m) SCHEDULE(STATIC)
        do i = 1, n_shells
            do l = 0, lmax
                lm_base = (l * (l + 1)) / 2 + 1
                e_sum = 0.0_real64
                do m = 0, l
                    if (calc_flux_local) then
                        e_sum = e_sum + 2.0_real64 * ( &
                        coeffs_filt(1, lm_base + m, i) * coeffs_N_filt(1, lm_base + m, i) + &
                        coeffs_filt(2, lm_base + m, i) * coeffs_N_filt(2, lm_base + m, i) )
                    else
                        e_sum = e_sum + coeffs_filt(1, lm_base + m, i)**2 &
                                + coeffs_filt(2, lm_base + m, i)**2
                    end if
                end do
                fse(l, i) = e_sum

            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine compute_filtered_spectral_energy

    ! 6. Computes the 3D volume spectrum for ONE velocity component by sweeping
    !    radial Gaussian filter scales (MSCG coarse-graining), accumulating on the fly into binner,
    !    and optionally accumulates scale-to-scale energy transfers into binner_TV.
    subroutine compute_volume_spectrum_comp(coeffs_comp, lmax, n_lm, n_shells, comp, binner, &
                                            coeffs_N_comp, binner_TV)
        integer,                 intent(in)    :: lmax, n_lm, n_shells, comp
        real(kind=real64),       intent(in)    :: coeffs_comp(2, n_lm, n_shells)
        type(spectral_binner_t), intent(inout) :: binner
        real(kind=real64),       intent(in),    optional :: coeffs_N_comp(2, n_lm, n_shells)
        type(spectral_binner_t), intent(inout), optional :: binner_TV

        ! Filter Scale Energy (fse) and Higher Scale Energy (hse)
        real(kind=real64), allocatable :: fse(:,:), hse(:,:)
        real(kind=real64), allocatable :: fst(:,:), hst(:,:)
        real(kind=real64), allocatable :: coeffs_filt(:,:,:), coeffs_N_filt(:,:,:)
        real(kind=real64) :: ell
        integer :: k
        logical :: calc_flux_local

        calc_flux_local = present(coeffs_N_comp) .and. present(binner_TV)

        allocate(coeffs_filt(2, n_lm, n_shells))
        allocate(hse(0:lmax, n_shells), fse(0:lmax, n_shells))
        if (calc_flux_local) then
            allocate(coeffs_N_filt(2, n_lm, n_shells))
            allocate(hst(0:lmax, n_shells), fst(0:lmax, n_shells))
        end if

        ! k = 0 means infinite filtering scale, i.e. simple radial averaging across all shells
        call compute_filtered_sht_coeffs(coeffs_comp, n_lm, n_shells, -1.0_real64, coeffs_filt)
        call compute_filtered_spectral_energy(coeffs_filt, lmax, n_lm, n_shells, fse)
        call binner%accumulate_volume_energy(fse, 0, comp)
        hse = fse
        if (calc_flux_local) then
            call compute_filtered_sht_coeffs(coeffs_N_comp, n_lm, n_shells, -1.0_real64, coeffs_N_filt)
            call compute_filtered_spectral_energy(coeffs_filt, lmax, n_lm, n_shells, fst, &
                                                  coeffs_N_filt)
            call binner_TV%accumulate_volume_energy(fst, 0, comp)
            hst = fst
        end if

        ! ---- Sweep over radial filter scales ----
        do k = 1, n_shells - 1
            ell = 0.5_real64 * real(n_shells, real64) / sqrt(real(k*(k+1), real64))

            call compute_filtered_sht_coeffs(coeffs_comp, n_lm, n_shells, ell, coeffs_filt)
            call compute_filtered_spectral_energy(coeffs_filt, lmax, n_lm, n_shells, fse)
            call binner%accumulate_volume_energy(fse - hse, k, comp)
            hse = fse
            if (calc_flux_local) then
                call compute_filtered_sht_coeffs(coeffs_N_comp, n_lm, n_shells, ell, coeffs_N_filt)
                call compute_filtered_spectral_energy(coeffs_filt, lmax, n_lm, n_shells, fst, &
                                                      coeffs_N_filt)
                call binner_TV%accumulate_volume_energy(fst - hst, k, comp)
                hst = fst
            end if

            if ((mod(k, n_shells/4) == 0) .and. (n_shells > 64)) write(*,'(A,I4,A,I4,A,I2)') &
                "      k = ", k, " * 2*pi/(", n_shells, "dr),   component ", comp
        end do   ! k

        ! Unfiltered spectral energy (k_max => ell = 0; pass raw coefficients directly)
        call compute_filtered_spectral_energy(coeffs_comp, lmax, n_lm, n_shells, fse)
        call binner%accumulate_volume_energy(fse - hse, n_shells, comp)
        if (calc_flux_local) then
            call compute_filtered_spectral_energy(coeffs_comp, lmax, n_lm, n_shells, fst, &
            coeffs_N_comp)
            call binner_TV%accumulate_volume_energy(fst - hst, n_shells, comp)
        end if
        
        deallocate(hse, fse, coeffs_filt)
        if (calc_flux_local) deallocate(hst, fst, coeffs_N_filt)
    end subroutine compute_volume_spectrum_comp

    ! 7. Wrapper computing both 2D horizontal and 3D volume spectra (and optional transfers)
    !    for a single spectral component, then freeing its coefficient buffers.
    subroutine process_component_spectra(coeffs_u, lmax, n_lm, n_shells, comp, &
                                         binner_H, binner_E, coeffs_N, binner_TH, binner_TV)
        integer,                 intent(in)    :: lmax, n_lm, n_shells, comp
        real(kind=real64), allocatable, intent(inout) :: coeffs_u(:,:,:)
        type(spectral_binner_t), intent(inout) :: binner_H, binner_E
        real(kind=real64), allocatable, intent(inout), optional :: coeffs_N(:,:,:)
        type(spectral_binner_t), intent(inout), optional :: binner_TH, binner_TV

        logical :: calc_flux_local

        write(*, '(A, I1, A)') "   -> Computing 2D Horizontal & 3D Volume Spectra (component ", comp, "/3)..."

        calc_flux_local = present(coeffs_N) .and. present(binner_TH) .and. present(binner_TV)
        if (calc_flux_local) calc_flux_local = allocated(coeffs_N)

        if (calc_flux_local) then
            call compute_horizontal_spectrum_comp(coeffs_u, lmax, n_lm, n_shells, binner_H, &
                                                  coeffs_N, binner_TH)
            call compute_volume_spectrum_comp(coeffs_u, lmax, n_lm, n_shells, comp, binner_E, &
                                              coeffs_N, binner_TV)
            deallocate(coeffs_N)
        else
            call compute_horizontal_spectrum_comp(coeffs_u, lmax, n_lm, n_shells, binner_H)
            call compute_volume_spectrum_comp(coeffs_u, lmax, n_lm, n_shells, comp, binner_E)
        end if
        deallocate(coeffs_u)
    end subroutine process_component_spectra


end module sht_processing_mod


program sht_3d_utility
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    use data_io_mod
    use binner_mod, only: spectral_binner_t
    use advection_mod, only: compute_advection_shell, shell_transfer
    use sht_processing_mod
    implicit none

    integer :: n1, n2, n3, lmax, n_lat, n_lon, n_lm
    integer :: arg_count, i, j, k, io_status, last_slash
    real(kind=real64) :: sim_time, pi, phi, sp

    character(len=256) :: input_path, coord_type, out_file, r_arg, arg
    integer :: arg_idx
    logical :: r_found, calc_flux
    integer, parameter :: out_unit = 11, vol_unit = 12, count_unit = 13, flux_unit = 14

    ! Velocity in spherical coordinates (vt=azimuthal/zonal, vr=radial, vp=polar/meridional)
    real(kind=real64), allocatable :: vt(:,:,:), vr(:,:,:), vp(:,:,:)

    ! Spectral arrays and binners (packed triangularly as (2, n_lm, n2))
    real(kind=real64), allocatable :: coeffs_r(:,:,:), coeffs_tor(:,:,:), coeffs_sph(:,:,:)
    type(spectral_binner_t) :: binner_H, binner_E
    real(kind=real64) :: rint, rext, dr
    integer :: n_bins

    ! Nonlinear advection and flux variables
    real(kind=real64), allocatable :: Nr_shell(:,:), Nt_shell(:,:), Np_shell(:,:)
    real(kind=real64), allocatable :: coeffs_Nr(:,:,:), coeffs_Ntor(:,:,:), coeffs_Nsph(:,:,:)
    real(kind=real64), allocatable :: Pi_tot(:), Pi_1(:), Pi_2(:)
    type(spectral_binner_t) :: binner_TH, binner_TV
    real(kind=real64) :: T_phys, T_binned_sum, Pi_max

    real(kind=real64) :: E_phys, d_theta, d_phi, r_j
    real(kind=real64) :: E_binned_sum, delta_log, w_hi

    ! 1. Parse Command Line
    arg_count = command_argument_count()
    if (arg_count >= 1) then
        call get_command_argument(1, input_path)
        if (trim(input_path) == '--help' .or. trim(input_path) == '-h') then
            write(*, '("Usage: ./sht_3d_utility <file_path> -r=<r_in> [-p|-c] [-f|--flux]")')
            write(*, '(" -r=<r_in>      : inner radius (required); r_ext = r_in + 1")')
            write(*, '(" -p (default)   = velocity field input is in spherical coordinates")')
            write(*, '(" -c             = velocity field input is in Cartesian coordinates")')
            write(*, '(" -f, --flux     = also compute horizontal and 3D volume energy fluxes")')
            stop 0
        end if
    end if

    if (arg_count < 2) then
        write(error_unit, '("Usage: ./sht_3d_utility <file_path> -r=<r_in> [-p|-c] [-f|--flux]")')
        write(error_unit, '(" -r=<r_in>      : inner radius (required); r_ext = r_in + 1")')
        write(error_unit, '(" -p (default)   = velocity field input is in spherical coordinates")')
        write(error_unit, '(" -c             = velocity field input is in Cartesian coordinates")')
        write(error_unit, '(" -f, --flux     = also compute horizontal and 3D volume energy fluxes")')
        write(error_unit, '("Example: ./sht_3d_utility dati/spbou001234.q -r=1.0 -c -f")')
        stop 1
    end if

    call get_command_argument(1, input_path)

    ! Default values
    coord_type = '-p'
    calc_flux = .false.
    r_found = .false.

    do arg_idx = 2, arg_count
        call get_command_argument(arg_idx, arg)
        if (arg(1:3) == '-r=') then
            read(arg(4:), *, iostat=io_status) rint
            if (io_status /= 0) then
                write(error_unit, '("Error: could not parse r_in from: ", A)') trim(arg)
                stop 1
            end if
            r_found = .true.
        else if (trim(arg) == '-p' .or. trim(arg) == '-c') then
            coord_type = arg
        else if (trim(arg) == '-f' .or. trim(arg) == '--flux') then
            calc_flux = .true.
        else
            write(error_unit, '("Warning: unknown argument skipped: ", A)') trim(arg)
        end if
    end do

    if (.not. r_found) then
        write(error_unit, '("Error: inner radius -r=<r_in> is required.")')
        stop 1
    end if

    rext = rint + 1.0_real64
    write(*, '(A, F10.4, A, F10.4)') "   Inner radius: ", rint, "   Outer radius: ", rext

    ! 2. Read Data
    call read_simulation_data(input_path, n1, n2, n3, sim_time, vt, vr, vp, coord_type)
    write(*, '(A, 3(I0, A), F10.4)') "   Dimensions: ", n1, "x", n2, "x", n3, " at time ", sim_time

    ! Compute radial spacing (uniform grid, shell thickness = 1)
    dr = (rext - rint) / real(n2 - 1, real64)

    ! Calculate total physical energy from the spatial grid
    pi = 2.0_real64 * asin(1.0_real64)
    E_phys = 0.0_real64
    d_phi = pi / real(n3-1, real64)
    d_theta = 2.0_real64 * pi / real(n1-1, real64)
    do k = 1, n3
        phi = real(k-1, real64) * pi / real(n3-1, real64)
        sp = sin(phi)
        do j = 1, n2
            r_j = rint + real(j-1, real64) * dr
            ! Only sum i=1 to n1-1 since i=n1 is a periodic copy of i=1
            do i = 1, n1-1
                E_phys = E_phys + 0.5_real64 * (vt(i,j,k)**2 + vr(i,j,k)**2 + vp(i,j,k)**2) * &
                         r_j**2 * dr * d_theta * d_phi * sp
            end do
        end do
    end do
    write(*,'(A, ES14.7)') "   Total Physical Energy (Grid sum) : ", E_phys

    ! 3. Setup SHT Parameters and Binners
    n_lat = n3 - mod(n3, 2)
    n_lon = n1 - 1
    lmax = n_lat / 2 - 1
    n_lm = ((lmax + 1) * (lmax + 2)) / 2

    ! Number of logarithmic bins (~12 per octave)
    n_bins = max(64, nint(12 * log(real(max(lmax, n2), real64)) / log(2.0_real64)))

    call binner_H%init(lmax, n2, rint, rext, n_bins)
    call binner_E%init(lmax, n2, rint, rext, n_bins)

    if (calc_flux) then
        allocate(Pi_tot(0:n_bins), Pi_1(0:n_bins), Pi_2(0:n_bins))
        call binner_TH%init(lmax, n2, rint, rext, n_bins)
        call binner_TV%init(lmax, n2, rint, rext, n_bins)
        T_phys = 0.0_real64
    end if

    ! ===== Stage 1: Radial component (scalar SHT, comp = 2) =====
    write(*, '(A)') "-> Stage 1: Processing radial component (SHT + MSCG)..."
    allocate(coeffs_r(2, n_lm, n2))
    if (calc_flux) then
        allocate(Nr_shell(n1, n3), coeffs_Nr(2, n_lm, n2))
    end if

    do j = 1, n2
        call compute_sht_coeffs(vr(1:n_lon, j, 1:n_lat), n_lat, n_lon, lmax, n_lm, &
                                coeffs_r(:, :, j))

        if (calc_flux) then
            r_j = rint + real(j - 1, real64) * dr
            call compute_advection_shell(vt, vr, vp, n1, n2, n3, rint, dr, j, 2, Nr_shell)
            T_phys = T_phys + shell_transfer(vr(:, j, :), Nr_shell, n1, n3, r_j, dr, d_theta, d_phi)
            call compute_sht_coeffs(Nr_shell(1:n_lon, 1:n_lat), n_lat, n_lon, lmax, n_lm, &
                                    coeffs_Nr(:, :, j))
        end if

        if (mod(j, 25) == 0) write(*,*) "      Progress: shell ", j, "/", n2
    end do

    if (calc_flux) then
        deallocate(Nr_shell)
    else
        deallocate(vr)
    end if

    call process_component_spectra(coeffs_r, lmax, n_lm, n2, 2, &
                                   binner_H, binner_E, coeffs_Nr, binner_TH, binner_TV)

    ! ===== Stage 2: Tangential components (comp = 1 and comp = 3) =====
    write(*, '(A)') "-> Stage 2: Processing tangential components (1/3 and 3/3, SHT)..."
    allocate(coeffs_tor(2, n_lm, n2), coeffs_sph(2, n_lm, n2))
    if (calc_flux) then
        allocate(Nt_shell(n1, n3), Np_shell(n1, n3))
        allocate(coeffs_Ntor(2, n_lm, n2), coeffs_Nsph(2, n_lm, n2))
    end if

    do j = 1, n2
        call compute_sht_coeffs(vt(1:n_lon, j, 1:n_lat), n_lat, n_lon, lmax, n_lm, &
                                coeffs_tor(:, :, j))
        call compute_sht_coeffs(vp(1:n_lon, j, 1:n_lat), n_lat, n_lon, lmax, n_lm, &
                                coeffs_sph(:, :, j))

        if (calc_flux) then
            r_j = rint + real(j - 1, real64) * dr
            call compute_advection_shell(vt, vr, vp, n1, n2, n3, rint, dr, j, 1, Nt_shell)
            call compute_advection_shell(vt, vr, vp, n1, n2, n3, rint, dr, j, 3, Np_shell)
            T_phys = T_phys + shell_transfer(vt(:, j, :), Nt_shell, n1, n3, r_j, dr, d_theta, d_phi) &
                            + shell_transfer(vp(:, j, :), Np_shell, n1, n3, r_j, dr, d_theta, d_phi)
            call compute_sht_coeffs(Nt_shell(1:n_lon, 1:n_lat), n_lat, n_lon, lmax, n_lm, &
                                    coeffs_Ntor(:, :, j))
            call compute_sht_coeffs(Np_shell(1:n_lon, 1:n_lat), n_lat, n_lon, lmax, n_lm, &
                                    coeffs_Nsph(:, :, j))
        end if

        if (mod(j, 25) == 0) write(*,*) "      Progress: shell ", j, "/", n2
    end do

    deallocate(vt, vp)
    if (calc_flux) then
        deallocate(vr, Nt_shell, Np_shell)
    end if

    call process_component_spectra(coeffs_tor, lmax, n_lm, n2, 1, &
                                   binner_H, binner_E, coeffs_Ntor, binner_TH, binner_TV)
    call process_component_spectra(coeffs_sph, lmax, n_lm, n2, 3, &
                                   binner_H, binner_E, coeffs_Nsph, binner_TH, binner_TV)

    ! 4. Output Horizontal Spectra
    call binner_H%normalize()

    last_slash = scan(trim(input_path), '/', back=.true.)
    if (last_slash > 0) then
        out_file = input_path(1:last_slash) // 'h_spec.dat'
    else
        out_file = 'h_spec.dat'
    end if
    write(*,*) "-> Saving Horizontal Spectra to: ", trim(out_file)
    open(unit=out_unit, file=trim(out_file), status='replace')
    write(out_unit, '(A)') "# k | E_tot_H(k) | E_zonal_H(k) | E_rest_H(k)"
    do i = 0, n_bins
        if (binner_H%val_tot(i) > 0.0_real64) then
            write(out_unit, '(ES15.7, ES15.7, ES15.7, ES15.7)') &
                binner_H%K_bins(i), binner_H%val_tot(i), binner_H%val_1(i), binner_H%val_2(i)
        end if
    end do
    close(out_unit)

    ! Output Horizontal Flux Spectra if requested
    if (calc_flux) then
        call binner_TH%normalize()
        call binner_TH%compute_flux(Pi_tot, Pi_1, Pi_2)

        if (last_slash > 0) then
            out_file = input_path(1:last_slash) // 'h_flux.dat'
        else
            out_file = 'h_flux.dat'
        end if
        write(*,*) "-> Saving Horizontal Flux Spectra to: ", trim(out_file)
        open(unit=flux_unit, file=trim(out_file), status='replace')
        write(flux_unit, '(A)') "# k | Pi_tot_H(k) | Pi_zonal_H(k) | Pi_rest_H(k) | T_tot_H(k) | T_zonal_H(k) | T_rest_H(k)"
        do i = 0, n_bins
            if (binner_TH%count_tot(i) > 0.0_real64) then
                write(flux_unit, '(7(ES15.7))') &
                    binner_TH%K_bins(i), Pi_tot(i), Pi_1(i), Pi_2(i), &
                    binner_TH%val_tot(i), binner_TH%val_1(i), binner_TH%val_2(i)
            end if
        end do
        close(flux_unit)
        call binner_TH%destroy()
    end if

    call binner_H%destroy()

    ! 5. Output 3D Volume Spectrum E(K)
    call binner_E%normalize()

    ! Parseval check: integrate binned spectrum as sum(E_K * dK)
    delta_log = log(binner_E%K_bins(2) / binner_E%K_bins(1))
    w_hi = exp(0.5_real64 * delta_log) - exp(-0.5_real64 * delta_log)
    
    E_binned_sum = binner_E%val_tot(0) ! DC mode doesn't have a width in the log bins
    do i = 1, n_bins
        E_binned_sum = E_binned_sum + binner_E%val_tot(i) * binner_E%K_bins(i) * w_hi
    end do

    write(*,'(A)')          "   --- Energy conservation check (Parseval thm.) ---"
    write(*,'(A, ES14.7)')  "       Physical energy (grid sum)   : ", E_phys
    write(*,'(A, ES14.7)')  "       Integral of binned spectrum  : ", E_binned_sum
    if (E_phys > 0.0_real64) then
        write(*,'(A, F8.4, A)') "       Relative difference: ", &
            abs(E_binned_sum - E_phys) / E_phys * 100.0_real64, " %"
    end if

    ! Write volume spectrum
    if (last_slash > 0) then
        out_file = input_path(1:last_slash) // 'v_spec.dat'
    else
        out_file = 'v_spec.dat'
    end if
    write(*,*) "-> Saving Volume Spectrum to: ", trim(out_file)
    open(unit=vol_unit, file=trim(out_file), status='replace')
    write(vol_unit, '(A)') "# k | E(K) | E_2d(k) | E_3d(k)"
    do i = 0, n_bins
        if (binner_E%count_tot(i) > 0.0_real64) then
            write(vol_unit, '(ES15.7, ES15.7, ES15.7, ES15.7)') &
                binner_E%K_bins(i), binner_E%val_tot(i), binner_E%val_1(i), binner_E%val_2(i)
        end if
    end do
    close(vol_unit)

    ! Write binning count files
    if (last_slash > 0) then
        out_file = input_path(1:last_slash) // 'binning.dat'
    else
        out_file = 'binning.dat'
    end if
    write(*,*) "-> Saving Binning Counts to: ", trim(out_file)
    open(unit=count_unit, file=trim(out_file), status='replace')
    write(count_unit, '(A)') "# k | count_tot | count_2d | count_3d"
    do i = 0, n_bins
        write(count_unit, '(ES15.7, ES15.7, ES15.7, ES15.7)') &
            binner_E%K_bins(i), binner_E%count_tot(i), binner_E%count_1(i), binner_E%count_2(i)
    end do
    close(count_unit)

    ! Output Volume Flux Spectra if requested
    if (calc_flux) then
        call binner_TV%normalize()
        call binner_TV%compute_flux(Pi_tot, Pi_1, Pi_2)

        T_binned_sum = -(Pi_tot(n_bins) - 0.5_real64 * binner_TV%val_tot(n_bins) * binner_TV%K_bins(n_bins) * w_hi)
        Pi_max = maxval(abs(Pi_tot))

        write(*,'(A)')          "   --- Nonlinear transfer conservation check ---"
        write(*,'(A, ES14.7)')  "       Physical transfer (grid sum) : ", T_phys
        write(*,'(A, ES14.7)')  "       Integral of binned transfer  : ", T_binned_sum
        write(*,'(A, ES14.7)')  "       Peak energy flux max|Pi(K)|  : ", Pi_max
        if (Pi_max > 0.0_real64) then
            write(*,'(A, F8.4, A)') "       Flux closure |Pi_end|/max|Pi|: ", &
                abs(T_binned_sum) / Pi_max * 100.0_real64, " %"
            write(*,'(A, F8.4, A)') "       Grid-spectral diff / max|Pi| : ", &
                abs(T_binned_sum - T_phys) / Pi_max * 100.0_real64, " %"
        end if

        if (last_slash > 0) then
            out_file = input_path(1:last_slash) // 'v_flux.dat'
        else
            out_file = 'v_flux.dat'
        end if
        write(*,*) "-> Saving Volume Flux Spectra to: ", trim(out_file)
        open(unit=flux_unit, file=trim(out_file), status='replace')
        write(flux_unit, '(A)') "# k | Pi(K) | Pi_2d(K) | Pi_3d(K) | T(K) | T_2d(K) | T_3d(K)"
        do i = 0, n_bins
            if (binner_TV%count_tot(i) > 0.0_real64) then
                write(flux_unit, '(7(ES15.7))') &
                    binner_TV%K_bins(i), Pi_tot(i), Pi_1(i), Pi_2(i), &
                    binner_TV%val_tot(i), binner_TV%val_1(i), binner_TV%val_2(i)
            end if
        end do
        close(flux_unit)
        call binner_TV%destroy()
        deallocate(Pi_tot, Pi_1, Pi_2)
    end if

    call binner_E%destroy()

    write(*,*) "-> Done."

end program sht_3d_utility
