program river_simulation
    use simulation_params  ! Global parameters module
    implicit none

    ! Variable declarations
    real(8) :: dt, time, tend, dtmax, CFL, dtk, durata_tot
    real(8) :: xL, xR, yL, yR
    real(8), dimension(IMAX) :: xb, dx_b, dx_int  ! `dx_b` now an array for non-uniform dx
    real(8), dimension(IMAX+1) :: xj, q, q0, qstar, uf, Af
    real(8), dimension(IMAX) :: A, Rh, eta, b_min, b_max, eta2d
    integer :: M, n_s, i, j

    ! ---------------------------------------------------------------
    ! Step 0: Read user-defined variables from input file
    ! ---------------------------------------------------------------
    call read_user_input("input.txt", xL, xR, yL, yR, dtmax, CFL, tend, M)

    ! ---------------------------------------------------------------
    ! Step 1: Read hydrograph data from file
    ! ---------------------------------------------------------------
    call read_hydrograph("data/hydrograph.dat", M, dtk, durata_tot)

    ! ---------------------------------------------------------------
    ! Step 2: Read or calculate non-uniform dx values
    ! ---------------------------------------------------------------
    call read_dx_values("data/dx_values.dat", dx_b, IMAX)
    ! This subroutine reads the `dx` values from a file and fills the `dx_b` array.

    ! ---------------------------------------------------------------
    ! Step 3: Set up computational domain and mesh based on non-uniform dx
    ! ---------------------------------------------------------------
    call define_mesh_with_nonconstant_dx(xb, xj, dx_b, IMAX, xL)
    ! This subroutine sets the mesh points `xb` and `xj` based on the non-uniform `dx_b` array.

    ! ---------------------------------------------------------------
    ! Step 4: Set up geometry and bed profile
    ! ---------------------------------------------------------------
    call compute_geometry(xb, eta, b_min, b_max, IMAX)

    ! ---------------------------------------------------------------
    ! Step 5: Initialize hydraulic variables
    ! ---------------------------------------------------------------
    call initialize_hydraulics(A, Rh, q, q0, IMAX)

    ! ---------------------------------------------------------------
    ! Step 6: Main time-stepping loop
    ! ---------------------------------------------------------------
    time = 0.0     ! Start time
    n_s = 0        ! Time step counter

    do while (time < tend .and. n_s < 1000000)
        ! Step 6.1: Compute velocity and time step based on non-uniform dx
        call compute_velocity_and_dt_with_nonconstant_dx(q, Af, uf, dt, dx_b, IMAX, CFL, dtmax)
        ! This subroutine computes velocities at each interface using non-uniform dx.

        ! Step 6.2: Update discharge (qstar)
        call update_discharge(qstar, q, A, dt, time, dtk, M)

        ! Step 6.3: Update water surface elevation and area
        call update_eta_and_area(eta, A, qstar, dt, dx_b, IMAX)

        ! Step 6.4: Update simulation time and counters
        time = time + dt
        n_s = n_s + 1
    end do

    ! ---------------------------------------------------------------
    ! Step 7: Output final results
    ! ---------------------------------------------------------------
    call output_results("simulations/output.dat", xb, b_min, eta, q, IMAX)

end program river_simulation

! ================================================================

! Subroutine to read dx values from a file
subroutine read_dx_values(filename, dx_b, IMAX)
    implicit none
    character(len=*), intent(in) :: filename
    real(8), dimension(IMAX), intent(out) :: dx_b
    integer :: i, ios

    open(unit=10, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening dx values file: ", filename
        stop
    endif

    do i = 1, IMAX
        read(10, *, iostat=ios) dx_b(i)
        if (ios /= 0) then
            print *, "Error reading dx value for index: ", i
            stop
        endif
    end do
    close(10)

end subroutine read_dx_values

! ================================================================

! Subroutine to define mesh using non-uniform dx
subroutine define_mesh_with_nonconstant_dx(xb, xj, dx_b, IMAX, xL)
    implicit none
    real(8), intent(out) :: xb(IMAX), xj(IMAX+1)
    real(8), intent(in) :: dx_b(IMAX)
    integer, intent(in) :: IMAX
    real(8) :: xL
    integer :: i

    ! First point is at xL
    xb(1) = xL
    xj(1) = xL - dx_b(1) / 2.0

    ! Loop to calculate the rest of the mesh points
    do i = 2, IMAX
        xb(i) = xb(i-1) + dx_b(i-1)
        xj(i) = xb(i-1) + dx_b(i-1) / 2.0
    end do

    ! Last j point
    xj(IMAX+1) = xb(IMAX) + dx_b(IMAX) / 2.0

end subroutine define_mesh_with_nonconstant_dx

! ================================================================

! Subroutine to compute velocity and time step using non-uniform dx
subroutine compute_velocity_and_dt_with_nonconstant_dx(q, Af, uf, dt, dx_b, IMAX, CFL, dtmax)
    implicit none
    real(8), dimension(IMAX+1), intent(in) :: q, Af
    real(8), dimension(IMAX+1), intent(out) :: uf
    real(8), dimension(IMAX), intent(in) :: dx_b
    real(8), intent(out) :: dt
    real(8), intent(in) :: CFL, dtmax
    integer, intent(in) :: IMAX
    integer :: i
    real(8) :: dx_min, dt_local

    ! Compute velocities at interfaces
    do i = 1, IMAX+1
        uf(i) = q(i) / Af(i)
    end do

    ! Compute the time step based on the minimum dx and velocity
    dx_min = minval(dx_b)
    dt_local = CFL * dx_min / maxval(abs(uf))

    ! Ensure the time step doesn't exceed the maximum allowed value
    dt = min(dt_local, dtmax)

end subroutine compute_velocity_and_dt_with_nonconstant_dx
