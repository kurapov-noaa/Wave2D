program wave2d_solver
    implicit none

    ! Grid Parameters (0 to Nx, 0 to Ny)
    integer, parameter :: Nx = 6000            ! Max grid index in x
    integer, parameter :: Ny = 12000            ! Max grid index in y
    real, parameter :: dx = 2.0e3               ! Grid spacing
    real, parameter :: dy = 2.0e3               ! Grid spacing
    
    ! The length of the periodic domain in each direction (max wave length) 
    real, parameter :: Lx = (Nx - 1) * dx     
    real, parameter :: Ly = (Ny - 1) * dy

    ! Derivative Inverse Multipliers
    real, parameter :: odx = 1.0 / dx    ! 1 / dx
    real, parameter :: ody = 1.0 / dy    ! 1 / dy
    real, parameter :: odx2 = 0.5 * odx  ! 1 / (2*dx)
    real, parameter :: ody2 = 0.5 * ody  ! 1 / (2*dy)

    ! Physical Parameters
    real, parameter :: Cx = 10.0              ! Coefficient Cx (m/s)
    real, parameter :: Cy =  5.0              ! Coefficient Cy (m/s)

    ! Time Stepping Parameters
    real, parameter :: dt = 40.0              ! Time step (s)
    integer, parameter :: nstep = 10000        ! Total number of time steps

    ! history file (place here, later will read from .in):
    character (len=256), parameter :: hisname = "wave2D_his.nc"

    ! Wave Parameters for Initial Condition
    real, parameter :: pi = 3.14159265
    real, parameter :: kx = 2.* 2.0 * pi / Lx     ! Wave number in x
    real, parameter :: ky = 1.* 2.0 * pi / Ly     ! Wave number in y

    integer, parameter :: NHIS = 10001

    real :: t, dudx, dudy
    integer :: it, i, j
    integer :: know, knew, kold               ! Index pointers

    ! dissipation parameter: 
    real, parameter :: Ak = 2000.

    !!!!!!!!!!!!!
    ! Define arrays: 
    real, allocatable     :: Fu(:,:)   ! u-points (x-interfaces)
    real, allocatable     :: Fv(:,:)   ! v-points (y-interfaces)
    real, allocatable     :: u(:,:,:), rhs(:,:,:)

    real :: tcpu1,tcpu2

    print '(A)', "----------------------------------------------------"
    print '(A)', " 2D Inviscid Wave Solver (Central Diff, Explicit Loops)"
    print '(A)', "----------------------------------------------------"
    print '(A, I0, A, I0)', " Domain Grid     : 0:", Nx, " x 0:", Ny
    print '(A, F10.1, A, F10.1)', " Domain Size (L) : Lx = ", Lx, " m, Ly = ", Ly, " m"
    print '(A, F8.1, A, F8.1)', " Grid Resolution : dx = ", dx, " m, dy = ", dy, " m"
    print '(A, F10.2, A)', " Coefficient Cx  : ", Cx, " m/s"
    print '(A, F10.2, A)', " Coefficient Cy  : ", Cy, " m/s"
    print '(A, F8.2, A)', " Time Step (dt)  : ", dt, " s"
    print '(A, I0)',       " Total Steps     : ", nstep
    print '(A)', "----------------------------------------------------"

    allocate( u(0:Nx, 0:Ny, 2) )
    allocate( rhs(0:Nx, 0:Ny, 2) )
    allocate( Fu(1:Nx, 0:Ny) )
    allocate( Fv(0:Nx, 1:Ny) )  


    ! -------------------------------------------------------------
    ! 1. Initial Conditions (Oblique Harmonic Wave)
    ! -------------------------------------------------------------
    knew = 1
    t = 0.0

    do j = 0, Ny
        do i = 0, Nx
            u(i, j, knew) = sin(kx * i * dx + ky * j * dy)
        end do
    end do

    ! Apply initial periodic BCs across the domain boundaries (explicit loops)
    do j = 1, Ny - 1
        u(0,  j, knew) = u(Nx-1, j, knew)
        u(Nx, j, knew) = u(1,    j, knew)
    end do

    do i = 1, Nx - 1
        u(i, 0,  knew) = u(i, Ny-1, knew)
        u(i, Ny, knew) = u(i, 1,    knew)
    end do

    u(0,  0,  knew) = u(Nx-1, Ny-1, knew)
    u(Nx, 0,  knew) = u(1,    Ny-1, knew)
    u(0,  Ny, knew) = u(Nx-1, 1,    knew)
    u(Nx, Ny, knew) = u(1,    1,    knew)

    call create_his(hisname,Nx+1,Ny+1)
    call write_his(hisname,u(:,:,knew),Nx+1,Ny+1,t)

    ! -------------------------------------------------------------
    ! 2. Main Time Loop
    ! -------------------------------------------------------------

    call cpu_time(tcpu1)

    do it = 1, nstep

!        t = dt*it

        know = knew
        knew = 3 - know
        kold = 3 - know  ! Points to rhs from previous step

        ! ---------------------------------------------------------
        ! Compute RHS for INTERIOR points only (1 to Nx-1, 1 to Ny-1)
        ! ---------------------------------------------------------

        ! Advection (also initialize rhs(:,:,know):
        do j = 1, Ny - 1
            do i = 1, Nx - 1

                dudx = (u(i+1, j, know) - u(i-1, j, know)) * odx2
                dudy = (u(i, j+1, know) - u(i, j-1, know)) * ody2

                rhs(i, j, know) = -Cx * dudx - Cy * dudy

            end do
        end do

        ! Add dissipation: 
        do j = 1, Ny - 1
         do i = 1, Nx
          Fu(i, j) = Ak * (u(i, j, know) - u(i-1, j, know)) * odx
         end do
        end do

        do j = 1, Ny
         do i = 1, Nx - 1
          Fv(i, j) = Ak * (u(i, j, know) - u(i, j-1, know)) * ody
         end do
        end do

        do j = 1, Ny - 1
         do i = 1, Nx - 1

          rhs(i, j, know) = rhs(i, j, know) + &
                          (Fu(i+1, j) - Fu(i, j)) * odx + &
                          (Fv(i, j+1) - Fv(i, j)) * ody
         end do
        end do

        ! ---------------------------------------------------------
        ! Update u(:,:,knew) for INTERIOR points using explicit loops
        ! ---------------------------------------------------------
        if (it == 1) then
            ! Step 1: Forward Euler
            do j = 1, Ny - 1
                do i = 1, Nx - 1
                    u(i, j, knew) = u(i, j, know) + dt * rhs(i, j, know)
                end do
            end do
        else
            ! Steps 2..N: Adams-Bashforth 2
            do j = 1, Ny - 1
                do i = 1, Nx - 1
                    u(i, j, knew) = u(i, j, know) + dt * (1.5 * rhs(i, j, know) - 0.5 * rhs(i, j, kold))
                end do
            end do
        end if

        ! ---------------------------------------------------------
        ! Apply Double Periodic Boundary Conditions via Explicit Loops
        ! ---------------------------------------------------------
        do j = 1, Ny - 1
            u(0,  j, knew) = u(Nx-1, j, knew)
            u(Nx, j, knew) = u(1,    j, knew)
        end do

        do i = 1, Nx - 1
            u(i, 0,  knew) = u(i, Ny-1, knew)
            u(i, Ny, knew) = u(i, 1,    knew)
        end do

        ! Corner Points
        u(0,  0,  knew) = u(Nx-1, Ny-1, knew)
        u(Nx, 0,  knew) = u(1,    Ny-1, knew)
        u(0,  Ny, knew) = u(Nx-1, 1,    knew)
        u(Nx, Ny, knew) = u(1,    1,    knew)


        if (mod(it,NHIS) == 0) then
         t = dt*it
         call write_his(hisname,u(:,:,knew),Nx+1,Ny+1,t)
        end if

    end do

    call cpu_time(tcpu2)
    write (*,*) 'tcpu2-tcpu1=',tcpu2-tcpu1

    deallocate(u,rhs,Fu,Fv)

end program wave2d_solver

