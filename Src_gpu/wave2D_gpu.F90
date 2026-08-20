program wave2d_solver_gpu
    use cudafor
    implicit none

    ! Grid Parameters (0 to Nx, 0 to Ny)
    integer, parameter :: Nx = 6000            ! Max grid index in x
    integer, parameter :: Ny = 12000            ! Max grid index in y
    real, parameter :: dx = 2.0e3              ! Grid spacing
    real, parameter :: dy = 2.0e3              ! Grid spacing
    
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

    real :: t  !, dudx, dudy
    integer :: it, i, j
    integer :: know, knew, kold               ! Index pointers

    ! dissipation parameter: 
    real, parameter :: Ak = 2000.

    !!!!!!!!!!!!!
    ! Define arrays: 
    real, dimension(1:Nx, 0:Ny)     :: Fu   ! u-points (x-interfaces)
    real, dimension(0:Nx, 1:Ny)     :: Fv   ! v-points (y-interfaces)
    real, dimension(0:Nx, 0:Ny, 2) :: u, rhs

    real, device :: Fu_d(1:Nx, 0:Ny)        ! u-points (x-interfaces)
    real, device :: Fv_d(0:Nx, 1:Ny)          ! v-points (y-interfaces)
    real, device :: u_d(0:Nx, 0:Ny, 2), rhs_d(0:Nx, 0:Ny, 2)

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

    Fu_d = 0.
    Fv_d = 0.
    u_d = u
    rhs_d = 0.

    call cpu_time(tcpu1)
    
    do it = 1, nstep

!        t = dt*it

        know = knew
        knew = 3 - know
        kold = 3 - know  ! Points to rhs from previous step

        ! ---------------------------------------------------------
        ! Compute RHS for INTERIOR points only (1 to Nx-1, 1 to Ny-1)
        ! ---------------------------------------------------------

        !- Advection (also initialize rhs(:,:,know):

        !$cuf kernel do (2) <<<*,*>>>        
        do j = 1, Ny - 1
            do i = 1, Nx - 1

!                dudx = (u(i+1, j, know) - u(i-1, j, know)) * odx2
!                dudy = (u(i, j+1, know) - u(i, j-1, know)) * ody2

                rhs_d(i, j, know) = - Cx * (u_d(i+1, j, know) - u_d(i-1, j, know)) * odx2  & 
                                    - Cy * (u_d(i, j+1, know) - u_d(i, j-1, know)) * ody2

            end do
        end do
        i = cudaDeviceSynchronize ()


        !- Add dissipation: 

        !$cuf kernel do (2) <<<*,*>>>
        do j = 1, Ny - 1
         do i = 1, Nx
          Fu_d(i, j) = Ak * (u_d(i, j, know) - u_d(i-1, j, know)) * odx
         end do
        end do
  
        !$cuf kernel do (2) <<<*,*>>>
        do j = 1, Ny
         do i = 1, Nx - 1
          Fv_d(i, j) = Ak * (u_d(i, j, know) - u_d(i, j-1, know)) * ody
         end do
        end do
        i = cudaDeviceSynchronize ()

        !$cuf kernel do (2) <<<*,*>>>
        do j = 1, Ny - 1
         do i = 1, Nx - 1

          rhs_d(i, j, know) = rhs_d(i, j, know) + &
                          (Fu_d(i+1, j) - Fu_d(i, j)) * odx + &
                          (Fv_d(i, j+1) - Fv_d(i, j)) * ody

         end do
        end do
        i = cudaDeviceSynchronize ()

        ! ---------------------------------------------------------
        ! Update u(:,:,knew) for INTERIOR points using explicit loops
        ! ---------------------------------------------------------
        if (it == 1) then
            ! Step 1: Forward Euler
          
            !$cuf kernel do (2) <<<*,*>>>
            do j = 1, Ny - 1
                do i = 1, Nx - 1
                    u_d(i, j, knew) = u_d(i, j, know) + dt * rhs_d(i, j, know)
                end do
            end do
!            i = cudaDeviceSynchronize ()

        else
            ! Steps 2..N: Adams-Bashforth 2

            !$cuf kernel do (2) <<<*,*>>>
            do j = 1, Ny - 1
                do i = 1, Nx - 1
                    u_d(i, j, knew) = u_d(i, j, know) +             &
                                       dt * (1.5 * rhs_d(i, j, know) - 0.5 * rhs_d(i, j, kold))
                end do
            end do
!            i = cudaDeviceSynchronize ()
        end if

        ! ---------------------------------------------------------
        ! Apply Double Periodic Boundary Conditions via Explicit Loops
        ! ---------------------------------------------------------

        !$cuf kernel do (1) <<<*,*>>>        
        do j = 1, Ny - 1
            u_d(0,  j, knew) = u_d(Nx-1, j, knew)
            u_d(Nx, j, knew) = u_d(1,    j, knew)
        end do

        !$cuf kernel do (1) <<<*,*>>>
        do i = 1, Nx - 1
            u_d(i, 0,  knew) = u_d(i, Ny-1, knew)
            u_d(i, Ny, knew) = u_d(i, 1,    knew)
        end do

        ! Corner Points
        u_d(0,  0,  knew) = u_d(Nx-1, Ny-1, knew)
        u_d(Nx, 0,  knew) = u_d(1,    Ny-1, knew)
        u_d(0,  Ny, knew) = u_d(Nx-1, 1,    knew)
        u_d(Nx, Ny, knew) = u_d(1,    1,    knew)

        if (mod(it,NHIS) == 0) then
         u(:,:,knew) = u_d(:,:,knew)
         t=dt*it
         call write_his(hisname,u(:,:,knew),Nx+1,Ny+1,t)
        end if

    end do
    
    call cpu_time(tcpu2)
    write (*,*) 'tcpu2-tcpu1=',tcpu2-tcpu1


end program wave2d_solver_gpu

