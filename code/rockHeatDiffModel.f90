PROGRAM CERN

IMPLICIT NONE

! Declaration of variables
INTEGER :: i, n, nnx, nstep
REAL (8) :: days, years, SeasonA, pi
REAL (8) :: total_time, time, dt, Tbc_left, Tbc_right, kappa, Lx, h
REAL (8), DIMENSION (:), ALLOCATABLE :: T, x, T_new
CHARACTER (30) :: filename

! State nnx
PRINT *, 'Give the total depth, LX [m] and spacing, h [m]'
READ *, Lx, h
PRINT *, 'Give the total time [days] and the time step size [days]'
READ *, total_time, dt
PRINT *, 'State the output filename (max 30 characters)'
READ *, filename

! Problem definition at t = 0
nnx = (Lx/h)+1

pi = 3.141593
days = 3600*24
years = days*365
total_time = total_time*days
dt = dt*days
nstep = total_time/dt

! Left = surface, right = at depth
Tbc_left = 10
Tbc_right = 10
SeasonA = 15
kappa = 1d-6

ALLOCATE (T(1:nnx), x(1:nnx), T_new(1:nnx))


    ! Calculate locations x along beam
        DO i = 1, nnx
            x(i) = (i-1)*h
        ENDDO
        
    ! Calculate temperature for all x, along beam at t = 0
        DO i = 1, nnx
            IF (i == 1) THEN
                T(i) = Tbc_left
            ELSE IF (i > 1) THEN
                T(i) = Tbc_right
            ENDIF
        ENDDO

    ! Print t = 0 scenario to file
    OPEN (unit = 10, file = 'T_init.txt')
        WRITE (10,*) '# position [m]     ', 'temp [degC]'
        DO i = 1, nnx
            WRITE (10,*) x(i), T(i)
        ENDDO
    CLOSE (10)


! Calculate temperature diffusion across the beam, for each position x
OPEN (unit = 11, file = ''//TRIM(filename)//'.txt')
WRITE (11,*) '# total time = ', time/days, '[days]'
WRITE (11,*) '# timestep = ', dt/days, '[days]'
WRITE (11,*) '# position [m]     ', 'time [days]     ', 'temp [degC]'

    ! Write initial conditions
    time = 0
    DO i = 1, nnx
        IF (MOD(i,50) == 1) THEN
            WRITE (11,*) x(i), time/days, T(i)
        ENDIF
    ENDDO

    WRITE (11,*) ''
    
    DO n = 2, nstep
        time = dt*(n-1)
    
        DO i = 2, (nnx-1)
            T_new(i) = T(i) + ((kappa*dt)/h**2)*(T(i+1) - 2*T(i) + T(i-1))
        ENDDO
    
        DO i = 1, nnx
            T_new(1) = SeasonA*sin(2*pi*time/years)+Tbc_left
            T_new(nnx) = Tbc_right
            IF (MOD(i,50) == 1) THEN
                WRITE (11,*) x(i), time/days, T_new(i)
            ENDIF
            T(i) = T_new(i)
        ENDDO
    
        WRITE (11,*) ''
    ENDDO
    
CLOSE (11)


END PROGRAM CERN

! Variables
! ------------
! nnx = number of points along Lx [-]
! nstep = number of time steps [-]
! time = time [s]
! dt = size of a time step [s]
! Tbc_left, Tb_right = temperature boundary condition at left and right, resp [degC]
! kappa = heat diffusivity, 1e-6 [m2/s]
! Lx = depth [m]
! h = position step [m]
! T = temperature at position x [degC]
! x = position along beam of length Lx [m]
! T_new = new temperature at location i
! T_old = temperature at location i in previous timestep



