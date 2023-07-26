  module fmot
    use output
    implicit none
    integer, parameter :: n_particles = 3
    real, parameter :: pi = 4.*atan(1.)

  contains

    subroutine main_loop
      use types

      integer :: i, j

      real, dimension(n_particles,3) :: x
      real, dimension(n_particles,3) :: v

      real, dimension(n_particles,6) :: state
      
      real :: box_size=2.
      real, dimension(2, 3) :: xlim

      real :: t, h
      integer :: write_output = 1


      ! f = ma = m dv/dt = m ddt(x)
      ! Runge-Kutta: dy/dt = g(t, y), g(t_0) = y_0
      ! g = f/m
      ! 

      t = 0.
      h = 0.1035

      call init_particles(x, v, speed=1., radius=1.)

      state(:,1:3) = x(:,1:3)
      state(:,4:6) = v(:,1:3)

      xlim(1, 1:3) = [-box_size,-box_size,-box_size] ! minimum x
      xlim(2, 1:3) = [box_size,box_size,box_size] ! maximum x

      call write_particles(state, 0, 0., n_particles)
      do i = 1, 100
        t = i*h
        state = state &
                + dRK4(state, t, h, pi_attracts_pj) !&
                !+ dRK4(state, t, h, constant_a_1)
        call set_walls(state, xlim)
        if (modulo(i, write_output) .eq. 0) then
          call write_particles(state, i, t, n_particles)
        end if
      end do
    end subroutine main_loop

    subroutine init_particles(x, v, speed, radius)
      real, dimension(n_particles,3) :: x
      real, dimension(n_particles,3) :: v
      real, intent(in) :: speed
      real, intent(in) :: radius
      real :: RM(3,3) = 0.
      integer :: i

      x = 0.
      v = 0.

      RM(1:3,1) = [cos(2.*pi/n_particles), -sin(2.*pi/n_particles), 0.]
      RM(1:3,2) = [sin(2.*pi/n_particles),  cos(2.*pi/n_particles), 0.]
      RM(1:3,3) = [0.                    , 0.                     , 1.]
               
      x(1, :) = [radius, 0.    , 0.]
      v(1, :) = [0.    , speed, 0.]
      do i = 2, n_particles
        x(i, :) = matmul(RM, x(i-1, :))
        v(i, :) = matmul(RM, v(i-1, :))
      end do

    end subroutine

    subroutine set_walls(state, xlim)
      real, dimension(n_particles,6) :: state
      real, dimension(2, 3), intent(in) :: xlim
      integer :: j

      do j = 1, n_particles
        if (state(j, 1) .le. xlim(1, 1) .and. state(j, 4) < 0.) then
          state(j, 4) = -state(j, 4)
        else if (state(j, 1) .ge. xlim(2, 1) .and. state(j, 4) > 0.) then
          state(j, 4) = -state(j, 4)
        else if (state(j, 2) .le. xlim(1, 2) .and. state(j, 5) < 0.) then
          state(j, 5) = -state(j, 5)
        else if (state(j, 2) .ge. xlim(2, 2) .and. state(j, 5) > 0.) then
          state(j, 5) = -state(j, 5)
        else if (state(j, 3) .le. xlim(1, 3) .and. state(j, 6) < 0.) then
          state(j, 6) = -state(j, 6)
        else if (state(j, 3) .ge. xlim(2, 3) .and. state(j, 6) > 0.) then
          state(j, 6) = -state(j, 6)
        end if
      end do
    end subroutine


    function dRK4(y1, t1, h, f) result (y2)
      real, dimension(:,:),                             intent(in) :: y1
      real,                                             intent(in) :: t1, h
      real, dimension(size(y1,1), size(y1,2))                :: y2
      interface
        function f(t, state) result(dstate)
          real,                                    intent(in) :: t
          real, dimension(:,:),                    intent(in) :: state
          real, dimension(size(state,1), size(state,2)) :: dstate
        end function
      end interface
      !------------------------------------------------
      real                      :: t2
      real, dimension(size(y1,1), size(y1,2)) :: k1, k2, k3, k4
      integer :: i
      ! Runge-Kutta: dy/dt = f(t, y), y(t_0) = y_0
      ! y_{n+1}=y_n + h/6*(k_1 + 2*k_2 + 2*k_3 + k_4)
      ! t_{n+1}=t_n + h
      ! 
      ! k_1 = f(t_n, y_n)
      ! k_2 = f(t_n+h/2, y_n+k_1*h/2)
      ! k_3 = f(t_n+h/2, y_n+k_2*h/2)
      ! k_4 = f(t_n+h, y_n+k_3*h)
      ! t2 = t1 + h
      k1 = f(t1       , y1           )
      k2 = f(t1 + h/2., y1 + k1*h/2.)
      k3 = f(t1 + h/2., y1 + k2*h/2.)
      k4 = f(t1 + h   , y1 + k3*h   )
      y2 = h/6.*(k1 + 2*k2 + 2*k3 + k4)
    end function

    function zero_force(t, state) result(dstate)
      real              , intent(in) :: t
      real, dimension(:,:), intent(in) :: state
      real, dimension(size(state,1), size(state,2))     :: dstate

      dstate(:,1:3) = state(:,4:6)
      dstate(1,4:6) = 0.
    end function zero_force

    function constant_a_1(t, state) result(dstate)
      real              , intent(in) :: t
      real, dimension(:,:), intent(in) :: state
      real, dimension(size(state,1), size(state,2))     :: dstate

      dstate(:,1:3) = state(:,4:6)
      dstate(1,4:6) = [0.,2.,0.]
    end function constant_a_1

    function towards_zero_x(t, state) result(dstate)
      real              , intent(in) :: t
      real, dimension(:,:), intent(in) :: state
      real, dimension(size(state,1), size(state,2))     :: dstate
      real, dimension(1,3) :: x

      x = state(:,1:3)
      dstate(:,1:3) = state(:,4:6)
      dstate(1,4:6) = -x(1,:)
    end function towards_zero_x

    function p1_attracts_p2(t, state) result(dstate)
      real              , intent(in) :: t
      real, dimension(:,:), intent(in) :: state
      real, dimension(size(state,1), size(state,2))     :: dstate
      real, dimension(3) :: x1, x2
      real :: r

      x1 = state(1,1:3)
      x2 = state(2,1:3)
      dstate(:,4:6) = 0.
      r = sqrt(sum((x1-x2)**2.))
      dstate(:,1:3) = state(:,4:6)
      dstate(1,4:6) = (x2-x1)/r**2.
      dstate(2,4:6) = (x1-x2)/r**2.
    end function p1_attracts_p2

    function pi_attracts_pj(t, state) result(dstate)
      real              , intent(in) :: t
      real, dimension(:,:), intent(in) :: state
      real, dimension(size(state,1), size(state,2))     :: dstate
      real, dimension(3) :: x1, x2
      real, dimension(3) :: a
      real :: r
      integer :: i, j

      dstate(:,1:3) = state(:,4:6)
      dstate(:,4:6) = 0.
      do i=1,n_particles
        do j=1,n_particles
          if (i .ne. j) then
            x1 = state(i,1:3)
            x2 = state(j,1:3)
            r = sqrt(sum((x1-x2)**2.))
            a = (x2-x1)/r**2.
            dstate(i,4:6) = dstate(i,4:6) + a
          end if
        end do
      end do
          
    end function pi_attracts_pj

end module fmot
