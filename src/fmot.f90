  module fmot
    implicit none
    integer, parameter :: n_particles = 3

  contains
    subroutine say_hello
      use types
      implicit none 

      type(point_mass) :: test_mass

      test_mass % coordinates(1) = 1.
      test_mass % coordinates(2) = 2.
      test_mass % coordinates(3) = 3.

      print *, "test_mass % coordinates", test_mass % coordinates
    end subroutine say_hello

    subroutine main_loop
      use types
      implicit none 

      integer :: i, j
      character(len=12) :: filename

      real, dimension(n_particles,3) :: x
      real, dimension(n_particles,3) :: v
      real, dimension(n_particles,3) :: a
      real, dimension(n_particles,3) :: f

      real, dimension(n_particles,6) :: state

      real :: t, h

      ! f = ma = m dv/dt = m ddt(x)
      ! Runge-Kutta: dy/dt = g(t, y), g(t_0) = y_0
      ! g = f/m
      ! 

    x = 0.
    x(1,:) = [1./sqrt(3.),1./2./sqrt(3.),0.]
    x(2,:) = -x(1,:)
    x(3,:) = [0.,-1.,0.]
    v = 0.
    v(1,:) = [-x(1,1), x(1,2), 0.]
    v(2,:) = [-x(1,1), -x(1,2), 0.]
    v(3,:) = [1.,0.,0.]
    a = 0.
    t = 0.
    h = 0.2

    state(:,1:3) = x(:,1:3)
    state(:,4:6) = v(:,1:3)

    do i = 1, 100
      t = i*h
      state = RK4(state, t, h, pi_attracts_pj)
      write(filename, '(A4,I4.4,A4)') "res_",i,".csv"
      open(unit=20, file=filename)

      write(20, '(A16)') "x1, x2, x3, t" 
      do j = 1, n_particles
        write(20, '(4(F8.4,A2))') state(j,1), ",", state(j,2), ",", state(j,3), ",", t
      end do

      close(20)
    end do
  end subroutine main_loop

  function RK4(y1, t1, h, f) result (y2)
    implicit none
    real, dimension(:,:), intent(in) :: y1
    real,                 intent(in) :: t1, h
    interface
        function f(t, x) result(y)
            implicit none
            real                , intent(in) :: t
            real, dimension(:,:), intent(in) :: x
            real, dimension(size(x(:,1)), size(x(1,:))) :: y
        end function f
    end interface
    real, dimension(size(y1(:,1)), size(y1(1,:))) :: y2
    !------------------------------------------------
    real                      :: t2
    real, dimension(size(y1(:,1)), size(y1(1,:))) :: k1, k2, k3, k4
      ! Runge-Kutta: dy/dt = f(t, y), y(t_0) = y_0
      ! y_{n+1}=y_n + h/6*(k_1 + 2*k_2 + 2*k_3 + k_4)
      ! t_{n+1}=t_n + h
      ! 
      ! k_1 = f(t_n, y_n)
      ! k_2 = f(t_n+h/2, y_n+k_1*h/2)
      ! k_3 = f(t_n+h/2, y_n+k_2*h/2)
      ! k_4 = f(t_n+h, y_n+k_3*h)
    t2 = t1 + h
    k1 = f(t1       , y1           )
    k2 = f(t1 + h/2., y1 + k1*h/2.)
    k3 = f(t1 + h/2., y1 + k2*h/2.)
    k4 = f(t1 + h   , y1 + k3*h   )
    y2 = y1 + h/6.*(k1 + 2*k2 + 2*k3 + k4)
  end function

  function constant_a_1(t, state) result(dstate)
    implicit none
    real              , intent(in) :: t
    real, dimension(:,:), intent(in) :: state
    real, dimension(size(state(:,1)), size(state(1,:)))     :: dstate

    dstate(:,1:3) = state(:,4:6)
    dstate(1,4:6) = [1.0,0.1,0.01]
  end function constant_a_1

  function towards_zero_x(t, state) result(dstate)
    implicit none
    real              , intent(in) :: t
    real, dimension(:,:), intent(in) :: state
    real, dimension(size(state(:,1)), size(state(1,:)))     :: dstate
    real, dimension(1,3) :: x

    x = state(:,1:3)
    dstate(:,1:3) = state(:,4:6)
    dstate(1,4:6) = -x(1,:)
  end function towards_zero_x

  function p1_attracts_p2(t, state) result(dstate)
    implicit none
    real              , intent(in) :: t
    real, dimension(:,:), intent(in) :: state
    real, dimension(size(state(:,1)), size(state(1,:)))     :: dstate
    real, dimension(3) :: x1, x2
    real :: r

    x1 = state(1,1:3)
    x2 = state(2,1:3)
    r = sqrt(sum((x1-x2)**2.))
    dstate(:,1:3) = state(:,4:6)
    dstate(1,4:6) = (x2-x1)/r**2.
    dstate(2,4:6) = (x1-x2)/r**2.
  end function p1_attracts_p2

  function pi_attracts_pj(t, state) result(dstate)
    implicit none
    real              , intent(in) :: t
    real, dimension(:,:), intent(in) :: state
    real, dimension(size(state(:,1)), size(state(1,:)))     :: dstate
    real, dimension(3) :: x1, x2
    real :: r
    integer :: i, j

    dstate(:,1:3) = state(:,4:6)
    do i=1,n_particles
      do j=1,n_particles
        if (i .ne. j) then
          x1 = state(i,1:3)
          x2 = state(j,1:3)
          r = sqrt(sum((x1-x2)**2.))
          dstate(i,4:6) = (x2-x1)/r**2.
        end if
      end do
    end do
        
  end function pi_attracts_pj
end module fmot
