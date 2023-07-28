  module fmot
    use initialization
    use output
    use integ
    use forces
    implicit none
    integer, parameter :: n_particles = 3

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

      call init_particles(x, v, speed=1., radius=1., n_particles=n_particles)
      x(1,3) = x(1,3)

      state(:,1:3) = x(:,1:3)
      state(:,4:6) = v(:,1:3)

      xlim(1, 1:3) = [-box_size,-box_size,-box_size] ! minimum x
      xlim(2, 1:3) = [box_size,box_size,box_size] ! maximum x

      call write_particles(state, 0, 0., n_particles)
      do i = 1, 100
        t = i*h
        state(:, 1:3) = state(:, 1:3) + dRK4(state(:, 4:6), t, h, velocity_to_delta_x, n_particles)
        state = state &
                + dRK4(state, t, h, pi_attracts_pj, n_particles) !&
!                + dRK4(state, t, h, medium_resistance, n_particles)
        call set_walls(state, xlim)
        if (modulo(i, write_output) .eq. 0) then
          call write_particles(state, i, t, n_particles)
        end if
      end do
    end subroutine main_loop

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


end module fmot
