  module fmot
    use initialization
    use output
    use integ
    use forces
    use types
    implicit none

  contains

    subroutine main_loop(model)
      type(model_t), intent(in)              :: model
      integer                                :: i, j
      real, dimension(model % n_particles,3) :: x
      real, dimension(model % n_particles,3) :: v
      real, dimension(model % n_particles,6) :: state
      real, dimension(2, 3)                  :: xlim
      real                                   :: t = 0.

      call init_particles(x, v, speed=1., radius=1., n_particles=model % n_particles)

      state(:,1:3) = x(:,1:3)
      state(:,4:6) = v(:,1:3)

      xlim(1, 1:3) = [-model % box_size,-model % box_size,-model % box_size] ! minimum x
      xlim(2, 1:3) = [ model % box_size, model % box_size, model % box_size] ! maximum x

      call write_particles(state, 0, 0., model % n_particles)
      do i = 1, 100
        t = i*model % h
        state(:, 1:3) = state(:, 1:3) + dRK4(state(:, 4:6), t, velocity_to_delta_x, model)
        state = state &
                + dRK4(state, t, pi_attracts_pj, model) !&
                !+ dRK4(state, t, medium_resistance, model)
        call set_walls(state, xlim)
        if (modulo(i, model % write_output) .eq. 0) then
          call write_particles(state, i, t, model % n_particles)
        end if
      end do
    end subroutine main_loop

    subroutine set_walls(state, xlim)
      real, dimension(:,:) :: state
      real, dimension(2, 3), intent(in) :: xlim
      integer :: j

      do j = 1, size(state, dim=1)
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
