module forces
  implicit none
  contains

    function zero_force(t, state, n_particles) result(dstate)
      real, intent(in)                              :: t
      real, dimension(:,:), intent(in)              :: state
      integer                                       :: n_particles
      real, dimension(size(state,1), size(state,2)) :: dstate

      dstate(:,1:3) = state(:,4:6)
      dstate(1,4:6) = 0.
    end function zero_force

    function constant_a_1(t, state, n_particles) result(dstate)
      real, intent(in)                              :: t
      real, dimension(:,:), intent(in)              :: state
      integer                                       :: n_particles
      real, dimension(size(state,1), size(state,2)) :: dstate

      dstate(:,1:3) = state(:,4:6)
      dstate(1,4:6) = [0.,2.,0.]
    end function constant_a_1

    function towards_zero_x(t, state, n_particles) result(dstate)
      real, intent(in)                              :: t
      real, dimension(:,:), intent(in)              :: state
      integer                                       :: n_particles
      real, dimension(size(state,1), size(state,2)) :: dstate
      real, dimension(1,3)                          :: x

      x = state(:,1:3)
      dstate(:,1:3) = state(:,4:6)
      dstate(1,4:6) = -x(1,:)
    end function towards_zero_x

    function p1_attracts_p2(t, state, n_particles) result(dstate)
      real, intent(in)                              :: t
      real, dimension(:,:), intent(in)              :: state
      integer                                       :: n_particles
      real, dimension(size(state,1), size(state,2)) :: dstate
      real, dimension(3)                            :: x1, x2
      real                                          :: r

      x1 = state(1,1:3)
      x2 = state(2,1:3)
      dstate(:,4:6) = 0.
      r = sqrt(sum((x1-x2)**2.))
      dstate(:,1:3) = state(:,4:6)
      dstate(1,4:6) = (x2-x1)/r**2.
      dstate(2,4:6) = (x1-x2)/r**2.
    end function p1_attracts_p2

    function velocity_to_delta_x(t, velocity, n_particles) result(dx)
      real, intent(in)                                    :: t
      real, dimension(:,:), intent(in)                    :: velocity
      integer                                             :: n_particles
      real, dimension(size(velocity,1), size(velocity,2)) :: dx

      dx(:,1:3) = velocity(:,1:3)
    end function velocity_to_delta_x

    function pi_attracts_pj(t, state, n_particles) result(dstate)
      real, intent(in)                              :: t
      real, dimension(:,:), intent(in)              :: state
      integer                                       :: n_particles
      real, dimension(size(state,1), size(state,2)) :: dstate
      real, dimension(3)                            :: x1, x2
      real, dimension(3)                            :: a
      real                                          :: r
      integer                                       :: i, j

      dstate(:,:) = 0.
      do i=1,n_particles
        do j=1,n_particles
          if (i .ne. j) then
            x1 = state(i,1:3)
            x2 = state(j,1:3)
            r = sqrt(sum((x1-x2)**2.))
            a = (x2-x1)/r**2.
            if (r<0.9) a=-2.*a
            dstate(i,4:6) = dstate(i,4:6) + a
          end if
        end do
      end do
          
    end function pi_attracts_pj

    function medium_resistance(t, state, n_particles) result(dstate)
      real, intent(in)                              :: t
      real, dimension(:,:), intent(in)              :: state
      integer                                       :: n_particles
      real, dimension(size(state,1), size(state,2)) :: dstate
      real, dimension(3)                            :: x1, x2
      real, dimension(3)                            :: a, v
      real                                          :: r
      integer                                       :: i

      dstate(:,:) = 0.
      do i=1,n_particles
        v = state(i,4:6)
        a = -10*v**2. * v/sqrt(sum(v**2.))
        dstate(i,4:6) = dstate(i,4:6) + a
      end do
          
    end function medium_resistance
end module forces
