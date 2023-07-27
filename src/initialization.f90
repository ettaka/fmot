module initialization
  implicit none
  real, parameter :: pi = 4.*atan(1.)
  contains

    subroutine init_particles(x, v, speed, radius, n_particles)
      real, dimension(n_particles,3) :: x
      real, dimension(n_particles,3) :: v
      real, intent(in) :: speed
      real, intent(in) :: radius
      integer, intent(in) :: n_particles
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

end module initialization
