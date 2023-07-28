module types
  implicit none

  type :: point_mass_t
    real, dimension(3) :: coordinates
  end type point_mass_t

  type :: model_t
    integer               :: n_particles = 3
    real                  :: box_size=2.
    real                  :: h = 0.1035
    integer               :: write_output = 1
  end type model_t
end module types
