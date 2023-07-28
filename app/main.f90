program main
  use fmot, only: main_loop
  use types, only: model_t
  implicit none

  type(model_t) :: model

  model % n_particles = 3
  model % box_size=2.
  model % h = 0.1035
  model % write_output = 1

  call main_loop(model)
end program main
