program main
  use fmot, only: main_loop
  use types, only: model_t
  implicit none

  type(model_t) :: model

  model % n_particles = 3

  call main_loop(model)
end program main
