program test_types
  use types
  implicit none

  type(point_mass) :: test_mass

  test_mass % coordinates(1) = 1.
  test_mass % coordinates(2) = 2.
  test_mass % coordinates(3) = 3.

  print *, "test_mass % coordinates", test_mass % coordinates
end program test_types
