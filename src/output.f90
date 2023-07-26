module output
  implicit none

  contains

    subroutine write_particles(state, i, t, n_particles)
      real, dimension(n_particles,6), intent(in) :: state
      integer, intent(in) :: i
      real, intent(in) :: t
      integer, intent(in) :: n_particles
      integer :: j
      character(len=12) :: filename

      write(filename, '(A4,I4.4,A4)') "res_",i,".csv"
      open(unit=20, file=filename)

      write(20, '(A16)') "x1, x2, x3, t" 
      do j = 1, n_particles
        write(20, '(4(F8.4,A2))') state(j,1), ",", state(j,2), ",", state(j,3), ",", t
      end do

      close(20)
    end subroutine

end module output
