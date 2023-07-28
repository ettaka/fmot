module integ
  use types
  implicit none

  contains
    function dRK4(y1, t1, f, model) result (dy)
      real, dimension(:,:), intent(in)        :: y1
      real, intent(in)                        :: t1
      type(model_t)                           :: model
      real, dimension(size(y1,1), size(y1,2)) :: dy
      interface
        function f(t, state, model) result(dstate)
          use types
          real, intent(in)                              :: t
          real, dimension(:,:), intent(in)              :: state
          real, dimension(size(state,1), size(state,2)) :: dstate
          type(model_t)                                 :: model
        end function
      end interface
      !------------------------------------------------
      real                                              :: h
      real, dimension(size(y1,1), size(y1,2))           :: k1, k2, k3, k4
      ! Runge-Kutta: dy/dt = f(t, y), y(t_0) = y_0
      ! y_{n+1}=y_n + h/6*(k_1 + 2*k_2 + 2*k_3 + k_4)
      ! t_{n+1}=t_n + h
      ! 
      ! k_1 = f(t_n, y_n)
      ! k_2 = f(t_n+h/2, y_n+k_1*h/2)
      ! k_3 = f(t_n+h/2, y_n+k_2*h/2)
      ! k_4 = f(t_n+h, y_n+k_3*h)
      ! t2 = t1 + h
      h = model % h
      k1 = f(t1       , y1          , model)
      k2 = f(t1 + h/2., y1 + k1*h/2., model)
      k3 = f(t1 + h/2., y1 + k2*h/2., model)
      k4 = f(t1 + h   , y1 + k3*h   , model)
      dy = h/6.*(k1 + 2*k2 + 2*k3 + k4)
    end function
end module integ
