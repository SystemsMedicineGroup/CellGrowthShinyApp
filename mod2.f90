module params_mod
  implicit none
  double precision :: r, K, kd, kt
end module params_mod

subroutine initcal(parms)
  use params_mod
  implicit none
  double precision, intent(in) :: parms(4)
  r = parms(1)
  K = parms(2)
  kd = parms(3)
  kt = parms(4)
end subroutine initcal


subroutine derivscal(neq, t, y, ydot, out, ip)
  use params_mod
  implicit none
  integer, intent(in) :: neq
  double precision, intent(in) :: t
  double precision, intent(in) :: y(neq)
  double precision, intent(out) :: ydot(neq)
  double precision, intent(out) :: out(*)
  integer, intent(in) :: ip(*)

  if (neq < 3) then
    print *, "Error: neq must be at least 3"
    stop
  endif

  ydot(1) = -kt*y(1)  
  ydot(2) = kt*y(1) + r * y(2) * (1.0d0 - y(2)/K) - kd * y(2)
  ydot(3) = kd * y(2)

  out(1) = y(1)  ! Q
  out(2) = y(2)  ! P
  out(3) = y(3)  ! D
end subroutine derivscal
