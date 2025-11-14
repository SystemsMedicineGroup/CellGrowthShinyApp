module params_mod
  implicit none
  double precision :: r, K, kd
end module params_mod

subroutine initcal(parms)
  use params_mod
  implicit none
  double precision, intent(in) :: parms(3)
  r = parms(1)
  K = parms(2)
  kd = parms(3)
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

  if (neq < 2) then
    print *, "Error: neq must be at least 3"
    stop
  endif

  ydot(1) = r * y(1) * (1.0d0 - y(1)/K) - kd * y(1)
  ydot(2) = kd * y(1)

  out(1) = y(1)  ! P
  out(2) = y(2)  ! D
end subroutine derivscal
