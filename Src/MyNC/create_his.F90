subroutine create_his(hisname,Nx,Ny)
  use netcdf
  implicit none

  integer, intent(in) :: nx, ny
  character(len=*), intent(in) :: hisname  ! <--- Assumed-length string

  integer :: ncid
  integer :: dim_xi, dim_eta
  integer :: dim_time

  integer :: var_time, var_u
  integer :: status

  character(len=8)  :: tmp_str
  character(len=30) :: att_value

  !----------------------------------------------------
  ! Create file
  !----------------------------------------------------
  status = nf90_create(trim(hisname), nf90_clobber, ncid)
  if (status /= nf90_noerr) call abort_nc(status)

  !----------------------------------------------------
  ! Dimensions
  !----------------------------------------------------
  status = nf90_def_dim(ncid, "xi_u", nx+1,  dim_xi)
  if (status /= nf90_noerr) call abort_nc(status)

  status = nf90_def_dim(ncid, "eta_u",  ny+1, dim_eta)
  if (status /= nf90_noerr) call abort_nc(status)

  status = nf90_def_dim(ncid, "ocean_time", nf90_unlimited, dim_time)
  if (status /= nf90_noerr) call abort_nc(status)

  !----------------------------------------------------
  ! Variables
  !----------------------------------------------------
  status = nf90_def_var(ncid, "ocean_time", nf90_real, &
                        (/ dim_time /), var_time)
  if (status /= nf90_noerr) call abort_nc(status)

  status = nf90_def_var(ncid, "u", nf90_real, &
                        (/ dim_xi, dim_eta, dim_time /), var_u)
  if (status /= nf90_noerr) call abort_nc(status)

  !----------------------------------------------------
  ! Attributes
  !----------------------------------------------------

  ! ocean_time attributes:
  status = nf90_put_att(ncid, var_time, "long_name", "time since start")
  if (status /= nf90_noerr) call abort_nc(status)

  status = nf90_put_att(ncid, var_time, "units", "seconds")
  if (status /= nf90_noerr) call abort_nc(status)


  ! others:
  status = nf90_put_att(ncid, var_u, "long_name", "u")
  if (status /= nf90_noerr) call abort_nc(status)


  !----------------------------------------------------
  ! End define mode and close
  !----------------------------------------------------
  status = nf90_enddef(ncid)
  if (status /= nf90_noerr) call abort_nc(status)

  status = nf90_close(ncid)
  if (status /= nf90_noerr) call abort_nc(status)

contains

  subroutine abort_nc(status)
    integer, intent(in) :: status
    write(*,*) "NetCDF error: ", trim(nf90_strerror(status))
    stop
  end subroutine abort_nc

end subroutine create_his
