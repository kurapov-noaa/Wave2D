subroutine write_his(hisname,u,nx,ny,t)

  use netcdf
  implicit none

  integer, intent(in) :: Nx, Ny
  ! Nx+1 and N+1 are provided as inputs for variables that start from 0
  real, intent(in) :: u(nx,ny)
  real, intent(in) :: t                  ! output time

  character (len=*), intent(in) :: hisname

  integer :: ncid
  integer :: dimid_time
  integer :: vid_time, vid_u
  integer :: nt, rec
  integer :: status

  !----------------------------------------------------
  ! Open file for writing
  !----------------------------------------------------
  status = nf90_open(trim(hisname), nf90_write, ncid)
  if (status /= nf90_noerr) call abort_nc(status)

  !----------------------------------------------------
  ! Determine append record index
  !----------------------------------------------------
  status = nf90_inq_dimid(ncid, "ocean_time", dimid_time)
  if (status /= nf90_noerr) call abort_nc(status)

  status = nf90_inquire_dimension(ncid, dimid_time, len=nt)
  if (status /= nf90_noerr) call abort_nc(status)

  rec = nt + 1     ! <-- append position

  !----------------------------------------------------
  ! Variable IDs
  !----------------------------------------------------
  status = nf90_inq_varid(ncid, "ocean_time", vid_time)
  if (status /= nf90_noerr) call abort_nc(status)

  status = nf90_inq_varid(ncid, "u", vid_u)
  if (status /= nf90_noerr) call abort_nc(status)

  !----------------------------------------------------
  ! Append ocean_time
  !----------------------------------------------------
  status = nf90_put_var(ncid, vid_time, t, start=(/rec/))
  if (status /= nf90_noerr) call abort_nc(status)

  !----------------------------------------------------
  ! Append fields at same record
  !----------------------------------------------------
  status = nf90_put_var(ncid, vid_u, u, start=(/1,1,rec/))
  if (status /= nf90_noerr) call abort_nc(status)

  !----------------------------------------------------
  ! Close file
  !----------------------------------------------------
  status = nf90_close(ncid)
  if (status /= nf90_noerr) call abort_nc(status)

contains

  subroutine abort_nc(status)
    integer, intent(in) :: status
    write(*,*) "NetCDF error: ", trim(nf90_strerror(status))
    stop
  end subroutine abort_nc

end subroutine write_his
