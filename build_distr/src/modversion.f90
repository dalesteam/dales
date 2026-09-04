! Variables in this file are filled in by CMake at build time,
! used to record the git version and hash.

module modversion

  implicit none

  character(80) :: git_version="v5.0.0-beta.1-247-g939508-dirty"

end module modversion
