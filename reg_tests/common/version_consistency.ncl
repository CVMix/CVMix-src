; The isfilepresent() routine is different in version 6.2.1 than it was in
; previous versions - it now only returns true if the file exists AND the file
; is of a valid type for addfile() [netCDF, HDF5, etc]. The new routine
; fileexists() replaces the old functionality... but fileexists() doesn't exist
; in older versions of NCL. So this library contains a wrapper to use
; fileexists() if it is part of NCL and isfilepresent() if it is not.

undef("check_for_file")
function check_for_file(filename)
begin
  if (isfunc("fileexists")) then
    return(fileexists(filename))
  end if
  return(isfilepresent(filename))

end


