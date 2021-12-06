undef("get1Dindex2")
function get1Dindex2 (cv[*],  cvWant[*])
local nWant, indWant, n, kn, kind, nkk
begin
  ncv     = dimsizes (cv)
  nWant   = dimsizes (cvWant)
  indWant = new (nWant*ncv, integer)
  kn = -1
  do n=0,nWant-1
     kind = ind( cv.eq.cvWant(n) )
     if (.not.ismissing(kind(0))) then
         nkk = dimsizes(kind)
         do km=0,nkk-1
            kn = kn+1
            indWant(kn) = kind(km)
         end do
         delete( nkk )
      end if
      delete(kind)
  end do

  return (indWant(0:kn))        
end

