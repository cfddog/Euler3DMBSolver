function vminmod(a1,a2)
real vminmod
real a1,a2
 if(a1*a2.le.0) then
       vminmod=0.
else
      if(abs(a1) .le. abs(a2)) then
       vminmod=a1
      else
       vminmod=a2
      endif
 endif
 return
 end function 