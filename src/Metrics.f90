!compute the metrics for all zones
! at the first stage, just for one zone,use visbal's method or direct method 
subroutine ComputeMetrics()
use Ctrlpm
use fieldpm
type(BLOCK_TYPE),pointer :: pBlk
real xxc,xet,xct,yxc,yet,yct,zxc,zet,zct
!now this is the non-SCMM method. this will be to do later. !todo: SCMM
do iblk=1,numblk
   pBlk=>compblock(iblk)
do i=pBlk%icmpst,pBlk%icmped
  do j=pBlk%jcmpst,pBlk%jcmped
    do k=pBlk%kcmpst,pBlk%kcmped
! xc direction metrics    
       if(i.eq.pBlk%ist) then
         xxc=pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i,j,k)
         yxc=pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i,j,k)
         zxc=pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i,j,k) 
       elseif(i.eq.pBlk%ied) then
         xxc=pBlk%xcoord(i,j,k)-pBlk%xcoord(i-1,j,k)
         yxc=pBlk%ycoord(i,j,k)-pBlk%ycoord(i-1,j,k)
         zxc=pBlk%zcoord(i,j,k)-pBlk%zcoord(i-1,j,k)
       elseif(i.eq.pBlk%ist+1 .or. i.eq.pBlk%ied-1) then
         xxc=0.5*(pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i-1,j,k))
         yxc=0.5*(pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i-1,j,k))
         zxc=0.5*(pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i-1,j,k))
        else
         xxc=0.5*(pBlk%xcoord(i+1,j,k) - pBlk%xcoord(i-1,j,k))*4./3  - 0.5*(pBlk%xcoord(i+2,j,k) - pBlk%xcoord(i-2,j,k))/6.
         yxc=0.5*(pBlk%ycoord(i+1,j,k) - pBlk%ycoord(i-1,j,k))*4./3. - 0.5*(pBlk%ycoord(i+2,j,k) - pBlk%ycoord(i-2,j,k))/6.
         zxc=0.5*(pBlk%zcoord(i+1,j,k) - pBlk%zcoord(i-1,j,k))*4./3. - 0.5*(pBlk%zcoord(i+2,j,k) - pBlk%zcoord(i-2,j,k))/6.         
       endif   
 ! et direction metrics
        if(j.eq.pBlk%jst) then
         xet=pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j,k)
         yet=pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j,k)
         zet=pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j,k) 
       elseif(j.eq.pBlk%jed) then
         xet=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j-1,k)
         yet=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j-1,k)
         zet=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j-1,k)
       elseif(j.eq.pBlk%jst+1 .or. j.eq.pBlk%jed-1) then
         xet=0.5*(pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j-1,k))
         yet=0.5*(pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j-1,k))
         zet=0.5*(pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j-1,k))
        else
         xet=0.5*(pBlk%xcoord(i,j+1,k) - pBlk%xcoord(i,j-1,k))*4./3  - 0.5*(pBlk%xcoord(i,j+2,k) - pBlk%xcoord(i,j-2,k))/6.
         yet=0.5*(pBlk%ycoord(i,j+1,k) - pBlk%ycoord(i,j-1,k))*4./3. - 0.5*(pBlk%ycoord(i,j+2,k) - pBlk%ycoord(i,j-2,k))/6.
         zet=0.5*(pBlk%zcoord(i,j+1,k) - pBlk%zcoord(i,j-1,k))*4./3. - 0.5*(pBlk%zcoord(i,j+2,k) - pBlk%zcoord(i,j-2,k))/6.         
       endif 
  ! ct diretion
       if(k.eq.pBlk%kst) then
         xct=pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k)
         yct=pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k)
         zct=pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k) 
       elseif(k.eq.pBlk%ked) then
         xct=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j,k-1)
         yct=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j,k-1)
         zct=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j,k-1)
       elseif(k.eq.pBlk%kst+1 .or. k.eq.pBlk%ked-1) then
         xct=0.5*(pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k-1))
         yct=0.5*(pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k-1))
         zct=0.5*(pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k-1))
        else
         xct=0.5*(pBlk%xcoord(i,j,k+1) - pBlk%xcoord(i,j,k-1))*4./3  - 0.5*(pBlk%xcoord(i,j,k+2) - pBlk%xcoord(i,j,k-2))/6.
         yct=0.5*(pBlk%ycoord(i,j,k+1) - pBlk%ycoord(i,j,k-1))*4./3. - 0.5*(pBlk%ycoord(i,j,k+2) - pBlk%ycoord(i,j,k-2))/6.
         zct=0.5*(pBlk%zcoord(i,j,k+1) - pBlk%zcoord(i,j,k-1))*4./3. - 0.5*(pBlk%zcoord(i,j,k+2) - pBlk%zcoord(i,j,k-2))/6.         
       endif  
    ! compute the metrics
         pBlk%xcx(i,j,k)=yet*zct-yct*zet
         pBlk%xcy(i,j,k)=zet*xct-zct*xet 
         pBlk%xcz(i,j,k)=xet*yct-xct*yet
         
         pBlk%etx(i,j,k)=yct*zxc-yxc*zct
         pBlk%ety(i,j,k)=zct*xxc-zxc*xct 
         pBlk%etz(i,j,k)=xct*yxc-xxc*yct
         
         pBlk%ctx(i,j,k)=yxc*zet-yet*zxc
         pBlk%cty(i,j,k)=zxc*xet-zet*xxc 
         pBlk%ctz(i,j,k)=xxc*yet-xet*yxc
         pBlk%dj(i,j,k)=xxc*(yet*zct-yct*zet)-xet*(yxc*zct-yct*zxc)+xct*(yxc*zet-yet*zxc)
         if( pBlk%dj(i,j,k) .lt. 0) then 
           print*, "the djac less than zero,check the grid please!"
           stop
         elseif( pBlk%dj(i,j,k) .eq. 0) then
           print*," the djac equals zero, specify it 1.e-18!"
           pBlk%dj(i,j,k)=1.e-20
         endif          
     enddo
   enddo
enddo 
!//i=icmpst buffer
  do i=pBlk%ist,pBlk%icmpst-1
    do j=pBlk%jcmpst,pBlk%jcmped
      do k=pBlk%kcmpst,pBlk%kcmped
  ! xc direction metrics    
         if(i.eq.pBlk%ist) then
           xxc=pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i,j,k)
           yxc=pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i,j,k)
           zxc=pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i,j,k) 
         elseif(i.eq.pBlk%ied) then
           xxc=pBlk%xcoord(i,j,k)-pBlk%xcoord(i-1,j,k)
           yxc=pBlk%ycoord(i,j,k)-pBlk%ycoord(i-1,j,k)
           zxc=pBlk%zcoord(i,j,k)-pBlk%zcoord(i-1,j,k)
         elseif(i.eq.pBlk%ist+1 .or. i.eq.pBlk%ied-1) then
           xxc=0.5*(pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i-1,j,k))
           yxc=0.5*(pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i-1,j,k))
           zxc=0.5*(pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i-1,j,k))
          else
           xxc=0.5*(pBlk%xcoord(i+1,j,k) - pBlk%xcoord(i-1,j,k))*4./3  - 0.5*(pBlk%xcoord(i+2,j,k) - pBlk%xcoord(i-2,j,k))/6.
           yxc=0.5*(pBlk%ycoord(i+1,j,k) - pBlk%ycoord(i-1,j,k))*4./3. - 0.5*(pBlk%ycoord(i+2,j,k) - pBlk%ycoord(i-2,j,k))/6.
           zxc=0.5*(pBlk%zcoord(i+1,j,k) - pBlk%zcoord(i-1,j,k))*4./3. - 0.5*(pBlk%zcoord(i+2,j,k) - pBlk%zcoord(i-2,j,k))/6.         
         endif   
   ! et direction metrics
          if(j.eq.pBlk%jcmpst) then
           xet=pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j,k)
           yet=pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j,k)
           zet=pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j,k) 
         elseif(j.eq.pBlk%jcmped) then
           xet=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j-1,k)
           yet=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j-1,k)
           zet=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j-1,k)
         elseif(j.eq.pBlk%jcmpst+1 .or. j.eq.pBlk%jcmped-1) then
           xet=0.5*(pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j-1,k))
           yet=0.5*(pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j-1,k))
           zet=0.5*(pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j-1,k))
          else
           xet=0.5*(pBlk%xcoord(i,j+1,k) - pBlk%xcoord(i,j-1,k))*4./3  - 0.5*(pBlk%xcoord(i,j+2,k) - pBlk%xcoord(i,j-2,k))/6.
           yet=0.5*(pBlk%ycoord(i,j+1,k) - pBlk%ycoord(i,j-1,k))*4./3. - 0.5*(pBlk%ycoord(i,j+2,k) - pBlk%ycoord(i,j-2,k))/6.
           zet=0.5*(pBlk%zcoord(i,j+1,k) - pBlk%zcoord(i,j-1,k))*4./3. - 0.5*(pBlk%zcoord(i,j+2,k) - pBlk%zcoord(i,j-2,k))/6.         
         endif 
    ! ct diretion
         if(k.eq.pBlk%kcmpst) then
           xct=pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k)
           yct=pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k)
           zct=pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k) 
         elseif(k.eq.pBlk%kcmped) then
           xct=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j,k-1)
           yct=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j,k-1)
           zct=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j,k-1)
         elseif(k.eq.pBlk%kcmpst+1 .or. k.eq.pBlk%kcmped-1) then
           xct=0.5*(pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k-1))
           yct=0.5*(pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k-1))
           zct=0.5*(pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k-1))
          else
           xct=0.5*(pBlk%xcoord(i,j,k+1) - pBlk%xcoord(i,j,k-1))*4./3  - 0.5*(pBlk%xcoord(i,j,k+2) - pBlk%xcoord(i,j,k-2))/6.
           yct=0.5*(pBlk%ycoord(i,j,k+1) - pBlk%ycoord(i,j,k-1))*4./3. - 0.5*(pBlk%ycoord(i,j,k+2) - pBlk%ycoord(i,j,k-2))/6.
           zct=0.5*(pBlk%zcoord(i,j,k+1) - pBlk%zcoord(i,j,k-1))*4./3. - 0.5*(pBlk%zcoord(i,j,k+2) - pBlk%zcoord(i,j,k-2))/6.         
         endif  
      ! compute the metrics
           pBlk%xcx(i,j,k)=yet*zct-yct*zet
           pBlk%xcy(i,j,k)=zet*xct-zct*xet 
           pBlk%xcz(i,j,k)=xet*yct-xct*yet
           
           pBlk%etx(i,j,k)=yct*zxc-yxc*zct
           pBlk%ety(i,j,k)=zct*xxc-zxc*xct 
           pBlk%etz(i,j,k)=xct*yxc-xxc*yct
           
           pBlk%ctx(i,j,k)=yxc*zet-yet*zxc
           pBlk%cty(i,j,k)=zxc*xet-zet*xxc 
           pBlk%ctz(i,j,k)=xxc*yet-xet*yxc
           pBlk%dj(i,j,k)=xxc*(yet*zct-yct*zet)-xet*(yxc*zct-yct*zxc)+xct*(yxc*zet-yet*zxc)
           if( pBlk%dj(i,j,k) .lt. 0) then 
             print*, "the djac less than zero,check the grid please!"
             stop
           elseif( pBlk%dj(i,j,k) .eq. 0) then
             print*," the djac equals zero, specify it 1.e-18!"
             pBlk%dj(i,j,k)=1.e-20
           endif          
       enddo
     enddo
  enddo
  !icmped buffer
  do i=pBlk%icmped+1,pBlk%ied
    do j=pBlk%jcmpst,pBlk%jcmped
      do k=pBlk%kcmpst,pBlk%kcmped
  ! xc direction metrics    
         if(i.eq.pBlk%ist) then
           xxc=pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i,j,k)
           yxc=pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i,j,k)
           zxc=pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i,j,k) 
         elseif(i.eq.pBlk%ied) then
           xxc=pBlk%xcoord(i,j,k)-pBlk%xcoord(i-1,j,k)
           yxc=pBlk%ycoord(i,j,k)-pBlk%ycoord(i-1,j,k)
           zxc=pBlk%zcoord(i,j,k)-pBlk%zcoord(i-1,j,k)
         elseif(i.eq.pBlk%ist+1 .or. i.eq.pBlk%ied-1) then
           xxc=0.5*(pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i-1,j,k))
           yxc=0.5*(pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i-1,j,k))
           zxc=0.5*(pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i-1,j,k))
          else
           xxc=0.5*(pBlk%xcoord(i+1,j,k) - pBlk%xcoord(i-1,j,k))*4./3  - 0.5*(pBlk%xcoord(i+2,j,k) - pBlk%xcoord(i-2,j,k))/6.
           yxc=0.5*(pBlk%ycoord(i+1,j,k) - pBlk%ycoord(i-1,j,k))*4./3. - 0.5*(pBlk%ycoord(i+2,j,k) - pBlk%ycoord(i-2,j,k))/6.
           zxc=0.5*(pBlk%zcoord(i+1,j,k) - pBlk%zcoord(i-1,j,k))*4./3. - 0.5*(pBlk%zcoord(i+2,j,k) - pBlk%zcoord(i-2,j,k))/6.         
         endif   
   ! et direction metrics
          if(j.eq.pBlk%jcmpst) then
           xet=pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j,k)
           yet=pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j,k)
           zet=pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j,k) 
         elseif(j.eq.pBlk%jcmped) then
           xet=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j-1,k)
           yet=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j-1,k)
           zet=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j-1,k)
         elseif(j.eq.pBlk%jcmpst+1 .or. j.eq.pBlk%jcmped-1) then
           xet=0.5*(pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j-1,k))
           yet=0.5*(pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j-1,k))
           zet=0.5*(pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j-1,k))
          else
           xet=0.5*(pBlk%xcoord(i,j+1,k) - pBlk%xcoord(i,j-1,k))*4./3  - 0.5*(pBlk%xcoord(i,j+2,k) - pBlk%xcoord(i,j-2,k))/6.
           yet=0.5*(pBlk%ycoord(i,j+1,k) - pBlk%ycoord(i,j-1,k))*4./3. - 0.5*(pBlk%ycoord(i,j+2,k) - pBlk%ycoord(i,j-2,k))/6.
           zet=0.5*(pBlk%zcoord(i,j+1,k) - pBlk%zcoord(i,j-1,k))*4./3. - 0.5*(pBlk%zcoord(i,j+2,k) - pBlk%zcoord(i,j-2,k))/6.         
         endif 
    ! ct diretion
         if(k.eq.pBlk%kcmpst) then
           xct=pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k)
           yct=pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k)
           zct=pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k) 
         elseif(k.eq.pBlk%kcmped) then
           xct=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j,k-1)
           yct=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j,k-1)
           zct=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j,k-1)
         elseif(k.eq.pBlk%kcmpst+1 .or. k.eq.pBlk%kcmped-1) then
           xct=0.5*(pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k-1))
           yct=0.5*(pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k-1))
           zct=0.5*(pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k-1))
          else
           xct=0.5*(pBlk%xcoord(i,j,k+1) - pBlk%xcoord(i,j,k-1))*4./3  - 0.5*(pBlk%xcoord(i,j,k+2) - pBlk%xcoord(i,j,k-2))/6.
           yct=0.5*(pBlk%ycoord(i,j,k+1) - pBlk%ycoord(i,j,k-1))*4./3. - 0.5*(pBlk%ycoord(i,j,k+2) - pBlk%ycoord(i,j,k-2))/6.
           zct=0.5*(pBlk%zcoord(i,j,k+1) - pBlk%zcoord(i,j,k-1))*4./3. - 0.5*(pBlk%zcoord(i,j,k+2) - pBlk%zcoord(i,j,k-2))/6.         
         endif  
      ! compute the metrics
           pBlk%xcx(i,j,k)=yet*zct-yct*zet
           pBlk%xcy(i,j,k)=zet*xct-zct*xet 
           pBlk%xcz(i,j,k)=xet*yct-xct*yet
           
           pBlk%etx(i,j,k)=yct*zxc-yxc*zct
           pBlk%ety(i,j,k)=zct*xxc-zxc*xct 
           pBlk%etz(i,j,k)=xct*yxc-xxc*yct
           
           pBlk%ctx(i,j,k)=yxc*zet-yet*zxc
           pBlk%cty(i,j,k)=zxc*xet-zet*xxc 
           pBlk%ctz(i,j,k)=xxc*yet-xet*yxc
           pBlk%dj(i,j,k)=xxc*(yet*zct-yct*zet)-xet*(yxc*zct-yct*zxc)+xct*(yxc*zet-yet*zxc)
           if( pBlk%dj(i,j,k) .lt. 0) then 
             print*, "the djac less than zero,check the grid please!"
             stop
           elseif( pBlk%dj(i,j,k) .eq. 0) then
             print*," the djac equals zero, specify it 1.e-18!"
             pBlk%dj(i,j,k)=1.e-20
           endif          
       enddo
     enddo
  enddo  
  !jcmpst buffer
  do i=pBlk%icmpst,pBlk%icmped
    do j=pBlk%jst,pBlk%jcmpst-1
      do k=pBlk%kcmpst,pBlk%kcmped
  ! xc direction metrics    
         if(i.eq.pBlk%icmpst) then
           xxc=pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i,j,k)
           yxc=pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i,j,k)
           zxc=pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i,j,k) 
         elseif(i.eq.pBlk%icmped) then
           xxc=pBlk%xcoord(i,j,k)-pBlk%xcoord(i-1,j,k)
           yxc=pBlk%ycoord(i,j,k)-pBlk%ycoord(i-1,j,k)
           zxc=pBlk%zcoord(i,j,k)-pBlk%zcoord(i-1,j,k)
         elseif(i.eq.pBlk%icmpst+1 .or. i.eq.pBlk%icmped-1) then
           xxc=0.5*(pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i-1,j,k))
           yxc=0.5*(pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i-1,j,k))
           zxc=0.5*(pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i-1,j,k))
          else
           xxc=0.5*(pBlk%xcoord(i+1,j,k) - pBlk%xcoord(i-1,j,k))*4./3  - 0.5*(pBlk%xcoord(i+2,j,k) - pBlk%xcoord(i-2,j,k))/6.
           yxc=0.5*(pBlk%ycoord(i+1,j,k) - pBlk%ycoord(i-1,j,k))*4./3. - 0.5*(pBlk%ycoord(i+2,j,k) - pBlk%ycoord(i-2,j,k))/6.
           zxc=0.5*(pBlk%zcoord(i+1,j,k) - pBlk%zcoord(i-1,j,k))*4./3. - 0.5*(pBlk%zcoord(i+2,j,k) - pBlk%zcoord(i-2,j,k))/6.         
         endif   
   ! et direction metrics
          if(j.eq.pBlk%jst) then
           xet=pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j,k)
           yet=pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j,k)
           zet=pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j,k) 
         elseif(j.eq.pBlk%jed) then
           xet=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j-1,k)
           yet=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j-1,k)
           zet=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j-1,k)
         elseif(j.eq.pBlk%jst+1 .or. j.eq.pBlk%jed-1) then
           xet=0.5*(pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j-1,k))
           yet=0.5*(pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j-1,k))
           zet=0.5*(pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j-1,k))
          else
           xet=0.5*(pBlk%xcoord(i,j+1,k) - pBlk%xcoord(i,j-1,k))*4./3  - 0.5*(pBlk%xcoord(i,j+2,k) - pBlk%xcoord(i,j-2,k))/6.
           yet=0.5*(pBlk%ycoord(i,j+1,k) - pBlk%ycoord(i,j-1,k))*4./3. - 0.5*(pBlk%ycoord(i,j+2,k) - pBlk%ycoord(i,j-2,k))/6.
           zet=0.5*(pBlk%zcoord(i,j+1,k) - pBlk%zcoord(i,j-1,k))*4./3. - 0.5*(pBlk%zcoord(i,j+2,k) - pBlk%zcoord(i,j-2,k))/6.         
         endif 
    ! ct diretion
         if(k.eq.pBlk%kcmpst) then
           xct=pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k)
           yct=pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k)
           zct=pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k) 
         elseif(k.eq.pBlk%kcmped) then
           xct=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j,k-1)
           yct=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j,k-1)
           zct=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j,k-1)
         elseif(k.eq.pBlk%kcmpst+1 .or. k.eq.pBlk%kcmped-1) then
           xct=0.5*(pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k-1))
           yct=0.5*(pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k-1))
           zct=0.5*(pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k-1))
          else
           xct=0.5*(pBlk%xcoord(i,j,k+1) - pBlk%xcoord(i,j,k-1))*4./3  - 0.5*(pBlk%xcoord(i,j,k+2) - pBlk%xcoord(i,j,k-2))/6.
           yct=0.5*(pBlk%ycoord(i,j,k+1) - pBlk%ycoord(i,j,k-1))*4./3. - 0.5*(pBlk%ycoord(i,j,k+2) - pBlk%ycoord(i,j,k-2))/6.
           zct=0.5*(pBlk%zcoord(i,j,k+1) - pBlk%zcoord(i,j,k-1))*4./3. - 0.5*(pBlk%zcoord(i,j,k+2) - pBlk%zcoord(i,j,k-2))/6.         
         endif  
      ! compute the metrics
           pBlk%xcx(i,j,k)=yet*zct-yct*zet
           pBlk%xcy(i,j,k)=zet*xct-zct*xet 
           pBlk%xcz(i,j,k)=xet*yct-xct*yet
           
           pBlk%etx(i,j,k)=yct*zxc-yxc*zct
           pBlk%ety(i,j,k)=zct*xxc-zxc*xct 
           pBlk%etz(i,j,k)=xct*yxc-xxc*yct
           
           pBlk%ctx(i,j,k)=yxc*zet-yet*zxc
           pBlk%cty(i,j,k)=zxc*xet-zet*xxc 
           pBlk%ctz(i,j,k)=xxc*yet-xet*yxc
           pBlk%dj(i,j,k)=xxc*(yet*zct-yct*zet)-xet*(yxc*zct-yct*zxc)+xct*(yxc*zet-yet*zxc)
           if( pBlk%dj(i,j,k) .lt. 0) then 
             print*, "the djac less than zero,check the grid please!"
             stop
           elseif( pBlk%dj(i,j,k) .eq. 0) then
             print*," the djac equals zero, specify it 1.e-18!"
             pBlk%dj(i,j,k)=1.e-20
           endif          
       enddo
     enddo
  enddo  
  !jcmped buffer
  do i=pBlk%icmpst,pBlk%icmped
    do j=pBlk%jcmped+1,pBlk%jed
      do k=pBlk%kcmpst,pBlk%kcmped
  ! xc direction metrics    
         if(i.eq.pBlk%icmpst) then
           xxc=pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i,j,k)
           yxc=pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i,j,k)
           zxc=pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i,j,k) 
         elseif(i.eq.pBlk%icmped) then
           xxc=pBlk%xcoord(i,j,k)-pBlk%xcoord(i-1,j,k)
           yxc=pBlk%ycoord(i,j,k)-pBlk%ycoord(i-1,j,k)
           zxc=pBlk%zcoord(i,j,k)-pBlk%zcoord(i-1,j,k)
         elseif(i.eq.pBlk%icmpst+1 .or. i.eq.pBlk%icmped-1) then
           xxc=0.5*(pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i-1,j,k))
           yxc=0.5*(pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i-1,j,k))
           zxc=0.5*(pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i-1,j,k))
          else
           xxc=0.5*(pBlk%xcoord(i+1,j,k) - pBlk%xcoord(i-1,j,k))*4./3  - 0.5*(pBlk%xcoord(i+2,j,k) - pBlk%xcoord(i-2,j,k))/6.
           yxc=0.5*(pBlk%ycoord(i+1,j,k) - pBlk%ycoord(i-1,j,k))*4./3. - 0.5*(pBlk%ycoord(i+2,j,k) - pBlk%ycoord(i-2,j,k))/6.
           zxc=0.5*(pBlk%zcoord(i+1,j,k) - pBlk%zcoord(i-1,j,k))*4./3. - 0.5*(pBlk%zcoord(i+2,j,k) - pBlk%zcoord(i-2,j,k))/6.         
         endif   
   ! et direction metrics
          if(j.eq.pBlk%jst) then
           xet=pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j,k)
           yet=pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j,k)
           zet=pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j,k) 
         elseif(j.eq.pBlk%jed) then
           xet=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j-1,k)
           yet=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j-1,k)
           zet=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j-1,k)
         elseif(j.eq.pBlk%jst+1 .or. j.eq.pBlk%jed-1) then
           xet=0.5*(pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j-1,k))
           yet=0.5*(pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j-1,k))
           zet=0.5*(pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j-1,k))
          else
           xet=0.5*(pBlk%xcoord(i,j+1,k) - pBlk%xcoord(i,j-1,k))*4./3  - 0.5*(pBlk%xcoord(i,j+2,k) - pBlk%xcoord(i,j-2,k))/6.
           yet=0.5*(pBlk%ycoord(i,j+1,k) - pBlk%ycoord(i,j-1,k))*4./3. - 0.5*(pBlk%ycoord(i,j+2,k) - pBlk%ycoord(i,j-2,k))/6.
           zet=0.5*(pBlk%zcoord(i,j+1,k) - pBlk%zcoord(i,j-1,k))*4./3. - 0.5*(pBlk%zcoord(i,j+2,k) - pBlk%zcoord(i,j-2,k))/6.         
         endif 
    ! ct diretion
         if(k.eq.pBlk%kcmpst) then
           xct=pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k)
           yct=pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k)
           zct=pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k) 
         elseif(k.eq.pBlk%kcmped) then
           xct=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j,k-1)
           yct=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j,k-1)
           zct=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j,k-1)
         elseif(k.eq.pBlk%kcmpst+1 .or. k.eq.pBlk%kcmped-1) then
           xct=0.5*(pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k-1))
           yct=0.5*(pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k-1))
           zct=0.5*(pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k-1))
          else
           xct=0.5*(pBlk%xcoord(i,j,k+1) - pBlk%xcoord(i,j,k-1))*4./3  - 0.5*(pBlk%xcoord(i,j,k+2) - pBlk%xcoord(i,j,k-2))/6.
           yct=0.5*(pBlk%ycoord(i,j,k+1) - pBlk%ycoord(i,j,k-1))*4./3. - 0.5*(pBlk%ycoord(i,j,k+2) - pBlk%ycoord(i,j,k-2))/6.
           zct=0.5*(pBlk%zcoord(i,j,k+1) - pBlk%zcoord(i,j,k-1))*4./3. - 0.5*(pBlk%zcoord(i,j,k+2) - pBlk%zcoord(i,j,k-2))/6.         
         endif  
      ! compute the metrics
           pBlk%xcx(i,j,k)=yet*zct-yct*zet
           pBlk%xcy(i,j,k)=zet*xct-zct*xet 
           pBlk%xcz(i,j,k)=xet*yct-xct*yet
           
           pBlk%etx(i,j,k)=yct*zxc-yxc*zct
           pBlk%ety(i,j,k)=zct*xxc-zxc*xct 
           pBlk%etz(i,j,k)=xct*yxc-xxc*yct
           
           pBlk%ctx(i,j,k)=yxc*zet-yet*zxc
           pBlk%cty(i,j,k)=zxc*xet-zet*xxc 
           pBlk%ctz(i,j,k)=xxc*yet-xet*yxc
           pBlk%dj(i,j,k)=xxc*(yet*zct-yct*zet)-xet*(yxc*zct-yct*zxc)+xct*(yxc*zet-yet*zxc)
           if( pBlk%dj(i,j,k) .lt. 0) then 
             print*, "the djac less than zero,check the grid please!"
             stop
           elseif( pBlk%dj(i,j,k) .eq. 0) then
             print*," the djac equals zero, specify it 1.e-18!"
             pBlk%dj(i,j,k)=1.e-20
           endif          
       enddo
     enddo
  enddo
  !kcmpst buffer
    do i=pBlk%icmpst,pBlk%icmped
    do j=pBlk%jcmpst,pBlk%jcmped
      do k=pBlk%kst,pBlk%kcmpst-1
  ! xc direction metrics    
         if(i.eq.pBlk%icmpst) then
           xxc=pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i,j,k)
           yxc=pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i,j,k)
           zxc=pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i,j,k) 
         elseif(i.eq.pBlk%icmped) then
           xxc=pBlk%xcoord(i,j,k)-pBlk%xcoord(i-1,j,k)
           yxc=pBlk%ycoord(i,j,k)-pBlk%ycoord(i-1,j,k)
           zxc=pBlk%zcoord(i,j,k)-pBlk%zcoord(i-1,j,k)
         elseif(i.eq.pBlk%icmpst+1 .or. i.eq.pBlk%icmped-1) then
           xxc=0.5*(pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i-1,j,k))
           yxc=0.5*(pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i-1,j,k))
           zxc=0.5*(pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i-1,j,k))
          else
           xxc=0.5*(pBlk%xcoord(i+1,j,k) - pBlk%xcoord(i-1,j,k))*4./3  - 0.5*(pBlk%xcoord(i+2,j,k) - pBlk%xcoord(i-2,j,k))/6.
           yxc=0.5*(pBlk%ycoord(i+1,j,k) - pBlk%ycoord(i-1,j,k))*4./3. - 0.5*(pBlk%ycoord(i+2,j,k) - pBlk%ycoord(i-2,j,k))/6.
           zxc=0.5*(pBlk%zcoord(i+1,j,k) - pBlk%zcoord(i-1,j,k))*4./3. - 0.5*(pBlk%zcoord(i+2,j,k) - pBlk%zcoord(i-2,j,k))/6.         
         endif   
   ! et direction metrics
          if(j.eq.pBlk%jcmpst) then
           xet=pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j,k)
           yet=pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j,k)
           zet=pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j,k) 
         elseif(j.eq.pBlk%jcmped) then
           xet=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j-1,k)
           yet=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j-1,k)
           zet=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j-1,k)
         elseif(j.eq.pBlk%jcmpst+1 .or. j.eq.pBlk%jcmped-1) then
           xet=0.5*(pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j-1,k))
           yet=0.5*(pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j-1,k))
           zet=0.5*(pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j-1,k))
          else
           xet=0.5*(pBlk%xcoord(i,j+1,k) - pBlk%xcoord(i,j-1,k))*4./3  - 0.5*(pBlk%xcoord(i,j+2,k) - pBlk%xcoord(i,j-2,k))/6.
           yet=0.5*(pBlk%ycoord(i,j+1,k) - pBlk%ycoord(i,j-1,k))*4./3. - 0.5*(pBlk%ycoord(i,j+2,k) - pBlk%ycoord(i,j-2,k))/6.
           zet=0.5*(pBlk%zcoord(i,j+1,k) - pBlk%zcoord(i,j-1,k))*4./3. - 0.5*(pBlk%zcoord(i,j+2,k) - pBlk%zcoord(i,j-2,k))/6.         
         endif 
    ! ct diretion
         if(k.eq.pBlk%kst) then
           xct=pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k)
           yct=pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k)
           zct=pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k) 
         elseif(k.eq.pBlk%ked) then
           xct=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j,k-1)
           yct=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j,k-1)
           zct=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j,k-1)
         elseif(k.eq.pBlk%kst+1 .or. k.eq.pBlk%ked-1) then
           xct=0.5*(pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k-1))
           yct=0.5*(pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k-1))
           zct=0.5*(pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k-1))
          else
           xct=0.5*(pBlk%xcoord(i,j,k+1) - pBlk%xcoord(i,j,k-1))*4./3  - 0.5*(pBlk%xcoord(i,j,k+2) - pBlk%xcoord(i,j,k-2))/6.
           yct=0.5*(pBlk%ycoord(i,j,k+1) - pBlk%ycoord(i,j,k-1))*4./3. - 0.5*(pBlk%ycoord(i,j,k+2) - pBlk%ycoord(i,j,k-2))/6.
           zct=0.5*(pBlk%zcoord(i,j,k+1) - pBlk%zcoord(i,j,k-1))*4./3. - 0.5*(pBlk%zcoord(i,j,k+2) - pBlk%zcoord(i,j,k-2))/6.         
         endif  
      ! compute the metrics
           pBlk%xcx(i,j,k)=yet*zct-yct*zet
           pBlk%xcy(i,j,k)=zet*xct-zct*xet 
           pBlk%xcz(i,j,k)=xet*yct-xct*yet
           
           pBlk%etx(i,j,k)=yct*zxc-yxc*zct
           pBlk%ety(i,j,k)=zct*xxc-zxc*xct 
           pBlk%etz(i,j,k)=xct*yxc-xxc*yct
           
           pBlk%ctx(i,j,k)=yxc*zet-yet*zxc
           pBlk%cty(i,j,k)=zxc*xet-zet*xxc 
           pBlk%ctz(i,j,k)=xxc*yet-xet*yxc
           pBlk%dj(i,j,k)=xxc*(yet*zct-yct*zet)-xet*(yxc*zct-yct*zxc)+xct*(yxc*zet-yet*zxc)
           if( pBlk%dj(i,j,k) .lt. 0) then 
             print*, "the djac less than zero,check the grid please!"
             stop
           elseif( pBlk%dj(i,j,k) .eq. 0) then
             print*," the djac equals zero, specify it 1.e-18!"
             pBlk%dj(i,j,k)=1.e-20
           endif          
       enddo
     enddo
  enddo
!kcmped buffer
    do i=pBlk%icmpst,pBlk%icmped
    do j=pBlk%jcmpst,pBlk%jcmped
      do k=pBlk%kcmped+1,pBlk%ked
  ! xc direction metrics    
         if(i.eq.pBlk%icmpst) then
           xxc=pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i,j,k)
           yxc=pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i,j,k)
           zxc=pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i,j,k) 
         elseif(i.eq.pBlk%icmped) then
           xxc=pBlk%xcoord(i,j,k)-pBlk%xcoord(i-1,j,k)
           yxc=pBlk%ycoord(i,j,k)-pBlk%ycoord(i-1,j,k)
           zxc=pBlk%zcoord(i,j,k)-pBlk%zcoord(i-1,j,k)
         elseif(i.eq.pBlk%icmpst+1 .or. i.eq.pBlk%icmped-1) then
           xxc=0.5*(pBlk%xcoord(i+1,j,k)-pBlk%xcoord(i-1,j,k))
           yxc=0.5*(pBlk%ycoord(i+1,j,k)-pBlk%ycoord(i-1,j,k))
           zxc=0.5*(pBlk%zcoord(i+1,j,k)-pBlk%zcoord(i-1,j,k))
          else
           xxc=0.5*(pBlk%xcoord(i+1,j,k) - pBlk%xcoord(i-1,j,k))*4./3  - 0.5*(pBlk%xcoord(i+2,j,k) - pBlk%xcoord(i-2,j,k))/6.
           yxc=0.5*(pBlk%ycoord(i+1,j,k) - pBlk%ycoord(i-1,j,k))*4./3. - 0.5*(pBlk%ycoord(i+2,j,k) - pBlk%ycoord(i-2,j,k))/6.
           zxc=0.5*(pBlk%zcoord(i+1,j,k) - pBlk%zcoord(i-1,j,k))*4./3. - 0.5*(pBlk%zcoord(i+2,j,k) - pBlk%zcoord(i-2,j,k))/6.         
         endif   
   ! et direction metrics
          if(j.eq.pBlk%jcmpst) then
           xet=pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j,k)
           yet=pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j,k)
           zet=pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j,k) 
         elseif(j.eq.pBlk%jcmped) then
           xet=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j-1,k)
           yet=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j-1,k)
           zet=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j-1,k)
         elseif(j.eq.pBlk%jcmpst+1 .or. j.eq.pBlk%jcmped-1) then
           xet=0.5*(pBlk%xcoord(i,j+1,k)-pBlk%xcoord(i,j-1,k))
           yet=0.5*(pBlk%ycoord(i,j+1,k)-pBlk%ycoord(i,j-1,k))
           zet=0.5*(pBlk%zcoord(i,j+1,k)-pBlk%zcoord(i,j-1,k))
          else
           xet=0.5*(pBlk%xcoord(i,j+1,k) - pBlk%xcoord(i,j-1,k))*4./3  - 0.5*(pBlk%xcoord(i,j+2,k) - pBlk%xcoord(i,j-2,k))/6.
           yet=0.5*(pBlk%ycoord(i,j+1,k) - pBlk%ycoord(i,j-1,k))*4./3. - 0.5*(pBlk%ycoord(i,j+2,k) - pBlk%ycoord(i,j-2,k))/6.
           zet=0.5*(pBlk%zcoord(i,j+1,k) - pBlk%zcoord(i,j-1,k))*4./3. - 0.5*(pBlk%zcoord(i,j+2,k) - pBlk%zcoord(i,j-2,k))/6.         
         endif 
    ! ct diretion
         if(k.eq.pBlk%kst) then
           xct=pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k)
           yct=pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k)
           zct=pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k) 
         elseif(k.eq.pBlk%ked) then
           xct=pBlk%xcoord(i,j,k)-pBlk%xcoord(i,j,k-1)
           yct=pBlk%ycoord(i,j,k)-pBlk%ycoord(i,j,k-1)
           zct=pBlk%zcoord(i,j,k)-pBlk%zcoord(i,j,k-1)
         elseif(k.eq.pBlk%kst+1 .or. k.eq.pBlk%ked-1) then
           xct=0.5*(pBlk%xcoord(i,j,k+1)-pBlk%xcoord(i,j,k-1))
           yct=0.5*(pBlk%ycoord(i,j,k+1)-pBlk%ycoord(i,j,k-1))
           zct=0.5*(pBlk%zcoord(i,j,k+1)-pBlk%zcoord(i,j,k-1))
          else
           xct=0.5*(pBlk%xcoord(i,j,k+1) - pBlk%xcoord(i,j,k-1))*4./3  - 0.5*(pBlk%xcoord(i,j,k+2) - pBlk%xcoord(i,j,k-2))/6.
           yct=0.5*(pBlk%ycoord(i,j,k+1) - pBlk%ycoord(i,j,k-1))*4./3. - 0.5*(pBlk%ycoord(i,j,k+2) - pBlk%ycoord(i,j,k-2))/6.
           zct=0.5*(pBlk%zcoord(i,j,k+1) - pBlk%zcoord(i,j,k-1))*4./3. - 0.5*(pBlk%zcoord(i,j,k+2) - pBlk%zcoord(i,j,k-2))/6.         
         endif  
      ! compute the metrics
           pBlk%xcx(i,j,k)=yet*zct-yct*zet
           pBlk%xcy(i,j,k)=zet*xct-zct*xet 
           pBlk%xcz(i,j,k)=xet*yct-xct*yet
           
           pBlk%etx(i,j,k)=yct*zxc-yxc*zct
           pBlk%ety(i,j,k)=zct*xxc-zxc*xct 
           pBlk%etz(i,j,k)=xct*yxc-xxc*yct
           
           pBlk%ctx(i,j,k)=yxc*zet-yet*zxc
           pBlk%cty(i,j,k)=zxc*xet-zet*xxc 
           pBlk%ctz(i,j,k)=xxc*yet-xet*yxc
           pBlk%dj(i,j,k)=xxc*(yet*zct-yct*zet)-xet*(yxc*zct-yct*zxc)+xct*(yxc*zet-yet*zxc)
           if( pBlk%dj(i,j,k) .lt. 0) then 
             print*, "the djac less than zero,check the grid please!"
             stop
           elseif( pBlk%dj(i,j,k) .eq. 0) then
             print*," the djac equals zero, specify it 1.e-18!"
             pBlk%dj(i,j,k)=1.e-20
           endif          
       enddo
     enddo
  enddo
  
enddo  
return
end subroutine