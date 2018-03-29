!compute the metrics for all zones
! at the first stage, just for one zone,use visbal's method or direct method 
subroutine ComputeMetrics()
use Ctrlpm
use fieldpm
use md_CompactSchm
type(BLOCK_TYPE),pointer :: pBlk
real,pointer,dimension(:,:,:) :: xxc,xet,xct,yxc,yet,yct,zxc,zet,zct
!now this is the non-SCMM method. this will be to do later. !todo: SCMM
do iblk=1,numblk
   pBlk=>compblock(iblk)
   allocate(xxc(pblk%icmpst:pblk%icmped,pblk%jcmpst:pblk%jcmped,pblk%kcmpst:pblk%kcmped))
   allocate(yxc(pblk%icmpst:pblk%icmped,pblk%jcmpst:pblk%jcmped,pblk%kcmpst:pblk%kcmped))
   allocate(zxc(pblk%icmpst:pblk%icmped,pblk%jcmpst:pblk%jcmped,pblk%kcmpst:pblk%kcmped))
   allocate(xet(pblk%icmpst:pblk%icmped,pblk%jcmpst:pblk%jcmped,pblk%kcmpst:pblk%kcmped))
   allocate(yet(pblk%icmpst:pblk%icmped,pblk%jcmpst:pblk%jcmped,pblk%kcmpst:pblk%kcmped))
   allocate(zet(pblk%icmpst:pblk%icmped,pblk%jcmpst:pblk%jcmped,pblk%kcmpst:pblk%kcmped))
   allocate(xct(pblk%icmpst:pblk%icmped,pblk%jcmpst:pblk%jcmped,pblk%kcmpst:pblk%kcmped))
   allocate(yct(pblk%icmpst:pblk%icmped,pblk%jcmpst:pblk%jcmped,pblk%kcmpst:pblk%kcmped))
   allocate(zct(pblk%icmpst:pblk%icmped,pblk%jcmpst:pblk%jcmped,pblk%kcmpst:pblk%kcmped))
   !x_xc,y_xc,z_xc
  do j=pBlk%jcmpst,pBlk%jcmped
    do k=pBlk%kcmpst,pBlk%kcmped
	     call Explicit_C4(pblk%icmpst,pblk%icmped,pBlk%xcoord(:,j,k),xxc(:,j,k))
		 call Explicit_C4(pblk%icmpst,pblk%icmped,pBlk%ycoord(:,j,k),yxc(:,j,k))
		 call Explicit_C4(pblk%icmpst,pblk%icmped,pBlk%zcoord(:,j,k),zxc(:,j,k))
	enddo
  enddo
 if(pBlk%kcmped-pBlk%kcmpst+1 .eq. 3) then
    do i=pBlk%icmpst,pBlk%icmped
    do k=pBlk%kcmpst,pBlk%kcmped
	   do j=pBlk%jcmpst,pBlk%jcmped
	     xet(i,j,k)=(pblk%xcoord(i,3,k)-pblk%xcoord(i,1,k))/2.
		 yet(i,j,k)=(pblk%ycoord(i,3,k)-pblk%ycoord(i,1,k))/2.
		 zet(i,j,k)=(pblk%zcoord(i,3,k)-pblk%zcoord(i,1,k))/2.
	   enddo
	enddo
    enddo
  else
  do i=pBlk%icmpst,pBlk%icmped
    do k=pBlk%kcmpst,pBlk%kcmped
	     call Explicit_C4(pblk%jcmpst,pblk%jcmped,pBlk%xcoord(i,:,k),xet(i,:,k))
		 call Explicit_C4(pblk%jcmpst,pblk%jcmped,pBlk%ycoord(i,:,k),yet(i,:,k))
		 call Explicit_C4(pblk%jcmpst,pblk%jcmped,pBlk%zcoord(i,:,k),zet(i,:,k))
	enddo
  enddo
  endif 
  if(pBlk%kcmped-pBlk%kcmpst+1 .eq. 3) then
    do i=pBlk%icmpst,pBlk%icmped
    do j=pBlk%jcmpst,pBlk%jcmped
	   do k=pBlk%kcmpst,pBlk%kcmped
	     xct(i,j,k)=(pblk%xcoord(i,j,3)-pblk%xcoord(i,j,1))/2.
		 yct(i,j,k)=(pblk%ycoord(i,j,3)-pblk%ycoord(i,j,1))/2.
		 zct(i,j,k)=(pblk%zcoord(i,j,3)-pblk%zcoord(i,j,1))/2.
	   enddo
	enddo
    enddo
  else
  do i=pBlk%icmpst,pBlk%icmped
    do j=pBlk%jcmpst,pBlk%jcmped
	     call Explicit_C4(pblk%jcmpst,pblk%jcmped,pBlk%xcoord(i,j,:),xct(i,j,:))
		 call Explicit_C4(pblk%jcmpst,pblk%jcmped,pBlk%ycoord(i,j,:),yct(i,j,:))
		 call Explicit_C4(pblk%jcmpst,pblk%jcmped,pBlk%zcoord(i,j,:),zct(i,j,:))
	enddo
   enddo
   endif   
    ! compute the metrics
    do i=pBlk%icmpst,pBlk%icmped
    do j=pBlk%jcmpst,pBlk%jcmped
    do k=pBlk%kcmpst,pBlk%kcmped	
       pBlk%xcx(i,j,k)=yet(i,j,k)*zct(i,j,k)-yct(i,j,k)*zet(i,j,k)
       pBlk%xcy(i,j,k)=zet(i,j,k)*xct(i,j,k)-zct(i,j,k)*xet(i,j,k) 
       pBlk%xcz(i,j,k)=xet(i,j,k)*yct(i,j,k)-xct(i,j,k)*yet(i,j,k)
                                                           
       pBlk%etx(i,j,k)=yct(i,j,k)*zxc(i,j,k)-yxc(i,j,k)*zct(i,j,k)
       pBlk%ety(i,j,k)=zct(i,j,k)*xxc(i,j,k)-zxc(i,j,k)*xct(i,j,k) 
       pBlk%etz(i,j,k)=xct(i,j,k)*yxc(i,j,k)-xxc(i,j,k)*yct(i,j,k)
                                                           
       pBlk%ctx(i,j,k)=yxc(i,j,k)*zet(i,j,k)-yet(i,j,k)*zxc(i,j,k)
       pBlk%cty(i,j,k)=zxc(i,j,k)*xet(i,j,k)-zet(i,j,k)*xxc(i,j,k) 
       pBlk%ctz(i,j,k)=xxc(i,j,k)*yet(i,j,k)-xet(i,j,k)*yxc(i,j,k)
       pBlk%dj(i,j,k)=xxc(i,j,k)*(yet(i,j,k)*zct(i,j,k)-yct(i,j,k)*zet(i,j,k))&
	                 -xet(i,j,k)*(yxc(i,j,k)*zct(i,j,k)-yct(i,j,k)*zxc(i,j,k))&
					 +xct(i,j,k)*(yxc(i,j,k)*zet(i,j,k)-yet(i,j,k)*zxc(i,j,k))
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
   !//deallocate the variables for the next block
   deallocate(xxc)
   deallocate(yxc)
   deallocate(zxc)
   deallocate(xet)
   deallocate(yet)
   deallocate(zet)
   deallocate(xct)
   deallocate(yct)
   deallocate(zct)
enddo
  
return
end subroutine