    subroutine RHSperblk_Invis()
    use fieldpm
    use ctrlpm
    
    integer no_blk
    type(BLOCK_TYPE),pointer :: pBlk
    integer iblkst,iblked,jblkst,jblked,kblkst,kblked
	do iblk=1,numblk
       pBlk=>compblock(no_blk)
       iblkst=pBlk%icmpst;iblked=pBlk%icmped
       jblkst=pBlk%jcmpst;jblked=pBlk%jcmped
       kblkst=pBlk%kcmpst;kblked=pBlk%kcmped
       if(iblked-iblkst+1 .gt. 3) then
       iinc=1;jinc=0;kinc=0
       nst=pBlk%ist
       ned=pBlk%ied
	   !the invicid flux
       do j=jblkst,jblked
           do k=kblkst,kblked
               call SplittingFlux(nst,ned,nst,j,k,iinc,jinc,kinc,no_blk,iflagSplitSchm)
               call ReconstrctFlux(nst,ned,nst,j,k,iinc,jinc,kinc,no_blk,iConScheme)
           enddo
       enddo
       endif
       if(jblked-jblkst+1 .gt. 3) then
       iinc=0;jinc=1;kinc=0
       nst=pBlk%jst
       ned=pBlk%jed
       do i=iblkst,iblked
           do k=kblkst,kblked
               call SplittingFlux(nst,ned,i,nst,k,iinc,jinc,kinc,no_blk,iflagSplitSchm)
               call ReconstrctFlux(nst,ned,i,nst,k,iinc,jinc,kinc,no_blk,iConScheme)
           enddo
       enddo
       endif
       if(jblked-jblkst+1 .gt. 3) then
       iinc=0;jinc=0;kinc=1
       nst=pBlk%kst
       ned=pBlk%ked
       do i=iblkst,iblked
           do j=jblkst,jblked
               call SplittingFlux(nst,ned,i,j,nst,iinc,jinc,kinc,no_blk,iflagSplitSchm)
               call ReconstrctFlux(nst,ned,i,j,nst,iinc,jinc,kinc,no_blk,iConScheme)
           enddo
       enddo
       endif
    enddo
    
    return
    end subroutine
    
    subroutine SplittingFlux(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk,iflagSplitSchm) 
    use fieldpm
    use ctrlpm
        if(iflagSplitSchm .eq. 0) then
            !empty,dont split flux,usually compact or central scheme
			call  Splitting_none(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
        else    
            call  Splitting_SW(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
        endif
    return
    end subroutine
    
    subroutine ReconstrctFlux(iblkst,iblked,icur,jcur,kcur,iinc,jinc,kinc,no_blk,iConScheme)
    use fieldpm  
    use ctrlpm
	
	    call Scheme_OptC4(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
		
	return
    end subroutine