    
    subroutine Scheme_NND2nd(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
    use fieldpm
    type(BLOCK_TYPE),pointer :: pBlk
    real,pointer,dimension(:,:,:,:) :: pRHS !pointer to RHS
    integer no_blk
    pblk=>compblock(no_blk)
    i=icur;j=jcur;k=kcur
    idimst=iinc*pblk%ist+jinc*pblk%jst+kinc*pblk%kst
    idimed=iinc*pblk%ied+jinc*pblk%jed+kinc*pblk%ked
	
	if(iinc .eq. 1) then
	  pRHS=>pblk%RHSi
	elseif(jinc .eq. 1) then
      pRHS=>pBlk%RHSj
    else
      pRHS=>pBLK%RHSk
    endif
	
    do icnt=nst,ned-1
       do i1=1,nvar
           if(icnt .eq. idimst) then
               hflux(i1,icnt)=fp(i1,icnt  )
           else
               dFp =fp(i1,icnt+1) - fp(i1,icnt  )
               dFpn=fp(i1,icnt  ) - fp(i1,icnt-1)
               hflux(i1,icnt)=fp(i1,icnt  )+0.5*vminmod(dFp,dFpn)
           endif
       enddo
    enddo
    do icnt=nst,ned-1
       do i1=1,nvar
           if(icnt .eq. idimed-1) then
               hflux(i1,icnt)=hflux(i1,icnt)+fn(i1,icnt+1)
           else
               dFn =fn(i1,icnt+1) - fn(i1,icnt  )
               dFnp=fn(i1,icnt+2) - fn(i1,icnt+1)
               hflux(i1,icnt)=hflux(i1,icnt)+fn(i1,icnt+1)-0.5*vminmod(dFn,dFnp)
           endif
        enddo
    enddo
    
    !rhs
    do icnt=nst,ned
        do i1=1,nvar
        if(icnt .eq. idimst) then
            pRHS(i1,i,j,k)=pblk%RHS(i1,i,j,k)+fn(i1,icnt+1)-fn(i1,icnt)
        elseif(icnt .eq. idimed) then
            pRHS(i1,i,j,k)=pblk%RHS(i1,i,j,k)+fp(i1,icnt)-fp(i1,icnt-1)
        else
            pRHS(i1,i,j,k)= pblk%RHS(i1,i,j,k)+hflux(i1,icnt)-hflux(i1,icnt-1) 
        endif 
        enddo
        i=i+iinc
        j=j+jinc
        k=k+kinc
    enddo
    
    return
    end subroutine
	
	subroutine Scheme_OptC4(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
    use fieldpm
	use md_CompactSchm
    type(BLOCK_TYPE),pointer :: pBlk
    real,pointer,dimension(:,:,:,:) :: pRHS !pointer to RHS
	real varin(nst:ned),varout(nst:ned)
    integer no_blk
    pblk=>compblock(no_blk)
    i=icur;j=jcur;k=kcur
    idimst=iinc*pblk%ist+jinc*pblk%jst+kinc*pblk%kst
    idimed=iinc*pblk%ied+jinc*pblk%jed+kinc*pblk%ked
	
	if(iinc .eq. 1) then
	  pRHS=>pblk%RHSi
	elseif(jinc .eq. 1) then
      pRHS=>pBlk%RHSj
    else
      pRHS=>pBLK%RHSk
    endif
	do i1=1,5
       varin=flux(i1,nst:ned)
       call CentralCompact_Opt4(nst,ned,varin,varout)
	   i=icur;j=jcur;k=kcur
	   do icnt=nst,ned
	      pRHS(i1,i,j,k)=varout(icnt)
		  i=i+iinc
		  j=j+jinc
		  k=k+kinc
	   enddo
	enddo
	
    return
    end subroutine
    
