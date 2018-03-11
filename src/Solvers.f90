    subroutine Solvers()
      use Fieldpm
      use Freepm
      use Ctrlpm
      use constant
      ! select the Time-advanced methods 
      do iter=1,Nmaxstep
        call CPU_TIME(t1)
        if(itimeschm.eq.1) then
         call Runge_Kutta3(iter)   
        endif 
        call CPU_TIME(t2)
        print*,"Step",iter,"time comsuming:",t2-t1
      enddo
      return
    end subroutine
    
   subroutine Runge_Kutta3(iter)
       use fieldpm
       use ctrlpm
       type(BLOCK_TYPE),pointer :: pBlk
       real RK3(3,3)         !the coefficients of RK method
       integer iter
      !specify  the coefficients of rk method
      RK3(1,1)=1.
      RK3(1,2)=0.
      RK3(1,3)=1.
      RK3(2,1)=0.75
      RK3(2,2)=0.25
      RK3(2,3)=0.25
      RK3(3,1)=1./3.
      RK3(3,2)=2./3.
      RK3(3,3)=2./3.

!       call compute_dt()   !todo: cmp dt
!       cmptime=cmptime+dt  !todo: cmptime needed to be initialization
       ddrho=0.0
       do iRk=1,3 !todo: 3 steps needed
          call ComputeRHS    
          do iblk=1,numblk
            pBlk=>compblock(iblk)
            do i=pBlk%icmpst,pBlk%icmped
                do j=pBlk%jcmpst,pBlk%jcmped
                    do k=pBlk%kcmpst,pBlk%kcmped
                        do ivar=1,nvar
                            pBlk%Q(ivar,i,j,k)=Rk3(irk,1)*pBlk%Qn(ivar,i,j,k)+&
                                               Rk3(irk,2)*pBlk%Q(ivar,i,j,k)-&
                                               dt*Rk3(irk,3)*pBlk%RHS(ivar,i,j,k)/pBlk%dj(i,j,k) 
                               ddrho=ddrho+sqrt((pBlk%Q(ivar,i,j,k)-pBlk%Qn(ivar,i,j,k))**2)
                        enddo
                    enddo
                enddo
            enddo

          enddo
          !update the primitive variables
           call QvarToPvar() 
          !the boundaries including the block-block interface
           call EnforceBC() 
       enddo 
       ! check converged?
       ddrho=ddrho/(pblk%ied-pblk%ist-1)/(pblk%jed-pblk%jst-1)/(pblk%ked-pblk%kst-1)/3.
       if(ddrho .le. eps_solver) then
         print*, "the solution has been converged,and output results."
         call OutputResults()
         stop
       elseif(ddrho .gt. 1.e20) then
         print*, "the solution has been diverged!stop!"
         stop
       else
         print*, iter,dt,dt*iter,ddrho  !not steady situation, continue to compute
         if(mod(iter,iperiod) .eq. 0) then
         print*, "output results at iter=",iter
         call OutputResults()
         endif
!update the conservation law  
         do iblk=1,numblk
             pBlk=>compblock(iblk)
             pblk%Qn=pblk%Q
         enddo
       endif 
       
!!the plot3D flow field
! open(10,file="output.dat",form="formatted")
!    write(10,*) numblk
!    do iblk=1,numblk
!        pBlk=>compblock(iblk)
!        write(10,*) pblk%maxi,pblk%maxj,pblk%maxk
!    enddo
!    do iblk=1,numblk
!        pBlk=>compblock(iblk)
!        write(10,*) 0.0,0.0,0.0,0.0
!        write(10,*) (((pBlk%RHS(1,i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
!        (((pBlk%RHS(2,i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
!        (((pBlk%RHS(3,i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
!        (((pBlk%RHS(4,i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
!        (((pBlk%RHS(5,i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped)
!    enddo
!    close(10)
!    pause
       
       return
   end subroutine
   
   subroutine ComputeRHS()
   use fieldpm
   type(BLOCK_TYPE),pointer :: pBlk
   type(BC_TYPE),pointer :: pbound
   do iblk=1,numblk
      call RHSperblk(iblk)  
   enddo   !iblk
   !exchange the RHS at the buffer
    do iblk=1,numblk
        pblk=>compblock(iblk)
        do  ibc=1,pblk%num_BC
            pbound=>pblk%BC_MSG(ibc)
            if(pbound%ibctype .lt. 0) then
                call ExchangeBlkInfo_RHS(pbound%ibst,pbound%ibed,pbound%jbst,pbound%jbed,pbound%kbst,pbound%kbed,iblk,&
                                         pbound%ibst_targ,pbound%ibed_targ,pbound%jbst_targ,pbound%jbed_targ,pbound%kbst_targ,&
                                         pbound%kbed_targ,pbound%no_targblk)
            endif
        enddo
    enddo

    return 
    end subroutine 
!   
!    
    subroutine RHSPerBLK(no_blk)
    use fieldpm
    use ctrlpm
    
    integer no_blk
    type(BLOCK_TYPE),pointer :: pBlk
    integer iblkst,iblked,jblkst,jblked,kblkst,kblked
    pBlk=>compblock(no_blk)
    iblkst=pBlk%icmpst;iblked=pBlk%icmped
    jblkst=pBlk%jcmpst;jblked=pBlk%jcmped
    kblkst=pBlk%kcmpst;kblked=pBlk%kcmped
    
    pBlk%RHS=0.0
    if(iblked-iblkst+1 .gt. 3) then
    !i-direction
    iinc=1;jinc=0;kinc=0
    nst=pBlk%ist
    ned=pBlk%ied
    do j=jblkst,jblked
        do k=kblkst,kblked
            call SplittingFlux(nst,ned,nst,j,k,iinc,jinc,kinc,no_blk,iflagSplitSchm)
            call ReconstrctFlux(nst,ned,nst,j,k,iinc,jinc,kinc,no_blk,iConScheme)
        enddo
    enddo
    endif

    !j-direction
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

    !k-direction
    if(kblked-kblkst+1 .gt. 3) then
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

    return
    end subroutine
    
    subroutine SplittingFlux(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk,iflagSplitSchm) 
    use fieldpm
    use ctrlpm
    
    if(iSplitSchm .eq. 1) then
        call  Splitting_SW(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
    elseif(iSplitSchm .eq. 2) then
        call  Splitting_GblMaxEign(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
    elseif(iSplitSchm .eq. 3) then
        call  Splitting_LclMaxEign(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
    endif
    return
    end subroutine
    
    subroutine ReconstrctFlux(iblkst,iblked,icur,jcur,kcur,iinc,jinc,kinc,no_blk,iConScheme)
    use fieldpm  
    use ctrlpm
    
        if(iConSchm .eq. 1) then 
            call Scheme_NND2nd(iblkst,iblked,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
        endif
        
    return
    end subroutine