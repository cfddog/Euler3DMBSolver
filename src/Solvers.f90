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
	  call OutputResults()
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
           !call EnforceBC() 
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
        !update the conservation law，continue the inner iteration 
         do iblk=1,numblk
             pBlk=>compblock(iblk)
             pblk%Qn=pblk%Q
         enddo
        endif 
       
       return
   end subroutine
   
   subroutine ComputeRHS()
   use fieldpm
   type(BLOCK_TYPE),pointer :: pBlk
   type(BC_TYPE),pointer :: pbound
   !///
   call RHSperblk_Invis()  
   !ExchangeRHS()
   call AddToRHS()
   !call RHSperblk_Vis()
   !ExchangeRHS()
   !call AddToRHS()
   call RHSperblk_Forcing()
   !ExchangeRHS() 
   call AddToRHS()
   !
    return 
    end subroutine 
	
    subroutine AddToRHS()
	use fieldpm
	
	return
	end subroutine 
   
