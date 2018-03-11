    
    subroutine Scheme_NND2nd(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
    use fieldpm
    type(BLOCK_TYPE),pointer :: pBlk
    integer no_blk
    pblk=>compblock(no_blk)
    i=icur;j=jcur;k=kcur
    idimst=iinc*pblk%ist+jinc*pblk%jst+kinc*pblk%kst
    idimed=iinc*pblk%ied+jinc*pblk%jed+kinc*pblk%ked
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
            pblk%RHS(i1,i,j,k)=pblk%RHS(i1,i,j,k)+fn(i1,icnt+1)-fn(i1,icnt)
        elseif(icnt .eq. idimed) then
            pblk%RHS(i1,i,j,k)=pblk%RHS(i1,i,j,k)+fp(i1,icnt)-fp(i1,icnt-1)
        else
            pblk%RHS(i1,i,j,k)= pblk%RHS(i1,i,j,k)+hflux(i1,icnt)-hflux(i1,icnt-1) 
        endif 
        enddo
        i=i+iinc
        j=j+jinc
        k=k+kinc
    enddo
    
    return
    end subroutine
    
    subroutine Scheme_stdWENO5(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
    use fieldpm
    use ctrlpm
    type(BLOCK_TYPE),pointer :: pBlk
    integer no_blk
    real    :: WENO_eps=1.e-6
    real    :: LM(5,5),RM(5,5),fvar(5,-3:3),hh(5)
    
    pblk=>compblock(no_blk)
    i=icur;j=jcur;k=kcur
    idimst=iinc*pblk%ist+jinc*pblk%jst+kinc*pblk%kst
    idimed=iinc*pblk%ied+jinc*pblk%jed+kinc*pblk%ked
    !weno linear weight
     CL0=0.1;CL1=0.6;CL2=0.3
     c01=1./3.;c02=-7./6.;c03=11./6.
     c11=-1./6.;c12=5./6.;c13=1./3.
     c21=1./3.;c22=5./6.;c23=-1./6.
     
    do icnt=nst,ned-1      
           if(icnt .eq. idimst) then
               do i1=1,nvar
                   hp(i1,icnt)=fp(i1,icnt)
               enddo
           elseif(icnt .eq.idimst+1) then
               do i1=1,nvar
                   dFp =fp(i1,icnt+1)-fp(i1,icnt)
                   dFpn=fp(i1,icnt)-fp(i1,icnt-1)
                   hp(i1,icnt)=fp(i1,icnt)+0.5*vminmod(dFp,dFpn)
               enddo
           elseif(icnt .eq. idimed-1) then
               do i1=1,nvar
                  dFp =fp(i1,icnt+1)-fp(i1,icnt  )
                  dFpn=fp(i1,icnt  )-fp(i1,icnt-1)
                  hp(i1,icnt)=fp(i1,icnt)+0.5*vminmod(dFp,dFpn)
               enddo
           else
              do i1=1,nvar
              if(iproj .eq. 1) then
                call  RLmaxtrix(no_blk,i,j,k,iinc,jinc,kinc,RM,LM)
                !projection
                 do ii=-2,2              
                     fvar(i1,ii)=0.0
                     do kk=1,5
                         fvar(i1,ii)=fvar(i1,ii)+LM(i1,kk)*fp(kk,icnt+ii)
                     enddo
                 enddo
              else
                  do ii=-2,2
                      fvar(i1,ii)=fp(i1,icnt+ii)
                  enddo
               endif
               h0=c01*fvar(i1,-2)+c02*fvar(i1,-1)+c03*fvar(i1,0  )
               h1=c11*fvar(i1,-1)+c12*fvar(i1, 0)+c13*fvar(i1,  1)
               h2=c21*fvar(i1,0 )+c22*fvar(i1, 1)+c23*fvar(i1,  2)
               vIS0=13./12.*(fvar(i1,-2)-2.*fvar(i1, -1)+   fvar(i1,0))**2+&
                    1./4.  *(fvar(i1,-2)-4.*fvar(i1, -1)+3.*fvar(i1,0))**2
               vIS1=13./12.*(fvar(i1,-1)-2.*fvar(i1,0  )+   fvar(i1,1))**2+&
                    1./4.  *(fvar(i1,-1)-fvar(i1,1))**2
               vIS2=13./12.*(fvar(i1,0  )-2.*fvar(i1, 1)+fvar(i1,  2))**2+&
                    1./4.  *(3.*fvar(i1,0)-4.*fvar(i1,1)+fvar(i1,  2))**2
               beta0=CL0/(vIS0+WENO_eps)**2
               beta1=CL1/(vIS1+WENO_eps)**2
               beta2=CL2/(vIS2+WENO_eps)**2
               w0=beta0/(beta0+beta1+beta2)
               w1=beta1/(beta0+beta1+beta2)
               w2=beta2/(beta0+beta1+beta2)
               hh(i1)=w0*h0+w1*h1+w2*h2
               enddo
               if(iproj .eq. 1) then
               !anti-projection
               do i1=1,nvar
                   hp(i1,icnt)=0.0
                 do kk=1,5
                     hp(i1,icnt)=hp(i1,icnt)+RM(i1,kk)*hh(kk)
                 enddo
               enddo
               else
                   do i1=1,nvar
                       hp(i1,icnt)=hh(i1)
                   enddo
               endif
           endif
       i=i+iinc
       j=j+jinc
       k=k+kinc
    enddo
    !//
    i=icur;j=jcur;k=kcur
    hflux(:,:)=0.0
    do icnt=nst,ned-1
           if(icnt .eq. idimst) then
               do i1=1,nvar
               dFn =fn(i1,icnt+1) - fn(i1,icnt)
               dFnp=fn(i1,icnt+2) - fn(i1,icnt+1)
               hn(i1,icnt)=fn(i1,icnt+1)-0.5*vminmod(dFn,dFnp)
               enddo
           elseif(icnt .eq. idimed-2) then
               do i1=1,nvar
               dFn =fn(i1,icnt+1) - fn(i1,icnt  )
               dFnp=fn(i1,icnt+2) - fn(i1,icnt+1)
               hn(i1,icnt)=fn(i1,icnt+1)-0.5*vminmod(dFn,dFnp)
               enddo
           elseif(icnt .eq. idimed-1) then
               do i1=1,nvar
               hn(i1,icnt)=fn(i1,icnt+1)
               enddo
           else
               do i1=1,nvar
               if(iproj .eq. 1) then
               call  RLmaxtrix(no_blk,i,j,k,iinc,jinc,kinc,RM,LM)
               !projection
                do ii=-1,3               
                    fvar(i1,ii)=0.0
                    do kk=1,5
                        fvar(i1,ii)=fvar(i1,ii)+LM(i1,kk)*fn(kk,icnt+ii)
                    enddo
                enddo
               else
                   do ii=-1,3
                       fvar(i1,ii)=fn(i1,icnt+ii)
                   enddo
               endif
               h0=c01*fvar(i1,3)+c02*fvar(i1,2)+c03*fvar(i1,+1)
               h1=c11*fvar(i1,2)+c12*fvar(i1,1)+c13*fvar(i1, 0)
               h2=c21*fvar(i1,1)+c22*fvar(i1,0) +c23*fvar(i1,-1)
               vIS0=13./12.*(fvar(i1,3)-2.*fvar(i1,2)+   fvar(i1,+1))**2+&
                    1./4.  *(fvar(i1,3)-4.*fvar(i1,2)+3.*fvar(i1,+1))**2
               vIS1=13./12.*(fvar(i1,2)-2.*fvar(i1,1)   +fvar(i1, 0))**2+&
                    1./4.  *(fvar(i1,2)-fvar(i1,0))**2
               vIS2=13./12.*(fvar(i1,1)-2.*fvar(i1,0)+fvar(i1,-1))**2+&
                    1./4.  *(3.*fvar(i1,1)-4.*fvar(i1,0)+fvar(i1,-1))**2
               beta0=CL0/(vIS0+WENO_eps)**2
               beta1=CL1/(vIS1+WENO_eps)**2
               beta2=CL2/(vIS2+WENO_eps)**2
               w0=beta0/(beta0+beta1+beta2)
               w1=beta1/(beta0+beta1+beta2)
               w2=beta2/(beta0+beta1+beta2)
               hh(i1)=w0*h0+w1*h1+w2*h2
               enddo
               if(iproj .eq. 1) then
               !anti-projection
               do i1=1,nvar
                   hn(i1,icnt)=0.0
                 do kk=1,5
                     hn(i1,icnt)=hn(i1,icnt)+RM(i1,kk)*hh(kk)
                 enddo
               enddo
               else
                   do i1=1,nvar
                       hn(i1,icnt)=hh(i1)
                   enddo
               endif
           endif
       i=i+iinc
       j=j+jinc
       k=k+kinc
    enddo               
    !rhs
    i=icur;j=jcur;k=kcur
    do icnt=nst,ned
        do i1=1,nvar
        if(icnt .eq. idimst) then
            pblk%RHS(i1,i,j,k)=pblk%RHS(i1,i,j,k)+fn(i1,icnt+1)-fn(i1,icnt)
        elseif(icnt .eq. idimed) then
            pblk%RHS(i1,i,j,k)=pblk%RHS(i1,i,j,k)+fp(i1,icnt)-fp(i1,icnt-1)
        else
            pblk%RHS(i1,i,j,k)=pblk%RHS(i1,i,j,k)+hp(i1,icnt)-hp(i1,icnt-1)+hn(i1,icnt)-hn(i1,icnt-1)
        endif  
        enddo
        i=i+iinc
        j=j+jinc
        k=k+kinc
    enddo
    !.........nothing
    return
    end subroutine
    !//Pirozzoli hybrid-compact-weno scheme
    subroutine HybridCompactWeno(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
    use fieldpm
    use ctrlpm
    type(BLOCK_TYPE),pointer :: pBlk
    integer no_blk
    real    :: rho_eps=0.05
    real    :: LM(5,5),RM(5,5),fvar(5,-3:3),hh(5)
    real    :: arr(nst-1:ned,3),brr(nst-1:ned),xrr(nst-1:ned) !solve a*x=b
    real    :: h
    i=icur
    j=jcur
    k=kcur
    do i1=1,5
       
    do icnt=nst-1,ned
       if(icnt .eq. nst-1) then
           arr(icnt,1)=0.;arr(icnt,2)=1.;arr(icnt,3)=0.
           brr(icnt)=25./12.*fp(i1,icnt+1)-23./12.*fp(i1,icnt+2)-13./12.*fp(i1,icnt+3)-1./4.*fp(i1,icnt+4)
       elseif(icnt .eq. nst) then
           arr(icnt,1)=0.;arr(icnt,2)=1.;arr(icnt,3)=0.
           brr(icnt)=1./4.*fp(i1,icnt)+13./12.*fp(i1,icnt+1)-5./12.*fp(i1,icnt+2)+1./12.*fp(i1,icnt+3)
       elseif(icnt .eq. ned) then
           arr(icnt,1)=0.;arr(icnt,2)=1.;arr(icnt,3)=0.
           brr(icnt)=25./12.*fp(i1,icnt)-23./12.*fp(i1,icnt-1)-13./12.*fp(i1,icnt-2)-1./4.*fp(i1,icnt-3)
       else
           arr(icnt,1)=3.;arr(icnt,2)=6.;arr(icnt,3)=1.
           brr(icnt)=1./3.*fp(i1,icnt-1)+19./3.*fp(i1,icnt)+10./3.*fp(i1,icnt+1)
       endif
       !NND-2 subroutine if oscillation ocurred
       beta_r=abs(rho(i+iinc,j+jinc,k+kinc)-rho(i,j,k))
       if(beta_r .le. rho_eps) then
          arr(icnt,1)=0.;arr(icnt,2)=1.;arr(icnt,3)=0.
          dFp =fp(i1,icnt+1) - fp(i1,icnt  )
          dFpn=fp(i1,icnt  ) - fp(i1,icnt-1)
          brr(icnt)=fp(i1,icnt  )+0.5*vminmod(dFp,dFpn) 
       endif
    enddo
    
    return
    end subroutine
    
