subroutine RHSperblk_Forcing()
    use fieldpm
    use ctrlpm
    integer no_blk
    type(BLOCK_TYPE),pointer :: pBlk
    integer iblkst,iblked,jblkst,jblked,kblkst,kblked
	do iblk=1,numblk
       pBlk=>compblock(iblk)
       iblkst=pBlk%icmpst;iblked=pBlk%icmped
       jblkst=pBlk%jcmpst;jblked=pBlk%jcmped
       kblkst=pBlk%kcmpst;kblked=pBlk%kcmped
       if(iblked-iblkst+1 .gt. 3) then
       iinc=1;jinc=0;kinc=0
       nst=pBlk%ist
       ned=pBlk%ied
	   do j=jblkst,jblked
	   do k=kblkst,kblked
	      call NlinArtDissipTerms(nst,ned,nst,j,k,iinc,jinc,kinc,iblk)
	   enddo
       enddo
       endif
       if(jblked-jblkst+1 .gt. 3) then
       iinc=0;jinc=1;kinc=0
       nst=pBlk%jst
       ned=pBlk%jed
	   do i=iblkst,iblked
	   do k=kblkst,kblked
	      call NlinArtDissipTerms(nst,ned,i,nst,k,iinc,jinc,kinc,iblk)
	   enddo
       enddo
       endif
       if(jblked-jblkst+1 .gt. 3) then
       iinc=0;jinc=0;kinc=1
       nst=pBlk%kst
       ned=pBlk%ked
	   do j=jblkst,jblked
	   do i=iblkst,iblked
	      call NlinArtDissipTerms(nst,ned,i,j,kst,iinc,jinc,kinc,iblk)
	   enddo
       enddo
       endif
    enddo
return
end subroutine

subroutine NlinArtDissipTerms(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
    use fieldpm
	use constant
    type(BLOCK_TYPE),pointer :: pBlk
    real,pointer,dimension(:,:,:,:) :: pRHS !pointer to RHS
	real,pointer,dimension(:,:,:) :: varx,vary,varz
    real df(nst:ned,5)
    integer no_blk
    b1=-0.1624382574577463
    b2= 0.07309131357825455
    b3=-0.01447042896399915
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
	if(iinc .eq. 1) then    
     varx=>pblk%xcx; vary=>pblk%xcy; varz=>pblk%xcz 
    elseif(jinc .eq. 1) then
     varx=>pblk%etx; vary=>pblk%ety; varz=>pblk%etz 
    else
     varx=>pblk%ctx; vary=>pblk%cty; varz=>pblk%ctz 
    endif
	i=icur
	j=jcur
	k=kcur
	   pmin=pblk%p(i,j,k)
	   pmax=pblk%p(i,j,k)
	   Uc=varx(i,j,k)*pblk%u(i,j,k)+vary(i,j,k)*pblk%v(i,j,k)+varz(i,j,k)*pblk%w(i,j,k)
       cL=sqrt(gama*pblk%p(i,j,k)/pblk%rho(i,j,k))*sqrt(varx(i,j,k)**2+vary(i,j,k)**2+varz(i,j,k)**2)
       vmod=sqrt(varx(i,j,k)**2+vary(i,j,k)**2+varz(i,j,k)**2)
	   vlmin=Uc+cL
	   vlmax=Uc+cL
       vlmin1=(abs(Uc)+cL)/vmod
	   vlmax1=(abs(Uc)+cL)/vmod
	do icnt=nst,ned
	   pmin=min(pblk%p(i,j,k),pmin)
	   pmax=max(pblk%p(i,j,k),pmax)
	   Uc=varx(i,j,k)*pblk%u(i,j,k)+vary(i,j,k)*pblk%v(i,j,k)+varz(i,j,k)*pblk%w(i,j,k)
       cL=sqrt(gama*pblk%p(i,j,k)/pblk%rho(i,j,k))*sqrt(varx(i,j,k)**2+vary(i,j,k)**2+varz(i,j,k)**2)
	   vmod=sqrt(varx(i,j,k)**2+vary(i,j,k)**2+varz(i,j,k)**2)
	   vlmin=min(vlmin,abs(Uc)+cL)
	   vlmax=max(vlmax,abs(Uc)+cL)
	   vlmin1=min(vlmin1,(abs(Uc)+cL)/vmod)
	   vlmax1=max(vlmax1,(abs(Uc)+cL)/vmod)
	   i=i+iinc
       j=j+jinc
       k=k+kinc
	enddo
	alpham=vlmax/(vlmin+1.e-7)
	betam=vlmax1/(vlmin1+1.E-7)
	sigm=pmax/pmin
	alpham1=(alpham+1)/(alpham-1)*tanh(alpham-1)
	betam1=(betam+1)/(betam-1)*tanh(betam-1)
	Rj=(alpham+betam)/2./alpham/betam
	fkj=(1./sigm**Rj)*(1+(sigm-1)*tanh(alpham/betam-1))*sqrt(alpham1*betam1)**(1+tanh(sigm-1))
	!/////
	i=icur
	j=jcur
	k=kcur
	do icnt=nst,ned-1
	   if(icnt .le. nst+1 .or. icnt .ge. ned-2) then
          i=i+iinc
          j=j+jinc
          k=k+kinc		  
	   else
	      Uc=varx(i,j,k)*pblk%u(i,j,k)+vary(i,j,k)*pblk%v(i,j,k)+varz(i,j,k)*pblk%w(i,j,k)
          c=sqrt(gama*pblk%p(i,j,k)/pblk%rho(i,j,k))*sqrt(varx(i,j,k)**2+vary(i,j,k)**2+varz(i,j,k)**2)
          vlmin=abs(Uc)+c
          vlmax=abs(Uc)+c
          do ii=-2,3
		     iL=i+ii*iinc;jL=j+ii*jinc;kL=j+kinc*ii
		     Uc=varx(iL,jL,kL)*pblk%u(iL,jL,kL)+vary(iL,jL,kL)*pblk%v(iL,jL,kL)+&
			    varz(iL,jL,kL)*pblk%w(iL,jL,kL)
             cL=sqrt(gama*pblk%p(iL,jL,kL)/pblk%rho(iL,jL,kL))*&
			   sqrt(varx(iL,jL,kL)**2+vary(iL,jL,kL)**2+varz(iL,jL,kL)**2)
		     vlmin=min(vlmin,abs(Uc)+cL)
			 vlmax=max(vlmax,abs(Uc)+cL)
          enddo
          if(icnt .eq. nst+2 .or. icnt .eq. ned-3) then
		    eps2=0.0
            do ii=-1,2
		       iL=i+ii*iinc;jL=j+ii*jinc;kL=j+kinc*ii
			   rptmp=pblk%p(iL-iinc,jL-jinc,kL-kinc)-2.*pblk%p(iL,jL,kL)+pblk%p(iL+iinc,jL+jinc,kL+kinc)
			   eps2=max(eps2,abs(rptmp)/(rptmp+1.e-7))
            enddo			
          else
		    eps2=0.0
            do ii=-2,3
		       iL=i+ii*iinc;jL=j+ii*jinc;kL=j+kinc*ii
			   rptmp=pblk%p(iL-iinc,jL-jinc,kL-kinc)-2.*pblk%p(iL,jL,kL)+pblk%p(iL+iinc,jL+jinc,kL+kinc)
			   eps2=max(eps2,abs(rptmp)/(rptmp+1.e-7))
            enddo		  
          endif
          eps2=fkj*eps2	  
		  eps4=max(fkj-eps2,0.0)
		  vlstn=max(vlmax-vlmin,0.01*(vlmax+vlmin))
		  do i1=1,5
             df(icnt,i1)=vlstn*(eps2*(pblk%Q(i1,i+iinc,j+jinc,k+kinc)-pblk%Q(i1,i,j,k))-&
			                    eps4*(b1*(pblk%Q(i1,i+iinc,j+jinc,k+kinc)-pblk%Q(i1,i,j,k))+&
								      b2*(pblk%Q(i1,i+2*iinc,j+2*jinc,k+2*kinc)-pblk%Q(i1,i-iinc,j-jinc,k-kinc))+&
									  b3*(pblk%Q(i1,i+3*iinc,j+3*jinc,k+3*kinc)-pblk%Q(i1,i-2*iinc,j-2*jinc,k-2*kinc))) )
		  enddo	
          i=i+iinc
          j=j+jinc
          k=k+kinc		  
       endif	   
	enddo
!/////////////////	
	i=icur
	j=jcur
	k=kcur
	do icnt=nst,ned
	   if(icnt .le. nst+2 .or. icnt .ge. ned-2) then
          do i1=1,5
            pRHS(i1,i,j,k)=0.0
          enddo
	      i=i+iinc
          j=j+jinc
          k=k+kinc 
	   else
          do i1=1,5	   
	         pRHS(i1,i,j,k)=-(df(icnt,i1)-df(icnt-1,i1))  !move the forcing terms to left hand.so "minus"
		  enddo
		  i=i+iinc
          j=j+jinc
          k=k+kinc 
	   endif	  
	enddo
	

return
end subroutine