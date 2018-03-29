   subroutine gridBCIO
      use fieldpm
      type(BLOCK_TYPE),pointer :: pblk
      type(BC_TYPE),pointer :: pbound
!read the dimensions
      open(10,file="xyz.fmt",form="formatted")
      read(10,*) numblk
      allocate(compblock(numblk))
      do iblk=1,numblk
          
          pblk=>compblock(iblk)
          pblk%no_blk=iblk
          read(10,*) imax,jmax,kmax
          pblk%icmpst=1
          pblk%icmped=imax
          pblk%jcmpst=1
          pblk%jcmped=jmax
          pblk%kcmpst=1
          pblk%kcmped=kmax 
          pblk%maxi=imax
          pblk%maxj=jmax
          pblk%maxk=kmax
      enddo
      close(10)
      !....
!read the boundary conditions fillout the bclist
      open(10,file="generic.inp",form="formatted")
      read(10,*)        !blank,solver type
      read(10,*)        !blank,total blk number
      do iblk=1,numblk
          pblk=>compblock(iblk)
          read(10,*)        !blank,the dimensions for each blk
          read(10,*)        !blank,the blk name
          read(10,*) pblk%num_BC
          allocate(pblk%BC_MSG(pblk%num_BC))
          do ibc=1,pblk%num_BC
              pbound=>pblk%BC_MSG(ibc)
              read(10,*) pbound%ibst, pbound%ibed,pbound%jbst,pbound%jbed,pbound%kbst,pbound%kbed,pbound%ibctype
              if(pbound%ibctype .lt. 0) then
              read(10,*) pbound%ibst_targ, pbound%ibed_targ,pbound%jbst_targ,pbound%jbed_targ,pbound%kbst_targ,&
                         pbound%kbed_targ,pbound%no_targblk    
              elseif(pbound%ibctype .eq. 0) then
              print*,"WARNING: No Specified Boundary at (no_blk,ibc):",iblk,ibc
              endif
              if(pbound%ibst .eq. pbound%ibed) then  !i-const plane
                 if(abs(pbound%ibst) .eq. 1) then
                     if(pbound%ibctype .ge. 0) then
                         pblk%ist=1
                     else
                         pblk%ist=1-len_buf
                     endif
                 elseif(abs(pbound%ibst) .eq. pblk%maxi) then
                     if(pbound%ibctype .ge. 0) then
                         pblk%ied=pblk%maxi
                     else
                         pblk%ied=pblk%maxi+len_buf
                     endif
                 endif       
              elseif(pbound%jbst .eq. pbound%jbed) then !j-const plane
                 if(abs(pbound%jbst) .eq. 1) then
                     if(pbound%ibctype .ge. 0) then
                         pblk%jst=1
                     else
                         pblk%jst=1-len_buf
                     endif
                 elseif(abs(pbound%jbst) .eq. pblk%maxj) then
                     if(pbound%ibctype .ge. 0) then
                         pblk%jed=pblk%maxj
                     else
                         pblk%jed=pblk%maxj+len_buf
                     endif
                 endif                 
              elseif(pbound%kbst .eq. pbound%kbed) then !k-const plane
                  if(abs(pbound%kbst) .eq. 1) then
                     if(pbound%ibctype .ge. 0) then
                         pblk%kst=1
                     else
                         pblk%kst=1-len_buf
                     endif
                 elseif(abs(pbound%kbst) .eq. pblk%maxk) then
                     if(pbound%ibctype .ge. 0) then
                         pblk%ked=pblk%maxk
                     else
                         pblk%ked=pblk%maxk+len_buf
                     endif
                 endif                 
              endif       
          enddo
      enddo
    close(10)
!got the max dimension ed and min dimension st
   MINDIMST=100000
   MAXDIMED=0
    do iblk=1,numblk
        pblk=>compblock(iblk)
        MINDIMST=min(MINDIMST,min(pblk%ist,pblk%jst,pblk%kst))
        MAXDIMED=max(MAXDIMED,max(pblk%ied,pblk%jed,pblk%ked))
    enddo
!read the grid and allocate the memory
      open(10,file="xyz.fmt",form="formatted")
      read(10,*) numblk
      do iblk=1,numblk
          read(10,*) idum,jdum,kdum
          pblk=>compblock(iblk)
          allocate(pblk%xcoord(pblk%ist:pblk%ied,pblk%jst:pblk%jed,pblk%kst:pblk%ked))
          allocate(pblk%ycoord(pblk%ist:pblk%ied,pblk%jst:pblk%jed,pblk%kst:pblk%ked))
          allocate(pblk%zcoord(pblk%ist:pblk%ied,pblk%jst:pblk%jed,pblk%kst:pblk%ked))
      enddo
      do iblk=1,numblk
          pblk=>compblock(iblk)
          read(10,*) (((pblk%xcoord(i,j,k),i=pblk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
                   (((pblk%ycoord(i,j,k),i=pblk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
                   (((pblk%zcoord(i,j,k),i=pblk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped)
      enddo
      close(10)
      
 !extended grid at interface
    do iblk=1,numblk
        pblk=>compblock(iblk)
        do  ibc=1,pblk%num_BC
            pbound=>pblk%BC_MSG(ibc)
            if(pbound%ibctype .lt. 0) then
                call ExchangeBlkInfo(pbound%ibst,pbound%ibed,pbound%jbst,pbound%jbed,pbound%kbst,pbound%kbed,iblk,&
                                     pbound%ibst_targ,pbound%ibed_targ,pbound%jbst_targ,pbound%jbed_targ,pbound%kbst_targ,&
                                     pbound%kbed_targ,pbound%no_targblk)
            endif
        enddo
    enddo
    
      open(10,file="xyz-TEST.fmt",form="formatted")
      WRITE(10,*) numblk
      do iblk=1,numblk 
          pblk=>compblock(iblk)
          WRITE(10,*) pblk%ied-pblk%ist+1,pblk%jed-pblk%jst+1,pblk%ked-pblk%kst+1
      enddo
      do iblk=1,numblk
          pblk=>compblock(iblk)
          write(10,*) (((pblk%xcoord(i,j,k),i=pblk%ist,pblk%ied),j=pblk%jst,pblk%jed),k=pblk%kst,pblk%ked),&
                    (((pblk%ycoord(i,j,k),i=pblk%ist,pblk%ied),j=pblk%jst,pblk%jed),k=pblk%kst,pblk%ked),&
                    (((pblk%zcoord(i,j,k),i=pblk%ist,pblk%ied),j=pblk%jst,pblk%jed),k=pblk%kst,pblk%ked)
      enddo
      close(10)
      
    return
    end subroutine
    
    subroutine ExchangeBlkInfo(ist_src,ied_src,jst_src,jed_src,kst_src,ked_src,noblk_src,&
                               ist_tar,ied_tar,jst_tar,jed_tar,kst_tar,ked_tar,noblk_tar)
    use fieldpm
    type(BLOCK_TYPE),pointer :: pblk1,pblk2
    pblk1=>compblock(noblk_src)
    pblk2=>compblock(noblk_tar)
    !decide three kinds of plane: const,negative,positive
        if(abs(ist_src) .eq. abs(ied_src)) then
            nconst=abs(ist_src)
            
            if(jed_src .lt. 0) then
              nvarn=abs(jed_src-jst_src)+1     
              nvarp=abs(ked_src-kst_src)+1
            else
              nvarn=abs(ked_src-kst_src)+1     
              nvarp=abs(jed_src-jst_src)+1
            endif  
        elseif(abs(jst_src) .eq. abs(jed_src)) then
            nconst=abs(jst_src)
            if(ied_src .lt. 0) then
              nvarn=abs(ied_src-ist_src)+1     
              nvarp=abs(ked_src-kst_src)+1
            else
              nvarn=abs(ked_src-kst_src)+1     
              nvarp=abs(ied_src-ist_src)+1
            endif  
        elseif(abs(kst_src) .eq. abs(ked_src)) then
             nconst=abs(kst_src)
            if(ied_src .lt. 0) then
              nvarn=abs(ied_src-ist_src)+1     
              nvarp=abs(jed_src-jst_src)+1
            else
              nvarn=abs(jed_src-jst_src)+1     
              nvarp=abs(ied_src-ist_src)+1
            endif
        endif

        do indx_c=1,len_buf+1    !len_buf increments in const direction
            do indx_n=1,nvarn
                do indx_p=1,nvarp
                    
        !decide the increment of block 1 and block 2 interface plane
                   call decide_increment(i1,j1,k1,ist_src,ied_src,jst_src,jed_src,kst_src,ked_src,indx_c,indx_n,indx_p, 1) 
                   call decide_increment(i2,j2,k2,ist_tar,ied_tar,jst_tar,jed_tar,kst_tar,ked_tar,indx_c,indx_n,indx_p,-1)                                 

                    !specify the values
                    pblk1%xcoord(i1,j1,k1)=pblk2%xcoord(i2,j2,k2)
                    pblk1%ycoord(i1,j1,k1)=pblk2%ycoord(i2,j2,k2)
                    pblk1%zcoord(i1,j1,k1)=pblk2%zcoord(i2,j2,k2)
                enddo
            enddo
        enddo
        
    return
    end subroutine
    
    subroutine decide_increment(inc_i,inc_j,inc_k,ist,ied,jst,jed,kst,ked,indx_c,indx_n,indx_p,iflag)
    !iflag : make some difference between src blk(1) and tar blk(-1)
    if(ist .eq. ied) then
        !i-plane is const plane
        if((ied .gt. 0 .and. ied .eq. 1) .or. (ied .lt. 0) .and. (ied .ne. 1) ) then
           inc_i=ist+(-1)*iflag*(indx_c-1)
        else  
           inc_i=ist+(+1)*iflag*(indx_c-1) 
        endif
        !negative plane?
        if((jed .lt. 0)) then !j-plane is negative plane
           if(jed-jst .gt. 0) then
               inc_j=jst+(indx_n-1)
           else
               inc_j=jst-(indx_n-1)
           endif  
           if(ked-kst .gt. 0) then
               inc_k=kst+(indx_p-1)
           else
               inc_k=kst-(indx_p-1)
           endif 
        else                  !k-plane is nagative plane
           if(jed-jst .gt. 0) then
               inc_j=jst+(indx_p-1)
           else
               inc_j=jst-(indx_p-1)
           endif  
           if(ked-kst .gt. 0) then
               inc_k=kst+(indx_n-1)
           else
               inc_k=kst-(indx_n-1)
           endif 
        endif 
    inc_i=inc_i
    inc_j=abs(inc_j)
    inc_k=abs(inc_k)
    endif
    
    if(jst .eq. jed) then
        !j-plane is const plane
        if((jed .gt. 0 .and. jed .eq. 1) .or. (jed .lt. 0) .and. (jed .ne. 1) ) then
           inc_j=jst+(-1)*iflag*(indx_c-1)
        else  
           inc_j=jst+(+1)*iflag*(indx_c-1)
        endif
        !negative plane?
        if((ied .lt. 0)) then !i-plane is negative plane
           if(ied-ist .gt. 0) then
               inc_i=ist+(indx_n-1)
           else
               inc_i=ist-(indx_n-1)
           endif  
           if(ked-kst .gt. 0) then
               inc_k=kst+(indx_p-1)
           else
               inc_k=kst-(indx_p-1)
           endif 
        else                  !k-plane is nagative plane
           if(ied-ist .gt. 0) then
               inc_i=ist+(indx_p-1)
           else
               inc_i=ist-(indx_p-1)
           endif  
           if(ked-kst .gt. 0) then
               inc_k=kst+(indx_n-1)
           else
               inc_k=kst-(indx_n-1)
           endif 
        endif  
    inc_j=inc_j
    inc_i=abs(inc_i)
    inc_k=abs(inc_k)
    endif
    if(kst .eq. ked) then
        !k-plane is const plane
        if((ked .gt. 0 .and. ked .eq. 1) .or. (ked .lt. 0) .and. (ked .ne. 1) ) then
           inc_k=kst+(-1)*iflag*(indx_c-1)
        else  
           inc_k=kst+(+1)*iflag*(indx_c-1)
        endif
        !negative plane?
        if((jed .lt. 0)) then !j-plane is negative plane
           if(jed-jst .gt. 0) then
               inc_j=jst+(indx_n-1)
           else
               inc_j=jst-(indx_n-1)
           endif  
           if(ied-ist .gt. 0) then
               inc_i=ist+(indx_p-1)
           else
               inc_i=ist-(indx_p-1)
           endif 
        else                  !i-plane is nagative plane
           if(jed-jst .gt. 0) then
               inc_j=jst+(indx_p-1)
           else
               inc_j=jst-(indx_p-1)
           endif  
           if(ied-ist .gt. 0) then
               inc_i=ist+(indx_n-1)
           else
               inc_i=ist-(indx_n-1)
           endif 
        endif 
    inc_k=inc_k
    inc_j=abs(inc_j)
    inc_i=abs(inc_i)
    endif    
    

    return
    end subroutine