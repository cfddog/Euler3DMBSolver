module md_filtering
integer ifilter
integer idir,jdir,kdir
integer :: ifinit=0
!//visbal-type filter right-coeffs
real :: af=0.49
real a20,a21,a40,a41,a42,a60,a61,a62,a63,a80,a81,a82,a84,&
     a100,a101,a102,a103,a104,a105
contains
!visbal-type filters 6th-order and 8th-order
subroutine setVisbalcoeffs()
    !//second order
    a20=1./2.+af
    a21=1./2.+af
    !fourth order
    a40=5./8.+6.*af/8.
    a41=1./2.+af
    a42=-1./8.+af/4.
    !sixth order
    a60=11./16.+5.*af/8.
    a61=(15+34*af)/32.
    a62=(-3.+6*af)/16
    a63=(1-2*af)/32.
    !eighth order
    a80=(93+70*af)/128.
    a81=(7+18*af)/16.
    a82=(-7.+14*af)/32.
    a83=(1-2*af)/16.
    a84=(-1+2*af)/128.
    !tenth order
    a100=(193+126*af)/256
    a101=(105+302*af)/256
    a102=(-15+30*af)/64
    a103=(45-90*af)/512
    a104=(-5+10*af)/256
    a105=(1-2*af)/512
    !//init ok
    ifinit=1
    return
end subroutine

subroutine FilterVisbalSixth(Nmax,varin)
    real varin(Nmax)
    real arr(Nmax,3),brr(Nmax),xrr(Nmax)
    if(ifinit .eq. 0) then
        call setVisbalcoeffs()
    endif
    !the matrix and right-terms
    do i=1,Nmax
        if(i .eq. 1 .or. i .eq. Nmax) then
            arr(i,1)=0.0;arr(i,2)=1.0;arr(i,3)=0.0
            brr(i)=varin(i)
        elseif(i .eq. 2 .or. i .eq. Nmax-1) then
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a20*varin(i)+a21*(varin(i-1)+varin(i+1))
        elseif(i .eq. 3 .or. i .eq. Nmax-2) then
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a40*varin(i)+a41*(varin(i-1)+varin(i+1))+a42*(varin(i-2)+varin(i+2))
        else
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a60*varin(i)+a61*(varin(i-1)+varin(i+1))+a62*(varin(i-2)+varin(i+2))+&
                   a63*(varin(i-3)+varin(i+3))
        endif          
    enddo 
    !//solve tri-band eqns
    call BandEqnSolve(Nmax,1,1,arr,brr,xrr)
    varin=xrr

    return
end subroutine

subroutine FilterVisbalEightth(Nmax,varin)
    real varin(Nmax)
    real arr(Nmax,3),brr(Nmax),xrr(Nmax)
    if(ifinit .eq. 0) then
        call setVisbalcoeffs()
    endif
    !the matrix and right-terms
    do i=1,Nmax
        if(i .eq. 1 .or. i .eq. Nmax) then
            arr(i,1)=0.0;arr(i,2)=1.0;arr(i,3)=0.0
            brr(i)=varin(i)
        elseif(i .eq. 2 .or. i .eq. Nmax-1) then
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a20*varin(i)+a21*(varin(i-1)+varin(i+1))
        elseif(i .eq. 3 .or. i .eq. Nmax-2) then
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a40*varin(i)+a41*(varin(i-1)+varin(i+1))+a42*(varin(i-2)+varin(i+2))
        elseif(i .eq. 4 .or. i .eq. Nmax-3) then
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a60*varin(i)+a61*(varin(i-1)+varin(i+1))+a62*(varin(i-2)+varin(i+2))+&
                                a63*(varin(i-3)+varin(i+3))
        else 
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a80*varin(i)+a81*(varin(i-1)+varin(i+1))+a82*(varin(i-2)+varin(i+2))+&
                                a83*(varin(i-3)+varin(i+3))+a84*(varin(i-4)+varin(i+4))                      
        endif          
    enddo 
    !//solve tri-band eqns
    call BandEqnSolve(Nmax,1,1,arr,brr,xrr)
    varin=xrr
    return
end subroutine

subroutine FilterVisbalTenth(Nmax,varin)
    real varin(Nmax)
    real arr(Nmax,3),brr(Nmax),xrr(Nmax)
    if(ifinit .eq. 0) then
        call setVisbalcoeffs()
    endif
    !the matrix and right-terms
    do i=1,Nmax
        if(i .eq. 1 .or. i .eq. Nmax) then
            arr(i,1)=0.0;arr(i,2)=1.0;arr(i,3)=0.0
            brr(i)=varin(i)
        elseif(i .eq. 2 .or. i .eq. Nmax-1) then
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a20*varin(i)+a21*(varin(i-1)+varin(i+1))
        elseif(i .eq. 3 .or. i .eq. Nmax-2) then
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a40*varin(i)+a41*(varin(i-1)+varin(i+1))+a42*(varin(i-2)+varin(i+2))
        elseif(i .eq. 4 .or. i .eq. Nmax-3) then
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a60*varin(i)+a61*(varin(i-1)+varin(i+1))+a62*(varin(i-2)+varin(i+2))+&
                                a63*(varin(i-3)+varin(i+3))
        elseif(i .eq. 5 .or. i .eq. Nmax-4) then 
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a80*varin(i)+a81*(varin(i-1)+varin(i+1))+a82*(varin(i-2)+varin(i+2))+&
                                a83*(varin(i-3)+varin(i+3))+a84*(varin(i-4)+varin(i+4)) 
        else 
            arr(i,1)=af;arr(i,2)=1.0;arr(i,3)=af 
            brr(i)=a100*varin(i)+a101*(varin(i-1)+varin(i+1))+a102*(varin(i-2)+varin(i+2))+&
                                 a103*(varin(i-3)+varin(i+3))+a104*(varin(i-4)+varin(i+4))+& 
                                 a105*(varin(i-5)+varin(i+5))                                
        endif          
    enddo 
    !//solve tri-band eqns
    call BandEqnSolve(Nmax,1,1,arr,brr,xrr)
    varin=xrr
    return
end subroutine

end module