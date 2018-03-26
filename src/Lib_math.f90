function vminmod(a1,a2)
real vminmod
real a1,a2
 if(a1*a2.le.0) then
       vminmod=0.
else
      if(abs(a1) .le. abs(a2)) then
       vminmod=a1
      else
       vminmod=a2
      endif
 endif
 return
end function 
    
subroutine BandEqnSolve(Np,m1,m2,arr,brr,xrr)
real arr(Np,m1+m2+1),brr(Np),xrr(Np)
real al(Np,m1+m2+1)
integer indx(Np)

call bandec(arr,Np,m1,m2,Np,m1+m2+1,al,m1+m2+1,indx,d) !LU decompostion of arr
call banbks(arr,Np,m1,m2,Np,m1+m2+1,al,m1+m2+1,indx,brr) !band-eqn solution

xrr=brr

return
end 
      SUBROUTINE bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
	  !% a(N,m1+m2+1),m1 the elem num below diag, m2 the elem num above diag,np=n,
	  !% mpl=m1+m2+1,al(Np,mpl) store the lower part of a, indx(n),d=+1/-1
      INTEGER m1,m2,mp,mpl,n,np,indx(n)
      REAL d,a(np,mp),al(np,mpl),TINY
      PARAMETER (TINY=1.e-20)
      INTEGER i,j,k,l,mm
      REAL dum
      mm=m1+m2+1
      if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in bandec'
      l=m1
      do 13 i=1,m1
        do 11 j=m1+2-i,mm
          a(i,j-l)=a(i,j)
11      continue
        l=l-1
        do 12 j=mm-l,mm
          a(i,j)=0.
12      continue
13    continue
      d=1.
      l=m1
      do 18 k=1,n
        dum=a(k,1)
        i=k
        if(l.lt.n)l=l+1
        do 14 j=k+1,l
          if(abs(a(j,1)).gt.abs(dum))then
            dum=a(j,1)
            i=j
          endif
14      continue
        indx(k)=i
        if(dum.eq.0.) a(k,1)=TINY
        if(i.ne.k)then
          d=-d
          do 15 j=1,mm
            dum=a(k,j)
            a(k,j)=a(i,j)
            a(i,j)=dum
15        continue
        endif
        do 17 i=k+1,l
          dum=a(i,1)/a(k,1)
          al(k,i-k)=dum
          do 16 j=2,mm
            a(i,j-1)=a(i,j)-dum*a(k,j)
16        continue
          a(i,mm)=0.
17      continue
18    continue
      return
    END
	
    SUBROUTINE banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
      INTEGER m1,m2,mp,mpl,n,np,indx(n)
      REAL a(np,mp),al(np,mpl),b(n)
      INTEGER i,k,l,mm
      REAL dum
      mm=m1+m2+1
      if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in banbks'
      l=m1
      do 12 k=1,n
        i=indx(k)
        if(i.ne.k)then
          dum=b(k)
          b(k)=b(i)
          b(i)=dum
        endif
        if(l.lt.n)l=l+1
        do 11 i=k+1,l
          b(i)=b(i)-al(k,i-k)*b(k)
11      continue
12    continue
      l=1
      do 14 i=n,1,-1
        dum=b(i)
        do 13 k=2,l
          dum=dum-a(i,k)*b(k+i-1)
13      continue
        b(i)=dum/a(i,1)
        if(l.lt.mm) l=l+1
14    continue
      return
      END