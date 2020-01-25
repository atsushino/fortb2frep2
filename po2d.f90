    program po2d



    implicit none
    integer, allocatable :: psi(:,:), omg(:,:) ,tmp(:,:)
    integer                 i, j, k, l
    integer                 Nmesh, Stop_itr, GZ_itr
    real                    Re, dt, eps, const
    real                    h, rhsp, rhso, bb, err


    write(*,*)'============== Grid generation ==================='
    write(*,*) 'input meshsize'
    read(*,*) Nmesh

    allocate(psi(Nmesh,Nmesh))
    allocate(omg(Nmesh,Nmesh))
    allocate(tmp(Nmesh,Nmesh))

    psi(:,:)=0.0
    omg(:,:)=0.0
    tmp(:,:)=0.0

    write(*,*)'Mesh size=',Nmesh,'*',Nmesh,'=',Nmesh**2

    h = 1.0/float(Nmesh-1)  !cellvertex

    write(*,*)'============== Grid generation: Done ============='


    Stop_itr=100

    !MAIN LOOP
    do l = 1, Stop_itr

      !========boundary condition for cuvic cavity

      !left/right
      do j = 1, Nmesh

        omg(1,j)      = -2.0*psi(2      ,j) /h**2.0
        omg(Nmesh,j)  = -2.0*psi(Nmesh-1,j) /h**2.0

      end do

      !bottom/top
      do i = 1, Nmesh

        omg(i,1)      = -2.0*psi(i,2)      /h**2.0
        omg(i,Nmesh)  = -2.0*(psi(i,Nmesh-1)+h)/h**2.0

      end do

      !========END boundary condition for cuvic cavity



      !calc psi
      GZ_itr = 100
      const = 1.0
      Re=0.20
      dt= 0.0025
      eps=0.00001
      do while (err.gt.eps)
        do j = 2, Nmesh-1
          do i = 2, Nmesh-1

            tmp(i,j) = psi(i,j)

            rhsp = (psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1))/4.0 &
                  +omg(i,j)*h**2.0/4.0-psi(i,j)

            psi(i,j) = psi(i,j)+const*rhs

          end do
        end do

        err=0.0
        do j = 2, Nmesh-1
          do i = 2, Nmesh-1
            bb=abs(psi(i,j)-tmp(i,j))
            if(bb.ge.err) err = bb
          end do
        end do

      end do



      !calc omg
      do k = 1, GZ_itr
        do j = 2, Nmesh-1
          do i = 2, Nmesh-1

            tmp(i,j) = omg(i,j)

            rhso = ((psi(i+1,j)-psi(i-1,j))*(omg(i,j+1)-omg(i,j-1)) &
                  -(psi(i,j+1)-psi(i,j-1))*(omg(i+1,j)-omg(i-1,j)))/4.0 &
                  +(omg(i+1,j)-omg(i-1,j)+omg(i,j+1)-omg(i,j-1) &
                  -4.0*omg(i,j))/Re

            omg(i,j) = omg(i,j)+dt*rhs/h**2

          end do
        end do
      end do



      write(*,*)psi,omg
      write(*,*)'==========================',l

      !close MAIN LOOP
    end do

    stop
    end
