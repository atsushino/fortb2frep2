    program po2d



    implicit none
    character(len=32)       Basename
    integer, allocatable :: psi(:,:), omg(:,:) ,tmp(:,:)
    integer(8)                 i, j, k, l, q, Nvert, FLAG
    integer(8)                 Ncell, Stop_itr
    real(8)                    Re, dt, eps, const
    real(8)                    h, rhsp, rhso, bb, err

    open(1, file='PARAM.dat')
    read(1,*)Basename, const, Re, dt, eps, Ncell, Stop_itr

    write(*,*)'============== Grid generation ==================='

    Nvert = Ncell + 1


    allocate(psi(Nvert,Nvert))
    allocate(omg(Nvert,Nvert))
    allocate(tmp(Nvert,Nvert))

    psi(:,:)=0.0
    omg(:,:)=0.0
    tmp(:,:)=0.0

    write(*,*)'Mesh size=',Ncell,'*',Ncell,'=',Ncell**2

    h = 1.0/float(Ncell)  !cellvertex

    write(*,*)'============== Grid generation: Done ============='




    !MAIN LOOP
    FLAG=0
    q=0
    l=0
    do while (FLAG == 0)


      !========boundary condition for cuvic cavity

      !left/right
      do j = 1, Nvert

        omg(1,j)      = -2.0*psi(2    ,j) /(h**2.0)
        omg(Nvert,j)  = -2.0*psi(Nvert,j) /(h**2.0)

      end do

      !bottom/top
      do i = 1, Nvert

        omg(i,1)      = -2.0*psi(i ,         2)/(h**2.0)
        omg(i,Nvert)  = -2.0*(psi(i,Nvert-1)+h)/(h**2.0)

      end do
      !========END boundary condition for cuvic cavity



      !calc psi

      do while (err.gt.eps)
        do j = 2, Nvert-1
          do i = 2, Nvert-1

            tmp(i,j) = psi(i,j)

            rhsp = (psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1))/4.0 &
                   +omg(i,j)*h**2.0/4.0-psi(i,j)

            psi(i,j) = psi(i,j)+const*rhsp
            !write(*,*)psi!!!!!!!!!!
          end do
        end do

        err=0.0
        do j = 2, Nvert-1
          do i = 2, Nvert-1

             bb=abs(psi(i,j)-tmp(i,j))

             if(bb.ge.err) then
              err = bb
            end if

          end do
        end do
      end do



      !calc omg
      do j = 2, Nvert-1
        do i = 2, Nvert-1

          tmp(i,j) = omg(i,j)

          rhso = ((psi(i+1,j)-psi(i-1,j))*(omg(i,j+1)-omg(i,j-1)) &
                 -(psi(i,j+1)-psi(i,j-1))*(omg(i+1,j)-omg(i-1,j)))/4.0 &
                 +(omg(i+1,j)+omg(i-1,j)+omg(i,j+1)+omg(i,j-1) &
                 -4.0*omg(i,j))/Re

          omg(i,j) = omg(i,j)+dt*rhso/(h**2)
          !write(*,*) omg!!!!!!!
        end do
      end do

      l=l+1
      write(*,*)'iteration=',l

      !close MAIN LOOP


      if(mod(l.Stop_itr) == Stop_itr) then
        write(*,*) "press 1=stop, else=continue"
        read(*,*) q
        if(q == 1) then
          FLAG=1
        end if
      end if



    end do

      write(*,*)psi
      write(*,*)
      write(*,*)omg
    stop
    end
