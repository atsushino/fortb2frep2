    program po2d



    implicit none
    character(len=32)          Basename
    real(8), allocatable ::    psi(:,:), omg(:,:) ,tmp(:,:)
    integer(8)                 i, j, k, l, m, q, Nvert, FLAG
    integer(8)                 Ncell, Stop_itr
    real(8)                    urel, Re, dt, eps, const, Calc_max
    real(8)                    h, rhsp, rhso, err, Err_max

    open(1, file='PARAM.dat',                     action="read")
    open(2, file='psi.txt',     status="replace", action="write")
    open(3, file='omg.txt',     status="replace", action="write")
    open(4, file='po2dlog.txt', status="replace", action="write")



    read(1,*) Basename, Ncell, urel, Re, const, dt, eps ,Stop_itr, Calc_max

    write(*,*)'============================================'
    write(*,*)' 2D PSI-OMEGA METHOD CFD ANALYSIS'
    write(*,*)' ver 1.1.0'
    write(*,*)' by ATSUSHI SHINOZUKA'
    write(*,*)' CASE: UNSTEDY FLOW IN CUBIC CAVITY'
    write(*,*)' METHOD: CENTRAL DIFFERENCE'
    write(*,*)'         GAUSS=SEIDEL METHOD'
    write(*,*)' Set parameters in PARAM.dat '
    write(*,*)' Check for update on GitHub:'
    write(*,*)'https://github.com/atsushino/fortb2frep2'
    write(*,*)'============================================'
    write(*,*)

    write(4,*)'============================================'
    write(4,*)' 2D PSI-OMEGA METHOD CFD ANALYSIS'
    write(4,*)' ver 1.1.0'
    write(4,*)' by ATSUSHI SHINOZUKA'
    write(4,*)' CASE: UNSTEDY FLOW IN CUBIC CAVITY'
    write(4,*)' METHOD: CENTRAL DIFFERENCE'
    write(4,*)'         GAUSS=SEIDEL METHOD'
    write(4,*)' Check for update on GitHub:'
    write(4,*)'https://github.com/atsushino/fortb2frep2'
    write(4,*)'============================================'
    write(4,*)



    write(4,*)'======= Grid generation ============='

    write(*,*)'======= Grid generation ============='
    write(4,*)'======= Grid generation ============='

    Nvert = Ncell + 1           !cellvertex

    allocate(psi(Nvert,Nvert))  !cellvertex
    allocate(omg(Nvert,Nvert))
    allocate(tmp(Nvert,Nvert))

    psi(:,:)=0.0
    omg(:,:)=0.0
    tmp(:,:)=0.0

    write(*,*)'Grid size=',Ncell,'*',Ncell,'=',Ncell**2
    write(4,*)'Grid size=',Ncell,'*',Ncell,'=',Ncell**2


    h = 1.0/float(Ncell)

    write(*,*)'======= Grid generation: Done ======='
    write(*,*)
    write(4,*)'======= Grid generation: Done ======='
    write(4,*)



    !MAIN LOOP
    write(*,*)'======= Psi-Omega Method ============'
    write(4,*)'======= Psi-Omega Method ============'

    FLAG=0
    q=0
    l=0
    m=0

    do while (FLAG == 0)
      l=l+1

      i=0.0
      j=0.0

      !======= boundary condition for cuvic cavity

      !left/right
      do j = 1, Nvert

        omg(1,j)      = -2.0*psi(2      ,j) /(h**2.0)
        omg(Nvert,j)  = -2.0*psi(Nvert-1,j) /(h**2.0)

      end do

      !bottom/top
      do i = 1, Nvert

        omg(i,1)      = -2.0*psi(i ,              2)/(h**2.0)
        omg(i,Nvert)  = -2.0*(psi(i,Nvert-1)+urel*h)/(h**2.0)

      end do
      !======= boundary condition for cuvic cavity: END


      !======= Calc psi and omega : Gauss-Seidel Method

      !calc psi

      m = 0             !reset counter GS_itr
      Err_max=eps*1.1   !enter DO loop below

      do while (FLAG==0 .and. Err_max.gt.eps)

        m = m + 1
        Err_max=0.0



        i=0.0
        j=0.0
        do j = 2, Nvert-1
          do i = 2, Nvert-1

            tmp(i,j) = psi(i,j)   !store current psi

            rhsp = (psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1))/4.0 &
                   +omg(i,j)*(h**2.0)/4.0-psi(i,j)

            psi(i,j) = psi(i,j)+const*rhsp  !calc new psi

          end do
        end do



        i=0.0
        j=0.0

        do j = 1, Nvert
          do i = 1, Nvert

            if(abs(psi(i,j)).ge.Calc_max .or. abs(omg(i,j)).ge.Calc_max) then
              FLAG=1   !ABRT calclatio to aboid overflow
              write(*,*)'Error: Check PARAM.dat'
              write(4,*)'Error: Check PARAM.dat'

              stop
            end if

          end do
        end do



        i=0.0
        j=0.0

        do j = 2, Nvert-1
          do i = 2, Nvert-1

            err=abs(psi(i,j)-tmp(i,j))

            if(err.ge.Err_max) then
              Err_max = err         !calc Err for GSitr
            end if

          end do
        end do



        write(*,*)'Itr',l,'GSitr',m,'Err_max',Err_max
        write(4,*)'Itr',l,'GSitr',m,'Err_max',Err_max

      end do



      !calc omg
      i=0.0
      j=0.0
      do j = 2, Nvert-1
        do i = 2, Nvert-1

          tmp(i,j) = omg(i,j)

          rhso = ((psi(i+1,j)-psi(i-1,j))*(omg(i,j+1)-omg(i,j-1)) &
                -(psi(i,j+1)-psi(i,j-1))*(omg(i+1,j)-omg(i-1,j)))/4.0 &
                +(omg(i+1,j)+omg(i-1,j)+omg(i,j+1)+omg(i,j-1) &
                -4.0*omg(i,j))/Re

          omg(i,j) = omg(i,j)+dt*rhso/(h**2)

        end do
      end do



    !close MAIN LOOP

    !  if(FLAG == 1) then

    !    write(*,*) 'psi'
    !    do j=1, Nvert
    !    	write(*,*) (psi(i,j),i=1, Nvert)
    !    end do

    !    write(*,*) 'omg'
    !    do j=1, Nvert
    !    	write(*,*) (omg(i,j),i=1, Nvert)
    !    end do



      if (mod(l,Stop_itr) == 0) then

        write(*,*) 'psi'
        do j=1, Nvert
        	write(*,*) (psi(i,j),i=1, Nvert)
        end do

        write(*,*) 'omg'
        do j=1, Nvert
        	write(*,*) (omg(i,j),i=1, Nvert)
        end do

        write(*,*) "press 1=stop, 2 =continue"
        read(*,*) q
        if(q == 1) then
          FLAG = 2

          do j=1, Nvert
            write(2,*) (psi(i,j),i=1, Nvert)
            write(3,*) (omg(i,j),i=1, Nvert)
          end do

          write(*,*)'======= Result Saved (psi.txt , omg.txt).'
          write(4,*)'======= Result Saved (psi.txt , omg.txt).'

        end if
      end if

    end do



    stop
    end program po2d
