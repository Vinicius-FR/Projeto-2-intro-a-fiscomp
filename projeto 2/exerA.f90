module functions

implicit none

contains

function f(n, h)
    real(8), intent(in) :: h
    integer, intent(in) :: n
    real(8) :: x, f

    x = 1.d0 / 3.d0 + n * h

    f = exp(4.d0 * x) * cos(x / 2.d0)

end function f

function df(x)
    real(8), intent(in) :: x
    real(8) :: df

    df = 4.d0 * exp(4.d0 * x) * cos(x / 2.d0) - (1.d0/2.d0) * exp(4.d0 * x) * sin(x / 2.d0)

end function df

function d2f(x)
    real(8), intent(in) :: x
    real(8) :: d2f

    d2f = (63.d0/4.d0) * exp(4.d0 * x) * cos(x / 2.d0) - 4.d0 * exp(4.d0 * x) * sin(x / 2.d0)

end function d2f

function d3f(x)
    real(8), intent(in) :: x
    real(8) :: d3f

    d3f = 61.d0 * exp(4.d0 * x) * cos(x / 2.d0) - (191.d0/8.d0) * exp(4.d0 * x) * sin(x / 2.d0)

end function d3f

end module functions

program exerA
use functions
implicit none

real(8), allocatable, dimension(:) :: h
real(8) :: dfrente_2, dtras_2, dsim_3, dsim_5, d2sim_3, d2sim_5, d3asim_5
real(8), parameter :: x = 1.d0 / 3.d0
integer :: i, n

open(1,file = 'tabA_in.dat',status = 'old',action = 'read')
read(1,*) n
allocate(h(n))

read(1,*) (h(i) , i = 1,n)

close(1)

open(2, file = 'tabA_out.dat', status='replace')

do i = 1, n

    dsim_3 = (f(1,h(i)) - f(-1,h(i))) / (2.d0 * h(i))

    dfrente_2 = (f(1,h(i)) - f(0,h(i))) / h(i)

    dtras_2 = (f(0,h(i)) - f(-1,h(i))) / h(i)

    !dsim_5 = (-f(2,h(i)) + 8 * f(1,h(i)) - 8 * f(-1,h(i) + f(-2,h(i)))) / (12 * h(i))

    d2sim_3 = (f(1,h(i)) - 2.d0 * f(0,h(i)) + f(-1,h(i))) / (h(i) ** 2.d0)

    d2sim_5 = (-f(2,h(i)) + 16.d0 * f(1,h(i)) - 30.d0 * f(0,h(i)) + 16.d0 * f(-1,h(i)) - f(-2,h(i))) / (12.d0 * h(i) ** 2.d0)

    d3asim_5 = (f(2,h(i)) - 2.d0 * f(1,h(i)) + 2.d0 * f(-1,h(i)) - f(-2,h(i))) / (2.d0 * h(i) ** 3.d0)

    write(2,*) h(i), dsim_3 - df(x), dfrente_2 - df(x), dtras_2 - df(x), d2sim_3 - d2f(x), &
        d2sim_5 - d2f(x), d3asim_5 - d3f(x)
end do

close(2)

print*, 'Baseando-se no valor de h que gera o menor desvio, tesm que os valores mais apropriados para cada caso são:'
print*, 'Derivada simétrica de 3 pontos: h = 0,000001 (12° valor)'
print*, 'Derivada para frente de 2 pontos: h = 0,00000001 (14° valor)'
print*, 'Derivada para trás de 2 pontos: h = 0,00000001 (14° valor)'
print*, 'Derivada segunda simétrica de 3 pontos: h = 0,00005 (10° valor)'
print*, 'Derivada segunda simétrica de 5 pontos: h = 0,0005 (8° valor)'
print*, 'Derivada terceira anti-simétrica de 5 pontos: h = 0,0005 (8° valor)'
print*, 'Por fim, nota-se que, conforme aumenta a ordem da derivada e o número de pontos, aumenta a relevância ',&
        'do erro computacional inerente, tornando-o cada vez mais significativo cada vez que h se aproxima da ',&
        'ordem de grandeza da precisão. Portanto nem sempre o menor h é o mais adequado.'

end program