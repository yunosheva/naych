clear();

h = 1
dt = 1
n = 10
m = 20
teta = 1/3
u = zeros(2,n,m)
tau = 1


//задаем начальную плотность
rho = zeros(n,m)
for i = 1 : n/2
    for j = 1 : m
        rho(i,j) = 1
    end
end
for i = n/2 : n
    for j = 1 : m
        rho(i,j) = 2
    end
end


// задаем равновесную функцию, sp = скалярное произведение
function F = f(sp,u2,w,rho)
    F = w*rho*(1 + sp/teta + (sp)^2/(2*teta^2) - u2/2/teta)
endfunction

// решеточное уравнение Больцмана
function F_1 = F(f_k,f_e)
    F_1 = f_k + (f_e - f_k)/tau
endfunction


// задаем вектор возможных скоростей частиц (c_0 = c(1) и т.д.)
c = zeros(9,2)
for k = 1 : 4
    c(k+1,1) = cos((k-1)*%pi/2)
    c(k+1,2) = sin((k-1)*%pi/2)
end
for k = 5 : 8
    c(k+1,1) = sqrt(2)*cos((k-1/2)*%pi/2)
    c(k+1,2) = sqrt(2)*sin((k-1/2)*%pi/2)
end

// задаем коэффициенты w_k (w_0 = w(1) и т.д.)
w = zeros(9)
w(1) = 4/9
for i = 2:5
    w(i) = 1/9
end
for i = 6:9
    w(i) = 1/36
end

// посчитаем начальные одночастичные функции распределения f_k 
// и равновесные функции f_e
f_e = zeros(9,n,m)
f_k = zeros(9,n,m)
buf = zeros(9,n,m)

for i = 1 : n
    for j = 1 : m
        for k = 1 : 9
            f_e(k,i,j) = f(c(k,1)*u(1,i,j) + c(k,2)*u(2,i,j), u(1,i,j)^2 + u(2,i,j)^2, w(k), rho(i,j))
        end
    end
end 
f_k = f_e

M = 50
for j = 1 : M
    // учитываем движение
    buf = f_k
    
    for j = 1 : m
        f_k(2,1,j) = buf(2,n,j)
    end
    for j = 1 : n
        f_k(3,j,1) = buf(3,j,m)
    end
    for j = 1 : m
        f_k(4,n,j) = buf(4,1,j)
    end
    for j = 1 : n
        f_k(5,j,m) = buf(5,j,1)
    end

    f_k(6,1,1) = buf(6,n,m)
    for j = 2 : m
        f_k(6,1,j) = buf(6,n,j-1)
    end
    for j = 2 : n
        f_k(6,j,1) = buf(6,j-1,m)
    end

    f_k(7,n,1) = buf(7,1,m)
    for j = 2 : m
        f_k(7,n,j) = buf(7,1,j-1)
    end
    for j = 1 : n-1
        f_k(7,j,1) = buf(7,j+1,m)
    end

    f_k(8,n,m) = buf(8,1,1)
    for j = 1 : m-1
        f_k(8,n,j) = buf(8,1,j+1)
    end
    for j = 1 : n-1
        f_k(8,j,m) = buf(8,j+1,1)
    end

    f_k(9,1,m) = buf(9,n,1)
    for j = 1 : m-1
        f_k(9,1,j) = buf(9,n,j+1)
    end
    for j = 2 : n
        f_k(9,j,m) = buf(9,j-1,1)
    end

    for i = 2 : n-1
        for j = 2 : m-1 
            f_k(2,i,j) = buf(2, i-1, j)
            f_k(3,i,j) = buf(3, i, j-1)
            f_k(4,i,j) = buf(4, i+1, j)
            f_k(5,i,j) = buf(5, i, j+1)
            f_k(6,i,j) = buf(6, i-1, j-1)
            f_k(7,i,j) = buf(7, i+1, j-1)
            f_k(8,i,j) = buf(8, i+1, j+1)
            f_k(9,i,j) = buf(9, i-1, j+1)
        end
    end
    sum = 0
    // посчитаем новую плотность
    for i = 1 : n
        for j = 1 : m
            rho(i,j) = f_k(1,i,j) + f_k(2,i,j) + f_k(3,i,j) + f_k(4,i,j) + f_k(5,i,j) + f_k(6,i,j) + f_k(7,i,j) + f_k(8,i,j) + f_k(9,i,j)
            sum = sum + rho(i,j)
        end
    end
    disp(sum)
    // посчитаем новую скорость вещества в узле u
    ux = 0
    uy = 0
    u = zeros(2,n,m)
    for i = 1 : n
        for j = 1 : m
            for k = 1 : 9
                u(1,i,j) = u(1,i,j) + f_k(k,i,j)*c(k,1)/rho(i,j)
                ux = ux + u(1,i,j)
                u(2,i,j) = u(2,i,j) + f_k(k,i,j)*c(k,2)/rho(i,j)
                uy = uy + u(2,i,j)
            end
        end
    end
    disp(ux, uy)

    // посчитаем новые равновесные функции
    for i = 1 : n
        for j = 1 : m
            for k = 1 : 9
                f_e(k,i,j) = f(c(k,1)*u(1,i,j) + c(k,2)*u(2,i,j), u(1,i,j)^2 + u(2,i,j)^2, w(k), rho(i,j))
            end
        end
    end 


    // сделаем одну итерацию
    
    for i = 1 : n
        for j = 1 : m
            for k = 1 : 9
                F_k(k,i,j) = F(f_k(k,i,j), f_e(k,i,j))
            end
        end
    end
    
    // задаем одночастичные функции
    f_k = F_k
end

rho_x = zeros(n)
for i = 1 : n
    rho_x(i) = rho(i,7)
end


x = [0 : 1 : n-1]
y = [0 : 1: m-1]
plot2d(x,rho_x)














