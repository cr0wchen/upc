%                              第一个程序模板

% for solving two-point boundary problem -d2u/dx2+xu=(1+x)sinx on [0,pi]
% with u(0)=u(pi)=0
% true solution u=sinx
    
    a = 0;
    b = pi;
    N = 1000;
    h = (b-a)/N;
    x = a+[1:N-1]*h;
    e = ones(N-1,1);
    A = spdiags([e,-2*e,e],[-1,0,1],N-1,N-1);
    A = -1*A/(h^2);
    fun1 = inline('x');
    matrix_q = diag(feval(fun1,x),0);
    fun2 = inline('(1+x).*sin(x)');
    vector_f = feval(fun2,x);
    vector_f = vector_f';
    u = (A+matrix_q)\vector_f;
    %plot(x,u)
    
    hold on;
    %plot(x,sin(x),'--r')
    
    plot(x,u'-sin(x))

    