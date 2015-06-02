global N = 5;    %Number of points in time
global M = 5;    %Number of points in space

t_max = 0.3;    %Bounds. t_min = x_min = 0
x_max = 1;

dt = t_max / N;     %delta t and delta x
dx = x_max / M;

%w = zeros(1:2*M*(N+1),1); We want to know this; Aw=b
b = zeros(1:2*M*(N+1),1);
A = spalloc(2*M*(N+1), 2*M*(N+1), 9*M*(N-1)+3*M); %use (y, x) or (row, col)

function result = upos(i, n)                        %Easy to use versions of get_uv_pos_periodic
    result = get_uv_pos_periodic(i, n, false);
endfunction

function result = vpos(i, n)                        %v(i*dx, n*dt) = w(vpos(i, n))
    result = get_uv_pos_periodic(i, n, true);
endfunction

function result = get_uv_pos_periodic(i, n, is_v)   %Notice peridic boundary conditions!
    global N;
    global M;
    i_hat = mod(i+M, M);
    result = n*2*M + i_hat*2 + 1;
    if is_v
        result += 1;
    endif
endfunction

function result = u_exact(x, t)     %Exact solution of PDE (model 1)
    result = 1/2*(-exp(2*pi^(3/2)*t)*sin(2*pi^(3/2)*t-2*pi*x)+exp(-2*pi^(3/2)*t)*sin(2*pi^(3/2)*t+2*pi*x));
endfunction

function result = v_exact(x, t)     %obtained using mathematica (derivative to t of u_exact)
    result = exp(-2*pi^(3/2)*t)*pi^(3/2)*(cos(2*pi*(sqrt(pi)*t+x))-sin(2*pi*(sqrt(pi)*t+x))-exp(4*pi^(3/2)*t)*(cos(2*pi^(3/2)*t-2*pi*x)+sin(2*pi^(3/2)*t-2*pi*x)));
endfunction

disp("Initialising matrix A...");

time_begin = time();

for i = 0:M-1                       %Fill the matrix

    for n = 1:N-1
        col = upos(i, n);               %Leapfrog
        A(upos(i, n+1) , col) =  1;
        A(upos(i, n-1) , col) = -1;
        A(vpos(i, n)   , col) = -2*dt;
        
        col = vpos(i, n);               %SC1 (+leapfrog)
        A(vpos(i, n+1) , col) =  1;
        A(vpos(i, n-1) , col) = -1;
        A(upos(i-2, n) , col) =    dt/(dx^3);
        A(upos(i-1, n) , col) = -2*dt/(dx^3);
        A(upos(i+1, n) , col) =  2*dt/(dx^3);
        A(upos(i+2, n) , col) =   -dt/(dx^3);
    endfor

    
    pos = upos(i, 0);                   %initial condition u
    A(pos, pos) = 1;
    b(pos) = sin(2*pi*dx*i);

    pos = vpos(i, 0);                   %initial condition v
    A(pos, pos) = 1;
    b(pos) = 0; %does not do anything, added for clarity

    pos = upos(i, N);                   %final condition u
    A(pos, pos) = 1;
    b(pos) = u_exact(dx*i, t_max);
    
    pos = vpos(i, N);                   %final condition v
    A(pos, pos) = 1;
    b(pos) = v_exact(dx*i, t_max);

endfor

time_end = time();

printf("Initialisation done! (%f seconds)\n", time_end - time_begin);

spy(A);
pause();
printf("Solving system of size %d\n", 2*M*(N+1));
time_begin = time();

w = A / b;                          %Solve the system!

time_end = time();
printf("Solved system! (%f seconds)\n", time_end - time_begin)

