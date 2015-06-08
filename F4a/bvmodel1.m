global N = 9;    %Number of points in time
global M = 100;    %Number of points in space

t_max = 0.3;    %Bounds. t_min = x_min = 0
x_max = 1;

dt = t_max / N;     %delta t and delta x
dx = x_max / M;

%w = zeros(1:2*M*(N+1),1); We want to know this; Aw=b
b = zeros(2*M*(N+1),1);
A = spalloc(2*M*(N+1), 2*M*(N+1), 9*M*(N-1)+3*M          +(M*N)+10); %use (y, x) or (row, col) (sparse!)
%A = zeros(2*M*(N+1),2*M*(N+1)); %use (y, x) or (row, col) (full!)

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
    result = exp(-2*pi^(3/2)*t)*pi^(3/2)*(
                    cos(2*pi*(sqrt(pi)*t+x)) - sin(2*pi*(sqrt(pi)*t+x))
                    - exp(4*pi^(3/2)*t)*(cos(2*pi^(3/2)*t-2*pi*x)+sin(2*pi^(3/2)*t-2*pi*x))
                );
endfunction

function save_u_v(w)
    global N
    global M
    filenameu = "u.dat";
    filenamev = "v.dat";
    fidu = fopen (filenameu, "w");
    fidv = fopen (filenamev, "w");
    for n = 0:N
        for i = 0:M-1
            if i != 0
                fprintf (fidu, " ");
                fprintf (fidv, " ");
            endif
            fprintf (fidu, "%f", w(upos(i, n)));
            fprintf (fidv, "%f", w(vpos(i, n)));
        endfor
        fprintf (fidu, "\n")
        fprintf (fidv, "\n")
    endfor
    fclose (fidu);
    fclose (fidv);
endfunction

disp("Initialising matrix A...");

time_begin = time();

for i = 0:M-1                       %Fill the matrix

    x = dx * i;

    for n = 1:N-1
        row = upos(i, n);               %Leapfrog
        A(row, upos(i, n+1)) =  1;
        A(row, upos(i, n-1)) = -1;
        A(row, vpos(i, n)  ) = -2*dt;
        
        row = vpos(i, n);               %SC1 (+leapfrog)
        A(row, vpos(i, n+1)) =  1;
        A(row, vpos(i, n-1)) = -1;
        A(row, upos(i-2, n)) =    dt/(dx^3);
        A(row, upos(i-1, n)) = -2*dt/(dx^3);
        A(row, upos(i+1, n)) =  2*dt/(dx^3);
        A(row, upos(i+2, n)) =   -dt/(dx^3);
    endfor


    pos = upos(i, 0);                   %initial condition u
    A(pos, pos) = 1;
    b(pos) = sin(2*pi*x);

    pos = vpos(i, 0);                   %initial condition v
    A(pos, pos) = 1;
    b(pos) = 0; %does not do anything, added for clarity

    pos = upos(i, N);                   %final condition u
    A(pos, pos) = 1;
    b(pos) = u_exact(x, t_max);

    pos = vpos(i, N);                   %final condition v
    A(pos, pos) = 1;
    b(pos) = v_exact(x, t_max);

endfor

time_end = time();

printf("Initialisation done! (%f seconds)\n", time_end - time_begin);

%spy(A);
%pause();
printf("Solving system of size %d\n", 2*M*(N+1));
time_begin = time();

w = A \ b;                          %Solve the system!

time_end = time();
printf("Solved system! (%f seconds)\n", time_end - time_begin)

save_u_v(w)
