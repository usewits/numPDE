global N = 9;    %Number of points in time
global M = 100;    %Number of points in space

t_min = 0;
t_max = 1;    %Bounds
x_min = -5;
x_max = 8;

dt = t_max / N;     %delta t and delta x
dx = (x_max-x_min) / (M+1);
%Let x_i = -5 + (i+1)*dx
%We have x_0 = -5 + dx, and x_M-1 = 8 - dx
%We can impose boundary conditions for x_-1=-5 and x_M=8

%w = zeros(1:2*M*(N+1),1); We want to know this; Aw=b
b = zeros(2*M*(N+1),1);
A = spalloc(2*M*(N+1), 2*M*(N+1), 9*M*(N-1)+3*M          +(M*N)+10); %use (y, x) or (row, col) (sparse!)
%A = zeros(2*M*(N+1),2*M*(N+1)); %use (y, x) or (row, col) (full!)

function result = upos(i, n)                        %Easy to use versions of get_uv_pos_periodic
    result = get_uv_pos_boundary(i, n, false);
endfunction

function result = vpos(i, n)                        %v(i*dx, n*dt) = w(vpos(i, n))
    result = get_uv_pos_boundary(i, n, true);
endfunction

function result = get_uv_pos_boundary(i, n, is_v)   %Notice the boundary conditions!
    global N;
    global M;
    if i < 0 || i >= M
        result = -1     %i out of bounds!
    endif
    result = n*2*M + i*2 + 1;
    if is_v
        result += 1;
    endif
endfunction

function result = u_initial(x)     %Initial value of u (model 2)
    result = exp(-3*(x-2)^2);
endfunction

function result = integrantv0(x)    %The integrant of v0(x). Was computed using mathematica.
    result = (6*exp(-3*(-2+x)^2)*(23+6*(-4+x)*x))/(sqrt(pi)*sqrt(tmp-x));
endfunction


function result = v_initial(x)     %Initial value of v (model 2)
    fun = @(x,c) (6*exp(-3*(-2+x)^2)*(23+6*(-4+x)*x))/(sqrt(pi)*sqrt(c-x));
    %TODO: replace this with a good value!
    tmp = x+0.1;
    result = quad(@(x)fun(x,tmp),0,x);
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

    x = x_min + (i+1)*dx;

    for n = 1:N-1
        row = upos(i, n);               %Leapfrog
        A(row, upos(i, n+1)) =  1;
        A(row, upos(i, n-1)) = -1;
        A(row, vpos(i, n)  ) = -2*dt;
        
        row = vpos(i, n);               %SC1 (+leapfrog)
        A(row, vpos(i, n+1)) =  1;
        A(row, vpos(i, n-1)) = -1;

%v(-5, t) = v(8, t) = 0 (chosen BC)
%u(-5, t) = u(8, t) = 0 (given BC)
        if i-2 >= 0
            A(row, upos(i-2, n)) =    dt/(dx^3);
        endif
        if i-1 >= 0
            A(row, upos(i-1, n)) = -2*dt/(dx^3);
        endif
        if i+1 <= M-1
            A(row, upos(i+1, n)) =  2*dt/(dx^3);
        endif
        if i+2 <= M-1
            A(row, upos(i+2, n)) =   -dt/(dx^3);
        endif
    endfor


    pos = upos(i, 0);                   %initial condition u (given)
    A(pos, pos) = 1;
    b(pos) = u_initial(x);

    pos = vpos(i, 0);                   %initial condition v (given)
    A(pos, pos) = 1;
    b(pos) = v_initial(x); 

    pos = upos(i, N);                   %final condition u; we choose this 0
    A(pos, pos) = 1;
    b(pos) = 0;

    pos = vpos(i, N);                   %final condition v; we choose this 0
    A(pos, pos) = 1;
    b(pos) = 0;

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
