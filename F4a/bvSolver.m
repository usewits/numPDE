global N = 201;    %Number of points in time
global M = 200;    %Number of points in space

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
A = spalloc(2*M*(N+1), 2*M*(N+1), 9*M*(N-1)+3*M          +5*(M*N)+10); %use (y, x) or (row, col) (sparse!)
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

function result = v_initial(x)     %Initial value of v (model 2)
    %Using the approximation provided on the webpage
    global M;
    k = 2*(M+1);             %Higher k is more accurate
    h = x/k;
    som = 0;
    for j=0:k-1
        som += (     u_initial( x - (j-1) * h )
                 - 2*u_initial( x -   j   * h )
                 +   u_initial( x - (j+1) * h ) ) * ( sqrt(j+1) - sqrt(j) );
    endfor
    result = 2 / (sqrt(pi)*h^(3/2)) * som;
endfunction

function result = u_final(i)    %Fitted functions found in mathematica
    result = 0.31314927132111237*exp(-3*(-2.486678585705514+0.025199383606127345*i)^2);
endfunction

function result = v_final(i)    %Fitted functions found in mathematica
    result = -0.36948628705439973*exp(-0.0021796983589234373*(-93.06891803678627+i)^2)+0.40559736585200146*exp(-0.003584576255722463*(-78.67776776214579+i)^2);
endfunction

function save_u_v(w)
    global N;
    global M;
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
    b(pos) = u_final(i);

    pos = vpos(i, N);                   %final condition v; we choose this 0
    A(pos, pos) = 1;
    b(pos) = v_final(i);

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
