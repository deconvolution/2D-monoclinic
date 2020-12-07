A = randi(10,300,300,300);
B = randi(10,300,300,300);
tic;C=A+B;toc;
A = gpuArray(A);
B = gpuArray(B);
tic;D=A+B;wait(gpuDevice);toc;