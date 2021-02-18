wt = wtlib(1);
A = zeros(128,128);
A(65,65) = 1;
% A(64,65) = 1;
% A(66,65) = 1;
% A(65,64) = 1;
% A(65,66) = 1;
% A(63,65) = 1;
% A(67,65) = 1;
% A(65,63) = 1;
% A(65,67) = 1;
Img1 = A;
Img2 = My_rotate_image_perf_test(A,pi/6-0.001);
Img2 = wt.rotateimagedeg(A,30);
%Img2 = My_rotate_image_perf(A,pi/4-0.001);
%Img2 = imrotate(A,45,'crop');
subplot(1,2,1)
imagesc(Img1);axis image
subplot(1,2,2)
imagesc(Img2);axis image