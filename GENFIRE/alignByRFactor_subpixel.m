function [R, new_x_center, new_y_center] = alignByRFactor_subpixel(obj,ref_img,img,window_half_size)

    radR = ceil(window_half_size);

    bestR = inf;
    for sx=-radR:radR
        for sy=-radR:radR
            Shiftimg = circshift(img,[sx sy]);
            R = sum(sum(abs(ref_img(obj.FSind)-Shiftimg(obj.FSind))));
            if R < bestR
                bestR = R;
                new_x_center = ncx - sx;
                new_y_center = ncy - sy;
            end
        end
    end

    
    newimg2 = circshift(img2,[dI1 dI2]);

    % subpixel alignment
    if Rres < 0.5

        DivFact = 5;

        NumIter = abs(floor(log(Rres)/log(DivFact)));
        curr_start1 = -0.5;
        curr_end1 = 0.5;
        curr_start2 = -0.5;
        curr_end2 = 0.5;


        for i=1:NumIter

            curr_boundary1arr = curr_start1:(curr_end1-curr_start1)/DivFact:curr_end1;
            curr_boundary2arr = curr_start2:(curr_end2-curr_start2)/DivFact:curr_end2;

            curr_ind1arr = (circshift(curr_boundary1arr,[0 -1]) + curr_boundary1arr)/2;
            curr_ind2arr = (circshift(curr_boundary2arr,[0 -1]) + curr_boundary2arr)/2;

            curr_ind1arr = curr_ind1arr(1:end-1);
            curr_ind2arr = curr_ind2arr(1:end-1);


            minRar = zeros(length(curr_ind1arr),length(curr_ind2arr));
            if doBG==1
                minBGar = zeros(length(curr_ind1arr),length(curr_ind2arr));
            end
            for j=1:length(curr_ind1arr)
                for k=1:length(curr_ind2arr)
                    
                    shift_img2 = real(My_FourierShift(newimg2,curr_ind1arr(j),curr_ind2arr(k)));
                    
                    if MaskDir == 2
                        comp_img1 = img1(:,Maskrange);
                        comp_img2 = shift_img2(:,Maskrange);
                    elseif MaskDir == 1
                        comp_img1 = img1(Maskrange,:);
                        comp_img2 = shift_img2(Maskrange,:);
                    end
                                 
                    comp_img2(comp_img1==0) = 0;
                    
                    if doBG==1
                        [cR, c_minBG] = fitBG_calcR_YY_PosProj(comp_img1,comp_img2,doNorm);
                    else
                        [cR] = calcR_YY(comp_img1,comp_img2,doNorm);
                    end


                    minRar(j,k) = cR;
                    if doBG==1
                        minBGar(j,k) = c_minBG;
                    end
                end
            end

            [~, minInd] = min(minRar(:));
            [minj, mink] = ind2sub(size(minRar),minInd);

            curr_start1 = curr_boundary1arr(minj);
            curr_end1 = curr_boundary1arr(minj+1);
            curr_start2 = curr_boundary2arr(mink);
            curr_end2 = curr_boundary2arr(mink+1);    
        end

        sd1 = curr_ind1arr(minj);
        sd2 = curr_ind2arr(mink);

        d1 = dI1+sd1;
        d2 = dI2+sd2;

        if doBG==1
            FinalminBG = minBGar(minj,mink);
        end

        minVal = min(minRar(:));
    else
        
        if doBG == 1
            FinalminBG = currMinBG;
        end
        minVal = currMinR;
        d1 = dI1;
        d2 = dI2;
    end
end





[R, new_x_center, new_y_center] = alignByRFactor(ref_img,img,window_half_size)

if all(size(ref_img)~=size(img))
    error('Images must be same size')
end

[dimx, dimy] = size(ref_img);
ncx = round((dimx+1)/2);
ncy = round((dimy+1)/2);

if nargin<3
    window_half_size = floor(dimx/10);
end

bestR = 1e30;
new_x_center = ncx;
new_y_center = ncy;
norm_factor = sum(abs(ref_img(:)));
for sx = -window_half_size:window_half_size
    for sy = -window_half_size:window_half_size
        R = sum(sum(abs(ref_img-circshift(img,[sx, sy])))) ./ norm_factor;
        if R < bestR
           bestR = R;
           new_x_center = ncx - sx;
           new_y_center = ncy - sy;
        end
    end
end

end


function img2 = My_FourierShift_full(img,dy,dx)

ny = size(img,1); nx = size(img,2);
[X Y] = meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
F = My_IFFTN(img);
Pfactor = exp(2*pi*1i*(dx*X/nx + dy*Y/ny));
img2 = My_FFTN(F.*Pfactor);

end

function R = calcR_YY(data1,data2,doNorm)

    data1 = data1(:);
    data2 = data2(:);
    
    if length(data1)~=length(data2)
        fprintf(1,'data length does not match!\n');
        R = -1;
        return
    end
    
    if doNorm
        lscale = dot(data1,data2)/norm(data2)^2;

        data2 = data2 * lscale;
    end
    
    R = sum(abs(abs(data1)-abs(data2)));
    
end


function [R, minBG] = fitBG_calcR_YY_PosProj(data1,data2,doNorm)

    data1 = data1(:);
    data2 = data2(:);
    
    if doNorm
        lscale = dot(data1,data2)/norm(data2)^2;

        data2 = data2 * lscale;
    end
    
    opts = optimset('Display','off');
    % fitting for x axis projection
    init_guess = 0;
    fun = @(p,xdata) sum(abs(abs(   xdata + (xdata>0)*1.*p -  ((xdata + (xdata>0)*1.*p)<0).*(xdata + (xdata>0)*1.*p -1e-6) )-abs(data1)));
    [p, fminres, ~, eflag] = lsqcurvefit(fun, init_guess, data2, 0, [], [], opts);
    fitresult = data2 + (data2>0)*p;  
    R = calcR_YY(data1,fitresult,0);    
    minBG = p;
end