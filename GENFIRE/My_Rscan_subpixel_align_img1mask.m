function [d1, d2, minVal] = My_Rscan_subpixel_align_img1mask(img1,img2,RscanRad,Rres)

    if Rres >= 0.5
        disp('Resolution should be less than 0.5');
        return;
    end

    if size(img1)~=size(img2)
        disp('two images should be in the same size!');
        return;
    end


    radR = ceil(RscanRad);

    currMinR = inf;
    for i=-radR:radR
        for j=-radR:radR
            shift_img2 = circshift(img2,[i j]);
            shift_img2(img1==0) = 0;
            currR = calcR_norm_YY(img1,shift_img2);
            if currR < currMinR
                dI1 = i;
                dI2 = j;
                currMinR = currR;
            end
        end
    end

    newimg2 = circshift(img2,[dI1 dI2]);


    DivFact = 4;

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
        for j=1:length(curr_ind1arr)
            for k=1:length(curr_ind2arr)
                shift_img2 = real(My_FourierShift(newimg2,curr_ind1arr(j),curr_ind2arr(k)));
                shift_img2(img1==0) = 0;
                minRar(j,k) = calcR_norm_YY(img1,shift_img2);
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
    
    minVal = min(minRar(:));

end


function img2 = My_FourierShift(img,dy,dx)

ny = size(img,1); nx = size(img,2);
[X Y] = meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
F = My_IFFTN(img);
Pfactor = exp(2*pi*1i*(dx*X/nx + dy*Y/ny));
img2 = My_FFTN(F.*Pfactor);

end

function R = calcR_norm_YY(data1,data2)

    data1 = data1(:);
    data2 = data2(:);
    
    if length(data1)~=length(data2)
        fprintf(1,'data length does not match!\n');
        R = -1;
        return
    end
    
    lscale = dot(data1,data2)/norm(data2)^2;
    
    data2 = data2 * lscale;
    
    R = sum(abs(abs(data1)-abs(data2)))/sum(abs(data1));
    
end