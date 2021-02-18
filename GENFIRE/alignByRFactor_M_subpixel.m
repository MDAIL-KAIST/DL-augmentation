function [R, new_x_center, new_y_center] = alignByRFactor_M_subpixel(ref_img,img,window_half_size,Rscanres,FSind)
    ncx = round((size(ref_img,1)+1)/2);
    ncy = round((size(ref_img,2)+1)/2);
    radR = ceil(window_half_size);

    bestR = inf;    
    for sx=-radR:radR
        for sy=-radR:radR
            Shiftimg = circshift(img,[sx sy]);
            R = sum(sum(abs(ref_img(FSind)-Shiftimg(FSind))));
            if R < bestR
                bestR = R;
                Dx = sx;
                Dy = sy;
            end
        end
    end
    

    % subpixel alignment
    if Rscanres < 0.5      
      newimg = circshift(img,[Dx Dy]);

        DivFact = 5;

        NumIter = abs(floor(log(Rscanres)/log(DivFact)));
        curr_start1 = -0.5;
        curr_end1 = 0.5;
        curr_start2 = -0.5;
        curr_end2 = 0.5;

        nx = size(img,1); ny = size(img,2);
        [Ygrid, Xgrid] = meshgrid(-ceil((nx-1)/2):floor((nx-1)/2),-ceil((ny-1)/2):floor((ny-1)/2));
        Xgrid = ifftshift(Xgrid);
        Ygrid = ifftshift(Ygrid);
        

        for i=1:NumIter

            curr_boundary1arr = curr_start1:(curr_end1-curr_start1)/DivFact:curr_end1;
            curr_boundary2arr = curr_start2:(curr_end2-curr_start2)/DivFact:curr_end2;

            curr_ind1arr = (circshift(curr_boundary1arr,[0 -1]) + curr_boundary1arr)/2;
            curr_ind2arr = (circshift(curr_boundary2arr,[0 -1]) + curr_boundary2arr)/2;

            curr_ind1arr = curr_ind1arr(1:end-1);
            curr_ind2arr = curr_ind2arr(1:end-1);

           
            bestR = inf; 
            for j=1:length(curr_ind1arr)
                for k=1:length(curr_ind2arr)
                    
                    shift_img2 = real(My_FourierShift_full(newimg, nx, ny, Xgrid, Ygrid, curr_ind1arr(j),curr_ind2arr(k)));
                    R = sum(sum(abs(ref_img(FSind)-shift_img2(FSind))));
                    if R < bestR
                        bestR = R;
                        minj = j;
                        mink = k;
                    end
                end
            end
            
            curr_start1 = curr_boundary1arr(minj);
            curr_end1 = curr_boundary1arr(minj+1);
            curr_start2 = curr_boundary2arr(mink);
            curr_end2 = curr_boundary2arr(mink+1);    
        end

        sdx = curr_ind1arr(minj);
        sdy = curr_ind2arr(mink);

        Dx_f = Dx+sdx;
        Dy_f = Dy+sdy;


    else
      Dx_f = Dx;
      Dy_f = Dy;
    end
    
    R = bestR / sum(ref_img(:));
    new_x_center = ncx + Dx_f;
    new_y_center = ncy + Dy_f;
end



function img2 = My_FourierShift_full(img, ny, nx, Ygrid, Xgrid, dy,dx)

img2 = fftshift(fftn(ifftn(ifftshift(img)).*exp(2*pi*1i*(dx*Xgrid/nx + dy*Ygrid/ny))));

end

