function [R, new_x_center, new_y_center] = alignByRFactor_M(ref_img,img,window_half_size,Rscanres,FSind)
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
    
    
    R = R / sum(ref_img(:));
    new_x_center = ncx - Dx;
    new_y_center = ncy - Dy;  
end