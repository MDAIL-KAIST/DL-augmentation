function intep_vol = Func_real_intepolation_OV_sampling(Vol,OV)
% J.Lee, KAIST, 2020

    volsize=length(Vol);
    %
    OV_low_limit = round((OV-1)/2);
    OV_high_limit = round((OV-2)/2);
    %
    [X,Y,Z] = meshgrid(1:volsize,1:volsize,1:volsize);
    OV_list = 1-OV_low_limit/OV:1/OV:volsize+OV_high_limit/OV;
    [Xq,Yq,Zq] = meshgrid(OV_list, OV_list, OV_list);
    
    %OV_vol = interp3(X,Y,Z,Vol,Xq,Yq,Zq,'linear');
    OV_vol = interp3(X,Y,Z,Vol,Xq,Yq,Zq,'cubic');
    OV_vol(isnan(OV_vol))=0;
    
    %binning
    intep_vol=imresize3(OV_vol,1/OV, 'box');

end