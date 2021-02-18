function CMline = get_Wolfgangcommonline_from_angle(PJ,Angle)   
    Angle = -1*Angle;
    wt = wtlib(1);

    a_90 = round(Angle/90);
    a_res = Angle - a_90*90;
    
    if a_90 ~=0
      cPJ = My_stripzero(PJ,floor((size(PJ)-1)/2)*2 + 1);
      rotcPJ = rot90(cPJ,a_90);
      cPJ = My_paddzero(rotcPJ,size(PJ));
      rotImg = wt.rotateimagedeg(cPJ,a_res);
    else
      rotImg = wt.rotateimagedeg(PJ,a_res);
      
      
    end
    CMline = squeeze(sum(rotImg,1));
    
    
end