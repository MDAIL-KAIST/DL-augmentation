function Rand_pos_list=Func_generating_pos_3D_rand(volsize, res, minimum_distance, density, getting_percent, scanning_unit) 
% J.Lee, KAIST, 2020

    %scanning_unit = 4;
    scanning_volume_size = volsize/scanning_unit;
    pre_i = 0;

    %
    i=1;
    array_list = [];
    array_list(:,1) = res*scanning_volume_size*rand(3,1);
    
    for i1 =1:scanning_unit
        for i2 =1:scanning_unit
            for i3 =1:scanning_unit

                unit_endflag = 0;
                while ~unit_endflag
                    generating_pos = res*scanning_volume_size*[i1-1;i2-1;i3-1]+ res*scanning_volume_size*rand(3,1);

                    if all(sum((array_list-generating_pos).^2).^0.5 > minimum_distance)
                        array_list(:,i)=generating_pos;
                        i=i+1;
                    end

                    curr_density = abs(i-pre_i) / (res*scanning_volume_size)^3;

                    if curr_density > getting_percent*density
                        unit_endflag = 1;
                        pre_i=i;
                    end
                end

                %[i1,i2,i3];
            end
        end
    end
    
    Rand_pos_list=array_list;


end
