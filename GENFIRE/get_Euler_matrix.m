function Matrix = get_Euler_matrix(Angles)

phi = Angles(1);
theta = Angles(2);
psi = Angles(3);

c_ph = cosd(phi);
s_ph = sind(phi);
c_th = cosd(theta);
s_th = sind(theta);
c_ps = cosd(psi);
s_ps = sind(psi);

c_ps_c_th = c_ps*c_th;
s_ps_c_ph = s_ps*c_ph;
s_ps_s_ph = s_ps*s_ph;

Matrix = [c_ps_c_th*c_ph-s_ps_s_ph , -s_ps_c_ph*c_th-c_ps*s_ph , s_th*c_ph;
          c_ps_c_th*s_ph+s_ps_c_ph , -s_ps_s_ph*c_th+c_ps*c_ph , s_th*s_ph;
         -c_ps*s_th                ,  s_ps*s_th                , c_th    ];

end