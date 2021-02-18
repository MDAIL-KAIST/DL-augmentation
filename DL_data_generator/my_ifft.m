function realout = my_ifft(k)
realout =fftshift((ifftn(ifftshift(k))));
end