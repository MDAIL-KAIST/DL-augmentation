%CentroSymmetricity
function obj = CentroSymmetricity(obj, FS, n1, n_k1, n2)
    obj.measuredK = single(zeros(n1,n2,n1));
    
    if mod(n2,2) == 1

      if mod(n1,2) == 1
          ftArray_half = reshape(FS,n_k1+1,n2,n1);        
          obj.measuredK(1:n_k1+1,:,:) = ftArray_half;

          ftArray_ahalf = ftArray_half(1:end-1,:,:);            
          ftArray_ahalf = fliplr(ftArray_ahalf);
          ftArray_ahalf = permute(flipud(permute(ftArray_ahalf, [3 1 2])),[2 3 1]);
          ftArray_ahalf = flipud(ftArray_ahalf);
          obj.measuredK(n_k1+2:end,:,:) = conj(ftArray_ahalf);

          clear ftArray_half ftArray_ahalf

      else
          ftArray_half_ext = reshape(FS,n_k1+1,n2,n1+1);   
          ftArray_half = ftArray_half_ext(:,:,1:end-1);
          obj.measuredK(1:n_k1+1,:,:) = ftArray_half;

          ftArray_ahalf_ext = ftArray_half_ext(2:end-1,:,:);            
          ftArray_ahalf_ext = fliplr(ftArray_ahalf_ext);
          ftArray_ahalf_ext = permute(flipud(permute(ftArray_ahalf_ext, [3 1 2])),[2 3 1]);
          ftArray_ahalf_ext = flipud(ftArray_ahalf_ext);
          obj.measuredK(n_k1+2:end,:,:) = conj(ftArray_ahalf_ext(:,:,1:end-1));

          clear ftArray_half ftArray_half_ext ftArray_ahalf_ext

      end

  else


      if mod(n1,2) == 1
          ftArray_half_ext = reshape(FS,n_k1+1,n2+1,n1);
          ftArray_half = ftArray_half_ext(:,1:end-1,:);
          obj.measuredK(1:n_k1+1,:,:) = ftArray_half;

          ftArray_ahalf_ext = ftArray_half_ext(1:end-1,:,:);
          ftArray_ahalf_ext = fliplr(ftArray_ahalf_ext);
          ftArray_ahalf_ext = permute(flipud(permute(ftArray_ahalf_ext, [3 1 2])),[2 3 1]);
          ftArray_ahalf_ext = flipud(ftArray_ahalf_ext);

          obj.measuredK(n_k1+2:end,:,:) = conj(ftArray_ahalf_ext(:,1:end-1,:));

          clear ftArray_half ftArray_half_ext ftArray_ahalf_ext



      else
          ftArray_half_ext = reshape(FS,n_k1+1,n2+1,n1+1);
          ftArray_half = ftArray_half_ext(:,1:end-1,1:end-1);
          obj.measuredK(1:n_k1+1,:,:) = ftArray_half;

          ftArray_ahalf_ext = ftArray_half_ext(2:end-1,:,:);   
          ftArray_ahalf_ext = fliplr(ftArray_ahalf_ext);
          ftArray_ahalf_ext = permute(flipud(permute(ftArray_ahalf_ext, [3 1 2])),[2 3 1]);
          ftArray_ahalf_ext = flipud(ftArray_ahalf_ext);

          obj.measuredK(n_k1+2:end,:,:) = conj(ftArray_ahalf_ext(:,1:end-1,1:end-1));

          clear ftArray_half ftArray_half_ext ftArray_ahalf_ext


      end


    end
end
        
        
        
        
        
   