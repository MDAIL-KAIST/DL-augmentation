% fillInFourierGrid_DFT_Class
% Y. Yang, UCLA Physics & Astronomy
% First version date: 2015. 04. 30.
% output parameter: rec (nx x ny x nx) IFFT of ftArray
%                   ftArray (nx x ny x nx) interpolated Fourier 3D array
%                   CWArray (nx x ny x nx) interpolated confidence weight
%                   SDArray (nx x ny x nx) interpolated distance
%
%
% input parameter: obj.InputAngles - measured projection, (nx x ny x P) array with P # of projections
%                  obj.InputAngles - 3 Euler angle array, (3 x P) array
%                  obj.interpolationCutoffDistance - threshold for acceptible distance
%                  obj.n1_oversampled (length of oversampled projection 1st dim)
%                  obj.n2_oversampled (length of oversampled projection 2nd dim) 
%                  obj.DFT_CentroSymmetricity - 0 for complex-valued reconstruction
%                                        1 for real-valued reconstruction,
%                                        centrosymmetry will be enforced to
%                                        save calculation time
%                  obj.DFT_doGPU - 0 for not using GPU
%                           1 for using GPU
%                  obj.doCTFcorrection - 0 for not doing CTF correction
%                                        1 for doing CTF correction
%                  obj.CTFparameters - CTF parameters if doing CTF correction
%                   
%
%
% Second version date: 2015. 7. 7. (YY)
% Change: 1. Now this code properly process arbitrary-sized input projections
%         and oversampling. nx_ori, ny_ori, nx, ny can be arbitrary positive
%         integer, and it does not matter if they are even or odd.
%         2. This code assumes that the pixels of original projections are
%         "square", i.e. the pixel resolution is the same in both
%         directions. This version of code does not work if pixels are rectangular.
%         3. Implemented spacial frequency dependent interpolation
%         threshold (turn on and off using ThreshType argument)
%
% Thrid version date: 2015. 7. 9. (YY)
% Change: 1. Modified dimensional convention (1st dim: x, 2nd dim: y, 
%                                            3rd dim: z)
%         2. Rotation matrix is now Z(phi)*X(theta)*Z(psi), instead of
%         previous version [Z(phi)*Y(theta)*Z(psi)].
%
%         3. Zero degree projection is in XY plane, not XZ plane as before.
%
% Fourth version date: 2016. 4. 11. (YY)
% Change: 1. Cleaned up the code a little bit
%
%         2. Made switch for centrosymmetricity, in case of complex
%         reconstruction
%
%         3. CTF correction
%
%         4. inline C function for speedup
%
% Sixth version date: 2016. 6. 26. (YY)
% Change: 1. Cleaned up the code a little bit
%
%         2. Made siwtches for CTF correction
%
%         3. Wiener filter CTF correction 
%
% Seventh version date: 2016. 8. 11. (YY)
% Change: 1. the output should be nx x ny x nx array (y is the rotation
%              axis)
%         2. fixed bug for correctly determining the cutoff sphere
%
% Eighth version date: 2016. 8. 24. (YY)
% Change: 1. C function disabled because it is actually slower
%         2. Implemented GPU capability
%
% Nineth version date: 2016. 8. 25. (YY)
% Change: 1. Made the function consistent with AJ's convention
%
% Class version date: 2016. 12. 18. (YY)
% Change: 1. Made the fuction for use in GENFIRE_Class
%         2. Separated CTF correction and centrosymmetry part as separate
%         functions
%         3. Removed confidence weight calculation and SD calculation
             
function obj = fillInFourierGrid_DFT(obj)

tic

% original projection dimensions
n1_ori = obj.Dim1;
n2_ori = obj.Dim2;
n1 = obj.n1_oversampled;
n2 = obj.n2_oversampled;

% if distance below minInvThresh, minInvThresh will be used
% this is to prevent division by zero
minInvThresh = 0.00001;

% initialize normal vectors and rotation matrices array
normVECs = zeros(size(obj.InputProjections,3),3);
rotMATs = zeros(3,3,size(obj.InputProjections,3));

phis = obj.InputAngles(:,1);
thetas = obj.InputAngles(:,2);
psis = obj.InputAngles(:,3);

% calculate rotation matrices and normal vectors from the rotation matrices
for i=1:size(obj.InputProjections,3)
    
    rotmat1 = MatrixQuaternionRot(obj.vector1,phis(i));
    rotmat2 = MatrixQuaternionRot(obj.vector2,thetas(i));   
    rotmat3 = MatrixQuaternionRot(obj.vector3,psis(i));

    rotMATs(:,:,i) =  (rotmat1*rotmat2*rotmat3);

    init_normvec = [0 0 1];
    normVECs(i,:) = squeeze(rotMATs(:,:,i))*init_normvec';
end


% initiate Fourier space indices
if obj.DFT_CentroSymmetricity
    k1 = (-1*floor(n1/2):1:0) ;
    n_k1 = floor(n1/2);    
    
    k2 = (-1*floor(n2/2):1:floor(n2/2)) ;  
    k3 = (-1*floor(n1/2):1:floor(n1/2)) ;    
    
else
    k1 = (-1*ceil((n1-1)/2):1:floor((n1-1)/2)) / n1;    
    k2 = (-1*ceil((n2-1)/2):1:floor((n2-1)/2)) / n2;    
    k3 = k1;
end


% Fourier grid
% in case of centrosymmetry, only half of k1 will be interpolated
% and centrosymmetricity will be enforced later
[~, K1, ~] = meshgrid(k2,k1,k3);

% initialize variables
FS = zeros(size(K1))*(-1); % Fourier points

Numpt = zeros(size(K1)); % array to store how many points found per Fourier point
invSumTotWeight = zeros(size(K1)); % array to store sum of weights of Fourier point

% initiate Fourier space indices
k1_ori = (-1*ceil((n1_ori-1)/2):1:floor((n1_ori-1)/2)) ;    
k2_ori = (-1*ceil((n2_ori-1)/2):1:floor((n2_ori-1)/2)) ;    

[K20, K10] = meshgrid(k2_ori,k1_ori);

if obj.DFT_doGPU
    K10G = gpuArray(K10(:));
    K20G = gpuArray(K20(:));
end
    
for p=1:size(obj.InputProjections,3)
    % current projection
    curr_proj = squeeze(obj.InputProjections(:,:,p));
        
    [K2, K1, K3] = meshgrid(k2,k1,k3);
    
    % obtain points-to-plane distance
    D = distancePointsPlane_YY([K1(:) K2(:) K3(:)]', normVECs(p,:));
    
    % find Fourier points within the threshold
    Dind = find(D < obj.interpolationCutoffDistance);
    
    %tic
    % rotate the plane to zero degree
    CP = closestpoint(normVECs(p,:)',0,[K1(Dind); K2(Dind); K3(Dind)]);

    %toc
    
    CP_plane = (squeeze(rotMATs(:,:,p)))\CP;
    
    clear KX KY KZ
    % picked closest point which is within the projection plain, x coordinate must be zero after rotation
    if sum(abs(CP_plane(3,:)) > 0.0001) > 0
        fprintf(1,'something wrong!\n');                
    end
    
    % consider Fourier points only within the resolution circle
    Gind = Dind(abs(CP_plane(1,:)) <= n1/2 & abs(CP_plane(2,:)) <= n2/2);  % good indices  
    G_CP_plane = CP_plane(:,abs(CP_plane(1,:)) <= n1/2 & abs(CP_plane(2,:)) <= n2/2 );  % good in-plane coordinates
    
    %determine the available memory in MB
    if obj.DFT_doGPU
        GPUinfo = gpuDevice();
        av_memory_size = round(GPUinfo.AvailableMemory/1000000);
    else
        
        % determine the available memory in MB
        if ispc % in case of Windows machine
            [~, SYSTEMVIEW] = memory;
            av_memory_size = SYSTEMVIEW.PhysicalMemory.Available / 1000000;
        elseif isunix && ~ismac  % in case of linux (or unix)
            [~,out]=system('cat /proc/meminfo | grep MemFree');
            av_memory_size=sscanf(out,'MemFree:          %d kB');
            av_memory_size = av_memory_size / 1000;
        else % in case of mac (I don't have mac now, to be implemented later)
            av_memory_size = 1000;
        end
    end

    memory_per_index = 40*length(curr_proj(:))/1000000;        
    
    % determine block size for vectorized calculation
    block_size = floor(av_memory_size/memory_per_index);
    if block_size < 1
        block_size = 1;
    end
    %block_size = 500;
    cutnum = floor(length(Gind)/block_size);
    cutrem = mod(length(Gind),block_size);
    
    if cutrem~=0
        cutloopnum = cutnum + 1;
    else
        cutloopnum = cutnum;
    end
    

    % loop over Fourier points within the threshold
    for i=1:cutloopnum 
        curr_indices = ((i-1)*block_size+1):(i*block_size);
        
        if i > cutnum
            curr_indices = (cutnum*block_size+1):length(Gind);
        end
         
        % CTF correction
        if obj.doCTFcorrection       
            % to be implemented as subfunction
            CTFcorr = ones(1,length(curr_indices));
        % no CTF correction
        else
            CTFcorr = ones(1,length(curr_indices));
        end               
        
        if obj.DFT_doGPU
            G_CP_plane1_GPU_n = gpuArray(G_CP_plane(1,curr_indices)/n1);
            G_CP_plane2_GPU_n = gpuArray(G_CP_plane(2,curr_indices)/n2);
            curr_proj_GPU = gpuArray(curr_proj(:));
        
        
            % DFT calculation
            FpointsG = sum(bsxfun(@times, curr_proj_GPU, exp(-1*1i*2*pi*(K10G*G_CP_plane1_GPU_n+K20G*G_CP_plane2_GPU_n))),1);

            Fpoints = gather(FpointsG);     
            Fpoints = CTFcorr.*Fpoints;
            
            clear G_CP_plane1_GPU_n G_CP_plane2_GPU_n curr_proj_GPU FpointsG     
        else
            Fpoints = CTFcorr.*sum(bsxfun(@times, curr_proj(:), exp(-1*1i*2*pi*(K10(:)*G_CP_plane(1,curr_indices)/n1+K20(:)*G_CP_plane(2,curr_indices)/n2))),1);
        end
        
        
        %weighted avearaging
        CIND = Gind(curr_indices);

        currDist = D(CIND);
                        
        currDist(currDist < minInvThresh) = minInvThresh; % if distance smaller than minInvThresh, put minInvThresh (to prevent divison by zero)
        currInvDist = 1./ currDist;           % inverse distance
        currTotWeight = currInvDist;

        % re-average inverse distance
        FS(CIND) = FS(CIND).* invSumTotWeight(CIND) + currTotWeight.*Fpoints;
        invSumTotWeight(CIND) = invSumTotWeight(CIND) + currTotWeight;
        FS(CIND) = FS(CIND) ./ invSumTotWeight(CIND);
        
        Numpt(CIND) = Numpt(CIND) + 1;     
        
        clear Fpoints
    end
    
    
    
   
end

clear Mindist invSumdist Numpt

% enforce centrosymmetricity
if obj.DFT_CentroSymmetricity 
    obj = obj.CentroSymmetricity(FS, n1, n_k1, n2);        
else
    obj.measuredK = reshape(FS,n1,n2,n1);
end

clear FS

obj.recIFFT = My_stripzero(real(my_ifft(obj.measuredK)),[obj.Dim1 obj.Dim2 obj.Dim1]);

timeTakenToFillInGrid = toc;
timeTakenToFillInGrid = round(10*timeTakenToFillInGrid)./10;
fprintf('GENFIRE: Fourier grid assembled in %.12g seconds.\n\n',timeTakenToFillInGrid);

end



function D = distancePointsPlane_YY(points, normvec)
    %distancePointsPlane_YY unsigned distances betwen 3D points and a plane
    % through origin
    %
    %   D = distancePointsPlane_YY(point, normvec)
    %   Returns the euclidean distance between points and a plane going through origin with normal vector normvec,
    %   given by: 
    %   points : (3 x n) array of 3D vectors
    %   normvec : (3 x 1) or (1 x 3) array of normal vector
    %   D     : (1 x n) vector  
 
    %
    %   ---------
    %   author : Y. Yang, UCLA Physics and Astronomy 
    %   created the 05/03/2015.
    %

    % normalized plane normal
    normvec = normvec(:) / norm(normvec);
    D = abs(points(1,:)*normvec(1) + points(2,:)*normvec(2) + points(3,:)*normvec(3));

end


function x = closestpoint(n, d, p)
    % n is the vector [A,B,C] that defines the plane
    % d is the distance of the plane from the origin
    % p is the point  [P,Q,R]
    if size(p,2) == 1
        v = (d - sum(p.*n)) / sum(n.*n);
        x = p + v * n;
    else
        nr = repmat(n,[1 size(p,2)]);
        v = (d - sum(p.*nr,1)) / sum(n.*n);
        x = p + repmat(v,[3 1]) .* nr;
    end
end


function dd = MatrixQuaternionRot(vector,theta)

theta = theta*pi/180;
vector = vector/sqrt(dot(vector,vector));
w = cos(theta/2); x = -sin(theta/2)*vector(1); y = -sin(theta/2)*vector(2); z = -sin(theta/2)*vector(3);
RotM = [1-2*y^2-2*z^2 2*x*y+2*w*z 2*x*z-2*w*y;
      2*x*y-2*w*z 1-2*x^2-2*z^2 2*y*z+2*w*x;
      2*x*z+2*w*y 2*y*z-2*w*x 1-2*x^2-2*y^2;];

dd = RotM;
end

