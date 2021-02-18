% WTip class with static methods for Fourier based image shif, shear, rotation
% up- and down sampling...
% Only verified for SQUARE images with EVEN number of pixels per dimension
% (very likely nc and nr, etc are randomly permutated in the code at the moment)
%
% Wolfgang Theis, August 2013
%
% Changes:
%   March 2014: doing integer part of shift with circshift
%   March 2014: converted all static to non-static methods



classdef wtlib < handle
    properties
        emdfilename =' ';
        tomosetup = [];
        tomo = [];
        tomo_set = [];
        tomofromsetup = [];
        mag = 1;
        data_dim = [];
        nrep = 1;
        ndefocus = 1;
        nrset =  0;
        image_dim = [];
        pix = 1024;
        res = 0;
        dstr = 'date';
        ta = [];
        philist =[];
        nsetofphilist = [];
        nrprojections = 0;
        isp = [];
        ta_set = [];
        isp_set = [];
        scanrotation = 6.5;
        scandistortion = eye(2);
        rotsummary = [];
        scalebarlength = 5;
        wto;
    end
    
    methods
        
        function obj = wtlib(emd)
            if (ischar(emd))
                obj.emdfilename = emd;
                datagroup = '/data/tomoraw/';
                obj.mag = h5readatt(emd,datagroup,'magnification');
                obj.data_dim = size(h5read(emd, strcat(datagroup, '/data'), [1 1 1 1 1],[Inf Inf Inf 1 1 ]));
                obj.nrep =  obj.data_dim(2);
                obj.ndefocus =  obj.data_dim(1);
                obj.nrset =  obj.data_dim(3);
                obj.image_dim = size( h5read(emd, strcat(datagroup, '/data'), [1 1 1 1 1],[1 1 1 Inf Inf]));
                obj.pix = obj.image_dim(end);
                obj.ta = h5read(emd, strcat(datagroup, '/tiltangles'), [1 1],[Inf Inf]);
                [tu, ia] = unique(obj.ta(:,1),'last');
                [tus,itu] = sort(tu);
                obj.philist = tus;
                obj.nsetofphilist = ia(itu);
                obj.nrprojections = size(obj.philist,1);
                obj.isp = h5read(emd, strcat(datagroup, '/imagesetparameters'), [1 1],[Inf Inf]);
                obj.ta_set = obj.ta(obj.nsetofphilist,:);
                obj.isp_set = obj.isp(obj.nsetofphilist,:);
                t = obj.isp(1,7);
                dn = t/86400.0 + datenum('01-Jan-1970');
                obj.dstr = datestr(double(dn));
                a = h5read(emd, '/data/tomoraw/dim1', [1],[2]);
                obj.res = a(2)-a(1);
            end
        end
        
        function tomo = tomosetupfromparams(obj, tomosetup, dir, alphafactor,alphaaxistobeamtilt, gammaaxistobeamtilt)
            obj.tomosetup = tomosetup;
            wto = wtiptomotilting();
            wto.alphafactor = alphafactor;
            wto.alphaaxistobeamtilt = alphaaxistobeamtilt; 
            wto.gammaaxistobeamtilt = gammaaxistobeamtilt;
            wto.CalcOffsetSlope( tomosetup(1,1), tomosetup(1,2), tomosetup(2,1),  tomosetup(2,2),  tomosetup(3,1),  tomosetup(3,2))
            wto.calcgeneralgammatiltangles(16,tomosetup(4,1), tomosetup(4,2), tomosetup(5,1),  tomosetup(5,2),dir)
            tomo =wto.tiltanglelist';
            obj.tomo = tomo;
            obj.tomo_set = tomo(obj.ta_set(:,1)+1,:);
            if (alphafactor == 1)
                if max(abs(obj.tomo_set(:,2)-obj.ta_set(:,3)))> 1e-1
                    warning('calculation of alpha angles does not match script''s');
                end
                if (max(abs(obj.tomo_set(:,3)-obj.ta_set(:,4)))> 1e-1)
                    warning('calculation of gamma angles does not match script''s');
                end
            end
            obj.wto = wto;
        end
        
        function tomo = tomosetupfromemd(obj, dir, alphafactor)
            tomosetup = h5read(obj.emdfilename, '/tomosetup/alignangles', [1 1],[Inf Inf]);
            tomo = obj.tomosetupfromparams(tomosetup, dir, alphafactor);
        end
        
        function sd = calcscandist(obj, imnr, vlen, lat_type)
            s = obj.gettomoimage(obj.emdfilename, imnr, round(obj.ndefocus/2),  round(obj.nrep/2));
            mfft = obj.calcmaskedfft(s, 0.15);
            vset = obj.findbraggminmaxk(mfft, 6, vlen-10, vlen+10);
            vnr1 = 1;
            vnr2 = 3;
            ang = vset(vnr2,4)-vset(vnr1,4);
            c_ang = round(ang/60)*60;
            v2c = [[cosd(c_ang), -sind(c_ang)]; [sind(c_ang), cosd(c_ang)]] *vset(vnr1,1:2)';
            vexpset = [vset(vnr1,1:2)',vset(vnr2,1:2)'];
            vcorset = [vset(vnr1,1:2)',v2c];
            m = vcorset*(vexpset^(-1));
            r  =atan2(m(2,1),m(1,1))
            mr = [[cos(r), sin(r)]; [-sin(r), cos(r)]];
            m = mr*m;
            sd = m/m(1,1);
            obj.scandistortion = sd;
        end
        
        function sc = image_correctdistortion(obj,s)
            s = obj.softenalledges(s, 40);
            sp = obj.addpadding(s);
            sp = obj.applytransf2(sp,obj.scandistortion);
            sc = obj.strippadding(sp);
        end
        
        function sc = fft_correctdistortion(obj,s)
            s = obj.softenalledges(s, 40);
            sp = obj.addpadding(s);
            sp = obj.applytransf2(sp,(obj.scandistortion^(-1))');
            sc = obj.strippadding(sp);
        end
        
        function vc = fftvectors_correctdistortion(obj,v)
            vc = obj.scandistortion*v;
        end
        
        function vc = fftvector_set_correctdistortion(obj,vs)
            vcp = (obj.scandistortion)*(vs(:,1:2)')/det(obj.scandistortion);
            nmax = size(vcp,2);
            vc = zeros(size(vs));
            vc(:,1) = vcp(1,:);
            vc(:,2) = vcp(2,:);
            vc(:,3) = sqrt(vc(:,1).*vc(:,1)+vc(:,2).*vc(:,2));
            vc(:,4) = atan2(vc(:,2),vc(:,1)) *180/pi;
            vc(:,5) = vs(:,5);
        end
        
        function calcrotation(obj, corrimage)
            rotsummary = zeros(size(obj.philist,1),2);
            for philistindex = 1:size(obj.philist,1)
                s = obj.gettomoimage(obj.emdfilename, obj.nsetofphilist(philistindex), round(obj.ndefocus/2),  round(obj.nrep/2));
                if corrimage
                    s = obj.image_correctdistortion(s);
                end
                rot = -obj.tomo_set(philistindex, 4);
                st = s(1:obj.pix/2, obj.pix/4:obj.pix/4+obj.pix/2-1);
                mfft = obj.calcmaskedfft(st, 0.45);
                %mfft= fftshift(fft2(st));
                vset = obj.findbraggmink(mfft, 10, 50);
                % find peak with rot+5.5 +-10 degree orientation
                if ~corrimage
                    vset = obj.fftvector_set_correctdistortion(vset);
                end
                knr = find(abs(vset(1:end,4)-rot-obj.scanrotation) <10, 1, 'first');
                rotangledeg = vset(knr,4)-90.0;
                rotsummary(philistindex,1:2) = [rot-90+obj.scanrotation, rotangledeg];
            end
            obj.rotsummary = rotsummary;
        end
        
        function recalc_rotation(obj)
            obj.rotsummary(:,1) =  -obj.tomo_set(:, 4)-90+obj.scanrotation;
        end
      
        function showplots(obj)
            figure(1)
            clf
            plot(obj.ta_set(:,2),obj.rotsummary-obj.scanrotation, 'o')
            xlabel('phi/degrees')
            ylabel('rot/degrees')
            title({'Measured and expected' , 'image rotation angle'})
            legend('Calculated', 'Measured');
            
            figure(2)
            clf
            da_set = obj.isp_set(:,4) - obj.ta_set(:,3);
            dg_set = obj.isp_set(:,5) - obj.ta_set(:,4);
            plot(da_set, 'r-')
            hold on
            plot(dg_set, 'b-')
            legend('\Delta \alpha', '\Delta \gamma');
            xlabel('projection number')
            ylabel('delta angle/degrees')
            title('Difference between actual and intended angles')
            hold off
            
            figure(3)
            clf
            % show close to 90 deg image
            [a ind] = min(abs(obj.ta_set(:,2)-90));
            [a ind] = min(abs(obj.philist-obj.ta_set(ind, 1)));
            n = obj.nsetofphilist(ind);
            image90 = obj.gettomoimage(obj.emdfilename, n, round(obj.ndefocus/2), round(obj.nrep/2));
            imshow(obj.addscalebar(image90, obj.res,5),[])
            title('Closest image to phi = 90')
            
            figure(4)
            clf
            image = obj.gettomoimage(obj.emdfilename, 1, round(obj.ndefocus/2), round(obj.nrep/2));
            imshow(obj.addscalebar(image, obj.res,5),[])
            title('first image')
            
            
        end
            
        function h = hanning(obj, n)
            h = (1:n)';
            h = 0.5*(1-cos(2*pi*h/n));
        end;
        
        function m = rotdegmat(obj, phideg)
            m = obj.rotmat(phideg*pi/180);
        end;
        
        function m = rotmat(obj, phi)
            m = [cos(phi) -sin(phi); sin(phi) cos(phi)];
        end;
        
        % utility function for input from EMD files
        function image = getimage(obj, ima)
            s = size(ima);
            image  = reshape(ima, s(end-1), s(end));
        end
        
        function im = gettomoimage(obj, emd, nrset, nrdefocus, nrrep)
            a = h5read(emd, '/data/tomoraw/data', [nrdefocus, nrrep, nrset 1 1],[1 1 1 Inf Inf]);
            im = obj.getimage(a);
        end
        
        function image = readfromh5(obj, emd, dataset, start, amount)
            a = h5read(emd, dataset, start, amount);
            image = getimage(a);
        end
        
        function im = cleartopffts(obj, image)
            f = fftshift(fft2(image));
            % to ensure ifft generates real data (no imaginary components)
            f(1,:) = 0;
            f(:,1) = 0;
            im = ifft2(ifftshift(f));
        end
        
        
        % back could be improved in quality and speed
        function back = getbackgroundvalue(obj, image)
            [nr nc] = size(image);
            backlist = sort(image(:));
            back = backlist(round(nr*nc/128.0));
        end
        
        function [back, st] = getbackground(obj, image, n)
            [nr nc] = size(image);
            ba = image(nr-n:nr,nc-n:nc);
            back = mean(ba(:));
            st = std(ba(:));
        end
        
        function sp =addpadding(obj, image)
            [nr, nc] = size(image);
            back = obj.getbackgroundvalue(image);
            sp = zeros(size(image)*2) + back;
            sp(nr-nr/2:nr-nr/2+nr-1,nc-nc/2:nc-nc/2+nc-1) = image;
        end
        
        function sp =addpaddingwithnoise(obj, image, n)
            [nr, nc] = size(image);
            [back, st] = obj.getbackground(image, n);
            if st>back
                st = back*0.9;
            end;
            sp = rand(size(image)*2)*st + back;
            sp(nr-nr/2:nr-nr/2+nr-1,nc-nc/2:nc-nc/2+nc-1) = image;
        end
        
        function si = strippadding(obj, image)
            [nr nc] = size(image);
            nr = nr/2;
            nc = nc/2;
            si = image(nr/2:nr/2+nr-1,nc/2:nc/2+nc-1);
        end
        
        function  f = calcmaskedfft(obj, image, pixrange)
            f = obj.calcoffsetmaskedfft(image, 0, 0, pixrange);
        end
        
        function [d, a, f, mask] =  calcplanespacingangle(obj, image, pcx, pcy, pixrange, dguess, angleguess)
            pcx =round(pcx);
            pcy =round(pcy);
            n = size(image,2);
            back = obj.getbackgroundvalue(image);
            mask = zeros(n);
            mask(round(n/2+pcy-pixrange*n):round(n/2+pcy+pixrange*n),round(n/2+pcx-pixrange*n):round(n/2+pcx+pixrange*n)) = obj.hanning(2*round(pixrange*n)+1) *obj.hanning(2*round(pixrange*n)+1)';
            f = (image-back).*mask;
            fs =  sum(f(:));
            ms = sum(mask(:));
            f = f-fs/ms*mask;
            fs =  sum(abs(f(:)))
            d = dguess;
            a = angleguess;
            v =1:1:n;
            [x, y] = meshgrid(v,v);
            
            g = sqrt(fs)/2*cos((cos(a)*y + sin(a)*x)/d*2*pi);
            g = g.*mask;
            f(:,1:500) =g(:,1:500);
            fs = sum(f(:))
            
        end
            
        function  [f mask] = calcoffsetmaskedfft(obj, image, pcx, pcy, pixrange)
            pcx =round(pcx);
            pcy =round(pcy);
            n = size(image,2);
            back = obj.getbackgroundvalue(image);
            mask = zeros(n);
            mask(round(n/2+pcy-pixrange*n):round(n/2+pcy+pixrange*n),round(n/2+pcx-pixrange*n):round(n/2+pcx+pixrange*n)) = obj.hanning(2*round(pixrange*n)+1) *obj.hanning(2*round(pixrange*n)+1)';
            f = fftshift(fft2((image-back).*mask));
        end
        
        function v = findbragg(obj, fftim, nmax)
            nr = size(fftim,2);
            af = abs(fftim);
            %find peaks
            p = af;
            p((circshift(af,[-1,0]) > af) | (circshift(af,[1,0]) > af) | (circshift(af,[0,-1])> af) | (circshift(af,[0,1]) > af)) = 0;
            [row,col,v] = find(p);
            [w, ind] = sort(v,'descend');
            v = zeros(nmax,5);
            %logint = zeros(nmax,1);
            for k = 1:nmax
                r = row(ind(k));
                c = col(ind(k));
                cf = c + (log(af(r,c+1)) - log(af(r,c-1)))/(-2*log(af(r,c+1)) + 4*  log(af(r,c)) -2* log(af(r,c-1)));
                rf = r + (log(af(r+1,c)) - log(af(r-1,c)))/(-2*log(af(r+1,c)) + 4*  log(af(r,c)) -2* log(af(r-1,c)));
                v(k,1) = cf-nr/2-1;
                v(k,2) = -(rf-nr/2-1);
                v(k,3) = sqrt(v(k,1)*v(k,1)+v(k,2)*v(k,2));
                v(k,4) = atan2(v(k,2),v(k,1)) *180/pi;
                v(k,5) = log(af(r,c));
            end
        end
        
        function v = findbraggmink(obj, fftim, nmax, mink)
            v = obj.findbraggminmaxk(fftim, nmax, mink, 10000);
        end
        
        function v = findbraggminmaxk(obj, fftim, nmax, mink, maxk)
            nr = size(fftim,2);
            af = abs(fftim);
            %find peaks
            p = af;
            p((circshift(af,[-1,0]) > af) | (circshift(af,[1,0]) > af) | (circshift(af,[0,-1])> af) | (circshift(af,[0,1]) > af)) = 0;
            [x,y] =meshgrid((-nr/2):(-nr/2+nr-1));
            p(x.*x+y.*y < mink*mink)=0;
            p(x.*x+y.*y > maxk*maxk)=0;
            p(x.*x+y.*y >= (nr-4)*(nr-4)/4)=0; % to facilitate c+-1, r+-1 below
            [row,col,v] = find(p);
            [w, ind] = sort(v,'descend');
            v = zeros(nmax,5);
            %logint = zeros(nmax,1);
            for k = 1:nmax
                r = row(ind(k));
                c = col(ind(k));
                cf = c + (log(af(r,c+1)) - log(af(r,c-1)))/(-2*log(af(r,c+1)) + 4*  log(af(r,c)) -2* log(af(r,c-1)));
                rf = r + (log(af(r+1,c)) - log(af(r-1,c)))/(-2*log(af(r+1,c)) + 4*  log(af(r,c)) -2* log(af(r-1,c)));
                v(k,1) = cf-nr/2-1;
                v(k,2) = -(rf-nr/2-1);
                v(k,3) = sqrt(v(k,1)*v(k,1)+v(k,2)*v(k,2));
                v(k,4) = atan2(v(k,2),v(k,1)) *180/pi;
                v(k,5) = log(af(r,c));
            end
        end
        
        function [sh f] = shearimage (obj, image, sheardim, shfactor)
            % only shear if shfactor not zero or tiny
            if abs(shfactor)>1e-6
                [nr nc] = size(image);
                n = size(image, sheardim);
                assert(mod(nc,2) == 0 && mod(nc,2)==0, 'shearimage is only implemented for EVEN sizes');
                [x, y] =meshgrid(-nr/2:(nr/2-1),-nc/2:(nc/2-1));
                f = fftshift(fft(image,[],sheardim),sheardim);
                f = f .* exp(-2.0i*pi/n .*x .*y * shfactor);
                % to ensure ifft generates real data (no imaginary components)
                f(1,:) = 0;
                f(:,1) = 0;
                sh = ifft(ifftshift(f,sheardim),[],sheardim);
                assert(max(abs(imag(sh(:))))<1e-6, 'shearimage: sheared image has imagenary components');
                % all imaginary components should already be zero within numerical accuracy
                sh = real(sh);
            else
                f = fftshift(fft2(image));
                sh=image;
            end
        end
        
        function [sh f] = shiftimage (obj, image, shift)
            intshift = round(shift);
            image = circshift(image, [intshift(1), intshift(2)]);
            shift = shift - intshift;
            % only shear if shift not zero or tiny
            if max(abs(shift))>1e-3
                [nr, nc] = size(image);
                assert(mod(nc,2) == 0 && mod(nc,2)==0, 'shiftimage is only implemented for EVEN sizes');
                [x, y] =meshgrid(((-nr/2):(nr/2-1)),((-nc/2):(nc/2-1)));
                % circshift for closest integer
                f = fftshift(fft2(image));
                f = f .* exp(-2.0i*pi/nr .* x*shift(2) -2.0i*pi/nc* y * shift(1));
                % to ensure ifft generates real data (no imaginary components)
                f(1,:) = 0;
                f(:,1) = 0;
                sh = ifft2(ifftshift(f));
                assert(max(abs(imag(sh(:))))<1e-6, 'shiftimage: shifted image has imaginary components');
                % all imaginary components should already be zero within numerical accuracy
                sh = real(sh);
            else
                f = fftshift(fft2(image));
                sh = image;
            end
        end
        
        function [a, b, c, d] = sheardecomposition(obj, m)
            mstart = m;
            d = det(m); % if <0 will fail below
            assert(d>0, 'function sheardecomposion requires det(m)>0');
            d = sqrt(d);
            m = m/d;
            if max(max(abs(m -eye(2)))) < 1e-6
                a = 0;
                b = 0;
                c = 0;
            else
                b = m(2,1); % will fail below if b = 0
                assert(not (b==0), 'function sheardecomposion requires m(2,1) != 0');
                a = (m(1,1)-1)/b;
                c = (m(2,2)-1)/b;
                % m = sa*sb*sc
            end
            assert (max(max(abs(mstart - d*[[1 a];[0 1]]* [[1 0];[b 1]]* [[1 c];[0 1]])))<1e-6, 'matrix decomposition error');
        end
        
        function sh = applytransf (obj, image, mat)
            [a, b, c, d] = obj.sheardecomposition(mat);
            sh = obj.shearimage (image, 1, c);
            sh = obj.shearimage (sh, 2, b);
            sh = obj.shearimage (sh, 1, a);
            sh = obj.scaleimage (sh, 1, d);
            sh = obj.scaleimage (sh, 2, d);
        end
        
        function [matA, matB, matC] = matdecompose(obj, mat)
            [lmat umat] = lu(mat);
            matA = diag(diag(umat));
            matB = umat/diag(diag(umat));
            matC = lmat;
            assert (max(max(abs(mat - matC*matB*matA)))<1e-6, 'matrix decomposition error');
        end
        
        function ima = applytransf2(obj, image, mat)
            sp = obj.cleartopffts(image);
            [matA, matB, matC] = obj.matdecompose(mat);
            sp =obj.scaleimage(sp,2,matA(1,1));
            sp =obj.scaleimage(sp,1,matA(2,2));
            sp = obj.shearimage (sp, 1, matB(1,2));
            ima = obj.shearimage (sp, 2, matC(2,1));
        end
        
        function sh = rotateimagedeg (obj, image, rot)
            [a, b, c, d] = obj.sheardecomposition(obj.rotdegmat(rot));
            assert(abs(d-1) <1e-6, 'function sheardegrot: rotation matrix should have det(m) == 1');
            sh = obj.shearimage (image, 1, c);
            sh = obj.shearimage (sh, 2, b);
            sh = obj.shearimage (sh, 1, a);
        end
        
        function sh = rotateimage (obj, image, rot)
            [a, b, c, d] = obj.sheardecomposition(obj.rotmat(rot));
            assert(abs(d-1) <1e-6, 'function sheardegrot: rotation matrix should have det(m) == 1');
            sh = obj.shearimage (image, 1, c);
            sh = obj.shearimage (sh, 2, b);
            sh = obj.shearimage (sh, 1, a);
        end
        
        function sc = scaleimage(obj, image, dim, scfactor)
            assert(scfactor > 0, 'scalepixels requires scalefactor > 0');
            if scfactor > 1.0
                sc = obj.shrinkpixels (image, dim, 1.0/scfactor);
            else
                if scfactor < 1.0
                    sc = obj.expandpixels(image, dim, 1.0/scfactor);
                else
                    sc = image;
                end
            end
        end
        
        function sc = expandpixels (obj, image, dim, scfactor)
            assert(scfactor >1, 'expandpixels requires scalefactor > 1');
            n = size(image,dim);
            if dim == 1
                myimage = image';
            else
                myimage = image;
            end
            nscaled = round(n*scfactor/2.0)*2;
            poff = nscaled/2-n/2;
            back = obj.getbackgroundvalue(image);
            padimage = zeros(nscaled,n)+ back;
            padimage((1:n)+poff, 1:n) = myimage;
            f = fftshift(fft2(padimage));
            fcroped = f(1+poff:n+poff,1:n);
            % set fcroped center row and column to zero to ensure real ifft
            fcroped(1,:) =0;
            fcroped(:,1) =0;
            sc =ifft2(ifftshift(fcroped));
            assert(max(abs(imag(sc(:))))<1e-6, 'expandpixels: pixelexpanded image has imaginary components');
            sc = real(sc);
            if dim == 1
                sc = sc';
            end
        end
        
        function sc = shrinkpixels (obj, image, dim, scfactor)
            assert(scfactor <1, 'shrinkpixels requires 0 < scalefactor < 1');
            n = size(image,dim);
            if dim == 1
                myimage = image';
            else
                myimage = image;
            end
            f = fftshift(fft2(myimage));
            f(1, :) = 0;
            f(:, 1) = 0;
            nscaled = round(n/scfactor/2.0)*2;
            poff = nscaled/2-n/2;
            padfft = zeros(nscaled,n);
            padfft((1:n)+poff, 1:n) = f;
            padfft(1, :) = 0; % not needed but here for simplicity and debugging
            padfft(:, 1) = 0;
            im = ifft2(ifftshift(padfft));
            sc = im(1+poff:n+poff,1:n);
            assert(max(abs(imag(sc(:))))<1e-6, 'shrinkpixels: pixelexpanded image has imaginary components');
            sc = real(sc);
            if dim == 1
                sc = sc';
            end
        end
        
        function image =split(obj, im1, im2)
            image = im1;
            [nx ny] = size(im1);
            image(1:nx/2,1:ny/2) = im2(1:nx/2,1:ny/2);
            image(nx/2:nx,ny/2:ny) = im2(nx/2:nx,ny/2:ny);
        end
        
        function [m, minv] = affinefrompair(obj, vset, vcorset)
            m = vcorset * (vset^(-1)); % vcorset = m*vset
            minv = vset * (vcorset^(-1)); % vset = minv*vcorset
        end
        
        function [mrsp, mrspinv] = realspaceaffinefrompair(obj, vset, vcorset)
            %m = vcorset * (vset^(-1)); % vcorset = m*vset
            %minv = vset * (vcorset^(-1)); % vset = minv*vcorset
            mrsp = vcorset^(-1)*vset;
            mrspinv = vset^(-1)*vcorset;
        end
        
        function ima = softenedge(obj, image, width, nr)
            % nr = 0:top, 1:left, ...
            im = rot90(image,nr);
            %assert (nr == 0, 'softenedge only implemented for top');
            m = ones(size(im));
            h = obj.hanning(2*width) *ones(1,size(im,1));
            m(1:width,:) = h(1:width,:);
            back = obj.getbackgroundvalue(im);
            ima = (im-back).*m +back;
            ima = rot90(ima,4-nr);
        end
        
        function ima = redopadding (obj, image, width)
            s = obj.strippadding(image);
            s = obj.softenedge(s, 40, 0);
            s = obj.softenedge(s, 40, 1);
            s = obj.softenedge(s, 40, 2);
            s = obj.softenedge(s, 40, 3);
            ima=obj.addpadding(s);
        end
        
        function ima = softenedgerect(obj, image, width, fromtop, back)
            m = ones(size(image));
            m(1:fromtop,:) = 0;
            h = obj.hanning(2*width) *ones(1,size(image,2));
            m(1+fromtop:width+fromtop,:) = h(1:width,:);
            ima = (image-back).*m +back;
        end
        
        function ima = softenalledges(obj, im, width)
            % nr = 0:top, 1:left, ...
            m = ones(size(im));
            h = obj.hanning(2*width) *ones(1,size(im,2));
            v = obj.hanning(2*width) *ones(1,size(im,1));
            w = size(m,2);
            m(1:width,:) = m(1:width,:) .* h(1:width,:);
            m(w:-1:w-width+1,:) =  m(w:-1:w-width+1,:) .* h(1:width,:);
            m(:, 1:width) = m(:, 1:width) .* v(1:width,:)';
            m(:, w:-1:w-width+1) =  m(:, w:-1:w-width+1)  .* v(1:width,:)';
            %m(1:width,:) = h(1:width,:);
            back = obj.getbackgroundvalue(im);
            ima = (im-back).*m +back;
        end
        
        function ima = multhanning(obj, image, side)
            n = size(image,1)/2;
            ha = zeros(2*n,1);
            ha(1:1.5*n) = obj.hanning(1.5*n);
            ha = ones(2*n,1)*ha';
            ha = rot90(ha,side);
            ima = image.*ha;
        end
        
        function ima = multslopedhanning(obj, image, side)
            n = size(image,1)/2;
            ha = zeros(2*n,1);
            ha(1:1.5*n) = obj.hanning(1.5*n).*power(2*n-(1:1.5*n)',2);
            ha = ones(2*n,1)*ha';
            ha = rot90(ha,side);
            mi = min(image(:));
            ima = (image-mi).*ha;
        end
       
        function ima = multslopedhanningdepth(obj, image, side, depth)
            d = round(depth);
            n = size(image,1)/2;
            ha = zeros(2*n,1);
            ha(1:d) = obj.hanning(d).*power(2*n-(1:d)',2);
            ha = ones(2*n,1)*ha';
            ha = rot90(ha,side);
            mi = min(image(:));
            ima = (image-mi).*ha;
        end
       
        function shfft(obj, image,base)
            imshow(log(abs(fftshift(fft2(image)))+base),[])
        end
        
        function s = testimage(obj, n, gp)
            assert(mod(n,2)==0, 'testimage: n has to be even');
            v = (-n/2:-n/2+n-1);
            [x, y] =  meshgrid(v,v);
            s =zeros(n);
            s((abs(x)<50)&(abs(y)<50)) = 100;
            s((abs(x)<20)&(abs(y)<10)) = 0;
            g = gausswin(gp);
            s = conv2(s,g*g', 'same');
        end
        
        function s = testlatticeimage(obj, n, na, gp)
            assert(mod(n,2)==0, 'testimage: n has to be even');
            v = (-n/2:-n/2+n-1);
            [x, y] =  meshgrid(v,v);
            s =zeros(n);
            s((abs(x)<n/4)&(abs(y)<n/4) & mod(x,na) == 0 & mod(y,na) == 0 ) = 100;
            g = gausswin(gp);
            s = conv2(s,g*g', 'same');
        end
        
        function s = testpointimage(obj, n, gp)
            assert(mod(n,2)==0, 'testimage: n has to be even');
            s =zeros(n);
            s(n/2+1,n/2+1) = 1000;
            g = gausswin(gp);
            s = conv2(s,g*g', 'same');
        end
        
        function showtests(obj)
            figure (1);
            im = obj.testimage(256,5);
            im = obj.cleartopffts(im);
            subplot(3,2,1);
            imshow(im,[]);
            title('original');
            subplot(3,2,2);
            imshow(obj.shiftimage(im,[10,50]),[]);
            title('shift [10,50]');
            subplot(3,2,3);
            imshow(obj.shearimage(im,1, 0.1),[]);
            title('shear dim = 1, factor = 0.1');
            subplot(3,2,4);
            imshow(obj.rotateimagedeg(im,15),[]);
            title('rotate 15 deg');
            subplot(3,2,5);
            imshow(obj.scalepixels(im,1,1.2),[]);
            title('scalepix dim = 1, factor = 1.2');
            subplot(3,2,6);
            imshow(obj.scalepixels(im,2,0.8),[]);
            title('scalepix dim = 2, factor = 0.8');
            
        end
        
        function im = addscalebar(obj, im, res, len)
            m =max(im(:));
            n = size(im);
            off = round(n(1)/40);
            p = round(len/res);
            im(n(1)-2*off: n(1)-off,n(2)-off-p:n(2)-off) = m;
        end
        
        function textlines = loadasciifile(obj,fname)
            fid = fopen(fname,'r');
            n = 1;
            tline = fgets(fid);
            while ischar(tline)
                textlines{n} = tline;
                n = n+1;
                tline = fgets(fid);
            end
            fclose(fid);
        end
       
        function ret = writeasciifile(obj,fname, textlines)
            fid = fopen(fname,'w');
            for n=1:size(textlines,2);
                fwrite(textlines{n});
            end
            fclose(fid);
        end

        function emdsummary(obj, emd, calcrot)
            %%
            datagroup ='/data/tomoraw';
            mag = h5readatt(emd,datagroup,'magnification')
            data_dim = size(h5read(emd, strcat(datagroup, '/data'), [1 1 1 1 1],[Inf Inf Inf 1 1 ]));
            nrep =  data_dim(2)
            ndefocus =  data_dim(1)
            nrset =  data_dim(3)
            image_dim = size( h5read(emd, strcat(datagroup, '/data'), [1 1 1 1 1],[1 1 1 Inf Inf]));
            pix = image_dim(end)
            ta = h5read(emd, strcat(datagroup, '/tiltangles'), [1 1],[Inf Inf]);
            [tu, ia] = unique(ta(:,1),'last');
            [tus,itu] = sort(tu);
            philist = tus;
            nsetofphilist = ia(itu);
            projections = size(philist,1)
            even_pr = philist(mod(philist,2) == 0);
            odd_pr = philist(mod(philist,2) == 1);
            
            % check angles
            isp = h5read(emd, strcat(datagroup, '/imagesetparameters'), [1 1],[Inf Inf]);
            ta_set = ta(nsetofphilist,:);
            isp_set = isp(nsetofphilist,:);
            da_set = isp_set(:,4) - ta_set(:,3);
            dg_set = isp_set(:,5) - ta_set(:,4);
            valid_set  = true(size(da_set));
            valid_set(abs(da_set)>3) = false;
            da_valid = da_set(valid_set,:);
            dg_valid = dg_set(valid_set,:);
            mean_da_valid = mean(da_valid)
            std_da_valid = std(da_valid)
            mean_dg_valid = mean(dg_valid)
            std_dg_valid = std(dg_valid)
            
            
            figure(1)
            clf
            plot(da_valid, 'r-')
            hold on
            plot(dg_valid, 'b-')
            legend('\Delta \alpha', '\Delta \gamma');
            xlabel('projection number')
            ylabel('delta angle/degrees')
            title('Difference between actual and intended angles')
            hold off
            
            % show first image
            image = obj.gettomoimage(emd, 1, round(ndefocus/2), round(nrep/2));
            imshow(image,[])
            % show close to 90 deg image
            [a ind] = min(abs(ta_set(:,2)-90))
            ta_set(ind, :)
            [a ind] = min(abs(philist-ta_set(ind, 1)))
            n = nsetofphilist(ind)
            image90 = obj.gettomoimage(emd, n, round(ndefocus/2), round(nrep/2));
            imshow(image90,[])
            
            %get expected image rotation
            
            %tomo_set = tomo(nsetofphilist,:);
            % unfortunately rotation axis angle was not written to early emd files - now
            % fixed.
            % calculate here from tomosetup:
            
            tomosetup = h5read(emd, '/tomosetup/alignangles', [1 1],[Inf Inf]);
            wto = wtiptomotilting();
            wto.CalcOffsetSlope( tomosetup(1,1), tomosetup(1,2), tomosetup(2,1),  tomosetup(2,2),  tomosetup(3,1),  tomosetup(3,2))
            wto.calcgeneralgammatiltangles(16,tomosetup(4,1), tomosetup(4,2), tomosetup(5,1),  tomosetup(5,2),1)
            tomo =wto.tiltanglelist';
            tomo_set = tomo(ta_set(:,1)+1,:);
            mean_tomota_a_diff = mean(tomo_set(:,2)-ta_set(:,3))
            std_tomota_a_diff = std(tomo_set(:,2)-ta_set(:,3)-  mean_tomota_a_diff )
            mean_tomota_g_diff = mean(tomo_set(:,3)-ta_set(:,4))
            std_tomota_g_diff = std(tomo_set(:,3)-ta_set(:,4)-  mean_tomota_g_diff )
            plot(ta_set(:,2),tomo_set(:,4), 'o')
            xlabel('phi/degrees')
            ylabel('rot/degrees')
            title('image rotation vs phi');
            
            if(calcrot)
                for philistindex = 1:size(philist,1)
                    s = obj.gettomoimage(emd, nsetofphilist(philistindex), round(ndefocus/2),  round(nrep/2));
                    rot = -tomo_set(philistindex, 4);
                    st = s(1:pix/2, pix/4:pix/4+pix/2-1);
                    mfft = obj.calcmaskedfft(st, 0.45);
                    %mfft= fftshift(fft2(st));
                    vset = obj.findbraggmink(mfft, 10, 50);
                    % find peak with rot+5.5 +-10 degree orientation
                    knr = find(abs(vset(1:end,4)-rot-5.5) <10, 1, 'first');
                    rotangledeg = vset(knr,4)-90.0;
                    rotsum(philistindex,1:2) = [rot-90, rotangledeg];
                end
                % put everything on a page for printing
                figure(2)
                clf
                set(2, 'PaperType', 'A4');
                subplot(3,2,3)
                set(text,'Interpreter', 'none')
                txt_emd = texlabel( emd,'literal');
                annotation('textbox', [.2 .88 .1 .1], 'String', emd);
                t = isp(1,7);
                dn = t/86400.0 + datenum('01-Jan-1970');
                dstr = datestr(double(dn));
                a = h5read(emd, '/data/tomoraw/dim1', [1],[2]);
                res = a(2)-a(1);
                txt = {['recorded: ', dstr], ['Mag =  ', num2str(mag/1e6), ' M'], ['pix =  ', num2str(pix)], ['resol = ',num2str(res), ' nm'], ['proj =  ', num2str(nrset)], ...
                    ['rep =  ', num2str(nrep)],  ['Nr defocus =  ', num2str(ndefocus)] };
                annotation('textbox', [.15 .8 .1 .1], 'String', txt);
                subplot(3,2,3)
                plot(rotsum(:,1),'b')
                hold on
                plot(rotsum(:,2)-5.8,'ro')
                xlabel('projection number')
                ylabel('image rotation angle/degrees')
                title({'Measured and expected' , 'image rotation angle'})
                legend('Calculated', 'Measured');
                subplot(3,2,4)
                plot(rotsum(:,2)-5.8-rotsum(:,1),'ro')
                xlabel('projection number')
                ylabel('\Delta image rotation angle/degrees')
                title({'Difference between measured and ','expected image rotation angle'})
                
                set(2, 'PaperType', 'A4');
                print([emd,'sum_p2.png'],'-dpng','-r200')
                
            end
            
            
            
            %%
            % put everything on a page for printing
            figure(1)
            set(1, 'PaperType', 'A4');
            clf
            subplot(3,2,2)
            set(text,'Interpreter', 'none')
            txt_emd = texlabel( emd,'literal');
            annotation('textbox', [.2 .88 .1 .1], 'String', emd);
            t = isp(1,7);
            dn = t/86400.0 + datenum('01-Jan-1970');
            dstr = datestr(double(dn));
            a = h5read(emd, '/data/tomoraw/dim1', [1],[2]);
            res = a(2)-a(1);
            txt = {['recorded: ', dstr], ['Mag =  ', num2str(mag/1e6), ' M'], ['pix =  ', num2str(pix)], ['resol = ',num2str(res), ' nm'], ['proj =  ', num2str(nrset)], ...
                ['rep =  ', num2str(nrep)],  ['Nr defocus =  ', num2str(ndefocus)], 'Scale bars = 5nm' };
            annotation('textbox', [.15 .8 .1 .1], 'String', txt);
            subplot(3,2,2)
            plot(da_valid, 'r-')
            hold on
            plot(dg_valid, 'b-')
            legend('\Delta \alpha', '\Delta \gamma');
            xlabel('projection number')
            ylabel('delta angle/degrees')
            title('Difference between actual and intended angles')
            hold off
            subplot(3,2,3)
            imshow(obj.addscalebar(image, res,5),[])
            title('Image at phi = 0')
            
            subplot(3,2,4)
            imshow(obj.addscalebar(image90, res,5),[])
            
            title('Closest image to phi = 90')
            
            pst = round(pix/4);
            pen = round(pix*3/4);
            
            subplot(3,2,5)
            imshow(obj.addscalebar(image(pst:pen,pst:pen), res,5),[])
            title('Image at phi = 0 (zoom)')
            subplot(3,2,6)
            imshow(obj.addscalebar(image90(pst:pen,pst:pen), res,5),[])
            title('Closest image to phi = 90 (zoom)')
            
            %set(0,'Interpreter', 'tex')
            set(1, 'PaperType', 'A4');
            print([emd,'sum.png'],'-dpng','-r200')
            
            %%
        end
    end
end
