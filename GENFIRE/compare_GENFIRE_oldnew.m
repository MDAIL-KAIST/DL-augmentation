for A=1:2
  for Pos=0:1
    for Sup=0:1
      for Moddd=1:3
        for Gridd=1:2
          
            filename_newResults = sprintf('./results/GENFIRE_rec_new_NRF_%d_%d_%d_%d_%d.mat',A,Pos,Sup,Moddd,Gridd);
            filename_oldResults = sprintf('./results/GENFIRE_rec_old_NRF_%d_%d_%d_%d_%d.mat',A,Pos,Sup,Moddd,Gridd);
            
            NEW = importdata(filename_newResults);
            OLD = importdata(filename_oldResults);
            
            %pause
            NEW_rec = NEW.reconstruction;
            OLD_rec = OLD.reconstruction;
            
            [FSCout, spatialFrequencyout] = FourierShellCorrelate(NEW_rec, OLD_rec, 25, 0.5); 
            plot(spatialFrequencyout,FSCout);pause
            
%             subplot(2,3,1)
%             imagesc(squeeze(sum(NEW_rec,1)));axis image
%             subplot(2,3,2)
%             imagesc(squeeze(sum(NEW_rec,2)));axis image
%             subplot(2,3,3)
%             imagesc(squeeze(sum(NEW_rec,3)));axis image
%             subplot(2,3,4)
%             imagesc(squeeze(sum(OLD_rec,1)));axis image
%             subplot(2,3,5)
%             imagesc(squeeze(sum(OLD_rec,2)));axis image
%             subplot(2,3,6)
%             imagesc(squeeze(sum(OLD_rec,3)));axis image
%             pause


            %B=sum(abs(NEW_rec(:)-OLD_rec(:))>1e-3)
            if B==0
%               A
%               Pos
%               Sup
%               Moddd
%              Gridd
            end
              
            
        end
      end
    end
  end
end