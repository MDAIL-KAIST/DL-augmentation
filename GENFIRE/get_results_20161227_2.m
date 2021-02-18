%% AllRanJutt5_RandStart_R_RealInterp
%NamePrefixString = 'RefineResult_FePt_Jit5_RS';
%NamePrefixString = 'RefineResult_FePt_Jit5MAInt1_RS';
% NamePrefixString = 'RefineResult_FePt_Jit5MA1_RS';
% NamePrefixString = 'RefineResult_FePt_Jit5MA2_RS';
% NamePrefixString = 'RefineResult_FePt_Jit5MAInt2_RS';
 NamePrefixString = 'RefineResult_FePt_Jit10_RS';



dispFigure = 1;
saveImage = 1;


NameCell = {'Zero','Puturb1','Puturb2','Puturb3'};
angle_filename   = '../data/Angles_AllRanJitt5_RI.mat';
given_angles = importdata(angle_filename);

 MeanErrPhi = zeros(1,length(NameCell));
 MeanErrThe = zeros(1,length(NameCell));
 MeanErrPsi = zeros(1,length(NameCell));
for Loopnum = 1:length(NameCell)

resultFileName = sprintf('./results/%s_%s.mat',NamePrefixString,NameCell{Loopnum});

RR = importdata(resultFileName);

 AngleEvolution = RR.REFINEMENT.AngleEvolution;
 
 true_angles = RR.true_angles;
    true_angles = reorient_Angles(true_angles,RR.REFINEMENT.RefineReferenceAngleInd,RR.REFINEMENT.RefineZeroCenterFlag,RR.REFINEMENT.RefineReferenceAngletoSet);

 refined_angles = RR.refined_angles;

 MeanErrPhi(Loopnum) = mean(abs(true_angles(:,1)-refined_angles(:,1)));
 MeanErrThe(Loopnum) = mean(abs(true_angles(:,2)-refined_angles(:,2)));
 MeanErrPsi(Loopnum) = mean(abs(true_angles(:,3)-refined_angles(:,3)));

 if dispFigure
figure(1)
subplot(1,3,1)
plot(true_angles(:,1),'bo');hold all
plot(refined_angles(:,1),'r.');hold off
title(sprintf('phi, avgErr = %5.2f', MeanErrPhi(Loopnum)));
legend('true','refined')

subplot(1,3,2)
plot(true_angles(:,3),'bo');hold all
plot(refined_angles(:,3),'r.');hold off
title(sprintf('psi, avgErr = %5.2f', MeanErrThe(Loopnum)));
legend('true','refined')

subplot(1,3,3)
plot(true_angles(:,2)-given_angles(:,2),'r-');hold on
plot(true_angles(:,2)-refined_angles(:,2),'b-');hold off
title(sprintf('theta, avgErr = %5.2f', MeanErrPsi(Loopnum)));
legend('given','refined')

if saveImage
print('-dpng','-r150',sprintf('./results/pngs/%s_Fig1_Loop%d.png',NamePrefixString,Loopnum));
end

figure(2)
Legend = {};
for i=1:size(AngleEvolution,3)
    Legend{end+1} = sprintf('Iter %d',i);
end

NameTagCell = {'phi','theta','psi'};
for i=1:3
    subplot(1,3,i)
    for j=1:size(AngleEvolution,3)
      plot(squeeze(AngleEvolution(:,i,j)));hold all;%pause
    end
    hold off;
    title(sprintf('%s Evolution',NameTagCell{i}));
    legend(Legend{:})
end

if saveImage
print('-dpng','-r150',sprintf('./results/pngs/%s_Fig2_Loop%d.png',NamePrefixString,Loopnum));
end

 end
%pause
end
% 
% mean(MeanErrPhi([1 2 4]))
% mean(MeanErrThe([1 2 4]))
% mean(MeanErrPsi([1 2 4]))

