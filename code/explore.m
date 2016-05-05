%% clean up

clc
clear all
close all

%% locations and files

projectDir = '/Users/rwhut/Documents/TU/MATLAB_work/graceCERN';

dataDir = [projectDir filesep 'data'];
figDir = [projectDir filesep 'fig'];

cernDataFile1 = [dataDir filesep 'circRing.2012.csv'];
cernDataFile2 = [dataDir filesep 'circRing.2015.csv'];
graceData = [dataDir filesep 'graceCERN25_100.txt'];

%% read files

% circumference data from CERN
fID = fopen(cernDataFile1);
dataCern1 = textscan(fID,'%u %q %f %f','Delimiter',',','headerLines',2); %%{dd-mm-yy_HH:MM:SS}D
fclose(fID);
fID = fopen(cernDataFile2);
dataCern2 = textscan(fID,'%u %q %f %f','Delimiter',',','headerLines',2); %%{dd-mm-yy_HH:MM:SS}D
fclose(fID);

cernTimeStamp=datenum([dataCern1{2};dataCern2{2}],'dd-mm-yy_HH:MM:SS');
cernCircumData=[dataCern1{3};dataCern2{3}];


%gravity data from 
fID = fopen(graceData);
dataGrace = textscan(fID,'%f %f %f %f %f','headerLines',7);
fclose(fID);
graceTimeStamp = datenum([dataGrace{1} zeros(size(dataGrace{1})) zeros(size(dataGrace{1})) ...
zeros(size(dataGrace{1})) zeros(size(dataGrace{1})) zeros(size(dataGrace{1}))]);
graceWaterData = dataGrace{2};

clear dataCern1 dataCern2 dataGrace fID;



%% initial plots for exploration of data

figure(1);
[ax,h1,h2]=plotyy(cernTimeStamp,cernCircumData,graceTimeStamp,graceWaterData);
set(h1,'LineStyle','none');
set(h1,'Marker','.')
set(h1,'LineStyle','none');
set(h2,'Marker','X')
labelLinks=get(ax(1),'YLabel');
set(labelLinks,'String','change in LHC circumference (monthly average)');
labelRechts=get(ax(2),'YLabel');
set(labelRechts,'String','monthly change in stored water (mm) from GRACE satallite');
set(gca,'XTickLabel',datestr(get(gca,'XTick'),'dd-mmm-yyyy'))
print(gcf,[figDir filesep 'rawData.eps'],'-depsc');
print(gcf,[figDir filesep 'rawData.jpg'],'-djpeg99');

%need to constrain GRACE data to CERN time domain.: also: doublecheck
%timezones.


%make raw data plots of the "cern" seasons
figure(2)
subplot(1,2,1)
[ax,h1,h2]=plotyy(cernTimeStamp,cernCircumData,graceTimeStamp,graceWaterData);
set(h1,'LineStyle','none');
set(h1,'Marker','.')
set(h1,'LineStyle','none');
set(h2,'Marker','X')
set(ax,'XLim',[datenum(2012,01,01) datenum(2012,12,31)]);
set(ax(2),'YLim',[-200 100]);
labelLinks=get(ax(1),'YLabel');
set(labelLinks,'String','change in LHC circumference (monthly average)');
set(gca,'XTickLabel',datestr(get(gca,'XTick'),'dd-mmm'))

subplot(1,2,2)
[ax,h1,h2]=plotyy(cernTimeStamp,cernCircumData,graceTimeStamp,graceWaterData);
set(h1,'LineStyle','none');
set(h1,'Marker','.')
set(h1,'LineStyle','none');
set(h2,'Marker','X')
set(ax,'XLim',[datenum(2015,01,01) datenum(2015,12,31)]);
set(ax(2),'YLim',[-200 100]);
labelRechts=get(ax(2),'YLabel');
set(labelRechts,'String','monthly change in stored water (mm) from GRACE satallite');
set(gca,'XTickLabel',datestr(get(gca,'XTick'),'dd-mmm'))

print(gcf,[figDir filesep 'Cern20122015DATA.eps'],'-depsc');
print(gcf,[figDir filesep 'Cern20122015DATA.jpg'],'-djpeg99');


%% data manipulation

%assuming that the GRACE data has broader support than the CERN data, but
%that the CERN data has more numerous we interpolate the CERN data to the
%GRACE timestamps
cernOnGraceTimeStamp=interp1(graceTimeStamp,graceTimeStamp,cernTimeStamp,'nearest');
cernFilter = not(isnan(cernOnGraceTimeStamp));

%check: show the percentage of data not connected
disp('percentage of cern data used:');
disp(sum(cernFilter)/length(cernFilter));

cernCircumDataOnGraceTimeStamp=cernCircumData(cernFilter);
cernOnGraceTimeStamp=cernOnGraceTimeStamp(cernFilter);

graceFilter = ismember(graceTimeStamp,cernOnGraceTimeStamp);
%check: show the percentage of data not connected
disp('percentage of grace data used:');
disp(sum(graceFilter)/length(graceFilter));

graceWaterDataFiltered=graceWaterData(graceFilter);
graceTimeStampFiltered=graceTimeStamp(graceFilter);


cernAvgData=[];
cernStdData=[];
cernDataBuffer=cernCircumDataOnGraceTimeStamp(1);


for cernDataCounter = 2:length(cernOnGraceTimeStamp);
    if cernOnGraceTimeStamp(cernDataCounter)==cernOnGraceTimeStamp(cernDataCounter-1)
        cernDataBuffer=[cernDataBuffer;cernCircumDataOnGraceTimeStamp(cernDataCounter)];
        
    else
        cernAvgData=[cernAvgData;mean(cernDataBuffer)];
        cernStdData=[cernStdData;std(cernDataBuffer)];
        cernDataBuffer=[];
    end %if cernOnGraceTimeStamp(cernDataCounter)==cernOnGraceTimeStamp(cernDataCounter-1)
end %for cernDataCounter = 1:length(cernOnGraceTimeStamp);

cernAvgData=[cernAvgData;mean(cernDataBuffer)];
cernStdData=[cernStdData;std(cernDataBuffer)];

[offset,helling]=kkmfit(graceWaterDataFiltered,cernAvgData,cernStdData);
yHat=offset+graceWaterDataFiltered*helling;
R2=1-(mean((yHat-cernAvgData).^2)/mean(cernAvgData.^2));

figure(3);
H=errorbar(graceWaterDataFiltered,cernAvgData,cernStdData);
set(H,'LineStyle','none')
set(H,'Marker','X')
hold on
plot(graceWaterDataFiltered,yHat,'r')
xlabel('monthly change in stored water (mm) from GRACE satallite');
ylabel('change in LHC circumference (monthly average)');
text(-60,0.9,['R^{2}: ' num2str(R2*100,2) '%']);

print(gcf,[figDir filesep 'cernVSGrace.eps'],'-depsc');
print(gcf,[figDir filesep 'cernVSGrace.jpg'],'-djpeg99');
