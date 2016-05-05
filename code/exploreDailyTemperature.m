%% clean up

clc
clear all
close all

addpath ~/Documents/MATLABLibs/regression/

%% locations and files

projectDir = '/Users/rwhut/Documents/TU/MATLAB_work/graceCERN';

dataDir = [projectDir filesep 'data'];
figDir = [projectDir filesep 'fig'];

cernDataFile1 = [dataDir filesep 'circRing.2012.csv'];
cernDataFile2 = [dataDir filesep 'circRing.2015.csv'];
meteoData = [dataDir filesep 'SZ000008440.dly.txt'];

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


%meteo data from Geneva airport
dataMeteo = textread(meteoData,'%s%*[^\n]','whitespace','none');
TAVG=[];
TAVGTimestamp=[];
PRCP=[];
PRCPTimestamp=[];
for dataMeteoLine = 1 : length(dataMeteo)
    if ((str2num(dataMeteo{dataMeteoLine}(12:15))==2012) | ...
            (str2num(dataMeteo{dataMeteoLine}(12:15))==2015))
        if strcmp(dataMeteo{dataMeteoLine}(18:21),'TAVG')
            for dayCounter = 1 : ((length(dataMeteo{dataMeteoLine})-21)/8)
                if (str2num(dataMeteo{dataMeteoLine}(((dayCounter-1)*8)+(1:5)+21))~=-9999)
                    TAVG=[TAVG,str2num(dataMeteo{dataMeteoLine}(((dayCounter-1)*8)+(1:5)+21))/10];
                    TAVGTimestamp=[TAVGTimestamp,datenum(str2num(dataMeteo{dataMeteoLine}(12:15)),...
                        str2num(dataMeteo{dataMeteoLine}(16:17)),dayCounter,12,0,0)];
                end
            end
        elseif strcmp(dataMeteo{dataMeteoLine}(18:21),'PRCP')
            for dayCounter = 1 : ((length(dataMeteo{dataMeteoLine})-21)/8)
                if (str2num(dataMeteo{dataMeteoLine}(((dayCounter-1)*8)+(1:5)+21))~=-9999)
                    PRCP=[PRCP,str2num(dataMeteo{dataMeteoLine}(((dayCounter-1)*8)+(1:5)+21))/10];
                    PRCPTimestamp=[PRCPTimestamp,datenum(str2num(dataMeteo{dataMeteoLine}(12:15)),...
                        str2num(dataMeteo{dataMeteoLine}(16:17)),dayCounter,12,0,0)];
                end
            end        
        end
    end
end %for dataMeteoLine = 1 : length(dataMeteo)



clear dataCern1 dataCern2 dataGrace fID;



%% initial plots for exploration of data

figure(1);
[ax,h1,h2]=plotyy(cernTimeStamp,cernCircumData,TAVGTimestamp,TAVG);
set(h1,'LineStyle','none');
set(h1,'Marker','.')
set(h2,'LineStyle','none');
set(h2,'Marker','X')
labelLinks=get(ax(1),'YLabel');
set(labelLinks,'String','change in LHC circumference (monthly average)');
labelRechts=get(ax(2),'YLabel');
set(labelRechts,'String','Temperature (C)');
set(gca,'XTickLabel',datestr(get(gca,'XTick'),'dd-mmm-yyyy'))
print(gcf,[figDir filesep 'rawDataTemperature.eps'],'-depsc');
print(gcf,[figDir filesep 'rawDataTemperature.jpg'],'-djpeg99');

%need to constrain GRACE data to CERN time domain.: also: doublecheck
%timezones.


%make raw data plots of the "cern" seasons
figure(2)
subplot(1,2,1)
[ax,h1,h2]=plotyy(cernTimeStamp,cernCircumData,TAVGTimestamp,TAVG);
set(h1,'LineStyle','none');
set(h1,'Marker','.')
set(h2,'LineStyle','none');
set(h2,'Marker','X')
set(ax,'XLim',[datenum(2012,01,01) datenum(2012,12,31)]);
%set(ax(2),'YLim',[-200 100]);
labelLinks=get(ax(1),'YLabel');
set(labelLinks,'String','change in LHC circumference');
set(gca,'XTickLabel',datestr(get(gca,'XTick'),'dd-mmm'))

subplot(1,2,2)
[ax,h1,h2]=plotyy(cernTimeStamp,cernCircumData,TAVGTimestamp,TAVG);
set(h1,'LineStyle','none');
set(h1,'Marker','.')
set(h2,'LineStyle','none');
set(h2,'Marker','X')
set(ax,'XLim',[datenum(2015,01,01) datenum(2015,12,31)]);
%set(ax(2),'YLim',[-200 100]);
labelRechts=get(ax(2),'YLabel');
set(labelRechts,'String','Temperature (C)');
set(gca,'XTickLabel',datestr(get(gca,'XTick'),'dd-mmm'))

print(gcf,[figDir filesep 'Cern20122015DATADailyTemperature.eps'],'-depsc');
print(gcf,[figDir filesep 'Cern20122015DATADailyTemperature.jpg'],'-djpeg99');


%% data manipulation

%since we have daily GRACE data, but not daily CERN data, we will look for
%the clostest GRACE point to every CERN point.

cernOnTAVGTimeStamp=interp1(TAVGTimestamp,TAVGTimestamp,cernTimeStamp,'nearest');
cernFilter = not(isnan(cernOnTAVGTimeStamp));

%check: show the percentage of data not connected
disp('percentage of cern data used:');
disp(sum(cernFilter)/length(cernFilter));

cernCircumDataOnTAVGTimeStamp=cernCircumData(cernFilter);
cernOnTAVGTimeStamp=cernOnTAVGTimeStamp(cernFilter);

TAVGFilter = ismember(TAVGTimestamp,cernOnTAVGTimeStamp);
%check: show the percentage of data not connected
disp('percentage of grace data used:');
disp(sum(TAVGFilter)/length(TAVGFilter));

TAVGFiltered=TAVG(TAVGFilter);
TAVGTimeStampFiltered=TAVGTimestamp(TAVGFilter);


cernAvgData=[];
cernStdData=[];
cernDataBuffer=cernCircumDataOnTAVGTimeStamp(1);


for cernDataCounter = 2:length(cernOnTAVGTimeStamp);
    if cernOnTAVGTimeStamp(cernDataCounter)==cernOnTAVGTimeStamp(cernDataCounter-1)
        cernDataBuffer=[cernDataBuffer;cernCircumDataOnTAVGTimeStamp(cernDataCounter)];
        
    else
        cernAvgData=[cernAvgData;mean(cernDataBuffer)];
        cernStdData=[cernStdData;std(cernDataBuffer)];
        cernDataBuffer=[cernCircumDataOnTAVGTimeStamp(cernDataCounter)];
    end %if cernOnGraceTimeStamp(cernDataCounter)==cernOnGraceTimeStamp(cernDataCounter-1)
end %for cernDataCounter = 1:length(cernOnGraceTimeStamp);

cernAvgData=[cernAvgData;mean(cernDataBuffer)];
cernStdData=[cernStdData;std(cernDataBuffer)];

%[offset,helling]=kkmfit(graceWaterDataFiltered,cernAvgData,cernStdData);
[offset,helling]=kkmfit(TAVGFiltered,cernAvgData);
yHat=offset+TAVGFiltered'*helling;
R2=1-(mean((yHat-cernAvgData).^2)/mean(cernAvgData.^2));

figure(3);
H=scatter(TAVGFiltered,cernAvgData);
%set(H,'LineStyle','none')
set(H,'Marker','X')
hold on
plot(TAVGFiltered,yHat,'r')
xlabel('Temperature at Geneva Airport (C)');
ylabel('change in LHC circumference (daily average in mm)');
text(3,0.45,['R^{2}: ' num2str(R2*100,2) '%']);

print(gcf,[figDir filesep 'cernVSGraceDailyTemperature.eps'],'-depsc');
print(gcf,[figDir filesep 'cernVSGraceDailyTemperature.jpg'],'-djpeg99');


