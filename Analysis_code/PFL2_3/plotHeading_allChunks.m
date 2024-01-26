clear
rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\raw_PFL2_3'; % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp
 % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp

date = input('Date? ','s');
cell_num = input('Cell? ','s');
cell_num = strcat('cell_',cell_num);
trial = input('Trial? ','s');
trial = strcat('trial_',trial);
fileName = fullfile(rootPath,date,cell_num,trial);

cd(fileName)
load('pro_trialData.mat');
load('pro_behaviourData.mat');
load('trialMeta.mat');

cd('C:\Code\EphysCode'); 

if iscell(processed_behaviourData)
    processed_behaviourData = processed_behaviourData{1};
end

if iscell(processed_trialData)
    processed_trialData = processed_trialData{1};
end
%%
%% separate idx at which fly is moving in CL vs OL
    
    startOL = 178.189*1000; % sec
    endOL = 203.249*1000;
    lengthCL = 179.99*1000;
    
% start = 1819.45;
% finish = 1844.47;
%     startOL  = find(processed_behaviourData.time == start); 
%     endOL = find(processed_behaviourData.time == finish); 

[CL, OL, CL_startStopIdx, OL_startStopIdx] = separateCL_OL(startOL, endOL, lengthCL, processed_behaviourData, processed_trialData);
%% CL_OL trial split into 1 min segements
% split 3 min continuous CL chunks into 1 min segments for more datapoints
count = 1;
CL_chunks = [];
for chunk = 1:length(CL_startStopIdx)
    distance = CL_startStopIdx(chunk,2) - CL_startStopIdx(chunk,1);
    step = floor(distance/3);
    l = length(CL_chunks); 
    CL_chunks(l+1,1) = CL_startStopIdx(chunk,1); 
    CL_chunks(l+1,2) = CL_startStopIdx(chunk,1) + (step);
    
    CL_chunks(l+2,1) = CL_startStopIdx(chunk,1) + (step)+1; 
    CL_chunks(l+2,2) = CL_startStopIdx(chunk,1) + (step)*2;
    
    CL_chunks(l+3,1) = CL_startStopIdx(chunk,1) + (step)*2+1; 
    CL_chunks(l+3,2) = CL_startStopIdx(chunk,2);
end

count = 1; 

%%

for chunk = 1:length(CL_chunks)
 
    start = CL_chunks(chunk, 1);
    finish = CL_chunks(chunk, 2);
    
    f = fieldnames(processed_behaviourData);
    for k=1:numel(f)
        chunk_behaviourData.(f{k}) = processed_behaviourData.(f{k})(start:finish); 
    end


    f = fieldnames(processed_trialData);
    for k=1:numel(f)
        chunk_trialData.(f{k}) = processed_trialData.(f{k})(start:finish); 
    end
    
    activity_values = chunk_trialData.smooth_Vm';
    x_values = chunk_behaviourData.angle';
    y_values = chunk_behaviourData.vel_for';

    

    x_edges = [-180:10:180]; %
    y_edges = [-2:.5:10];
    
    [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);

    fig = figure(1); clf;
           
            imagesc(flip(heatmap_Vm))
            ylabel('Vf mm/s')
            xt = get(gca, 'XTick');
            xtnew = linspace(min(xt), max(xt), 7);                             
            xtlbl = linspace(-135, 165, numel(xtnew));                  
            set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
            xlabel('angle')
            colorbar
            title('Vm')
            yt = get(gca, 'YTick');
            ytnew = linspace(min(yt),max(yt),6);
            ytlbl = linspace(8.5, -1.5, numel(ytnew));
            set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl)
            set(gcf,'color','w')
            
       fileName = ['Z:\Dropbox (HMS)\Wilson_Lab_Data\Analysis\labMeeting102021\heatmaps_070121\','map_',num2str(chunk),'.jpeg'];  
        saveas(fig,fileName)
end 
    