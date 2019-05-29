clc
close all
clearvars

%% Visualization of the detection result from NPD

% load data
data = hzhu_read_all([],'detail');
data_x = hzhu_read_all([],'data_x');
data_y = hzhu_read_all([],'data_y');

data_cell = struct2cell(data);
data_cell_x = struct2cell(data_x);
data_cell_y = struct2cell(data_y);

t = -15:1:35;
dm = 0.95;

sampling_rate = 1000;

%% Plot results

N = data.n;

for i = 2:N+1
    data_cell{i}(:,3) = data_cell{i}(:,3)+1;
    data_cell{i}(:,7) = data_cell{i}(:,7)+1;
    figure
    subplot(2,1,1)
    hold all    
    leg = [];
    for j = 1:length(data_cell{i}(:,1))
        plot(data_cell{i}(j,3)+t,...
            h_s(t,data_cell{i}(j,4),data_cell{i}(j,5),data_cell{i}(j,6),dm)+data_cell_x{i}(data_cell{i}(j,3)),...
            'linewidth',1.5);
        leg{j} = ['saccadic event ',num2str(j)];
    end
    scatter(data_cell{i}(:,3),data_cell_x{i}(data_cell{i}(:,3)+1),30,'ro','LineWidth',0.8)
    
    theta1_x = data_cell{i}(:,4);
    theta2_x = data_cell{i}(:,5);
    theta3_x = data_cell{i}(:,6);
    
    tmp = (theta3_x-1)./theta3_x;
    d = theta2_x.*(-log(0.5-0.5*dm)).^(1./theta3_x)-theta2_x.*(-log(0.5+0.5*dm)).^(1./theta3_x);
    
    vmax_x = theta1_x.*theta3_x./theta2_x.*exp(-tmp).*tmp.^tmp;
    t1 = theta2_x.*(-log(0.5+0.5*dm)).^(1./theta3_x);
    tv = theta2_x.*tmp.^(1.0./theta3_x);
    
    scatter(round(data_cell{i}(:,3)+d),data_cell_x{i}(round(data_cell{i}(:,3)+d)),50,'bx','LineWidth',0.8)
    
    scatter(round(data_cell{i}(:,3)+tv-t1),data_cell_x{i}(round(data_cell{i}(:,3)+tv-t1)),50,'*','LineWidth',0.8)
    
    for p = 1:length(vmax_x)
        text(round(data_cell{i}(p,3)+tv(p)-t1(p)),data_cell_x{i}(round(data_cell{i}(p,3)+tv(p)-t1(p))),...
            ['\leftarrow ',num2str(vmax_x(p)*sampling_rate,'%.1f'),'^o/s'])
    end
    
    plot(data_cell_x{i})
    
    leg{j+1} = 'start point';
    leg{j+2} = 'end point';
    leg{j+3} = 'max velocity';
    leg{j+4} = 'gaze signal (X)';
    legend(leg,'Location','best')
    xlabel('Discrete time index k')
    ylabel('X position (degree)')
    
    subplot(2,1,2)
    hold all
    
    leg = [];
    for j = 1:length(data_cell{i}(:,1))
        plot(data_cell{i}(j,7)+t,...
            h_s(t,data_cell{i}(j,8),data_cell{i}(j,9),data_cell{i}(j,10),dm)+data_cell_y{i}(data_cell{i}(j,7)),...
            'linewidth',2);
        leg{j} = ['saccadic event ',num2str(j)];
    end
    scatter(data_cell{i}(:,7),data_cell_y{i}(data_cell{i}(:,7)),30,'ro','LineWidth',0.8)
    
    theta1_y = data_cell{i}(:,8);
    theta2_y = data_cell{i}(:,9);
    theta3_y = data_cell{i}(:,10);
    
    tmp = (theta3_y-1)./theta3_y;
    d = theta2_y.*(-log(0.5-0.5*dm)).^(1./theta3_y)-theta2_y.*(-log(0.5+0.5*dm)).^(1./theta3_y);
    
    vmax_y = theta1_y.*theta3_y./theta2_y.*exp(-tmp).*tmp.^tmp;
    t1 = theta2_y.*(-log(0.5+0.5*dm)).^(1./theta3_y);
    tv = theta2_y.*tmp.^(1.0./theta3_y);
    
    scatter(round(data_cell{i}(:,7)+d),data_cell_y{i}(round(data_cell{i}(:,7)+d)),50,'bx','LineWidth',0.8)
    
    scatter(round(data_cell{i}(:,7)+tv-t1),data_cell_y{i}(round(data_cell{i}(:,7)+tv-t1)),50,'*','LineWidth',0.8)
    
    for p = 1:length(vmax_y)
        text(round(data_cell{i}(p,7)+tv(p)-t1(p)),data_cell_y{i}(round(data_cell{i}(p,7)+tv(p)-t1(p))),...
            ['\leftarrow ',num2str(vmax_y(p)*1000,'%.1f'),'^o/s'])
    end
    
    plot(data_cell_y{i})

    leg{j+1} = 'start point';
    leg{j+2} = 'end point';
    leg{j+3} = 'max velocity';
    leg{j+4} = 'gaze signal (Y)';
    legend(leg,'Location','best')
    xlabel('Discrete time index k')
    ylabel('Y position (degree)')
    
end